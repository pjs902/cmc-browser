import glob
import requests
import gzip
import tarfile
import os
import time
import h5py
import numpy as np
import pandas as pd

import cmctoolkit as ck

__all__ = ["CMCBrowser"]


class CMCBrowser:
    def __init__(self, ss_dir="/home/peter/research/CMC-obs/snapshots"):
        """
        Create a CMCBrowser object to load and interact with CMC snapshots.

        Parameters
        ----------
        ss_dir : str
            The directory containing the CMC snapshots.
        """

        self.ss_dir = ss_dir
        self.models_list = glob.glob(ss_dir + "/*.*")
        # remove directory from model name
        self.models_list = [model.split("/")[-1] for model in self.models_list]

        # initialize empty dictionary to store loaded snapshots
        self.loaded_snapshots = {}

        # enumerate the snapshots for each model
        self.model_snapshots = {}
        for model in self.models_list:
            # check first if there is a tar.gz snapshot
            snap_fn = glob.glob(f"{self.ss_dir}/{model}/initial.*.dat.gz")
            if len(snap_fn) == 0:
                # check instead for the hdf5 snapshot
                snap_fn = glob.glob(f"{self.ss_dir}/{model}/*.window.snapshots.h5")

                # enumerate the datasets in the hdf5 file
                with h5py.File(snap_fn[0], "r") as f:
                    snap_fn = [ss for ss in f.keys()]

            # save the snapshot filenames
            self.model_snapshots[model] = snap_fn

        # remove directory from snapshot name
        for model in self.models_list:
            self.model_snapshots[model] = [
                ss.split("/")[-1] for ss in self.model_snapshots[model]
            ]

    def list_models(self):
        """
        Print a list of all models in snapshot directory.
        """
        print("Found the following models:")
        for model in self.models_list:
            print(model)

    def list_snapshots(self, model_name):
        """
        Print a list of all snapshots for a given model.
        """
        print("Found the following snapshots:")
        for ss in self.model_snapshots[model_name]:
            print(ss)

    def load_snapshot(
        self,
        model_name,
        distance=5.0,
        ss_name="king.window.snapshots.h5",
        mode="h5",
        *,
        h5_key=None,
    ):
        """
        Load a snapshot into the CMCBrowser object.

        Parameters
        ----------
        model_name : str
            The name of the model to load.
        distance : float
            The distance to use for the snapshot.
        ss_name : str
            The name of the snapshot to load. Either a hdf5 file or an dat.gz file.
        mode : str
            The mode to use for loading the snapshot. Must be one of ["h5", "dat.gz"].
        h5_key : str
            The key to use for the hdf5 snapshot. If None, will use the last key in the file. Only used if mode is "h5".
        """

        # parse the output prefix from the snapshot name
        prefix = ss_name.split(".")[0]

        # need to parse metallicity
        Z = float(model_name.split("Z")[-1])

        print(
            f"Loading Model: {model_name}, snapshot: {ss_name}, Z={Z}, distance={distance}"
        )

        if mode == "dat.gz":
            snap = ck.Snapshot(
                fname=f"{self.ss_dir}/{model_name}/{ss_name}",
                conv=f"{self.ss_dir}/{model_name}/{prefix}.conv.sh",
                z=Z,
                dist=distance,
            )
            # compatibility with old column names
            snap.data["m_MSUN"] = snap.data["m[MSUN]"]
            snap.data["m0_MSUN"] = snap.data["m0[MSUN]"]
            snap.data["m1_MSUN"] = snap.data["m1[MSUN]"]
            snap.name = f"{model_name}/{ss_name}"
        elif mode == "h5":
            snap = ck.Snapshot(
                fname=f"{self.ss_dir}/{model_name}/{ss_name}",
                conv=f"{self.ss_dir}/{model_name}/{prefix}.conv.sh",
                snapshot_name=h5_key,
                z=Z,
                dist=distance,
            )
            # create a fresh copy of the data to avoid memory fragmentation
            snap.data = snap.data.copy()

            # create a new column with the old mass column name for compatibility
            snap.data["m[MSUN]"] = snap.data["m_MSUN"]
            snap.data["m0[MSUN]"] = snap.data["m0_MSUN"]
            snap.data["m1[MSUN]"] = snap.data["m1_MSUN"]
            snap.name = f"{model_name}/{h5_key}"

        # add some useful info
        snap.mass = snap.data["m_MSUN"].sum()
        snap.FeH = np.log10(snap.z / 0.02)

        # calculate the central escape velocity, using the gravitational potential at the center of the cluster and the tidal boundary
        log_file = f"{self.ss_dir}/{model_name}/{prefix}.esc.dat"

        # load the log file with pandas, delim is just a space
        esc = pd.read_csv(log_file, delim_whitespace=True)

        esc["t[Myr]"] = esc["#2:t"] * snap.unitdict["myr"]

        # conversions to physical units, following https://github.com/tomas-cabrera/hvss-bsco/blob/main/src/scripts/cmc_single_clusters_vesc.py
        nb_kms = 1e-5 * snap.unitdict["cm"] / snap.unitdict["nb_s"]
        esc["#10:phi_rtidal"] *= nb_kms**2
        esc["#11:phi_zero"] *= nb_kms**2

        esc["vesc"] = np.sqrt(2 * (esc["#10:phi_rtidal"] - esc["#11:phi_zero"]))

        snap.vesc_initial = esc[esc["t[Myr]"] < 20]["vesc"].mean()
        snap.vesc_final = esc["vesc"][-5000:].mean()
        del esc

        # BH info
        snap.bh_masses = pd.concat(
            [
                snap.data.loc[(snap.data["startype"] == 14)]["m_MSUN"],
                snap.data.loc[(snap.data["bin_startype0"] == 14)]["m0_MSUN"],
                snap.data.loc[(snap.data["bin_startype1"] == 14)]["m1_MSUN"],
            ],
            axis=0,
        ).to_list()
        snap.M_BH = np.sum(snap.bh_masses)
        snap.N_BH = len(snap.bh_masses)

        snap.bh_radii = pd.concat(
            [
                snap.data.loc[(snap.data["startype"] == 14)]["d[PC]"],
                snap.data.loc[(snap.data["bin_startype0"] == 14)]["d[PC]"],
                snap.data.loc[(snap.data["bin_startype1"] == 14)]["d[PC]"],
            ],
            axis=0,
        ).to_list()

        # TODO: could eventually add this sort of stuff for WDs or anything else we're interested in

        # half mass radius
        snap.rh = snap.calculate_renclosed(enclosed_frac=0.5, qty="mass")

        if mode == "dat.gz":
            self.loaded_snapshots[f"{model_name}/{ss_name}"] = snap
        elif mode == "h5":
            self.loaded_snapshots[f"{model_name}/{h5_key}"] = snap

    def download_new_model(self, N, rv, rg, Z):  # pragma: no cover
        """
        Download new models from the CMC website. Will not download if the model already exists.

        Parameters
        ----------
        N : str
            The number of particles in the model. Must be one of ["1.6e6", "2e5", "3.2e6", "4e5", "8e5"]

        rv : str
            The virial radius of the model. Must be one of ["0.5", "1", "2", "4"]

        rg : str
            The galactocentric radius of the model. Must be one of ["2", "20", "8"]

        Z : str
            The metallicity of the model. Must be one of ["0.02", "0.002", "0.0002"]

        Notes
        -----
        The largest models with N=3.2e6 are only computed for a subset of the possible parameters.

        """
        base_url = (
            "https://cmc.ciera.northwestern.edu/download-cluster/download-cluster/"
        )

        possible_Ns = ["1.6e6", "2e5", "3.2e6", "4e5", "8e5"]
        possible_rvs = ["0.5", "1", "2", "4"]
        possible_rgs = ["2", "20", "8"]
        possible_Zs = ["0.02", "0.002", "0.0002"]

        if str(N) not in possible_Ns:
            raise ValueError(f"Invalid N: {N}, must be one of {possible_Ns}")
        if str(rv) not in possible_rvs:
            raise ValueError(f"Invalid rv: {rv}, must be one of {possible_rvs}")
        if str(rg) not in possible_rgs:
            raise ValueError(f"Invalid rg: {rg}, must be one of {possible_rgs}")
        if str(Z) not in possible_Zs:
            raise ValueError(f"Invalid Z: {Z}, must be one of {possible_Zs}")

        # check if model already exists
        model_name = f"N{str(N)}_rv{str(rv)}_rg{str(rg)}_Z{str(Z)}"
        if model_name in self.models_list:
            raise ValueError(f"Model {model_name} already exists!")

        # record time
        start_time = time.time()

        # construct the url
        url = f"{base_url}?number_of_objects=N{str(N)}&virial_radius=rv{str(rv)}&galactocentric_distance=rg{str(rg)}&metallicity=Z{str(Z)}"

        # construct the model name
        model_name = f"N{str(N)}_rv{str(rv)}_rg{str(rg)}_Z{str(Z)}"

        # download the file
        print(f"Downloading model: {model_name}")
        r = requests.get(url, allow_redirects=True, timeout=10, verify=False)

        # check response
        if r.status_code != 200:
            raise ValueError(f"Download failed with status code: {r.status_code}")

        # write the file
        with open(f"{self.ss_dir}/{model_name}.tar.gz", "wb") as f:
            f.write(r.content)

        # extract the file
        print(f"Unzipping model: {model_name}")
        with gzip.open(f"{self.ss_dir}/{model_name}.tar.gz", "rb") as f_in:
            with open(f"{self.ss_dir}/{model_name}.tar", "wb") as f_out:
                f_out.write(f_in.read())

        # untar the file into the subdirectory
        print(f"Untarring model: {model_name}")
        with tarfile.open(f"{self.ss_dir}/{model_name}.tar", "r") as tar:
            tar.extractall(path=f"{self.ss_dir}/{model_name}")

        # remove the tar and zipped tar files
        print("Cleaning up")
        os.remove(f"{self.ss_dir}/{model_name}.tar.gz")
        os.remove(f"{self.ss_dir}/{model_name}.tar")

        # add the new model to the list of models
        self.models_list.append(model_name)

        # add snapshots to the list of snapshots
        self.model_snapshots[model_name] = glob.glob(
            f"{self.ss_dir}/{model_name}/initial.*.dat.gz"
        )

        # print done and time taken
        print(f"Done! Took {time.time() - start_time} seconds")
