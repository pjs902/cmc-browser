import glob
import requests
import gzip
import tarfile
import os
import time
import h5py
import numpy as np
import pandas as pd
import scipy as sp

import cmctoolkit as ck

__all__ = ["CMCBrowser"]


class CMCBrowser:
    def __init__(self, ss_dir="/home/peter/research/CMC-validation/snapshots", model_prefix="king"):
        """
        Create a CMCBrowser object to load and interact with CMC snapshots.

        Parameters
        ----------
        ss_dir : str
            The directory containing the CMC snapshots.
        model_prefix : str (optional) [default="king"]
            The prefix to use for the model name. This is used to identify the hdf5 snapshots.
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
            # check first if there is a tar.gz snapshot, this indicates we're loading from the
            # public grid where there are no window snapshots and everything is in the tar.gz format
            tar_fns = glob.glob(f"{self.ss_dir}/{model}/initial.*.dat.gz")
            if len(tar_fns) == 0:
                # check instead for the hdf5 snapshot
                snap_fn_window = glob.glob(f"{self.ss_dir}/{model}/{model_prefix}.window.snapshots.h5")

                # check also regular hdf5 snapshots
                snap_fn_reg = glob.glob(f"{self.ss_dir}/{model}/{model_prefix}.snapshots.h5")


                # if no window snapshots, set to None
                if len(snap_fn_window) == 0:
                    snap_fn_window = None
                else:
                    # enumerate the datasets in the hdf5 file
                    with h5py.File(snap_fn_window[0], "r") as f:
                        snap_fn_window = list(f.keys())

                # do the same for the regular snapshots
                if len(snap_fn_reg) == 0:
                    snap_fn_reg = None
                else:
                    with h5py.File(snap_fn_reg[0], "r") as f:
                        snap_fn_reg = list(f.keys())
            else:
                snap_fn_window = None
                snap_fn_reg = tar_fns

            # merge the two lists
            self.model_snapshots[model] = {}
            self.model_snapshots[model]["regular"] = snap_fn_reg
            self.model_snapshots[model]["window"] = snap_fn_window

        # remove directory from snapshot name
        for model in self.models_list:
            if self.model_snapshots[model]["regular"] is None:
                continue
            self.model_snapshots[model]["regular"] = [
                ss.split("/")[-1] for ss in self.model_snapshots[model]["regular"]
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
        print("Regular Snapshots:")
        for ss in self.model_snapshots[model_name]["regular"]:
            print(ss)
        if self.model_snapshots[model_name]["window"] is not None:
            print("Window Snapshots:")
            for ss in self.model_snapshots[model_name]["window"]:
                print(ss)

    def load_snapshot(
        self,
        model_name,
        distance=5.0,
        ss_name="king.window.snapshots.h5",
        mode="h5",
        *,
        h5_key=None,
        strict=True,
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
            The key to use for the hdf5 snapshot. If None, will use the last key in the file. Only
            used if mode is "h5".
        strict : bool
            If True, will raise an error if quantities like the tidal radius or core radius are not
            able to be interpolated. If False, will set these quantities to np.nan.
        """

        # parse the output prefix from the snapshot name
        prefix = ss_name.split(".")[0]

        # need to parse metallicity
        Z = float(model_name.split("Z")[-1])

        print(f"Loading Model: {model_name}, snapshot: {ss_name}, Z={Z}, distance={distance}")

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
            snap.data["luminosity_LSUN"] = snap.data["luminosity[LSUN]"]
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
            if "m[MSUN]" in snap.data.columns:
                snap.data["m_MSUN"] = snap.data["m[MSUN]"]
                snap.data["m0_MSUN"] = snap.data["m0[MSUN]"]
                snap.data["m1_MSUN"] = snap.data["m1[MSUN]"]
                snap.data["luminosity_LSUN"] = snap.data["luminosity[LSUN]"]
            elif "m_MSUN" in snap.data.columns:
                snap.data["m[MSUN]"] = snap.data["m_MSUN"]
                snap.data["m0[MSUN]"] = snap.data["m0_MSUN"]
                snap.data["m1[MSUN]"] = snap.data["m1_MSUN"]
                snap.data["luminosity[LSUN]"] = snap.data["luminosity_LSUN"]
            snap.name = f"{model_name}/{h5_key}"

        # add some useful info
        snap.mass = snap.data["m_MSUN"].sum()
        snap.FeH = np.log10(snap.z / 0.02)

        # calculate the central escape velocity, using the gravitational potential at the center of the cluster and the tidal boundary
        esc_log_file = f"{self.ss_dir}/{model_name}/{prefix}.esc.dat"

        # load the log file with pandas, delim is just a space
        esc = pd.read_csv(esc_log_file, sep="\s+")

        esc["t[Myr]"] = esc["#2:t"] * snap.unitdict["myr"]

        # conversions to physical units, following https://github.com/tomas-cabrera/hvss-bsco/blob/main/src/scripts/cmc_single_clusters_vesc.py
        nb_kms = 1e-5 * snap.unitdict["cm"] / snap.unitdict["nb_s"]
        esc["#10:phi_rtidal"] *= nb_kms**2
        esc["#11:phi_zero"] *= nb_kms**2

        esc["vesc"] = np.sqrt(2 * (esc["#10:phi_rtidal"] - esc["#11:phi_zero"]))

        snap.vesc_initial = esc[esc["t[Myr]"] < 20]["vesc"].mean()
        snap.vesc_final = esc["vesc"][-5000:].mean()

        # while we have the escape logs open, get the tidal radius at the current time
        # need to interpolate to get the tidal radius at the current time

        try:
            # fill outside values with the first and last value
            rtidal_interp = sp.interpolate.interp1d(
                esc["t[Myr]"],
                esc["#9:Rtidal"],
                kind="linear",
                bounds_error=False,
                fill_value=(esc["#9:Rtidal"].iloc[0], esc["#9:Rtidal"].iloc[-1]),
            )
            snap.rtidal = float(rtidal_interp(snap.age * 1000)) * snap.unitdict["pc"]
        except ValueError as e:
            if strict:
                msg = f"Unable to interpolate tidal radius for model {model_name} at time {snap.age}, error: {e}"
                raise ValueError(msg) from None
            else:
                print(f"Unable to interpolate tidal radius for model {model_name} at time {snap.age}, error: {e}")
                snap.rtidal = np.nan

        # convert to pc
        snap.rtidal *= snap.unitdict["pc"]

        # initial cluster mass
        dyn_log_file = f"{self.ss_dir}/{model_name}/{prefix}.dyn.dat"

        # load the log file with pandas, delim is just a space
        dyn = pd.read_csv(dyn_log_file, skiprows=1, sep="\s+")

        # get the initial mass
        snap.initial_mass = dyn["#5:M"][0] * snap.unitdict["msun"]

        # while we have the dynamics logs open, get the core radius at the current time
        # need to interpolate to get the core radius at the current time

        dyn["t[Myr]"] = dyn["#1:t"] * snap.unitdict["myr"]

        try:
            rcore_interp = sp.interpolate.interp1d(
                dyn["t[Myr]"],
                dyn["#8:r_c"],
                kind="linear",
                bounds_error=False,
                fill_value=(dyn["#8:r_c"].iloc[0], dyn["#8:r_c"].iloc[-1]),
            )
            snap.rcore = float(rcore_interp(snap.age * 1000)) * snap.unitdict["pc"]
        except ValueError:
            if strict:
                msg = f"Unable to interpolate core radius for model {model_name} at time {snap.age}"
                raise ValueError(msg) from None
            else:
                snap.rcore = np.nan

        # convert to pc
        snap.rcore *= snap.unitdict["pc"]

        # save the mass, half mass radius, tidal radius, and core radius over time to the snapshot object

        snap.evolutionary_quantities = {}
        snap.evolutionary_quantities["time_Gyr"] = dyn["t[Myr]"].to_numpy() / 1000
        snap.evolutionary_quantities["cluster_mass_MSUN"] = dyn["#5:M"].to_numpy() * snap.unitdict["msun"]
        snap.evolutionary_quantities["rh_pc"] = dyn["#21:r_h"].to_numpy() * snap.unitdict["pc"]
        snap.evolutionary_quantities["rcore_pc"] = dyn["#8:r_c"].to_numpy() * snap.unitdict["pc"]

        # grab the central density too
        snap.evolutionary_quantities["rho0_MSUN_pc3"] = (
            dyn["#22:rho_0"].to_numpy() * snap.unitdict["msun"] / snap.unitdict["pc"] ** 3
        )

        # also just save the current central density, interpolated from the log file
        rho0_interp = sp.interpolate.interp1d(
            dyn["t[Myr]"],
            dyn["#22:rho_0"],
            kind="linear",
            bounds_error=False,
            fill_value=(dyn["#22:rho_0"].iloc[0], dyn["#22:rho_0"].iloc[-1]),
        )
        rho0_MSUN_pc3 = float(rho0_interp(snap.age * 1000))
        snap.rho0_MSUN_pc3 = rho0_MSUN_pc3 * snap.unitdict["msun"] / snap.unitdict["pc"] ** 3

        # do the same thing with the BH logs, we want to know the number of BHs over time

        bh_log_file = f"{self.ss_dir}/{model_name}/{prefix}.bh.dat"

        # load the log file with pandas, delim is just a space
        bh = pd.read_csv(bh_log_file, sep="\s+")

        # add in the time column
        bh["t[Myr]"] = bh["#2:TotalTime"] * snap.unitdict["myr"]

        # save the number of BHs over time to the snapshot object
        snap.evolutionary_quantities["n_bh"] = bh["#3:Nbh,tot"].to_numpy()

        # save the bh timestamps too
        snap.evolutionary_quantities["bh_time_Gyr"] = bh["t[Myr]"].to_numpy() / 1000

        # delete the pandas dataframes to free up memory
        del esc
        del dyn
        del bh

        # half mass radius
        r = snap.data["r"]
        M_cdf = np.cumsum(snap.data["m[MSUN]"]) / np.sum(snap.data["m[MSUN]"])
        try:
            # try with scipy interp1d first
            rh = sp.interpolate.interp1d(y=r, x=M_cdf, kind="linear")(0.5)
            rh = float(rh)
        except ValueError:
            # if that fails, try with numpy interp
            try:
                rh = np.interp(0.5, M_cdf, r)
            except ValueError as e:
                # if that fails, raise the error
                if strict:
                    msg = (
                        f"Unable to interpolate half mass radius for model {model_name} at time {snap.age}, error: {e}"
                    )
                    raise ValueError(msg) from None
                else:
                    print(
                        f"Unable to interpolate half mass radius for model {model_name} at time {snap.age}, error: {e}"
                    )
                    rh = np.nan

        # convert to pc
        rh *= snap.unitdict["pc"]
        snap.rh = float(rh)

        # calculate the half mass relaxation time, here using eq 7.108 from Binney and Tremaine
        N = len(snap.data["m[MSUN]"])
        G = 4.3e-3  # pc^3 Msun^-1 Myr^-2
        trh = (0.17 * N) / (np.log(0.1 * N)) * np.sqrt(rh**3 / (G * snap.mass))
        # convert to Gyr
        snap.Trh = trh / 1000

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
        base_url = "https://cmc.ciera.northwestern.edu/download-cluster/download-cluster/"

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
        self.model_snapshots[model_name] = glob.glob(f"{self.ss_dir}/{model_name}/initial.*.dat.gz")

        # print done and time taken
        print(f"Done! Took {time.time() - start_time} seconds")
