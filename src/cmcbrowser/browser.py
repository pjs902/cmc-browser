import glob
import requests
import gzip
import tarfile
import os
import time

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
            self.model_snapshots[model] = glob.glob(
                f"{self.ss_dir}/{model}/initial.*.dat.gz"
            )
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
        self, model_name, ss_name="initial.snap0000.dat.gz", distance=5.0
    ):
        """
        Load a snapshot into the CMCBrowser object.

        Parameters
        ----------
        ss_name : str
            The name of the snapshot to load.
        distance : float
            The distance to use for the snapshot.
        """
        # need to parse metallicity
        Z = float(model_name.split("Z")[-1])
        print(
            f"Loading Model: {model_name}, snapshot: {ss_name}, Z={Z}, distance={distance}"
        )
        snap = ck.Snapshot(
            fname=f"{self.ss_dir}/{model_name}/{ss_name}",
            conv=f"{self.ss_dir}/{model_name}/initial.conv.sh",
            z=Z,
            dist=distance,
        )
        self.loaded_snapshots[f"{model_name}/{ss_name}"] = snap

    def download_new_model(self, N, rv, rg, Z):
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
