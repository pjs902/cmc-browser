import glob
import sys

cmctk_path = "/home/peter/research/cmctoolkit"
sys.path.append(cmctk_path)
import cmctoolkit as ck


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
            self.model_snapshots[model] = glob.glob(f"{self.ss_dir}/{model}/initial.*.dat.gz")
        # remove directory from snapshot name
        for model in self.models_list:
            self.model_snapshots[model] = [ss.split("/")[-1] for ss in self.model_snapshots[model]]

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


    def load_snapshot(self, model_name, ss_name="initial.snap0000.dat.gz", distance=5.0):
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
        print(f"Loading Model: {model_name}, snapshot: {ss_name}, Z={Z}, distance={distance}")
        snap = ck.Snapshot(
            fname=f"{self.ss_dir}/{model_name}/{ss_name}", conv=f"{self.ss_dir}/{model_name}/initial.conv.sh", z=Z, dist=distance
        )
        self.loaded_snapshots[f"{model_name}/{ss_name}"] = snap

