## cmc-browser

This small library has two main purposes, managing a local grid of CMC models and handling the
loading of CMC snapshots.

On the managing side of things, this library will index all models and snapshots in a given
directory and also allows for the automated downloading of new CMC models from the web server.

Here's an example of a simple workflow, indexing some local models, listing the snapshots of a model
and downloading a new model:

```python
from cmcbrowser import CMCBrowser

# initialize cmcbrowser in a specific directory
b = CMCBrowser(ss_dir="/home/peter/research/CMC-grid/grid/")


b.list_models()
```

```shell
Found the following models:
N4e5_rv1_rg20_Z0.0002
N4e5_rv1_rg8_Z0.02
N1.6e6_rv0.5_rg8_Z0.002
N2e5_rv0.5_rg2_Z0.02
N1.6e6_rv2_rg20_Z0.02
N1.6e6_rv0.5_rg20_Z0.002
N1.6e6_rv1_rg2_Z0.0002
N4e5_rv1_rg8_Z0.0002
N2e5_rv0.5_rg20_Z0.02
N4e5_rv1_rg8_Z0.002
N4e5_rv2_rg20_Z0.02
N1.6e6_rv2_rg20_Z0.002
```

```python
b.list_snapshots(model_name="N2e5_rv0.5_rg2_Z0.02")
```

```shell
Found the following snapshots:
initial.snap0489.dat.gz
initial.snap0036.dat.gz
initial.snap0000.dat.gz
```

On the loading side of things, this library augments the Snapshot loading functionality provided by
[cmctoolkit](https://github.com/NicholasRui/cmctoolkit).

Most of this functionality is convenience, like parsing the metallicity from the model name or
collecting useful quantities like the initial cluster mass, half-mass radius, escape velocity and
core radius. It also allows for easy loading of snapshots in either the old tar.gz format or the new
HDF5 format.

```python
# load in a snapshot, place it at 5 kpc
b.load_snapshot(
    model_name="N2e5_rv0.5_rg20_Z0.02",
    h5_key="3(t=9.5000169Gyr)",
    distance=5.0,
    mode="h5",
)
```

This returns a regular `cmctoolkit` `Snapshot` object, just with some of that extra info attached
which means that we can directly use the methods implemented in `cmctoolkit` like adding photometry
or make a 2D projection:

```python
# select a loaded snapshot
snap = b.loaded_snapshots["N2e5_rv0.5_rg20_Z0.02/3(t=9.5000169Gyr)"]

# load in the list of filters we want to calculate photometry for
filttable = ck.load_filtertable("/home/peter/research/cmctoolkit/filt_index.txt")

# add photometry
snap.add_photometry(filttable)

# add in projected 2d postions and velocities
snap.make_2d_projection()
```

Finally, downloading new models:

```python
b.download_new_model(N="4e5", rv="1", rg="20", Z="0.0002")
```

Which will construct the appropriate URL, send the request to the CMC web server, download the
model, unpack it, add it the existing grid, index the snapshots and clean up any leftover files.


### To Do
- Fix/clean up the handling of window vs regular snapshots


## See also

[`cmc-obs`](https://github.com/pjs902/cmc-obs): Another small library which takes in a `Snapshot` object from `cmc-browser` and extracts a series of mock observations, designed to be as realistic as possible. 
