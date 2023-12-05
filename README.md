# Boosted Higgs to bb with Run 3 PF scouting

# Setup the correct environment

### LPC

Follow [these intrsuctions](https://github.com/CoffeaTeam/lpcjobqueue) to install LPC DASK executor.

### LXPLUS

```
source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
```

# Analysis flow

### Optimise tagger discriminant score by maximising significance

1. Create a 2D histogram of (a) tagger discriminant score and (b) pT.

```
python submit_dask_[lxplus/lpc]_ddb_score.py
```
2. Find the tagger discriminant score that maximises the significance with [notebooks/ddb_score.ipynb](notebooks/ddb_score.ipynb)

### Run cutflow given the optimised tagger discriminant score

1. Create multi-D histograms containing cutflow

```
python submit_dask_[lxplus/lpc]_cutflow.py
```

2. Use notebooks such as [notebooks/qcd.ipynb](notebooks/.ipynb) and [notebooks/cutflow.ipynb](notebooks/cutflow.ipynb) to visualise the cutflow.

### Create ROOT histograms to use with combine

1. Use notebook [notebooks/roothist.ipynb]

Next, follow the instructions in [combine/README.MD](combine/README.MD) to create the correct environment and how to run Combine.
