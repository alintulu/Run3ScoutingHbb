# Table of contents

- [Setup the correct environment](#setup-the-correct-environment)
- [Run the code](#run-the-code)

# Setup the correct environment

The following instructions were executed on lxplus7. 

1. Setup CMSSW. These commands were taken from the [Combine tutorial](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part5/longexercise/#getting-started) at November 29th 2023.

```bash
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit

cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.0.0
scramv1 b clean; scramv1 b
```

2. Install rhalphalib

```bash
cd $CMSSW_BASE/src/
python3 -m pip install --user https://github.com/nsmith-/rhalphalib/archive/master.zip
```

# Run the code

Start by making the workspace:

```bash
python3 make_cards.py
source make_workspace.sh
```

Creates a workspace with 1 POI: rggF, in `output/testmodel_2022/model_combined.root`

2D plots of the QCD MC pass/fail transfer factor are automatically saved to plots. To draw plots of this fit in each pT bin, do

```
root -l draw_PFratio_QCDMC.C
```
