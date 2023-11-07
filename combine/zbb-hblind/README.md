_The code was developed by Adelina Lintuluoto, starting from the [hbb-fitcode](https://github.com/jennetd/hbb-fitcode/tree/master) repository by Jennet Dickinson._

# Analysis workflow

Running the analysis from scratch requires the following steps:

1. Making of input cards
2. Making of Combine workspace
3. Computing the expected and observed significance
4. Computing the observed signal strength

Everything needed to run each step is stored in this directory, together with their outputs.

## Making the input cards

Make sure that the directory contains a file called `signalregion.root`, which is the ROOT histogram created from the step detailed [here](https://github.com/alintulu/Run3ScoutingHbb/tree/lxplus?tab=readme-ov-file#converting-the-hist-histogram-to-root-histogram).

```bash
python3 make_cards.py
```

This creates a directory called `output`, which contains the Combine input cards.

## Make the Combine workspace

```bash
source make_workspace.sh
```

This creates a workspace with 1 parameter of interest (POI): rZbb, in `output/testmodel_2022/model_combined.root`

2D plots of the QCD MC pass/fail transfer factor are automatically saved to plots. To draw plots of this fit in each pT bin, do:

```
root -l draw_PFratio_QCDMC.C
```

## Compute the expected and observed significance

```bash
source exp_significance.sh
source obs_significance.sh
```

## Compute the observed signal strength

```bash
source obs_shapes.sh
```

The previous step may take a while. When it is finished and a file named `fitDiagnosticsTest.root` is created, run:

```bash
combine output/testModel_2022/model_combined.root -m 125 -M MultiDimFit --saveWorkspace -n zbb.postfit

combine higgsCombinehzbb.postfit.MultiDimFit.mH125.root -M MultiDimFit -n zbb.total --algo grid --snapshotName MultiDimFit --setParameterRanges r=0,4

combine higgsCombinehzbb.postfit.MultiDimFit.mH125.root -M MultiDimFit --algo grid --snapshotName MultiDimFit --setParameterRanges r=0,4  --freezeParameters allConstrainedNuisances -n zbb.freeze_all
```

When it is finished, files named `higgsCombinezbb.total.MultiDimFit.mH120.root ` and `higgsCombinezbb.freeze_all.MultiDimFit.mH120.root` are created. 

Finally, create the likelihood plot with:

```bash
./../../../../CombineHarvester/CombineTools/scripts/plot1DScan.py higgsCombinezbb.total.MultiDimFit.mH120.root --main-label "Total uncertainty" --others higgsCombinezbb.freeze_all.MultiDimFit.mH120.root:"Statistical uncertanity":6 --output breakdown --y-max 10 --y-cut 40 --breakdown "sys.,stat." --POI rZbb --logo-sub Private
```

An example plot can be found at [plots/likelihood_obs.pdf](plots/likelihood_obs.pdf).
