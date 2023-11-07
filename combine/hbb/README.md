_The code was developed by Adelina Lintuluoto, starting from the [hbb-fitcode](https://github.com/jennetd/hbb-fitcode/tree/master) repository by Jennet Dickinson._

# Analysis workflow

Running the analysis from scratch requires the following steps:

1. Making of input cards
2. Making of Combine workspace
3. Computing the expected significance
4. Computing the expected signal strength

Everything needed to run each step is stored in this directory, together with their outputs.

**Since this directory is a work-in-progress**, all Combine commands are ran blinded. This is achieved by adding `-t -1` to the arguments.

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

This creates a workspace with 1 parameter of interest (POI): rggF, in `output/testmodel_2022/model_combined.root`

2D plots of the QCD MC pass/fail transfer factor are automatically saved to plots. To draw plots of this fit in each pT bin, do:

```
root -l draw_PFratio_QCDMC.C
```

## Compute the expected significance

```bash
source exp_significance.sh
```

## Compute the expected signal strength

```bash
source exp_shapes.sh
```

The previous step may take a while. When it is finished and a file named `fitDiagnosticsTest.roo` is created, run:

```bash
source exp_mu.sh
```

When it is finished, a file named `higgsCombinerggF.MultiDimFit.mH125.root` is created. 

Finally, create the likelihood plot with:

```bash
root -l draw_likelihood.C
```

An example plot can be found at [plots/likelihood_exp.pdf](plots/likelihood_exp.pdf).
