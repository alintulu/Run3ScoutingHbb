_The code was developed by Adelina Lintuluoto, starting from the [hbb-fitcode](https://github.com/jennetd/hbb-fitcode/tree/master) repository by Jennet Dickinson._

# Directory structure

Please note, this repository is a "work-in-progess". It uses the [Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) framework and [rhalphalib](https://github.com/nsmith-/rhalphalib) python package. The 2 folders contain:

- [hbb](hbb), the code used for the analysis where $H \rightarrow b\bar{b}$ is signal.
- [zbb-hblind](zbb-hblind), the code used for the analysis where $Z \rightarrow b\bar{b}$ is signal and Higgs mass window is blinded.

Each folder contains a README with more information specific to each analysis. However, first, **setup the correct environment** by following the below instructions.
  
# Setup the correct environment

The following instructions were executed on lxplus7. 

1. Setup CMSSW. These commands were taken from the [Combine tutorial](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part5/longexercise/#getting-started) at November 29th 2023.

```bash
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv

bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/main/CombineTools/scripts/sparse-checkout-https.sh)

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
