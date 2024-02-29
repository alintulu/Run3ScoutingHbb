_The code was developed by Adelina Lintuluoto, starting from the [hbb-coffea](https://github.com/jennetd/hbb-coffea/tree/master) repository by Jennet Dickinson._

This README documents:

1. [Directory structure](#directory-structure), the contents of this Github repository
2. [Analysis workflow](#analysis-workflow), details on how to run the analysis from scratch
3. [Open items](#open-items), list of items that should be addressed

# Directory structure

Please note, this repository is a "work-in-progess". It contains code developed with Coffea 0.7.21 for the analysis of Higgs bosons produced at high transverse momentum and their subsequent decay into bottom quark-antiquark pairs using Run 3 scouting data. The 4 folders contain:

- [combine](combine), the [Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) code used for the statistical analysis
- [data](data), the input data
- [notebooks](notebooks), the notebooks used for plotting
- [processors](processors), the Coffea 0.7 processors used to create histograms 

# Analysis workflow

Running the analysis from scratch requires the following steps:

1. [Event selection](#event-selection)
   - Optimising the event selection by stuyding the cutflow
   - Optimising the double b-tagging discriminant score
2. [Creating the final histogram](#creating-the-final-histogram)
   - Creating a hist histogram
   - Converting the hist histogram to ROOT histogram
3. [Statistical analysis](combine)

Each step will now be discussed in detail, **but first, setup the correct environment:**

1. Log into lxplus
2. Make sure you have a valid grid proxy
3. Source an environment containing Coffea 0.7.21

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh
```

## Event selection

### Optimising the event selection by stuyding the cutflow

In order to optimse the event selection, a multi-dimensional histogram is created with [processors/cutflow.py](processors/cutflow.py). 

Selection steps are added to the [PackedSelection](https://coffeateam.github.io/coffea/api/coffea.analysis_tools.PackedSelection.html#packedselection) object, for example as follows:

```python
selection.add('minjetkin',
      (candidatejet.pt >= 300)
      & (candidatejet.pt < 1200)
      & (abs(candidatejet[self._mass]) >= 40)
      & (abs(candidatejet[self._mass]) < 201)
      & (abs(candidatejet.eta) < 2.5)
  )
```

This example can be found [here](https://github.com/alintulu/Run3ScoutingHbb/blob/lxplus/processors/cutflow.py#L128-L136) in the code, and shows how the minimal kinematic jet selection is added.

In order to actually apply the selection before filling the histogram, you are required to add it to the following [list](https://github.com/alintulu/Run3ScoutingHbb/blob/lxplus/processors/cutflow.py#L205).

```python
regions = {
      'signal': ['trigger','minjetkin','jetid','met','noleptons'],
    }
```

The order of the selections in this list dictates the order in which they are applied.

Next, each selection step is consecutively applied to the events, filling the histogram before each new selection. 

The histogram variable `cut` corresponds to the index of each selection step in the above list: `0` corresponds to no selection, `1` corresponds to the 1st selection applied, `2` corresponds to the 2nd selection applied on top of the 1st one, etc.

#### Submission

To create the histogram using DASK, submit the jobs with  [submit_dask_lxplus_cutflow.py](submit_dask_lxplus_cutflow.py).

Before submission, create a grid proxy located in your AFS user space (in the command below, **exchange the AFS path to your path**):

```bash
voms-proxy-init --rfc --voms cms:/cms/dcms -valid 192:00 --out /afs/cern.ch/user/a/adlintul/private/gridproxy.pem
```

The path to the grid proxy is then inserted into the submission file [like this](https://github.com/alintulu/Run3ScoutingHbb/blob/main/submit_dask_lxplus_cutflow.py#L18):

```python
proxy_path = "/afs/cern.ch/user/a/adlintul/private/gridproxy.pem"
os.environ['X509_USER_PROXY'] = proxy_path
if os.path.isfile(os.environ['X509_USER_PROXY']):
    print("Found proxy at {}".format(os.environ['X509_USER_PROXY']))
else:
    print("os.environ['X509_USER_PROXY'] ",os.environ['X509_USER_PROXY'])
os.environ['X509_CERT_DIR'] = '/cvmfs/cms.cern.ch/grid/etc/grid-security/certificates'
os.environ['X509_VOMS_DIR'] = '/cvmfs/cms.cern.ch/grid/etc/grid-security/vomsdir'
os.environ['X509_USER_CERT'] = proxy_path

env_extra = [
        'export XRD_RUNFORKHANDLER=1',
        'export X509_USER_PROXY={}'.format(proxy_path),
        'export X509_CERT_DIR={}'.format(os.environ["X509_CERT_DIR"]),            
    ]
```

Finally, submit with:

```bash
python submit_dask_lxplus_cutflow.py
```

#### Reweighting the MC histograms

To reweight the MC histogram according to luminosity and cross section, use the following notebook: [notebooks/group.ipynb](notebooks/group.ipynb).

#### Plotting

To create tables and plots to evaluate the cutflow, use the following notebook: [notebooks/cutflow.ipynb](notebooks/cutflow.ipynb).

### Optimising the double b-tagging discriminant score

To facilitate background estimation, signal candidate jets are divided into a signal and control region. These regions are identified by being either above or below a double b-tagging discriminant value, and are chosen to optimise sensitivity to the signal by maximising the number-counting significance.

A multi-dimensional histogram is created with [processors/ddb_score.py](processors/ddb_score.py). In the script, the same selections as detailed in [processors/cutflow.py](processors/cutflow.py) are added.

#### Submission

To create the histogram using DASK, submit the jobs with  [submit_dask_lxplus_ddb_score.py](submit_dask_lxplus_ddb_score.py).

Make sure you have created the grid proxy as detailed in [Optimising the event selection by stuyding the cutflow - Submission](#optimising-the-event-selection-by-stuyding-the-cutflow), then submit with:

```
python submit_dask_lxplus_ddb_score.py
```

#### Reweighting the MC histograms

To reweight the MC histogram according to luminosity and cross section, use the following notebook: [notebooks/group.ipynb](notebooks/group.ipynb).

#### Plotting

To find the optimised value, use the following notebook: [notebooks/ddb_score.ipynb](notebooks/ddb_score.ipynb).

The value is found by maximising the significance:

$Z = \sqrt{2 \cdot (N_{S} +  N_{B}) \cdot \text{ln}(1 + N_{S} /  N_{B}) - 2  N_{S}}$,

where $N_{S}$ and $N_{B}$ is the number of signal and background events, respectively.

## Creating the final histogram

### Creating a hist histogram

Now that the event selection is optimised, and the signal and control regions are selected, it's time to fill the final histogram.

A multi-dimensional histogram is created with [processors/hist.py](processors/hist.py). In the script, the same selections as detailed in [processors/cutflow.py](processors/cutflow.py) are added.

Up- and down-varied histograms of the JES and JER systematics are created as [follows](https://github.com/alintulu/Run3ScoutingHbb/blob/master/processors/hist.py#L85-L95):

```python
shifts = [({"FatJet": fatjets}, None)]
if self._systematics:
   shifts = [
       ({"ScoutingFatJet": fatjets}, None),
       ({"ScoutingFatJet": fatjets.JES_jes.up}, "JESUp"),
       ({"ScoutingFatJet": fatjets.JES_jes.down}, "JESDown"),
       ({"ScoutingFatJet": fatjets.JER.up}, "JERUp"),
       ({"ScoutingFatJet": fatjets.JER.down}, "JERDown"),
   ]
       
return processor.accumulate(self.process_shift(update(events, collections), name) for collections, name in shifts)
```

The process that fills the final histogram is repeated several times, each time using a nominal, up- or down-varied verison of the large-radius jets.

#### Submission

To create the histogram using DASK, submit the jobs with  [submit_dask_lxplus_hist.py](submit_dask_lxplus_hist.py).

Make sure you have created the grid proxy as detailed in [Optimising the event selection by stuyding the cutflow - Submission](#optimising-the-event-selection-by-stuyding-the-cutflow), then submi with:

```
python submit_dask_lxplus_hist.py
```

#### Reweighting the MC histograms

To reweight the MC histogram according to luminosity and cross section, use the following notebook: [notebooks/group.ipynb](notebooks/group.ipynb).

### Converting the hist histogram to ROOT histogram

In order to input the histogram to Combine, the hist histogram needs to be converted to a ROOT histogram. To do so, use the following notebook: [notebooks/roothist.ipynb](notebooks/roothist.ipynb).

# Open items

#### The following items should be addressed:

1. Run 2 JES and JER SF and SF associated uncertanities are applied [here](https://github.com/alintulu/Run3ScoutingHbb/blob/main/processors/correction.py#L83-L88). These were added to provide an estimate of the JES and JER systematics, however, they should be exchanged for Run 3 (preferably scouting specific) SFs.
2. MET corrections for data are provided by Connor Houghton, and are applied [here](https://github.com/alintulu/Run3ScoutingHbb/blob/main/processors/hist.py#L184-L189). Corrections for MC should also be applied (my understanding is that Connor is working on this). In additon, only 2022 corrections exist at the moment so they are applied also for data recorded in 2023. When 2023 corrections exist, makes sure to add them.
3. Inclusion of more systematic uncertanities. For the implementation, Jennet's file [here](https://github.com/jennetd/hbb-fitcode/blob/master/main-vbf/prefit/make_cards.py#L254) is a good example to follow.

#### The following items could be addressed (but are not necessary to my knowledge):

1. If an AK4 jet flavour tagger exist, one may add the following selection: veto any event with a b-tagged AK4 jet in the opposite hemisphere of the signal candidate jet. This rejects TTbar contamination. An idea of how it might be implemented can be found [here](https://github.com/jennetd/hbb-coffea/blob/master/boostedhiggs/vbfprocessor.py#L278-L280).
2. If scouting taus exist, reject events with taus in the same way as events with muons or electrons are rejected [here](https://github.com/alintulu/Run3ScoutingHbb/blob/main/processors/hist.py#L218).


