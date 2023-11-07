import numpy as np
import awkward as ak
import os
import uproot
import pickle
import json
import gzip
import cloudpickle
from coffea.lumi_tools import LumiData, LumiList
from coffea.nanoevents.methods import vector
from coffea.lookup_tools import extractor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty

def n2ddt_shift(fatjets, ddt_map):
    path_ddt_map = os.path.join(
        os.path.dirname(__file__),
        f"../data/n2b1/{ddt_map}",
    )
    ddt_map = pickle.load(open(path_ddt_map,'rb'))

    return ddt_map(fatjets.qcdrho, fatjets.pt)

def correct_met(met, pv, year, era):
    path_met_corr = os.path.join(
        os.path.dirname(__file__),
        "../data/met/PFScouting_METXY_corrections.json",
    )
    with open(path_met_corr, "r") as json_file:
        met_corr_dict = json.load(json_file)
        
    def linear_function(slope, intercept, npm):
        return(slope * npm + intercept)
    
    correction_met = ak.zip(
        {
            "x": linear_function(
                met_corr_dict[year][era]['x']['m'],
                met_corr_dict[year][era]['x']['c'],
                ak.num(pv),
            ),
            "y": linear_function(
                met_corr_dict[year][era]['y']['m'],
                met_corr_dict[year][era]['y']['c'],
                ak.num(pv),
            ),
            "z": np.ones(len(met))
        },
        with_name="ThreeVector",
        behavior=vector.behavior,
    )
    
    raw_met = ak.zip(
        {
            "pt": met.pt,
            "eta": np.ones(len(met)),
            "phi": met.phi,
            "mass": np.ones(len(met)),
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )
    
    new_met = (raw_met - correction_met)
    
    return ak.zip(
        {
            "pt": new_met.pt,
            "eta": np.ones(len(met)),
            "phi": new_met.phi,
            "mass": np.ones(len(met)),
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )

path_jec = os.path.join(
    os.path.dirname(__file__),
    "../data/jec",
)

ext = extractor()
ext.add_weight_sets([
    f"* * {path_jec}/jennet/Summer19UL18_V5_MC_L2Residual_AK8PFPuppi.txt",
    f"* * {path_jec}/jennet/Summer19UL18_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
    f"* * {path_jec}/jennet/Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.jr.txt",
    f"* * {path_jec}/jennet/Summer19UL18_JRV2_MC_SF_AK8PFPuppi.jersf.txt",
])
ext.finalize()
evaluator = ext.make_evaluator()

jec_stack_names = evaluator.keys()
jec_inputs = {name: evaluator[name] for name in jec_stack_names}
jec_stack = JECStack(jec_inputs)

name_map = jec_stack.blank_name_map
name_map['JetPt'] = 'pt'
name_map['JetMass'] = 'mass'
name_map['JetEta'] = 'eta'
name_map['JetA'] = 'area'
name_map['ptRaw'] = 'pt_raw'
name_map['massRaw'] = 'mass_raw'
name_map['Rho'] = 'event_rho'
name_map['ptGenJet'] = 'pt_gen'

jet_factory = CorrectedJetsFactory(name_map, jec_stack)

#with gzip.open(path_jec + '/jennet/jec_compiled.pkl.gz') as fin:
#    jmestuff = cloudpickle.load(fin)
#
#jet_factory = jmestuff["fatjet_factory"]["2018mc"]
    
def apply_jec(jets, rho_name, events):

    jets["pt_raw"] = jets["pt"]
    jets["mass_raw"] = jets["mass"]
    jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["event_rho"] = ak.broadcast_arrays(events[rho_name], jets.pt)[0]
    corrected_jets = jet_factory.build(jets, lazy_cache=events.caches[0])

    return corrected_jets
