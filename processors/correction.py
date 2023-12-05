import numpy as np
import awkward as ak
import os
import uproot
import pickle
import json
from coffea.lumi_tools import LumiData, LumiList
from coffea.nanoevents.methods import vector

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
