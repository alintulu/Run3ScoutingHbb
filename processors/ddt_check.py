import awkward as ak
import os, sys
import subprocess
import json
import uproot
from coffea.lookup_tools.lookup_base import lookup_base
import numpy as np
from coffea import processor, util
from hist import Hist
import hist
from coffea.analysis_tools import Weights, PackedSelection
from collections import defaultdict

from processors.helper import (
    add_pileup_weight,
    pn_disc,
)

from processors.correction import (
    n2ddt_shift,
)

class DDTCheckProcessor(processor.ProcessorABC):
    def __init__(self, do_jetid=True, mass='particleNet_mass'):
        self._do_jetid = do_jetid
        self._mass = mass
        self._ddt_map = "ddt_msd_map.pkl" if "soft" in mass else "ddt_map.pkl"
        
    @property
    def accumulator(self):
        return {
            "sumw": defaultdict(float),
            "events": defaultdict(int),
            "h": (
                    Hist.new.Reg(
                        100, 20, 200, name="mass", label=r"Mass (GeV)"
                    #).Reg(
                    #    100, 300, 1200, name="pt", label=r"$p_T$ (GeV)"
                    ).Reg(
                        100, 0, 1, name="n2ddt", label=r"$n^2_1$"
                    ).Reg(
                        100, 0.8, 1, name="disc", label=r"$H\rightarrow b\bar{b}$ vs QCD discriminator"
                    ).StrCategory(
                         [], name="dataset", label="Dataset", growth=True
                    ).Weight()
                ),
        }
           
        
    def process(self, events):
        
        events = events[
                (events.Pileup.nPU < 100)
        ]
        
        output = self.accumulator
        dataset = events.metadata['dataset']
                
        weights = Weights(len(events), storeIndividual=True)
        
        output['events'][dataset] += len(events)
        output['sumw'][dataset] += ak.sum(events.genWeight)
        weights.add('genweight', events.genWeight)
        
        add_pileup_weight(events)
        weights.add('pileup', events['weight_pileup'])
            
        if len(events) == 0:
            return output
                
        fatjets = events.ScoutingFatJet                                      
                                                                     
        if self._do_jetid:
            fatjets = fatjets[
                (fatjets.neHEF < 0.9)
                & (fatjets.neEmEF < 0.9)
                & (fatjets.muEmEF < 0.8)
                & (fatjets.chHEF > 0.01)
                & (fatjets.nCh > 0)
                & (fatjets.chEmEF < 0.8)
                & (fatjets.pt > 600)
            ]
        fatjets["qcdrho"] = 2 * np.log(fatjets[self._mass] / fatjets.pt)
        fatjets['n2ddt'] = fatjets.n2b1 - n2ddt_shift(fatjets, self._ddt_map)
        fatjets["pn_Hbb"] = pn_disc(
            fatjets.particleNet_prob_Hbb,
            fatjets.particleNet_prob_QCD
        )
        fatjet = ak.firsts(fatjets)
        
        def normalise(val, cut=None):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar
        
        output['h'].fill(
            dataset = dataset,
            n2ddt = normalise(fatjet.n2ddt),
            #pt = normalise(fatjet.pt),
            mass = normalise(fatjet[self._mass]),
            disc = normalise(fatjet.pn_Hbb),
            weight = weights.weight(),
        )

        return output
    
    def postprocess(self, accumulator):
        return accumulator
