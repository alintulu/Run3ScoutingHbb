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
)

class DDTProcessor(processor.ProcessorABC):
    def __init__(self, do_jetid=True, mass='particleNet_mass'):
        self._do_jetid = do_jetid
        self._mass = mass
        
    @property
    def accumulator(self):
        return {
            "sumw": defaultdict(float),
            "events": defaultdict(int),
            "h": (
                    Hist.new.Reg(
                        100, -8, -0.5, name="rho", label=r"$\rho=ln(m^2_{reg}/p_T^2)$"
                    ).Reg(
                        100, 250, 1350, name="pt", label=r"$p_T$ (GeV)"
                    ).Reg(
                        1000, 0, 1, name="n2b1", label=r"$n^2_1$"
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
        fatjets['weight'] = ak.broadcast_arrays(weights.weight(), fatjets.pt)[0]                                                            
                                                                     
        if self._do_jetid:
            fatjets = fatjets[
                (fatjets.neHEF < 0.9)
                & (fatjets.neEmEF < 0.9)
                & (fatjets.muEmEF < 0.8)
                & (fatjets.chHEF > 0.01)
                & (fatjets.nCh > 0)
                & (fatjets.chEmEF < 0.8)
            ]
        fatjets["qcdrho"] = 2 * np.log(fatjets[self._mass] / fatjets.pt)
#         fatjet = ak.firsts(fatjets)
        
        def normalise(val, cut=None):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar
        
#         output['h'].fill(
#             dataset=dataset,
#             n2b1 = normalise(fatjet.n2b1),
#             pt = normalise(fatjet.pt),
#             rho = normalise(fatjet["qcdrho"]),
#             weight=weights.weight(),
#         )

        output['h'].fill(
            dataset=dataset,
            n2b1 = normalise(ak.flatten(fatjets.n2b1)),
            pt = normalise(ak.flatten(fatjets.pt)),
            rho = normalise(ak.flatten(fatjets["qcdrho"])),
            weight = normalise(ak.flatten(fatjets["weight"])),
        )

        return output
    
    def postprocess(self, accumulator):
        return accumulator
