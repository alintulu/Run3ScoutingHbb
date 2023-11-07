import awkward as ak
import matplotlib.pyplot as plt
import os, sys
import subprocess
import json
import uproot
from coffea.nanoevents import NanoEventsFactory, ScoutingNanoAODSchema
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
    def __init__(self, do_jetid=True):
        self._do_jetid = do_jetid
        
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
            ]
        fatjet = ak.firsts(fatjets)
        fatjet["qcdrho"] = 2 * np.log(fatjet.particleNet_mass / fatjet.pt)
        
        def normalise(val, cut=None):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar
        
        output['h'].fill(
            dataset=dataset,
            n2b1 = normalise(fatjet.n2b1),
            pt = normalise(fatjet.pt),
            rho = normalise(fatjet["qcdrho"]),
            weight=weights.weight(),
        )

        return output
    
    def postprocess(self, accumulator):
        return accumulator
