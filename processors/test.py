import awkward as ak
import matplotlib.pyplot as plt
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
    n2ddt_shift,
)

class TestProcessor(processor.ProcessorABC):
    def __init__(self, do_jetid=True):
        self._do_jetid = do_jetid
        
    @property
    def accumulator(self):
        return {
            "sumw": defaultdict(float),
            "events": defaultdict(int),
            "cutflow": (
                    Hist.new.Reg(
                        10, 250, 1350, name="pt", label=r"$p_T$ (GeV)"
                    ).IntCategory(
                        [], name="cut", label="Cut Idx", growth=True
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
            ]
        fatjet = ak.firsts(fatjets)
        fatjet["qcdrho"] = 2 * np.log(fatjet.particleNet_mass / fatjet.pt)
        fatjet['n2ddt'] = fatjet.n2b1 - n2ddt_shift(fatjet)
        
        def normalise(val, cut=None):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar

        selection = PackedSelection()
        selection.add('n2ddt', (fatjet.n2ddt < 0.))
        
        regions = {
            'signal': ['n2ddt'],
        } 
        
        for region, cuts in regions.items():
            if region == "noselection":
                continue
            allcuts = set([])
            cut = selection.all(*allcuts)
            weight = weights.weight()[cut]
            
            output['cutflow'].fill(
                pt=normalise(fatjet.pt, cut),
                cut=0,
                weight=weight,
            )
            
            for i, cut in enumerate(cuts):
                allcuts.add(cut)
                cut = selection.all(*allcuts)
                weight = weights.weight()[cut]
                
                output['cutflow'].fill(
                    pt=normalise(fatjet.pt, cut),
                    cut=i+1,
                    weight=weight,
                )

        return output
    
    def postprocess(self, accumulator):
        return accumulator
