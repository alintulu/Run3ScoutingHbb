import awkward as ak
import os, sys
import subprocess
import json
import uproot
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
from coffea.lookup_tools.lookup_base import lookup_base
import numpy as np
from coffea import processor, util
import coffea.hist as hist
from coffea.analysis_tools import Weights, PackedSelection

from processors.helper import (
    add_pileup_weight,
)

class DDTCoffeaProcessor(processor.ProcessorABC):
    def __init__(self, do_jetid=True, mass='particleNet_mass'):
        self._do_jetid = do_jetid
        self._mass = mass
        self._accumulator = processor.dict_accumulator(
            {
                "sumw": processor.defaultdict_accumulator(float),
                "jets": hist.Hist("Events",
                    hist.Cat("sample", "Sample"),
                    hist.Bin("rho", r"$\rho = \ln (m^2 / p_T^2)$", 100, -8, -0.5),
                    hist.Bin("pt", r"$p_T$ (GeV)", 100, 250, 1350),
                    hist.Bin("n2b1", r"n2b1", 1000, 0, 1),
                ),
            }
        )

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        
        events = events[
                (events.Pileup.nPU < 100)
        ]
        
        output = self.accumulator.identity()
        dataset = events.metadata['dataset']
                
        weights = Weights(len(events), storeIndividual=True)
        weights.add('genweight', events.genWeight)

        fatjets = events.ScoutingFatJet                                      
        fatjets['weight'] = ak.broadcast_arrays(weights.weight(), fatjets.pt)[0]                                                            
                                                                     
        fatjets = fatjets[
            (fatjets.neHEF < 0.9)
            & (fatjets.neEmEF < 0.9)
            & (fatjets.muEmEF < 0.8)
            & (fatjets.chHEF > 0.01)
            & (fatjets.nCh > 0)
            & (fatjets.chEmEF < 0.8)
        ]
            
        fatjets["qcdrho"] = 2 * np.log(fatjets[self._mass] / fatjets.pt)
        
        output["sumw"][dataset] += ak.sum(events.genWeight)
        
        def normalise(val, cut=None):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar

        output['jets'].fill(
            sample = dataset,
            n2b1 = normalise(ak.flatten(fatjets.n2b1)),
            pt = normalise(ak.flatten(fatjets.pt)),
            rho = normalise(ak.flatten(fatjets["qcdrho"])),
            weight = normalise(ak.flatten(fatjets["weight"])),
        )
        
        return output

    def postprocess(self, accumulator):
        return accumulator