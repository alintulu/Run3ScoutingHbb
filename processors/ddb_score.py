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
from coffea.lumi_tools import LumiList
from processors.lumiacc import LumiAccumulator

from processors.helper import (
    add_pileup_weight,
    getBosons,
    bosonFlavour,
    pn_disc,
)

from processors.correction import (
    n2ddt_shift,
    correct_met,
)

class DDBScoreProcessor(processor.ProcessorABC):
    def __init__(self, year="2022", jet_arbitration='pt', systematics=False, mass='particleNet_mass'):
        self._year = year
        self._jet_arbitration = jet_arbitration
        self._tightMatch = False
        self._systematics = systematics
        self._mass = mass
        self._ddt_map = "ddt_msd_map.pkl" if "soft" in mass else "ddt_map.pkl"
        
    @property
    def accumulator(self):
        return {
            "lumilist": {},
            "sumw": defaultdict(float),
            "events": defaultdict(int),
            "hist": (
                    Hist.new.Reg(
                        1000, 0.9, 1, name="disc", label=r"$H\rightarrow b\bar{b}$ vs QCD discriminant"
                    ).Var(
                        [40, 110, 139, 201], name="mass", label=r"Mass (GeV)"
                    ).Var(
                        [300, 350, 400, 450, 500, 550, 600, 675, 800, 1200], name="pt", label=r"$p_T$ (GeV)"
#                     ).IntCategory(
#                         [], name="genflav", label="Gen flavour", growth=True
                    ).StrCategory(
                        [], name="dataset", label="Dataset", growth=True
                    ).Weight()
                ),
        }
           
    def process(self, events):
        
        output = self.accumulator
        dataset = events.metadata['dataset']
        
        isRealData = not hasattr(events, "genWeight")
        isQCDMC = 'QCD' in dataset
        
        if not isRealData:
            events = events[
                    (events.Pileup.nPU < 100)
            ]
                
        selection = PackedSelection()
        weights = Weights(len(events), storeIndividual=True)
        
        output['events'][dataset] += len(events)
        
        if isRealData:
            lumi_list = LumiAccumulator(events.run, events.luminosityBlock, auto_unique=True)
            output['lumilist'] = {dataset : lumi_list}
        
        if not isRealData:
            output['sumw'][dataset] += ak.sum(events.genWeight)
            add_pileup_weight(events)
            weights.add('pileup', events['weight_pileup'])
            
        if len(events) == 0:
            return output
        
        fatjets = events.ScoutingFatJet
        fatjets["qcdrho"] = 2 * np.log(fatjets[self._mass] / fatjets.pt)
        fatjets["pn_Hbb"] = pn_disc(
            fatjets.particleNet_prob_Hbb,
            fatjets.particleNet_prob_QCD
        )
        fatjets['n2ddt'] = fatjets.n2b1 - n2ddt_shift(fatjets, self._ddt_map)
        
        jets = events.ScoutingJet
        jets = jets[
            (jets.neHEF < 0.9)
            & (jets.neEmEF < 0.9)
            & (jets.muEmEF < 0.8)
            & (jets.chHEF > 0.01)
            & (jets.nCh > 0)
            & (jets.chEmEF < 0.8)
        ]
        jets["pn_b"] = pn_disc(
            jets.particleNet_prob_b,
            jets.particleNet_prob_g
        )
        
        candidatejet = fatjets[:, :2]
        if self._jet_arbitration == 'pt':
            candidatejet = ak.firsts(candidatejet)
        elif self._jet_arbitration == 'n2':
            candidatejet = ak.firsts(candidatejet[ak.argmin(candidatejet.n2ddt, axis=1, keepdims=True)])
        elif self._jet_arbitration == 'ddb':
            candidatejet = ak.firsts(candidatejet[ak.argmax(candidatejet.pn_Hbb, axis=1, keepdims=True)])
        else:
            raise RuntimeError("Unknown candidate jet arbitration")
                
        if isRealData:
            selection.add("trigger", (events.L1["SingleJet180"] | events.L1["HTT360er"]))
        else:
            selection.add('trigger', np.ones(len(events), dtype='bool'))
            
        selection.add('minjetkin',
            (candidatejet.pt >= 300)
            & (candidatejet.pt < 1200)
#             & (candidatejet.qcdrho < -1.7)
#             & (candidatejet.qcdrho > -6.0)
            & (abs(candidatejet[self._mass]) >= 40)
            & (abs(candidatejet[self._mass]) < 201)
            & (abs(candidatejet.eta) < 2.5)
        )
        
        selection.add('jetid',
            (candidatejet.neHEF < 0.9)
            & (candidatejet.neEmEF < 0.9)
            & (candidatejet.muEmEF < 0.8)
            & (candidatejet.chHEF > 0.01)
            & (candidatejet.nCh > 0)
            & (candidatejet.chEmEF < 0.8)
         )
        
        selection.add('n2ddt', (candidatejet.n2ddt < 0.))
        
        if isRealData:
            era = dataset.split(self._year, 1)[1][0]
            corrected_met = correct_met(events.ScoutingMET, events.ScoutingPrimaryVertex, self._year, era)
            selection.add('met', corrected_met.pt < 140.)
        else:
            selection.add('met', events.ScoutingMET.pt < 140.)
        
        goodmuon = (
            (events.ScoutingMuon.pt > 10)
            & (abs(events.ScoutingMuon.eta) < 2.4)
            & (abs(events.ScoutingMuon.trk_dxy) < 0.2)
            & (abs(events.ScoutingMuon.trackIso) < 0.15)
            & (abs(events.ScoutingMuon.trk_dz) < 0.5)
            & (events.ScoutingMuon.normchi2 < 10)
            & (events.ScoutingMuon.nValidRecoMuonHits > 0)
            & (events.ScoutingMuon.nRecoMuonMatchedStations > 1)
            & (events.ScoutingMuon.nValidPixelHits > 0)
            & (events.ScoutingMuon.nTrackerLayersWithMeasurement > 5)
        )
 
        goodelectron = (
            (events.ScoutingElectron.pt > 10)
            & (abs(events.ScoutingElectron.eta) < 2.4)
            & (events.ScoutingElectron.sigmaIetaIeta < 0.0103)
            & (abs(events.ScoutingElectron.dPhiIn) < 0.127)
            & (abs(events.ScoutingElectron.dEtaIn) < 0.00481)
            & (abs(events.ScoutingElectron.ooEMOop) < 0.0966)
            & (events.ScoutingElectron.missingHits <= 1)
        )
        
        nmuons = ak.sum(goodmuon, axis=1)
        nelectrons = ak.sum(goodelectron, axis=1)
        leadingmuon = ak.firsts(events.ScoutingMuon[goodmuon])
        
        selection.add('noleptons', (nmuons == 0) & (nelectrons == 0))
        
        selection.add('muonDphiAK8', abs(leadingmuon.delta_phi(candidatejet)) > 2*np.pi/3)
        
        if isRealData :
            genflavour = ak.zeros_like(candidatejet.pt)
        else:
            weights.add('genweight', events.genWeight)

            bosons = getBosons(events.GenPart)
            matchedBoson = candidatejet.nearest(bosons, axis=None, threshold=0.8)
            if self._tightMatch:
                match_mask = ((candidatejet.pt - matchedBoson.pt)/matchedBoson.pt < 0.5) & ((candidatejet[self._mass] - matchedBoson.mass)/matchedBoson.mass < 0.3)
                selmatchedBoson = ak.mask(matchedBoson, match_mask)
                genflavour = bosonFlavour(selmatchedBoson)
            else:
                genflavour = bosonFlavour(matchedBoson)
            
        mass_matched = candidatejet[self._mass] * (genflavour > 0) + candidatejet[self._mass] * (genflavour == 0)
            
        regions = {
#             'signal': ['n2ddt'],
            'signal': ['trigger','minjetkin','jetid','met','noleptons'],
        }
        
        def normalise(val, cut):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar

        def normalise(val, cut):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar

        for region, cuts in regions.items():
            if region == "noselection":
                continue
            cut = selection.all(*cuts)
            weight = weights.weight()[cut]

            output['hist'].fill(
                mass=normalise(mass_matched, cut),
                dataset=dataset,
                pt=normalise(candidatejet.pt, cut),
#                 genflav=normalise(genflavour, cut),
                disc=normalise(candidatejet.pn_Hbb, cut),
                weight=weight,
            )

        return output
    
    def postprocess(self, accumulator):
        return accumulator
