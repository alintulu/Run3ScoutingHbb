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
    apply_jec,
)

def update(events, collections):
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
    for name, value in collections.items():
        out = ak.with_field(out, value, name)
    return out

class HistProcessor(processor.ProcessorABC):
    def __init__(self, year="2023", jet_arbitration='pt', systematics=True, mass='particleNet_mass'):
        self._year = year
        self._jet_arbitration = jet_arbitration
        self._tightMatch = False
        self._systematics = systematics
        self._mass = mass
        self._ddt_map = "ddt_msd_map.pkl" if "soft" in mass else "ddt_map.pkl"
        self._pt_debug = "pt"
        
    @property
    def accumulator(self):
        return {
            "lumilist": {},
            "sumw": defaultdict(float),
            "events": defaultdict(int),
            "hist": (
                    Hist.new.Reg(
                        23, 40, 201, name="mass", label=r"Mass"
                    ).Var(
                        [300, 350, 400, 450, 500, 550, 600, 675, 800, 1200], name="pt", label=r"$p_T$ (GeV)"
                    ).Var(
                        [-0.1, 0.8167194, 0.95448214, 0.9707, 0.9782, 0.9859, 0.9864132, 0.9945, 0.9962, 0.997, 0.9984, 0.9988, 0.9991, 0.9994, 1.1], name="disc", label=r"$H\rightarrow b\bar{b}$ vs QCD discriminator",
                    ).IntCategory(
                        [], name="genflav", label="Gen flavour", growth=True
                    ).StrCategory(
                        [], name="dataset", label="Dataset", growth=True
                    ).StrCategory(
                        [], name="systematic", label="Systematic", growth=True
                    ).Weight()
                ),
        }
    
    def process(self, events):
        isRealData = not hasattr(events, "genWeight")
        isQCDMC = 'QCD' in events.metadata['dataset']

        if isRealData or isQCDMC:
            self._pt_debug = "pt"
            # Nominal JEC are already applied in data
            return self.process_shift(events, None)

        if np.sum(ak.num(events.ScoutingFatJet, axis=1)) < 1:
            return self.process_shift(events, None)

        jec_cache = {}

        fatjets = apply_jec(events.ScoutingFatJet, "ScoutingRho", events)
                
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
           
    def process_shift(self, events, shift_name):
        
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
        
        if isRealData:
            lumi_list = LumiAccumulator(events.run, events.luminosityBlock, auto_unique=True)
            output['lumilist'] = {dataset : lumi_list}
        
        if not isRealData:
            if shift_name is None:
                output['sumw'][dataset] += ak.sum(events.genWeight)
            add_pileup_weight(events)
            weights.add('pileup', events['weight_pileup'])
            
        if len(events) == 0:
            return output
               
        fatjets = events.ScoutingFatJet
        fatjets["qcdrho"] = 2 * np.log(fatjets[self._mass] / fatjets[self._pt_debug])
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
            (candidatejet[self._pt_debug] >= 300)
            & (candidatejet[self._pt_debug] < 1200)
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
            genflavour = ak.zeros_like(candidatejet[self._pt_debug])
        else:
            weights.add('genweight', events.genWeight)

            bosons = getBosons(events.GenPart)
            matchedBoson = candidatejet.nearest(bosons, axis=None, threshold=0.8)
            if self._tightMatch:
                match_mask = ((candidatejet[self._pt_debug] - matchedBoson.pt)/matchedBoson.pt < 0.5) & ((candidatejet[self._mass] - matchedBoson.mass)/matchedBoson.mass < 0.3)
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
            
        if shift_name is None:
            systematics = [None] + list(weights.variations)
        else:
            systematics = [shift_name]
            
        def fill(region, systematic, wmod=None):
            selections = regions[region]
            cut = selection.all(*selections)
            sname = 'nominal' if systematic is None else systematic
            if wmod is None:
                if systematic in weights.variations:
                    weight = weights.weight(modifier=systematic)[cut]
                else:
                    weight = weights.weight()[cut]
            else:
                weight = weights.weight()[cut] * wmod[cut]
            
            output['hist'].fill(
                dataset=dataset,
                systematic=sname,
                mass=normalise(mass_matched, cut),
                pt=normalise(candidatejet[self._pt_debug], cut),
                disc=normalise(candidatejet.pn_Hbb, cut),
                genflav=normalise(genflavour, cut),
                weight=weight,
            )
            
        for region in regions:
            for systematic in systematics:
                if isRealData and systematic is not None:
                    continue
                fill(region, systematic)

        return output
    
    def postprocess(self, accumulator):
        return accumulator
