#!/usr/bin/env python
import os, sys
import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from os import system, environ
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
import numpy as np
#import "$CMSSW_BASE/src/AnalysisTools/QuickRefold/interface/TObjectContainer.h"

class qcdjetptSFProducer(Module): 
    def __init__(self, era):
	self.writeHistFile=True
	self.metBranchName="MET"
	self.debug = True
	# if era == "2016":
        # 	self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2016/qcdJetptCorr.root"
	# elif era == "2017":
        # 	self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2017/qcdJetptCorr.root"
	# elif era == "2018":
        # 	self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2018/qcdJetptCorr.root"
        # self.correctionhistLM="ptnorm"
        # self.correctionhistHM="ptnorm1"
        # self.correctionhistLM1="ptnorm2"
        # self.correctionhistHM1="ptnorm3"
        if era == "2016":
            self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2016/qcdJetRespTailCorr.root"
        elif era == "2017":
            self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2017/qcdJetRespTailCorr.root"
        elif era == "2018":
            self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2018/qcdJetRespTailCorr.root"
        self.correctionhist="RespTailCorr"
        self.correctionhistHM="RespTailCorrHM"

    def beginJob(self,histFile=None,histDirName=None):
   	pass

    def endJob(self):
	pass 

    def loadhisto(self,filename,hname):
        file =ROOT.TFile.Open(filename)
        hist_ = file.Get(hname)
        hist_.SetDirectory(0)
        return hist_

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        # self.out.branch("qcdJetPtWeightLMwq","F")
        # self.out.branch("qcdJetPtWeightLM","F")
        # self.out.branch("qcdJetPtWeightHMwq","F")
        # self.out.branch("qcdJetPtWeightHM","F")

        # self.corrhistLM1=self.loadhisto(self.correctionfile,self.correctionhistLM1)
        # self.corrhistHM1=self.loadhisto(self.correctionfile,self.correctionhistHM1)
        # self.corrhistLM=self.loadhisto(self.correctionfile,self.correctionhistLM)
        # self.corrhistHM=self.loadhisto(self.correctionfile,self.correctionhistHM)
        self.out.branch("qcdRespTailWeight","F")
        self.out.branch("qcdRespTailWeightHM","F")
        self.corrhist=self.loadhisto(self.correctionfile,self.correctionhist)
        self.corrhistHM=self.loadhisto(self.correctionfile,self.correctionhistHM)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def calculateNLeptons(self, eles, muons, isks, taus):
        countEle = sum([e.Stop0l for e in eles])
        countMu  = sum([m.Stop0l for m in muons])
        countIsk = sum([i.Stop0l for i in isks])
        countTau = sum([t.Stop0l for t in taus])
        return countEle, countMu, countIsk, countTau

    def PassEventFilter(self, flags):
       # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2016_data
        passEventFilter = None
       
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2016 ~~~~~
        # if self.era == "2016":
        #     ## Common filters
        #     passEventFilter = flags.goodVertices and flags.HBHENoiseFilter and \
        #         flags.HBHENoiseIsoFilter and flags.EcalDeadCellTriggerPrimitiveFilter \
        #         and flags.BadPFMuonFilter 
                      
                
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2017 ~~~~~
        #if self.era == "2017" or self.era == "2018":
            ## Common filters
            ## Missing the latest ecalBadCalibReducedMINIAODFilter, not in MiniAOD
            ## But still using the old ecalBadCalibFilter from MiniAOD
        passEventFilter = flags.goodVertices and flags.HBHENoiseFilter and flags.HBHENoiseIsoFilter and flags.EcalDeadCellTriggerPrimitiveFilter and flags.BadPFMuonFilter and flags.ecalBadCalibFilter  
            ## Only data
        
        return passEventFilter
    def getQCDRespTailCorrector(self, jets, genJets, met):
        MM = -1
        ind = -1
        resp = -1
       
        for iG in range(0,len(genJets)):
                if iG == 4: break
                if genJets[iG].pt == 0: break
                fpt = -1
                for rJ in xrange(len(jets)):
                        if not jets[rJ].Stop0l: continue
                        if jets[rJ].genJetIdx != iG:continue
                        fpt = jets[rJ].pt
                        break
                if fpt < 0: fpt = 9.5 #for the cases with no reco jet due to being below thresh                                                                                        
                if abs(fpt - genJets[iG].pt) > MM:
                        ind = iG
                        resp =float(fpt)/float(genJets[iG].pt)
                        MM = abs(float(fpt - genJets[iG].pt))
       
        if ind >= 0:
                mmResp = resp
                mmInd = ind
       
        else:
                mmResp = -1
                mmInd = -1
       

        #mmout = [mmInd, mmResp, mmFlv]                                                                                                                                                
        return mmResp


    def PassJetID(self, jets):
       
        jetIDs = [j.jetId & 0b010 for j in jets if j.pt > 30]
        return (0 not in jetIDs)


    def PassNjets(self, jets):
        countJets = sum([j.Stop0l for j in jets])
        return countJets >= 2

    def GetJetSortedIdx(self, jets):
        ptlist = []
        etalist = []
        dphiMET = []
        for j in jets:
            if math.fabs(j.eta) > 4.7 or j.pt < 30:
                pass
            else:
               ptlist.append(-j.pt)
               etalist.append(math.fabs(j.eta))
               dphiMET.append(j.dPhiMET)
        sortIdx = np.lexsort((etalist, ptlist))
        return sortIdx, [dphiMET[j] for j in sortIdx]
    
    def PassdPhi(self, sortedPhi, dPhiCuts, invertdPhi =False):
        if invertdPhi:
            return any( a < b for a, b in zip(sortedPhi, dPhiCuts))
        else:
            return all( a > b for a, b in zip(sortedPhi, dPhiCuts))
        
    def PassdPhiVal(self, sortedPhi, dPhiCutsLow, dPhiCutsHigh):
        return all( (a < b and b < c) for a, b, c in zip(dPhiCutsLow, sortedPhi, dPhiCutsHigh))
    def analyze(self, event):
        jets      = Collection(event, 	"Jet")
        met       = Object(event,     "MET")
        genjets   = Collection(event,   "GenJet")
        stop0l    = Object(event,     "Stop0l")
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        isotracks = Collection(event, "IsoTrack")
        taus      = Collection(event, "Tau")
        countEle, countMu, countIsk, countTauPOG = self.calculateNLeptons(electrons, muons, isotracks, taus)
        PassElecVeto   = countEle == 0
        PassMuonVeto   = countMu == 0
        PassIsoTrkVeto = countIsk == 0
        PassTauVeto    = countTauPOG == 0
        PassLeptonVeto  = PassElecVeto and PassMuonVeto and PassIsoTrkVeto and PassTauVeto
        sortedIdx, sortedPhi = self.GetJetSortedIdx(jets)        
        jet1pt = jets[0].pt
        caloMET   = Object(event, "CaloMET")
        flags     = Object(event,     "Flag")
        PassJetID       = self.PassJetID(jets)
        PassEventFilter = self.PassEventFilter(flags)
        PassCaloMETRatio= (met.pt / caloMET.pt ) < 5 if caloMET.pt > 0 else True
        PassNjets       = self.PassNjets(jets)
        PassMET         = met.pt >= 250
        PassHT          = stop0l.HT >= 300
        PassdPhiQCD     = self.PassdPhi(sortedPhi, [0.1, 0.1, 0.1], invertdPhi =True)
        PassQCDCR       = PassEventFilter and PassJetID and PassLeptonVeto and PassNjets and PassMET and PassHT and PassdPhiQCD
        qcdjetptweightLM =1.0
        qcdjetptweightLMwq=1.0
        qcdjetptweightHM=1.0
        qcdjetptweightHMwq=1.0
        #iG =jets[0].genJetIdx
        qcdresptailweight=1.0
        qcdresptailweightHM=1.0
        mmResp = self.getQCDRespTailCorrector(jets, genjets, met)
        if PassQCDCR and stop0l.nTop == 0 and stop0l.nW == 0 and stop0l.nResolved == 0 and stop0l.Mtb < 175 and stop0l.ISRJetPt >= 200 and stop0l.METSig > 10 and PassCaloMETRatio :
            #print 'test1'
            #qcdjetptweightLMwq=self.corrhistLM1.GetBinContent(self.corrhistLM1.GetXaxis().FindBin(jet1pt))
            #qcdjetptweightLM=self.corrhistLM.GetBinContent(self.corrhistLM.GetXaxis().FindBin(jet1pt))        
            qcdresptailweight = self.corrhist.GetBinContent(self.corrhist.GetXaxis().FindBin(mmResp))
       
        #if qcdjetptweightLMwq==0: qcdjetptweightLMwq =1
        #if qcdjetptweightLM==0: qcdjetptweightLM =1
        if qcdresptailweight ==0: qcdresptailweight =1
        self.out.fillBranch("qcdRespTailWeight", qcdresptailweight)
        #self.out.fillBranch("qcdJetPtWeightLMwq", qcdjetptweightLMwq)
        #self.out.fillBranch("qcdJetPtWeightLM", qcdjetptweightLM)
        if PassQCDCR and stop0l.nJets >= 5 and stop0l.nbtags >= 1 and PassCaloMETRatio:
            #qcdjetptweightHMwq=self.corrhistHM1.GetBinContent(self.corrhistHM1.GetXaxis().FindBin(jet1pt))
            #qcdjetptweightHM=self.corrhistHM.GetBinContent(self.corrhistHM.GetXaxis().FindBin(jet1pt))
            qcdresptailweightHM = self.corrhistHM.GetBinContent(self.corrhistHM.GetXaxis().FindBin(mmResp))
        if qcdresptailweightHM ==0: qcdresptailweightHM =1
        #if qcdjetptweightHMwq==0: qcdjetptweightHMwq =1
        #if qcdjetptweightHM==0: qcdjetptweightHM =1
        self.out.fillBranch("qcdRespTailWeightHM", qcdresptailweightHM)
        #self.out.fillBranch("qcdJetPtWeightHMwq", qcdjetptweightHMwq)
        #self.out.fillBranch("qcdJetPtWeightHM", qcdjetptweightHM)
        return True
 
