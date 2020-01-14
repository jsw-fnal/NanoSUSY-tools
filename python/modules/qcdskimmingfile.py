#!/usr/bin/env python
import os, sys
import ROOT
import math
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from functools import reduce
import operator

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoSUSYTools.modules.Stop0lObjectsProducer import DeepCSVMediumWP, DeepCSVLooseWP
#import "$CMSSW_BASE/src/AnalysisTools/QuickRefold/interface/TObjectContainer.h"


DeepCSVMediumWP ={
    "2016" : 0.6324,
    "2017" : 0.4941,
    "2018" : 0.4184
}

class qcdskimmingfile(Module): 

    def __init__(self, era,isData = False):
## WP from Hui's study https://indico.cern.ch/event/780000/contributions/3248659/attachments/1768782/2873013/Search_bin_study_with_combine_tools_v13.pdf## Updated WP from https://indico.cern.ch/event/840827/contributions/3527925/attachments/1895214/3126510/DeepAK8_Top_W_SFs_2017_JMAR_PK.pdf
	self.writeHistFile=True
        self.era = era
        self.isData = isData
	self.metBranchName="MET"
        if self.era == "2017":
            self.metBranchName = "METFixEE2017"
        self.nBootstraps = 50
        self.minAK8TopMass = 105
        self.maxAK8TopMass = 210
        self.DeepAK8TopPt  = 400.0 # Using 400 Pt cut 
        ## Mistag 0.5% WP, using 2017 WP as 2018
        self.DeepAK8TopWP  = {
            "2016" : 0.937,
            "2017" : 0.895,
            "2018" : 0.895
        }

        ## Mistag 0.5% WP, using 2017 WP as 2018
        self.minAK8WMass   = 65
        self.maxAK8WMass   = 105
        self.DeepAK8WPt    = 200.0
        self.DeepAK8WWP    = {
            "2016" : 0.973,
            "2017" : 0.991,
            "2018" : 0.991
        }

        self.DeepResolveWP = 0.92
        self.etaMax        = 2.0
        self.bJetEtaMax    = 2.4
        self.resAK4bTagWP  = DeepCSVMediumWP[era]
        self.dR2AK4Subjet  = 0.4*0.4
        self.era = era 
       # if self.era == "2017":
           # self.metBranchName = "METFixEE2017"   
    def beginJob(self,histFile=None,histDirName=None):
   	pass
    def endJob(self):
	pass 

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nJetPass", "I")
        self.out.branch("JetPass_pt","F",lenVar="nJetPass") 
        self.out.branch("JetPass_eta","F",lenVar="nJetPass") 
        self.out.branch("JetPass_phi", "F",lenVar="nJetPass")
        self.out.branch("JetPass_mass", "F",lenVar="nJetPass")
        self.out.branch("Stop0l_HT30",                          "F")
        self.out.branch("Pass_HT30",                            "O")
        #self.out.branch("Pass_NJets30",                         "O")
        self.out.branch("Stop0l_nJets30",                       "I")
        self.out.branch("Pass_dPhiMET30",                       "O")
        self.out.branch("Pass_dPhiMETLowDM30",                  "O")
        self.out.branch("Pass_dPhiMETMedDM30",                  "O")
        self.out.branch("Pass_dPhiMETHighDM30",                 "O")
        self.out.branch("Pass_dPhiQCD30",                         "O")
        self.out.branch("Pass_dPhiQCDSF30",                       "O")
        self.out.branch("Stop0l_nJets30_2p2",                     "I")
        self.out.branch("Stop0l_nJets20_2p2",                     "I")
        if not self.isData:
            self.out.branch("topTagWeight","F")
            self.out.branch("wTagWeight","F")
            self.out.branch("restopTagWeight","F")
            self.out.branch("softbTagWeight","F")
            self.out.branch("wSFErr", "F")
            self.out.branch("topSFErr", "F")
            self.out.branch("softbSFErr","F")
    

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def jetResFunction(self, jets, genjets):
	res = jets.pt/genjets.pt
	return res
    def SelJets(self, jet):
        if jet.pt < 20 or math.fabs(jet.eta) > 2.4 :
            return False
        return True

    def SelJets2p2(self, jet):
         if jet.pt < 20 or math.fabs(jet.eta) > 2.2 :
             return False
         return True
    def SelJets302p2(self, jet):
         if jet.pt < 30 or math.fabs(jet.eta) > 2.2 :
             return False
         return True

    def PassNjets(self, jets, jetpt = 20.):
        countJets = sum([j.Stop0l for j in jets if j.pt >= jetpt])
        return countJets
    def PassdPhiVal(self, sortedPhi, dPhiCutsLow, dPhiCutsHigh):
        return all( (a < b and b < c) for a, b, c in zip(dPhiCutsLow, sortedPhi, dPhiCutsHigh))

    def CalHT(self, jets, jetpt=20):
        HT = sum([j.pt for i, j in enumerate(jets) if (self.Jet_Stop0l[i] and j.pt > jetpt)])
        return HT
    def GetJetSortedIdx(self, jets, jetpt = 20, jeteta = 4.7):
        ptlist = []
        etalist = []
        dphiMET = []
        for j in jets:
            if math.fabs(j.eta) > jeteta or j.pt < jetpt:
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
    def ScaleFactorErrFatjet(self, obj, objType):
        sf = 1
        sfErr = 0
        for s in obj:
            if not (s.Stop0l == objType): continue
            sf *= s.SF
            sfErr += ((s.SFerr)**2)*((s.SF)**2)
            
        return sf, math.sqrt(sfErr)
        
    def ScaleFactorErrSoftB(self, obj):
        sf = 1
        sfErr = 0
        for s in obj:
            if not s.Stop0l: continue
            sf *= s.SF
            sfErr += ((s.SFerr)**2)*((s.SF)**2)
            
        return sf, math.sqrt(sfErr)

    def parselist(self, array, which):
        listx_, listy_ = which, len(array)
        out = []
        for i in xrange(listy_):
            out.append(array[i][listx_])
        return out
    def analyze(self, event):
	jets      = Collection(event, "Jet")
	#Npv   = Collection(event, "PV")
        if self.isData == False:
            genjets     = Collection(event, "GenJet")
	met       = Object(event,     self.metBranchName)
        stop0l = Object(event, "Stop0l")
        fatjets  = Collection(event, "FatJet")
        resolves  = Collection(event, "ResolvedTopCandidate")
        restop  = Collection(event, "ResolvedTopCandidate")
        res  = Collection(event, "ResolvedTop", lenVar="nResolvedTopCandidate")
        SB  = Collection(event, "SB")

        #for n in Npv:    
         #   npv = n.npvsGood
        # npv=event.PV_npvsGood      
        # #print 'npv: ',npv
        # i_npvweight=[]
        # i_npvweight.append(192.076)
        # i_npvweight.append(2.70278)
        #..................      
        njets = len(jets)
        sortedPhi = self.GetJetSortedIdx(jets)
        self.Jet_Stop0l = map(self.SelJets, jets)
        self.Jet_Stop0l2p2 = map(self.SelJets2p2,jets)
        self.Jet_Stop0l302p2 = map(self.SelJets302p2,jets)
        jet_pass =[]
        for i in xrange(len(jets)):
               jet_ =[]
               j=jets[i]
               if self.Jet_Stop0l[i] == 1:
                   jet_.append(j.pt)
                   jet_.append(j.eta)
                   jet_.append(j.phi)
                   jet_.append(j.mass)
                   jet_pass.append(jet_)
                   break
        PassHT30             = self.CalHT(jets, 30.) >= 300
        HT30                 = self.CalHT(jets, 30.)
        PassNjets30          = self.PassNjets(jets, 30.) >= 2
        nJets30              = self.PassNjets(jets, 30.)
        sortedIdx, sortedPhi = self.GetJetSortedIdx(jets, 30.)
        PassdPhiLowDM30      = self.PassdPhi(sortedPhi, [0.5, 0.15, 0.15])
        PassdPhiMedDM30      = self.PassdPhiVal(sortedPhi, [0.15, 0.15, 0.15], [0.5, 4., 4.]) #Variable for LowDM Validation bins
        PassdPhiHighDM30     = self.PassdPhi(sortedPhi, [0.5, 0.5, 0.5, 0.5])
        PassdPhiQCD30        = self.PassdPhi(sortedPhi, [0.1, 0.1, 0.1], invertdPhi =True)
        PassdPhiQCDSF30      = self.PassdPhi(sortedPhi, [0.1, 0.1], invertdPhi =True)
        if not self.isData:
            i_topTagWeight,i_topSFErr = self.ScaleFactorErrFatjet(fatjets, 1)
            i_wTagWeight,i_wSFErr = self.ScaleFactorErrFatjet(fatjets, 2) 
            i_restopTagWeight = reduce(operator.mul, (restop[rt].sf for rt in xrange(len(restop)) if res[rt].Stop0l), 1)
            i_softbTagWeight, i_softbSFErr    = self.ScaleFactorErrSoftB(SB)
            self.out.fillBranch("topTagWeight",i_topTagWeight)
            self.out.fillBranch("wTagWeight",i_wTagWeight)
            self.out.fillBranch("restopTagWeight",i_restopTagWeight)
            self.out.fillBranch("softbTagWeight",i_softbTagWeight)
            self.out.fillBranch("topSFErr",i_topSFErr)
            self.out.fillBranch("wSFErr",i_wSFErr)
            self.out.fillBranch("softbSFErr",i_softbSFErr)
        self.out.fillBranch("Stop0l_nJets20_2p2",sum(self.Jet_Stop0l2p2))
        self.out.fillBranch("Stop0l_nJets30_2p2",sum(self.Jet_Stop0l302p2)  )
        self.out.fillBranch("nJetPass", len(jet_pass))
        self.out.fillBranch("JetPass_pt", self.parselist(jet_pass, 0))
        self.out.fillBranch("JetPass_eta", self.parselist(jet_pass, 1))
        self.out.fillBranch("JetPass_phi", self.parselist(jet_pass, 2))
        self.out.fillBranch("JetPass_mass", self.parselist(jet_pass,3))
        self.out.fillBranch("Stop0l_HT30",              HT30)
        self.out.fillBranch("Pass_HT30",                PassHT30)
        #self.out.fillBranch("Pass_NJets30",             PassNjets30)
        self.out.fillBranch("Stop0l_nJets30",           nJets30)
        self.out.fillBranch("Pass_dPhiMET30",           PassdPhiLowDM30)
        self.out.fillBranch("Pass_dPhiMETLowDM30",      PassdPhiLowDM30)
        self.out.fillBranch("Pass_dPhiMETMedDM30",      PassdPhiMedDM30)
        self.out.fillBranch("Pass_dPhiMETHighDM30",     PassdPhiHighDM30)
        self.out.fillBranch("Pass_dPhiQCD30",         PassdPhiQCD30)
        self.out.fillBranch("Pass_dPhiQCDSF30",       PassdPhiQCDSF30)

       #print 'no of primary vertex',npv
        #if npv >=3 and npv <=120:  
            #self.out.fillBranch("npvweight",i_npvweight[npv-3])
            #print 'no of primary vertex',npv
            #print 'weight value is ',i_npvweight[npv-3]
        #else: 
         #   self.out.fillBranch("npvweight",1)

        return True
