import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np
from functools import reduce
import operator

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

DeepCSVLooseWP = {
    "2016" : 0.2217,
    "2017" : 0.1522,
    "2018" : 0.1241
}

DeepCSVMediumWP ={
    "2016" : 0.6324,
    "2017" : 0.4941,
    "2018" : 0.4184
}

CSVv2MediumWP = {
    "2016" : 0.8484,
    "2017" : 0.8838,
    "2018" : 0.8838  # Not recommended, use 2017 as temp
}


class LLObjectsProducer(Module):
    def __init__(self, era, isData = False):
        self.era = era
	self.isData = isData
        self.metBranchName = "MET"
        # EE noise mitigation in PF MET
        # https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1865.html
        if self.era == "2017":
            self.metBranchName = "METFixEE2017"

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("Stop0l_nbtags_Loose",  		"I")
	self.out.branch("Stop0l_MtLepMET", 			"F")
	self.out.branch("Stop0l_nVetoElecMuon", 		"I")
	self.out.branch("Stop0l_nVetoMuon",                 "I")
        self.out.branch("Stop0l_nVetoElec",                 "I")
        self.out.branch("Stop0l_noMuonJet",			"O")
	self.out.branch("Pass_dPhiQCD",				"O")
	self.out.branch("Pass_dPhiQCDSF",			"O")
	self.out.branch("Stop0l_dPhiISRMET",			"F")
	if not self.isData :
		self.out.branch("ElecSF",			"F")
		self.out.branch("ElecSFErr",		"F")
		self.out.branch("MuSF",			"F")
		self.out.branch("MuSFErr",			"F")
		self.out.branch("TaSF",			"F")
                print 'what problem is here'
         	# self.out.branch("TauSFUp",			"F")
		# self.out.branch("TauSFDown",			"F")
		# self.out.branch("WSF",				"F")
		# self.out.branch("WSFErr",			"F")
		# self.out.branch("TopSF",			"F")
		# self.out.branch("TopSFErr",			"F")
		# self.out.branch("restopSF",			"F")
		# self.out.branch("SoftBSF",			"F")
		# self.out.branch("SoftBSFErr",			"F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelMtlepMET(self, ele, muon, isks, met):
	mt = 0.0
	for l in ele:
		if l.Stop0l: mt += math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi))))
	for l in muon:
		if l.Stop0l: mt += math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi))))
	for l in isks:
		if l.Stop0l: mt += math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi))))
	return mt

    def SelJets(self, jet):
        if jet.pt < 20 or math.fabs(jet.eta) > 2.4 :
            return False
        return True

    def SelBtagJets(self, jet):
        global DeepCSVLooseWP
        if jet.btagDeepB >= DeepCSVLooseWP[self.era]:
            return True
        return False

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

    def SelNoMuon(self, jets, met):
	noMuonJet = True
	for j in jets:
		if j.pt > 200 and j.muEF > 0.5 and abs(deltaPhi(j.phi, met.phi)) > (math.pi - 0.4):
			noMuonJet = False
	return noMuonJet

    def ScaleFactorErrElectron(self, obj):
	sf = 1
	sfErr = 0
	for s in obj:
		if not s.Stop0l: continue
	        sf *= s.VetoSF
                sfErr += ((s.VetoSFErr)**2)

	return sf, math.sqrt(sfErr)

    def ScaleFactorErrMuon(self, obj):
	sf = 1
	sfErr = 0
	for s in obj:
		if not s.Stop0l: continue
		sf *= s.MediumSF
                sfErr += ((s.MediumSFErr)**2)
     	return sf, math.sqrt(sfErr)

    def ScaleFactorErrTau(self, obj):
	sf = 1
	sfUp = 0
	sfDown = 0
	for s in obj:
		if not s.Stop0l: continue
		sf *= s.MediumSF
		sfUp += ((s.MediumSF_Up)**2)*((s.MediumSF)**2)
		sfDown += ((s.MediumSF_Down)**2)*((s.MediumSF)**2)

	return sf, math.sqrt(sfUp), math.sqrt(sfDown)

    def ScaleFactorErrFatjet(self, obj, objType):
	sf = 1
	sfErr = 0
	for s in obj:
		if not (s.Stop0l == objType): continue
		sf *= s.SF
		sfErr += ((s.SFerr)**2)*((s.SF)**2)

	return sf, math.sqrt(sfErr)

    def ScaleFactorErrResTop(self, obj, objType):
	sf = 1
	for s, rT in zip(obj, objType):
		if not rT.Stop0l: continue
		sf *= s.sf

	return sf

    def ScaleFactorErrSoftB(self, obj):
	sf = 1
	sfErr = 0
	for s in obj:
		if not s.Stop0l: continue
		sf *= s.SF
		sfErr += ((s.SFerr)**2)*((s.SF)**2)

	return sf, math.sqrt(sfErr)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        isotracks = Collection(event, "IsoTrack")
	jets      = Collection(event, "Jet")
	met       = Object(event, self.metBranchName)
	taus	  = Collection(event, "Tau")
	stop0l    = Object(event, "Stop0l")
	fatjets   = Collection(event, "FatJet")
	restop	  = Collection(event, "ResolvedTopCandidate")
	res	  = Collection(event, "ResolvedTop", lenVar="nResolvedTopCandidate")
	SB	  = Collection(event, "SB")

        ## Selecting objects
	self.Jet_Stop0l      = map(self.SelJets, jets)
	local_BJet_Stop0l    = map(self.SelBtagJets, jets)
        self.BJet_Stop0l     = [a and b for a, b in zip(self.Jet_Stop0l, local_BJet_Stop0l )]
	mt 		     = sum([ e.MtW for e in electrons if e.Stop0l ] + [ m.MtW for m in muons if m.Stop0l ])
	countEle	     = sum([e.Stop0l for e in electrons])
	countMuon	     = sum([m.Stop0l for m in muons])
	noMuonJet	     = self.SelNoMuon(jets, met)
	sortedIdx, sortedPhi = self.GetJetSortedIdx(jets)
	PassdPhiQCD          = self.PassdPhi(sortedPhi, [0.1, 0.1, 0.1], invertdPhi =True)
	PassdPhiQCDSF        = self.PassdPhi(sortedPhi, [0.1, 0.1], invertdPhi =True)
	dphiISRMet	     = abs(deltaPhi(fatjets[stop0l.ISRJetIdx].phi, met.phi)) if stop0l.ISRJetIdx >= 0 else -1

	if not self.isData :
		electronSF, electronSFErr = self.ScaleFactorErrElectron(electrons)
		muonSF, muonSFErr    = self.ScaleFactorErrMuon(muons)
		tauSF, tauSFUp, tauSFDown = self.ScaleFactorErrTau(taus)
		## type top = 1, W = 2, else 0
		WSF, WSFErr	     = self.ScaleFactorErrFatjet(fatjets, 2)
		topSF, topSFErr	     = self.ScaleFactorErrFatjet(fatjets, 1)
		resSF		     = self.ScaleFactorErrResTop(restop, res)
		softSF, softSFErr    = self.ScaleFactorErrSoftB(SB)

        ### Store output
	self.out.fillBranch("Stop0l_nbtags_Loose",   	sum(self.BJet_Stop0l))
	self.out.fillBranch("Stop0l_MtLepMET",  	mt)
	self.out.fillBranch("Stop0l_nVetoElecMuon", 	countEle + countMuon)
        self.out.fillBranch("Stop0l_nVetoElec",     countEle )
        self.out.fillBranch("Stop0l_nVetoMuon",     countMuon)
	self.out.fillBranch("Stop0l_noMuonJet",		noMuonJet)
	self.out.fillBranch("Pass_dPhiQCD",		PassdPhiQCD)
	self.out.fillBranch("Pass_dPhiQCDSF",		PassdPhiQCDSF)
	self.out.fillBranch("Stop0l_dPhiISRMET",	dphiISRMet)
	
	if not self.isData :
		self.out.fillBranch("ElecSF",	electronSF)
		self.out.fillBranch("ElecSFErr",	electronSFErr)
		self.out.fillBranch("MuSF",		muonSF)
		self.out.fillBranch("MuSFErr",	muonSFErr)
		self.out.fillBranch("TaSF",		tauSF)
		#self.out.fillBranch("TauSFUp",		tauSFUp)
		#self.out.fillBranch("TauSFDown",	tauSFDown)
		# self.out.fillBranch("WSF",		WSF)
		# self.out.fillBranch("WSFErr",		WSFErr)
		# self.out.fillBranch("TopSF",		topSF)
		# self.out.fillBranch("TopSFErr",		topSFErr)
		# self.out.fillBranch("restopSF",		resSF)
		# self.out.fillBranch("SoftBSF",		softSF)
		# self.out.fillBranch("SoftBSFErr",	softSFErr)
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
