#!/usr/bin/env python
import os, sys
import ROOT
import math
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from os import system, environ
import xgboost as xgb
import logging

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

class XGBHelper:
    def __init__(self, model_file, var_list):
        self.bst = xgb.Booster(params={'nthread': 1}, model_file=model_file)
        self.var_list = var_list
        logging.info('Load XGBoost model %s, input variables:\n  %s' % (model_file, str(var_list)))

    def eval(self, inputs):
        dmat = xgb.DMatrix(np.array([[inputs[k] for k in self.var_list]]), feature_names=self.var_list)
        return self.bst.predict(dmat)[0]


class tauMVA(Module):
    def __init__(self):
	self.writeHistFile=True
	self.metBranchName = "MET"
	#still trying to find appropriate cut, but the this is the best training model
	self.tauMVADisc = 0.71
	self.bdt_file = environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_003000_maxdepth10.model"
	self.bdt_vars = ['pt', 'abseta', 'chiso0p1', 'chiso0p2', 'chiso0p3', 'chiso0p4', 'totiso0p1', 'totiso0p2', 'totiso0p3', 'totiso0p4', 'neartrkdr', 'contjetdr', 'contjetcsv']
	self.xgb = XGBHelper(self.bdt_file, self.bdt_vars)

    def beginJob(self,histFile=None,histDirName=None):
   	pass
    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("mt", 			"F", lenVar="nPFcand")
	self.out.branch("BDT_Output", 		"F", lenVar="nPFcand")
	self.out.branch("TauMVA_Stop0l",        "I")
	self.out.branch("nPFcand_save", 	"I")
	self.out.branch("PFCand_pt",		"F", lenVar="nPFcand_save")   
	self.out.branch("PFCand_abseta",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_chiso0p1",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_chiso0p2",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_chiso0p3",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_chiso0p4",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_totiso0p1",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_totiso0p2",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_totiso0p3",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_totiso0p4",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_neartrkdr",	"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_jetdr",		"F", lenVar="nPFcand_save")
	self.out.branch("PFCand_jetcsv",	"F", lenVar="nPFcand_save")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelTauMVA(self, mva):
	if mva > self.tauMVADisc:
		return True
	else:
		return False

    def getNearPhotonIndex(self, pfc, pfcands):
	minPhotonPt = 0.5
	maxPhotonDR = 0.2
	photonInd = -1
	maxPhotonPT = 0.0
	
	for ic in xrange(len(pfcands)):
		c = pfcands[ic]
		if(c.nearphopt < minPhotonPt): continue
		dr = deltaR(c.eta, c.phi, pfc.nearphoeta, pfc.nearphophi)
		if(dr > maxPhotonDR): continue;
		if(c.nearphopt > maxPhotonPT):
			maxPhotonPT = c.pt;
			photonInd = ic;
	
	return photonInd;

    def transverseMass(self, visible, invisible):
	cosDPhi   = np.cos( deltaPhi(visible.Phi(), invisible.phi) );
	return np.sqrt( 2 * visible.Pt() * invisible.pt * (1 - cosDPhi) );
        
    def computeMT(self, pfc, met, pfcands):
	photonInd = self.getNearPhotonIndex(pfc, pfcands);
	candP4 = ROOT.TLorentzVector()
	candP4.SetPtEtaPhiM(pfc.pt, pfc.eta, pfc.phi, pfc.mass)
	if(photonInd > -1): 
		pfcand_buff = ROOT.TLorentzVector()
		pfcand_buff.SetPtEtaPhiM(pfcands[photonInd].pt, pfcands[photonInd].eta, pfcands[photonInd].phi, pfcands[photonInd].mass)
		candP4+=pfcand_buff;
	return self.transverseMass(candP4, met);

    def CreateTaus(self, mva, pfc):
	mva_store = []
	if mva[self.bdt_vars[0]] > 10.0 and abs(mva[self.bdt_vars[1]]) < 2.4 and abs(pfc.dz) < 0.2:
		mva_store = mva.values()
		print mva_store
	return mva_store[0], mva_store[1], mva_store[2], mva_store[3], mva_store[4], mva_store[5], mva_store[6], mva_store[7], mva_store[8], mva_store[9], mva_store[10], mva_store[11], mva_store[12]

    def analyze(self, event):
        ## Getting objects
	met	  = Object(event, self.metBranchName)
	jets	  = Collection(event, "Jet")
	pfcand    = Collection(event, "PFcand")
	eventNum  = event.event

        pfchargedhads = []
	mva = {}
	mva_ = []
	mt_ = []
	MVASave_pt = []
	MVASave_eta = []
	MVASave_chiso0p1 = []
	MVASave_chiso0p2 = []
	MVASave_chiso0p3 = []
	MVASave_chiso0p4 = []
	MVASave_totiso0p1 = []
	MVASave_totiso0p2 = []
	MVASave_totiso0p3 = []
	MVASave_totiso0p4 = []
	MVASave_neartrkdr = []
	MVASave_contjetdr = []
	MVASave_contjetcsv = []
	for pfc in pfcand:
		mva_buff = 0.0
		mt = 0.0
		if(pfc.pt > 10.0 and abs(pfc.eta) < 2.4 and abs(pfc.dz) < 0.2):
			mt = self.computeMT(pfc, met, pfcand)

			pt 	     = min(pfc.pt,float(300.0))
			abseta       = min(abs(pfc.eta), float(2.4))
			chiso0p1     = min(pfc.chiso0p1,float(700.0))
			chiso0p2     = min(pfc.chiso0p2,float(700.0))
			chiso0p3     = min(pfc.chiso0p3,float(700.0))
			chiso0p4     = min(pfc.chiso0p4,float(700.0))
			totiso0p1    = min(pfc.totiso0p1,float(700.0))
			totiso0p2    = min(pfc.totiso0p2,float(700.0))
			totiso0p3    = min(pfc.totiso0p3,float(700.0))
			totiso0p4    = min(pfc.totiso0p4,float(700.0))
			neartrkdr    = pfc.nearestTrkDR
			jetmatch     = (pfc.contJetIndex > -1) and (jets[pfc.contJetIndex].pt >= 20.0) and (abs(jets[pfc.contJetIndex].eta) < 2.4)
			jetdr        = deltaR(jets[pfc.contJetIndex].eta, jets[pfc.contJetIndex].phi, pfc.eta, pfc.phi) if jetmatch else -1.0
			jetcsv       = jets[pfc.contJetIndex].btagDeepB if jetmatch else -1.0
			
			contjetdr  = min(float(0.4), jetdr)
			if(contjetdr < 0.0): contjetdr = 0.0
			contjetcsv =  jetcsv
			if(contjetcsv < 0.0): contjetcsv = 0.0
	
			mva = {self.bdt_vars[0]: pt, 
			       self.bdt_vars[1]: abseta,
			       self.bdt_vars[2]: chiso0p1, 
			       self.bdt_vars[3]: chiso0p2, 
			       self.bdt_vars[4]: chiso0p3, 
			       self.bdt_vars[5]: chiso0p4, 
			       self.bdt_vars[6]: totiso0p1, 
			       self.bdt_vars[7]: totiso0p2, 
			       self.bdt_vars[8]: totiso0p3, 
			       self.bdt_vars[9]: totiso0p4, 
			       self.bdt_vars[10]: neartrkdr, 
			       self.bdt_vars[11]: contjetdr, 
			       self.bdt_vars[12]: contjetcsv}
			pt, eta, chiso0p1, chiso0p2, chiso0p3, chiso0p4, totiso0p1, totiso0p2, totiso0p3, totiso0p4, neartrkdr, contjetdr, contjetcsv = self.CreateTaus(mva, pfc)
			mva_buff = self.xgb.eval(mva)
		MVASave_pt.append(pt)
		MVASave_eta.append(eta)
		MVASave_chiso0p1.append(chiso0p1)
		MVASave_chiso0p2.append(chiso0p2)
		MVASave_chiso0p3.append(chiso0p3)
		MVASave_chiso0p4.append(chiso0p4)
		MVASave_totiso0p1.append(totiso0p1)
		MVASave_totiso0p2.append(totiso0p2)
		MVASave_totiso0p3.append(totiso0p3)
		MVASave_totiso0p4.append(totiso0p4)
		MVASave_neartrkdr.append(neartrkdr)
		MVASave_contjetdr.append(contjetdr)
		MVASave_contjetcsv.append(contjetcsv)
		mt_.append(mt)
		mva_.append(mva_buff)
	self.TauMVA_Stop0l = map(self.SelTauMVA, mva_)

	#print "mva output: ", mva_
	self.out.fillBranch("mt", 		mt_)
        self.out.fillBranch("BDT_Output", 	mva_)
	self.out.fillBranch("TauMVA_Stop0l", sum(self.TauMVA_Stop0l))
	self.out.fillBranch("PFCand_pt"		, MVASave_pt)   
	self.out.fillBranch("PFCand_abseta"	, MVASave_eta)
	self.out.fillBranch("PFCand_chiso0p1"	, MVASave_chiso0p1)
	self.out.fillBranch("PFCand_chiso0p2"	, MVASave_chiso0p2)
	self.out.fillBranch("PFCand_chiso0p3"	, MVASave_chiso0p3)
	self.out.fillBranch("PFCand_chiso0p4"	, MVASave_chiso0p4)
	self.out.fillBranch("PFCand_totiso0p1"	, MVASave_totiso0p1)
	self.out.fillBranch("PFCand_totiso0p2"	, MVASave_totiso0p2)
	self.out.fillBranch("PFCand_totiso0p3"	, MVASave_totiso0p3)
	self.out.fillBranch("PFCand_totiso0p4"	, MVASave_totiso0p4)
	self.out.fillBranch("PFCand_neartrkdr"	, MVASave_neartrkdr)
	self.out.fillBranch("PFCand_jetdr"	, MVASave_contjetdr)
	self.out.fillBranch("PFCand_jetcsv"	, MVASave_contjetcsv)
		
	return True
