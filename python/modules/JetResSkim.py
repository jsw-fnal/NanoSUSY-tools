#!/usr/bin/env python                                                                                                                                   
import os, sys
import ROOT
import math
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol


class JetResSkim(Module):
    def __init__(self,era):
	pass
    def beginJob(self):
        pass
    def endjob(self):
        pass
    def beginFile(self,inputFile,outputFile,inputTree,wrappedOutputTree):
	self.out = wrappedOutputTree
	self.out.branch("weight"	,"F")
	self.out.branch("genjetpt"	,"F")
	self.out.branch("genjeteta"	,"F")
	self.out.branch("recojetpt"	,"F")
	self.out.branch("genjetrank"	,"I")
	self.out.branch("flavor"	,"I")
	self.out.branch("rechempt"	,"F")
	self.out.branch("genhempt"	,"F")

    def analyze(self, event):
	jets        = Collection(event, "Jet")
	genjets     = Collection(event, "GenJet")
	weight      = event.Stop0l_evtWeight
	eventNum    = event.event
	PassFilter  = event.Pass_EventFilter
	PassJetID   = event.Pass_JetID

	if PassFilter and PassJetID:
		for gJ in xrange(len(genjets)):
			gJet = genjets[gJ]
			if gJet.pt < 30: continue  
			
			rJet = 0
			for iR in xrange(len(jets)) :
				if jets[iR].genJetIdx != gJ: continue
				rJet = jets[iR]
				break
			
			self.out.fillBranch("weight",	weight)
			self.out.fillBranch("genjetpt", 	gJet.pt)
			self.out.fillBranch("genjeteta", 	gJet.eta)
			self.out.fillBranch("recojetpt", 	rJet.pt if rJet != 0 else 9.5)
			self.out.fillBranch("genjetrank", 	min(gJ, 250))
			self.out.fillBranch("flavor",	gJet.partonFlavour)
		
			if(gJet.eta > -2.8 and gJet.eta < -1.6 and gJet.phi >-1.37 and gJet.phi < -1.07):
				self.out.fillBranch("genhempt", gJet.pt)
			else:
				self.out.fillBranch("genhempt", 0)
			if rJet != 0:
				if(rJet.eta > -2.8 and rJet.eta < -1.6 and rJet.phi >-1.37 and rJet.phi < -1.07): 
					self.out.fillBranch("rechempt",rJet.pt)
				else:
					self.out.fillBranch("rechempt", 0)
			else:
				self.out.fillBranch("rechempt", 0)

			self.out.fill()	
	return True  
