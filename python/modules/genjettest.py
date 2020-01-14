#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
#import "$CMSSW_BASE/src/AnalysisTools/QuickRefold/interface/TObjectContainer.h"

class genjettest(Module): 
    def __init__(self):
	self.writeHistFile=True
	self.metBranchName="MET"

    h1 = ROOT.TH1F("h1"," genjets", 100, 0, 1000)
    h2 = ROOT.TH1F("h2"," recojets", 100, 0, 1000)

    def analyze1(self, event):
        jets      = Collection(event, "Jet")
        genjets   = Collection(event, "GenJet")
        met       = Object(event,     self.metBranchName)
        for gJ in genjets:
            #print 'pt of genjets',gJ.pt
            self.h1.Fill(gJ.pt)
            for rJ in jets:
                if(deltaR(gJ.eta,gJ.phi,rJ.eta,rJ.phi)< 0.2):
                   #print 'pt of recojets',rJ.pt
                   self.h2.Fill(rJ.pt)
        return True
    c1 = ROOT.TCanvas('c1','The Ntuple canvas',800,800)
    c1.cd()
    print '******'
    h1.Draw("")
    h2.Draw("same")
    c1.SaveAs('genreco_jet_map.png')
   
   
  
    



    
