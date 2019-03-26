import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class UpdateEvtWeight(Module):
    def __init__(self, isData, CrossSection, nEvent):
        self.isData = isData
        self.xs = CrossSection
        self.nEvent = nEvent

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.isData:
            infostr = "(Lumi=%f)" % self.xs
        else:
            infostr = "(CrossSection=%f, nEvent=%f)" % (self.xs, self.nEvent)
        self.out.branch("Stop0l_evtWeight",         "F", title="Storing cross section/nEvent for MC, lumi for Data" + infostr)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        sign = 1
        if not self.isData:
            initgenWeight = getattr(event, "genWeight")
            sign          = 1 if initgenWeight > 0 else -1

        neweight = self.xs/self.nEvent * sign

        ### Store output
        self.out.fillBranch("Stop0l_evtWeight",        neweight)
        return True
