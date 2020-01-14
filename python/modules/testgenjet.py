#!/usr/bin/env python                                                    
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from ROOT import TCanvas, TPad, TFile


file = ROOT.TFile.Open('../ttbar_v1_tree.root')

c1 = ROOT.TCanvas('c1','The Ntuple canvas',800,800)
h1 = ROOT.TH1F("h1"," genjets", 100, 0, 600)
h2 = ROOT.TH1F("h2"," recojets", 100, 0, 600)
h3 = ROOT.TH1F("h3","genIdxpt",100,0,600)
tree = file.Get('Events')

for event in tree:
 genptvec = event.GenJet_pt
 genetavec = event.GenJet_eta
 genphivec = event.GenJet_phi
 jetptvec = event.Jet_pt
 jetetavec = event.Jet_eta
 jetphivec = event.Jet_phi
 genjetidvec = event.Jet_genJetIdx

 for j in xrange(0, len(genptvec)):
   h1.Fill(genptvec[j])
   for r in xrange(0,len(jetptvec)):
     if deltaR(genetavec[j],genphivec[j],jetetavec[r],jetphivec[r])< 0.2 :
       h2.Fill(jetptvec[r])
 for i in xrange(0,len(genjetidvec)):
   print 'genjet Index value is',genjetidvec[i],'length of genjet is',len(genptvec) 
  # print 'gen eta' ,genetavec[i],'gen phi ',genphivec[i],'jet eta',jetetavec[i],'jet phi',jetphivec[i]
   if genjetidvec[i] >= len(genptvec): continue
   if genjetidvec[i]==-1:continue
   print 'genjet eta:',genetavec[genjetidvec[i]],'genjet phi:',genphivec[genjetidvec[i]],'jet eta:',jetetavec[i],'jet phi:',jetphivec[i]
   h3.Fill(genptvec[genjetidvec[i]])
   
c1.cd()
pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.Draw()
pad1.cd()
pad1.SetLogy()
h1.Draw("")
h2.Draw("same")
h3.Draw("same")
h1.SetLineColor(4)
h2.SetLineColor(6)
h3.SetLineColor(8)
h1.SetLineWidth(3)
h2.SetLineWidth(3)
h3.SetLineWidth(3)
h1.SetTitle("genpt and recopt mapping;pt;# events")
h1.SetLabelSize(15,"x")
h1.SetLabelFont(43,"x")
l1=ROOT.TLegend(.6,.7,.8,.9,"")
l1.AddEntry(h1,"genjet pt","l")
l1.AddEntry(h2,"recojet pt","l")
l1.AddEntry(h3,"genIdx pt","l")
l1.Draw("same")
c1.cd()
pad2 = TPad("pad2", "pad2", 0, 0.001, 1, 0.3)
pad2.Draw()
pad2.cd()
hr1 = h1.Clone("hr1")
hr1.Sumw2()
hr1.Divide(h2)
hr1.SetMarkerStyle(21)
hr1.SetTitle("ratio")
hr1.SetMarkerColor(8)
hr1.SetAxisRange(0,2,"y")
hr1.SetTitleFont(75,"y")
hr1.SetLabelSize(0.08,"y")
hr1.SetStats(0)
hr1.Draw("ep")
line1 = ROOT.TLine(0,1,600,1)
line1.SetLineColor(1)
line1.Draw("same")
c1.Update()
c1.Print("ptcheck.png")
# pad1 = TPad("pad3", "pad3", 0, 0.3, 1, 1.0);
# pad1->Draw();
# pad1->cd();
# pad1->SetLogy();
# h1->Draw("");
# h2->Draw("same");
# h1->SetLineColor(kRed);
# h2->SetLineColor(kGreen);
# h1->SetLineWidth(3);
# h2->SetLineWidth(3);
# h1->SetTitle("jet response curve (second leading jet);recojet / genjet; # of events");
# leg1 = new TLegend(.6,.7,.8,.9,"");
# leg1->SetFillColor(0);
# leg1->AddEntry(hres1H,"HEM pt","l");
# leg1->AddEntry(hres1nH,"NoHEM pt","l");
# leg1->Draw("same");
# c2->cd();
#     TPad *pad4 = new TPad("pad4", "pad4", 0, 0.01, 1, 0.3);
#     pad4->Draw();
#     pad4->cd();
#     TH1D *hratio2 = (TH1D*)hres1H->Clone("hratio2");
#     hratio2->Divide(hres1nH);
#     hratio2->SetMarkerStyle(3);
#     hratio2->SetMarkerColor(kRed);
#     hratio2->SetAxisRange(0,2,"y");
#     hratio2->SetLabelSize(0.07,"Y");
#     hratio2->Draw("same");
#     TLine *line2 = new TLine(0,1,2,1);
#     line2->SetLineColor(kBlack);
#     line2->Draw("same");
#     c2->SaveAs("HEMres1_pt.png");

