import ROOT
import sys
import numpy
import argparse
from math import sqrt, fabs, sin, exp, log10
from ROOT import TFile, TTree, TChain, TBranch, TH1D, TH1I, TH1F, TH2F, TGraphAsymmErrors, TEfficiency, Math
from ROOT import TLorentzVector
from ROOT.Math import LorentzVector

parser = argparse.ArgumentParser(description="Meausre the trigger efficiency")

parser.add_argument("-f", "--file", default="", help="The path to the output tree of displaced jets analysis")
parser.add_argument("-d", "--data", default=False, action="store_true", help="CMSSW_7_4_X or CMSSW_7_6_X")
parser.add_argument("-o", "--output", default="", help="The name of the output file")
parser.add_argument("-w", "--weight", default="1.0",  help="The reweight factor of the sample")
args = parser.parse_args()

filename = args.file
weight = float(args.weight)
tfile = TFile.Open(filename)
tree = tfile.Get("TriggerNtuple/tree")

outname = args.output
fout = TFile("{0}.root".format(outname), "recreate")


eff_CalojetPt40 = TEfficiency("eff_CalojetPt40", "Filter Efficiency for Online Calo Jet", 120, 0, 120)
eff_CalojetEta  = TEfficiency("eff_CalojetEta", "Filter Efficiency v.s. Eta for Online Calo Jet", 100, -2.5, 2.5) 	
eff_PromptTrack = TEfficiency("eff_PromptTrack", "Filter Efficiency for Prompt Track Requirement", 20, -0.5, 19.5)
eff_DisplacedTrack = TEfficiency("eff_DisplacedTrack", "Filter Efficiency for Displaced Track Requirement", 20, -0.5, 19.5) 
#eff_hlt_HT_disp_430 = TEfficiency("eff_hlt_HT_dis_430", "Trigger Efficiency of HT430_displaceddijet40", 200, 100, 1000)

#eff_HT_disp_430_3p0 = TEfficiency("eff_HT_dis_430_3p0", "Trigger Efficiency of HT430_displaceddijet40_3p0", 200, 100, 1000)

nentries = tree.GetEntries()

passed_evt = 0
for entry in range(nentries):
   
    tree.GetEntry(entry)
    #if entry%1000==1:

    #HT430_filter = tree.HTfilter_fired
    if not tree.HTfilter_fired:
        continue
    if not tree.CaloHT>200:
        continue
    
#    ncalojet = tree.calojet_pt.size()
    ncalojet = min(tree.calojet_pt.size(), tree.calojet_eta.size())
    for icalojet in range(ncalojet):
        if tree.calojet_pt.at(icalojet)>60:
            eff_CalojetEta.FillWeighted(bool(tree.calojet_tagged.at(icalojet)), weight, tree.calojet_eta.at(icalojet))
        #if fabs(tree.calojet_eta.at(icalojet))<2.0:
        eff_CalojetPt40.FillWeighted(bool(tree.calojet_tagged.at(icalojet)), weight, tree.calojet_pt.at(icalojet))

    njet = min(tree.offline_jetpt.size(), tree.offline_jeteta.size())
    for ijet in range(njet):
        if not (tree.offline_jetpt.at(ijet)>60 and fabs(tree.offline_jeteta.at(ijet))<2.0):
            continue
        eff_PromptTrack.FillWeighted(bool(tree.PromptTagged.at(ijet)), weight, tree.nPromptTracks_1000.at(ijet))
        if tree.PromptTagged.at(ijet):
            eff_DisplacedTrack.FillWeighted(bool(tree.DisplacedTagged.at(ijet)), weight, tree.nDisplacedTracks_500.at(ijet)) 
    #recalculate the offline caloHT

fout.Write()
fout.Close()

