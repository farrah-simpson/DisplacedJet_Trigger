import sys, os
import argparse

sys.path.append('/Users/jingyuluo/Downloads/root_dir/lib')

import ROOT

from ROOT import TFile, TCanvas, TH1F
from ROOT import gStyle
from ROOT import kGreen, kRed, kBlue, kBlack
import ctypes

gStyle.SetOptStat(ROOT.kFALSE)
parser=argparse.ArgumentParser()

parser.add_argument("-d", "--data1", default="test_Run2023D.root", help="The rootfile for 2023 data") 
parser.add_argument("-d2", "--data2", default="DisplacedJets_Run2018All_Trigger.root", help="The rootfile for 2018 data")
parser.add_argument("-m", "--mc", default="Trigger_OffSel_QCD_HT_Merged_New_92X.root", help="The rootfile for Run2 MC")

args = parser.parse_args()

dataname = args.data1
#data18name = args.data2
mcname  = args.mc
if "Run2023C" in dataname: run = "Run2023C"
elif "Run2023D" in dataname: run = "Run2023D"

can = TCanvas("can", "can", 800, 800)
can.cd()
can.SetTickx()
can.SetTicky()
#can.SetLogy()
#can.SetGrid(1,1)
can.Divide(1,2)

pad1_1 = can.cd(1)
pad1_1.SetPad(0.0, 0.3, 1.0, 1.0)   # xlow, ylow, xup, yup
pad1_1.SetBottomMargin(0.02)        # reduce x-axis space

pad2_1 = can.cd(2)
pad2_1.SetPad(0.0, 0.0, 1.0, 0.3)
pad2_1.SetTopMargin(0.05)
pad2_1.SetBottomMargin(0.3)

datafile = TFile(dataname)
#data18file = TFile(data18name)
mcfile = TFile(mcname)

eff_jetpt = datafile.Get("eff_CalojetPt40")
#eff_jetpt_2018 = data18file.Get("eff_CalojetPt40")
eff_jetpt_mc = mcfile.Get("eff_CalojetPt40")
pad1_1.cd()

gra_jetpt = eff_jetpt.CreateGraph()
#gra_jetpt_2018 = eff_jetpt_2018.CreateGraph()
gra_jetpt_mc = eff_jetpt_mc.CreateGraph()
gra_jetpt.SetTitle("")
gra_jetpt.SetMarkerColor(kBlack)
gra_jetpt.SetMarkerStyle(ROOT.kFullCircle)
gra_jetpt.SetMarkerSize(0.75)
#gra_jetpt_2018.SetMarkerColor(kBlue)
#gra_jetpt_2018.SetMarkerStyle(ROOT.kFullCircle)
#gra_jetpt_2018.SetMarkerSize(0.75)
#gra_jetpt_2018.SetLineColor(kBlue)
gra_jetpt_mc.SetMarkerColor(kRed)
gra_jetpt_mc.SetMarkerStyle(ROOT.kFullCircle)
gra_jetpt_mc.SetMarkerSize(0.75)
gra_jetpt_mc.SetLineColor(kRed)
#gra_jetpt.GetXaxis().SetTitle("Offline PF jet p_{T} [GeV]")
#gra_PromptTrk.GetYaxis().SetRangeUser(0, 1.0)
gra_jetpt.GetYaxis().SetTitle("Online Single Jet Tag Efficiency")
gra_jetpt.Draw("APE")
#gra_jetpt_2018.Draw("PE SAME")
gra_jetpt_mc.Draw("PE SAME")

gra_jetpt.GetXaxis().SetLabelSize(0)
gra_jetpt.GetXaxis().SetTitleSize(0)

text=ROOT.TLatex(0.65, 0.92, "2023 ")
text.SetNDC()
text.SetTextFont(62)
text.SetTextSize(0.05)
text2=ROOT.TLatex(0.15, 0.92, "CMS #bf{#scale[0.75]{#it{Preliminary}}}")
text2.SetNDC()
text2.SetTextSize(0.05)
text2.SetTextFont(62)
text.Draw("SAME")
text2.Draw("SAME")
leg=ROOT.TLegend(0.3,0.65, 0.55,0.8)
leg.SetLineColor(0)
leg.SetFillColorAlpha(0,0)

leg.AddEntry(gra_jetpt, "Run2023D", "pe")
#leg.AddEntry(gra_jetpt_2018, "Run2018", "pe")
leg.AddEntry(gra_jetpt_mc, "W+Jets selection MC", "pe")
#leg.AddEntry(gra_PromptTrk_mc, "MC (QCD MultiJet)", "pe")
leg.Draw("SAME")

# --- Build ratio_jetpt graph (data / MC) with errors ---
n_data = gra_jetpt.GetN()
n_mc   = gra_jetpt_mc.GetN()
npoints = min(n_data, n_mc)

ratio_jetpt = ROOT.TGraphErrors(npoints)
ratio_jetpt.SetName("ratio_jetpt")

for i in range(npoints):
    # Prepare containers for GetPoint
    x_d = ctypes.c_double()
    y_d = ctypes.c_double()
    x_m = ctypes.c_double()
    y_m = ctypes.c_double()

    gra_jetpt.GetPoint(i, x_d, y_d)
    gra_jetpt_mc.GetPoint(i, x_m, y_m)

    xd_val, yd_val = x_d.value, y_d.value
    xm_val, ym_val = x_m.value, y_m.value

    # Errors
    ey_d = gra_jetpt.GetErrorY(i)
    ey_m = gra_jetpt_mc.GetErrorY(i)

    if ym_val > 0 and yd_val > 0:
        R = yd_val / ym_val
        err = R * ((ey_d / yd_val)**2 + (ey_m / ym_val)**2)**0.5
        ratio_jetpt.SetPoint(i, xd_val, R)
        ratio_jetpt.SetPointError(i, 0.0, err)
    else:
        ratio_jetpt.SetPoint(i, xd_val, 0.0)
        ratio_jetpt.SetPointError(i, 0.0, 0.0)
# --- Style ratio_jetpt plot ---
pad2_1.cd()
ratio_jetpt.SetTitle("")
ratio_jetpt.GetYaxis().SetTitle("Data/MC")
ratio_jetpt.GetYaxis().SetRangeUser(0.0, 3.0)
ratio_jetpt.GetYaxis().SetNdivisions(505)
ratio_jetpt.GetYaxis().SetTitleSize(0.12)
ratio_jetpt.GetYaxis().SetTitleOffset(0.5)
ratio_jetpt.GetYaxis().SetLabelSize(0.1)

ratio_jetpt.GetXaxis().SetTitle("Offline PF jet p_{T} [GeV]")
ratio_jetpt.GetXaxis().SetTitleSize(0.12)
ratio_jetpt.GetXaxis().SetLabelSize(0.1)
ratio_jetpt.GetXaxis().SetRangeUser(0.0,140.0)

ratio_jetpt.SetMarkerStyle(ROOT.kFullCircle)
ratio_jetpt.SetMarkerSize(0.75)
ratio_jetpt.Draw("APE")

# Draw reference line at ratio_jetpt=1
line = ROOT.TLine(0.0, 1.0, 130.0, 1.0)
line.SetLineStyle(2)
line.Draw("SAME")

can.Update()


can.SaveAs(run+"-trigger_CalojetPt40_Eff.png")
can.SaveAs(run+"-trigger_CalojetPt40_Eff.pdf")
#can.SaveAs("Run3_CalojetPt40_Eff.eps")
#can.SaveAs("Run3_CalojetPt40_Eff.C")


can2 = ROOT.TCanvas("can2", "PromptTrack", 800, 800)
can2.SetTickx()
can2.SetTicky()
can2.Divide(1,2)

# Adjust pads (top for efficiency, bottom for ratio_prompttrk)
pad1_2 = can2.cd(1)
pad1_2.SetPad(0.0, 0.3, 1.0, 1.0)   # xlow, ylow, xup, yup
pad1_2.SetBottomMargin(0.02)        # reduce x-axis space

pad2_2 = can2.cd(2)
pad2_2.SetPad(0.0, 0.0, 1.0, 0.3)
pad2_2.SetTopMargin(0.05)
pad2_2.SetBottomMargin(0.3)

eff_PromptTrk = datafile.Get("eff_PromptTrack")
#eff_PromptTrk_2018 = data18file.Get("eff_PromptTrack")
eff_PromptTrk_mc = mcfile.Get("eff_PromptTrack")
pad1_2.cd()

gra_PromptTrk = eff_PromptTrk.CreateGraph()
#gra_PromptTrk_2018 = eff_PromptTrk_2018.CreateGraph()
gra_PromptTrk_mc = eff_PromptTrk_mc.CreateGraph()
gra_PromptTrk.SetTitle("")
gra_PromptTrk.SetMarkerColor(kBlack)
gra_PromptTrk.SetMarkerStyle(ROOT.kFullCircle)
gra_PromptTrk.SetMarkerSize(0.75)
#gra_PromptTrk_2018.SetMarkerColor(kBlue)
#gra_PromptTrk_2018.SetMarkerStyle(ROOT.kFullCircle)
#gra_PromptTrk_2018.SetMarkerSize(0.75)
#gra_PromptTrk_2018.SetLineColor(kBlue)
gra_PromptTrk_mc.SetMarkerColor(kRed)
gra_PromptTrk_mc.SetMarkerStyle(ROOT.kFullCircle)
gra_PromptTrk_mc.SetMarkerSize(0.75)
gra_PromptTrk_mc.SetLineColor(kRed)
#gra_PromptTrk.GetXaxis().SetTitle("Number of Offline Prompt Tracks")
gra_PromptTrk.GetYaxis().SetRangeUser(0, 1.0)
gra_PromptTrk.GetYaxis().SetTitle("Online Single Jet Tag Efficiency")
gra_PromptTrk.Draw("APE")
#gra_PromptTrk_2018.Draw("PE SAME")
gra_PromptTrk_mc.Draw("PE SAME")

gra_PromptTrk.GetXaxis().SetLabelSize(0)
gra_PromptTrk.GetXaxis().SetTitleSize(0)

text=ROOT.TLatex(0.65, 0.92, "2023 ")
text.SetNDC()
text.SetTextFont(62)
text.SetTextSize(0.05)
text2=ROOT.TLatex(0.15, 0.92, "CMS #bf{#scale[0.75]{#it{Preliminary}}}")
text2.SetNDC()
text2.SetTextSize(0.05)
text2.SetTextFont(62)
text.Draw("SAME")
text2.Draw("SAME")
leg=ROOT.TLegend(0.5,0.70, 0.85,0.85)
leg.SetLineColor(0)
leg.SetFillColorAlpha(0,0)

leg.AddEntry(gra_PromptTrk, run, "pe")
#leg.AddEntry(gra_PromptTrk_2018, "Run2018", "pe")
leg.AddEntry(gra_PromptTrk_mc, "W+Jets selection MC", "pe")
leg.Draw("SAME")
# --- Build ratio_prompttrk graph (data / MC) with errors ---
n_data = gra_PromptTrk.GetN()
n_mc   = gra_PromptTrk_mc.GetN()
npoints = min(n_data, n_mc)

ratio_prompttrk = ROOT.TGraphErrors(npoints)
ratio_prompttrk.SetName("ratio_prompttrk")

for i in range(npoints):
    # Prepare containers for GetPoint
    x_d = ctypes.c_double()
    y_d = ctypes.c_double()
    x_m = ctypes.c_double()
    y_m = ctypes.c_double()

    gra_PromptTrk.GetPoint(i, x_d, y_d)
    gra_PromptTrk_mc.GetPoint(i, x_m, y_m)

    xd_val, yd_val = x_d.value, y_d.value
    xm_val, ym_val = x_m.value, y_m.value

    # Errors
    ey_d = gra_PromptTrk.GetErrorY(i)
    ey_m = gra_PromptTrk_mc.GetErrorY(i)

    if ym_val > 0 and yd_val > 0:
        R = yd_val / ym_val
        err = R * ((ey_d / yd_val)**2 + (ey_m / ym_val)**2)**0.5
        ratio_prompttrk.SetPoint(i, xd_val, R)
        ratio_prompttrk.SetPointError(i, 0.0, err)
    else:
        ratio_prompttrk.SetPoint(i, xd_val, 0.0)
        ratio_prompttrk.SetPointError(i, 0.0, 0.0)
# --- Style ratio_prompttrk plot ---
pad2_2.cd()
ratio_prompttrk.SetTitle("")
ratio_prompttrk.GetYaxis().SetTitle("Data/MC")
ratio_prompttrk.GetYaxis().SetRangeUser(0.0, 3.0)
ratio_prompttrk.GetYaxis().SetNdivisions(505)
ratio_prompttrk.GetYaxis().SetTitleSize(0.12)
ratio_prompttrk.GetYaxis().SetTitleOffset(0.5)
ratio_prompttrk.GetYaxis().SetLabelSize(0.1)

ratio_prompttrk.GetXaxis().SetTitle("Number of Offline Prompt Tracks")
ratio_prompttrk.GetXaxis().SetTitleSize(0.12)
ratio_prompttrk.GetXaxis().SetLabelSize(0.1)
ratio_prompttrk.GetXaxis().SetRangeUser(-0.5,20.0)

ratio_prompttrk.SetMarkerStyle(ROOT.kFullCircle)
ratio_prompttrk.SetMarkerSize(0.75)
ratio_prompttrk.Draw("APE")

# Draw reference line at ratio_prompttrk=1
line = ROOT.TLine(0., 1.0, 20.0, 1.0)
line.SetLineStyle(2)
line.Draw("SAME")

can2.Update()


can2.SaveAs(run+"-trigger_PromptTrack_Eff.png")
can2.SaveAs(run+"-trigger_PromptTrack_Eff.pdf")
#can.SaveAs("Run3_PromptTrack_Eff.eps")
#can.SaveAs("Run3_PromptTrack_Eff.C")

can3 = ROOT.TCanvas("can2", "PromptTrack", 800, 800)
can3.SetTickx()
can3.SetTicky()
can3.Divide(1,2)

# Adjust pads (top for efficiency, bottom for ratio_displacedtrk)
pad1_3 = can3.cd(1)
pad1_3.SetPad(0.0, 0.3, 1.0, 1.0)   # xlow, ylow, xup, yup
pad1_3.SetBottomMargin(0.02)        # reduce x-axis space

pad2_3 = can3.cd(2)
pad2_3.SetPad(0.0, 0.0, 1.0, 0.3)
pad2_3.SetTopMargin(0.05)
pad2_3.SetBottomMargin(0.3)

eff_DisplacedTrk = datafile.Get("eff_DisplacedTrack")
#eff_DisplacedTrk_2018 = data18file.Get("eff_DisplacedTrack")
eff_DisplacedTrk_mc = mcfile.Get("eff_DisplacedTrack")

pad1_3.cd()
gra_DisplacedTrk = eff_DisplacedTrk.CreateGraph()
#gra_DisplacedTrk_2018 = eff_DisplacedTrk_2018.CreateGraph()
gra_DisplacedTrk_mc = eff_DisplacedTrk_mc.CreateGraph()
gra_DisplacedTrk.SetTitle("")
gra_DisplacedTrk.SetMarkerColor(kBlack)
gra_DisplacedTrk.SetMarkerStyle(ROOT.kFullCircle)
gra_DisplacedTrk.SetMarkerSize(0.75)
#gra_DisplacedTrk_2018.SetMarkerColor(kBlue)
#gra_DisplacedTrk_2018.SetMarkerStyle(ROOT.kFullCircle)
#gra_DisplacedTrk_2018.SetMarkerSize(0.75)
#gra_DisplacedTrk_2018.SetLineColor(kBlue)
gra_DisplacedTrk_mc.SetMarkerColor(kRed)
gra_DisplacedTrk_mc.SetMarkerStyle(ROOT.kFullCircle)
gra_DisplacedTrk_mc.SetMarkerSize(0.75)
gra_DisplacedTrk_mc.SetLineColor(kRed)
#gra_DisplacedTrk.GetXaxis().SetTitle("Number of Offline Displaced Tracks")
gra_DisplacedTrk.GetYaxis().SetRangeUser(0, 1.5)
gra_DisplacedTrk.GetYaxis().SetTitle("Online Single Jet Tag Efficiency")
gra_DisplacedTrk.Draw("APE")
#gra_DisplacedTrk_2018.Draw("PE SAME")
gra_DisplacedTrk_mc.Draw("PE SAME")

gra_DisplacedTrk.GetXaxis().SetLabelSize(0)
gra_DisplacedTrk.GetXaxis().SetTitleSize(0)

text=ROOT.TLatex(0.65, 0.92, "2023 ")
text.SetNDC()
text.SetTextFont(62)
text.SetTextSize(0.05)
text2=ROOT.TLatex(0.15, 0.92, "CMS #bf{#scale[0.75]{#it{Preliminary}}}")
text2.SetNDC()
text2.SetTextSize(0.05)
text2.SetTextFont(62)
text.Draw("SAME")
text2.Draw("SAME")
leg=ROOT.TLegend(0.5,0.70, 0.85,0.85)
leg.SetLineColor(0)
leg.SetFillColorAlpha(0,0)

leg.AddEntry(gra_DisplacedTrk, run, "pe")
#leg.AddEntry(gra_DisplacedTrk_2018, "Run2018", "pe")
leg.AddEntry(gra_DisplacedTrk_mc, "W+Jets selection MC", "pe")
leg.Draw("SAME")

# --- Build ratio_displacedtrk graph (data / MC) with errors ---
n_data = gra_DisplacedTrk.GetN()
n_mc   = gra_DisplacedTrk_mc.GetN()
npoints = min(n_data, n_mc)

ratio_displacedtrk = ROOT.TGraphErrors(npoints)
ratio_displacedtrk.SetName("ratio_displacedtrk")

for i in range(npoints):
    # Prepare containers for GetPoint
    x_d = ctypes.c_double()
    y_d = ctypes.c_double()
    x_m = ctypes.c_double()
    y_m = ctypes.c_double()

    gra_DisplacedTrk.GetPoint(i, x_d, y_d)
    gra_DisplacedTrk_mc.GetPoint(i, x_m, y_m)

    xd_val, yd_val = x_d.value, y_d.value
    xm_val, ym_val = x_m.value, y_m.value

    # Errors
    ey_d = gra_DisplacedTrk.GetErrorY(i)
    ey_m = gra_DisplacedTrk_mc.GetErrorY(i)

    if ym_val > 0 and yd_val > 0:
        R = yd_val / ym_val
        err = R * ((ey_d / yd_val)**2 + (ey_m / ym_val)**2)**0.5
        ratio_displacedtrk.SetPoint(i, xd_val, R)
        ratio_displacedtrk.SetPointError(i, 0.0, err)
    else:
        ratio_displacedtrk.SetPoint(i, xd_val, 0.0)
        ratio_displacedtrk.SetPointError(i, 0.0, 0.0)
# --- Style ratio_displacedtrk plot ---
pad2_3.cd()
ratio_displacedtrk.SetTitle("")
ratio_displacedtrk.GetYaxis().SetTitle("Data/MC")
ratio_displacedtrk.GetYaxis().SetRangeUser(0.0, 3.0)
ratio_displacedtrk.GetYaxis().SetNdivisions(505)
ratio_displacedtrk.GetYaxis().SetTitleSize(0.12)
ratio_displacedtrk.GetYaxis().SetTitleOffset(0.5)
ratio_displacedtrk.GetYaxis().SetLabelSize(0.1)

ratio_displacedtrk.GetXaxis().SetTitle("Number of Offline Displaced Tracks")
ratio_displacedtrk.GetXaxis().SetTitleSize(0.12)
ratio_displacedtrk.GetXaxis().SetLabelSize(0.1)
ratio_displacedtrk.GetXaxis().SetRangeUser(-0.5,12.5)

ratio_displacedtrk.SetMarkerStyle(ROOT.kFullCircle)
ratio_displacedtrk.SetMarkerSize(0.75)
ratio_displacedtrk.Draw("APE")

# Draw reference line at ratio_displacedtrk=1
line = ROOT.TLine(0., 1.0, 12.5, 1.0)
line.SetLineStyle(2)
line.Draw("SAME")

can3.Update()

can3.SaveAs(run+"-trigger_DisplacedTrack_Eff.png")
can3.SaveAs(run+"-trigger_DisplacedTrack_Eff.pdf")
#can.SaveAs("Run3_DisplacedTrack_Eff.eps")
#can.SaveAs("Run3_DisplacedTrack_Eff.C")

threshold = 40.0

# --- Open one file for all SFs ---
with open(run+"scale_factors_PtrkShortSig5.txt", "w") as fout:
    fout.write("# Scale factors (Data/MC)\n")
    fout.write("# Each section corresponds to one observable\n\n")

    # --- Jet pt scale factors ---
    fout.write("[Offline PF jet pT (GeV)]\n")
    fout.write("# pT\tRatio(Data/MC)\tStat. error\n")
    for i in range(ratio_jetpt.GetN()):
        x = ctypes.c_double()
        y = ctypes.c_double()
        ratio_jetpt.GetPoint(i, x, y)
        if x.value < threshold: continue
        fout.write(f"{x.value:.1f}\t{y.value:.3f}\t{ratio_jetpt.GetEY()[i]:.3f}\n")
    fout.write("\n")

    # --- Prompt track scale factors ---
    fout.write("[Number of Offline Prompt Tracks]\n")
    fout.write("# Ntrk\tRatio(Data/MC)\tStat. error\n")
    for i in range(ratio_prompttrk.GetN()):
        x = ctypes.c_double()
        y = ctypes.c_double()
        ratio_prompttrk.GetPoint(i, x, y)
        fout.write(f"{x.value:.1f}\t{y.value:.3f}\t{ratio_prompttrk.GetEY()[i]:.3f}\n")
    fout.write("\n")
