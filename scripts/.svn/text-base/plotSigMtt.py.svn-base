#!/usr/bin/env python
from os import path

from ROOT import *

Internal = True
Internal = False

gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

cacc = TCanvas("cacc", "", 800, 600);
cacc.SetBottomMargin(0.2)
cacc.SetLeftMargin(0.14)
l = TLegend(0.6,0.50,0.87,0.75)
l.SetBorderSize(0)

cn10 = TChain("mini")
cn20 = TChain("mini")
cn25 = TChain("mini")
cn30 = TChain("mini")

hn_zp10 = TH1F("hn_zp10", "", 36, 0.0, 3.6)
hn_zp20 = TH1F("hn_zp20", "", 36, 0.0, 3.6)
hn_zp25 = TH1F("hn_zp25", "", 36, 0.0, 3.6)
hn_zp30 = TH1F("hn_zp30", "", 36, 0.0, 3.6)

def setStyle(h, col, sty):
    h.GetYaxis().SetRangeUser(0, 0.6);
    h.GetYaxis().SetTitle("Arbitrary Units")
    h.GetXaxis().SetTitle("m_{t#bar{t}}^{reco} [TeV]");
    h.GetXaxis().SetTitleOffset(1.0);
    h.GetYaxis().SetTitleOffset(1.0);
    h.GetXaxis().SetLabelSize(0.05);
    h.GetXaxis().SetTitleSize(0.05);
    h.GetYaxis().SetLabelSize(0.05);
    h.GetYaxis().SetTitleSize(0.05);
    h.SetLineWidth(3);
    h.SetLineStyle(sty);
    h.SetLineColor(col);

setStyle(hn_zp10, kBlue,  4)
setStyle(hn_zp20, kGreen+2,  5)
setStyle(hn_zp25, kPink+2,  5)
setStyle(hn_zp30, kBlack,  2)

direc = "ttres_mc15b_splitb_bugfix2"
cn10.Add(direc+"/boosted_*Zprime1000.root")
cn20.Add(direc+"/boosted_*Zprime2000.root")
cn25.Add(direc+"/boosted_*Zprime2500.root")
cn30.Add(direc+"/boosted_*Zprime3000.root")

massList = {}
massList["zp10"] = [cn10, hn_zp10]
massList["zp20"] = [cn20, hn_zp20]
massList["zp25"] = [cn25, hn_zp25]
massList["zp30"] = [cn30, hn_zp30]

for m in massList:
  cn = massList[m][0]
  hn = massList[m][1]
  for i in xrange(cn.GetEntries()):
    cn.GetEntry(i)
    if cn.syst != '':
        continue
    hn.Fill(cn.mtt*1e-3, cn.weight)

for hist in [hn_zp10, hn_zp20, hn_zp25, hn_zp30]:
    s = hist.Integral(-1, 9999)
    hist.Scale(1.0/s)

hn_zp10.GetYaxis().SetRangeUser(0, 0.6)
hn_zp10.Draw("hist")
hn_zp20.Draw("hist same")
hn_zp25.Draw("hist same")
hn_zp30.Draw("hist same")

l.AddEntry(hn_zp10, "m(Z')=1.0TeV", "L")
l.AddEntry(hn_zp20, "m(Z')=2.0TeV", "L")
l.AddEntry(hn_zp25, "m(Z')=2.5TeV", "L")
l.AddEntry(hn_zp30, "m(Z')=3.0TeV", "L")
l.Draw()

gPad.RedrawAxis()

def stampATLAS(text, x, y):
  t = TLatex()
  t.SetNDC()
  t.SetTextFont(72)
  t.SetTextColor(1)
  delx = 0.115*696*gPad.GetWh()/(472*gPad.GetWw()) + 0.05
  t.SetTextSize(0.06)
  t.DrawLatex(x,y,"ATLAS")
  t.SetTextFont(42)
  t.SetTextSize(0.06)
  t.DrawLatex(x+delx, y, text)

def stampLumiText(lumi, x, y, text, size):
  t = TLatex()
  t.SetNDC()
  t.SetTextFont(42)
  t.SetTextColor(1)
  t.SetTextSize(size)
  #t.DrawLatex(x,y,"#int L dt = "+str(lumi)+" fb^{-1}, "+text)
  t.DrawLatex(x,y, text+", "+str(lumi)+" fb^{-1}")

def stampText(x, y, text, size):
  t = TLatex()
  t.SetNDC()
  t.SetTextFont(42)
  t.SetTextColor(1)
  t.SetTextSize(size)
  t.DrawLatex(x,y, text)

if Internal:
    stampATLAS("Simulation Internal", 0.18, 0.83)
else:
    stampATLAS("Simulation", 0.18, 0.83)
stampText(0.18, 0.75, "#sqrt{s} = 13 TeV", 0.05)

suf = ""
if not Internal:
    suf = "_ATLASSimul"
for i in [".eps", ".png", ".pdf"]:
    cacc.SaveAs("signalmtt"+suf+i)

  
