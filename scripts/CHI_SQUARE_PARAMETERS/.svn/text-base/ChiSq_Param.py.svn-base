#--- Takes as input the output from './read ' of AnaTtresSL.cxx -------#
from ROOT import *
gStyle.SetPalette(1)
gStyle.SetPadColor(0)
gStyle.SetCanvasColor(0)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleOffset(.45)
gStyle.SetTitleFont(42)
gStyle.SetPaintTextFormat("5.1f")

import math

def GetFunction(hist, low, high) :

    sig = 2.

    hist.GetXaxis().SetRangeUser(low,high)

    f = TF1("func", "gaus", low, high)

    hist.Fit("func","","", low, high)
    f = hist.GetFunction("func")
    #f.SetRange(f.GetParameter(0) - sig*f.GetParameter(1), f.GetParameter(0) + sig*f.GetParameter(1))

    f1 = TF1("func1", "gaus", f.GetParameter(0) - sig*f.GetParameter(1), f.GetParameter(0) + sig*f.GetParameter(1))
    hist.Fit("func1","","", f.GetParameter(0) - sig*f.GetParameter(1), f.GetParameter(0) + sig*f.GetParameter(1))
    f1 = hist.GetFunction("func1")
    #f1.SetRange(f1.GetParameter(0) - sig*f1.GetParameter(1), f1.GetParameter(0) + sig*f1.GetParameter(1))

    f2 = TF1("func2", "gaus", f1.GetParameter(0) - sig*f1.GetParameter(1), f1.GetParameter(0) + sig*f1.GetParameter(1))
    hist.Fit("func2","","", f1.GetParameter(0) - sig*f1.GetParameter(1), f1.GetParameter(0) + sig*f1.GetParameter(1))
    f2 = hist.GetFunction("func2")
    #f2.SetRange(f2.GetParameter(0) - sig*f2.GetParameter(1), f2.GetParameter(0) + sig*f2.GetParameter(1))

    f2.SetRange(low, high)
    f2.SetLineWidth(2)

    return f2

## end of GetFunction(hist) :

def Draw(h, f) :
    #h.Scale(1./h.Integral())
    f1 = f.Clone('NewFunc')
    h.SetLineColor(2)
    h.SetMarkerColor(2)
    h.SetMarkerStyle(20)
    name = h.GetName()
    c = TCanvas('canv_'+name, 'canv_'+name, 600, 600)

    c.cd()
    h.Draw()
    f.SetLineColor(1)
    f.SetLineStyle(7)
    f.Draw('same')
    f1.SetLineColor(1)
    f1.SetLineWidth(3)
    f1.SetRange(f.GetParameter(1) - f.GetParameter(2), f.GetParameter(1) + f.GetParameter(2))
    f1.Draw('same')

    pave = TPaveText(0.61,0.80,0.89,0.90,"NDC");
    mu = round(f.GetParameter(1),2)
    sigma = round(f.GetParameter(2),2)
    n_mu = '#mu = %f' % mu
    pave.AddText(n_mu)
    n_sig = '#sigma = %f' % sigma
    pave.AddText(n_sig)
    pave.SetFillColor(0);
    pave.SetLineColor(0);
    pave.SetBorderSize(0);
    pave.SetFillStyle(0);
    pave.SetTextFont(62);
    pave.SetTextSize(0.04);
    pave.SetTextAlign(12);
    pave.Draw("same");

    c.SaveAs('canv_'+name+'.png')

## end of Draw(h, f) :

f_e  = TFile.Open("resolved_e_o.root","READ")
f_mu = TFile.Open("resolved_mu_o.root","READ")

h_WH_e  = f_e.Get('W_Hadronic')
h_WH_mu = f_mu.Get('W_Hadronic')

h_TH_e  = f_e.Get('T_Hadronic')
h_TH_mu = f_mu.Get('T_Hadronic')

h_TL_e  = f_e.Get('T_Leptonic')
h_TL_mu = f_mu.Get('T_Leptonic')

h_PT_e  = f_e.Get('PT_Diff')
h_PT_mu = f_mu.Get('PT_Diff')

h_WH_e.Add(h_WH_mu)
h_TH_e.Add(h_TH_mu)
h_TL_e.Add(h_TL_mu)
h_PT_e.Add(h_PT_mu)

h_WH_e.GetXaxis().SetTitle('M_{jj} [GeV]')
h_TH_e.GetXaxis().SetTitle('M_{jjb} - M_{jj} [GeV]')
h_TL_e.GetXaxis().SetTitle('M_{jl#nu} [GeV]')

h_PT_e.Rebin(10)
#h_WH_e.Rebin(10)
f_WH = GetFunction(h_WH_e, 50, 110)
f_TH = GetFunction(h_TH_e, 60, 110)
f_TL = GetFunction(h_TL_e, 100, 250)
f_PT = GetFunction(h_PT_e, -100, 100)


Draw(h_TL_e, f_TL)
Draw(h_TH_e, f_TH)
Draw(h_WH_e, f_WH)
Draw(h_PT_e, f_PT)
#print "MaxBin = %i" % (maxBin)
