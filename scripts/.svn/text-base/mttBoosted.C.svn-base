
#include "TopEventReconstructionTools/Root/TtresNeutrinoBuilder.cxx"
#include "TopEventReconstructionTools/TopEventReconstructionTools/TtresNeutrinoBuilder.h"

#include "TLorentzVector.h"
#include "TChain.h"
#include "TH1F.h"

#include "nominal.h"
#include "nominal.C"

#include <vector>
#include <map>
#include <string>
#include <sstream>

#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

map<string, TH1F *> h[2];

void addH(int idx, const string &n, const string &t, int nbins, float a, float b) {
  std::stringstream ss;
  ss << n << "_" << idx;
  h[idx][n.c_str()] = new TH1F(ss.str().c_str(), t.c_str(), nbins, a, b);
}

void draw(const string &s) {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  TCanvas *canv = new TCanvas("canv", "", 800, 600);
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(h[0][s], "Standard");
  leg->AddEntry(h[1][s], "VR");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  h[0][s]->SetLineColor(kBlue);
  h[1][s]->SetLineColor(kRed);
  h[0][s]->Draw("hist");
  h[1][s]->Draw("hist same");
  leg->Draw();
  stringstream ss;
  ss << "mttBoosted_" << s << ".eps";
  canv->SaveAs(ss.str().c_str());
  delete canv;
  delete leg;
}

void mttBoosted() {

  TtresNeutrinoBuilder m_neutrinoBuilder("MeV");

  TChain c1("nominal");
  c1.Add("output_e4.root");

  TChain c2("nominal");
  c2.Add("output_vr.root");

  TChain *clist[] = {&c1, &c2};

  addH(0, "mtt", "; m_{tt} ; Entries", 70, 0, 7000);
  addH(1, "mtt", "; m_{tt} ; Entries", 70, 0, 7000);

  for (int nc = 0; nc < 2; ++nc) {
    TChain *c = clist[nc];

    nominal nom(c);
    for (int z = 0; z < c->GetEntries(); ++z) {
      c->GetEntry(z);

      if (!(nom.bejets || nom.bmujets)) continue;

      int l = -1;
      for (int k = 0; k < nom.ljet_good->size(); ++k) {
        if (nom.ljet_good->at(k) == 1) {
          l = k;
          break;
        }
      }

      int s = -1;
      for (int k = 0; k < nom.jet_closeToLepton->size(); ++k) {
        if (nom.jet_closeToLepton->at(k) == 1) {
          s = k;
          break;
        }
      }


      if (l == -1) continue;
      if (s == -1) continue;

      TLorentzVector ljet(0,0,0,0);
      ljet.SetPtEtaPhiE(nom.ljet_pt->at(l), nom.ljet_eta->at(l), nom.ljet_phi->at(l), nom.ljet_e->at(l));

      TLorentzVector sjet(0,0,0,0);
      sjet.SetPtEtaPhiE(nom.jet_pt->at(s), nom.jet_eta->at(s), nom.jet_phi->at(s), nom.jet_e->at(s));

      bool electron = false;
      TLorentzVector lepton(0,0,0,0);
      if (nom.el_pt->size() == 0 && nom.mu_pt->size() == 1) {
        electron = false;
        lepton.SetPtEtaPhiE(nom.mu_pt->at(0), nom.mu_eta->at(0), nom.mu_phi->at(0), nom.mu_e->at(0));
      } else if (nom.el_pt->size() == 1 && nom.mu_pt->size() == 0) {
        electron = true;
        lepton.SetPtEtaPhiE(nom.el_pt->at(0), nom.el_eta->at(0), nom.el_phi->at(0), nom.el_e->at(0));
      } else continue;

      
      TLorentzVector tt(0,0,0,0);
      tt += ljet;
      tt += lepton;
      tt += sjet;
      std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&lepton, nom.met_met, nom.met_phi, true);   
      TLorentzVector nu(0,0,0,0);
      if (vec_nu.size() > 0) {
        nu = *(vec_nu[0]);
        for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
        vec_nu.clear();
      }
      tt += nu;

      h[nc]["mtt"]->Fill(tt.M()*1e-3, nom.mcWeight);
    } // all events in this file
  } // for each mode

  draw("mtt");
}

