/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresHtaumu.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"


AnaTtresHtaumu::AnaTtresHtaumu(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
  : Analysis(filename, systList), m_electron(electron), m_boosted(boosted),
    m_neutrinoBuilder("MeV"), m_chi2("MeV") {

  m_chi2.Init(TtresChi2::DATA2015_MC15C);

  m_hSvc.m_tree->Branch("truemtt",    &_tree_truemtt);
  m_hSvc.m_tree->Branch("mtt",    &_tree_mtt);
  m_hSvc.m_tree->Branch("weight", &_tree_weight);
  m_hSvc.m_tree->Branch("cat",    &_tree_cat);
  m_hSvc.m_tree->Branch("syst",   &_tree_syst);


  double varBin1[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 340, 380, 450, 500};
  int varBinN1 = sizeof(varBin1)/sizeof(double) - 1;
  double varBin2[] = {300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 540, 580, 620, 660, 700, 800, 1e3, 1.2e3, 1.5e3};
  int varBinN2 = sizeof(varBin2)/sizeof(double) - 1;
  double varBin3[] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500};
  int varBinN3 = sizeof(varBin3)/sizeof(double) - 1;
  double varBin4[] = {80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500};
  int varBinN4 = sizeof(varBin4)/sizeof(double) - 1;
  double varBin5[] = {0, 80, 160, 240, 320, 400, 480, 560,640,720,800,920,1040,1160,1280,1400,1550,1700,2000,2300,2600,2900,3200,3600,4100,4600,5100,6000};
  int varBinN5 = sizeof(varBin5)/sizeof(double) - 1;
  double varBin6[] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 540, 580, 620, 660, 700, 800, 1e3, 1.2e3, 1.5e3};
  int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;

  //Plots for QCD validation
  if (m_electron){
     if(m_boosted){
        double varBin6[6] = {30, 40, 60, 120, 400, 700};
	double varBin7[5]  = {0., 0.4, 0.6, 1.0, 1.5}; 
	int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
	int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;	
	m_hSvc.create1DVar("lepPt_effBins", "; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
	m_hSvc.create1DVar("closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);
	
     }else{
        double varBin6[8] = {30, 35, 40, 50, 60, 120, 400, 700};
	double varBin7[7]  = {0., 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
	int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
	int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;	
	m_hSvc.create1DVar("lepPt_effBins", "; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
	m_hSvc.create1DVar("closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);
	
     }//m_boosted

  }else{
     if(m_boosted){
        double varBin6[7] = {25, 30, 40, 50, 100, 400, 700};
	double varBin7[6] = {0., 0.4, 0.6, 0.8, 1.0, 1.5};
	int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
	int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;	
	m_hSvc.create1DVar("lepPt_effBins", "; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
	m_hSvc.create1DVar("closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);

     }else{
        double varBin6[9] = {25, 30, 35, 40, 50, 70, 100, 400, 700};
	double varBin7[7] = {0., 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
	int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
	int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;
	m_hSvc.create1DVar("lepPt_effBins", "; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
	m_hSvc.create1DVar("closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);
	
     }//m_boosted
  }//m_electron

  m_hSvc.create1D("higgsMass", "; mass  [GeV]; Events", 20, 0, 200);
  m_hSvc.create1D("MET", "; met  [GeV]; Events", 30, 0, 150);
  m_hSvc.create1D("METPhi", "; #phi of met  ; Events", 32, -3.2, 3.2);
  m_hSvc.create1D("muonPt", "; muon p_{T}  [GeV]; Events", 30, 0, 150);
  m_hSvc.create1D("muonEta", "; muon #eta; Events", 25, -2.5, 2.5);
  m_hSvc.create1D("muonPhi", "; muon #phi; Events", 32, -3.2, 3.2);
  m_hSvc.create1D("electronPt", "; electron pt  [GeV]; Events", 30, 0, 150);
  m_hSvc.create1D("electronEta", "; electron #eta; Events", 25, -2.5, 2.5);
  m_hSvc.create1D("electronPhi", "; electron #phi; Events", 32, -3.2, 3.2);
}

AnaTtresHtaumu::~AnaTtresHtaumu() {
}

void AnaTtresHtaumu::run(const Event &evt, double weight, const std::string &s, int is2016run) {
  std::cout << "AnaTtresHtaumu::run" << std::endl;
  //std::cout << "bmujets/met110=" << evt.passes("bmujets_2016") << "/" << evt.passes("bmujetsmet110") << std::endl;
  // check channel
  //
  char name[200];
  std::string suffix = s;
  HistogramService *h = &m_hSvc;

  for(int im=0;im<evt.muon().size();im++){
    TLorentzVector mu = evt.muon()[im].mom();
    float pt = mu.Perp()*1e-3;
    float eta = mu.Eta();
    float phi = mu.Phi();
    std::cout << "muon pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl;
  }

  for(int ie=0;ie<evt.electron().size();ie++){
    TLorentzVector el = evt.electron()[ie].mom();
    float pt = el.Perp()*1e-3;
    float eta = el.Eta();
    float phi = el.Phi();
    std::cout << "electron pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl;
  }
  
  if(evt.muon().size() && evt.electron().size()){
    TLorentzVector mu = evt.muon()[0].mom();
    TLorentzVector el = evt.electron()[0].mom();
    float met = evt.met().Perp();
    float metphi = evt.met().Phi();
    TLorentzVector nu;
    nu.SetPtEtaPhiM(met,el.Eta(),el.Phi(),0);
    float mass = (mu+el+nu).M();
    std::cout << "mass=" << mass << std::endl;
    h->h1D("higgsMass", "", suffix)->Fill(mass*1e-3);
    h->h1D("MET", "", suffix)->Fill(met*1e-3);
    h->h1D("METPhi", "", suffix)->Fill(metphi);
    h->h1D("muonPt", "", suffix)->Fill(mu.Perp()*1e-3);
    h->h1D("muonEta", "", suffix)->Fill(mu.Eta());
    h->h1D("muonPhi", "", suffix)->Fill(mu.Phi());
    h->h1D("electronPt", "", suffix)->Fill(el.Perp()*1e-3);
    h->h1D("electronEta", "", suffix)->Fill(el.Eta());
    h->h1D("electronPhi", "", suffix)->Fill(el.Phi());
  }

}    
