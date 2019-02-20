/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
//#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnalysisValidation.h"
#include "TopNtupleAnalysis/AnaTtresCheckTrigger.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"
#include "TopNtupleAnalysis/UserFunktions.h"


AnaTtresCheckTrigger::AnaTtresCheckTrigger(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
  : AnalysisValidation(filename, systList), m_electron(electron), m_boosted(boosted),
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

    int varBinN10  = 100;
    double varBin10[101];
    for(int j=0;j<101;j++) varBin10[j] = -1 + 0.02*j; 

    //Plots for QCD validation
    if (m_electron){
      if(m_boosted){
        double varBin6[6] = {30, 40, 60, 120, 400, 700};
        double varBin7[5]  = {0., 0.4, 0.6, 1.0, 1.5}; 
        int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
        int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;	

      }else{
        double varBin6[8] = {30, 35, 40, 50, 60, 120, 400, 700};
        double varBin7[7]  = {0., 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
        int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
        int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;	

      }//m_boosted

    }else{
      if(m_boosted){
        double varBin6[7] = {25, 30, 40, 50, 100, 400, 700};
        double varBin7[6] = {0., 0.4, 0.6, 0.8, 1.0, 1.5};
        int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
        int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;	

      }else{
        double varBin6[9] = {25, 30, 35, 40, 50, 70, 100, 400, 700};
        double varBin7[7] = {0., 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
        int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
        int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;


      }//m_boosted
    }//m_electron

    double varBin8[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500,3000};
    int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;
    double varBin9[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000};
    int varBinN9 = sizeof(varBin9)/sizeof(double) - 1;

    m_hSvc.create1D("mu", "; #mu; Events", 50,0,50);

    //check trigger
    m_hSvc.create1D("wpt", "; p_{T,W} [GeV]; Events", 100, 0, 1000);
    m_hSvc.create1D("wptPassMet110", "; p_{T,W} [GeV]; Events", 100, 0, 1000);
    m_hSvc.create1D("wptPassPufitMet110", "; p_{T,W} [GeV]; Events", 100, 0, 1000);

    m_hSvc.create1D("diMuonMass", "; Di-muon Mass [GeV]; Events", 40, 0, 200);
    m_hSvc.create1D("muonPtBarrel", "; muon p_{T} [GeV]; Events", 50,0,100);  
    m_hSvc.create1D("muonPtEndcap", "; muon p_{T} [GeV]; Events", 50,0,100);  
    m_hSvc.create1D("muonPtBarrelPassTrigger", "; muon p_{T} [GeV]; Events", 50,0,100);  
    m_hSvc.create1D("muonPtEndcapPassTrigger", "; muon p_{T} [GeV]; Events", 50,0,100);  
    m_hSvc.create1D("muonPhiBarrel", "; muon #phi; Events", 32,-3.2,3.2);  
    m_hSvc.create1D("muonPhiEndcap", "; muon #phi; Events", 32,-3.2,3.2);  
    m_hSvc.create1D("muonPhiBarrelPassTrigger", "; muon #phi; Events", 32,-3.2,3.2);  
    m_hSvc.create1D("muonPhiEndcapPassTrigger", "; muon #phi; Events", 32,-3.2,3.2);  
    m_hSvc.create1D("muonEta", "; muon #eta; Events", 50,-2.5,2.5);  
    m_hSvc.create1D("muonEtaPassTrigger", "; muon #eta; Events", 50,-2.5,2.5);  

    m_hSvc.create1D("diElectronMass", "; Di-electron Mass [GeV]; Events", 40, 0, 200);
    m_hSvc.create1D("electronPt", "; electron p_{T} [GeV]; Events", 50,0,100);  
    m_hSvc.create1D("electronPtPassTrigger", "; electron p_{T} [GeV]; Events", 50,0,100);  
    m_hSvc.create1D("electronEta", "; electron #eta; Events", 50,-2.5,2.5);  
    m_hSvc.create1D("electronEtaPassTrigger", "; electron #eta; Events", 50,-2.5,2.5);  
    m_hSvc.create1D("electronPhi", "; electron #phi; Events", 32,-3.2,3.2);  
    m_hSvc.create1D("electronPhiPassTrigger", "; electron #phi; Events", 32,-3.2,3.2);  
  }

AnaTtresCheckTrigger::~AnaTtresCheckTrigger() {
}

void AnaTtresCheckTrigger::run(const Event &evt, double weight, const std::string &s, int mode, TSpline* m_spline) {

  HistogramService *h = &m_hSvc;
  bool isTight = false;
  std::string suffix = s;

  //MET trigger
  bool pass_met_trigger = evt.trigger("HLT_xe110_mht_L1XE50"); 
  bool pass_pufit_met_trigger = evt.trigger("HLT_xe110_pufit_L1XE50"); 
  //std::cout << "pass met trigger=" << pass_met_trigger << std::endl;
  //std::cout << "pass pufit met trigger=" << pass_pufit_met_trigger << std::endl;
  //Muon trigger
  bool pass_muon_trigger = evt.trigger("HLT_mu26_ivarmedium") || evt.trigger("HLT_mu50"); 
  //std::cout << "pass muon trigger=" << pass_muon_trigger << std::endl;
  //Electron trigger
  bool pass_electron_trigger = evt.trigger("HLT_e26_lhtight_nod0_ivarloose") || evt.trigger("HLT_e60_lhmedium_nod0") || evt.trigger("HLT_e140_lhloose_nod0"); 
  //std::cout << "pass electron trigger=" << pass_electron_trigger << std::endl;

  //Number of muon
  int nMuon30 = 0, nMuon25=0, muon_id=-1;
  for(int im=0;im<evt.muon().size();im++){
    TLorentzVector mu = evt.muon()[im].mom();
    //bool isTight = evt.muon()[im].isTight();
    //std::cout << "muon pt/eta/phi=" << mu.Perp() << "/" << mu.Eta() << "/" << mu.Phi() << std::endl;
    //if(!isTight) continue;
    if(mu.Perp()*1e-3 > 25) nMuon25++;
    if(mu.Perp()*1e-3 > 30) {
      if(nMuon30==0) muon_id=im;
      nMuon30++;
    }
  }

  //Number of electron
  int nElectron30=0,nElectron25=0,electron_id=-1;
  for(int ie=0;ie<evt.electron().size();ie++){
    TLorentzVector el = evt.electron()[ie].mom();
    //std::cout << "electron pt/eta/phi=" << el.Perp() << "/" << el.Eta() << "/" << el.Phi() << std::endl;
    //bool isTight = evt.muon()[ie].isTight();
    //if(!isTight) continue;
    if(el.Perp()*1e-3 > 25) nElectron25++;
    if(el.Perp()*1e-3 > 30) {
      if(nElectron30==0) electron_id=ie;
      nElectron30++;
    }
  }

  bool pass_muon_selection = evt.passes("check_muon_trigger");
  bool pass_electron_selection = evt.passes("check_electron_trigger");
  bool pass_met_selection = evt.passes("check_met_trigger");
  //std::cout << "pass selection muon/electron/met=" << pass_muon_selection << "/" << pass_electron_selection << "/" << pass_met_selection << std::endl;


  if(pass_met_selection){
    TLorentzVector lead_mu = evt.muon()[muon_id].mom();
    TLorentzVector vmet = evt.met();
    float wpt = (vmet+lead_mu).Perp()*1e-3;
    h->h1D("wpt", "", suffix)->Fill(wpt, weight);
    if(pass_met_trigger) h->h1D("wptPassMet110", "", suffix)->Fill(wpt, weight);
    if(pass_pufit_met_trigger) h->h1D("wptPassPufitMet110", "", suffix)->Fill(wpt, weight);
  }
  
  if(pass_muon_selection){
    TLorentzVector lead_mu = evt.muon()[muon_id].mom();
    for(int im=0;im<evt.muon().size();im++){
      if(im==muon_id) continue;
      TLorentzVector mu = evt.muon()[im].mom();
      //bool isTight = evt.muon()[im].isTight();
      float pt  = mu.Perp()*1e-3;
      float eta = mu.Eta();
      float phi = mu.Phi();
      float mass = (lead_mu+mu).M()*1e-3;
      bool pass = evt.muon()[im].HLT_mu26_ivarmedium() || evt.muon()[im].HLT_mu50();
      //std::cout << "muon pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl;
      //std::cout << "isTight=" << isTight << std::endl;
      //std::cout << "di-muon mass=" << mass << std::endl; 
      h->h1D("diMuonMass", "", suffix)->Fill(mass, weight);
      if(mass>80 && mass<100){
        if(pt>30)h->h1D("muonEta", "", suffix)->Fill(eta, weight);
        if(pt>30 && pass)h->h1D("muonEtaPassTrigger", "", suffix)->Fill(eta, weight);
        if(fabs(eta)<1.05) {//Barrel
          h->h1D("muonPtBarrel", "", suffix)->Fill(pt, weight);
          if(pt>30) h->h1D("muonPhiBarrel", "", suffix)->Fill(phi, weight);
          if(pass){
            h->h1D("muonPtBarrelPassTrigger", "", suffix)->Fill(pt, weight);
            if(pt>30)h->h1D("muonPhiBarrelPassTrigger", "", suffix)->Fill(phi, weight);
          }
        } else {//Endcap
          h->h1D("muonPtEndcap", "", suffix)->Fill(pt, weight);
          if(pt>30) h->h1D("muonPhiEndcap", "", suffix)->Fill(phi, weight);
          if(pass){
            h->h1D("muonPtEndcapPassTrigger", "", suffix)->Fill(pt, weight);
            if(pt>30)h->h1D("muonPhiEndcapPassTrigger", "", suffix)->Fill(phi, weight);
          }
        }
      }
    }
  }

  if(pass_electron_selection){
    TLorentzVector lead_mu = evt.electron()[electron_id].mom();
    for(int ie=0;ie<evt.electron().size();ie++){
      if(ie==electron_id) continue;
      TLorentzVector mu = evt.electron()[ie].mom();
      //bool isTight = evt.electron()[ie].isTight();
      float pt  = mu.Perp()*1e-3;
      float eta = mu.Eta();
      float phi = mu.Phi();
      float mass = (lead_mu+mu).M()*1e-3;
      bool pass = evt.electron()[ie].HLT_e26_lhtight_nod0_ivarloose() || evt.electron()[ie].HLT_e60_lhmedium_nod0() || 
        evt.electron()[ie].HLT_e140_lhloose_nod0();
      //std::cout << "electron pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl;
      //std::cout << "isTight=" << isTight << std::endl;
      //std::cout << "di-electron mass=" << mass << std::endl; 
      h->h1D("diElectronMass", "", suffix)->Fill(mass, weight);
      if(mass>80 && mass<100){
        h->h1D("electronPt", "", suffix)->Fill(pt, weight);
        if(pt>30) h->h1D("electronPhi", "", suffix)->Fill(phi, weight);
        if(pt>30)h->h1D("electronEta", "", suffix)->Fill(eta, weight);
        if(pt>30 && pass)h->h1D("electronEtaPassTrigger", "", suffix)->Fill(eta, weight);
        if(pass){
          h->h1D("electronPtPassTrigger", "", suffix)->Fill(pt, weight);
          if(pt>30)h->h1D("electronPhiPassTrigger", "", suffix)->Fill(phi, weight);
        }
      }
    }
  }

}
