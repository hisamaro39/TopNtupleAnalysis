/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
//#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnalysisValidation.h"
#include "TopNtupleAnalysis/AnaTtresMetBtag.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"


AnaTtresMetBtag::AnaTtresMetBtag(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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
        //m_hSvc.create1DVar("lepPt_effBins", "; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
        //m_hSvc.create1DVar("closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);

      }else{
        double varBin6[8] = {30, 35, 40, 50, 60, 120, 400, 700};
        double varBin7[7]  = {0., 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
        int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
        int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;	
        //m_hSvc.create1DVar("lepPt_effBins", "; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
        //m_hSvc.create1DVar("closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);

      }//m_boosted

    }else{
      if(m_boosted){
        double varBin6[7] = {25, 30, 40, 50, 100, 400, 700};
        double varBin7[6] = {0., 0.4, 0.6, 0.8, 1.0, 1.5};
        int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
        int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;	
        //m_hSvc.create1DVar("lepPt_effBins", "; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
        //m_hSvc.create1DVar("closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);

      }else{
        double varBin6[9] = {25, 30, 35, 40, 50, 70, 100, 400, 700};
        double varBin7[7] = {0., 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
        int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
        int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;
        //m_hSvc.create1DVar("lepPt_effBins", "; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
        //m_hSvc.create1DVar("closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);


      }//m_boosted
    }//m_electron

    double varBin8[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500,3000};
    int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;
    double varBin9[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000};
    int varBinN9 = sizeof(varBin9)/sizeof(double) - 1;

    //m_hSvc.create1D("yields", "; One ; Events", 1, 0.5, 1.5);

    m_hSvc.create1D("lepPt", "; lepton p_{T} [GeV]; Events", 50,0,1000);  
    m_hSvc.create1D("lepEta", "; lepton #eta; Events", 25,-2.5,2.5);  
    m_hSvc.create1D("lepPhi", "; lepton #phi; Events", 32,-3.2,3.2);  
    m_hSvc.create1D("nLepton", "; # of lepton; Events", 5,0,5);  

    m_hSvc.create1D("smallJetPt", "; small jet p_{T} [GeV]; Events", 50, 0, 1000);
    m_hSvc.create1D("smallJetEta", "; small jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("smallJetPhi", "; small jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("smallJetMass", "; small jet mass [GeV] ; Events", 30, 0, 300);
    m_hSvc.create1D("leadSmallJetPt", "; leading small jet p_{T} [GeV]; Events", 50, 0, 2000);
    m_hSvc.create1D("leadSmallJetEta", "; leading small jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("leadSmallJetPhi", "; leading small jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("leadSmallJetMass", "; leading small jet mass [GeV] ; Events", 30, 0, 300);
    m_hSvc.create1D("secondLeadSmallJetPt", "; 2nd leading small jet p_{T} [GeV]; Events", 50, 0, 1000);
    m_hSvc.create1D("secondLeadSmallJetEta", "; 2nd leading small jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("secondLeadSmallJetPhi", "; 2nd leading small jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("secondLeadSmallJetMass", "; 2nd leading small jet mass [GeV] ; Events", 30, 0, 300);
    m_hSvc.create1D("nSmallJet", "; # of small jets ; Events", 15, 0, 15);

    m_hSvc.create1D("largeJetPt", "; large jet p_{T} [GeV] ; Events", 50, 0, 2000);
    m_hSvc.create1D("largeJetEta", "; large jet #eta ; Events", 20, -2., 2.);
    m_hSvc.create1D("largeJetPhi", "; large jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("largeJetMass", "; large jet mass [GeV] ; Events", 30, 0, 300);
    m_hSvc.create1D("largeJetTau32", "; large jet tau32  ; Events", 20, 0, 1);
    m_hSvc.create1D("nLargeJet", "; # of large jets ; Events", 10, 0, 10);

    m_hSvc.create1D("trackJetPt", "; track jet p_{T} [GeV]; Events", 25, 0, 500);
    m_hSvc.create1D("trackJetEta", "; track jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("trackJetPhi", "; track jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("trackJetMv2c10", "; track jet Mv2c10 ; Events", 20, -1, 1);
    m_hSvc.create1D("bmatchedTrackJetMv2c10", "; track jet Mv2c10 ; Events", 20, -1, 1);
    m_hSvc.create1D("notBmatchedTrackJetMv2c10", "; track jet Mv2c10 ; Events", 20, -1, 1);
    m_hSvc.create1D("nTrackJet", "; # of track jets ; Events", 15, 0, 15);  
    m_hSvc.create1D("nBtaggedTrackJet", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nBtaggedBmatchedTrackJet", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nBtaggedBmatchedTrackJetLeptonic", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nBtaggedBmatchedTrackJetHadronic", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nBtaggedBmatchedTrackJet2b", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nBtaggedBmatchedTrackJet1bLeptonic", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nBtaggedNotBmatchedTrackJet", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nFlatBtaggedTrackJet", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nFlatBtaggedBmatchedTrackJet", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nFlatBtaggedBmatchedTrackJetLeptonic", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nFlatBtaggedBmatchedTrackJetHadronic", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nFlatBtaggedBmatchedTrackJet2b", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nFlatBtaggedBmatchedTrackJet1bLeptonic", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nFlatBtaggedNotBmatchedTrackJet", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nBmatchedTrackJet", "; # of track jets ; Events", 5, 0, 5);  
    m_hSvc.create1D("nBmatchedTrackJetLeptonic", "; # of track jets ; Events", 5, 0, 5);  
    m_hSvc.create1D("nBmatchedTrackJetHadronic", "; # of track jets ; Events", 5, 0, 5);  
    m_hSvc.create1D("nNotBmatchedTrackJet", "; # of track jets ; Events", 15, 0, 15);  
    m_hSvc.create1D("nFlatBtaggedTrackJet", "; # of track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("btaggedEventType", "; type ; Events", 5, 0, 5);  
    m_hSvc.create1D("bmatchedEventType", "; type ; Events", 5, 0, 5);  

    m_hSvc.create1D("initialTrackJetPt", "; Track jet p_{T} [GeV]; Events", 25, 0, 500);
    m_hSvc.create1D("initialTrackJetEta", "; Track jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("initialTrackJetPhi", "; Track jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("initialTrackJetMv2c10", "; Track jet Mv2c10 ; Events", 20, -1, 1);
    m_hSvc.create1D("drTrackJetMuon", "; #DeltaR_{mu,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drTrackJetElectron", "; #DeltaR_{e,tjet}; Events", 30, 0, 3);
    m_hSvc.create2D("muonPtVsDrTrackJetMuon", "; muon p_[T];#DeltaR_{mu,tjet}",50,0,1000, 60, 0, 3);
    m_hSvc.create2D("eletrconPtVsDrTrackJetElectron", ";electron p_{T}; #DeltaR_{e,tjet}",50,0,1000, 60, 0, 3);
    m_hSvc.create2D("muonPtVsDrInitialTrackJetMuon", "; muon p_[T];#DeltaR_{mu,tjet}",50,0,1000, 60, 0, 3);
    m_hSvc.create2D("muonPtVsDrBmatchedInitialTrackJetMuon", "; muon p_[T];#DeltaR_{mu,tjet}",50,0,1000, 60, 0, 3);
    m_hSvc.create2D("muonPtVsDrNotBmatchedInitialTrackJetMuon", "; muon p_[T];#DeltaR_{mu,tjet}",50,0,1000, 60, 0, 3);
    m_hSvc.create2D("eletrconPtVsDrInitialTrackJetElectron", ";electron p_{T}; #DeltaR_{e,tjet}",50,0,1000, 60, 0, 3);
    m_hSvc.create2D("eletrconPtVsDrBmatchedInitialTrackJetElectron", ";electron p_{T}; #DeltaR_{e,tjet}",50,0,1000, 60, 0, 3);
    m_hSvc.create2D("eletrconPtVsDrNotBmatchedInitialTrackJetElectron", ";electron p_{T}; #DeltaR_{e,tjet}",50,0,1000, 60, 0, 3);
    m_hSvc.create1D("drInitialTrackJetMuon", "; #DeltaR_{mu,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drInitialTrackJetElectron", "; #DeltaR_{e,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drInitialBtaggedTrackJetMuon", "; #DeltaR_{mu,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drInitialBtaggedTrackJetElectron", "; #DeltaR_{e,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("nInitialTrackJet", "; # of track jets ; Events", 15, 0, 15);  
    m_hSvc.create1D("nInitialBtaggedTrackJet", "; # of track jets ; Events", 10, 0, 10);  

    m_hSvc.create1D("met", "; missing E_[T] [GeV]; Events", 25, 0, 1000);
    m_hSvc.create1D("metPhi", "; missing E_{T} #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("mwt", "; M_{T,W} [GeV]; Events", 25, 0, 250);
    m_hSvc.create1D("wpt", "; p_{T,W} [GeV]; Events", 25, 0, 2000);

    m_hSvc.create1D("deltaR_lepSmallJet", "; #DeltaR_{lep,sjet}; Events", 20, 0, 2);
    m_hSvc.create1D("deltaR_largeJetSmallJet", "; #DeltaR_{ljet,sjet}; Events", 50, 1, 6);
    m_hSvc.create1D("deltaPhi_largeJetLepton", "; #Delta#phi_{ljet,lep}; Events", 24, 2, 3.2);

    std::string channel[6] = {"All","FullHad","FullLep","SemiLepTau","SemiLepMu","SemiLepEl"};
    for(int ch=0;ch<6;ch++) {
      m_hSvc.create1D(Form("cutFlow%s",channel[ch].c_str()), ";Number of events; step", 20,0,20);
      m_hSvc.create1D(Form("mtt%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mlt%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mht%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; Events", 50,0,500);
    }
    m_hSvc.create1DVar("mtt", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
  }

AnaTtresMetBtag::~AnaTtresMetBtag() {
}

//void AnaTtresMetBtag::run(const Event &evt, double weight, const std::string &s, int analysis_mode, float totalEvents) {
void AnaTtresMetBtag::run(const Event &evt, double weight, const std::string &s, int analysis_mode, TSpline* m_spline) {
  //std::cout << "analysis mode = " << analysis_mode << std::endl;

  int w1h_pdgId = evt.MC_w1h_pdgId();
  int w2h_pdgId = evt.MC_w2h_pdgId();
  int w1l_pdgId = evt.MC_w1l_pdgId();
  int w2l_pdgId = evt.MC_w2l_pdgId();
  int Wdecay1_from_t_pdgId = evt.MC_Wdecay1_from_t_pdgId();
  int Wdecay2_from_t_pdgId = evt.MC_Wdecay2_from_t_pdgId();
  int Wdecay1_from_tbar_pdgId = evt.MC_Wdecay1_from_tbar_pdgId();
  int Wdecay2_from_tbar_pdgId = evt.MC_Wdecay2_from_tbar_pdgId();
  bool isLeptonicT = (abs(Wdecay1_from_t_pdgId)<7)? false : true;
  bool isLeptonicTbar = (abs(Wdecay1_from_tbar_pdgId)<7)? false : true;
  TLorentzVector b_from_t,b_from_tbar,t,tbar,W_from_t,W_from_tbar;
  b_from_t = evt.MC_b_from_t();
  b_from_tbar = evt.MC_b_from_tbar();
  W_from_t = evt.MC_Wdecay1_from_t() + evt.MC_Wdecay2_from_t(); 
  W_from_tbar = evt.MC_Wdecay1_from_tbar() + evt.MC_Wdecay2_from_tbar(); 
  t = evt.MC_t();
  tbar = evt.MC_tbar();
  float t_pt = t.Perp()*1e-3, t_eta = t.Eta(), t_phi = t.Phi();
  float tbar_pt = tbar.Perp()*1e-3, tbar_eta = tbar.Eta(), tbar_phi = tbar.Phi();
  float b_from_t_pt = b_from_t.Perp()*1e-3, b_from_t_eta = b_from_t.Eta(), b_from_t_phi = b_from_t.Phi();
  float b_from_tbar_pt = b_from_tbar.Perp()*1e-3, b_from_tbar_eta = b_from_tbar.Eta(), b_from_tbar_phi = b_from_tbar.Phi();

  bool isFullHad = (!isLeptonicT && !isLeptonicTbar);
  bool isFullLep = (isLeptonicT && isLeptonicTbar);
  bool isSemiLep = (!isFullHad && !isFullLep)? true : false;
  bool containTau = (Wdecay2_from_t_pdgId==-15 || Wdecay1_from_tbar_pdgId==15)? true : false;
  bool containMu = (Wdecay2_from_t_pdgId==-13 || Wdecay1_from_tbar_pdgId==13)? true : false;
  bool containEl = (Wdecay2_from_t_pdgId==-11 || Wdecay1_from_tbar_pdgId==11)? true : false;

  //if(isFullHad || isFullLep || containTau) return;

  std::string channel;
  TLorentzVector semilepmu;
  if(isFullHad) channel="FullHad";
  else if(isFullLep) channel="FullLep";
  else if(containTau) channel="SemiLepTau";
  else if(containMu) {
    channel="SemiLepMu";
    if(Wdecay2_from_t_pdgId==-13) semilepmu = evt.MC_Wdecay2_from_t(); 
    else if(Wdecay1_from_tbar_pdgId==13) semilepmu = evt.MC_Wdecay1_from_tbar(); 
  }
  else if(containEl) channel="SemiLepEl";
  else std::cout << "Unknow channel" << std::endl;
  //std::cout << "channel is " << channel << std::endl;

  HistogramService *h = &m_hSvc;

  bool isTight = false;

  std::string suffix = s;

  //MET trigger
  bool pass_met_trigger = evt.trigger("HLT_xe120_mht"); 
  //std::cout << "pass_met_trigger=" << pass_met_trigger << std::endl;
  //Muon trigger
  bool pass_muon_trigger = evt.trigger("HLT_mu26_ivarmedium") || evt.trigger("HLT_mu50"); 
  //std::cout << "pass muon trigger=" << pass_muon_trigger << std::endl;

  //Fixed & Flat B-tag
  int nTrackBtagJets=0,nTrackFlatBtagJets=0;
  float temp_thr=-1;
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    float pt = evt.tjet()[jidx].mom().Pt()*1e-3;
    float eta = evt.tjet()[jidx].mom().Eta();
    float phi = evt.tjet()[jidx].mom().Phi();
    if(pt<300) temp_thr = 0.6455;
    else temp_thr = -0.3;
    bool is_trackBtag = evt.tjet()[jidx].btag_mv2c10_70_trk();
    bool is_trackFlatBtag = evt.tjet()[jidx].mv2c10() > temp_thr;
    if (is_trackBtag && pt > 10 && fabs(eta) < 2.5 && evt.tjet()[jidx].numConstituents() >= 2) nTrackBtagJets++; 
    if (is_trackFlatBtag && pt > 10 && fabs(eta) < 2.5 && evt.tjet()[jidx].numConstituents() >= 2) nTrackFlatBtagJets++; 
  }
  bool pass_fixed_btag = nTrackBtagJets > 0 ;
  bool pass_flat_btag = nTrackFlatBtagJets > 0 ;
  //std::cout << "pass b-tag fixed/flat=" << pass_fixed_btag << "/" << pass_flat_btag << std::endl;
  //std::cout << "nTrackBtagJets=" << nTrackBtagJets << std::endl;

  //Calo B-tag
  int nCaloBtagJets=0,nCaloHybBtagJets=0;
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    float pt = evt.jet()[jidx].mom().Perp()*1e-3;
    float mv2c10 = evt.jet()[jidx].mv2c10();
    bool is_caloBtag = mv2c10>0.824;
    bool is_hybCaloBtag = mv2c10>m_spline->Eval(pt);
    if(is_caloBtag) nCaloBtagJets++;
    if(is_hybCaloBtag) nCaloHybBtagJets++;
  }
  bool pass_calo_btag = nCaloBtagJets > 0 ;
  bool pass_hyb_calo_btag = nCaloHybBtagJets > 0 ;

  //b-tag for no OR track jet
  int nBtaggedNoORTrackJets=0,nFlatBtaggedNoORTrackJets=0;
  for(int itj=0;itj<evt.nocalibtjet().size();itj++){
    TLorentzVector nocalibtjet = evt.nocalibtjet()[itj].mom();
    if(nocalibtjet.Perp()*1e-3<300) temp_thr = 0.6455;
    else temp_thr = -0.3;
    if (nocalibtjet.Perp()*1e-3 < 10 || fabs(nocalibtjet.Eta()) > 2.5 || evt.nocalibtjet()[itj].numConstituents() < 2) continue; 
    if (evt.nocalibtjet()[itj].mv2c10()>0.6455) nBtaggedNoORTrackJets++; 
    if (evt.nocalibtjet()[itj].mv2c10()>temp_thr) nFlatBtaggedNoORTrackJets++; 
  }
  bool pass_noor_btag = nBtaggedNoORTrackJets > 0 ;
  bool pass_noor_flat_btag = nFlatBtaggedNoORTrackJets > 0 ;
  //std::cout << "nBtaggedNoORTrackJets=" << nBtaggedNoORTrackJets << std::endl;

  //wpt & muon trigger matching
  float tmp_wpt=0;
  bool pass_muon_trigmatch=false;
  if(evt.muon().size()!=0){
    TLorentzVector mu = evt.muon()[0].mom();
    TLorentzVector met = evt.met();
    tmp_wpt = (mu+met).Perp()*1e-3;
    pass_muon_trigmatch = evt.muon()[0].HLT_mu26_ivarmedium() || evt.muon()[0].HLT_mu50();
  }
  //std::cout << "wpt=" << tmp_wpt << std::endl;
  //std::cout << "pass muon trigger matching=" << pass_muon_trigmatch << std::endl;

  if(analysis_mode==0){//default
    if(m_electron && !evt.passes("bejets_2016")) return;
    if(!m_electron && !evt.passes("bmujets_2016")) return;
  }
  else if (analysis_mode==1){//use met trigger
    if(m_electron && !evt.passes("bejets_2016")) return;
    if(!m_electron && !evt.passes("bmujets_notrigger_nobtag")) return;
    if(!m_electron){//met trigger is used for muon channel only
      if(tmp_wpt<220){ //use muon trigger
        if(!pass_muon_trigger || !pass_muon_trigmatch) return;
      } else{ //use met trigger
        if(!pass_met_trigger) return;
      }
    }
    if(!pass_fixed_btag) return;
  }
  else if (analysis_mode==2){//use flat b-tag
    if(m_electron && !evt.passes("bejets_nobtag")) return;
    if(!m_electron && !evt.passes("bmujets_notrigger_nobtag")) return;
    if(!m_electron)
      if(!pass_muon_trigger || !pass_muon_trigmatch) return;
    if(!pass_flat_btag) return;
  }
  else if (analysis_mode==3){//use met trigger & flat b-tag
    if(m_electron && !evt.passes("bejets_nobtag")) return;
    if(!m_electron && !evt.passes("bmujets_notrigger_nobtag")) return;
    if(!m_electron){//met trigger is used for muon channel only
      if(tmp_wpt<220){ //use muon trigger
        if(!pass_muon_trigger || !pass_muon_trigmatch) return;
      } else{ //use met trigger
        if(!pass_met_trigger) return;
      }
    }
    if(!pass_flat_btag) return;
  }
  else if (analysis_mode==4){//NO OR for track jet
    //std::cout << "pass e/mu=" << evt.passes("bejets_nobtag") << "/" << evt.passes("bmujets_notrigger_nobtag") << std::endl;
    if(m_electron && !evt.passes("bejets_nobtag")) return;
    if(!m_electron && !evt.passes("bmujets_notrigger_nobtag")) return;
    if(!m_electron){
      if(!pass_muon_trigger || !pass_muon_trigmatch) return;
    }
    if(!pass_noor_btag) return;
  }
  else if (analysis_mode==5){//NO OR for track jet & flat b-tag
    //std::cout << "pass e/mu=" << evt.passes("bejets_nobtag") << "/" << evt.passes("bmujets_notrigger_nobtag") << std::endl;
    if(m_electron && !evt.passes("bejets_nobtag")) return;
    if(!m_electron && !evt.passes("bmujets_notrigger_nobtag")) return;
    if(!m_electron){
      if(!pass_muon_trigger || !pass_muon_trigmatch) return;
    }
    if(!pass_noor_flat_btag) return;
  }
  else if (analysis_mode==6){//use MET trigger & NO OR for track jet & flat b-tag
    //std::cout << "pass e/mu=" << evt.passes("bejets_nobtag") << "/" << evt.passes("bmujets_notrigger_nobtag") << std::endl;
    if(m_electron && !evt.passes("bejets_nobtag")) return;
    if(!m_electron && !evt.passes("bmujets_notrigger_nobtag")) return;
    if(!m_electron){//met trigger is used for muon channel only
      if(tmp_wpt<220){ //use muon trigger
        if(!pass_muon_trigger || !pass_muon_trigmatch) return;
      } else{ //use met trigger
        if(!pass_met_trigger) return;
      }
    }
    if(!pass_noor_flat_btag) return;
  }
  else if (analysis_mode==7){//validation
    //std::cout << "pass e/mu=" << evt.passes("bejets_nobtag") << "/" << evt.passes("bmujets_notrigger_nobtag") << std::endl;
    if(m_electron && !evt.passes("bejets_nobtag")) return;
    if(!m_electron && !evt.passes("bmujets_notrigger_nobtag")) return;
    if(!m_electron){
      if(!pass_muon_trigger || !pass_muon_trigmatch) return;
    }
    //if(!pass_fixed_btag) return;
  }
  else if (analysis_mode==8 || analysis_mode==9){//use calo b-tag
    if(m_electron && !evt.passes("bejets_nobtag")) return;
    if(!m_electron && !evt.passes("bmujets_notrigger_nobtag")) return;
    if(!m_electron)
      if(!pass_muon_trigger || !pass_muon_trigmatch) return;
    if(analysis_mode==8 && !pass_calo_btag) return;
    if(analysis_mode==9 && !pass_hyb_calo_btag) return;
  }
  //std::cout << "make plot!!" << std::endl;

  h->h1D("cutFlowAll", "", suffix)->Fill(0);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(0);

  TLorentzVector l;
  if(m_electron) {//electron
    l = evt.electron()[0].mom();
    h->h1D("nLepton", "", suffix)->Fill(evt.electron().size(), weight);
  }
  else {//muon
    l = evt.muon()[0].mom();
    h->h1D("nLepton", "", suffix)->Fill(evt.muon().size(), weight);
  }
  h->h1D("lepPt", "", suffix)->Fill(l.Perp()*1e-3, weight);
  h->h1D("lepEta", "", suffix)->Fill(l.Eta(), weight);
  h->h1D("lepPhi", "", suffix)->Fill(l.Phi(), weight);

  TLorentzVector vmet = evt.met();
  float mwt = sqrt(2. * l.Perp() * vmet.Perp() * (1. - cos(vmet.DeltaPhi(l)) ));
  float wpt = (vmet+l).Perp();
  h->h1D("met", "", suffix)->Fill(vmet.Perp()*1e-3, weight);
  h->h1D("mwt", "", suffix)->Fill(mwt*1e-3, weight);
  h->h1D("wpt", "", suffix)->Fill(wpt*1e-3, weight);
  h->h1D("metPhi", "", suffix)->Fill(vmet.Phi(), weight);
  h->h1D("metPlusMwt", "", suffix)->Fill((vmet.Perp()+mwt)*1e-3, weight);

  size_t close_idx = 0;
  for (; close_idx < evt.jet().size(); ++close_idx){
    TLorentzVector sjet = evt.jet()[close_idx].mom();
    if (evt.jet()[close_idx].closeToLepton())
      break;
  }
  const TLorentzVector &sj = evt.jet()[close_idx].mom();
  const TLorentzVector &lead_sj = evt.jet()[0].mom();
  const TLorentzVector &second_lead_sj = evt.jet()[1].mom();
  h->h1D("nSmallJet", "", suffix)->Fill(evt.jet().size(), weight);
  h->h1D("smallJetPt", "", suffix)->Fill(sj.Perp()*1e-3, weight);
  h->h1D("smallJetEta", "", suffix)->Fill(sj.Eta(), weight);
  h->h1D("smallJetPhi", "", suffix)->Fill(sj.Phi(), weight);
  h->h1D("smallJetMass", "", suffix)->Fill(sj.M()*1e-3, weight);
  h->h1D("leadSmallJetPt", "", suffix)->Fill(lead_sj.Perp()*1e-3, weight);
  h->h1D("leadSmallJetEta", "", suffix)->Fill(lead_sj.Eta(), weight);
  h->h1D("leadSmallJetPhi", "", suffix)->Fill(lead_sj.Phi(), weight);
  h->h1D("leadSmallJetMass", "", suffix)->Fill(lead_sj.M()*1e-3, weight);
  h->h1D("secondLeadSmallJetPt", "", suffix)->Fill(second_lead_sj.Perp()*1e-3, weight);
  h->h1D("secondLeadSmallJetEta", "", suffix)->Fill(second_lead_sj.Eta(), weight);
  h->h1D("secondLeadSmallJetPhi", "", suffix)->Fill(second_lead_sj.Phi(), weight);
  h->h1D("secondLeadSmallJetMass", "", suffix)->Fill(second_lead_sj.M()*1e-3, weight);

  size_t goodljet_idx = 0;
  /*for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx) {
    TLorentzVector ljet = evt.largeJet()[goodljet_idx].mom();
    std::cout << "ljet pt/eta/phi=" << ljet.Perp()*1e-3 << "/" << ljet.Eta() << "/" << ljet.Phi() << std::endl;
    std::cout << "good=" << evt.largeJet()[goodljet_idx].good() << std::endl;
    if (evt.largeJet()[goodljet_idx].good())
      break;
  }*/

  //temporal because of missing top tag information
  for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx) {
    const TLorentzVector &lj = evt.largeJet()[goodljet_idx].mom();
    if(lj.Perp()*1e-3 < 300 || fabs(lj.Eta())>2.0) continue;
    float delta_phi_large_jet_lepton = lj.DeltaPhi(l);
    float delta_r_large_jet_small_jet = lj.DeltaR(sj);
    //std::cout << "dphi_lj_lep/dr_lj_sj=" << delta_phi_large_jet_lepton << "/" << delta_r_large_jet_small_jet << std::endl;
    if(fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) break;
  }
  
  const TLorentzVector &ljt = evt.largeJet()[goodljet_idx].mom();
  h->h1D("nLargeJet", "", suffix)->Fill(evt.largeJet().size(), weight);
  h->h1D("largeJetPt", "", suffix)->Fill(ljt.Perp()*1e-3, weight);
  h->h1D("largeJetEta", "", suffix)->Fill(ljt.Eta(), weight);
  h->h1D("largeJetPhi", "", suffix)->Fill(ljt.Phi(), weight);
  h->h1D("largeJetMass", "", suffix)->Fill(ljt.M()*1e-3, weight);
  //h->h1D("largeJetTau32", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau32_wta"), weight);

  h->h1D("deltaR_lepSmallJet", "",suffix)->Fill(l.DeltaR(sj),weight);
  h->h1D("deltaR_largeJetSmallJet", "",suffix)->Fill(ljt.DeltaR(sj),weight);
  h->h1D("deltaPhi_largeJetLepton", "",suffix)->Fill(fabs(ljt.DeltaPhi(l)),weight);

  //track jet after object selection
  int nBmatchedTrackJetLep=0,nNotBmatchedTrackJet=0,nBtaggedBmatchedTrackJetLep=0;
  int nFlatBtaggedBmatchedTrackJetLep=0,nFlatBtaggedNotBmatchedTrackJet=0;
  int nBmatchedTrackJetHad=0,nBtaggedBmatchedTrackJetHad=0;
  int nFlatBtaggedBmatchedTrackJetHad=0,nBtaggedNotBmatchedTrackJet=0;
  for(int itj=0;itj<evt.tjet().size();itj++){
    TLorentzVector tjet = evt.tjet()[itj].mom();
    float thr = (tjet.Perp()*1e-3>300)? -0.3 : 0.6455;
    h->h1D("trackJetPt", "", suffix)->Fill(tjet.Perp()*1e-3, weight);
    h->h1D("trackJetEta", "", suffix)->Fill(tjet.Eta(), weight);
    h->h1D("trackJetPhi", "", suffix)->Fill(tjet.Phi(), weight);
    h->h1D("trackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
    if(isLeptonicT){
      float deta = tjet.Eta() - b_from_t.Eta();
      float dphi = acos(cos(tjet.Phi() - b_from_t.Phi()));
      float dr = sqrt(deta*deta+dphi*dphi);
      if(dr<0.2) {
        nBmatchedTrackJetLep++;
        if (evt.tjet()[itj].mv2c10()>0.6455) nBtaggedBmatchedTrackJetLep++; 
        if (evt.tjet()[itj].mv2c10()>thr) nFlatBtaggedBmatchedTrackJetLep++; 
        h->h1D("bmatchedTrackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
      }
      else {
        nNotBmatchedTrackJet++;
        if (evt.tjet()[itj].mv2c10()>0.6455) nBtaggedNotBmatchedTrackJet++; 
        if (evt.tjet()[itj].mv2c10()>thr) nFlatBtaggedNotBmatchedTrackJet++; 
        h->h1D("notBmatchedTrackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
      }
    }
    else{
      float deta = tjet.Eta() - b_from_t.Eta();
      float dphi = acos(cos(tjet.Phi() - b_from_t.Phi()));
      float dr = sqrt(deta*deta+dphi*dphi);
      if(dr<0.2) {
        nBmatchedTrackJetHad++;
        if (evt.tjet()[itj].mv2c10()>0.6455) nBtaggedBmatchedTrackJetHad++; 
        if (evt.tjet()[itj].mv2c10()>thr) nFlatBtaggedBmatchedTrackJetHad++; 
        h->h1D("bmatchedTrackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
      }
      else {
        nNotBmatchedTrackJet++;
        if (evt.tjet()[itj].mv2c10()>0.6455) nBtaggedNotBmatchedTrackJet++; 
        if (evt.tjet()[itj].mv2c10()>thr) nFlatBtaggedNotBmatchedTrackJet++; 
        h->h1D("notBmatchedTrackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
      }
    }
    if(isLeptonicTbar){
      float deta = tjet.Eta() - b_from_tbar.Eta();
      float dphi = acos(cos(tjet.Phi() - b_from_tbar.Phi()));
      float dr = sqrt(deta*deta+dphi*dphi);
      if(dr<0.2) {
        nBmatchedTrackJetLep++;
        if (evt.tjet()[itj].mv2c10()>0.6455) nBtaggedBmatchedTrackJetLep++; 
        if (evt.tjet()[itj].mv2c10()>thr) nFlatBtaggedBmatchedTrackJetLep++; 
        h->h1D("bmatchedTrackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
      }
      else {
        nNotBmatchedTrackJet++;
        if (evt.tjet()[itj].mv2c10()>0.6455) nBtaggedNotBmatchedTrackJet++; 
        if (evt.tjet()[itj].mv2c10()>thr) nFlatBtaggedNotBmatchedTrackJet++; 
        h->h1D("notBmatchedTrackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
      }
    }
    else{
      float deta = tjet.Eta() - b_from_tbar.Eta();
      float dphi = acos(cos(tjet.Phi() - b_from_tbar.Phi()));
      float dr = sqrt(deta*deta+dphi*dphi);
      if(dr<0.2) {
        nBmatchedTrackJetHad++;
        if (evt.tjet()[itj].mv2c10()>0.6455) nBtaggedBmatchedTrackJetHad++; 
        if (evt.tjet()[itj].mv2c10()>thr) nFlatBtaggedBmatchedTrackJetHad++; 
        h->h1D("bmatchedTrackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
      }
      else {
        nNotBmatchedTrackJet++;
        if (evt.tjet()[itj].mv2c10()>0.6455) nBtaggedNotBmatchedTrackJet++; 
        if (evt.tjet()[itj].mv2c10()>thr) nFlatBtaggedNotBmatchedTrackJet++; 
        h->h1D("notBmatchedTrackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
      }
    }
    for(int im=0;im<evt.muon().size();im++){
      TLorentzVector mu = evt.muon()[im].mom();
      float deta = mu.Eta() - tjet.Eta();
      float dphi = acos(cos(mu.Phi() - tjet.Phi() ));
      float dr = sqrt(deta*deta+dphi*dphi);
      h->h1D("drTrackJetMuon", "",suffix)->Fill(dr,weight);
      h->h2D("muonPtVsDrTrackJetMuon", "",suffix)->Fill(mu.Perp()*1e-3,dr);
    }
    for(int ie=0;ie<evt.electron().size();ie++){
      TLorentzVector el = evt.electron()[ie].mom();
      float deta = el.Eta() - tjet.Eta();
      float dphi = acos(cos(el.Phi() - tjet.Phi() ));
      float dr = sqrt(deta*deta+dphi*dphi);
      h->h1D("drTrackJetElectron", "",suffix)->Fill(dr,weight);
      h->h2D("electronPtVsDrTrackJetElectron", "",suffix)->Fill(el.Perp()*1e-3,dr);
    }
  }
  h->h1D("nTrackJet", "", suffix)->Fill(evt.tjet().size(), weight);
  h->h1D("nBtaggedTrackJet", "", suffix)->Fill(nTrackBtagJets, weight);
  h->h1D("nFlatBtaggedTrackJet", "", suffix)->Fill(nTrackFlatBtagJets, weight);
  h->h1D("nBmatchedTrackJet", "", suffix)->Fill(nBmatchedTrackJetLep+nBmatchedTrackJetHad, weight);
  h->h1D("nBtaggedBmatchedTrackJet", "", suffix)->Fill(nBtaggedBmatchedTrackJetLep+nBtaggedBmatchedTrackJetHad, weight);
  h->h1D("nBtaggedNotBmatchedTrackJet", "", suffix)->Fill(nBtaggedNotBmatchedTrackJet, weight);
  h->h1D("nFlatBtaggedBmatchedTrackJet", "", suffix)->Fill(nFlatBtaggedBmatchedTrackJetLep+nFlatBtaggedBmatchedTrackJetHad, weight);
  h->h1D("nFlatBtaggedNotBmatchedTrackJet", "", suffix)->Fill(nFlatBtaggedNotBmatchedTrackJet, weight);
  h->h1D("nNotBmatchedTrackJet", "", suffix)->Fill(nNotBmatchedTrackJet, weight);
  h->h1D("nBmatchedTrackJetLeptonic", "", suffix)->Fill(nBmatchedTrackJetLep, weight);
  h->h1D("nBtaggedBmatchedTrackJetLeptonic", "", suffix)->Fill(nBtaggedBmatchedTrackJetLep, weight);
  h->h1D("nFlatBtaggedBmatchedTrackJetLeptonic", "", suffix)->Fill(nFlatBtaggedBmatchedTrackJetLep, weight);
  h->h1D("nBmatchedTrackJetHadronic", "", suffix)->Fill(nBmatchedTrackJetHad, weight);
  h->h1D("nBtaggedBmatchedTrackJetHadronic", "", suffix)->Fill(nBtaggedBmatchedTrackJetHad, weight);
  h->h1D("nFlatBtaggedBmatchedTrackJetHadronic", "", suffix)->Fill(nFlatBtaggedBmatchedTrackJetHad, weight);
  if(nBtaggedBmatchedTrackJetLep==0 && nBtaggedBmatchedTrackJetHad==0) h->h1D("btaggedEventType", "", suffix)->Fill(0., weight);
  else if(nBtaggedBmatchedTrackJetLep==0 && nBtaggedBmatchedTrackJetHad>0) h->h1D("btaggedEventType", "", suffix)->Fill(1., weight);
  else if(nBtaggedBmatchedTrackJetLep>0 && nBtaggedBmatchedTrackJetHad==0) h->h1D("btaggedEventType", "", suffix)->Fill(2., weight);
  else if(nBtaggedBmatchedTrackJetLep>0 && nBtaggedBmatchedTrackJetHad>0) h->h1D("btaggedEventType", "", suffix)->Fill(3., weight);
  else h->h1D("btaggedEventType", "", suffix)->Fill(4., weight);
  if(nBmatchedTrackJetLep==0 && nBmatchedTrackJetHad==0) h->h1D("bmatchedEventType", "", suffix)->Fill(0., weight);
  else if(nBmatchedTrackJetLep==0 && nBmatchedTrackJetHad>0) {
    h->h1D("bmatchedEventType", "", suffix)->Fill(1., weight);
    h->h1D("nBtaggedBmatchedTrackJet1bLeptonic", "", suffix)->Fill(nBtaggedBmatchedTrackJetLep+nBtaggedBmatchedTrackJetHad, weight);
    h->h1D("nFlatBtaggedBmatchedTrackJet1bLeptonic", "", suffix)->Fill(nFlatBtaggedBmatchedTrackJetLep+nFlatBtaggedBmatchedTrackJetHad, weight);
  }
  else if(nBmatchedTrackJetLep>0 && nBmatchedTrackJetHad==0) h->h1D("bmatchedEventType", "", suffix)->Fill(2., weight);
  else if(nBmatchedTrackJetLep>0 && nBmatchedTrackJetHad>0) {
    h->h1D("bmatchedEventType", "", suffix)->Fill(3., weight);
    h->h1D("nBtaggedBmatchedTrackJet2b", "", suffix)->Fill(nBtaggedBmatchedTrackJetLep+nBtaggedBmatchedTrackJetHad, weight);
    h->h1D("nFlatBtaggedBmatchedTrackJet2b", "", suffix)->Fill(nFlatBtaggedBmatchedTrackJetLep+nFlatBtaggedBmatchedTrackJetHad, weight);
  }
  else h->h1D("bmatchedEventType", "", suffix)->Fill(4., weight);

  //track jet before object selection
  int nInitialTrackJet=0,nInitialTrackBtagJets=0;
  for(int itj=0;itj<evt.nocalibtjet().size();itj++){
    TLorentzVector nocalibtjet = evt.nocalibtjet()[itj].mom();
    if (nocalibtjet.Perp()*1e-3 < 10 || fabs(nocalibtjet.Eta()) > 2.5 || evt.nocalibtjet()[itj].numConstituents() < 2) continue; 
    h->h1D("initialTrackJetPt", "", suffix)->Fill(nocalibtjet.Perp()*1e-3, weight);
    h->h1D("initialTrackJetEta", "", suffix)->Fill(nocalibtjet.Eta(), weight);
    h->h1D("initialTrackJetPhi", "", suffix)->Fill(nocalibtjet.Phi(), weight);
    h->h1D("initialTrackJetMv2c10", "", suffix)->Fill(evt.nocalibtjet()[itj].mv2c10(), weight);
    nInitialTrackJet++;
    if (evt.nocalibtjet()[itj].mv2c10()>0.6455) nInitialTrackBtagJets++; 
    for(int im=0;im<evt.muon().size();im++){
      TLorentzVector mu = evt.muon()[im].mom();
      float deta = mu.Eta() - nocalibtjet.Eta();
      float dphi = acos(cos(mu.Phi() - nocalibtjet.Phi() ));
      float dr = sqrt(deta*deta+dphi*dphi);
      h->h1D("drInitialTrackJetMuon", "",suffix)->Fill(dr,weight);
      h->h2D("muonPtVsDrInitialTrackJetMuon", "",suffix)->Fill(mu.Perp()*1e-3,dr);
      if (evt.nocalibtjet()[itj].mv2c10()>0.6455) h->h1D("drInitialBtaggedTrackJetMuon", "",suffix)->Fill(dr,weight);
      if(isLeptonicT){
        float bdeta = nocalibtjet.Eta() - b_from_t.Eta();
        float bdphi = acos(cos(nocalibtjet.Phi() - b_from_t.Phi()));
        float bdr = sqrt(bdeta*bdeta+bdphi*bdphi);
        if(bdr<0.2) h->h2D("muonPtVsDrBmatchedInitialTrackJetMuon", "",suffix)->Fill(mu.Perp()*1e-3,dr);
        else h->h2D("muonPtVsDrNotBmatchedInitialTrackJetMuon", "",suffix)->Fill(mu.Perp()*1e-3,dr);
      }
      if(isLeptonicTbar){
        float bdeta = nocalibtjet.Eta() - b_from_tbar.Eta();
        float bdphi = acos(cos(nocalibtjet.Phi() - b_from_tbar.Phi()));
        float bdr = sqrt(bdeta*bdeta+bdphi*bdphi);
        if(bdr<0.2) h->h2D("muonPtVsDrBmatchedInitialTrackJetMuon", "",suffix)->Fill(mu.Perp()*1e-3,dr);
        else h->h2D("muonPtVsDrNotBmatchedInitialTrackJetMuon", "",suffix)->Fill(mu.Perp()*1e-3,dr);
      }
    }
    for(int ie=0;ie<evt.electron().size();ie++){
      TLorentzVector el = evt.electron()[ie].mom();
      float deta = el.Eta() - nocalibtjet.Eta();
      float dphi = acos(cos(el.Phi() - nocalibtjet.Phi() ));
      float dr = sqrt(deta*deta+dphi*dphi);
      h->h1D("drInitialTrackJetElectron", "",suffix)->Fill(dr,weight);
      h->h2D("electronPtVsDrInitialTrackJetElectron", "",suffix)->Fill(el.Perp()*1e-3,dr);
      if (evt.nocalibtjet()[itj].mv2c10()>0.6455) h->h1D("drInitialBtaggedTrackJetElectron", "",suffix)->Fill(dr,weight);
      if(isLeptonicT){
        float bdeta = nocalibtjet.Eta() - b_from_t.Eta();
        float bdphi = acos(cos(nocalibtjet.Phi() - b_from_t.Phi()));
        float bdr = sqrt(bdeta*bdeta+bdphi*bdphi);
        if(bdr<0.2) h->h2D("muonPtVsDrBmatchedInitialTrackJetMuon", "",suffix)->Fill(el.Perp()*1e-3,dr);
        else h->h2D("muonPtVsDrNotBmatchedInitialTrackJetMuon", "",suffix)->Fill(el.Perp()*1e-3,dr);
      }
      if(isLeptonicTbar){
        float bdeta = nocalibtjet.Eta() - b_from_tbar.Eta();
        float bdphi = acos(cos(nocalibtjet.Phi() - b_from_tbar.Phi()));
        float bdr = sqrt(bdeta*bdeta+bdphi*bdphi);
        if(bdr<0.2) h->h2D("electronPtVsDrBmatchedInitialTrackJetElectron", "",suffix)->Fill(el.Perp()*1e-3,dr);
        else h->h2D("electronPtVsDrNotBmatchedInitialTrackJetElectron", "",suffix)->Fill(el.Perp()*1e-3,dr);
      }
    }
  }
  h->h1D("nInitialTrackJet", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nInitialBtaggedTrackJet", "", suffix)->Fill(nInitialTrackBtagJets, weight);

  std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), true);   
  TLorentzVector nu(0,0,0,0);
  if (vec_nu.size() > 0) {
    nu = *(vec_nu[0]);
    for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
    vec_nu.clear();
  }
  float mtt = (ljt+sj+nu+l).M();
  float mlt = (sj+nu+l).M();
  float mht = ljt.M();
  float ljet_pt = ljt.Perp();
  float ljet_eta = ljt.Eta();
  float ljet_phi = ljt.Phi();
  h->h1D("mtt", "", suffix)->Fill(mtt*1e-3, weight);
  h->h1D("mttAll", "", suffix)->Fill(mtt*1e-3, weight);
  h->h1D("mltAll", "", suffix)->Fill(mlt*1e-3, weight);
  h->h1D("mhtAll", "", suffix)->Fill(mht*1e-3, weight);
  h->h1D(Form("mtt%s",channel.c_str()), "", suffix)->Fill(mtt*1e-3, weight);
  h->h1D(Form("mlt%s",channel.c_str()), "", suffix)->Fill(mlt*1e-3, weight);
  h->h1D(Form("mht%s",channel.c_str()), "", suffix)->Fill(mht*1e-3, weight);
}
