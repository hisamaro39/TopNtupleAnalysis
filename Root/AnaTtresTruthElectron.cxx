/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresTruthElectron.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"


AnaTtresTruthElectron::AnaTtresTruthElectron(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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

    double varBin8[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500,3000};
    int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;
    double varBin9[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000};
    int varBinN9 = sizeof(varBin9)/sizeof(double) - 1;

    m_hSvc.create1D("yields", "; One ; Events", 1, 0.5, 1.5);

    m_hSvc.create1DVar("lepPt", "; lepton p_{T} [GeV]; Events", varBinN1, varBin1);  
    m_hSvc.create1D("nJets", "; number of jets ; Events", 10, -0.5, 9.5);
    m_hSvc.create1D("nLargeJets", "; number of large jets ; Events", 10, -0.5, 9.5);
    m_hSvc.create2D("lepPt_vs_drMuonJet",    "; Pt of lept [GeV]; #DeltaR_{#mu,jet}", 100, 25, 525,100,0,2);

    m_hSvc.create1D("nTrkBtagJets", "; number of b-tagged track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nCaloBtagJets", "; number of b-tagged calo jets ; Events", 10, 0, 10);  

    m_hSvc.create1DVar("MET", "; missing E_{T} [GeV]; Events", varBinN1, varBin1);
    m_hSvc.create1DVar("MET_plus_MWT", "; missing E_{T} + m^{W}_{T} [GeV]; Events", varBinN1, varBin1);
    m_hSvc.create1D("MWT", "; W transverse mass [GeV]; Events", 20, 0, 200);

    m_hSvc.create1DVar("ltop_pt_jet",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_jet",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_jet_pt",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_jet_passCaloBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_jet_passCaloBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_jet_pt_passCaloBtag",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create2DVar("bMatchedTjetPt_vs_btagScore",    ";track jet pt[GeV]; btag score", varBinN9,varBin9,varBinN10,varBin10);
    m_hSvc.create1DVar("ltop_pt_tjet",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_tjet",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_tjet_pt",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_tjet_passTrackBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_tjet_passTrackBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_tjet_pt_passTrackBtag",    "; Pt of b from leptonic top matched track jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_tjet_pass85TrackBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_tjet_pass85TrackBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_tjet_pt_pass85TrackBtag",    "; Pt of b from leptonic top matched track jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_tjet_passFlatTrackBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_tjet_passFlatTrackBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_tjet_pt_passFlatTrackBtag",    "; Pt of b from leptonic top matched track jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_jet",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_jet",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_jet_pt",    "; Pt of b from hadronic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_jet_passCaloBtag",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_jet_passCaloBtag",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_jet_pt_passCaloBtag",    "; Pt of b from hadronic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_tjet",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_tjet",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_tjet_pt",    "; Pt of b from hadronic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_tjet_passTrackBtag",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_tjet_passTrackBtag",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_tjet_pt_passTrackBtag",    "; Pt of b from hadronic top matched track jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_tjet_pass85TrackBtag",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_tjet_pass85TrackBtag",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_tjet_pt_pass85TrackBtag",    "; Pt of b from hadronic top matched track jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_tjet_passFlatTrackBtag",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_tjet_passFlatTrackBtag",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_tjet_pt_passFlatTrackBtag",    "; Pt of b from hadronic top matched track jet [GeV]; Events", varBinN9, varBin9);

    m_hSvc.create1D("trackBtagBDT", "; track b-tag BDT ; Events", 100, -1, 1);
    m_hSvc.create1D("caloBtagBDT", "; calo b-tag BDT ; Events", 100, -1, 1);

    m_hSvc.create1D("small_jet_pt", "; small jet p_{T} [GeV]; Events", 50, 0, 500);
    m_hSvc.create1D("small_jet_eta", "; small jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("small_jet_phi", "; small jet #phi ; Events", 32, -3.2, 3.2);

    m_hSvc.create1DVar("largeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    m_hSvc.create1DVar("leadingLargeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    m_hSvc.create1D("largeJetEta", "; large jet #eta ; Events", 20, -2., 2.);
    m_hSvc.create1D("leadingLargeJetEta", "; large jet #eta ; Events", 20, -2., 2.);
    m_hSvc.create1DVar("TopMatchedLargeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    m_hSvc.create1DVar("TopMatchedLargeJetPt_passSTT", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    m_hSvc.create1DVar("allLargeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    m_hSvc.create1DVar("allLargeJetPt_passSTT", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    m_hSvc.create1DVar("drLargeJetTruthTop", "; #DeltaR_{lj,truth-top} ; Events", 100,0,1);
    m_hSvc.create1DVar("drJetHB", "; #DeltaR_{j,hb} ; Events", 100,0,1);
    m_hSvc.create1DVar("drTjetHB", "; #DeltaR_{tj,hb} ; Events", 100,0,1);

    m_hSvc.create1D("initialTrackJetPt", "; Track jet p_{T} [GeV]; Events", 25, 0, 500);
    m_hSvc.create1D("initialTrackJetEta", "; Track jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("initialTrackJetPhi", "; Track jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("initialTrackJetMv2c10", "; Track jet Mv2c10 ; Events", 20, -1, 1);
    m_hSvc.create1D("drTrackJetMuon", "; #DeltaR_{mu,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drTrackJetElectron", "; #DeltaR_{e,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drInitialTrackJetMuon", "; #DeltaR_{mu,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drInitialTrackJetElectron", "; #DeltaR_{e,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drInitialBtaggedTrackJetMuon", "; #DeltaR_{mu,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("drInitialBtaggedTrackJetElectron", "; #DeltaR_{e,tjet}; Events", 30, 0, 3);
    m_hSvc.create1D("nInitialBtaggedTrackJet", "; # of track jets ; Events", 10, 0, 10);  

    for (int i=0;i<14;i++){
      m_hSvc.create1D(Form("nInitialTrackJet_cut%d",i), "; # of track jets ; Events", 20, 0, 20);  
      m_hSvc.create1D(Form("nTrackJet_cut%d",i), "; # of track jets ; Events", 20, 0, 20);  
    }

    std::string bdt_thr[3] = {"loose","middle","tight"};
    for(int i=0;i<3;i++){
      m_hSvc.create1DVar(Form("allLargeJetPt_passBDT_%s",bdt_thr[i].c_str()), "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    }
    for(int lpt=0;lpt<4;lpt++){
      m_hSvc.create1D(Form("BDT_TOPtag_allljet_pt%d_%d",lpt*500,500+lpt*500),    "; BDT TOP tag; ", 100, -1, 1);
      m_hSvc.create1D(Form("BDT_TOPtag_passljet_pt%d_%d",lpt*500,500+lpt*500),    "; BDT TOP tag; ", 100, -1, 1);
    }

    std::string channel[5] = {"FullHad","FullLep","SemiLepTau","SemiLepMu","SemiLepEl"};
    for(int ch=0;ch<5;ch++) {
      m_hSvc.create1D(Form("cutFlow%s",channel[ch].c_str()), ";Number of events; step", 20,0,20);
      m_hSvc.create1D(Form("nMuon%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1D(Form("nMuonFinal%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1D(Form("nMuonNoIso%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1DVar(Form("mtt%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
    }
    m_hSvc.create1D("cutFlowAll", ";Number of events; step", 20,0,20);
    m_hSvc.create1D("nMuonAll", ";# #mu;Number of events", 5,0,5);
    m_hSvc.create1D("nMuonFinalAll", ";# #mu;Number of events", 5,0,5);
    m_hSvc.create1D("nMuonNoIsoAll", "# #mu;Number of events", 5,0,5);
    m_hSvc.create1D("muonPt", ";Number of events; p_{T,#mu}", 50,0,500);
    m_hSvc.create1D("initialMuonPt", ";Number of events; p_{T,#mu}", 50,0,500);
    m_hSvc.create1D("leadingInitialMuonPt", ";Number of events; p_{T,#mu}", 100,0,500);
    m_hSvc.create1D("electronPt", ";Number of events; p_{T,e}", 50,0,500);
    m_hSvc.create1D("nElectron", ";Number of events; # e", 5,0,5);
    m_hSvc.create1D("muon_ptvarcone30_over_pt", ";#mu ptvarcone30 / pt", 100,0,0.1);
    m_hSvc.create1DVar("mttAll", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
    m_hSvc.create1D("truthMuonPt", ";Number of events; p_{T,#mu}", 50,0,500);
  }

AnaTtresTruthElectron::~AnaTtresTruthElectron() {
}

void AnaTtresTruthElectron::run(const Event &evt, double weight, const std::string &s, int is2016run) {
  //std::cout << "AnaTtresTruthElectron::run" << std::endl;

  int w1h_pdgId = evt.MC_w1h_pdgId();
  int w2h_pdgId = evt.MC_w2h_pdgId();
  int w1l_pdgId = evt.MC_w1l_pdgId();
  int w2l_pdgId = evt.MC_w2l_pdgId();
  //std::cout << "w1h_pdgId=" << w1h_pdgId << std::endl;
  //std::cout << "w2h_pdgId=" << w2h_pdgId << std::endl;
  //std::cout << "w1l_pdgId=" << w1l_pdgId << std::endl;
  //std::cout << "w2l_pdgId=" << w2l_pdgId << std::endl;
  int Wdecay1_from_t_pdgId = evt.MC_Wdecay1_from_t_pdgId();
  int Wdecay2_from_t_pdgId = evt.MC_Wdecay2_from_t_pdgId();
  int Wdecay1_from_tbar_pdgId = evt.MC_Wdecay1_from_tbar_pdgId();
  int Wdecay2_from_tbar_pdgId = evt.MC_Wdecay2_from_tbar_pdgId();
  //std::cout << "Wdecay1_from_t_pdgId=" << Wdecay1_from_t_pdgId << std::endl;
  //std::cout << "Wdecay2_from_t_pdgId=" << Wdecay2_from_t_pdgId << std::endl;
  //std::cout << "Wdecay1_from_tbar_pdgId=" << Wdecay1_from_tbar_pdgId << std::endl;
  //std::cout << "Wdecay2_from_tbar_pdgId=" << Wdecay2_from_tbar_pdgId << std::endl;
  bool isLeptonicT = (abs(Wdecay1_from_t_pdgId)<7)? false : true;
  bool isLeptonicTbar = (abs(Wdecay1_from_tbar_pdgId)<7)? false : true;
  //std::cout << "isLeptonic top/tbar=" << isLeptonicT << "/" << isLeptonicTbar << std::endl;
  TLorentzVector b_from_t,b_from_tbar,t,tbar;
  b_from_t = evt.MC_b_from_t();
  b_from_tbar = evt.MC_b_from_tbar();
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

  bool trig1(0); 
  bool trig2(0); 
  bool trig3(0);
  bool trig4(0);
  bool trig5(0);
  bool trig6(0);

  bool isTight = false;

  std::string suffix = s;


  //track jet before object selection
  int nInitialTrackJet=0,nInitialTrackBtagJets=0;
  std::cout << "# of initial track jet is " << evt.nocalibtjet().size() << std::endl;
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
      if (evt.nocalibtjet()[itj].mv2c10()>0.6455) h->h1D("drInitialBtaggedTrackJetMuon", "",suffix)->Fill(dr,weight);
    }
    for(int ie=0;ie<evt.electron().size();ie++){
      TLorentzVector el = evt.electron()[ie].mom();
      float deta = el.Eta() - nocalibtjet.Eta();
      float dphi = acos(cos(el.Phi() - nocalibtjet.Phi() ));
      float dr = sqrt(deta*deta+dphi*dphi);
      h->h1D("drInitialTrackJetElectron", "",suffix)->Fill(dr,weight);
      if (evt.nocalibtjet()[itj].mv2c10()>0.6455) h->h1D("drInitialBtaggedTrackJetElectron", "",suffix)->Fill(dr,weight);
    }
  }
  h->h1D("nInitialTrackJet", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nInitialBtaggedTrackJet", "", suffix)->Fill(nInitialTrackBtagJets, weight);

  if(!(evt.passes("no_cut"))) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(0);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(0);
  h->h1D("nInitialTrackJet_cut0", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut0", "", suffix)->Fill(evt.tjet().size(), weight);

  //Trigger selection
  bool pass_e26_lhtight_nod0_ivarloose = evt.trigger("HLT_e26_lhtight_nod0_ivarloose"); 
  bool pass_e60_lhmedium_nod0 = evt.trigger("HLT_e60_lhmedium_nod0"); 
  bool pass_e140_lhloose_nod0 = evt.trigger("HLT_e140_lhloose_nod0"); 
  std::cout << "pass_e26_lhtight_nod0_ivarloose=" << pass_e26_lhtight_nod0_ivarloose << std::endl;
  std::cout << "pass_e60_lhmedium_nod0=" << pass_e60_lhmedium_nod0 << std::endl;
  std::cout << "pass_e140_lhloose_nod0=" << pass_e140_lhloose_nod0 << std::endl;
  //if(!pass_e26_lhtight_nod0_ivarloose && !pass_e60_lhmedium_nod0 && !pass_e140_lhloose_nod0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(1);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(1);
  h->h1D("nInitialTrackJet_cut1", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut1", "", suffix)->Fill(evt.tjet().size(), weight);
  //std::cout << "Pass Trigger selection" << std::endl;

  //Number of electron
  int nElectron = 0,nElectron25=0;
  for(int ie=0;ie<evt.electron().size();ie++){
    //std::cout << "electron " << ie << std::endl;
    TLorentzVector el = evt.electron()[ie].mom();
    if(ie==0) h->h1D("electronPt", "", suffix)->Fill(el.Perp()*1e-3, weight);
    if(el.Perp()*1e-3 > 30) nElectron++;
    if(el.Perp()*1e-3 > 25) nElectron25++;
    bool isTight = evt.electron()[ie].isTightPP();
    //std::cout << "pt/isTight=" << "/" << el.Perp()*1e-3 << "/" << isTight << std::endl;
  }
  h->h1D("nElectron", "", suffix)->Fill(nElectron, weight);
  if(nElectron==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(2);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(2);
  //std::cout << "Pass #Electron == 0" << std::endl;
  h->h1D("nInitialTrackJet_cut2", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut2", "", suffix)->Fill(evt.tjet().size(), weight);

  if(nElectron25>1) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(3);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(3);
  h->h1D("nInitialTrackJet_cut3", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut3", "", suffix)->Fill(evt.tjet().size(), weight);

  //Number of muon
  int nMuon = 0;
  for(int im=0;im<evt.muon().size();im++){
    TLorentzVector mu = evt.muon()[im].mom();
    bool isTight = evt.muon()[im].isTight();
    h->h1D("muonPt", "", suffix)->Fill(mu.Perp()*1e-3, weight);
    if(mu.Perp()*1e-3 < 25) continue;
    nMuon++;
  }
  //Number of muon == 0
  if(nMuon>0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(4);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(4);
  h->h1D("nInitialTrackJet_cut4", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut4", "", suffix)->Fill(evt.tjet().size(), weight);

  //Trigger matching
  trig1 = evt.electron()[0].HLT_e26_lhtight_nod0_ivarloose();
  trig2 = evt.electron()[0].HLT_e60_lhmedium_nod0();
  trig3 = evt.electron()[0].HLT_e140_lhloose_nod0();
  if(!trig1 && !trig2 && !trig3) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(5);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(5);
  //std::cout << "Pass Trigger Matching" << std::endl;
  h->h1D("nInitialTrackJet_cut5", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut5", "", suffix)->Fill(evt.tjet().size(), weight);

  //Jet Cleaning
  //std::cout << "jet clean is " << evt.passes("jet_clean") << std::endl;
  //if(!(evt.passes("jet_clean"))) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(6);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(6);
  //std::cout << "Pass Jet Clean" << std::endl;
  h->h1D("nInitialTrackJet_cut6", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut6", "", suffix)->Fill(evt.tjet().size(), weight);

  //MET & MET+MWT 
  TLorentzVector l = evt.electron()[0].mom();
  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3;
  float MET = evt.met().Perp()*1e-3;
  float MET_plus_mWt = mWt + MET;
  //std::cout << "MET/MWT=" << MET << "/" << mWt << std::endl;
  h->h1D("MET", "", suffix)->Fill(MET, weight);
  h->h1D("MET_plus_MWT", "", suffix)->Fill(MET_plus_mWt, weight);
  h->h1D("MWT", "", suffix)->Fill(mWt, weight);
  if(MET<20) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(7);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(7);
  h->h1D("nInitialTrackJet_cut7", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut7", "", suffix)->Fill(evt.tjet().size(), weight);
  //std::cout << "Pass MET" << std::endl;
  if(MET_plus_mWt<60) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(8);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(8);
  //std::cout << "Pass MET+mWt" << std::endl;
  h->h1D("nInitialTrackJet_cut8", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut8", "", suffix)->Fill(evt.tjet().size(), weight);

  //Number of small jets
  //std::cout << "Number of small jet is " << evt.jet().size() << std::endl;
  int nJet=0;
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    float jet_pt = evt.jet()[jidx].mom().Perp()*1e-3;
    float jet_eta = evt.jet()[jidx].mom().Eta();
    float jet_phi = evt.jet()[jidx].mom().Phi();
    h->h1D("small_jet_pt", "", suffix)->Fill(jet_pt, weight);
    h->h1D("small_jet_eta", "", suffix)->Fill(jet_eta, weight);
    h->h1D("small_jet_phi", "", suffix)->Fill(jet_phi, weight);
    if(jet_pt>25) nJet++;
  }
  h->h1D("nJets", "", suffix)->Fill(nJet, weight);
  if(nJet==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(9);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(9);
  //std::cout << "Pass #Small Jet" << std::endl;
  h->h1D("nInitialTrackJet_cut9", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut9", "", suffix)->Fill(evt.tjet().size(), weight);

  //Jet close to lepton
  size_t close_idx = 0;
  for (; close_idx < evt.jet().size(); ++close_idx){    
    float dphi = l.DeltaPhi(evt.jet()[close_idx].mom());
    float dy = l.Rapidity() - evt.jet()[close_idx].mom().Rapidity();
    float deltaR_tmp = std::sqrt(dy*dy + dphi*dphi);
    if(evt.jet()[close_idx].mom().Pt()*1e-3 > 25 && deltaR_tmp < 1.5) break;
    //if (evt.jet()[close_idx].closeToLepton()) break;
  }//for     
  if(close_idx==evt.jet().size()) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(10);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(10);
  std::cout << "Pass Jet Close to Lepton" << std::endl;
  h->h1D("nInitialTrackJet_cut10", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut10", "", suffix)->Fill(evt.tjet().size(), weight);

  //Number of large jet
  int nLargeJet=0;
  int ljetid=-1;
  std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
  for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
    const TLorentzVector &lj = evt.largeJet()[ljid].mom();
    bool good = evt.largeJet()[ljid].good();//default
    std::cout << "large jet pt/eta/good=" << lj.Perp()*1e-3 << "/" << lj.Eta() << "/" << good << std::endl;
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0 && good) {
      ljetid=ljid;
      nLargeJet++;
    }
    h->h1D("largeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    h->h1D("largeJetEta", "", suffix)->Fill(lj.Eta(), weight);
    if(ljid==0) {
      h->h1D("leadingLargeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      h->h1D("leadingLargeJetEta", "", suffix)->Fill(lj.Eta(), weight);
    }
  }
  h->h1D("nLargeJets", "", suffix)->Fill(nLargeJet, weight);
  if(nLargeJet==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(11);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(11);
  h->h1D("nInitialTrackJet_cut11", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut11", "", suffix)->Fill(evt.tjet().size(), weight);
  //std::cout << "Pass #Large Jet" << std::endl;

  //Angular cut
  int nGoodJets=0;
  TLorentzVector sj = evt.jet()[close_idx].mom();
  for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
    const TLorentzVector &lj = evt.largeJet()[ljid].mom();
    bool good = evt.largeJet()[ljid].good();//default
    if(lj.Perp()*1e-3 < 300 || fabs(lj.Eta())>2.0 || !good) continue;
    float delta_phi_large_jet_lepton = l.DeltaPhi(lj);
    float delta_r_large_jet_small_jet = lj.DeltaR(sj);
    //std::cout << "large jet pt/eta/good=" << lj.Perp()*1e-3 << "/" << lj.Eta() << "/" << good << std::endl;
    //std::cout << "dphi_lj_lep/dr_lj_sj=" << delta_phi_large_jet_lepton << "/" << delta_r_large_jet_small_jet << std::endl;
    if(fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) nGoodJets++;
  }
  //std::cout << "nGoodJets=" << nGoodJets << std::endl;
  if(nGoodJets==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(12);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(12);
  //std::cout << "Pass Angular Cut" << std::endl;
  h->h1D("nInitialTrackJet_cut12", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut12", "", suffix)->Fill(evt.tjet().size(), weight);

  //Track B-tag
  int nTrackBtagJets=0,nTrackFlatBtagJets=0;
  float temp_thr=-1;
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    float pt = evt.tjet()[jidx].mom().Pt()*1e-3;
    float eta = evt.tjet()[jidx].mom().Eta();
    float phi = evt.tjet()[jidx].mom().Phi();
    if(pt<300) temp_thr = 0.6455;
    else temp_thr = -0.3;
    //std::cout << "tjet pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl; 
    //std::cout << "b-tag value = " << evt.tjet()[jidx].mv2c10() << std::endl;
    bool is_trackBtag = evt.tjet()[jidx].btag_mv2c10_70_trk();
    bool is_trackFlatBtag = evt.tjet()[jidx].mv2c10() > temp_thr;
    if (is_trackBtag && evt.tjet()[jidx].mom().Pt() > 10e3 && 
        fabs(evt.tjet()[jidx].mom().Eta()) < 2.5 && evt.tjet()[jidx].numConstituents() >= 2) nTrackBtagJets++; 
    if (is_trackFlatBtag && evt.tjet()[jidx].mom().Pt() > 10e3 && 
        fabs(evt.tjet()[jidx].mom().Eta()) < 2.5 && evt.tjet()[jidx].numConstituents() >= 2) nTrackFlatBtagJets++; 
  }
  //std::cout << "Number of track b-tagged jets is " << nTrackBtagJets << std::endl;
  h->h1D("nTrkBtagJets", "", suffix)->Fill(nTrackBtagJets, weight);  

  //Calo B-tag
  int nCaloBtagJets=0;
  for (int jid=0; jid < evt.jet().size(); ++jid){    
    bool is_caloBtag = evt.jet()[jid].btag_mv2c20_70();
    if(is_caloBtag) nCaloBtagJets++;
  }
  //std::cout << "Number of calo b-tagged jets is " << nCaloBtagJets << std::endl;
  h->h1D("nCaloBtagJets", "", suffix)->Fill(nCaloBtagJets, weight);  

  if(nTrackBtagJets==0) return;
  //if(nTrackFlatBtagJets==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(13);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(13);
  //std::cout << "Pass B-tag" << std::endl;
  h->h1D("nInitialTrackJet_cut13", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut13", "", suffix)->Fill(evt.tjet().size(), weight);

  std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), true);   
  TLorentzVector nu(0,0,0,0);
  if (vec_nu.size() > 0) {
    nu = *(vec_nu[0]);
    for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
    vec_nu.clear();
  }
  const TLorentzVector &ljt = evt.largeJet()[ljetid].mom();
  float mtt = (ljt+sj+nu+l).M();
  float mlt = (sj+nu+l).M();
  h->h1D("mttAll", "", suffix)->Fill(mtt*1e-3, weight);
  h->h1D(Form("mtt%s",channel.c_str()), "", suffix)->Fill(mtt*1e-3, weight);

}
