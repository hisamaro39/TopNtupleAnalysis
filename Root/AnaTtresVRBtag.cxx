/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresVRBtag.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"


AnaTtresVRBtag::AnaTtresVRBtag(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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

    //double varBin8[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500,3000};
    double varBin8[] = {0, 20,40,60,80,100,120,140,160,180,200,250, 300, 350, 400, 450, 500};
    int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;
    double varBin9[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000};
    //double varBin9[] = {0, 20,40,60,80,100,120,140,160,180,200,250, 300, 350, 400, 450, 500};
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
    m_hSvc.create1DVar("ltop_pt_tjet",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_tjet",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_tjet_pt",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_vrtjet",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_vrtjet",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_vrtjet_pt",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_tjet_passTrackBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_tjet_passTrackBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_tjet_pt_passTrackBtag",    "; Pt of b from leptonic top matched track jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_vrtjet_passVRTrackBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_vrtjet_passVRTrackBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_vrtjet_pt_passVRTrackBtag",    "; Pt of b from leptonic top matched track jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_jet",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_jet",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_jet_pt",    "; Pt of b from hadronic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_jet_passCaloBtag",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_jet_passCaloBtag",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_jet_pt_passCaloBtag",    "; Pt of b from hadronic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_tjet",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_tjet",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_tjet_pt",    "; Pt of b from hadronic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_vrtjet",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_vrtjet",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_vrtjet_pt",    "; Pt of b from hadronic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_tjet_passTrackBtag",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_tjet_passTrackBtag",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_tjet_pt_passTrackBtag",    "; Pt of b from hadronic top matched track jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("htop_pt_vrtjet_passVRTrackBtag",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_vrtjet_passVRTrackBtag",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_vrtjet_pt_passVRTrackBtag",    "; Pt of b from hadronic top matched track jet [GeV]; Events", varBinN9, varBin9);

    m_hSvc.create1D("bmatchedCaloBtagBDT", "; b matched calo b-tag BDT ; Events", 100, -1, 1);
    m_hSvc.create1D("bmatchedTrackBtagBDT", "; b matched track b-tag BDT ; Events", 100, -1, 1);
    m_hSvc.create1D("bmatchedVRTrackBtagBDT", "; b matched VR track b-tag BDT ; Events", 100, -1, 1);
    m_hSvc.create2D("bmatchedTrackPtVsBtagBDT", "; b matched track p_{T} [GeV]; b matched track b-tag BDT", 100,0,500, 100, -1, 1);
    m_hSvc.create2D("bmatchedVRTrackPtVsBtagBDT", "; b matched VR track p_{T} [GeV]; b matched VR track b-tag BDT", 100,0,500, 100, -1, 1);
    m_hSvc.create1D("TrackBtagBDT", "; track b-tag BDT ; Events", 100, -1, 1);
    m_hSvc.create1D("CaloBtagBDT", "; calo b-tag BDT ; Events", 100, -1, 1);
    m_hSvc.create1D("VRTrackBtagBDT", "; VR track b-tag BDT ; Events", 100, -1, 1);

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
    m_hSvc.create1DVar("drVRTjetHB", "; #DeltaR_{vrtj,hb} ; Events", 100,0,1);
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

    m_hSvc.create1DVar("lbMatchedTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("hbMatchedTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("lbMatchedVRTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("hbMatchedVRTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("lbMatchedTrackJetPtPassTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("hbMatchedTrackJetPtPassTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("lbMatchedVRTrackJetPtPassVRTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("hbMatchedVRTrackJetPtPassVRTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN9, varBin9);
    for(int i=0;i<varBinN9;i++){
      m_hSvc.create1D(Form("lbMatchedTrackJetMv2c10_pt%d_%d",varBin9[i],varBin9[i+1]), "; track jet Mv2c10; Events", 20,-1,1);
      m_hSvc.create1D(Form("hbMatchedTrackJetMv2c10_pt%d_%d",varBin9[i],varBin9[i+1]), "; track jet Mv2c10; Events", 20,-1,1);
    }
  }

AnaTtresVRBtag::~AnaTtresVRBtag() {
}

void AnaTtresVRBtag::run(const Event &evt, double weight, const std::string &s, int is2016run) {
  //std::cout << "AnaTtresVRBtag::run" << std::endl;

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
  else if(containMu) channel="SemiLepMu";
  else if(containEl) channel="SemiLepEl";
  else std::cout << "Unknow channel" << std::endl;
  //std::cout << "channel is " << channel << std::endl;

  TLorentzVector b_from_htop,b_from_ltop,htop,ltop;
  if(channel=="SemiLepMu" || channel=="SemiLepEl"){
    if(Wdecay2_from_t_pdgId==-13 || Wdecay2_from_t_pdgId==-11){//leptonic top
      b_from_ltop = b_from_t;
      ltop = t;
      b_from_htop = b_from_tbar;
      htop = tbar;
    }
    else {//hadronic top
      b_from_ltop = b_from_tbar;
      ltop = tbar;
      b_from_htop = b_from_t;
      htop = t;
    }
  }
  //std::cout << "pt of ltop/htop=" << ltop.Perp() << "/" << htop.Perp() << std::endl;

  HistogramService *h = &m_hSvc;

  std::string suffix = s;

  h->h1D("cutFlowAll", "", suffix)->Fill(0);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(0);

  if(!evt.passes("bmujets_VRbtag") && !evt.passes("bejets_VRbtag"))return;
  h->h1D("cutFlowAll", "", suffix)->Fill(1);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(1);

   float min_dr_hb_tjet=100,min_dr_lb_tjet=100;
   int hb_id_tjet=-1,lb_id_tjet=-1;
   for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
     float pt = evt.tjet()[jidx].mom().Perp()*1e-3;
     float eta = evt.tjet()[jidx].mom().Eta();
     float phi = evt.tjet()[jidx].mom().Phi();
     float dr_hb = TuDoAtlas::calc_delta_r(eta,b_from_htop.Eta(),phi,b_from_htop.Phi());
     float dr_lb = TuDoAtlas::calc_delta_r(eta,b_from_ltop.Eta(),phi,b_from_ltop.Phi());
     if(dr_lb<min_dr_lb_tjet && dr_lb<0.2){
       min_dr_lb_tjet = dr_lb;
       lb_id_tjet = jidx;
     }
     if(dr_hb<min_dr_hb_tjet && dr_hb<0.2){
       min_dr_hb_tjet = dr_hb;
       hb_id_tjet = jidx;
     }
   }
   //std::cout << "tjet id lb/hb=" << lb_id_tjet << "/" << hb_id_tjet << std::endl;

   //std::cout << "Number of VR track jet is " << evt.vrtjet().size() << std::endl;
   int numVRTrackJets=0;
   float min_dr_hb_vrtjet=100,min_dr_lb_vrtjet=100;
   int hb_id_vrtjet=-1,lb_id_vrtjet=-1;
   for (size_t jidx = 0; jidx < evt.vrtjet().size(); ++jidx){
     float pt = evt.vrtjet()[jidx].mom().Perp()*1e-3;
     float eta = evt.vrtjet()[jidx].mom().Eta();
     float phi = evt.vrtjet()[jidx].mom().Phi();
     int numConstituents = evt.vrtjet()[jidx].numConstituents();
     if(pt<10 || fabs(eta)>2.5 || numConstituents<2) continue;//object selection for track jets
     //OR with lepton
     bool or_muon=false,or_electron=false;
     for(int im=0;im<evt.muon().size();im++){
       TLorentzVector mu = evt.muon()[im].mom();
       float dr = TuDoAtlas::calc_delta_r(eta,mu.Eta(),phi,mu.Phi());
       if(dr<0.2) or_muon=true;
     }
     for(int ie=0;ie<evt.electron().size();ie++){
       TLorentzVector el = evt.electron()[ie].mom();
       float dr = TuDoAtlas::calc_delta_r(eta,el.Eta(),phi,el.Phi());
       if(dr<0.2) or_electron=true;
     }
     if(or_muon || or_electron) continue;
     float dr_hb = TuDoAtlas::calc_delta_r(eta,b_from_htop.Eta(),phi,b_from_htop.Phi());
     float dr_lb = TuDoAtlas::calc_delta_r(eta,b_from_ltop.Eta(),phi,b_from_ltop.Phi());
     if(dr_lb<min_dr_lb_vrtjet && dr_lb<0.2){
       min_dr_lb_vrtjet = dr_lb;
       lb_id_vrtjet = jidx;
     }
     if(dr_hb<min_dr_hb_vrtjet && dr_hb<0.2){
       min_dr_hb_vrtjet = dr_hb;
       hb_id_vrtjet = jidx;
     }
     numVRTrackJets++;
   }
   //std::cout << "vrtjet id lb/hb=" << lb_id_vrtjet << "/" << hb_id_vrtjet << std::endl;
   //std::cout << "Number of VR track jets finaly is " << numVRTrackJets << std::endl;

    double varBin9[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000};
    int varBinN9 = sizeof(varBin9)/sizeof(double) - 1;
   //track jet b-tag
   if(lb_id_tjet!=-1){
     h->h1D("lbMatchedTrackJetPt", "", suffix)->Fill(evt.tjet()[lb_id_tjet].mom().Perp()*1e-3, weight);
     if(evt.tjet()[lb_id_tjet].mv2c10()>0.6455) h->h1D("lbMatchedTrackJetPtPassTrackBtag", "", suffix)->Fill(evt.tjet()[lb_id_tjet].mom().Perp()*1e-3, weight);
     for(int i=0;i<varBinN9;i++){
       if(evt.tjet()[lb_id_tjet].mom().Perp()*1e-3>varBin9[i] && evt.tjet()[lb_id_tjet].mom().Perp()*1e-3<varBin9[i+1])
         h->h1D(Form("lbMatchedTrackJetMv2c10_pt%d_%d",varBin9[i],varBin9[i+1]), "", suffix)->Fill(evt.tjet()[lb_id_tjet].mv2c10(), weight);
     }
   }
   if(hb_id_tjet!=-1){
     h->h1D("hbMatchedTrackJetPt", "", suffix)->Fill(evt.tjet()[hb_id_tjet].mom().Perp()*1e-3, weight);
     if(evt.tjet()[hb_id_tjet].mv2c10()>0.6455) h->h1D("hbMatchedTrackJetPtPassTrackBtag", "", suffix)->Fill(evt.tjet()[hb_id_tjet].mom().Perp()*1e-3, weight);
     for(int i=0;i<varBinN9;i++){
       if(evt.tjet()[hb_id_tjet].mom().Perp()*1e-3>varBin9[i] && evt.tjet()[hb_id_tjet].mom().Perp()*1e-3<varBin9[i+1])
         h->h1D(Form("hbMatchedTrackJetMv2c10_pt%d_%d",varBin9[i],varBin9[i+1]), "", suffix)->Fill(evt.tjet()[hb_id_tjet].mv2c10(), weight);
     }
   }

   //VR track jet b-tag
   if(lb_id_vrtjet!=-1){
     h->h1D("lbMatchedVRTrackJetPt", "", suffix)->Fill(evt.vrtjet()[lb_id_vrtjet].mom().Perp()*1e-3, weight);
     if(evt.vrtjet()[lb_id_vrtjet].mv2c10()>0.8492) h->h1D("lbMatchedVRTrackJetPtPassVRTrackBtag", "", suffix)->Fill(evt.vrtjet()[lb_id_vrtjet].mom().Perp()*1e-3, weight);
     for(int i=0;i<varBinN9;i++){
       if(evt.vrtjet()[lb_id_vrtjet].mom().Perp()*1e-3>varBin9[i] && evt.vrtjet()[lb_id_vrtjet].mom().Perp()*1e-3<varBin9[i+1])
         h->h1D(Form("lbMatchedVRTrackJetMv2c10_pt%d_%d",varBin9[i],varBin9[i+1]), "", suffix)->Fill(evt.vrtjet()[lb_id_vrtjet].mv2c10(), weight);
     }
   }
   if(hb_id_vrtjet!=-1){
     h->h1D("hbMatchedVRTrackJetPt", "", suffix)->Fill(evt.vrtjet()[hb_id_vrtjet].mom().Perp()*1e-3, weight);
     if(evt.vrtjet()[hb_id_vrtjet].mv2c10()>0.8492) h->h1D("hbMatchedVRTrackJetPtPassVRTrackBtag", "", suffix)->Fill(evt.vrtjet()[hb_id_vrtjet].mom().Perp()*1e-3, weight);
     for(int i=0;i<varBinN9;i++){
       if(evt.vrtjet()[hb_id_vrtjet].mom().Perp()*1e-3>varBin9[i] && evt.vrtjet()[hb_id_vrtjet].mom().Perp()*1e-3<varBin9[i+1])
         h->h1D(Form("hbMatchedVRTrackJetMv2c10_pt%d_%d",varBin9[i],varBin9[i+1]), "", suffix)->Fill(evt.vrtjet()[hb_id_vrtjet].mv2c10(), weight);
     }
   }

  /*
  //truth b & jet matching
  int bt_sj_idx=-1, btbar_sj_idx=-1;
  float min_dr_sj_bt = 999., min_dr_sj_btbar = 999.;
  //std::cout << "Nubmer of small jet is " << evt.jet().size() << std::endl;
  for (size_t bidx = 0; bidx < evt.jet().size(); ++bidx){
    float sj_eta = evt.jet()[bidx].mom().Eta();
    float sj_phi = evt.jet()[bidx].mom().Phi();
    float drt = TuDoAtlas::calc_delta_r(sj_eta,b_from_t_eta,sj_phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(sj_eta,b_from_tbar_eta,sj_phi,b_from_tbar_phi);
    float pt = evt.jet()[bidx].mom().Pt();
    float eta = evt.jet()[bidx].mom().Eta();
    float phi = evt.jet()[bidx].mom().Phi();
    int isBtag = evt.jet()[bidx].btag_mv2c20_70();
    //std::cout << "small jet pt/eta/phi/btag=" << pt << "/" << eta << "/" << phi << "/" << isBtag << std::endl;
    //std::cout << "btag value=" << evt.jet()[bidx].mv2c20() << std::endl;
    if(!isLeptonicT) h->h1D("drJetHB", "", suffix)->Fill(drt, weight);  
    if(!isLeptonicTbar) h->h1D("drJetHB", "", suffix)->Fill(drtbar, weight);  
    if(pt*1e-3>25) h->h1D("CaloBtagBDT", "", suffix)->Fill(evt.jet()[bidx].mv2c20());
    if(drt<0.4 && drt<min_dr_sj_bt) {
      min_dr_sj_bt = drt;
      bt_sj_idx=bidx;
    }
    if(drtbar<0.4 && drtbar<min_dr_sj_btbar) {
      min_dr_sj_btbar = drtbar;
      btbar_sj_idx=bidx;
    }
  }

  //truth b & tjet matching
  std::cout << "Number of track jet is " << evt.tjet().size() << std::endl;
  std::vector<int> bt_tj_idx, btbar_tj_idx;
  bt_tj_idx.clear();btbar_tj_idx.clear();
  float min_dr_tj_bt = 999., min_dr_tj_btbar = 999.;
  int min_bt_tj_idx=-1,min_btbar_tj_idx=-1;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
    float pt = evt.tjet()[bidx].mom().Pt();
    float eta = evt.tjet()[bidx].mom().Eta();
    float phi = evt.tjet()[bidx].mom().Phi();
    int numConstituents = evt.tjet()[bidx].numConstituents();
    if(pt*1e-3<10 || fabs(eta)>2.5 || numConstituents<2) continue;
    float drt = TuDoAtlas::calc_delta_r(eta,b_from_t_eta,phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(eta,b_from_tbar_eta,phi,b_from_tbar_phi);
    if(!isLeptonicT) h->h1D("drTjetHB", "", suffix)->Fill(drt, weight);  
    if(!isLeptonicTbar) h->h1D("drTjetHB", "", suffix)->Fill(drtbar, weight);  
    if(pt*1e-3>25) h->h1D("TrackBtagBDT", "", suffix)->Fill(evt.tjet()[bidx].mv2c10());
    if(drt<0.4) {
      bt_tj_idx.push_back(bidx);
      if(drt<min_dr_tj_bt){
        min_dr_tj_bt=drt;
        min_bt_tj_idx=bidx;
      }
    }
    if(drtbar<0.4) {
      btbar_tj_idx.push_back(bidx);
      if(drtbar<min_dr_tj_btbar){
        min_dr_tj_btbar=drtbar;
        min_btbar_tj_idx=bidx;
      }
    }
  }

  //truth b & vrtjet matching
  std::cout << "Number of VR track jet is " << evt.vrtjet().size() << std::endl;
  std::vector<int> bt_vrtj_idx, btbar_vrtj_idx;
  bt_vrtj_idx.clear();btbar_vrtj_idx.clear();
  float min_dr_vrtj_bt = 999., min_dr_vrtj_btbar = 999.;
  int num_close_tb=0,num_close_tbarb=0;
  int min_bt_vrtj_idx=-1,min_btbar_vrtj_idx=-1;
  int numVRTrackJets=0;
  for (size_t bidx = 0; bidx < evt.vrtjet().size(); ++bidx){
    float pt = evt.vrtjet()[bidx].mom().Pt();
    float eta = evt.vrtjet()[bidx].mom().Eta();
    int numConstituents = evt.vrtjet()[bidx].numConstituents();
    if(pt*1e-3<10 || fabs(eta)>2.5 || numConstituents<2) continue;
    float phi = evt.vrtjet()[bidx].mom().Phi();
    float drt = TuDoAtlas::calc_delta_r(eta,b_from_t_eta,phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(eta,b_from_tbar_eta,phi,b_from_tbar_phi);
    if(!isLeptonicT) h->h1D("drVRTjetHB", "", suffix)->Fill(drt, weight);  
    if(!isLeptonicTbar) h->h1D("drVRTjetHB", "", suffix)->Fill(drtbar, weight);  
    if(pt*1e-3>25) h->h1D("VRTrackBtagBDT", "", suffix)->Fill(evt.vrtjet()[bidx].mv2c10());
    if(drt<0.4) {
      num_close_tb++;
      bt_vrtj_idx.push_back(bidx);
      if(drt<min_dr_vrtj_bt){
        min_dr_vrtj_bt=drt;
        min_bt_vrtj_idx=bidx;
      }
    }
    if(drtbar<0.4) {
      num_close_tbarb++;
      btbar_vrtj_idx.push_back(bidx);
      if(drtbar<min_dr_vrtj_btbar){
        min_dr_vrtj_btbar=drtbar;
        min_btbar_vrtj_idx=bidx;
      }
    }
    numVRTrackJets++;
  }
  //std::cout << "num close to tb/tbarb=" << num_close_tb << "/" << num_close_tbarb << std::endl;
  std::cout << "Number of good VR track jet is " << numVRTrackJets << std::endl;

  //calo jet b-tag efficiency
  if(bt_sj_idx!=-1){
    h->h1D("bmatchedCaloBtagBDT", "", suffix)->Fill(evt.jet()[bt_sj_idx].mv2c20());
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_jet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_ltop_pt_jet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_ltop_matched_jet_pt", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[bt_sj_idx].btag_mv2c20_70()){
        //std::cout << "ltop jet  btagged!!" << std::endl;
        h->h1D("ltop_pt_jet_passCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
      //else std::cout << "ltop jet not btagged!!" << std::endl;
    }
    else{//hadronic top
      h->h1D("htop_pt_jet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_jet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_htop_matched_jet_pt", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[bt_sj_idx].btag_mv2c20_70()){
        //std::cout << "htop jet btagged!!" << std::endl;
        h->h1D("htop_pt_jet_passCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
      //else std::cout << "htop jet not btagged!!" << std::endl;
    }
  }
  if(btbar_sj_idx!=-1){
    h->h1D("bmatchedCaloBtagBDT", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mv2c20());
    if(isLeptonicTbar){//leptonic antitop
      h->h1D("ltop_pt_jet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_ltop_pt_jet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_ltop_matched_jet_pt", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[btbar_sj_idx].btag_mv2c20_70()){
        //std::cout << "ltop jet btagged!!" << std::endl;
        h->h1D("ltop_pt_jet_passCaloBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_ltop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      }
      //else std::cout << "ltop jet not btagged!!" << std::endl;
    }
    else{//hadronic antitop
      h->h1D("htop_pt_jet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_htop_pt_jet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_htop_matched_jet_pt", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[btbar_sj_idx].btag_mv2c20_70()){
        //std::cout << "htop jet btagged!!" << std::endl;
        h->h1D("htop_pt_jet_passCaloBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_htop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      }
      //else std::cout << "htop jet not btagged!!" << std::endl;
    }
  }

  //Track jet b-tag efficiency (include)
  if(bt_tj_idx.size()){//top
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_tjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_ltop_pt_tjet", "", suffix)->Fill(b_from_t_pt);
      bool isBtag=false;
      for(int b=0;b<bt_tj_idx.size();b++)
        if(evt.tjet()[bt_tj_idx.at(b)].mv2c10()>0.6455) isBtag=true;
      if(isBtag) {
        h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_t_pt);
      }
    }
    else {//hadronic top
      h->h1D("htop_pt_tjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_t_pt);
      bool isBtag=false;
      for(int b=0;b<bt_tj_idx.size();b++)
        if(evt.tjet()[bt_tj_idx.at(b)].mv2c10()>0.6455) isBtag=true;
      if(isBtag) {
        h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_t_pt);
      }
    }
  }
  if(btbar_tj_idx.size()){//anti top
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_tjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_ltop_pt_tjet", "", suffix)->Fill(b_from_tbar_pt);
      bool isBtag=false;
      for(int b=0;b<btbar_tj_idx.size();b++)
        if(evt.tjet()[btbar_tj_idx.at(b)].mv2c10()>0.6455) isBtag=true;
      if(isBtag) {
        h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
      }
    }
    else {//hadronic top
      h->h1D("htop_pt_tjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_tbar_pt);
      bool isBtag=false;
      for(int b=0;b<btbar_tj_idx.size();b++)
        if(evt.tjet()[btbar_tj_idx.at(b)].mv2c10()>0.6455) isBtag=true;
      if(isBtag) {
        h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
      }
    }
  }

  //Track jet b-tag efficiency (match)
  if(min_bt_tj_idx!=-1){//top
    float pt = evt.tjet()[min_bt_tj_idx].mom().Perp()*1e-3;
    bool isBtag = evt.tjet()[min_bt_tj_idx].mv2c10()>0.6455;
    h->h1D("bmatchedTrackBtagBDT", "", suffix)->Fill(evt.tjet()[min_bt_tj_idx].mv2c10());
    h->h2D("bmatchedTrackPtVsBtagBDT", "", suffix)->Fill(pt,evt.tjet()[min_bt_tj_idx].mv2c10());
    if(isLeptonicT){//leptonic top
      h->h1D("b_from_ltop_matched_tjet_pt", "", suffix)->Fill(pt);
      if(isBtag) h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(pt);
    }
    else {//hadronic top
      h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(pt);
      if(isBtag) h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(pt);
    }
  }
  if(min_btbar_tj_idx!=-1){//top
    float pt = evt.tjet()[min_btbar_tj_idx].mom().Perp()*1e-3;
    bool isBtag = evt.tjet()[min_btbar_tj_idx].mv2c10()>0.6455;
    h->h1D("bmatchedTrackBtagBDT", "", suffix)->Fill(evt.tjet()[min_btbar_tj_idx].mv2c10());
    h->h2D("bmatchedTrackPtVsBtagBDT", "", suffix)->Fill(pt,evt.tjet()[min_btbar_tj_idx].mv2c10());
    if(isLeptonicT){//leptonic top
      h->h1D("b_from_ltop_matched_tjet_pt", "", suffix)->Fill(pt);
      if(isBtag) h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(pt);
    }
    else {//hadronic top
      h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(pt);
      if(isBtag) h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(pt);
    }
  }

  //VR track jet b-tag efficiency (include)
  if(bt_vrtj_idx.size()){//top
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_vrtjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_ltop_pt_vrtjet", "", suffix)->Fill(b_from_t_pt);
      bool isBtag=false;
      for(int b=0;b<bt_vrtj_idx.size();b++)
        if(evt.vrtjet()[bt_vrtj_idx.at(b)].mv2c10()>0.8492) isBtag=true;
      if(isBtag) {
        h->h1D("ltop_pt_vrtjet_passVRTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_vrtjet_passVRTrackBtag", "", suffix)->Fill(b_from_t_pt);
      }
    }
    else {//hadronic top
      h->h1D("htop_pt_vrtjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_vrtjet", "", suffix)->Fill(b_from_t_pt);
      bool isBtag=false;
      for(int b=0;b<bt_vrtj_idx.size();b++)
        if(evt.vrtjet()[bt_vrtj_idx.at(b)].mv2c10()>0.8492) isBtag=true;
      if(isBtag) {
        h->h1D("htop_pt_vrtjet_passVRTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_vrtjet_passVRTrackBtag", "", suffix)->Fill(b_from_t_pt);
      }
    }
  }
  if(btbar_vrtj_idx.size()){//anti top
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_vrtjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_ltop_pt_vrtjet", "", suffix)->Fill(b_from_tbar_pt);
      bool isBtag=false;
      for(int b=0;b<btbar_vrtj_idx.size();b++)
        if(evt.vrtjet()[btbar_vrtj_idx.at(b)].mv2c10()>0.8492) isBtag=true;
      if(isBtag) {
        h->h1D("ltop_pt_vrtjet_passVRTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_vrtjet_passVRTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
      }
    }
    else {//hadronic top
      h->h1D("htop_pt_vrtjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_htop_pt_vrtjet", "", suffix)->Fill(b_from_tbar_pt);
      bool isBtag=false;
      for(int b=0;b<btbar_vrtj_idx.size();b++)
        if(evt.vrtjet()[btbar_vrtj_idx.at(b)].mv2c10()>0.8492) isBtag=true;
      if(isBtag) {
        h->h1D("htop_pt_vrtjet_passVRTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_vrtjet_passVRTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
      }
    }
  }

  //VR track jet b-tag efficiency (match)
  if(min_bt_vrtj_idx!=-1){//top
    float pt = evt.vrtjet()[min_bt_vrtj_idx].mom().Perp()*1e-3;
    bool isBtag = evt.vrtjet()[min_bt_vrtj_idx].mv2c10()>0.8492;
    h->h1D("bmatchedVRTrackBtagBDT", "", suffix)->Fill(evt.vrtjet()[min_bt_vrtj_idx].mv2c10());
    h->h2D("bmatchedVRTrackPtVsBtagBDT", "", suffix)->Fill(pt,evt.vrtjet()[min_bt_vrtj_idx].mv2c10());
    if(isLeptonicT){//leptonic top
      h->h1D("b_from_ltop_matched_vrtjet_pt", "", suffix)->Fill(pt);
      if(isBtag) h->h1D("b_from_ltop_matched_vrtjet_pt_passVRTrackBtag", "", suffix)->Fill(pt);
    }
    else {//hadronic top
      h->h1D("b_from_htop_matched_vrtjet_pt", "", suffix)->Fill(pt);
      if(isBtag) h->h1D("b_from_htop_matched_vrtjet_pt_passVRTrackBtag", "", suffix)->Fill(pt);
    }
  }
  if(min_btbar_vrtj_idx!=-1){//anti top
    float pt = evt.vrtjet()[min_btbar_vrtj_idx].mom().Perp()*1e-3;
    bool isBtag = evt.vrtjet()[min_btbar_vrtj_idx].mv2c10()>0.8492;
    h->h1D("bmatchedVRTrackBtagBDT", "", suffix)->Fill(evt.vrtjet()[min_btbar_vrtj_idx].mv2c10());
    h->h2D("bmatchedVRTrackPtVsBtagBDT", "", suffix)->Fill(pt,evt.vrtjet()[min_btbar_vrtj_idx].mv2c10());
    if(isLeptonicT){//leptonic top
      h->h1D("b_from_ltop_matched_vrtjet_pt", "", suffix)->Fill(pt);
      if(isBtag) h->h1D("b_from_ltop_matched_vrtjet_pt_passVRTrackBtag", "", suffix)->Fill(pt);
    }
    else {//hadronic top
      h->h1D("b_from_htop_matched_vrtjet_pt", "", suffix)->Fill(pt);
      if(isBtag) h->h1D("b_from_htop_matched_vrtjet_pt_passVRTrackBtag", "", suffix)->Fill(pt);
    }
  }
*/
}
