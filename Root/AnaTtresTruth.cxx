/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresTruth.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"


AnaTtresTruth::AnaTtresTruth(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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


    m_hSvc.create1DVar("ltop_pt_jet",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_jet",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_jet_pt",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_jet_passCaloBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_jet_passCaloBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_jet_pt_passCaloBtag",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_jet_passFlatCaloBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_jet_passFlatCaloBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_jet_pt_passFlatCaloBtag",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
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
    m_hSvc.create1DVar("htop_pt_jet_passFlatCaloBtag",    "; Pt of truth hadronic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_htop_pt_jet_passFlatCaloBtag",    "; Pt of truth b from hadronic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_htop_matched_jet_pt_passFlatCaloBtag",    "; Pt of b from hadronic top matched jet [GeV]; Events", varBinN9, varBin9);
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
    m_hSvc.create1D("Mv2c10_trackb_lowpt", "; track b-tag BDT ; Events",40, -1, 1);
    m_hSvc.create1D("Mv2c10_trackb_highpt", "; track b-tag BDT ; Events", 40, -1, 1);

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
    std::string bdt_thr[3] = {"loose","middle","tight"};
    for(int i=0;i<3;i++){
      m_hSvc.create1DVar(Form("allLargeJetPt_passBDT_%s",bdt_thr[i].c_str()), "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    }
    for(int lpt=0;lpt<4;lpt++){
      m_hSvc.create1D(Form("BDT_TOPtag_allljet_pt%d_%d",lpt*500,500+lpt*500),    "; BDT TOP tag; ", 100, -1, 1);
      m_hSvc.create1D(Form("BDT_TOPtag_passljet_pt%d_%d",lpt*500,500+lpt*500),    "; BDT TOP tag; ", 100, -1, 1);
    }

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
    m_hSvc.create1D("nInitialTrackJetFinal", "; # of track jets ; Events", 15, 0, 15);  

    for (int i=0;i<14;i++){
      m_hSvc.create1D(Form("nInitialTrackJet_cut%d",i), "; # of track jets ; Events", 20, 0, 20);  
      m_hSvc.create1D(Form("nTrackJet_cut%d",i), "; # of track jets ; Events", 20, 0, 20);  
    }

    std::string channel[5] = {"FullHad","FullLep","SemiLepTau","SemiLepMu","SemiLepEl"};
    for(int ch=0;ch<5;ch++) {
      m_hSvc.create1D(Form("cutFlow%s",channel[ch].c_str()), ";Number of events; step", 20,0,20);
      m_hSvc.create1D(Form("nMuon%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1D(Form("nMuonFinal%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1D(Form("nMuonNoIso%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1D(Form("mtt%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mttNew%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mttTruthReal%s",channel[ch].c_str()), "; truth mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mttTruthImag%s",channel[ch].c_str()), "; truth mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mlt%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mltTruthReal%s",channel[ch].c_str()), "; truth mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mltTruthImag%s",channel[ch].c_str()), "; truth mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mltNew%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mht%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mhtGoodMtt%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mhtBadMtt%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mttRealSolSmall%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mttRealSolBig%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mttImagSolRotate%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mttImagSolRealPart%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mltRealSolSmall%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mltRealSolBig%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mltImagSolRotate%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mltImagSolRealPart%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("nupzRealSolSmall%s",channel[ch].c_str()), "; calculated neutrino p_{Z} [GeV]; Events", 50,0,1000);
      m_hSvc.create1D(Form("nupzRealSolBig%s",channel[ch].c_str()), "; calculated neutrino p_{Z} [GeV]; Events", 50,0,1000);
      m_hSvc.create1D(Form("nupzImagSolRotate%s",channel[ch].c_str()), "; calculated neutrino p_{Z} [GeV]; Events", 50,0,1000);
      m_hSvc.create1D(Form("nupzImagSolRealPart%s",channel[ch].c_str()), "; calculated neutrino p_{Z} [GeV]; Events", 50,0,1000);
      m_hSvc.create1D(Form("nupzTruthReal%s",channel[ch].c_str()), "; truth neutrino p_{Z} [GeV]; Events", 50,0,1000);
      m_hSvc.create1D(Form("nupzTruthImag%s",channel[ch].c_str()), "; truth neutrino p_{Z} [GeV]; Events", 50,0,1000);
      m_hSvc.create1D(Form("resMttRealSolSmall%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; Events", 100,-1,1);
      m_hSvc.create1D(Form("resMttRealSolBig%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; Events", 100,-1,1);
      m_hSvc.create1D(Form("resMttImagSolRotate%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; Events", 100,-1,1);
      m_hSvc.create1D(Form("resMttImagSolRealPart%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; Events", 100,-1,1);
      m_hSvc.create1D(Form("resMltRealSolSmall%s",channel[ch].c_str()), "; residual for mass of the leptonic top; Events", 100,-3,3);
      m_hSvc.create1D(Form("resMltRealSolBig%s",channel[ch].c_str()), "; residual for mass of the leptonic top; Events", 100,-3,3);
      m_hSvc.create1D(Form("resMltImagSolRotate%s",channel[ch].c_str()), "; residual for mass of the leptonic top; Events", 100,-3,3);
      m_hSvc.create1D(Form("resMltImagSolRealPart%s",channel[ch].c_str()), "; residual for mass of the leptonic top; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzRealSolSmall%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzRealSolBig%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzImagSolRotate%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzImagSolRealPart%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzLowRealSolSmall%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzLowRealSolBig%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzLowImagSolRotate%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzLowImagSolRealPart%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzHighRealSolSmall%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzHighRealSolBig%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzHighImagSolRotate%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNupzHighImagSolRealPart%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNuetaRealSolSmall%s",channel[ch].c_str()), "; residual for neutrino #eta; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNuetaRealSolBig%s",channel[ch].c_str()), "; residual for neutrino #eta; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNuetaImagSolRotate%s",channel[ch].c_str()), "; residual for neutrino #eta; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNuetaImagSolRealPart%s",channel[ch].c_str()), "; residual for neutrino #eta; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNuphiRealSolSmall%s",channel[ch].c_str()), "; residual for neutrino #phi; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNuphiRealSolBig%s",channel[ch].c_str()), "; residual for neutrino #phi; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNuphiImagSolRotate%s",channel[ch].c_str()), "; residual for neutrino #phi; Events", 100,-3,3);
      m_hSvc.create1D(Form("resNuphiImagSolRealPart%s",channel[ch].c_str()), "; residual for neutrino #phi; Events", 100,-3,3);
      m_hSvc.create2D(Form("resNupzImagSolRotateVsMetAlpha%s",channel[ch].c_str()), "; residual for neutrino p_{Z}; MET #alpha", 100,-3,3,100,-3,3);
      m_hSvc.create2D(Form("resMttVsTrueMtt%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; mass of true top-antitop system [GeV]", 100,-1,1,60,0,6000);
      m_hSvc.create2D(Form("resMttVsMtt%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; mass of top-antitop system [GeV]", 100,-1,1,60,0,6000);
      m_hSvc.create2D(Form("resMttVsMlt%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; mass of leptonic top [GeV]", 100,-1,1,50,0,500);
      m_hSvc.create2D(Form("resMttVsMht%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; mass of hadronic top [GeV]", 100,-1,1,50,0,500);
      m_hSvc.create2D(Form("resMttVsResMht%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; residual for mass of hadronic top", 100,-1,1,100,-1,1);
      m_hSvc.create2D(Form("mhtVsTrueMht%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; mass of true hadronic top [GeV]", 50,0,500,50,0,500);
      m_hSvc.create1D(Form("resMtt%s",channel[ch].c_str()), "; residual for mass of the top-antitop system; Events", 100,-1,1);
      m_hSvc.create1D(Form("resMlt%s",channel[ch].c_str()), "; residual for mass of the leptonic top; Events", 100,-3,3);
      m_hSvc.create1D(Form("resMht%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetPt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetPtGoodMtt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetPtBadMtt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetEta%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetEtaGoodMtt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetEtaBadMtt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetPhi%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetPhiGoodMtt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resLargeJetPhiBadMtt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resMhtGoodMtt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("resMhtBadMtt%s",channel[ch].c_str()), "; residual for mass of the hadronic top; Events", 100,-1,1);
      m_hSvc.create1D(Form("drTrueWTrueBHad%s",channel[ch].c_str()), "; #DeltaR_{true W, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drTrueWTrueBHadGoodMtt%s",channel[ch].c_str()), "; #DeltaR_{true W, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drTrueWTrueBHadBadMtt%s",channel[ch].c_str()), "; #DeltaR_{true W, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueBHad%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueBHadGoodMtt%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueBHadBadMtt%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueHtop%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true hadronic top}; Events", 50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueHtopGoodMtt%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true hadronic top}; Events", 50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueHtopBadMtt%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true hadronic top}; Events", 50,0,5);
      m_hSvc.create1D(Form("drSmallJetTrueBfromLtop%s",channel[ch].c_str()), "; #DeltaR_{small jet, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drSmallJetTrueBfromLtopGoodMtt%s",channel[ch].c_str()), "; #DeltaR_{small jet, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drSmallJetTrueBfromLtopBadMtt%s",channel[ch].c_str()), "; #DeltaR_{small jet, true b}; Events", 50,0,5);
      m_hSvc.create1D(Form("drMuonTrueMuon%s",channel[ch].c_str()), "; #DeltaR_{muon, true muon}; Events", 50,0,5);
      m_hSvc.create1D(Form("drMuonTrueMuonGoodMtt%s",channel[ch].c_str()), "; #DeltaR_{muon, true muon}; Events", 50,0,5);
      m_hSvc.create1D(Form("drMuonTrueMuonBadMtt%s",channel[ch].c_str()), "; #DeltaR_{muon, true muon}; Events", 50,0,5);
      m_hSvc.create2D(Form("mhtVsDrTrueWTrueBHad%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; #DeltaR_{true W, true b}", 50,0,500,50,0,5);
      m_hSvc.create2D(Form("mhtVsDrLargeJetTrueBHad%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; #DeltaR_{large-R jet, true b}", 50,0,500,50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueBfromLtop%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true b from ltop}; Events", 50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueBfromLtopGoodMtt%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true b from ltop}; Events", 50,0,5);
      m_hSvc.create1D(Form("drLargeJetTrueBfromLtopBadMtt%s",channel[ch].c_str()), "; #DeltaR_{large-R jet, true b from ltop}; Events", 50,0,5);
      m_hSvc.create1D(Form("muonPt_%s",channel[ch].c_str()), ";p_{T,#mu}; Events", 50,0,500);
      m_hSvc.create1D(Form("muonEta_%s",channel[ch].c_str()), ";muon #eta; Events", 25,-2.5,2.5);
      m_hSvc.create1D(Form("muonPhi_%s",channel[ch].c_str()), ";muon #phi; Events", 32,-3.2,3.2);
      m_hSvc.create1D(Form("muonPtPassMuTrig_%s",channel[ch].c_str()), ";p_{T,#mu}; Events", 50,0,500);
      m_hSvc.create1D(Form("muonEtaPassMuTrig_%s",channel[ch].c_str()), ";muon #eta; Events", 25,-2.5,2.5);
      m_hSvc.create1D(Form("muonPhiPassMuTrig_%s",channel[ch].c_str()), ";muon #phi; Events", 32,-3.2,3.2);
      m_hSvc.create1D(Form("initialMuonPt_%s",channel[ch].c_str()), ";p_{T,#mu}; Events", 50,0,500);
      m_hSvc.create1D(Form("initialMuonEta_%s",channel[ch].c_str()), ";muon #eta; Events", 25,-2.5,2.5);
      m_hSvc.create1D(Form("initialMuonPhi_%s",channel[ch].c_str()), ";muon #phi; Events", 32,-3.2,3.2);
      m_hSvc.create1D(Form("initialMuonPtPassMuTrig_%s",channel[ch].c_str()), ";p_{T,#mu}; Events", 50,0,500);
      m_hSvc.create1D(Form("initialMuonEtaPassMuTrig_%s",channel[ch].c_str()), ";muon #eta; Events", 25,-2.5,2.5);
      m_hSvc.create1D(Form("initialMuonPhiPassMuTrig_%s",channel[ch].c_str()), ";muon #phi; Events", 32,-3.2,3.2);
      m_hSvc.create1D(Form("met_%s",channel[ch].c_str()), "; missing E_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("met_plus_mwt_%s",channel[ch].c_str()), "; missing E_{T} + m^{W}_{T} [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mwt_%s",channel[ch].c_str()), "; W transverse mass [GeV]; Events", 20, 0, 200);
      m_hSvc.create1D(Form("wpt_%s",channel[ch].c_str()), "; W transverse p_{T} [GeV]; Events", 30, 0, 300);
      m_hSvc.create1D(Form("metPassMetTrig_%s",channel[ch].c_str()), "; missing E_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("met_plus_mwtPassMetTrig_%s",channel[ch].c_str()), "; missing E_{T} + m^{W}_{T} [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mwtPassMetTrig_%s",channel[ch].c_str()), "; W transverse mass [GeV]; Events", 20, 0, 200);
      m_hSvc.create1D(Form("wptPassMetTrig_%s",channel[ch].c_str()), "; W transverse p_{T} [GeV]; Events", 30, 0, 300);
      m_hSvc.create1D(Form("wptPassMuonTrig_%s",channel[ch].c_str()), "; W transverse p_{T} [GeV]; Events", 30, 0, 300);
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

AnaTtresTruth::~AnaTtresTruth() {
}

void AnaTtresTruth::run(const Event &evt, double weight, const std::string &s, int is2016run) {
  //std::cout << "AnaTtresTruth::run" << std::endl;

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

  if(!(evt.passes("no_cut"))) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(0);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(0);

  //track jet before object selection
  int nInitialTrackJet=0,nInitialTrackBtagJets=0;
  std::cout << "# of initial track jet is " << evt.nocalibtjet().size() << std::endl;
  std::cout << "# of initial muon is " << evt.nocalibmu().size() << std::endl;
  for(int itj=0;itj<evt.nocalibtjet().size();itj++){
    TLorentzVector nocalibtjet = evt.nocalibtjet()[itj].mom();
    for(int im=0;im<evt.nocalibmu().size();im++){
      TLorentzVector nocalibmu = evt.nocalibmu()[im].mom();
      float deta = nocalibtjet.Eta() - nocalibmu.Eta();
      float dphi = acos(cos(nocalibtjet.Phi() - nocalibmu.Phi()));
      float dr = sqrt(deta*deta+dphi*dphi);
      h->h1D("drTrackJetMuon", "",suffix)->Fill(dr,weight);
    }
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
  h->h1D("nInitialTrackJet_cut0", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut0", "", suffix)->Fill(evt.tjet().size(), weight);
  h->h1D("nInitialBtaggedTrackJet", "", suffix)->Fill(nInitialTrackBtagJets, weight);

  //muon trigger efficeincy
  for(int im=0;im<evt.inimu().size();im++){
    float pt = evt.inimu()[im].pt();
    float eta = evt.inimu()[im].eta();
    float phi = evt.inimu()[im].phi();
    float ptvarcone30 = evt.inimu()[im].ptvarcone30();
    float iso = ptvarcone30 / pt;
    int quality = evt.inimu()[im].quality();
    int accept = evt.inimu()[im].accept();
    float d0sig = evt.inimu()[im].d0sig();
    float z0sintheta = evt.inimu()[im].z0sintheta();
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float deta = eta - evt.jet()[jet_idx].mom().Eta(); 
      float dphi = acos(cos(phi - evt.jet()[jet_idx].mom().Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    //if(pt*1e-3<25) continue;
    if(fabs(eta)>2.5) continue;
    if(iso>0.06) continue;
    if(quality>1) continue;
    if(!accept) continue;
    if(d0sig>3) continue;
    if(z0sintheta>0.5) continue;
    if(closejl_deltaR<0.04+10./(pt*1e-3)) continue;
    h->h1D(Form("initialMuonPt_%s",channel.c_str()), "", suffix)->Fill(pt*1e-3,weight);
    h->h1D(Form("initialMuonEta_%s",channel.c_str()), "", suffix)->Fill(eta,weight);
    h->h1D(Form("initialMuonPhi_%s",channel.c_str()), "", suffix)->Fill(phi,weight);
    if(evt.inimu()[im].trigger_match()){
      h->h1D(Form("initialMuonPtPassMuTrig_%s",channel.c_str()), "", suffix)->Fill(pt*1e-3,weight);
      h->h1D(Form("initialMuonEtaPassMuTrig_%s",channel.c_str()), "", suffix)->Fill(eta,weight);
      h->h1D(Form("initialMuonPhiPassMuTrig_%s",channel.c_str()), "", suffix)->Fill(phi,weight);
    }
  }

  //Trigger selection
  bool pass_mu26_ivarmedium = evt.trigger("HLT_mu26_ivarmedium"); 
  bool pass_mu50 = evt.trigger("HLT_mu50"); 
  bool pass_met = evt.trigger("HLT_xe120_mht");
  //if(!pass_mu26_ivarmedium && !pass_mu50) return;
  //if(!pass_met) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(1);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(1);
  //std::cout << "Pass Trigger selection" << std::endl;
  h->h1D("nInitialTrackJet_cut1", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut1", "", suffix)->Fill(evt.tjet().size(), weight);

  //Number of muon
  int nMuon = 0, nMuon25=0, muon_id=-1;
  for(int im=0;im<evt.muon().size();im++){
    //std::cout << "muon " << im << std::endl;
    TLorentzVector mu = evt.muon()[im].mom();
    bool isTight = evt.muon()[im].isTight();
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float dphi = mu.DeltaPhi(evt.jet()[jet_idx].mom());
      float dy = mu.Rapidity() - evt.jet()[jet_idx].mom().Rapidity();
      float deltaR_tmp = std::sqrt(dy*dy + dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    float iso  = evt.muon()[im].ptvarcone30() / evt.muon()[im].mom().Pt();
    h->h1D("muon_ptvarcone30_over_pt", "", suffix)->Fill(iso, weight);
    //std::cout << "pt/isTight=" << "/" << mu.Perp()*1e-3 << "/" << isTight << std::endl;
    //std::cout << "ptvarcone30/iso=" << evt.muon()[im].ptvarcone30() << "/" << iso << std::endl;
    //if( closejl_deltaR < 0.04 + 10./(mu.Perp()*1e-3) ) continue;
    h->h2D("lepPt_vs_drMuonJet", "", suffix)->Fill(mu.Perp()*1e-3,closejl_deltaR, weight);
    h->h1D("muonPt", "", suffix)->Fill(mu.Perp()*1e-3, weight);
    //std::cout << "pass OR" << std::endl;
    //if(!isTight && iso>0.06) continue;
    //if(iso>0.06) continue;
    if(mu.Perp()*1e-3 > 25) nMuon25++;
    if(mu.Perp()*1e-3 > 30) {
      if(nMuon==0) muon_id=im;
      nMuon++;
    }
  }

  //Number of muon > 0
  h->h1D("nMuonAll", "", suffix)->Fill(nMuon, weight);
  h->h1D(Form("nMuon%s",channel.c_str()), "", suffix)->Fill(nMuon, weight);
  if(nMuon==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(2);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(2);
  //std::cout << "Pass #Muon >= 1" << std::endl;
  h->h1D("nInitialTrackJet_cut2", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut2", "", suffix)->Fill(evt.tjet().size(), weight);

  //Number of muon == 1
  if(nMuon25>1) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(3);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(3);
  //std::cout << "Pass #Muon == 1" << std::endl;
  h->h1D("nInitialTrackJet_cut3", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut3", "", suffix)->Fill(evt.tjet().size(), weight);

  //Initial muon
  int nIniMuon=0;
  float maxpt=0;
  for(int imu=0;imu<evt.inimu().size();imu++){
    float pt = evt.inimu()[imu].pt()*1e-3;
    float eta = evt.inimu()[imu].eta();
    float phi = evt.inimu()[imu].phi();
    std::cout << "initial muon pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl;
    h->h1D("initialMuonPt", "", suffix)->Fill(pt, weight);
    if(pt<25) continue;
    if(pt>maxpt) maxpt=pt;
    nIniMuon++;
  }
  h->h1D(Form("nMuonNoIso%s",channel.c_str()), "", suffix)->Fill(nIniMuon, weight);
  if(channel=="SemiLepMu"){
    if(nIniMuon!=0) h->h1D("leadingInitialMuonPt", "", suffix)->Fill(maxpt, weight);
    h->h1D("truthMuonPt", "", suffix)->Fill(semilepmu.Perp()*1e-3, weight);
  }

  //Number of electron == 0
  int nElectron = 0;
  for(int ie=0;ie<evt.electron().size();ie++){
    //std::cout << "electron " << ie << std::endl;
    TLorentzVector el = evt.electron()[ie].mom();
    if(ie==0) h->h1D("electronPt", "", suffix)->Fill(el.Perp()*1e-3, weight);
    if(el.Perp()*1e-3 > 25) nElectron++;
    bool isTight = evt.electron()[ie].isTightPP();
    //std::cout << "pt/isTight=" << "/" << el.Perp()*1e-3 << "/" << isTight << std::endl;
  }
  h->h1D("nElectron", "", suffix)->Fill(nElectron, weight);
  if(nElectron!=0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(4);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(4);
  //std::cout << "Pass #Electron == 0" << std::endl;
  h->h1D("nInitialTrackJet_cut4", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut4", "", suffix)->Fill(evt.tjet().size(), weight);

  //Trigger matching
  trig1 = evt.muon()[muon_id].HLT_mu26_ivarmedium();
  trig2 = evt.muon()[muon_id].HLT_mu50();
  std::cout << "trig match mu26/mu50=" << trig1 << "/" << trig2 << std::endl;
  if(!trig1 && !trig2) return;
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
  TLorentzVector l = evt.muon()[muon_id].mom();
  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3;
  float MET = evt.met().Perp()*1e-3;
  float MET_plus_mWt = mWt + MET;
  float wpt = (evt.met()+l).Perp()*1e-3;
  //std::cout << "MET/MWT=" << MET << "/" << mWt << std::endl;
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

  //if(wpt<220) return;

  if(wpt<210) {
    //if(!pass_mu26_ivarmedium && !pass_mu50) return;
  }
  else {
    //if(!pass_met) return;
  }

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
  std::cout << "Pass #Small Jet" << std::endl;
  std::cout << "trig match v2 mu26/mu50=" << trig1 << "/" << trig2 << std::endl;
  h->h1D("nInitialTrackJet_cut9", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut9", "", suffix)->Fill(evt.tjet().size(), weight);

  //Jet close to lepton
  size_t close_idx = 0;
  for (; close_idx < evt.jet().size(); ++close_idx){    
    float dphi = l.DeltaPhi(evt.jet()[close_idx].mom());
    float dy = l.Rapidity() - evt.jet()[close_idx].mom().Rapidity();
    float deltaR_tmp = std::sqrt(dy*dy + dphi*dphi);
    if(evt.jet()[close_idx].mom().Pt()*1e-3 > 25 && deltaR_tmp < 1.5) break;
  }//for     
  if(close_idx==evt.jet().size()) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(10);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(10);
  std::cout << "Pass Jet Close to Lepton" << std::endl;
  std::cout << "trig match v3 mu26/mu50=" << trig1 << "/" << trig2 << std::endl;
  h->h1D("nInitialTrackJet_cut10", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut10", "", suffix)->Fill(evt.tjet().size(), weight);

  //Number of large jet
  int nLargeJet=0,nLargeJetPt=0;
  std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
  for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
    const TLorentzVector &lj = evt.largeJet()[ljid].mom();
    bool good = evt.largeJet()[ljid].good();//default
    std::cout << "large jet pt/eta/good=" << lj.Perp()*1e-3 << "/" << lj.Eta() << "/" << good << std::endl;
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0) nLargeJetPt++;
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0 && good) nLargeJet++;
    h->h1D("largeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    h->h1D("largeJetEta", "", suffix)->Fill(lj.Eta(), weight);
    if(ljid==0) {
      h->h1D("leadingLargeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      h->h1D("leadingLargeJetEta", "", suffix)->Fill(lj.Eta(), weight);
    }
  }
  h->h1D("nLargeJets", "", suffix)->Fill(nLargeJet, weight);
  if(nLargeJet==0) return;
  //if(nLargeJetPt==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(11);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(11);
  std::cout << "Pass #Large Jet" << std::endl;
  std::cout << "trig match v4 mu26/mu50=" << trig1 << "/" << trig2 << std::endl;
  h->h1D("nInitialTrackJet_cut11", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut11", "", suffix)->Fill(evt.tjet().size(), weight);

  //Angular cut
  int nGoodJets=0;
  int ljetid=-1;
  TLorentzVector sj = evt.jet()[close_idx].mom();
  for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
    const TLorentzVector &lj = evt.largeJet()[ljid].mom();
    bool good = evt.largeJet()[ljid].good();//default
    if(lj.Perp()*1e-3 < 300 || fabs(lj.Eta())>2.0 || !good) continue;
    float delta_phi_large_jet_lepton = l.DeltaPhi(lj);
    float delta_r_large_jet_small_jet = lj.DeltaR(sj);
    //std::cout << "large jet pt/eta/good=" << lj.Perp()*1e-3 << "/" << lj.Eta() << "/" << good << std::endl;
    //std::cout << "dphi_lj_lep/dr_lj_sj=" << delta_phi_large_jet_lepton << "/" << delta_r_large_jet_small_jet << std::endl;
    if(fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) {
      if(nGoodJets==0) ljetid=ljid;
      nGoodJets++;
    }
  }
  //std::cout << "nGoodJets=" << nGoodJets << std::endl;
  if(nGoodJets==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(12);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(12);
  std::cout << "Pass Angular Cut" << std::endl;
  std::cout << "trig match v5 mu26/mu50=" << trig1 << "/" << trig2 << std::endl;
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
 // std::cout << "Number of track b-tagged jets is " << nTrackBtagJets << std::endl;
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
  std::cout << "Pass B-tag" << std::endl;
  std::cout << "trig match v6 mu26/mu50=" << trig1 << "/" << trig2 << std::endl;
  h->h1D("nInitialTrackJet_cut13", "", suffix)->Fill(nInitialTrackJet, weight);
  h->h1D("nTrackJet_cut13", "", suffix)->Fill(evt.tjet().size(), weight);

  h->h1D("nMuonFinalAll", "", suffix)->Fill(nMuon, weight);
  h->h1D(Form("nMuonFinal%s",channel.c_str()), "", suffix)->Fill(nMuon, weight);

  std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), true);   
  TLorentzVector nu(0,0,0,0);
  std::cout << "vec_nu.size=" << vec_nu.size() << std::endl;
  if (vec_nu.size() > 0) {
    nu = *(vec_nu[0]);
    for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
    vec_nu.clear();
  }
  const TLorentzVector &ljt = evt.largeJet()[ljetid].mom();
  float mtt = (ljt+sj+nu+l).M();
  float mlt = (sj+nu+l).M();
  float mht = ljt.M();
  float ljet_pt = ljt.Perp();
  float ljet_eta = ljt.Eta();
  float ljet_phi = ljt.Phi();
  //if(mlt*1e-3>300) return;
  h->h1D("mttAll", "", suffix)->Fill(mtt*1e-3, weight);
  h->h1D(Form("mtt%s",channel.c_str()), "", suffix)->Fill(mtt*1e-3, weight);
  h->h1D(Form("mlt%s",channel.c_str()), "", suffix)->Fill(mlt*1e-3, weight);
  h->h1D(Form("mht%s",channel.c_str()), "", suffix)->Fill(mht*1e-3, weight);
  h->h1D("cutFlowAll", "", suffix)->Fill(14);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(14);

  float mtt_truth = (t+tbar).M();
  float mlt_truth = (isLeptonicT)? t.M() : tbar.M();
  float mht_truth = (isLeptonicT)? tbar.M() : t.M();
  float ht_pt_truth = (isLeptonicT)? tbar.Perp() : t.Perp();
  float ht_eta_truth = (isLeptonicT)? tbar.Eta() : t.Eta();
  float ht_phi_truth = (isLeptonicT)? tbar.Phi() : t.Phi();
  float nupz_truth = (isLeptonicT)? evt.MC_Wdecay1_from_t().Pz() : evt.MC_Wdecay2_from_tbar().Pz();  
  float nueta_truth = (isLeptonicT)? evt.MC_Wdecay1_from_t().Eta() : evt.MC_Wdecay2_from_tbar().Eta();  
  float nuphi_truth = (isLeptonicT)? evt.MC_Wdecay1_from_t().Phi() : evt.MC_Wdecay2_from_tbar().Phi();  
  float dr_Wb_had_truth = (isLeptonicT)? W_from_tbar.DeltaR(b_from_tbar) : W_from_t.DeltaR(b_from_t);
  float dr_ljet_trueb_had = (isLeptonicT)? ljt.DeltaR(b_from_tbar) : ljt.DeltaR(b_from_t);
  float dr_ljet_true_htop = (isLeptonicT)? ljt.DeltaR(tbar) : ljt.DeltaR(t);
  float dr_ljet_true_b_from_ltop = (isLeptonicT)? ljt.DeltaR(b_from_t) : ljt.DeltaR(b_from_tbar);
  float dr_sjet_true_b_from_ltop = (isLeptonicT)? sj.DeltaR(b_from_t) : sj.DeltaR(b_from_tbar);
  float dr_muon_true_muon = (isLeptonicT)? l.DeltaR(evt.MC_Wdecay2_from_t()) : l.DeltaR(evt.MC_Wdecay1_from_tbar());
  float res_mlt = (mlt_truth - mlt)/mlt_truth;
  float res_mht = (mht_truth - mht)/mht_truth;
  float res_mtt = (mtt_truth - mtt)/mtt_truth;
  float res_ht_pt = (ht_pt_truth - ljet_pt)/ht_pt_truth;
  float res_ht_eta = (ht_eta_truth - ljet_eta)/ht_eta_truth;
  float res_ht_phi = (ht_phi_truth - ljet_phi)/ht_phi_truth;
  h->h1D(Form("drTrueWTrueBHad%s",channel.c_str()), "", suffix)->Fill(dr_Wb_had_truth, weight);
  h->h1D(Form("drLargeJetTrueBHad%s",channel.c_str()), "", suffix)->Fill(dr_ljet_trueb_had, weight);
  h->h1D(Form("drLargeJetTrueHtop%s",channel.c_str()), "", suffix)->Fill(dr_ljet_true_htop, weight);
  h->h1D(Form("drLargeJetTrueBfromLtop%s",channel.c_str()), "", suffix)->Fill(dr_ljet_true_b_from_ltop, weight);
  h->h1D(Form("drSmallJetTrueBfromLtop%s",channel.c_str()), "", suffix)->Fill(dr_sjet_true_b_from_ltop, weight);
  h->h1D(Form("drMuonTrueMuon%s",channel.c_str()), "", suffix)->Fill(dr_muon_true_muon, weight);
  h->h2D(Form("mhtVsDrTrueWTrueBHad%s",channel.c_str()), "", suffix)->Fill(mht*1e-3,dr_Wb_had_truth);
  h->h2D(Form("mhtVsDrLargeJetTrueBHad%s",channel.c_str()), "", suffix)->Fill(mht*1e-3,dr_ljet_trueb_had);
  h->h1D(Form("resMtt%s",channel.c_str()), "", suffix)->Fill(res_mtt, weight);
  h->h1D(Form("resMlt%s",channel.c_str()), "", suffix)->Fill(res_mlt, weight);
  h->h1D(Form("resMht%s",channel.c_str()), "", suffix)->Fill(res_mht, weight);
  h->h1D(Form("resLargeJetPt%s",channel.c_str()), "", suffix)->Fill(res_ht_pt, weight);
  h->h1D(Form("resLargeJetEta%s",channel.c_str()), "", suffix)->Fill(res_ht_eta, weight);
  h->h1D(Form("resLargeJetPhi%s",channel.c_str()), "", suffix)->Fill(res_ht_phi, weight);
  h->h2D(Form("mhtVsTrueMht%s",channel.c_str()), "", suffix)->Fill(mht*1e-3,mht_truth*1e-3);
  h->h2D(Form("resMttVsTrueMtt%s",channel.c_str()), "", suffix)->Fill(res_mtt,mtt_truth*1e-3);
  h->h2D(Form("resMttVsMtt%s",channel.c_str()), "", suffix)->Fill(res_mtt,mtt*1e-3);
  h->h2D(Form("resMttVsMlt%s",channel.c_str()), "", suffix)->Fill(res_mtt,mlt*1e-3);
  h->h2D(Form("resMttVsMht%s",channel.c_str()), "", suffix)->Fill(res_mtt,mht*1e-3);
  h->h2D(Form("resMttVsResMht%s",channel.c_str()), "", suffix)->Fill(res_mtt,res_mht);
  if(res_mtt>0.6) {
    h->h1D(Form("mhtBadMtt%s",channel.c_str()), "", suffix)->Fill(mht*1e-3, weight);
    h->h1D(Form("resLargeJetPtBadMtt%s",channel.c_str()), "", suffix)->Fill(res_ht_pt, weight);
    h->h1D(Form("resLargeJetEtaBadMtt%s",channel.c_str()), "", suffix)->Fill(res_ht_eta, weight);
    h->h1D(Form("resLargeJetPhiBadMtt%s",channel.c_str()), "", suffix)->Fill(res_ht_phi, weight);
    h->h1D(Form("resMhtBadMtt%s",channel.c_str()), "", suffix)->Fill(res_mht, weight);
    h->h1D(Form("drTrueWTrueBHadBadMtt%s",channel.c_str()), "", suffix)->Fill(dr_Wb_had_truth, weight);
    h->h1D(Form("drLargeJetTrueBHadBadMtt%s",channel.c_str()), "", suffix)->Fill(dr_ljet_trueb_had, weight);
    h->h1D(Form("drLargeJetTrueHtopBadMtt%s",channel.c_str()), "", suffix)->Fill(dr_ljet_true_htop, weight);
    h->h1D(Form("drSmallJetTrueBfromLtopBadMtt%s",channel.c_str()), "", suffix)->Fill(dr_sjet_true_b_from_ltop, weight);
    h->h1D(Form("drMuonTrueMuonBadMtt%s",channel.c_str()), "", suffix)->Fill(dr_muon_true_muon, weight);
    h->h1D(Form("drLargeJetTrueBfromLtopBadMtt%s",channel.c_str()), "", suffix)->Fill(dr_ljet_true_b_from_ltop, weight);
  }
  else {
    h->h1D(Form("mhtGoodMtt%s",channel.c_str()), "", suffix)->Fill(mht*1e-3, weight);
    h->h1D(Form("resLargeJetPtGoodMtt%s",channel.c_str()), "", suffix)->Fill(res_ht_pt, weight);
    h->h1D(Form("resLargeJetEtaGoodMtt%s",channel.c_str()), "", suffix)->Fill(res_ht_eta, weight);
    h->h1D(Form("resLargeJetPhiGoodMtt%s",channel.c_str()), "", suffix)->Fill(res_ht_phi, weight);
    h->h1D(Form("resMhtGoodMtt%s",channel.c_str()), "", suffix)->Fill(res_mht, weight);
    h->h1D(Form("drTrueWTrueBHadGoodMtt%s",channel.c_str()), "", suffix)->Fill(dr_Wb_had_truth, weight);
    h->h1D(Form("drLargeJetTrueBHadGoodMtt%s",channel.c_str()), "", suffix)->Fill(dr_ljet_trueb_had, weight);
    h->h1D(Form("drLargeJetTrueHtopGoodMtt%s",channel.c_str()), "", suffix)->Fill(dr_ljet_true_htop, weight);
    h->h1D(Form("drSmallJetTrueBfromLtopGoodMtt%s",channel.c_str()), "", suffix)->Fill(dr_sjet_true_b_from_ltop, weight);
    h->h1D(Form("drMuonTrueMuonGoodMtt%s",channel.c_str()), "", suffix)->Fill(dr_muon_true_muon, weight);
    h->h1D(Form("drLargeJetTrueBfromLtopGoodMtt%s",channel.c_str()), "", suffix)->Fill(dr_ljet_true_b_from_ltop, weight);
  }

  float met_alpha = m_neutrinoBuilder.fitAlpha(&l, evt.met().Perp(), evt.met().Phi());   
  bool hasRealSol = m_neutrinoBuilder.hasRealSolution(&l, evt.met().Perp(), evt.met().Phi());   
  std::cout << "hasRealSol=" << hasRealSol << std::endl;
  std::vector<TLorentzVector*> vec_nu_rotate_small = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), true);   
  std::vector<TLorentzVector*> vec_nu_rotate_big = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), false);   
  std::vector<TLorentzVector*> vec_nu_real_small = m_neutrinoBuilder.candidatesFromWMass_RealPart(&l, evt.met().Perp(), evt.met().Phi(), true);   
  std::vector<TLorentzVector*> vec_nu_real_big = m_neutrinoBuilder.candidatesFromWMass_RealPart(&l, evt.met().Perp(), evt.met().Phi(), false);   
  std::cout << "neutrion pz rotate_small=" << vec_nu_rotate_small[0]->Pz() << std::endl;
  std::cout << "neutrion pz rotate_big=" << vec_nu_rotate_big[0]->Pz() << std::endl;
  std::cout << "neutrion pz real_small=" << vec_nu_real_small[0]->Pz() << std::endl;
  std::cout << "neutrion pz real_big=" << vec_nu_real_big[0]->Pz() << std::endl;
  float mtt_rotate_small = (ljt+sj+l+*(vec_nu_rotate_small[0])).M();
  float mtt_rotate_big = (ljt+sj+l+*(vec_nu_rotate_big[0])).M();
  float mtt_realpart = (ljt+sj+l+*(vec_nu_real_small[0])).M();
  float mlt_rotate_small = (sj+l+*(vec_nu_rotate_small[0])).M();
  float mlt_rotate_big = (sj+l+*(vec_nu_rotate_big[0])).M();
  float mlt_realpart = (sj+l+*(vec_nu_real_small[0])).M();
  float res_mtt_rotate_small = (mtt_truth-mtt_rotate_small)/mtt_truth;
  float res_mtt_rotate_big = (mtt_truth-mtt_rotate_big)/mtt_truth;
  float res_mtt_realpart = (mtt_truth-mtt_realpart)/mtt_truth;
  float res_mlt_rotate_small = (mlt_truth-mlt_rotate_small)/mlt_truth;
  float res_mlt_rotate_big = (mlt_truth-mlt_rotate_big)/mlt_truth;
  float res_mlt_realpart = (mlt_truth-mlt_realpart)/mlt_truth;
  float res_nupz_rotate_small = (nupz_truth - vec_nu_rotate_small[0]->Pz())/fabs(nupz_truth);
  float res_nupz_rotate_big = (nupz_truth - vec_nu_rotate_big[0]->Pz())/fabs(nupz_truth);
  float res_nupz_realpart = (nupz_truth - vec_nu_real_small[0]->Pz())/fabs(nupz_truth);
  float res_nueta_rotate_small = (nueta_truth - vec_nu_rotate_small[0]->Eta())/fabs(nueta_truth);
  float res_nueta_rotate_big = (nueta_truth - vec_nu_rotate_big[0]->Eta())/fabs(nueta_truth);
  float res_nueta_realpart = (nueta_truth - vec_nu_real_small[0]->Eta())/fabs(nueta_truth);
  float res_nuphi_rotate_small = acos(cos(nuphi_truth - vec_nu_rotate_small[0]->Phi()))/fabs(nuphi_truth);
  float res_nuphi_rotate_big = acos(cos(nuphi_truth - vec_nu_rotate_big[0]->Phi()))/fabs(nuphi_truth);
  float res_nuphi_realpart = acos(cos(nuphi_truth - vec_nu_real_small[0]->Phi()))/fabs(nuphi_truth);
  if(hasRealSol) {
    h->h1D(Form("mttTruthReal%s",channel.c_str()), "", suffix)->Fill(mtt_truth*1e-3, weight);
    h->h1D(Form("mltTruthReal%s",channel.c_str()), "", suffix)->Fill(mlt_truth*1e-3, weight);
    h->h1D(Form("nupzTruthReal%s",channel.c_str()), "", suffix)->Fill(fabs(nupz_truth)*1e-3, weight);
    h->h1D(Form("mttRealSolSmall%s",channel.c_str()), "", suffix)->Fill(mtt_rotate_small*1e-3, weight);
    h->h1D(Form("mttRealSolBig%s",channel.c_str()), "", suffix)->Fill(mtt_rotate_big*1e-3, weight);
    h->h1D(Form("mltRealSolSmall%s",channel.c_str()), "", suffix)->Fill(mlt_rotate_small*1e-3, weight);
    h->h1D(Form("mltRealSolBig%s",channel.c_str()), "", suffix)->Fill(mlt_rotate_big*1e-3, weight);
    h->h1D(Form("nupzRealSolSmall%s",channel.c_str()), "", suffix)->Fill(fabs(vec_nu_rotate_small[0]->Pz())*1e-3, weight);
    h->h1D(Form("nupzRealSolBig%s",channel.c_str()), "", suffix)->Fill(fabs(vec_nu_rotate_big[0]->Pz())*1e-3, weight);
    h->h1D(Form("resMttRealSolSmall%s",channel.c_str()), "", suffix)->Fill(res_mtt_rotate_small, weight);
    h->h1D(Form("resMttRealSolBig%s",channel.c_str()), "", suffix)->Fill(res_mtt_rotate_big, weight);
    h->h1D(Form("resMltRealSolSmall%s",channel.c_str()), "", suffix)->Fill(res_mlt_rotate_small, weight);
    h->h1D(Form("resMltRealSolBig%s",channel.c_str()), "", suffix)->Fill(res_mlt_rotate_big, weight);
    h->h1D(Form("resNupzRealSolSmall%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_small, weight);
    h->h1D(Form("resNupzRealSolBig%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_big, weight);
    if(fabs(nupz_truth)*1e-3<200){
      h->h1D(Form("resNupzLowRealSolSmall%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_small, weight);
      h->h1D(Form("resNupzLowRealSolBig%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_big, weight);
    } else{
      h->h1D(Form("resNupzHighRealSolSmall%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_small, weight);
      h->h1D(Form("resNupzHighRealSolBig%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_big, weight);
    }
    h->h1D(Form("resNuetaRealSolSmall%s",channel.c_str()), "", suffix)->Fill(res_nueta_rotate_small, weight);
    h->h1D(Form("resNuetaRealSolBig%s",channel.c_str()), "", suffix)->Fill(res_nueta_rotate_big, weight);
    h->h1D(Form("resNuphiRealSolSmall%s",channel.c_str()), "", suffix)->Fill(res_nuphi_rotate_small, weight);
    h->h1D(Form("resNuphiRealSolBig%s",channel.c_str()), "", suffix)->Fill(res_nuphi_rotate_big, weight);
    h->h1D(Form("mttNew%s",channel.c_str()), "", suffix)->Fill(mtt_rotate_big*1e-3, weight);
    h->h1D(Form("mltNew%s",channel.c_str()), "", suffix)->Fill(mlt_rotate_big*1e-3, weight);
  }
  else{
    if(fabs(res_nupz_rotate_small-1)<0.01){
      std::cout << "met alpha=" << met_alpha << std::endl;
      std::cout << "truth nu pz/eta/phi=" << nupz_truth << "/" << nueta_truth << "/" << nuphi_truth << std::endl;
      std::cout << "nupz rotate/realpart=" << vec_nu_rotate_small[0]->Pz() << "/" << vec_nu_real_small[0]->Pz() << std::endl; 
      std::cout << "nueta rotate/realpart=" << vec_nu_rotate_small[0]->Eta() << "/" << vec_nu_real_small[0]->Eta() << std::endl; 
      std::cout << "nuphi rotate/realpart=" << vec_nu_rotate_small[0]->Phi() << "/" << vec_nu_real_small[0]->Phi() << std::endl; 
    }
    h->h1D(Form("mttTruthImag%s",channel.c_str()), "", suffix)->Fill(mtt_truth*1e-3, weight);
    h->h1D(Form("mltTruthImag%s",channel.c_str()), "", suffix)->Fill(mlt_truth*1e-3, weight);
    h->h1D(Form("nupzTruthImag%s",channel.c_str()), "", suffix)->Fill(fabs(nupz_truth)*1e-3, weight);
    h->h1D(Form("mttImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(mtt_realpart*1e-3, weight);
    h->h1D(Form("mttImagSolRotate%s",channel.c_str()), "", suffix)->Fill(mtt_rotate_small*1e-3, weight);
    h->h1D(Form("mltImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(mlt_realpart*1e-3, weight);
    h->h1D(Form("mltImagSolRotate%s",channel.c_str()), "", suffix)->Fill(mlt_rotate_small*1e-3, weight);
    h->h1D(Form("nupzImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(fabs(vec_nu_real_small[0]->Pz())*1e-3, weight);
    h->h1D(Form("nupzImagSolRotate%s",channel.c_str()), "", suffix)->Fill(fabs(vec_nu_rotate_small[0]->Pz())*1e-3, weight);
    h->h1D(Form("resMttImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(res_mtt_realpart, weight);
    h->h1D(Form("resMttImagSolRotate%s",channel.c_str()), "", suffix)->Fill(res_mtt_rotate_small, weight);
    h->h1D(Form("resMltImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(res_mlt_realpart, weight);
    h->h1D(Form("resMltImagSolRotate%s",channel.c_str()), "", suffix)->Fill(res_mlt_rotate_small, weight);
    h->h1D(Form("resNupzImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(res_nupz_realpart, weight);
    h->h1D(Form("resNupzImagSolRotate%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_small, weight);
    h->h2D(Form("resNupzImagSolRotateVsMetAlpha%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_small, met_alpha,weight);
    if(fabs(nupz_truth)*1e-3<200){
      h->h1D(Form("resNupzLowImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(res_nupz_realpart, weight);
      h->h1D(Form("resNupzLowImagSolRotate%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_small, weight);
    } else{
      h->h1D(Form("resNupzHighImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(res_nupz_realpart, weight);
      h->h1D(Form("resNupzHighImagSolRotate%s",channel.c_str()), "", suffix)->Fill(res_nupz_rotate_small, weight);
    }
    h->h1D(Form("resNuetaImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(res_nueta_realpart, weight);
    h->h1D(Form("resNuetaImagSolRotate%s",channel.c_str()), "", suffix)->Fill(res_nueta_rotate_small, weight);
    h->h1D(Form("resNuphiImagSolRealPart%s",channel.c_str()), "", suffix)->Fill(res_nuphi_realpart, weight);
    h->h1D(Form("resNuphiImagSolRotate%s",channel.c_str()), "", suffix)->Fill(res_nuphi_rotate_small, weight);
    h->h1D(Form("mttNew%s",channel.c_str()), "", suffix)->Fill(mtt_rotate_small*1e-3, weight);
    h->h1D(Form("mltNew%s",channel.c_str()), "", suffix)->Fill(mlt_rotate_small*1e-3, weight);
  }

  //truth jet matching
  int bt_sj_idx=-1, btbar_sj_idx=-1;
  float min_dr_sj_bt = 999., min_dr_sj_btbar = 999.;
  std::cout << "Number of small jet is " << evt.jet().size() << std::endl;
  for (size_t bidx = 0; bidx < evt.jet().size(); ++bidx){
    float sj_eta = evt.jet()[bidx].mom().Eta();
    float sj_phi = evt.jet()[bidx].mom().Phi();
    float drt = TuDoAtlas::calc_delta_r(sj_eta,b_from_t_eta,sj_phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(sj_eta,b_from_tbar_eta,sj_phi,b_from_tbar_phi);
    float pt = evt.jet()[bidx].mom().Pt();
    float eta = evt.jet()[bidx].mom().Eta();
    float phi = evt.jet()[bidx].mom().Phi();
    int isBtag = evt.jet()[bidx].btag_mv2c20_70();
    std::cout << "small jet pt/eta/phi/btag=" << pt << "/" << eta << "/" << phi << "/" << isBtag << std::endl;
    //std::cout << "btag value=" << evt.jet()[bidx].mv2c20() << std::endl;
    if(!isLeptonicT) h->h1D("drJetHB", "", suffix)->Fill(drt, weight);  
    if(!isLeptonicTbar) h->h1D("drJetHB", "", suffix)->Fill(drtbar, weight);  
    if(drt<0.2 && drt<min_dr_sj_bt) {
      min_dr_sj_bt = drt;
      bt_sj_idx=bidx;
    }
    if(drtbar<0.2 && drtbar<min_dr_sj_btbar) {
      min_dr_sj_btbar = drtbar;
      btbar_sj_idx=bidx;
    }
  }

  //truth tjet matching
  int bt_tj_idx=-1, btbar_tj_idx=-1;
  float min_dr_tj_bt = 999., min_dr_tj_btbar = 999.;
  //std::cout << "Number of small tjet is " << evt.tjet().size() << std::endl;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
    float pt = evt.tjet()[bidx].mom().Pt();
    float eta = evt.tjet()[bidx].mom().Eta();
    float phi = evt.tjet()[bidx].mom().Phi();
    float drt = TuDoAtlas::calc_delta_r(eta,b_from_t_eta,phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(eta,b_from_tbar_eta,phi,b_from_tbar_phi);
    int isBtag = evt.tjet()[bidx].btag_mv2c10_70_trk();
    //std::cout << "track jet pt/eta/phi/btag=" << pt << "/" << eta << "/" << phi << "/" << isBtag << std::endl;
    //std::cout << "btag value=" << evt.tjet()[bidx].mv2c10() << std::endl;
    if(!isLeptonicT) h->h1D("drTjetHB", "", suffix)->Fill(drt, weight);  
    if(!isLeptonicTbar) h->h1D("drTjetHB", "", suffix)->Fill(drtbar, weight);  
    if(drt<0.2 && drt<min_dr_tj_bt) {
      min_dr_tj_bt = drt;
      bt_tj_idx=bidx;
    }
    if(drtbar<0.2 && drtbar<min_dr_tj_btbar) {
      min_dr_tj_btbar = drtbar;
      btbar_tj_idx=bidx;
    }
  }

  //calo jet b-tag efficiency
  std::cout << "bt matched small jet id=" << bt_sj_idx << std::endl;
  //std::cout << "btbar matched small jet id=" << btbar_sj_idx << std::endl;
  if(bt_sj_idx!=-1){
    h->h1D("caloBtagBDT", "", suffix)->Fill(evt.jet()[bt_sj_idx].mv2c20());
    h->h1D("caloBtagBDT", "", suffix)->Fill(evt.jet()[bt_sj_idx].mv2c20());
    //std::cout << "calo btag value bt/btbar=" << evt.jet()[bt_sj_idx].mv2c20() << "/" << evt.jet()[btbar_sj_idx].mv2c20() << std::endl;
    //std::cout << "track btag value bt/btbar=" << evt.tjet()[bt_tj_idx].mv2c10() << "/" << evt.tjet()[btbar_tj_idx].mv2c10() << std::endl;
    if(isLeptonicT){//leptonic top
      std::cout << "jet is_flatbtag=" << evt.jet()[bt_sj_idx].is_flatbtag() << std::endl;
      h->h1D("ltop_pt_jet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_ltop_pt_jet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_ltop_matched_jet_pt", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[bt_sj_idx].btag_mv2c20_70()){
        h->h1D("ltop_pt_jet_passCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
      if(evt.jet()[bt_sj_idx].is_flatbtag()){
        h->h1D("ltop_pt_jet_passFlatCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_jet_passFlatCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_jet_pt_passFlatCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
    }
    else{//hadronic top
      std::cout << "jet is_flatbtag=" << evt.jet()[bt_sj_idx].is_flatbtag() << std::endl;
      h->h1D("htop_pt_jet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_jet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_htop_matched_jet_pt", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[bt_sj_idx].btag_mv2c20_70()){
        h->h1D("htop_pt_jet_passCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
      if(evt.jet()[bt_sj_idx].is_flatbtag()){
        h->h1D("htop_pt_jet_passFlatCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_jet_passFlatCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_jet_pt_passFlatCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
    }
  }
  if(btbar_sj_idx!=-1){
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

  float flat_btag_thr=-1;
  //track tjet b-tag efficiency
  //std::cout << "bt matched small tjet id=" << bt_tj_idx << std::endl;
  //std::cout << "btbar matched small tjet id=" << btbar_tj_idx << std::endl;
  if(bt_tj_idx!=-1 && bt_sj_idx!=-1){ 
    float tjet_pt = evt.tjet()[bt_tj_idx].mom().Perp()*1e-3;
    if(tjet_pt<300) flat_btag_thr = 0.6455;
    else flat_btag_thr = -0.3;
    if(evt.tjet()[bt_tj_idx].btag_mv2c10_70_trk() && evt.tjet()[bt_tj_idx].mv2c10() <0.6455) std::cout << "strange!!!!" << std::endl;
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_tjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_ltop_pt_tjet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_ltop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      h->h2D("bMatchedTjetPt_vs_btagScore", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3,evt.tjet()[bt_tj_idx].mv2c10());
      if(tjet_pt<300) h->h1D("Mv2c10_trackb_lowpt", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mv2c10());
      else h->h1D("Mv2c10_trackb_highpt", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mv2c10());
      std::cout << "mv2c10/thr=" << evt.tjet()[bt_tj_idx].mv2c10() << "/" << flat_btag_thr << std::endl;
      if(evt.tjet()[bt_tj_idx].btag_mv2c10_70_trk()){
        //std::cout << "ltop tjet btagged!!" << std::endl;
        h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
      if(evt.tjet()[bt_tj_idx].mv2c10()>flat_btag_thr){
        h->h1D("ltop_pt_tjet_passFlatTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_tjet_passFlatTrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_passFlatTrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
      if(evt.tjet()[bt_tj_idx].mv2c10()>-0.1416){//80%
        h->h1D("ltop_pt_tjet_pass85TrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_tjet_pass85TrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_pass85TrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }

      //else std::cout << "ltop tjet not btagged!!" << std::endl;
    }
    else{//hadronic top
      h->h1D("htop_pt_tjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      h->h2D("bMatchedTjetPt_vs_btagScore", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3,evt.tjet()[bt_tj_idx].mv2c10());
      if(evt.tjet()[bt_tj_idx].btag_mv2c10_70_trk()){
        //std::cout << "htop tjet btagged!!" << std::endl;
        h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
      if(evt.tjet()[bt_tj_idx].mv2c10()>flat_btag_thr){
        h->h1D("htop_pt_tjet_passFlatTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_tjet_passFlatTrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_tjet_pt_passFlatTrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
      if(evt.tjet()[bt_tj_idx].mv2c10()>-0.1416){
        h->h1D("htop_pt_tjet_pass85TrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_tjet_pass85TrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_tjet_pt_pass85TrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
      //else std::cout << "htop tjet not btagged!!" << std::endl;
    }
  }
  if(btbar_tj_idx!=-1 && btbar_sj_idx!=-1){
    float tjet_pt = evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3;
    if(tjet_pt<300) flat_btag_thr = 0.6455;
    else flat_btag_thr = -0.3;
    if(isLeptonicTbar){//leptonic top
      h->h1D("ltop_pt_tjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_ltop_pt_tjet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_ltop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      if(tjet_pt<300) h->h1D("Mv2c10_trackb_lowpt", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mv2c10());
      else h->h1D("Mv2c10_trackb_highpt", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mv2c10());
      if(evt.tjet()[btbar_tj_idx].btag_mv2c10_70_trk()){
        //std::cout << "ltop tjet btagged!!" << std::endl;
        h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
      if(evt.tjet()[btbar_tj_idx].mv2c10()>flat_btag_thr){
        //std::cout << "ltop tjet btagged!!" << std::endl;
        h->h1D("ltop_pt_tjet_passFlatTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_tjet_passFlatTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_passFlatTrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
      if(evt.tjet()[btbar_tj_idx].mv2c10()>-0.1416){//85%
        h->h1D("ltop_pt_tjet_pass85TrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_tjet_pass85TrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_pass85TrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
      //else std::cout << "ltop tjet not btagged!!" << std::endl;
    }
    else{//hadronic top
      h->h1D("htop_pt_tjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      if(evt.tjet()[btbar_tj_idx].btag_mv2c10_70_trk()){
        //std::cout << "htop tjet btagged!!" << std::endl;
        h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
      if(evt.tjet()[btbar_tj_idx].mv2c10()>flat_btag_thr){
        //std::cout << "ltop tjet btagged!!" << std::endl;
        h->h1D("htop_pt_tjet_passFlatTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_tjet_passFlatTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_htop_matched_tjet_pt_passFlatTrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
      if(evt.tjet()[btbar_tj_idx].mv2c10()>-0.1416){//85%
        h->h1D("htop_pt_tjet_pass85TrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_tjet_pass85TrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_htop_matched_tjet_pt_pass85TrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
      //else std::cout << "htop tjet not btagged!!" << std::endl;
    }
  }

  //top tagging
  //std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
  int t_matched_idx=-1;
  float min_dr_t=999.;
  for (int ij=0; ij < evt.largeJet().size(); ++ij) {
    const TLorentzVector &lj = evt.largeJet()[ij].mom();
    float ljet_eta = lj.Eta();
    float ljet_phi = lj.Phi();
    float deta = t_eta-ljet_eta;
    float dphi = acos(cos(t_phi-ljet_phi));
    float dr = deta*deta+dphi*dphi; 
    h->h1D("drLargeJetTruthTop", "", suffix)->Fill(dr, weight);
    bool good_def = evt.largeJet()[ij].good();//default
    //bool good_bdt = (evt.largeJet()[ij].BDT_TOPtag()>bdtthr)? true : false;
    float delta_phi_large_jet_lepton = l.DeltaPhi(lj);
    //float delta_r_large_jet_small_jet = lj.DeltaR(sj);
    //bool pass_ljet = (lj.Perp()*1e-3>300 && fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) ?true:false;
    h->h1D("allLargeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    for(int lpt=0;lpt<4;lpt++){
      if(lj.Perp()*1e-3>500*lpt && lj.Perp()*1e-3<500+500*lpt){
        h->h1D(Form("BDT_TOPtag_allljet_pt%d_%d",500*lpt,500+500*lpt),"",suffix)->Fill(evt.largeJet()[ij].BDT_TOPtag(),weight);
        //if(pass_ljet) h->h1D(Form("BDT_TOPtag_passljet_pt%d_%d",500*lpt,500+500*lpt),"",suffix)->Fill(evt.largeJet()[ij].BDT_TOPtag(),weight);
      }
    }
    //std::cout << "top tag def/new=" << good_def << "/" << evt.largeJet()[ij].BDT_TOPtag() << std::endl;
    if(good_def) h->h1D("allLargeJetPt_passSTT", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    if(evt.largeJet()[ij].BDT_TOPtag()>-0.5) h->h1D("allLargeJetPt_passBDT_loose", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    if(evt.largeJet()[ij].BDT_TOPtag()>0) h->h1D("allLargeJetPt_passBDT_middle", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    if(evt.largeJet()[ij].BDT_TOPtag()>0.5) h->h1D("allLargeJetPt_passBDT_tight", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    //std::cout << "dr=" << dr << std::endl;
    if(dr<0.2 && dr<min_dr_t){
      min_dr_t = dr;
      t_matched_idx=ij;
    }
  }
  if(t_matched_idx!=-1){
    const TLorentzVector &matched_t = evt.largeJet()[t_matched_idx].mom();
    h->h1D("TopMatchedLargeJetPt", "", suffix)->Fill(matched_t.Perp()*1e-3, weight);
    if(evt.largeJet()[t_matched_idx].good()) h->h1D("TopMatchedLargeJetPt_passSTT", "", suffix)->Fill(matched_t.Perp()*1e-3, weight);
    if(evt.largeJet()[t_matched_idx].BDT_TOPtag()>-0.5) h->h1D("TopMatchedLargeJetPt_passBDT_loose", "", suffix)->Fill(matched_t.Perp()*1e-3, weight);
    if(evt.largeJet()[t_matched_idx].BDT_TOPtag()>0) h->h1D("TopMatchedLargeJetPt_passBDT_middle", "", suffix)->Fill(matched_t.Perp()*1e-3, weight);
    if(evt.largeJet()[t_matched_idx].BDT_TOPtag()>0.5) h->h1D("TopMatchedLargeJetPt_passBDT_tight", "", suffix)->Fill(matched_t.Perp()*1e-3, weight);
  }

  //Trigger efficiency
  for(int im=0;im<evt.muon().size();im++){
    //std::cout << "muon " << im << std::endl;
    TLorentzVector mu = evt.muon()[im].mom();
    if(mu.Perp()*1e-3<30) continue;
    bool pass_trigmu = evt.muon()[im].HLT_mu26_ivarmedium() || evt.muon()[im].HLT_mu50();
    std::cout << "pass_trigmu=" << pass_trigmu << std::endl;
    h->h1D(Form("muonPt_%s",channel.c_str()), "", suffix)->Fill(mu.Perp()*1e-3, weight);
    h->h1D(Form("muonEta_%s",channel.c_str()), "", suffix)->Fill(mu.Eta(), weight);
    h->h1D(Form("muonPhi_%s",channel.c_str()), "", suffix)->Fill(mu.Phi(), weight);
    if(pass_trigmu) {
      h->h1D(Form("muonPtPassMuTrig_%s",channel.c_str()), "", suffix)->Fill(mu.Perp()*1e-3, weight);
      h->h1D(Form("muonEtaPassMuTrig_%s",channel.c_str()), "", suffix)->Fill(mu.Eta(), weight);
      h->h1D(Form("muonPhiPassMuTrig_%s",channel.c_str()), "", suffix)->Fill(mu.Phi(), weight);
    }
  }

  h->h1D(Form("met_%s",channel.c_str()), "", suffix)->Fill(MET);
  h->h1D(Form("met_plus_mwt_%s",channel.c_str()), "", suffix)->Fill(MET_plus_mWt);
  h->h1D(Form("mwt_%s",channel.c_str()), "", suffix)->Fill(mWt);
  h->h1D(Form("wpt_%s",channel.c_str()), "", suffix)->Fill(wpt);
  if(evt.trigger("HLT_xe120_mht")){
    h->h1D(Form("metPassMetTrig_%s",channel.c_str()), "", suffix)->Fill(MET);
    h->h1D(Form("met_plus_mwtPassMetTrig_%s",channel.c_str()), "", suffix)->Fill(MET_plus_mWt);
    h->h1D(Form("mwtPassMetTrig_%s",channel.c_str()), "", suffix)->Fill(mWt);
    h->h1D(Form("wptPassMetTrig_%s",channel.c_str()), "", suffix)->Fill(wpt);
  }
  if((pass_mu26_ivarmedium || pass_mu50) && (trig1 || trig2)){
    h->h1D(Form("wptPassMuonTrig_%s",channel.c_str()), "", suffix)->Fill(wpt);
  }
}
