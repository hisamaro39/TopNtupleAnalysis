/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
//#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnalysisMVATopTag.h"
#include "TopNtupleAnalysis/AnaTtresMVATopTag.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"


AnaTtresMVATopTag::AnaTtresMVATopTag(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
  : AnalysisMVATopTag(filename, systList), m_electron(electron), m_boosted(boosted),
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

    double varBin8[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500,3000};
    int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;
    double varBin9[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000};
    int varBinN9 = sizeof(varBin9)/sizeof(double) - 1;

    m_hSvc.create1D("yields", "; One ; Events", 1, 0.5, 1.5);

    //m_hSvc.create1DVar("lepPt", "; lepton p_{T} [GeV]; Events", varBinN1, varBin1);  
    m_hSvc.create1D("lepPt",    "; Pt of lept [GeV]; Events", 100, 25, 525);
    m_hSvc.create1D("lepPt_genm",    "; Pt of lept [GeV]; Events", 100, 25, 525);
    m_hSvc.create1D("lepPt_fromb",    "; Pt of lept [GeV]; Events", 100, 25, 525);
    m_hSvc.create1D("lepEta", "; lepton #eta ; Events", 20, -2.5, 2.5);
    m_hSvc.create1D("lepPhi", "; lepton #phi [rd] ; Events", 32, -3.2, 3.2);
    m_hSvc.create1DVar("leadJetPt", "; leading Jet p_{T} [GeV]; Events", varBinN1, varBin1);  
    m_hSvc.create1D("nJets", "; number of jets ; Events", 10, -0.5, 9.5);

    m_hSvc.create1D("W_Hadronic", "; M_{j_{1}j_{2}} [GeV] ; Events", 1000, 0, 1000);
    m_hSvc.create1D("T_Hadronic", "; M_{j_{1}j_{2}j_{bh}} [GeV] ; Events", 1000, 0, 1000);
    m_hSvc.create1D("T_Leptonic", "; M_{j_{bl}l#nu} [GeV] ; Events", 1000, 0, 1000);
    m_hSvc.create1D("PT_Diff", "; #Delta p_{T} [GeV] ; Events", 1000, -500, 500);

    m_hSvc.create1DVar("leadbJetPt", "; leading b-jet p_{T} [GeV]; Events", varBinN1, varBin1);
    m_hSvc.create1DVar("leadTrkbJetPt", "; leading b-jet p_{T} [GeV]; Events", varBinN1, varBin1);
    m_hSvc.create1D("nBtagJets", "; number of b-tagged jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nTrkBtagJets", "; number of b-tagged track jets ; Events", 10, 0, 10);  

    m_hSvc.create1DVar("MET", "; missing E_{T} [GeV]; Events", varBinN1, varBin1);
    m_hSvc.create1D("MET_phi", "; missing E_{T} #phi [rd] ; Events", 32, -3.2, 3.2);

    m_hSvc.create1D("wpt", "; W p_{T} [GeV] ; Events", 100, 0, 1000);

    m_hSvc.create1D("DeltaPhiLargeJetLepton", "; #Delta#phi_{J,l} ; Events", 32, 0, 3.2);
    m_hSvc.create1D("DeltaRLargeJetSmallJet", "; #DeltaR_{J,j} ; Events", 100, 0, 3);

    m_hSvc.create1D("DeltaR_Leading", "; #Delta R Between Leading Jet and Track", 100, 0, 3.2);
    m_hSvc.create1D("DeltaR_Inclusive", "; #Delta R Between Inclusive Jet and Track", 100, 0, 3.2);

    m_hSvc.create1DVar("bjmatched_jet_pt", "; p_{T} of R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);  
    m_hSvc.create1DVar("bjmatched_btagged_jet_pt", "; p_{T} of R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);
    m_hSvc.create1DVar("Nonbjmatched_jet_pt", "; p_{T} of R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);   
    m_hSvc.create1DVar("Nonbjmatched_btagged_jet_pt", "; p_{T} of R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);

    m_hSvc.create1DVar("ltop_pt_jet",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_jet",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_jet_pt",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_jet_passCaloBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_jet_passCaloBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_jet_pt_passCaloBtag",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_tjet",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_tjet",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_tjet_pt",    "; Pt of b from leptonic top matched jet [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("ltop_pt_tjet_passTrackBtag",    "; Pt of truth leptonic top [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("b_from_ltop_pt_tjet_passTrackBtag",    "; Pt of truth b from leptonic top [GeV]; Events", varBinN9, varBin9);
    m_hSvc.create1DVar("b_from_ltop_matched_tjet_pt_passTrackBtag",    "; Pt of b from leptonic top matched track jet [GeV]; Events", varBinN9, varBin9);
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

    m_hSvc.create1D("trackBtagBDT", "; track b-tag BDT ; Events", 100, -1, 1);
    m_hSvc.create1D("caloBtagBDT", "; calo b-tag BDT ; Events", 100, -1, 1);

    char name[200];
    m_hSvc.create1DP("Profile_DelR", "; #Delta R; Events", 31, 0.2, 0.5, 0., 1.); 
    for(int idr = 20; idr <=50; idr++)
    {
      sprintf(name, "bjmatched_jet_pt_DelR%i",idr);
      m_hSvc.create1DVar(name, "; p_{T} of R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);
      //    sprintf(name, "Profile_DelR%i",idr);
      //    m_hSvc.create1DP(name, "; #Delta R; Events", 31, 0.2, 0.5, 0., 1.);
      sprintf(name, "bjmatched_btagged_jet_pt_DelR%i",idr);
      m_hSvc.create1DVar(name, "; p_{T} of R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);
      sprintf(name, "Nonbjmatched_jet_pt_DelR%i",idr);
      m_hSvc.create1DVar(name, "; p_{T} of R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);
      sprintf(name, "Nonbjmatched_btagged_jet_pt_DelR%i",idr);  
      m_hSvc.create1DVar(name, "; p_{T} of R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);
    } // for(int idr = 20; idr <=50; idr++)

    m_hSvc.create1D("mwt", "; W transverse mass [GeV]; Events", 20, 0, 200);
    m_hSvc.create1D("mu", "; <mu>; Events", 100, 0, 100);
    m_hSvc.create1D("mu_original", "; <mu_origianl>; Events", 100, 0, 100);
    m_hSvc.create1D("vtxz", ";Z position of truth primary vertex; Events", 40, -400, 400);
    m_hSvc.create1D("npv", "; npv; Events", 50, 0, 50);

    m_hSvc.create1D("closejl_minDeltaR", "; min #Delta R(lep, jet); Events", 50, 0, 5);
    m_hSvc.create1DVar("closejl_pt", "; Pt of closest jet to lep [GeV]; Events", varBinN1, varBin1);

    m_hSvc.create1D("weight_leptSF", "; QCD weights; Events", 200, 0, 2);
    m_hSvc.create1D("weight", "; QCD weights; Events", 2000, -100, 100);
    m_hSvc.create2D("weight_leptPt", ";lept Pt (GeV); QCD weights",100, 25, 525, 2000, -100, 100);

    m_hSvc.create1D("jet0_m", "; mass of the leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
    m_hSvc.create1D("jet1_m", "; mass of the sub-leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
    m_hSvc.create1D("jet2_m", "; mass of the third leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
    m_hSvc.create1D("jet3_m", "; mass of the fourth leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
    m_hSvc.create1D("jet4_m", "; mass of the fifth leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
    m_hSvc.create1D("jet5_m", "; mass of the sixth leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 

    m_hSvc.create1DVar("jet0_pt", "; p_{T} of the leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);
    m_hSvc.create1DVar("jet1_pt", "; p_{T} of the sub-leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
    m_hSvc.create1DVar("jet2_pt", "; p_{T} of the third leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
    m_hSvc.create1DVar("jet3_pt", "; p_{T} of the fourth leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
    m_hSvc.create1DVar("jet4_pt", "; p_{T} of the fifth leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
    m_hSvc.create1DVar("jet5_pt", "; p_{T} of the sixth leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 

    m_hSvc.create1D("jet0_eta", "; #eta of the leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
    m_hSvc.create1D("jet1_eta", "; #eta of the sub-leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
    m_hSvc.create1D("jet2_eta", "; #eta of the third leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
    m_hSvc.create1D("jet3_eta", "; #eta of the fourth leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
    m_hSvc.create1D("jet4_eta", "; #eta of the fifth leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
    m_hSvc.create1D("jet5_eta", "; #eta of the sixth leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 

    m_hSvc.create1D("jet0_phi", "; #phi of the leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
    m_hSvc.create1D("jet1_phi", "; #phi of the sub-leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
    m_hSvc.create1D("jet2_phi", "; #phi of the third leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
    m_hSvc.create1D("jet3_phi", "; #phi of the fourth leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
    m_hSvc.create1D("jet4_phi", "; #phi of the fifth leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
    m_hSvc.create1D("jet5_phi", "; #phi of the sixth leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 

    m_hSvc.create1DVar("tjet0_pt", "; p_{T} of the leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3);
    m_hSvc.create1DVar("tjet1_pt", "; p_{T} of the sub-leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
    m_hSvc.create1DVar("tjet2_pt", "; p_{T} of the third leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
    m_hSvc.create1DVar("tjet3_pt", "; p_{T} of the fourth leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
    m_hSvc.create1DVar("tjet4_pt", "; p_{T} of the fifth leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
    m_hSvc.create1DVar("tjet5_pt", "; p_{T} of the sixth leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 

    std::string bdtlh[2] = {"BDTlow","BDThigh"};
    if (m_boosted) {
      m_hSvc.create1DVar("closeJetPt", "; selected Jet p_{T} [GeV] ; Events", varBinN1, varBin1);
      m_hSvc.create1DVar("largeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
      m_hSvc.create1DVar("largeJetPt_passSTT", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
      m_hSvc.create1DVar("TopMatchedLargeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
      m_hSvc.create1DVar("TopMatchedLargeJetPt_passSTT", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
      m_hSvc.create1DVar("allLargeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
      m_hSvc.create1DVar("allLargeJetPt_passSTT", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
      m_hSvc.create1DVar("drLargeJetTruthTop", "; #Delta_{lj,truth-top} ; Events", 100,0,1);
      m_hSvc.create1D("largeJetC2", "; C2 ; Events", 100,0, 1);
      m_hSvc.create1D("largeJetD2", "; D2 ; Events", 100,0, 25);
      for(int j=0;j<2;j++){
        m_hSvc.create1DVar(Form("largeJetPt_%s",bdtlh[j].c_str()), "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
        m_hSvc.create1D(Form("largeJetM_%s",bdtlh[j].c_str()), "; large jet mass [GeV] ; Events", 30, 0, 300);
        m_hSvc.create1D(Form("largeJetEta_%s",bdtlh[j].c_str()), "; large jet #eta ; Events", 20, -2., 2.);
        m_hSvc.create1D(Form("largeJetPhi_%s",bdtlh[j].c_str()), "; large jet #phi [rd] ; Events", 32, -3.2, 3.2);
        m_hSvc.create1D(Form("largeJetSd12_%s",bdtlh[j].c_str()), "; large jet #sqrt{d_{12}} [GeV] ; Events", 20, 0, 200);
      }
      m_hSvc.create1D("largeJetSd12", "; large jet #sqrt{d_{12}} [GeV] ; Events", 20, 0, 200);
      m_hSvc.create1DVar("mtlep_boo", "; leptonic top mass [GeV] ; Events", varBinN4, varBin4);
      std::string bdt_thr[3] = {"loose","middle","tight"};
      for(int i=0;i<3;i++){
        m_hSvc.create1DVar(Form("largeJetPt_passBDT_%s",bdt_thr[i].c_str()), "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
        m_hSvc.create1DVar(Form("TopMatchedLargeJetPt_passBDT_%s",bdt_thr[i].c_str()), "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
        m_hSvc.create1DVar(Form("allLargeJetPt_passBDT_%d",bdt_thr[i].c_str()), "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
      }
    } else {
      m_hSvc.create1DVar("mtlep_res", "; leptonic top mass [GeV]; Events", varBinN4, varBin4);
      m_hSvc.create1DVar("mthad_res", "; hadronic top mass [GeV]; Events", varBinN4, varBin4);
      m_hSvc.create1DVar("mwhad_res", "; hadronic W boson mass [GeV]; Events", 40, 0, 400);
      m_hSvc.create1D("chi2", "; log(#chi^{2}) ; Events", 50, -3, 7);
    }

    m_hSvc.create1DVar("mtt", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
    m_hSvc.create1DVar("trueMtt", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);

    m_hSvc.create1D("largeJet_tau32", "; #tau_{32}; Events", 20, 0, 1);
    m_hSvc.create1D("largeJet_tau21", "; #tau_{21}; Events", 20, 0, 1);
    for(int j=0;j<2;j++){
      m_hSvc.create1D(Form("largeJet_tau32_wta_%s",bdtlh[j].c_str()), "; #tau_{32} wta; Events", 20, 0, 1);
      m_hSvc.create1D(Form("largeJet_tau21_wta_%s",bdtlh[j].c_str()), "; #tau_{21} wta; Events", 20, 0, 1);
    }

    m_hSvc.create1D("BDT_TOPtag_ljet_pt800_1000",    "; BDT TOP tag; ", 100, -1, 1);
    m_hSvc.create1D("BDT_TOPtag_ljet_pt500",    "; BDT TOP tag; ", 100, -1, 1);
    for(int lpt=0;lpt<4;lpt++){
      m_hSvc.create1D(Form("BDT_TOPtag_ljet_pt%d_%d",lpt*500,500+lpt*500),    "; BDT TOP tag; ", 100, -1, 1);
      m_hSvc.create1D(Form("BDT_TOPtag_allljet_pt%d_%d",lpt*500,500+lpt*500),    "; BDT TOP tag; ", 100, -1, 1);
      m_hSvc.create1D(Form("BDT_TOPtag_passljet_pt%d_%d",lpt*500,500+lpt*500),    "; BDT TOP tag; ", 100, -1, 1);
    }

    m_hSvc.create1D("minDrLjetTruthID",    "; id; ", 20, -10, 10);
    m_hSvc.create1D("minDrLjetTruth",    "; #DeltaR; ", 30, 0, 3);
  }

AnaTtresMVATopTag::~AnaTtresMVATopTag() {
}

void AnaTtresMVATopTag::run(const Event &evt, double weight, const std::string &s, int is2016run, float bdtthr) {

  // check channel
  //
  char name[200];

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  HistogramService *h = &m_hSvc;
  std::string suffix = s;
  for (int ij=0; ij < evt.largeJet().size(); ++ij) {
    const TLorentzVector &lj = evt.largeJet()[ij].mom();
    for(int lpt=0;lpt<4;lpt++){
      if(lj.Perp()*1e-3>500*lpt && lj.Perp()*1e-3<500+500*lpt){
        h->h1D(Form("BDT_TOPtag_allljet_pt%d_%d",500*lpt,500+500*lpt),"",suffix)->Fill(evt.largeJet()[ij].BDT_TOPtag(),weight);
      }
    }
  }

  if (m_boosted){
    if (!(evt.passes("bmujets_nobtag_notoptag") || evt.passes("bejets_nobtag_notoptag"))) 
      return;
  }

  if (!m_boosted){
    if (!(evt.passes("rejets_2016") || evt.passes("rmujets_2016") || evt.passes("rmujetsmet110")))
      return;
  }

  if (!m_boosted)	if(evt.jet().size()<4)	return;


  //HistogramService *h = &m_hSvc;

  bool trig1(0); 
  bool trig2(0); 
  bool trig3(0);
  bool trig4(0);
  bool trig5(0);
  bool trig6(0);

  bool isTight = false;

  TLorentzVector l;
  float mWt=0., MET=0., WPt=0.;
  if (m_electron) {
    l = evt.electron()[0].mom();

    isTight = evt.electron()[0].isTightPP();

    //Electron trigers
    trig1 = evt.electron()[0].HLT_e24_lhmedium_L1EM20VH();  
    trig2 = evt.electron()[0].HLT_e60_lhmedium();
    trig3 = evt.electron()[0].HLT_e120_lhloose();
    trig4 = evt.electron()[0].HLT_e26_lhtight_nod0_ivarloose();
    trig5 = evt.electron()[0].HLT_e60_lhmedium_nod0();
    trig6 = evt.electron()[0].HLT_e140_lhloose_nod0();

    bool trig_2015 = trig1 || trig2 || trig3; 
    bool trig_2016 = trig4 || trig5 || trig6;

    if(is2016run){
      if(!trig_2016) return;
    }else{
      if(!trig_2015) return;
    }

    if(!isTight) return;

  } else {

    l = evt.muon()[0].mom();
    TLorentzVector met = evt.met();
    mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3;
    MET = evt.met().Perp()*1e-3;
    WPt = sqrt(pow(MET,2) + pow(l.Perp()*1e-3,2));

    isTight = evt.muon()[0].isTight();

    //Muon trigers
    trig1 = evt.muon()[0].HLT_mu20_iloose_L1MU15();
    trig2 = evt.muon()[0].HLT_mu26_ivarmedium() && evt.trigger("HLT_mu26_ivarmedium");
    trig3 = evt.muon()[0].HLT_mu50() && evt.trigger("HLT_mu50");
    //MET triggers
    trig4 = evt.trigger("HLT_xe110_mht_L1XE60");

    bool trig_muon = trig2 || trig3 ;
    bool trig_met = trig4;

    if(!trig_muon) return;

    if(!isTight) return;

  }//m_electron

  //selection for dR_Jj dphi_Jl
  size_t close_idx = 0;
  for (; close_idx < evt.jet().size(); ++close_idx)
    if (evt.jet()[close_idx].closeToLepton())
      break;
  const TLorentzVector &sj = evt.jet()[close_idx].mom();
  bool pass_ljet_dr=false;
  size_t ljet_idx = 0;
  if(evt.largeJet().size()==0) return;
  for (; ljet_idx < evt.largeJet().size(); ++ljet_idx) {
    const TLorentzVector &lj = evt.largeJet()[ljet_idx].mom();
    float delta_phi_large_jet_lepton = l.DeltaPhi(lj);
    float delta_r_large_jet_small_jet = lj.DeltaR(sj);
    bool good_def = evt.largeJet()[ljet_idx].good();//default
    bool good_bdt = (evt.largeJet()[ljet_idx].BDT_TOPtag()>bdtthr)? true : false;
    if(lj.Perp()*1e-3>300 && fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) break;
  }
  if(ljet_idx==evt.largeJet().size()) return;

  TLorentzVector lgen;

  if(fabs(evt.MC_Wdecay1_from_t_pdgId()) == 11 || fabs(evt.MC_Wdecay1_from_t_pdgId()) == 13 || fabs(evt.MC_Wdecay1_from_t_pdgId()) == 15)
    lgen = evt.MC_Wdecay1_from_t();

  else if(fabs(evt.MC_Wdecay2_from_t_pdgId()) == 11 || fabs(evt.MC_Wdecay2_from_t_pdgId()) == 13 || fabs(evt.MC_Wdecay2_from_t_pdgId()) == 15)
    lgen = evt.MC_Wdecay2_from_t();

  else if(fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 11 || fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 13  || fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 15)
    lgen = evt.MC_Wdecay1_from_tbar();

  else if(fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 11 || fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 13 || fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 15)
    lgen = evt.MC_Wdecay2_from_tbar();

  bool isGenM = l.DeltaR(lgen) < 0.2;

  TLorentzVector bgen = evt.MC_b_from_t();
  TLorentzVector b_bar_gen = evt.MC_b_from_tbar();

  bool isFromB = l.DeltaR(bgen) < 0.5 || l.DeltaR(b_bar_gen) < 0.5;

  float sd0(99);

  if (m_electron) {
    isTight = evt.electron()[0].isTightPP();
    sd0    = evt.electron()[0].sd0();
  } else{
    isTight = evt.muon()[0].isTight();
    sd0    = evt.muon()[0].sd0();
  }//m_electron

  //std::string suffix = s;

  h->h1D("weight_leptSF", "", suffix)->Fill(evt.weight_leptonSF());
  h->h1D("weight", "", suffix)->Fill(weight);
  h->h2D("weight_leptPt", "", suffix)->Fill(l.Perp()*1e-3, weight);

  h->h1D("lepPt", "", suffix)->Fill(l.Perp()*1e-3, weight);
  if(isGenM)
    h->h1D("lepPt_genm", "", suffix)->Fill(l.Perp()*1e-3, weight);
  else if(isFromB)
    h->h1D("lepPt_fromb", "", suffix)->Fill(l.Perp()*1e-3, weight);
  h->h1D("lepPt_effBins", "", suffix)->Fill(l.Perp()*1e-3, weight);
  h->h1D("lepEta", "", suffix)->Fill(l.Eta(), weight);
  h->h1D("lepPhi", "", suffix)->Fill(l.Phi(), weight);

  h->h1D("MET_phi", "", suffix)->Fill(evt.met().Phi(), weight);

  const TLorentzVector &j = evt.jet()[0].mom();
  h->h1D("leadJetPt", "", suffix)->Fill(j.Perp()*1e-3, weight);

  // for now
  int nJets = evt.jet().size(); //njets 
  h->h1D("nJets", "", suffix)->Fill(nJets, weight);

  int w1h_pdgId = evt.MC_w1h_pdgId();
  int w2h_pdgId = evt.MC_w2h_pdgId();
  int w1l_pdgId = evt.MC_w1l_pdgId();
  int w2l_pdgId = evt.MC_w2l_pdgId();
  std::cout << "w1h_pdgId=" << w1h_pdgId << std::endl;
  std::cout << "w2h_pdgId=" << w2h_pdgId << std::endl;
  std::cout << "w1l_pdgId=" << w1l_pdgId << std::endl;
  std::cout << "w2l_pdgId=" << w2l_pdgId << std::endl;
  int Wdecay1_from_t_pdgId = evt.MC_Wdecay1_from_t_pdgId();
  int Wdecay2_from_t_pdgId = evt.MC_Wdecay2_from_t_pdgId();
  int Wdecay1_from_tbar_pdgId = evt.MC_Wdecay1_from_tbar_pdgId();
  int Wdecay2_from_tbar_pdgId = evt.MC_Wdecay2_from_tbar_pdgId();
  std::cout << "Wdecay1_from_t_pdgId=" << Wdecay1_from_t_pdgId << std::endl;
  std::cout << "Wdecay2_from_t_pdgId=" << Wdecay2_from_t_pdgId << std::endl;
  std::cout << "Wdecay1_from_tbar_pdgId=" << Wdecay1_from_tbar_pdgId << std::endl;
  std::cout << "Wdecay2_from_tbar_pdgId=" << Wdecay2_from_tbar_pdgId << std::endl;
  bool isLeptonicT = (abs(Wdecay1_from_t_pdgId)<7)? true : false;
  bool isLeptonicTbar = (abs(Wdecay1_from_tbar_pdgId)<7)? true : false;
  std::cout << "isLeptonic top/tbar=" << isLeptonicT << "/" << isLeptonicTbar << std::endl;

  TLorentzVector b_from_t,b_from_tbar,t,tbar;
  b_from_t = evt.MC_b_from_t();
  b_from_tbar = evt.MC_b_from_tbar();
  t = evt.MC_t();
  tbar = evt.MC_tbar();
  float t_pt = t.Perp()*1e-3, t_eta = t.Eta(), t_phi = t.Phi();
  float tbar_pt = tbar.Perp()*1e-3, tbar_eta = tbar.Eta(), tbar_phi = tbar.Phi();
  float b_from_t_pt = b_from_t.Perp()*1e-3, b_from_t_eta = b_from_t.Eta(), b_from_t_phi = b_from_t.Phi();
  float b_from_tbar_pt = b_from_tbar.Perp()*1e-3, b_from_tbar_eta = b_from_tbar.Eta(), b_from_tbar_phi = b_from_tbar.Phi();

  //truth matching
  int bt_sj_idx=-1, btbar_sj_idx=-1;
  float min_dr_sj_bt = 999., min_dr_sj_btbar = 999.;
  std::cout << "Nubmer of small jet is " << evt.jet().size() << std::endl;
  for (size_t bidx = 0; bidx < evt.jet().size(); ++bidx){
    float sj_eta = evt.jet()[bidx].mom().Eta();
    float sj_phi = evt.jet()[bidx].mom().Phi();
    float drt = TuDoAtlas::calc_delta_r(sj_eta,b_from_t_eta,sj_phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(sj_eta,b_from_tbar_eta,sj_phi,b_from_tbar_phi);
    if(drt<0.2 && drt<min_dr_sj_bt) {
      min_dr_sj_bt = drt;
      bt_sj_idx=bidx;
    }
    if(drtbar<0.2 && drtbar<min_dr_sj_btbar) {
      min_dr_sj_btbar = drtbar;
      btbar_sj_idx=bidx;
    }
  }

  int bt_tj_idx=-1, btbar_tj_idx=-1;
  float min_dr_tj_bt = 999., min_dr_tj_btbar = 999.;
  std::cout << "Nubmer of small tjet is " << evt.tjet().size() << std::endl;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
    float tj_eta = evt.tjet()[bidx].mom().Eta();
    float tj_phi = evt.tjet()[bidx].mom().Phi();
    float drt = TuDoAtlas::calc_delta_r(tj_eta,b_from_t_eta,tj_phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(tj_eta,b_from_tbar_eta,tj_phi,b_from_tbar_phi);
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
  std::cout << "btbar matched small jet id=" << btbar_sj_idx << std::endl;
  if(bt_sj_idx!=-1 && bt_tj_idx!=-1){
    h->h1D("trackBtagBDT", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mv2c10());
    h->h1D("trackBtagBDT", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mv2c10());
    h->h1D("caloBtagBDT", "", suffix)->Fill(evt.jet()[bt_sj_idx].mv2c20());
    h->h1D("caloBtagBDT", "", suffix)->Fill(evt.jet()[bt_sj_idx].mv2c20());
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_jet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_ltop_pt_jet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_ltop_matched_jet_pt", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[bt_sj_idx].btag_mv2c20_70()){
        std::cout << "btagged!!" << std::endl;
        h->h1D("ltop_pt_jet_passCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
    }
    else{//hadronic top
      h->h1D("htop_pt_jet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_jet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_htop_matched_jet_pt", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[bt_sj_idx].btag_mv2c20_70()){
        std::cout << "btagged!!" << std::endl;
        h->h1D("htop_pt_jet_passCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
    }
  }
  if(btbar_sj_idx!=-1 && btbar_tj_idx!=-1){
    if(isLeptonicTbar){//leptonic top
      h->h1D("ltop_pt_jet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_ltop_pt_jet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_ltop_matched_jet_pt", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[btbar_sj_idx].btag_mv2c20_70()){
        std::cout << "btagged!!" << std::endl;
        h->h1D("ltop_pt_jet_passCaloBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_ltop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      }
    }
    else{//hadronic top
      h->h1D("htop_pt_jet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_htop_pt_jet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_htop_matched_jet_pt", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[btbar_sj_idx].btag_mv2c20_70()){
        std::cout << "btagged!!" << std::endl;
        h->h1D("htop_pt_jet_passCaloBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_htop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      }
    }
  }


  //track tjet b-tag efficiency
  std::cout << "bt matched small tjet id=" << bt_tj_idx << std::endl;
  std::cout << "btbar matched small tjet id=" << btbar_tj_idx << std::endl;
  if(bt_tj_idx!=-1 && bt_sj_idx!=-1){
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_tjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_ltop_pt_tjet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_ltop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      if(evt.tjet()[bt_tj_idx].btag_mv2c10_70_trk()){
        std::cout << "btagged!!" << std::endl;
        h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
    }
    else{//hadronic top
      h->h1D("htop_pt_tjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      if(evt.tjet()[bt_tj_idx].btag_mv2c10_70_trk()){
        std::cout << "btagged!!" << std::endl;
        h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
    }
  }
  if(btbar_tj_idx!=-1 && btbar_sj_idx!=-1){
    if(isLeptonicTbar){//leptonic top
      h->h1D("ltop_pt_tjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_ltop_pt_tjet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_ltop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      if(evt.tjet()[btbar_tj_idx].btag_mv2c10_70_trk()){
        std::cout << "btagged!!" << std::endl;
        h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
    }
    else{//hadronic top
      h->h1D("htop_pt_tjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      if(evt.tjet()[btbar_tj_idx].btag_mv2c10_70_trk()){
        std::cout << "btagged!!" << std::endl;
        h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
    }
  }

  /*float ltop_pt = evt.MC_tl().Perp()*1e-3;
  float ltop_eta = evt.MC_tl().Eta();
  float ltop_phi = evt.MC_tl().Phi();
  //std::cout << "truth ltop eta/phi=" << ltop_eta << "/" << ltop_phi << std::endl;
  TLorentzVector b_from_lt = evt.MC_bl();
  float b_from_lt_pt = b_from_lt.Perp()*1e-3;
  float b_from_lt_eta = b_from_lt.Eta();
  float b_from_lt_phi = b_from_lt.Phi();
  std::cout << "truth b from ltop pt/eta/phi=" << b_from_lt_pt << "/" << b_from_lt_eta << "/" << b_from_lt_phi << std::endl;
  float htop_pt = evt.MC_th().Perp()*1e-3;
  float htop_eta = evt.MC_th().Eta();
  float htop_phi = evt.MC_th().Phi();
  //std::cout << "truth htop eta/phi=" << htop_eta << "/" << htop_phi << std::endl;
  TLorentzVector b_from_ht = evt.MC_bh();
  float b_from_ht_pt = b_from_ht.Perp()*1e-3;
  float b_from_ht_eta = b_from_ht.Eta();
  float b_from_ht_phi = b_from_ht.Phi();
  std::cout << "truth b from htop pt/eta/phi=" << b_from_ht_pt << "/" << b_from_ht_eta << "/" << b_from_ht_phi << std::endl;
  int nBtagged = 0; //nB-tagged jets 
  int lb_sj_idx=-1, hb_sj_idx=-1;
  float min_dr_sj_lb = 999., min_dr_sj_hb = 999.;
  std::cout << "Nubmer of small jet is " << evt.tjet().size() << std::endl;
  for (size_t bidx = 0; bidx < evt.jet().size(); ++bidx){
    float sj_eta = evt.jet()[bidx].mom().Eta();
    float sj_phi = evt.jet()[bidx].mom().Phi();
    float hdr = TuDoAtlas::calc_delta_r(sj_eta,b_from_ht_eta,sj_phi,b_from_ht_phi);
    float ldr = TuDoAtlas::calc_delta_r(sj_eta,b_from_lt_eta,sj_phi,b_from_lt_phi);
    if(ldr<0.2 && ldr<min_dr_sj_lb) {
      min_dr_sj_lb = ldr;
      lb_sj_idx=bidx;
    }
    if(hdr<0.2 && hdr<min_dr_sj_hb) {
      min_dr_sj_hb = hdr;
      hb_sj_idx=bidx;
    }
    if (evt.jet()[bidx].btag_mv2c20_70()){
      if(nBtagged==0)h->h1D("leadbJetPt", "", suffix)->Fill(evt.jet()[bidx].mom().Perp()*1e-3, weight);
      nBtagged += 1;
    }
  }
  h->h1D("nBtagJets", "", suffix)->Fill(nBtagged, weight);
  //std::cout << "b matched small jet id=" << lb_sj_idx << std::endl;
  if(lb_sj_idx!=-1){
    std::cout << "lb matched small jet pt/eta/phi=" 
      << evt.jet()[lb_sj_idx].mom().Perp()*1e-3 << "/" << evt.jet()[lb_sj_idx].mom().Eta() << "/" << evt.jet()[lb_sj_idx].mom().Phi() << std::endl; 
    std::cout << "caloBtag mv2c20=" << evt.jet()[lb_sj_idx].mv2c20() << std::endl;
    h->h1D("ltop_pt_jet", "", suffix)->Fill(ltop_pt);
    h->h1D("b_from_ltop_pt_jet", "", suffix)->Fill(b_from_lt_pt);
    h->h1D("b_from_ltop_matched_jet_pt", "", suffix)->Fill(evt.jet()[lb_sj_idx].mom().Perp()*1e-3);
    if(evt.jet()[lb_sj_idx].btag_mv2c20_70()){
      std::cout << "btagged!!" << std::endl;
      h->h1D("ltop_pt_jet_passCaloBtag", "", suffix)->Fill(ltop_pt);
      h->h1D("b_from_ltop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_lt_pt);
      h->h1D("b_from_ltop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[lb_sj_idx].mom().Perp()*1e-3);
    }
  }
  if(hb_sj_idx!=-1){
    std::cout << "hb matched small jet pt/eta/phi=" 
      << evt.jet()[hb_sj_idx].mom().Perp()*1e-3 << "/" << evt.jet()[hb_sj_idx].mom().Eta() << "/" << evt.jet()[hb_sj_idx].mom().Phi() << std::endl; 
    std::cout << "caloBtag mv2c20=" << evt.jet()[hb_sj_idx].mv2c20() << std::endl;
    h->h1D("htop_pt_jet", "", suffix)->Fill(htop_pt);
    h->h1D("b_from_htop_pt_jet", "", suffix)->Fill(b_from_ht_pt);
    h->h1D("b_from_htop_matched_jet_pt", "", suffix)->Fill(evt.jet()[hb_sj_idx].mom().Perp()*1e-3);
    if(evt.jet()[hb_sj_idx].btag_mv2c20_60()){
      std::cout << "btagged!!" << std::endl;
      h->h1D("htop_pt_jet_passCaloBtag", "", suffix)->Fill(htop_pt);
      h->h1D("b_from_htop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_ht_pt);
      h->h1D("b_from_htop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[hb_sj_idx].mom().Perp()*1e-3);
    }
  }*/

  std::vector<float> tjetPt_vector;
  /*
  int nTrkBtagged = 0; //nTrkB-tagged jets 
  int lb_tj_idx=-1,hb_tj_idx=-1;
  float min_dr_tj_lb = 999.,min_dr_tj_hb=999.;
  std::cout << "Nubmer of track jet is " << evt.tjet().size() << std::endl;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx) {
    float tj_eta = evt.tjet()[bidx].mom().Eta();
    float tj_phi = evt.tjet()[bidx].mom().Phi();
    float ldr = TuDoAtlas::calc_delta_r(tj_eta,b_from_lt_eta,tj_phi,b_from_lt_phi);
    float hdr = TuDoAtlas::calc_delta_r(tj_eta,b_from_ht_eta,tj_phi,b_from_ht_phi);
    bool isBtagged = evt.tjet()[bidx].btag_mv2c10_70_trk();
    std::cout << "id/isBtagged=" << bidx << "/" << isBtagged << std::endl;
    if(ldr<0.2 && ldr<min_dr_tj_lb) {
      min_dr_tj_lb = ldr;
      lb_tj_idx=bidx;
    }
    if(hdr<0.2 && hdr<min_dr_tj_hb) {
      min_dr_tj_hb = hdr;
      hb_tj_idx=bidx;
    }
    tjetPt_vector.push_back(evt.tjet()[bidx].mom().Perp());
    if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk()){
      if(nTrkBtagged==0)h->h1D("leadTrkbJetPt", "", suffix)->Fill(evt.tjet()[bidx].mom().Perp()*1e-3, weight);
      nTrkBtagged += 1;
    }
  }
  h->h1D("nTrkBtagJets", "", suffix)->Fill(nTrkBtagged, weight);
  std::cout << "lb_tj_idx=" << lb_tj_idx << std::endl;
  std::cout << "hb_tj_idx=" << hb_tj_idx << std::endl;
  if(lb_tj_idx!=-1 && evt.tjet()[lb_tj_idx].pass_trk()){
    std::cout << "lb matched track jet pt/eta/phi=" 
      << evt.tjet()[lb_tj_idx].mom().Perp()*1e-3 << "/" << evt.tjet()[lb_tj_idx].mom().Eta() << "/" << evt.tjet()[lb_tj_idx].mom().Phi() << std::endl; 
    std::cout << "trackBtag mv2c10=" << evt.tjet()[lb_tj_idx].mv2c10() << std::endl;
    h->h1D("ltop_pt_tjet", "", suffix)->Fill(ltop_pt);
    h->h1D("b_from_ltop_pt_tjet", "", suffix)->Fill(b_from_lt_pt);
    h->h1D("b_from_ltop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[lb_tj_idx].mom().Perp()*1e-3);
    if(evt.tjet()[lb_tj_idx].btag_mv2c10_70_trk()){
      std::cout << "btagged!!" << std::endl;
      h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(ltop_pt);
      h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_lt_pt);
      h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[lb_tj_idx].mom().Perp()*1e-3);
    }
  }
  if(hb_tj_idx!=-1 && evt.tjet()[hb_tj_idx].pass_trk()){
    std::cout << "hb matched track jet pt/eta/phi=" 
      << evt.tjet()[hb_tj_idx].mom().Perp()*1e-3 << "/" << evt.tjet()[hb_tj_idx].mom().Eta() << "/" << evt.tjet()[hb_tj_idx].mom().Phi() << std::endl; 
    std::cout << "trackBtag mv2c10=" << evt.tjet()[hb_tj_idx].mv2c10() << std::endl;
    h->h1D("htop_pt_tjet", "", suffix)->Fill(htop_pt);
    h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_ht_pt);
    h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[hb_tj_idx].mom().Perp()*1e-3);
    if(evt.tjet()[hb_tj_idx].btag_mv2c10_70_trk()){
      std::cout << "btagged!!" << std::endl;
      h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(htop_pt);
      h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_ht_pt);
      h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[hb_tj_idx].mom().Perp()*1e-3);
    }
  }
  */

  return;

  int maxNtjet = (evt.tjet().size()<6) ? evt.tjet().size() : 6;
  for (int i = 0; i < maxNtjet; ++i){ 
    std::stringstream ss;
    ss << i;
    std::string nameJet_pt = "tjet" + ss.str()+"_pt";
    h->h1D(nameJet_pt, "", suffix)->Fill(tjetPt_vector[i]*1e-3, weight);
  }

  // Jet kinematics  

  std::vector<float> jetMass_vector;
  std::vector<float> jetPt_vector;
  std::vector<float> jetEta_vector;
  std::vector<float> jetPhi_vector;

  jetMass_vector.resize(evt.jet().size());
  jetPt_vector.resize(evt.jet().size());
  jetEta_vector.resize(evt.jet().size());
  jetPhi_vector.resize(evt.jet().size());  

  size_t iJet = 0;
  for (; iJet < evt.jet().size(); ++iJet){
    const TLorentzVector &jet_p4 = evt.jet()[iJet].mom();
    jetMass_vector[iJet] =  evt.jet()[iJet].mom().M();
    jetPt_vector[iJet]   =  evt.jet()[iJet].mom().Pt();
    jetEta_vector[iJet]  =  evt.jet()[iJet].mom().Eta();
    jetPhi_vector[iJet]  =  evt.jet()[iJet].mom().Phi();    
  }//for

  int maxNjet = (evt.jet().size()<6) ? evt.jet().size() : 6;
  for (int i = 0; i < maxNjet; ++i){  
    std::stringstream ss;
    ss << i;

    std::string nameJet_m = "jet" + ss.str()+"_m";
    h->h1D(nameJet_m, "", suffix)->Fill(jetMass_vector[i]*1e-3, weight);

    std::string nameJet_pt = "jet" + ss.str()+"_pt";
    h->h1D(nameJet_pt, "", suffix)->Fill(jetPt_vector[i]*1e-3, weight);

    std::string nameJet_eta = "jet" + ss.str()+"_eta";
    h->h1D(nameJet_eta, "", suffix)->Fill(jetEta_vector[i], weight);

    std::string nameJet_phi = "jet" + ss.str()+"_phi";
    h->h1D(nameJet_phi, "", suffix)->Fill(jetPhi_vector[i], weight);    

  }//for

  float mtt = -1;

  //missing et
  h->h1D("MET", "", suffix)->Fill(evt.met().Perp()*1e-3, weight);

  //transverse W mass  
  h->h1D("mwt", "", suffix)->Fill(sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(l.Phi() - evt.met().Phi())))*1e-3, weight); 

  //transverse W pt  
  h->h1D("wpt", "", suffix)->Fill(WPt, weight); 

  h->h1D("yields", "", suffix)->Fill(1, weight);

  //mu
  h->h1D("mu", "", suffix)->Fill(evt.mu()*1.16, weight); 
  h->h1D("mu_original", "", suffix)->Fill(evt.mu_original(), weight); 

  //npv
  h->h1D("npv", "", suffix)->Fill(evt.npv(), weight);

  //z prosition of primary vertex
  h->h1D("vtxz", "", suffix)->Fill(evt.vtxz(), weight);

  //deltaR between lepton and the closest narrow jet
  float closejl_deltaR  = 99;
  float deltaR_tmp      = 99;
  int closejl_idx       = -1;

  size_t jet_idx = 0;
  for (; jet_idx < evt.jet().size(); ++jet_idx){    
    deltaR_tmp = 99;
    float dphi = l.DeltaPhi(evt.jet()[jet_idx].mom());
    float dy = l.Rapidity() - evt.jet()[jet_idx].mom().Rapidity();
    deltaR_tmp = std::sqrt(dy*dy + dphi*dphi);
    if (deltaR_tmp < closejl_deltaR){
      closejl_deltaR = deltaR_tmp;
      closejl_idx = jet_idx;
    }   
  }//for     

  float closejl_pt = -1;
  if (closejl_idx>0)    closejl_pt = evt.jet()[closejl_idx].mom().Perp();  
  h->h1D("closejl_minDeltaR", "", suffix)->Fill(closejl_deltaR, weight);
  h->h1D("closejl_minDeltaR_effBins", "", suffix)->Fill(closejl_deltaR, weight);
  h->h1D("closejl_pt", "", suffix)->Fill(closejl_pt*1e-3, weight);

  h->h1D("trueMtt", "", suffix)->Fill(evt.MC_ttbar_beforeFSR().M()*1e-3, weight);
  _tree_truemtt = evt.MC_ttbar_beforeFSR().M()*1e-3;
  if (m_boosted && (is2016run && (evt.passes("bejets_nobtag_notoptag") || evt.passes("bmujets_nobtag_notoptag") ))) {

    float truth_htop_eta = evt.MC_th().Eta();
    float truth_htop_phi = evt.MC_th().Phi();
    //std::cout << "truth htop eta/phi=" << truth_htop_eta << "/" << truth_htop_phi << std::endl;

    float truth_top_eta = evt.MC_t().Eta();
    float truth_top_phi = evt.MC_t().Phi();
    //std::cout << "truth top eta/phi=" << truth_htop_eta << "/" << truth_htop_phi << std::endl;

    size_t close_idx = 0;
    for (; close_idx < evt.jet().size(); ++close_idx)
      if (evt.jet()[close_idx].closeToLepton())
        break;
    const TLorentzVector &sj = evt.jet()[close_idx].mom();
    h->h1D("closeJetPt", "", suffix)->Fill(sj.Perp()*1e-3, weight);

    //std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
    for (int ij=0; ij < evt.largeJet().size(); ++ij) {
      const TLorentzVector &lj = evt.largeJet()[ij].mom();
      float ljet_eta = lj.Eta();
      float ljet_phi = lj.Phi();
      float deta = truth_htop_eta-ljet_eta;
      float dphi = acos(cos(truth_htop_phi-ljet_phi));
      float dr = deta*deta+dphi*dphi; 
      h->h1D("drLargeJetTruthTop", "", suffix)->Fill(dr, weight);
      bool good_def = evt.largeJet()[ij].good();//default
      bool good_bdt = (evt.largeJet()[ij].BDT_TOPtag()>bdtthr)? true : false;
      float delta_phi_large_jet_lepton = l.DeltaPhi(lj);
      float delta_r_large_jet_small_jet = lj.DeltaR(sj);
      bool pass_ljet = (lj.Perp()*1e-3>300 && fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) ?true:false;
      h->h1D("allLargeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      for(int lpt=0;lpt<4;lpt++){
        if(lj.Perp()*1e-3>500*lpt && lj.Perp()*1e-3<500+500*lpt){
          h->h1D(Form("BDT_TOPtag_allljet_pt%d_%d",500*lpt,500+500*lpt),"",suffix)->Fill(evt.largeJet()[ij].BDT_TOPtag(),weight);
          //if(pass_ljet) h->h1D(Form("BDT_TOPtag_passljet_pt%d_%d",500*lpt,500+500*lpt),"",suffix)->Fill(evt.largeJet()[ij].BDT_TOPtag(),weight);
        }
      }
      if(good_def) h->h1D("allLargeJetPt_passSTT", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      if(evt.largeJet()[ij].BDT_TOPtag()>-0.5) h->h1D("allLargeJetPt_passBDT_loose", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      if(evt.largeJet()[ij].BDT_TOPtag()>0) h->h1D("allLargeJetPt_passBDT_middle", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      if(evt.largeJet()[ij].BDT_TOPtag()>0.5) h->h1D("allLargeJetPt_passBDT_tight", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      //std::cout << "dr=" << dr << std::endl;
      if(dr>0.2) continue;
      //std::cout << "bdt score=" << evt.largeJet()[ij].BDT_TOPtag() << std::endl;
      h->h1D("TopMatchedLargeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      //std::cout << "truth matched top good_def/good_bdt=" << good_def << "/" << good_bdt << std::endl;
      if(good_def) h->h1D("TopMatchedLargeJetPt_passSTT", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      if(evt.largeJet()[ij].BDT_TOPtag()>-0.5) h->h1D("TopMatchedLargeJetPt_passBDT_loose", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      if(evt.largeJet()[ij].BDT_TOPtag()>0) h->h1D("TopMatchedLargeJetPt_passBDT_middle", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      if(evt.largeJet()[ij].BDT_TOPtag()>0.5) h->h1D("TopMatchedLargeJetPt_passBDT_tight", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    }


    size_t goodljet_idx = 0;
    for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx) {
      const TLorentzVector &lj = evt.largeJet()[goodljet_idx].mom();
      float delta_phi_large_jet_lepton = l.DeltaPhi(lj);
      float delta_r_large_jet_small_jet = lj.DeltaR(sj);
      bool good_def = evt.largeJet()[goodljet_idx].good();//default
      bool good_bdt = (evt.largeJet()[goodljet_idx].BDT_TOPtag()>bdtthr)? true : false;
      if(lj.Perp()*1e-3>300 && fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) break;
    }
    if (goodljet_idx==evt.largeJet().size()) return; 

    //goodljet_idx=0;
    const TLorentzVector &lj = evt.largeJet()[goodljet_idx].mom();
    float delta_phi_large_jet_lepton = l.DeltaPhi(lj);
    float delta_r_large_jet_small_jet = lj.DeltaR(sj);

    h->h1D("DeltaPhiLargeJetLepton","",suffix)->Fill(delta_phi_large_jet_lepton,weight);
    h->h1D("DeltaRLargeJetSmallJet","",suffix)->Fill(delta_r_large_jet_small_jet,weight);

    h->h1D("largeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    h->h1D("largeJetC2", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("C2"), weight);
    h->h1D("largeJetD2", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("D2"), weight);
    if(evt.largeJet()[goodljet_idx].good()) h->h1D("largeJetPt_passSTT", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    if(evt.largeJet()[goodljet_idx].BDT_TOPtag()>-0.5) h->h1D("largeJetPt_passBDT_loose", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    if(evt.largeJet()[goodljet_idx].BDT_TOPtag()>0) h->h1D("largeJetPt_passBDT_middle", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    if(evt.largeJet()[goodljet_idx].BDT_TOPtag()>0.5) h->h1D("largeJetPt_passBDT_tight", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    std::string bdtlh[2] = {"BDTlow","BDThigh"};
    int lh = (evt.largeJet()[goodljet_idx].BDT_TOPtag()>0)?1:0;
    h->h1D(Form("largeJetPt_%s",bdtlh[lh].c_str()), "", suffix)->Fill(lj.Perp()*1e-3, weight);
    h->h1D(Form("largeJetM_%s",bdtlh[lh].c_str()), "", suffix)->Fill(lj.M()*1e-3, weight);
    h->h1D(Form("largeJetEta_%s",bdtlh[lh].c_str()), "", suffix)->Fill(lj.Eta(), weight);
    h->h1D(Form("largeJetPhi_%s",bdtlh[lh].c_str()), "", suffix)->Fill(lj.Phi(), weight);
    h->h1D(Form("largeJetSd12_%s",bdtlh[lh].c_str()), "", suffix)->Fill(evt.largeJet()[goodljet_idx].split12()*1e-3, weight);

    h->h1D("largeJet_tau32", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau32"), weight);
    h->h1D(Form("largeJet_tau32_wta_%s",bdtlh[lh].c_str()), "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau32_wta"), weight);

    h->h1D("largeJet_tau21", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau21"), weight);
    h->h1D(Form("largeJet_tau21_wta_%s",bdtlh[lh].c_str()), "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau21_wta"), weight);
    if(lj.Perp()*1e-3>500)
      h->h1D("BDT_TOPtag_ljet_pt500","",suffix)->Fill(evt.largeJet()[goodljet_idx].BDT_TOPtag(),weight);
    if(lj.Perp()*1e-3>800 && lj.Perp()*1e-3<1000)
      h->h1D("BDT_TOPtag_ljet_pt800_1000","",suffix)->Fill(evt.largeJet()[goodljet_idx].BDT_TOPtag(),weight);
    for(int lpt=0;lpt<4;lpt++){
      if(lj.Perp()*1e-3>500*lpt && lj.Perp()*1e-3<500+500*lpt)
        h->h1D(Form("BDT_TOPtag_ljet_pt%d_%d",500*lpt,500+500*lpt),"",suffix)->Fill(evt.largeJet()[goodljet_idx].BDT_TOPtag(),weight);
    }

    if(lj.Perp()*1e-3>1000 && lj.Perp()*1e-3<1500){
      //std::cout << "ljet index=" << goodljet_idx << std::endl;
      //std::cout << "BDT/weight=" << evt.largeJet()[goodljet_idx].BDT_TOPtag() << "/" << weight << std::endl;
      //std::cout << "ljet pt/eta/phi=" << lj.Perp()*1e-3 << "/" << lj.Eta() << "/" << lj.Phi() << std::endl;
    }

    // recalc. mtt
    // lepton = l
    // large-R jet = hadronic top = lj
    // selected jet = leptonic top's b-jet = sj
    // neutrino px, py = met
    std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), true);   
    TLorentzVector nu(0,0,0,0);
    if (vec_nu.size() > 0) {
      nu = *(vec_nu[0]);
      for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
      vec_nu.clear();
    }

    if (evt.largeJet().size()!=0)    mtt = (lj+sj+nu+l).M();
    h->h1D("mtt", "", suffix)->Fill(mtt*1e-3, weight);
    h->h1D("mtlep_boo", "", suffix)->Fill((sj+nu+l).M()*1e-3, weight);

    //std::cout << "Number of truth particle is " << evt.truth().size() << std::endl;
    float min_dr=999.;
    int min_id=-1;
    for (uint32_t it=0; it < evt.truth().size(); ++it) {
      float px = evt.truth()[it].px();
      float py = evt.truth()[it].py();
      float pz = evt.truth()[it].pz();
      float m = evt.truth()[it].m();
      float e = evt.truth()[it].e();
      int id = evt.truth()[it].id();
      int status = evt.truth()[it].status();
      TLorentzVector particle;
      particle.SetPxPyPzE(px,py,pz,e);
      float pt = particle.Pt();
      float eta = particle.Eta();
      float phi = particle.Phi();
      if(pt<1) continue;
      if(abs(id)>6) continue;
      //std::cout << "particle id/status/pt/eta/phi=" << id << "/" << status << "/" << pt << "/" << eta << "/" << phi << std::endl;
      float ljet_pt = lj.Perp()*1e-3;
      float ljet_eta = lj.Eta();
      float ljet_phi = lj.Phi();
      float deta = eta - ljet_eta;
      float dphi = acos(cos(phi - ljet_phi));
      float dr = sqrt(deta*deta+dphi*dphi);
      if(dr<min_dr){
        min_dr = dr;
        min_id = id;
      }
    }
    h->h1D("minDrLjetTruthID", "", suffix)->Fill(min_id, weight);
    h->h1D("minDrLjetTruth", "", suffix)->Fill(min_dr, weight);

    _tree_mtt = mtt*1e-3;
    _tree_weight = weight;
    _tree_cat = -1;
    if (((is2016run && evt.passes("bejets_nobtag_notoptag")) )  && m_boosted && m_electron) _tree_cat = 0;
    if ((((is2016run && evt.passes("bmujets_nobtag_notoptag")))) && m_boosted && !m_electron) _tree_cat = 1;
    _tree_syst = s;
    h->m_tree->Fill();

  }else if (!m_boosted && ((is2016run && evt.passes("rejets_2016")) || ((is2016run && evt.passes("rmujets_2016")) || evt.passes("rmujetsmet110")))) {

    // inputs 
    // LEPTON --> TLorentzVector for your lepton
    // vjets -->  std::vector<TLorentzVector*> for the jets
    // vjets_btagged --> std::vector<bool> to say if the jets are btagged or not
    // met --> TLorentzVector for your MET

    // outputs, they will be filled by the TTBarLeptonJetsBuilder_chi2
    int  igj3, igj4; // index for the Whad
    int igb3, igb4; // index for the b's
    int  ign1;  // index for the neutrino (because chi2 can test both pz solution)
    double chi2ming1, chi2ming1H, chi2ming1L;

    for(double idr=20.0; idr<=50.0; idr++)
    {
      TLorentzVector bl = evt.MC_bl();
      TLorentzVector bh = evt.MC_bh();
      double N_BMatch(0.), N_BMatch_BTag(0.), N_nBMatch(0.), N_nBMatch_BTag(0.);
      for (size_t z = 0; z < evt.jet().size(); ++z)
      {
        TLorentzVector tmpAk4Jet = evt.jet()[z].mom();
        bool is_bmatch = tmpAk4Jet.DeltaR(bl) < 0.4 || tmpAk4Jet.DeltaR(bh) < 0.4;
        bool is_btag(false);

        for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx)
        {
          TLorentzVector tmpTJet = evt.tjet()[bidx].mom();
          if(tmpAk4Jet.DeltaR(tmpTJet) < idr/100.) {
            if(evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
              is_btag = true;
          } // if(tmpAk4Jet.DeltaR(tmpTJet) < idr/100.)
        } // for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx)

        //if(is_bmatch && is_btag) std::cout<<"btag found"<<std::endl;
        if(is_bmatch) N_BMatch++;
        if(!is_bmatch) N_nBMatch++;
        if(is_bmatch && is_btag) N_BMatch_BTag++;
        if(!is_bmatch && is_btag) N_nBMatch_BTag++;
      } // for (size_t z = 0; z < evt.jet().size(); ++z)

      if(N_BMatch == 0 || N_nBMatch == 0) continue;
      double e_b = N_BMatch_BTag/N_BMatch;
      double e_bbar = N_nBMatch_BTag/N_nBMatch;

      h->p1D("Profile_DelR", "", suffix)->Fill(idr/100., e_b*(1-e_bbar));
    } // or(double idr=20.0; idr<=50.0; idr++)


    TLorentzVector a_w1h( evt.MA_w1h() ), a_w2h( evt.MA_w2h() ), a_bh( evt.MA_bh() );
    TLorentzVector a_w1l( evt.MA_w1l() ), a_w2l( evt.MA_w2l() ), a_bl( evt.MA_bl() );
    h->h1D("W_Hadronic", "", suffix)->Fill((a_w1h + a_w2h).M()*1e-3);
    h->h1D("T_Hadronic", "", suffix)->Fill((a_w1h + a_w2h + a_bh).M()*1e-3 - (a_w1h + a_w2h).M()*1e-3);
    h->h1D("T_Leptonic", "", suffix)->Fill((a_w1l + a_w2l + a_bl).M()*1e-3);
    h->h1D("PT_Diff", "", suffix)->Fill(((a_w1h+a_w2h+a_bh).Pt() - (a_w1l+a_w2l+a_bl).Pt() )*1e-3);

    std::vector<TLorentzVector *> vjets;
    std::vector<bool> vjets_btagged;
    for (size_t z = 0; z < evt.jet().size(); ++z) {
      vjets.push_back(new TLorentzVector(0,0,0,0));
      vjets[z]->SetPtEtaPhiE(evt.jet()[z].mom().Perp(), evt.jet()[z].mom().Eta(), evt.jet()[z].mom().Phi(), evt.jet()[z].mom().E());
      // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BTagingxAODEDM
      // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks
      //vjets_btagged.push_back(evt.jet()[z].btag_mv2c20_60());

      TLorentzVector tmpAk4Jet = evt.jet()[z].mom();
      bool is_btagged(false), is_bmatched(false);
      std::vector<bool> is_btagged_delR(31, false);

      TLorentzVector bl = evt.MC_bl();
      TLorentzVector bh = evt.MC_bh();

      if(tmpAk4Jet.DeltaR(bl) < 0.4 || tmpAk4Jet.DeltaR(bh) < 0.4)
        is_bmatched = true;

      for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx)
      {
        TLorentzVector tmpTJet = evt.tjet()[bidx].mom();
        if(z==0 && bidx==0)
          h->h1D("DeltaR_Leading", "", suffix)->Fill(tmpAk4Jet.DeltaR(tmpTJet));
        h->h1D("DeltaR_Inclusive", "", suffix)->Fill(tmpAk4Jet.DeltaR(tmpTJet));
        if(tmpAk4Jet.DeltaR(tmpTJet) <= 0.4){ 
          if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
            is_btagged = true;
          // break;
        }

        for(double idr=20.0; idr<=50.0; idr++)
        {
          if(tmpAk4Jet.DeltaR(tmpTJet) < idr/100.){
            if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk()) 
              is_btagged_delR[int(idr)-20] = true;
            // break;
          }
          //std::cout<<"index = "<<int(idr)-20<<", UAngle = "<<idr/100.<<", angle = "<<tmpAk4Jet.DeltaR(tmpTJet)<<
          //", bool = "<<is_btagged_delR[int(idr)-20]<<std::endl;
        } // for(double idr=20.0; idr<=50.0; idr++)
      } // for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx)

      vjets_btagged.push_back(is_btagged);
      //std::cout<<"Jet pT = "<<tmpAk4Jet.Pt()<<std::endl;
      if(is_bmatched)
        h->h1D("bjmatched_jet_pt", "", suffix)->Fill(tmpAk4Jet.Pt()*1e-3);
      if(is_bmatched && is_btagged)
        h->h1D("bjmatched_btagged_jet_pt", "", suffix)->Fill(tmpAk4Jet.Pt()*1e-3);
      if(!is_bmatched)
        h->h1D("Nonbjmatched_jet_pt", "", suffix)->Fill(tmpAk4Jet.Pt()*1e-3);
      if(!is_bmatched && is_btagged)
        h->h1D("Nonbjmatched_btagged_jet_pt", "", suffix)->Fill(tmpAk4Jet.Pt()*1e-3);

      //char name[200];
      for(int idr = 20; idr <=50; idr++)
      {
        sprintf(name, "bjmatched_jet_pt_DelR%i",idr);
        if(is_bmatched)
          h->h1D(name, "", suffix)->Fill(tmpAk4Jet.Pt()*1e-3);

        sprintf(name, "bjmatched_btagged_jet_pt_DelR%i",idr);
        if(is_bmatched && is_btagged_delR[idr-20])
          h->h1D(name, "", suffix)->Fill(tmpAk4Jet.Pt()*1e-3);

        sprintf(name, "Nonbjmatched_jet_pt_DelR%i",idr);
        if(!is_bmatched)
          h->h1D(name, "", suffix)->Fill(tmpAk4Jet.Pt()*1e-3); 

        sprintf(name, "Nonbjmatched_btagged_jet_pt_DelR%i",idr);
        if(!is_bmatched && is_btagged_delR[idr-20])
          h->h1D(name, "", suffix)->Fill(tmpAk4Jet.Pt()*1e-3);
      } // for(int idr = 20; idr <=50; idr++)

    }
    TLorentzVector met = evt.met();
    bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L); 

    float chi2Value = 1000000; // log10(1000000) = 6
    float mwh = -1;
    float mtl = -1;
    float mth = -1;

    if (status){
      chi2Value = chi2ming1;
      mwh = m_chi2.getResult_Mwh();
      mtl = m_chi2.getResult_Mtl();
      mth = m_chi2.getResult_Mth();
      mtt = m_chi2.getResult_Mtt();
    }

    for (size_t z = 0; z < vjets.size(); ++z) {
      delete vjets[z];
    }
    vjets.clear();
    vjets_btagged.clear();

    //Fill histograms
    h->h1D("mtt", "", suffix)->Fill(mtt*1e-3, weight);
    h->h1D("mtlep_res", "", suffix)->Fill(mtl*1e-3, weight);
    h->h1D("mthad_res", "", suffix)->Fill(mth*1e-3, weight);
    h->h1D("mwhad_res", "", suffix)->Fill(mwh*1e-3, weight);
    h->h1D("chi2", "", suffix)->Fill(log10(chi2Value), weight);

    _tree_mtt = mtt*1e-3;
    _tree_weight = weight;
    _tree_syst = s;
    _tree_cat = -1;
    if (is2016run && evt.passes("rejets_2016") && !m_boosted && m_electron) _tree_cat = 2;
    if (is2016run && (evt.passes("rmujets_2016") || evt.passes("rmujetsmet110")) && !m_boosted && !m_electron) _tree_cat = 3;
    h->m_tree->Fill();
  }
}

