/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresTruthResolved.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"


AnaTtresTruthResolved::AnaTtresTruthResolved(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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
    
    m_hSvc.create1D("Chi2All", "; log_{10}#chi^{2} ; Events", 60, -1, 5);  

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

    m_hSvc.create1D("small_jet_pt", "; small jet p_{T} [GeV]; Events", 50, 0, 500);
    m_hSvc.create1D("small_jet_eta", "; small jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("small_jet_phi", "; small jet #phi ; Events", 32, -3.2, 3.2);

    m_hSvc.create1DVar("largeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN6, varBin6);
    m_hSvc.create1D("largeJetEta", "; large jet #eta ; Events", 20, -2., 2.);
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

    m_hSvc.create1D("cutFlow", ";Number of events; step", 20,0,20);
    m_hSvc.create1D("nMuon", ";Number of events; # #mu", 5,0,5);
    m_hSvc.create1D("muonPt", ";Number of events; p_{T,#mu}", 50,0,500);
    m_hSvc.create1D("electronPt", ";Number of events; p_{T,e}", 50,0,500);
    m_hSvc.create1D("nElectron", ";Number of events; # e", 5,0,5);
    m_hSvc.create1D("muon_ptvarcone30_over_pt", ";#mu ptvarcone30 / pt", 100,0,0.1);
  }

AnaTtresTruthResolved::~AnaTtresTruthResolved() {
}

void AnaTtresTruthResolved::run(const Event &evt, double weight, const std::string &s, int is2016run) {
  std::cout << "AnaTtresTruthResolved::run" << std::endl;

  if(!(evt.passes("no_cut"))) return;

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

  HistogramService *h = &m_hSvc;

  bool trig1(0); 
  bool trig2(0); 
  bool trig3(0);
  bool trig4(0);
  bool trig5(0);
  bool trig6(0);

  bool isTight = false;

  std::string suffix = s;
  h->h1D("cutFlow", "", suffix)->Fill(0);

  //Trigger selection
  bool pass_mu26_ivarmedium = evt.trigger("HLT_mu26_ivarmedium"); 
  bool pass_mu50 = evt.trigger("HLT_mu50"); 
  if(!pass_mu26_ivarmedium && !pass_mu50) return;
  h->h1D("cutFlow", "", suffix)->Fill(1);
  std::cout << "Pass Trigger selection" << std::endl;

  //Number of muon
  int nMuon_25=0, nMuon_30=0, muon_id=-1;
  for(int im=0;im<evt.muon().size();im++){
    std::cout << "muon " << im << std::endl;
    TLorentzVector mu = evt.muon()[im].mom();
    bool isTight = evt.muon()[im].isTight();
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float dphi = mu.DeltaPhi(evt.jet()[jet_idx].mom());
      float dy = mu.Rapidity() - evt.jet()[jet_idx].mom().Rapidity();
      float deltaR_tmp = std::sqrt(dy*dy + dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    h->h2D("lepPt_vs_drMuonJet", "", suffix)->Fill(mu.Perp()*1e-3,closejl_deltaR, weight);
    h->h1D("muonPt", "", suffix)->Fill(mu.Perp()*1e-3, weight);
    float iso  = evt.muon()[im].ptvarcone30() / evt.muon()[im].mom().Pt();
    h->h1D("muon_ptvarcone30_over_pt", "", suffix)->Fill(iso, weight);
    std::cout << "pt/isTight=" << "/" << mu.Perp()*1e-3 << "/" << isTight << std::endl;
    std::cout << "ptvarcone30/iso=" << evt.muon()[im].ptvarcone30() << "/" << iso << std::endl;
    if( closejl_deltaR < 0.04 + 10./(mu.Perp()*1e-3) ) continue;
    std::cout << "pass OR" << std::endl;
    //if(!isTight && iso>0.06) continue;
    if(iso>0.06) continue;
    if(mu.Perp()*1e-3 > 30){
      if(nMuon_30==0) muon_id=im;
      nMuon_30++;
    }
    if(mu.Perp()*1e-3>25) nMuon_25++;
  }

  //Number of muon > 0
  h->h1D("nMuon", "", suffix)->Fill(nMuon_30, weight);
  if(nMuon_30==0) return;
  h->h1D("cutFlow", "", suffix)->Fill(2);
  std::cout << "Pass #Muon >= 1" << std::endl;

  //Number of muon == 1
  if(nMuon_25>1) return;
  h->h1D("cutFlow", "", suffix)->Fill(3);
  std::cout << "Pass #Muon == 1" << std::endl;

  //Number of electron == 0
  int nElectron = 0;
  for(int ie=0;ie<evt.electron().size();ie++){
    std::cout << "electron " << ie << std::endl;
    TLorentzVector el = evt.electron()[ie].mom();
    if(ie==0) h->h1D("electronPt", "", suffix)->Fill(el.Perp()*1e-3, weight);
    if(el.Perp()*1e-3 > 25) nElectron++;
    bool isTight = evt.electron()[ie].isTightPP();
    std::cout << "pt/isTight=" << "/" << el.Perp()*1e-3 << "/" << isTight << std::endl;
  }
  h->h1D("nElectron", "", suffix)->Fill(nElectron, weight);
  if(nElectron!=0) return;
  h->h1D("cutFlow", "", suffix)->Fill(4);
  std::cout << "Pass #Electron == 0" << std::endl;

  //Trigger matching
  trig1 = evt.muon()[muon_id].HLT_mu26_ivarmedium();
  trig2 = evt.muon()[muon_id].HLT_mu50();
  if(!trig1 && !trig2) return;
  h->h1D("cutFlow", "", suffix)->Fill(5);
  std::cout << "Pass Trigger Matching" << std::endl;

  //Jet Cleaning
  std::cout << "jet clean is " << evt.passes("jet_clean") << std::endl;
  if(!(evt.passes("jet_clean"))) return;
  h->h1D("cutFlow", "", suffix)->Fill(6);
  std::cout << "Pass Jet Clean" << std::endl;

  //MET & MET+MWT 
  std::cout << "muon_id=" << muon_id << std::endl;
  TLorentzVector l = evt.muon()[muon_id].mom();
  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3;
  float MET = evt.met().Perp()*1e-3;
  float MET_plus_mWt = mWt + MET;
  std::cout << "MET/MWT=" << MET << "/" << mWt << std::endl;
  h->h1D("MET", "", suffix)->Fill(MET, weight);
  h->h1D("MET_plus_MWT", "", suffix)->Fill(MET_plus_mWt, weight);
  h->h1D("MWT", "", suffix)->Fill(mWt, weight);
  if(MET<20) return;
  h->h1D("cutFlow", "", suffix)->Fill(7);
  std::cout << "Pass MET" << std::endl;
  if(MET_plus_mWt<60) return;
  h->h1D("cutFlow", "", suffix)->Fill(8);
  std::cout << "Pass MET+mWt" << std::endl;

  //Number of small jets
  std::cout << "Number of small jet is " << evt.jet().size() << std::endl;
  int nJet=0;
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    float jet_pt = evt.jet()[jidx].mom().Perp()*1e-3;
    float jet_eta = evt.jet()[jidx].mom().Eta();
    float jet_phi = evt.jet()[jidx].mom().Phi();
    if(jet_pt>25) nJet++;
    h->h1D("small_jet_pt", "", suffix)->Fill(jet_pt, weight);
    h->h1D("small_jet_eta", "", suffix)->Fill(jet_eta, weight);
    h->h1D("small_jet_phi", "", suffix)->Fill(jet_phi, weight);
  }
  h->h1D("nJets", "", suffix)->Fill(nJet, weight);
  if(nJet==0) return;
  h->h1D("cutFlow", "", suffix)->Fill(9);
  
  if(nJet==1) return;
  h->h1D("cutFlow", "", suffix)->Fill(10);
  
  if(nJet==2) return;
  h->h1D("cutFlow", "", suffix)->Fill(11);

  if(nJet==3) return;
  h->h1D("cutFlow", "", suffix)->Fill(12);
  
  std::cout << "Pass #Small Jet" << std::endl;

  //Track B-tag
  int nTrackBtagJets=0;
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    float pt = evt.tjet()[jidx].mom().Pt()*1e-3;
    float eta = evt.tjet()[jidx].mom().Eta();
    float phi = evt.tjet()[jidx].mom().Phi();
    std::cout << "tjet pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl; 
    std::cout << "b-tag value = " << evt.tjet()[jidx].mv2c10() << std::endl;
    bool is_trackBtag = evt.tjet()[jidx].btag_mv2c10_70_trk();
    if (is_trackBtag && evt.tjet()[jidx].mom().Pt() > 10e3 && 
        fabs(evt.tjet()[jidx].mom().Eta()) < 2.5 && evt.tjet()[jidx].numConstituents() >= 2) nTrackBtagJets++; 
  }
  std::cout << "Number of track b-tagged jets is " << nTrackBtagJets << std::endl;
  h->h1D("nTrkBtagJets", "", suffix)->Fill(nTrackBtagJets, weight);  

  //Calo B-tag
  int nCaloBtagJets=0;
  std::vector<TLorentzVector*> vjets;
  std::vector<bool> vjets_btagged;
  vjets.clear();vjets_btagged.clear();
  for (int jid=0; jid < evt.jet().size(); ++jid){    
    bool is_caloBtag = evt.jet()[jid].btag_mv2c20_70();
    if(is_caloBtag) nCaloBtagJets++;
    vjets.push_back(new TLorentzVector(0,0,0,0));
    vjets[jid]->SetPtEtaPhiE(evt.jet()[jid].mom().Perp(), evt.jet()[jid].mom().Eta(), evt.jet()[jid].mom().Phi(), evt.jet()[jid].mom().E());
    std::cout << "calc chi2 jet pt/eta/phi=" << evt.jet()[jid].mom().Perp() << "/" <<  evt.jet()[jid].mom().Eta() << "/" << evt.jet()[jid].mom().Phi() << std::endl;
    bool is_btagged(false);
    for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx)
    {
      TLorentzVector tmpTJet;
      tmpTJet.SetPtEtaPhiE(evt.tjet()[bidx].mom().Perp(),evt.tjet()[bidx].mom().Eta(),evt.tjet()[bidx].mom().Phi(),evt.tjet()[bidx].mom().E());
      if(tmpTJet.DeltaR(*vjets[jid]) <= 0.4) {
        std::cout << "matched tjet pt/eta/phi=" << evt.tjet()[bidx].mom().Perp() << "/" << evt.tjet()[bidx].mom().Eta() << "/" << evt.tjet()[bidx].mom().Phi() << std::endl;
        std::cout << "btag score=" << evt.tjet()[bidx].mv2c10() << std::endl;
        if(evt.tjet()[bidx].mom().Pt() > 10e3 && fabs(evt.tjet()[bidx].mom().Eta())<2.5 &&
            evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].numConstituents()>=2)  is_btagged = true ;
        //if(evt.jet()[jid].mv2c10() > 0.6455) is_btagged = true;
        break; 
      } 
    } 
    std::cout << "isBtag=" << is_btagged << std::endl;
    vjets_btagged.push_back(is_btagged);
  }
  std::cout << "Number of calo b-tagged jets is " << nCaloBtagJets << std::endl;
  h->h1D("nCaloBtagJets", "", suffix)->Fill(nCaloBtagJets, weight);  

  if(nTrackBtagJets==0) return;
  h->h1D("cutFlow", "", suffix)->Fill(13);
  std::cout << "Pass B-tag" << std::endl;

  //TLorentzVector met = evt.met();
  TLorentzVector met;
  met.SetPtEtaPhiM(evt.met().Perp(), 0, evt.met().Phi(), 0);
  int igj3, igj4, igb3, igb4, ign1;
  double chi2ming1, chi2ming1H, chi2ming1L;
  bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L); 
  if(status){
    std::cout << "chi2ming1/chi2ming1H/chi2ming1L=" << chi2ming1 << "/" << chi2ming1H << "/" << chi2ming1L << std::endl;
  }

  std::cout << "chi2_all=" << evt.chi2_all() << std::endl;
  //float log10chi2 = (evt.chi2_all()>0)? log10(evt.chi2_all()) : 100;
  float log10chi2 = (chi2ming1>0)? log10(chi2ming1) : 100;
  std::cout << "log10chi2=" << log10chi2 << std::endl;
  h->h1D("Chi2All", "", suffix)->Fill(log10chi2, weight);  
  if(log10chi2>0.9) return;
  h->h1D("cutFlow", "", suffix)->Fill(14);
  std::cout << "Pass Chi2" << std::endl;

  //truth jet matching
  int bt_sj_idx=-1, btbar_sj_idx=-1;
  float min_dr_sj_bt = 999., min_dr_sj_btbar = 999.;
  std::cout << "Nubmer of small jet is " << evt.jet().size() << std::endl;
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
    std::cout << "btag value=" << evt.jet()[bidx].mv2c20() << std::endl;
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
  std::cout << "Number of small tjet is " << evt.tjet().size() << std::endl;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
    float pt = evt.tjet()[bidx].mom().Pt();
    float eta = evt.tjet()[bidx].mom().Eta();
    float phi = evt.tjet()[bidx].mom().Phi();
    float drt = TuDoAtlas::calc_delta_r(eta,b_from_t_eta,phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(eta,b_from_tbar_eta,phi,b_from_tbar_phi);
    int isBtag = evt.tjet()[bidx].btag_mv2c10_70_trk();
    std::cout << "track jet pt/eta/phi/btag=" << pt << "/" << eta << "/" << phi << "/" << isBtag << std::endl;
    std::cout << "btag value=" << evt.tjet()[bidx].mv2c10() << std::endl;
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
  std::cout << "btbar matched small jet id=" << btbar_sj_idx << std::endl;
  if(bt_sj_idx!=-1 && bt_tj_idx!=-1){
    h->h1D("trackBtagBDT", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mv2c10());
    h->h1D("trackBtagBDT", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mv2c10());
    h->h1D("caloBtagBDT", "", suffix)->Fill(evt.jet()[bt_sj_idx].mv2c20());
    h->h1D("caloBtagBDT", "", suffix)->Fill(evt.jet()[bt_sj_idx].mv2c20());
    std::cout << "calo btag value bt/btbar=" << evt.jet()[bt_sj_idx].mv2c20() << "/" << evt.jet()[btbar_sj_idx].mv2c20() << std::endl;
    std::cout << "track btag value bt/btbar=" << evt.tjet()[bt_tj_idx].mv2c10() << "/" << evt.tjet()[btbar_tj_idx].mv2c10() << std::endl;
    if(isLeptonicT){//leptonic top
      h->h1D("ltop_pt_jet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_ltop_pt_jet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_ltop_matched_jet_pt", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[bt_sj_idx].btag_mv2c20_70()){
        std::cout << "ltop jet  btagged!!" << std::endl;
        h->h1D("ltop_pt_jet_passCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
      else std::cout << "ltop jet not btagged!!" << std::endl;
    }
    else{//hadronic top
      h->h1D("htop_pt_jet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_jet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_htop_matched_jet_pt", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[bt_sj_idx].btag_mv2c20_70()){
        std::cout << "htop jet btagged!!" << std::endl;
        h->h1D("htop_pt_jet_passCaloBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[bt_sj_idx].mom().Perp()*1e-3);
      }
      else std::cout << "htop jet not btagged!!" << std::endl;
    }
  }
  if(btbar_sj_idx!=-1 && btbar_tj_idx!=-1){
    if(isLeptonicTbar){//leptonic antitop
      h->h1D("ltop_pt_jet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_ltop_pt_jet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_ltop_matched_jet_pt", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[btbar_sj_idx].btag_mv2c20_70()){
        std::cout << "ltop jet btagged!!" << std::endl;
        h->h1D("ltop_pt_jet_passCaloBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_ltop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      }
      else std::cout << "ltop jet not btagged!!" << std::endl;
    }
    else{//hadronic antitop
      h->h1D("htop_pt_jet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_htop_pt_jet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_htop_matched_jet_pt", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      if(evt.jet()[btbar_sj_idx].btag_mv2c20_70()){
        std::cout << "htop jet btagged!!" << std::endl;
        h->h1D("htop_pt_jet_passCaloBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_jet_passCaloBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_htop_matched_jet_pt_passCaloBtag", "", suffix)->Fill(evt.jet()[btbar_sj_idx].mom().Perp()*1e-3);
      }
      else std::cout << "htop jet not btagged!!" << std::endl;
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
        std::cout << "ltop tjet btagged!!" << std::endl;
        h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
      else std::cout << "ltop tjet not btagged!!" << std::endl;
    }
    else{//hadronic top
      h->h1D("htop_pt_tjet", "", suffix)->Fill(t_pt);
      h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_t_pt);
      h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      if(evt.tjet()[bt_tj_idx].btag_mv2c10_70_trk()){
        std::cout << "htop tjet btagged!!" << std::endl;
        h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(t_pt);
        h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_t_pt);
        h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[bt_tj_idx].mom().Perp()*1e-3);
      }
      else std::cout << "htop tjet not btagged!!" << std::endl;
    }
  }
  if(btbar_tj_idx!=-1 && btbar_sj_idx!=-1){
    if(isLeptonicTbar){//leptonic top
      h->h1D("ltop_pt_tjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_ltop_pt_tjet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_ltop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      if(evt.tjet()[btbar_tj_idx].btag_mv2c10_70_trk()){
        std::cout << "ltop tjet btagged!!" << std::endl;
        h->h1D("ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_ltop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_ltop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
      else std::cout << "ltop tjet not btagged!!" << std::endl;
    }
    else{//hadronic top
      h->h1D("htop_pt_tjet", "", suffix)->Fill(tbar_pt);
      h->h1D("b_from_htop_pt_tjet", "", suffix)->Fill(b_from_tbar_pt);
      h->h1D("b_from_htop_matched_tjet_pt", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      if(evt.tjet()[btbar_tj_idx].btag_mv2c10_70_trk()){
        std::cout << "htop tjet btagged!!" << std::endl;
        h->h1D("htop_pt_tjet_passTrackBtag", "", suffix)->Fill(tbar_pt);
        h->h1D("b_from_htop_pt_tjet_passTrackBtag", "", suffix)->Fill(b_from_tbar_pt);
        h->h1D("b_from_htop_matched_tjet_pt_passTrackBtag", "", suffix)->Fill(evt.tjet()[btbar_tj_idx].mom().Perp()*1e-3);
      }
      else std::cout << "htop tjet not btagged!!" << std::endl;
    }
  }

  //top tagging
  std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
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
    std::cout << "top tag def/new=" << good_def << "/" << evt.largeJet()[ij].BDT_TOPtag() << std::endl;
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

}
