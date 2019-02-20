/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtres21ResoEl.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"


AnaTtres21ResoEl::AnaTtres21ResoEl(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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

    m_hSvc.create1D("nMuon", ";Number of events; # #mu", 5,0,5);
    m_hSvc.create1D("muonPt", ";Number of events; p_{T,#mu}", 50,0,500);
    m_hSvc.create1D("electronPt", ";Number of events; p_{T,e}", 50,0,500);
    m_hSvc.create1D("nElectron", ";Number of events; # e", 5,0,5);
    m_hSvc.create1D("muon_ptvarcone30_over_pt", ";#mu ptvarcone30 / pt", 100,0,0.1);

    std::string channel[6] = {"All","FullHad","FullLep","SemiLepTau","SemiLepMu","SemiLepEl"};
    for(int ch=0;ch<6;ch++){
      m_hSvc.create1D(Form("cutFlow_%s",channel[ch].c_str()), ";Number of events; step", 20,0,20);
    }
  }

AnaTtres21ResoEl::~AnaTtres21ResoEl() {
}

void AnaTtres21ResoEl::run(const Event &evt, double weight, const std::string &s, int is2016run) {
  std::cout << "AnaTtres21ResoEl::run" << std::endl;

  if(!(evt.passes("no_cut"))) return;

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

  bool isFullHad = (!isLeptonicT && !isLeptonicTbar);
  bool isFullLep = (isLeptonicT && isLeptonicTbar);
  bool isSemiLep = (!isFullHad && !isFullLep)? true : false;
  bool containTau = (Wdecay2_from_t_pdgId==-15 || Wdecay1_from_tbar_pdgId==15)? true : false;
  bool containMu = (Wdecay2_from_t_pdgId==-13 || Wdecay1_from_tbar_pdgId==13)? true : false;
  bool containEl = (Wdecay2_from_t_pdgId==-11 || Wdecay1_from_tbar_pdgId==11)? true : false;

  std::string channel;
  if(isFullHad) channel="FullHad";
  else if(isFullLep) channel="FullLep";
  else if(containTau) channel="SemiLepTau";
  else if(containMu) channel="SemiLepMu";
  else if(containEl) channel="SemiLepEl";
  //std::cout << "channel is " << channel << std::endl;

  HistogramService *h = &m_hSvc;

  bool isTight = false;

  std::string suffix = s;
  h->h1D("cutFlow_All", "", suffix)->Fill(0);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(0);

  //Trigger selection
  bool passtrig_el1 = evt.trigger("HLT_e26_lhtight_nod0_ivarloose"); 
  bool passtrig_el2 = evt.trigger("HLT_e60_lhmedium_nod0"); 
  bool passtrig_el3 = evt.trigger("HLT_e140_lhloose_nod0"); 
  if(!passtrig_el1 && !passtrig_el2 && !passtrig_el3) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(1);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(1);
  std::cout << "Pass Trigger selection" << std::endl;

  //no selection muon
  int numNocalibMuon=0,numNocalibMuon25=0,muon_id_nocalib=-1;
  for(int im=0;im<evt.nocalibmu().size();im++){
    TLorentzVector mu = evt.nocalibmu()[im].mom();
    float pt = mu.Perp();
    float eta = mu.Eta();
    float phi = mu.Phi();
    float ptvarcone30 = evt.nocalibmu()[im].ptvarcone30();
    float iso = ptvarcone30 / pt;
    int quality = evt.nocalibmu()[im].quality();
    int accept = evt.nocalibmu()[im].accept();
    float d0sig = evt.nocalibmu()[im].sd0();
    float z0sintheta = evt.nocalibmu()[im].Dz0();
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float deta = eta - evt.jet()[jet_idx].mom().Eta(); 
      float dphi = acos(cos(phi - evt.jet()[jet_idx].mom().Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    if(fabs(eta)>2.5) continue;
    if(iso>0.06) continue;
    if(quality>1) continue;
    if(!accept) continue;
    if(d0sig>3) continue;
    if(z0sintheta>0.5) continue;
    if(closejl_deltaR<0.04+10./(pt*1e-3)) continue;
    if(pt*1e-3>25) numNocalibMuon25++;
    if(pt*1e-3<30) continue;
    if(numNocalibMuon==0) muon_id_nocalib=im;
    numNocalibMuon++;
  }

  //Number of electron == 0
  int nElectron = 0,nElectron25=0,el_id=-1;
  for(int ie=0;ie<evt.electron().size();ie++){
    std::cout << "electron " << ie << std::endl;
    TLorentzVector el = evt.electron()[ie].mom();
    if(ie==0) h->h1D("electronPt", "", suffix)->Fill(el.Perp()*1e-3, weight);
    if(el.Perp()*1e-3 > 25) nElectron25++;
    if(el.Perp()*1e-3 > 30) {
      if(nElectron==0) el_id=ie;
      nElectron++;
    }
    bool isTight = evt.electron()[ie].isTightPP();
    std::cout << "pt/isTight=" << "/" << el.Perp()*1e-3 << "/" << isTight << std::endl;
  }
  h->h1D("nElectron", "", suffix)->Fill(nElectron, weight);
  if(nElectron==0) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(2);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(2);
  std::cout << "Pass #Electron > 0" << std::endl;

  if(nElectron>1) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(3);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(3);
  std::cout << "Pass #Electron == 1" << std::endl;

  if(numNocalibMuon25>0) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(4);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(4);
  std::cout << "Pass #Muon == 0" << std::endl;

  //Trigger matching
  bool trig1 = evt.electron()[el_id].HLT_e26_lhtight_nod0_ivarloose();
  bool trig2 = evt.electron()[el_id].HLT_e60_lhmedium_nod0();
  bool trig3 = evt.electron()[el_id].HLT_e140_lhloose_nod0();
  std::cout << "matching trig1/trig2/trig3=" << trig1 << "/" << trig2 << "/" << trig3 << std::endl;
  if(!trig1 && !trig2 && !trig3) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(5);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(5);
  std::cout << "Pass Trigger Matching" << std::endl;

  //Jet Cleaning
  //std::cout << "jet clean is " << evt.passes("jet_clean") << std::endl;
  //if(!(evt.passes("jet_clean"))) return;
  //h->h1D("cutFlow", "", suffix)->Fill(6);
  //std::cout << "Pass Jet Clean" << std::endl;

  //MET & MET+MWT 
  TLorentzVector l = evt.electron()[el_id].mom();
  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3;
  float MET = evt.met().Perp()*1e-3;
  float MET_plus_mWt = mWt + MET;
  std::cout << "MET/MWT=" << MET << "/" << mWt << std::endl;
  h->h1D("MET", "", suffix)->Fill(MET, weight);
  h->h1D("MET_plus_MWT", "", suffix)->Fill(MET_plus_mWt, weight);
  h->h1D("MWT", "", suffix)->Fill(mWt, weight);
  if(MET<20) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(6);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(6);
  std::cout << "Pass MET" << std::endl;
  if(MET_plus_mWt<60) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(7);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(7);
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
  h->h1D("cutFlow_All", "", suffix)->Fill(8);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(8);
  
  if(nJet==1) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(9);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(9);
  
  if(nJet==2) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(10);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(10);

  if(nJet==3) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(11);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(11);
  
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

  //No selection track jet
  int numTrackJets=0,numBtaggedTrackJets=0;
  for (size_t jidx = 0; jidx < evt.nocalibtjet().size(); ++jidx){
    TLorentzVector tj = evt.nocalibtjet()[jidx].mom();
    if(tj.Perp()*1e-3<10) continue;
    if(fabs(tj.Eta())>2.4) continue;
    if(evt.nocalibtjet()[jidx].numConstituents()<2) continue;
    float closeje_deltaR=999.;
    for (unsigned int ie=0; ie < evt.electron().size(); ++ie){    
      TLorentzVector el = evt.electron()[ie].mom();
      float deta = el.Eta() - tj.Eta(); 
      float dphi = acos(cos(el.Phi() - tj.Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closeje_deltaR) closeje_deltaR = deltaR_tmp;
    }//for     
    if(closeje_deltaR<0.2) continue;
    float closejmu_deltaR=999.,closejmu_deltaR_def=999.;
    for(int im=0;im<evt.nocalibmu().size();im++){
      TLorentzVector mu = evt.nocalibmu()[im].mom();
      float closejl_deltaR=999.;
      for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
        float deta = mu.Eta() - evt.jet()[jet_idx].mom().Eta(); 
        float dphi = acos(cos(mu.Phi() - evt.jet()[jet_idx].mom().Phi())); 
        float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
        if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
      }//for     
      if(mu.Perp()*1e-3<30) continue;
      if(fabs(mu.Eta())>2.5) continue;
      if(evt.nocalibmu()[im].ptvarcone30()/mu.Perp()>0.06) continue;
      if(evt.nocalibmu()[im].quality()>1) continue;
      if(!evt.nocalibmu()[im].accept()) continue;
      if(evt.nocalibmu()[im].sd0()>3) continue;
      if(evt.nocalibmu()[im].Dz0()>0.5) continue;
      float deta = mu.Eta() - tj.Eta(); 
      float dphi = acos(cos(mu.Phi() - tj.Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if(closejl_deltaR>0.04+10./(mu.Perp()*1e-3)) {//OR for ttbar resonance
        if (deltaR_tmp < closejmu_deltaR) closejmu_deltaR = deltaR_tmp;
      }
      if(closejl_deltaR>0.4) {//default OR 
        if (deltaR_tmp < closejmu_deltaR_def) closejmu_deltaR_def = deltaR_tmp;
      }
    }//for     
    numTrackJets++;
    if(evt.nocalibtjet()[jidx].btag_mv2c10_70_trk()) numBtaggedTrackJets++;
  }

  //if(nTrackBtagJets==0) return;
  if(numBtaggedTrackJets==0) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(12);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(12);
  std::cout << "Pass B-tag" << std::endl;

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
  h->h1D("cutFlow_All", "", suffix)->Fill(13);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(13);
  std::cout << "Pass Chi2" << std::endl;

}
