/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresFullHad.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"
#include "TopNtupleAnalysis/KinematicUtils.h"


AnaTtresFullHad::AnaTtresFullHad(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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
    m_hSvc.create1D("nGoodLargeJets", "; number of good large-R jets ; Events", 10, -0.5, 9.5);
    m_hSvc.create1D("nGoodLargeJets50", "; number of good large-R jets ; Events", 10, -0.5, 9.5);

    m_hSvc.create1D("nTrkBtagJets", "; number of b-tagged track jets ; Events", 10, 0, 10);  
    m_hSvc.create1D("nCaloBtagJets", "; number of b-tagged calo jets ; Events", 10, 0, 10);  

    m_hSvc.create1DVar("MET", "; missing E_{T} [GeV]; Events", varBinN1, varBin1);
    m_hSvc.create1DVar("MET_plus_MWT", "; missing E_{T} + m^{W}_{T} [GeV]; Events", varBinN1, varBin1);
    m_hSvc.create1D("MWT", "; W transverse mass [GeV]; Events", 20, 0, 200);

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

    std::string channel[5] = {"FullHad","FullLep","SemiLepTau","SemiLepMu","SemiLepEl"};
    std::string type[3] = {"FromW","FromB","FromOther"};
    m_hSvc.create1D("trueMuonPtToB", "; p_{T,#mu} to b [GeV];", 100,0,100);
    m_hSvc.create1D("trueMuonPtToW", "; p_{T,#mu} to W [GeV];", 100,0,100);
    for(int ch=0;ch<5;ch++) {
      m_hSvc.create1D(Form("cutFlow%s",channel[ch].c_str()), ";Number of events; step", 20,0,20);
      m_hSvc.create1D(Form("nMuon%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1D(Form("nMuonFinal%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1D(Form("nMuonNoIso%s",channel[ch].c_str()), ";# #mu;Number of events", 5,0,5);
      m_hSvc.create1DVar(Form("mtt%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
      for(int t=0;t<3;t++){
        m_hSvc.create1D(Form("initialMuon%s_pt_%s",type[t].c_str(),channel[ch].c_str()), "; p_{T,#mu};", 50,0,500);
        m_hSvc.create1D(Form("initialMuonPtToJet%s_%s",type[t].c_str(),channel[ch].c_str()), "; p_{T,#mu} to jet [GeV];", 100,0,200);
        m_hSvc.create1D(Form("initialMuonPtToW%s_%s",type[t].c_str(),channel[ch].c_str()), "; p_{T,#mu} to W [GeV];", 100,0,200);
        m_hSvc.create2D(Form("initialMuonPtToJet_vs_ToW%s_%s",type[t].c_str(),channel[ch].c_str()), "; p_{T,#mu} to jet [GeV];p_{T,#mu} to W [GeV]", 100,0,200,100,0,200);
        m_hSvc.create1D(Form("initialMuon%s_quality_%s",type[t].c_str(),channel[ch].c_str()), "; quality;", 5,0,5);
        m_hSvc.create1D(Form("leadingInitialMuon%s_pt_%s",type[t].c_str(),channel[ch].c_str()), ";p_{T,#mu};", 50,0,500);
        m_hSvc.create1D(Form("initialMuon_sd0_%s_%s",type[t].c_str(),channel[ch].c_str()), ";sd0;Number of events ", 200,0,20);
        m_hSvc.create1D(Form("initialMuon_dz0_%s_%s",type[t].c_str(),channel[ch].c_str()), ";dz0;Number of events ", 100,0,2);
        m_hSvc.create1D(Form("leadingInitialMuon_sd0_%s_%s",type[t].c_str(),channel[ch].c_str()), ";sd0;Number of events; ", 200,0,20);
        m_hSvc.create1D(Form("leadingInitialMuon_dz0_%s_%s",type[t].c_str(),channel[ch].c_str()), ";dz0;Number of events; ", 100,0,2);
        m_hSvc.create1D(Form("leadingInitialMuon%s_quality_%s",type[t].c_str(),channel[ch].c_str()), "; quality;", 5,0,5);
        m_hSvc.create1D(Form("initialMuonPtvarcone30OverPt_%s_%s",type[t].c_str(),channel[ch].c_str()), ";#mu ptvarcone30 / pt", 100,0,1);
        m_hSvc.create2D(Form("lepPt_vs_drMuonJet_%s_%s",type[t].c_str(),channel[ch].c_str()),    "; Pt of lept [GeV]; #DeltaR_{#mu,jet}", 100, 0, 500,100,0,2);
        m_hSvc.create2D(Form("initialMuon_sd0_vs_dz0_%s_%s",type[t].c_str(),channel[ch].c_str()), ";sd0;dz0 ", 100,0,10,100,0,0.2);
        m_hSvc.create1DVar(Form("mtt%s_%s",type[t].c_str(),channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
      }
    }
    m_hSvc.create1D("cutFlowAll", ";Number of events; step", 20,0,20);
    m_hSvc.create1D("nMuonAll", ";# #mu;Number of events", 5,0,5);
    m_hSvc.create1D("nMuonFinalAll", ";# #mu;Number of events", 5,0,5);
    m_hSvc.create1D("nMuonNoIsoAll", "# #mu;Number of events", 5,0,5);
    m_hSvc.create1D("muonPt", ";Number of events; p_{T,#mu}", 50,0,500);
    m_hSvc.create1D("truthMuonPt", ";Number of events; p_{T,#mu}", 50,0,500);
    m_hSvc.create1D("drTruthMuonInitialMuon", ";#DeltaR;Number of events", 100,0,6);
    m_hSvc.create1D("drTruthMuonInitialMuon_final", ";#DeltaR;Number of events", 100,0,6);
    m_hSvc.create1D("initialMuonPt", ";Number of events; p_{T,#mu}", 50,0,500);
    m_hSvc.create1D("leadingInitialMuonPt", ";Number of events; p_{T,#mu}", 50,0,500);
    m_hSvc.create1D("electronPt", ";Number of events; p_{T,e}", 50,0,500);
    m_hSvc.create1D("nElectron", ";Number of events; # e", 5,0,5);
    m_hSvc.create1D("muon_ptvarcone30_over_pt", ";#mu ptvarcone30 / pt", 100,0,1);
    m_hSvc.create1D("initial_muon_ptvarcone30_over_pt", ";#mu ptvarcone30 / pt", 100,0,1);
    m_hSvc.create2D("muonPt_vs_ptvarcone30_over_pt",    "; #mu pt [GeV]; #mu ptvarcone30 / pt", 50, 0, 500,100,0,1);
    m_hSvc.create1DVar("mttAll", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
    m_hSvc.create1D("muon_sd0", "sd0;Number of events; ", 100,0,10);
    m_hSvc.create1D("muon_dz0", "dz0;Number of events; ", 100,0,10);
    m_hSvc.create2D("muon_sd0_vs_dz0", "sd0;dz0; ", 100,0,10,100,0,10);
    m_hSvc.create2D("initialMuon_sd0_vs_dz0", "sd0;dz0; ", 100,0,10,100,0,10);
    m_hSvc.create1D("nInitialMuonFinal", ";# of mu;Number of events", 5,0,5);


  }

AnaTtresFullHad::~AnaTtresFullHad() {
}

void AnaTtresFullHad::run(const Event &evt, double weight, const std::string &s, int is2016run) {

  std::cout << "AnaTtresFullHad::run" << std::endl;
  HistogramService *h = &m_hSvc;
  std::string suffix = s;

  std::string type[3] = {"FromW","FromB","FromOther"};

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
  std::cout << "Wdecay1_from_t_pdgId=" << Wdecay1_from_t_pdgId << std::endl;
  std::cout << "Wdecay2_from_t_pdgId=" << Wdecay2_from_t_pdgId << std::endl;
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
  if(isLeptonicT){
    std::cout << "leptonic top" << std::endl;
    TLorentzVector W_from_t = evt.MC_Wdecay1_from_t() + evt.MC_Wdecay2_from_t();
    float dphi_muW = evt.MC_Wdecay2_from_t().DeltaPhi(W_from_t); 
    std::cout << "muon pt/eta/phi=" << evt.MC_Wdecay2_from_t().Perp() << "/" << evt.MC_Wdecay2_from_t().Eta() << "/" << evt.MC_Wdecay2_from_t().Phi() << std::endl;
    float mupt_from_W = evt.MC_Wdecay2_from_t().P() * sin(dphi_muW);
    std::cout << "mupt_from_W = " << mupt_from_W*1e-3 << std::endl;
    h->h1D("trueMuonPtToW", "", suffix)->Fill(fabs(mupt_from_W)*1e-3,weight);
  }
  else {
    std::cout << "hadronic top" << std::endl;
    float b_pt = b_from_t.Perp();
    float b_eta = b_from_t.Eta();
    float b_phi = b_from_t.Phi();
    std::cout << "b from t pt/eta/phi=" << b_pt << "/" << b_eta << "/" << b_phi << std::endl;
    float min_dr=999.;
    float min_mupt_from_b=0.;
    for(int t=0;t<evt.truth().size();t++){
      TLorentzVector truemu;
      float px = evt.truth()[t].px();
      float py = evt.truth()[t].py();
      float pz = evt.truth()[t].pz();
      float e = evt.truth()[t].e();
      truemu.SetPxPyPzE(px,py,pz,e);
      float pt = truemu.Perp();
      float eta = truemu.Eta();
      float phi = truemu.Phi();
      float deta = b_eta - eta;
      float dphi = acos(cos(b_phi - phi));
      float dr = sqrt(deta*deta+dphi*dphi);
      int type = evt.truth()[t].muon_type();
      float dphi_mub = truemu.DeltaPhi(b_from_t);
      float mupt_from_b = truemu.P() * sin(dphi_mub);
      std::cout << "muon pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl; 
      std::cout << "dr/type=" << dr << "/" << type << std::endl;
      if(dr<min_dr && type==1){
        min_dr=dr;
        min_mupt_from_b = mupt_from_b; 
      }
    } 
    if(min_dr<0.2){
      std::cout << "mupt from b = " << min_mupt_from_b << std::endl;
      h->h1D("trueMuonPtToB", "", suffix)->Fill(fabs(min_mupt_from_b),weight);
    }
  }

  bool isFullHad = (!isLeptonicT && !isLeptonicTbar);
  bool isFullLep = (isLeptonicT && isLeptonicTbar);
  bool isSemiLep = (!isFullHad && !isFullLep)? true : false;
  bool containTau = (Wdecay2_from_t_pdgId==-15 || Wdecay1_from_tbar_pdgId==15)? true : false;
  bool containMu = (Wdecay2_from_t_pdgId==-13 || Wdecay1_from_tbar_pdgId==13)? true : false;
  bool containEl = (Wdecay2_from_t_pdgId==-11 || Wdecay1_from_tbar_pdgId==11)? true : false;

  TLorentzVector semilepmu;
  std::string channel;
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
  std::cout << "channel is " << channel << std::endl;
  std::cout << "true muon from t pt/eta/phi=" << semilepmu.Perp()*1e-3 << "/" << semilepmu.Eta() << "/" << semilepmu.Phi() << std::endl;

  bool trig1(0); 
  bool trig2(0); 
  bool trig3(0);
  bool trig4(0);
  bool trig5(0);
  bool trig6(0);

  bool isTight = false;

  if(!(evt.passes("no_cut"))) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(0);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(0);

  //Trigger selection
  bool pass_mu26_ivarmedium = evt.trigger("HLT_mu26_ivarmedium"); 
  bool pass_mu50 = evt.trigger("HLT_mu50"); 
  bool pass_met = evt.trigger("HLT_xe120_mht"); 
  bool pass_ljet = evt.trigger("HLT_j420_a10r_L1J100"); 
  //std::cout << "pass_ljet=" << pass_ljet << std::endl;
  //std::cout << "pass_met=" << pass_met << std::endl;
  //if(!pass_mu26_ivarmedium && !pass_mu50) return;
  //if(pass_mu26_ivarmedium || pass_mu50) return;
  //if(!pass_met) return;
  if(!pass_ljet) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(1);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(1);
  //std::cout << "Pass Trigger selection" << std::endl;

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
  h->h1D("cutFlowAll", "", suffix)->Fill(2);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(2);
  std::cout << "Pass #Electron == 0" << std::endl;

  //Number of muon
  int nMuon = 0, muon_id=-1;
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
    float dz0 = evt.muon()[im].Dz0();
    float sd0 = evt.muon()[im].sd0();
    //std::cout << "isolated muon pt/eta/phi=" << mu.Perp()*1e-3 << "/" << mu.Eta() << "/" << mu.Phi() << std::endl;
    //std::cout << "dz0=" << dz0 << std::endl;
    //std::cout << "sd0=" << sd0 << std::endl;
    //std::cout << "iso=" << iso << std::endl;
    //std::cout << "isTight/ptvarcone30/iso=" << isTight << "/" << evt.muon()[im].ptvarcone30() << "/" << iso << std::endl;
    //if( closejl_deltaR < 0.04 + 10./(mu.Perp()*1e-3) ) continue;
    //std::cout << "pass OR" << std::endl;
    //if(!isTight && iso>0.06) continue;
    //if(iso>0.06) continue;
    if(mu.Perp()*1e-3 < 25) continue;
    h->h2D("lepPt_vs_drMuonJet", "", suffix)->Fill(mu.Perp()*1e-3,closejl_deltaR, weight);
    h->h1D("muonPt", "", suffix)->Fill(mu.Perp()*1e-3, weight);
    h->h1D("muon_ptvarcone30_over_pt", "", suffix)->Fill(iso, weight);
    h->h1D("muon_sd0", "", suffix)->Fill(fabs(sd0), weight);
    h->h1D("muon_dz0", "", suffix)->Fill(fabs(dz0), weight);
    h->h2D("muon_sd0_vs_dz0", "", suffix)->Fill(fabs(sd0),fabs(dz0), weight);
    if(nMuon==0) muon_id=im;
    nMuon++;
  }

  //Number of muon > 0
  h->h1D("nMuonAll", "", suffix)->Fill(nMuon, weight);
  h->h1D(Form("nMuon%s",channel.c_str()), "", suffix)->Fill(nMuon, weight);
  if(nMuon!=0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(3);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(3);
  //std::cout << "Pass #Muon >= 1" << std::endl;

  //Initial muon
  int nIniMuon=0;
  int muid=-1;
  float maxpt=0.;
  for(int imu=0;imu<evt.inimu().size();imu++){
    float pt = evt.inimu()[imu].pt()*1e-3;
    float eta = evt.inimu()[imu].eta();
    float phi = evt.inimu()[imu].phi();
    float ptvarcone30 = evt.inimu()[imu].ptvarcone30();
    float iso = ptvarcone30 / evt.inimu()[imu].pt();
    int quality = evt.inimu()[imu].quality();
    int accept = evt.inimu()[imu].accept();
    float d0sig = evt.inimu()[imu].d0sig();
    float z0sintheta = evt.inimu()[imu].z0sintheta();
    //std::cout << "initial muon pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl;
    //std::cout << "iso=" << iso << std::endl;
    //std::cout << "quality=" << quality << std::endl;
    //std::cout << "d0sig/z0sintheta=" << d0sig << "/" << z0sintheta << std::endl;
    std::cout << "accept=" << accept << std::endl;
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float deta = evt.inimu()[imu].eta() - evt.jet()[jet_idx].mom().Eta(); 
      float dphi = acos(cos(evt.inimu()[imu].phi() - evt.jet()[jet_idx].mom().Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    //std::cout << "closejl_deltaR/thr=" << closejl_deltaR << "/" << 0.04 + 10./(evt.inimu()[imu].pt()*1e-3) << std::endl;
    if(closejl_deltaR>0.04+10./(evt.inimu()[imu].pt()*1e-3))
      h->h1D("initial_muon_ptvarcone30_over_pt", "", suffix)->Fill(iso, weight);
    if(pt<25) continue;
    if(quality>1) continue;
    if(accept==0) continue;
    h->h1D("initialMuonPt", "", suffix)->Fill(pt, weight);
    h->h2D("initialMuon_sd0_vs_dz0", "", suffix)->Fill(fabs(d0sig),fabs(z0sintheta), weight);
    if(fabs(d0sig)>2) continue;
    if(fabs(z0sintheta)>0.1) continue;
    if(pt>maxpt){
      maxpt=pt;
      muid=imu;
    }
    nIniMuon++;
  }
  //std::cout << "leading initial muon pt/iso=" << maxpt << "/" << evt.inimu()[muid].ptvarcone30()/evt.inimu()[muid].pt() << std::endl;
  std::cout << "leading initial muon pt/eta/phi=" << maxpt << "/" << evt.inimu()[muid].eta() << "/" << evt.inimu()[muid].phi() << std::endl;
  if(channel=="SemiLepMu"){
    float deta = evt.inimu()[muid].eta() - semilepmu.Eta(); 
    float dphi = acos(cos(evt.inimu()[muid].phi() - semilepmu.Phi())); 
    float dr = sqrt(deta*deta+dphi*dphi);
    //std::cout << "dr truth&initial=" << dr << std::endl;
    h->h1D("drTruthMuonInitialMuon", "", suffix)->Fill(dr, weight);
    h->h1D("leadingInitialMuonPt", "", suffix)->Fill(maxpt, weight);
    h->h1D("truthMuonPt", "", suffix)->Fill(semilepmu.Perp()*1e-3, weight);
    h->h2D("muonPt_vs_ptvarcone30_over_pt", "", suffix)->Fill(maxpt, evt.inimu()[muid].ptvarcone30()/evt.inimu()[muid].pt(), weight);
  }
  h->h1D(Form("nMuonNoIso%s",channel.c_str()), "", suffix)->Fill(nIniMuon, weight);
  if(nIniMuon>0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(4);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(4);
  for(int t=0;t<evt.truth().size();t++){
    if(abs(evt.truth()[t].id())==13) {
      float px = evt.truth()[t].px();
      float py = evt.truth()[t].py();
      float pt = sqrt(px*px+py*py);
      std::cout << "truth muon!!" << std::endl;
      std::cout << "pt=" << pt << std::endl;
    }
  }

  //std::cout << "size of truth particle is " << evt.truth().size() << std::endl;
  int leading_muon_type=-1;
  for(int m=0;m<evt.inimu().size();m++){
    if(evt.inimu()[m].pt()*1e-3<25) continue;
    if(evt.inimu()[m].quality()>1) continue;
    if(evt.inimu()[m].accept()==0) continue;
    if(fabs(evt.inimu()[m].d0sig())>2) continue;
    if(fabs(evt.inimu()[m].z0sintheta())>0.1) continue;
    float min_dr=999;
    int min_truth_id=-1;
    for(int t=0;t<evt.truth().size();t++){
      float px = evt.truth()[t].px();
      float py = evt.truth()[t].py();
      float pt = sqrt(px*px+py*py);
      float pz = evt.truth()[t].pz();
      float e = evt.truth()[t].e();
      int pdgid = evt.truth()[t].id();
      int status = evt.truth()[t].status();
      int muon_type = evt.truth()[t].muon_type();//0:fromW 1:fromB 2:fromOther
      if(muon_type==-1) continue;
      TLorentzVector particle;
      particle.SetPxPyPzE(px,py,pz,e);
      if(pt<1) continue;
      float eta = particle.Eta(); 
      float phi = particle.Phi(); 
      float deta = eta - evt.inimu()[m].eta();
      float dphi = acos(cos(phi - evt.inimu()[m].phi()));
      float dr = sqrt(deta*deta+dphi*dphi);
      float deta_mu = eta - semilepmu.Eta();
      float dphi_mu = acos(cos(phi - semilepmu.Phi()));
      float dpt_mu = fabs(pt) - fabs(semilepmu.Pt()*1e-3);
      if(dr<0.2 && dr<min_dr) {
        min_dr=dr;
        min_truth_id=t;
      }
    }
    if(min_truth_id!=-1){
      float px = evt.truth()[min_truth_id].px();
      float py = evt.truth()[min_truth_id].py();
      float pz = evt.truth()[min_truth_id].pz();
      float e = evt.truth()[min_truth_id].e();
      TLorentzVector ma_particle;
      ma_particle.SetPxPyPzE(px,py,pz,e);
      float pt = ma_particle.Pt();
      int pdgid = evt.truth()[min_truth_id].id();
      int status = evt.truth()[min_truth_id].status();
      int muon_type = evt.truth()[min_truth_id].muon_type();
      float eta = ma_particle.Eta(); 
      float phi = ma_particle.Phi(); 
      float deta_mu = eta - semilepmu.Eta();
      float dphi_mu = acos(cos(phi - semilepmu.Phi()));
      float dpt_mu = fabs(pt) - fabs(semilepmu.Pt()*1e-3);
      float ptvarcone30 = evt.inimu()[m].ptvarcone30();
      float iso = ptvarcone30 / evt.inimu()[m].pt();
      int accept = evt.inimu()[m].accept();
      float closejl_deltaR=999.;
      TLorentzVector closest_jet;
      for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
        float jet_pt = evt.jet()[jet_idx].mom().Perp();
        float jet_eta = evt.jet()[jet_idx].mom().Eta();
        float jet_phi = evt.jet()[jet_idx].mom().Phi();
        float jet_mass = evt.jet()[jet_idx].mom().M();
        std::cout << "jet pt/eta/phi=" << jet_pt << "/" << jet_eta << "/" << jet_phi << std::endl;
        float deta = evt.inimu()[m].eta() - jet_eta; 
        float dphi = acos(cos(evt.inimu()[m].phi() - jet_phi)); 
        float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
        if (deltaR_tmp < closejl_deltaR) {
          closejl_deltaR = deltaR_tmp;
          closest_jet.SetPtEtaPhiM(jet_pt,jet_eta,jet_phi,jet_mass);
        }
      }//for     
      //if(closejl_deltaR < 0.04+10./(evt.inimu()[m].pt()*1e-3)) continue;
      std::cout << "########################################" << std::endl;
      std::cout << "Matched to truth muon!!" << std::endl;
      std::cout << "iso=" << iso << std::endl;
      std::cout << "muon_type is " << muon_type << std::endl;
      std::cout << "pdgid/status/pt/eta/phi=" << pdgid << "/" << status << "/" << pt << "/" << eta << "/" << phi << std::endl;
      std::cout << "dif true mu_from_t and current true mu pt/eta/phi=" << dpt_mu << "/" << deta_mu << "/" << dphi_mu << std::endl;
      std::cout << "closest_jet p/pt/eta/phi=" << closest_jet.P() << "/" << closest_jet.Pt() << "/" << closest_jet.Eta() << "/" << closest_jet.Phi() << std::endl;
      std::cout << "accept=" << accept << std::endl;
      //if( closejl_deltaR < 0.04 + 10./(evt.inimu()[m].pt()*1e-3) ) continue;
      h->h1D(Form("initialMuon%s_pt_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(evt.inimu()[m].pt()*1e-3, weight);
      h->h1D(Form("initialMuon%s_quality_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(evt.inimu()[m].quality(), weight);
      h->h1D(Form("initialMuon_sd0_%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(fabs(evt.inimu()[m].d0sig()), weight);
      h->h1D(Form("initialMuon_dz0_%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(fabs(evt.inimu()[m].z0sintheta()), weight);
      h->h2D(Form("initialMuon_sd0_vs_dz0_%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(fabs(evt.inimu()[m].d0sig()),fabs(evt.inimu()[m].z0sintheta()), weight);
      h->h1D(Form("initialMuonPtvarcone30OverPt_%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(iso, weight);
      TLorentzVector initial_muon,static_muon;
      initial_muon.SetPtEtaPhiM(evt.inimu()[m].pt(),evt.inimu()[m].eta(),evt.inimu()[m].phi(),105.6);
      static_muon = initial_muon - closest_jet;
      float dphi_jet_muon = initial_muon.DeltaPhi(closest_jet); 
      std::cout << "dphi_jet_muon=" << dphi_jet_muon << std::endl;
      float initial_muon_p = initial_muon.P();
      float muon_pt_to_jet = fabs(initial_muon_p*sin(dphi_jet_muon));
      std::cout << "muon pt_to_jet=" << muon_pt_to_jet << std::endl;
      float PX,PY,PZ;
      float closest_jet_theta = closest_jet.Theta();
      float closest_jet_phi = closest_jet.Phi();
      float closest_jet_px = closest_jet.Px();
      float closest_jet_py = closest_jet.Py();
      float closest_jet_pz = closest_jet.Pz();
      float closest_jet_p = closest_jet.P();
      std::cout << "real px/py/pz=" << closest_jet_px << "/" << closest_jet_py << "/" << closest_jet_pz << std::endl; 
      //std::cout << "calc px/py/pz=" << closest_jet_p*cos(closest_jet_theta)*cos(closest_jet_phi) << "/" 
      //<< closest_jet_p*cos(closest_jet_theta)*sin(closest_jet_phi) << "/" << closest_jet_p*sin(closest_jet_theta) << std::endl; 
      std::cout << "calc px/py/pz=" << closest_jet_p*sin(closest_jet_theta)*cos(closest_jet_phi) << "/" 
        << closest_jet_p*sin(closest_jet_theta)*sin(closest_jet_phi) << "/" << closest_jet_p*cos(closest_jet_theta) << std::endl; 
      KinematicUtils::rotateXaxis(PX, PY, PZ, closest_jet_px, closest_jet_py, closest_jet_pz, closest_jet_theta, closest_jet_phi);
      std::cout << "after rotation px/py/pz=" << PX << "/" << PY << "/" << PZ << std::endl;
      std::cout << "static_muon pt/P/E=" << static_muon.Perp() << "/" << static_muon.P() << "/" << static_muon.E()  << std::endl;;
      std::cout << "closejl_deltaR/thr=" << closejl_deltaR << "/" << 0.04 + 10./(evt.inimu()[m].pt()*1e-3) << std::endl;
      std::cout << "no-isolated muon pt/eta/phi=" 
        << evt.inimu()[m].pt() << "/" << evt.inimu()[m].eta() << "/" << evt.inimu()[m].phi() << std::endl;
      std::cout << "quality=" << evt.inimu()[m].quality() << std::endl;
      std::cout << "d0sig/z0sintheta=" << evt.inimu()[m].d0sig() << "/" << evt.inimu()[m].z0sintheta() << std::endl;
      h->h2D(Form("lepPt_vs_drMuonJet_%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(evt.inimu()[m].pt()*1e-3, closejl_deltaR);
      h->h1D(Form("initialMuonPtToJet%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(muon_pt_to_jet*1e-3, weight);

      std::vector<TLorentzVector*> vec_nu_tmp = m_neutrinoBuilder.candidatesFromWMass_Rotation(&initial_muon, evt.met().Perp(), evt.met().Phi(), true);   
      TLorentzVector nu_tmp(0,0,0,0);
      if (vec_nu_tmp.size() > 0) {
        nu_tmp = *(vec_nu_tmp[0]);
        for (size_t z = 0; z < vec_nu_tmp.size(); ++z) delete vec_nu_tmp[z];
        vec_nu_tmp.clear();
      }
      TLorentzVector W_tmp = initial_muon + nu_tmp;
      std::cout << "W pt=" << W_tmp.Perp() << std::endl;
      float dphi_W_muon = initial_muon.DeltaPhi(W_tmp); 
      std::cout << "dphi_W_muon=" << dphi_W_muon << std::endl;
      float muon_pt_to_W = fabs(initial_muon_p*sin(dphi_W_muon));
      std::cout << "muon pt_to_W=" << muon_pt_to_W << std::endl;
      h->h1D(Form("initialMuonPtToW%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(muon_pt_to_W*1e-3, weight);
      h->h2D(Form("initialMuonPtToJet_vs_ToW%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(muon_pt_to_jet*1e-3,muon_pt_to_W*1e-3, weight);
      if(m==muid){//leading
        leading_muon_type = muon_type;
        h->h1D(Form("leadingInitialMuon%s_pt_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(evt.inimu()[m].pt()*1e-3, weight);
        h->h1D(Form("leadingInitialMuon_sd0_%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(fabs(evt.inimu()[m].d0sig()), weight);
        h->h1D(Form("leadingInitialMuon_dz0_%s_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(fabs(evt.inimu()[m].z0sintheta()), weight);
        h->h1D(Form("leadingInitialMuon%s_quality_%s",type[muon_type].c_str(),channel.c_str()), "", suffix)->Fill(evt.inimu()[m].quality(), weight);
      }
    }
    else {
      std::cout << "NOT mathed to muon!!" << std::endl;
      if(m==muid) leading_muon_type = 2;
    }
  }

  //Number of large jet
  int nLargeJet=0;
  int nGoodLargeJets50=0;
  int ljetid=-1;
  //std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
  for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
    const TLorentzVector &lj = evt.largeJet()[ljid].mom();
    bool good = evt.largeJet()[ljid].good();//default
    bool good50 = evt.largeJet()[ljid].good50();//50% WP
    //std::cout << "large jet pt/eta/good=" << lj.Perp()*1e-3 << "/" << lj.Eta() << "/" << good << std::endl;
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0 && good) {
      nLargeJet++;
      if(ljetid=-1) ljetid=ljid;
    }
    if(lj.Perp()*1e-3 > 300 && good50) nGoodLargeJets50++;
    h->h1D("largeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    h->h1D("largeJetEta", "", suffix)->Fill(lj.Eta(), weight);
    if(ljid==0) {
      h->h1D("leadingLargeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
      h->h1D("leadingLargeJetEta", "", suffix)->Fill(lj.Eta(), weight);
    }
  }
  std::cout << "nGoodLargeJets50=" << nGoodLargeJets50 << std::endl;
  h->h1D("nLargeJets", "", suffix)->Fill(nLargeJet, weight);
  if(nGoodLargeJets50<2) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(5);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(5);
  std::cout << "Pass #Large Jet" << std::endl;

  //Jet Cleaning
  //std::cout << "jet clean is " << evt.passes("jet_clean") << std::endl;
  if(!(evt.passes("jet_clean"))) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(6);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(6);
  std::cout << "Pass Jet Clean" << std::endl;

  //Track B-tag
  int nTrackBtagJets=0;
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    float pt = evt.tjet()[jidx].mom().Pt()*1e-3;
    float eta = evt.tjet()[jidx].mom().Eta();
    float phi = evt.tjet()[jidx].mom().Phi();
    //std::cout << "tjet pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl; 
    //std::cout << "b-tag value = " << evt.tjet()[jidx].mv2c10() << std::endl;
    bool is_trackBtag = evt.tjet()[jidx].btag_mv2c10_70_trk();
    if (is_trackBtag && evt.tjet()[jidx].mom().Pt() > 10e3 && 
        fabs(evt.tjet()[jidx].mom().Eta()) < 2.5 && evt.tjet()[jidx].numConstituents() >= 2) nTrackBtagJets++; 
  }
  //std::cout << "Number of track b-tagged jets is " << nTrackBtagJets << std::endl;
  h->h1D("nTrkBtagJets", "", suffix)->Fill(nTrackBtagJets, weight);  

  if(nTrackBtagJets==0) return;
  h->h1D("cutFlowAll", "", suffix)->Fill(7);
  h->h1D(Form("cutFlow%s",channel.c_str()), "", suffix)->Fill(7);
  std::cout << "Pass trk B-tag" << std::endl;

  h->h1D("nInitialMuonFinal", "", suffix)->Fill(nIniMuon,weight);

  /*std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), true);   
  TLorentzVector nu(0,0,0,0);
  if (vec_nu.size() > 0) {
    nu = *(vec_nu[0]);
    for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
    vec_nu.clear();
  }
  const TLorentzVector &ljt = evt.largeJet()[ljetid].mom();
  float mtt = (ljt+sj+nu+l).M();
  h->h1D("mttAll", "", suffix)->Fill(mtt*1e-3, weight);
  h->h1D(Form("mtt%s",channel.c_str()), "", suffix)->Fill(mtt*1e-3, weight);
  std::cout << "leading_muon_type=" << leading_muon_type << std::endl;
  h->h1D(Form("mtt%s_%s",type[leading_muon_type].c_str(),channel.c_str()), "", suffix)->Fill(mtt*1e-3, weight);
  if(channel=="SemiLepMu") {
    float deta = evt.inimu()[muid].eta() - semilepmu.Eta(); 
    float dphi = acos(cos(evt.inimu()[muid].phi() - semilepmu.Phi())); 
    float dr = sqrt(deta*deta+dphi*dphi);
  }
  std::cout << "pass_ljet=" << pass_ljet << std::endl;
  std::cout << "final nGoodLargeJets50=" << nGoodLargeJets50 << std::endl;
  if(channel=="FullHad") h->h1D("nGoodLargeJets50", "", suffix)->Fill(nGoodLargeJets50,weight);
  */

}

