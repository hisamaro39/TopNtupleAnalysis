/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtres21.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"
#include "TopNtupleAnalysis/KinematicUtils.h"


AnaTtres21::AnaTtres21(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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
    float pt_range[] = { 250000.000,325000.000,375000.000,425000.000,475000.000,525000.000,575000.000,625000.000,675000.000,725000.000,775000.000,850000.000,950000.000,1100000.000,1300000.000,1680000.000, 1e10 } ;
    int num_pt_range = sizeof(pt_range)/sizeof(float);

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

    std::string channel[6] = {"All","FullHad","FullLep","SemiLepTau","SemiLepMu","SemiLepEl"};
    std::string region[2] = {"Barrel","Endcap"};

      m_hSvc.create1D("EventNumber", "; Event number; Events", 500,10000000,400000000);
    for(int ch=0;ch<6;ch++){
      ////////////////////plot for validation
      m_hSvc.create1D(Form("NPV_%s",channel[ch].c_str()), "; # of PV; Events", 25,0,25);
      m_hSvc.create1D(Form("VtxZ_%s",channel[ch].c_str()), "; vtx z [mm]; Events", 50,-100,100);
      m_hSvc.create1D(Form("mu_%s",channel[ch].c_str()), "; #mu; Events", 50,0,50);
      
      m_hSvc.create1D(Form("MET_%s",channel[ch].c_str()), "; missing E_{T} [GeV]; Events", 100,0,1000);
      m_hSvc.create1D(Form("METPhi_%s",channel[ch].c_str()), "; missing E_{T} #phi ; Events", 32,-3.2,3.2);

      m_hSvc.create1D(Form("smallJetPt_%s",channel[ch].c_str()), "; small jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("smallJetEta_%s",channel[ch].c_str()), "; small jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("smallJetPhi_%s",channel[ch].c_str()), "; small jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("smallJetMass_%s",channel[ch].c_str()), "; small jet mass [GeV]; Events", 30, 0, 300);
      m_hSvc.create1D(Form("smallJetEtaHighPt_%s",channel[ch].c_str()), "; small jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("smallJetPhiHighPt_%s",channel[ch].c_str()), "; small jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("smallJetMassHighPt_%s",channel[ch].c_str()), "; small jet mass [GeV]; Events", 30, 0, 300);
      m_hSvc.create1D(Form("smallJetPtBoost_%s",channel[ch].c_str()), "; small jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("smallJetEtaBoost_%s",channel[ch].c_str()), "; small jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("smallJetPhiBoost_%s",channel[ch].c_str()), "; small jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("smallJetMassBoost_%s",channel[ch].c_str()), "; small jet mass [GeV]; Events", 30, 0, 300);
      m_hSvc.create1D(Form("nSmallJet_%s",channel[ch].c_str()), ";# of small jet;", 10,0,10);
      m_hSvc.create1D(Form("nBtaggedSmallJet_%s",channel[ch].c_str()), ";# of b-tagged small jet;", 10,0,10);
      m_hSvc.create1D(Form("nTrackJet_%s",channel[ch].c_str()), ";# of track jet;", 20,0,20);
      m_hSvc.create1D(Form("nNocalibTrackJet_%s",channel[ch].c_str()), ";# of track jet;", 20,0,20);
      m_hSvc.create1D(Form("nBtaggedTrackJet_%s",channel[ch].c_str()), ";# of b-tagged track jet;", 6,0,6);
      m_hSvc.create1D(Form("nBtaggedNocalibTrackJet_%s",channel[ch].c_str()), ";# of b-tagged track jet;", 6,0,6);
      m_hSvc.create1D(Form("smallJetJVTInRange_%s",channel[ch].c_str()), "; small jet JVT in range; Events", 20, 0, 1);
      m_hSvc.create1D(Form("smallJetJVTOutRange_%s",channel[ch].c_str()), "; small jet JVT out range; Events", 20, 0, 1);
      m_hSvc.create1D(Form("leadingSmallJetPt_%s",channel[ch].c_str()), "; leading small jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("leadingSmallJetEta_%s",channel[ch].c_str()), "; leading small jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("leadingSmallJetPhi_%s",channel[ch].c_str()), "; leading small jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("leadingSmallJetMass_%s",channel[ch].c_str()), "; leading small jet mass [GeV]; Events", 30, 0, 300);

      m_hSvc.create1D(Form("nocalibSmallJetPt_%s",channel[ch].c_str()), "; small jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("nocalibSmallJetPtNoJVT_%s",channel[ch].c_str()), "; small jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("nocalibSmallJetEta_%s",channel[ch].c_str()), "; small jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("nocalibSmallJetPhi_%s",channel[ch].c_str()), "; small jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("nocalibSmallJetMass_%s",channel[ch].c_str()), "; small jet mass [GeV]; Events", 30, 0, 300);
      m_hSvc.create1D(Form("nocalibSmallJetMassNoJVT_%s",channel[ch].c_str()), "; small jet mass [GeV]; Events", 30, 0, 300);
      m_hSvc.create1D(Form("nNocalibSmallJet_%s",channel[ch].c_str()), ";# of small jet;", 10,0,10);
      m_hSvc.create1D(Form("nocalibSmallJetJVTInRange_%s",channel[ch].c_str()), "; small jet JVT in range; Events", 20, 0, 1);
      m_hSvc.create1D(Form("nocalibSmallJetJVTOutRange_%s",channel[ch].c_str()), "; small jet JVT out range; Events", 20, 0, 1);
      m_hSvc.create1D(Form("nocalibSmallJetMinNoCalibElectronDr_%s",channel[ch].c_str()), "; min #DeltaR no selection jet & no selection electron;", 50,0,2);
      m_hSvc.create1D(Form("leadingNocalibSmallJetPt_%s",channel[ch].c_str()), "; leading no selection small jet p_{T} [GeV]; Events", 50, 0, 500);

      m_hSvc.create1D(Form("trackJetPt_%s",channel[ch].c_str()), "; track jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("trackJetEta_%s",channel[ch].c_str()), "; track jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("trackJetPhi_%s",channel[ch].c_str()), "; track jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("trackJetMv2c10_%s",channel[ch].c_str()), "; MV2c10 value of track jet; Events", 40, -1, 1);
      m_hSvc.create2D(Form("trackJetPtVsMv2c10_%s",channel[ch].c_str()), "track jet p_{T} [GeV]; MV2c10 value of track jet", 50, 0, 500, 40, -1, 1);
      m_hSvc.create1D(Form("trackJetNumConstituents_%s",channel[ch].c_str()), "; # of tracks in jet; Events", 20, 0, 20);
      m_hSvc.create1D(Form("leadingTrackJetPt_%s",channel[ch].c_str()), "; leading track jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("leadingTrackJetEta_%s",channel[ch].c_str()), "; leading track jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("leadingTrackJetPhi_%s",channel[ch].c_str()), "; leading track jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("btaggedTrackJetPt_%s",channel[ch].c_str()), "; b-tagged track jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("btaggedTrackJetEta_%s",channel[ch].c_str()), "; b-tagged track jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("btaggedTrackJetPhi_%s",channel[ch].c_str()), "; b-tagged track jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("trackBjetPt_%s",channel[ch].c_str()), "; b matched track jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("trackBjetEta_%s",channel[ch].c_str()), "; b matched track jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("trackBjetPhi_%s",channel[ch].c_str()), "; b matched track jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("trackBjetMv2c10_%s",channel[ch].c_str()), "; MV2c10 value of track jet; Events", 40, -1, 1);
      m_hSvc.create1D(Form("trackLightJetMv2c10_%s",channel[ch].c_str()), "; MV2c10 value of track jet; Events", 40, -1, 1);
      m_hSvc.create1D(Form("trackBjetNumConstituents_%s",channel[ch].c_str()), "; # of tracks in jet; Events", 20, 0, 20);
      m_hSvc.create1D(Form("btaggedTrackBjetPt_%s",channel[ch].c_str()), "; b matched track jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("btaggedTrackBjetEta_%s",channel[ch].c_str()), "; b matched track jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("btaggedTrackBjetPhi_%s",channel[ch].c_str()), "; b matched track jet #phi ; Events", 32, -3.2, 3.2);

      m_hSvc.create1D(Form("nocalibTrackJetPt_%s",channel[ch].c_str()), "; track jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("nocalibTrackJetEta_%s",channel[ch].c_str()), "; track jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("nocalibTrackJetPhi_%s",channel[ch].c_str()), "; track jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("nocalibTrackJetMv2c10_%s",channel[ch].c_str()), "; MV2c10 value of track jet; Events", 40, -1, 1);
      m_hSvc.create1D(Form("nocalibTrackJetNumConstituents_%s",channel[ch].c_str()), "; # of tracks in jet; Events", 20, 0, 20);
      m_hSvc.create1D(Form("nocalibTrackJetPt_def_%s",channel[ch].c_str()), "; track jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("nocalibTrackJetEta_def_%s",channel[ch].c_str()), "; track jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("nocalibTrackJetPhi_def_%s",channel[ch].c_str()), "; track jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("btaggedNocalibTrackJetPt_%s",channel[ch].c_str()), "; b-tagged track jet p_{T} [GeV]; Events", 50, 0, 500);
      m_hSvc.create1D(Form("btaggedNocalibTrackJetEta_%s",channel[ch].c_str()), "; b-tagged track jet #eta ; Events", 30,-3,3);
      m_hSvc.create1D(Form("btaggedNocalibTrackJetPhi_%s",channel[ch].c_str()), "; b-tagged track jet #phi ; Events", 32, -3.2, 3.2);

      m_hSvc.create1D(Form("nocalibLargeJetPt_%s",channel[ch].c_str()), "; large jet p_{T} [GeV]; Events", 500, 0, 500);
      m_hSvc.create1D(Form("nocalibLargeJetEta_%s",channel[ch].c_str()), "; large jet #eta ; Events", 25,-2.5,2.5);
      m_hSvc.create1D(Form("nocalibLargeJetPhi_%s",channel[ch].c_str()), "; large jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("nocalibLargeJetMass_%s",channel[ch].c_str()), "; large jet mass [GeV]; Events", 30, 0, 300);
      m_hSvc.create1D(Form("nNocalibLargeJet_%s",channel[ch].c_str()), ";# of large jet;", 10,0,10);
      m_hSvc.create1D(Form("leadingNocalibLargeJetPt_%s",channel[ch].c_str()), "; leading no selection large jet p_{T} [GeV]; Events", 100, 0, 2000);

      m_hSvc.create1D(Form("largeJetPt_%s",channel[ch].c_str()), "; large jet p_{T} [GeV] ; Events", 100,0,2000);
      m_hSvc.create1D(Form("largeJetEta_%s",channel[ch].c_str()), "; large jet #eta ; Events", 25, -2.5, 2.5);
      m_hSvc.create1D(Form("largeJetPhi_%s",channel[ch].c_str()), "; large jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("largeJetMass_%s",channel[ch].c_str()), "; large jet mass [GeV] ; Events", 30, 0, 300);
      m_hSvc.create1D(Form("largeJetEtaHighPt_%s",channel[ch].c_str()), "; large jet #eta ; Events", 25, -2.5, 2.5);
      m_hSvc.create1D(Form("largeJetPhiHighPt_%s",channel[ch].c_str()), "; large jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("largeJetMassHighPt_%s",channel[ch].c_str()), "; large jet mass [GeV] ; Events", 30, 0, 300);
      m_hSvc.create1D(Form("largeJetPtBoost_%s",channel[ch].c_str()), "; large jet p_{T} [GeV] ; Events", 100,0,2000);
      m_hSvc.create1D(Form("largeJetEtaBoost_%s",channel[ch].c_str()), "; large jet #eta ; Events", 25, -2.5, 2.5);
      m_hSvc.create1D(Form("largeJetPhiBoost_%s",channel[ch].c_str()), "; large jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("largeJetMassBoost_%s",channel[ch].c_str()), "; large jet mass [GeV] ; Events", 30, 0, 300);
      m_hSvc.create1D(Form("nLargeJet_%s",channel[ch].c_str()), ";# of large jet;", 7,0,7);
      m_hSvc.create1D(Form("nGoodLargeJet_%s",channel[ch].c_str()), ";# of good large jet;", 7,0,7);
      m_hSvc.create1D(Form("largeJetTau32_%s",channel[ch].c_str()), "; large jet tau32 ; Events", 100,0,1);
      m_hSvc.create1D(Form("leadingLargeJetPt_%s",channel[ch].c_str()), "; leading large jet p_{T} [GeV] ; Events", 100,0,2000);
      m_hSvc.create1D(Form("leadingLargeJetEta_%s",channel[ch].c_str()), "; leading large jet #eta ; Events", 25, -2.5, 2.5);
      m_hSvc.create1D(Form("leadingLargeJetPhi_%s",channel[ch].c_str()), "; leading large jet #phi ; Events", 32, -3.2, 3.2);
      m_hSvc.create1D(Form("leadingLargeJetMass_%s",channel[ch].c_str()), "; leading large jet mass [GeV] ; Events", 30, 0, 300);

      m_hSvc.create1D(Form("muonPt_%s",channel[ch].c_str()), "; muon p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("muonEta_%s",channel[ch].c_str()), "; muon #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("muonPhi_%s",channel[ch].c_str()), ";muon #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("nMuon_%s",channel[ch].c_str()), ";# of #mu;", 4,0,4);
      m_hSvc.create1D(Form("muonIsoValue_%s",channel[ch].c_str()), "; ptvarcone30/pt;", 20,0,0.2);
      m_hSvc.create1D(Form("muonSd0_%s",channel[ch].c_str()), "; sd0;", 50,0,5);
      m_hSvc.create1D(Form("muonDz0_%s",channel[ch].c_str()), "; dz0 [mm];", 20,0,2);
      m_hSvc.create1D(Form("muonQuality_%s",channel[ch].c_str()), "; quality;", 5,0,5);
      m_hSvc.create1D(Form("muonAccept_%s",channel[ch].c_str()), "; accept;", 2,0,2);
      m_hSvc.create1D(Form("muonMinJetDr_%s",channel[ch].c_str()), "; min #DeltaR muon jet;", 50,0,2);
      m_hSvc.create2D(Form("muonMinJetDrVsMuonPt_%s",channel[ch].c_str()), "; muon p_{T} [GeV];min #DeltaR muon & jet", 50,0,500,50,0,2);
      m_hSvc.create1D(Form("leadingMuonPt_%s",channel[ch].c_str()), "; leading muon p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("leadingMuonEta_%s",channel[ch].c_str()), "; leading muon #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("leadingMuonPhi_%s",channel[ch].c_str()), ";leading muon #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("muonPtResidual_%s",channel[ch].c_str()), "; residual muon p_{T} ;", 100,-0.4,0.4);

      m_hSvc.create1D(Form("nocalibMuonPt_%s",channel[ch].c_str()), "; muon p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("nocalibMuonEta_%s",channel[ch].c_str()), "; muon #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("nocalibMuonPhi_%s",channel[ch].c_str()), "; muon #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("nocalibMuonEtaHighPt_%s",channel[ch].c_str()), "; muon #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("nocalibMuonPhiHighPt_%s",channel[ch].c_str()), "; muon #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("nocalibMuonPtTightIso_%s",channel[ch].c_str()), "; muon p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("nocalibMuonEtaTightIso_%s",channel[ch].c_str()), "; muon #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("nocalibMuonPhiTightIso_%s",channel[ch].c_str()), "; muon #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("nocalibMuonPtBoost_%s",channel[ch].c_str()), "; muon p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("nocalibMuonEtaBoost_%s",channel[ch].c_str()), "; muon #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("nocalibMuonPhiBoost_%s",channel[ch].c_str()), "; muon #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("nNocalibMuon_%s",channel[ch].c_str()), ";# of muon;", 4,0,4);
      m_hSvc.create1D(Form("nocalibMuonIsoValue_%s",channel[ch].c_str()), "; ptvarcone30/pt;", 20,0,0.2);
      m_hSvc.create1D(Form("nocalibMuonSd0_%s",channel[ch].c_str()), "; sd0;", 50,0,5);
      m_hSvc.create1D(Form("nocalibMuonDz0_%s",channel[ch].c_str()), "; dz0 [mm];", 20,0,2);
      m_hSvc.create1D(Form("nocalibMuonQuality_%s",channel[ch].c_str()), "; quality;", 5,0,5);
      m_hSvc.create1D(Form("nocalibMuonAccept_%s",channel[ch].c_str()), "; accept;", 2,0,2);
      m_hSvc.create1D(Form("nocalibMuonMinJetDr_%s",channel[ch].c_str()), "; min #DeltaR nocalibMuon jet;", 50,0,2);
      m_hSvc.create2D(Form("nocalibMuonMinJetDrVsMuonPt_%s",channel[ch].c_str()), "; muon p_{T} [GeV];min #DeltaR muon & jet", 50,0,500,50,0,2);
      m_hSvc.create1D(Form("leadingNocalibMuonPt_%s",channel[ch].c_str()), "; leading muon p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("leadingNocalibMuonEta_%s",channel[ch].c_str()), "; leading muon #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("leadingNocalibMuonPhi_%s",channel[ch].c_str()), "; leading muon #phi;", 32,-3.2,3.2);

      for(int r=0;r<2;r++){
        m_hSvc.create1D(Form("muonPt_%s_%s",region[r].c_str(),channel[ch].c_str()), "; muon p_{T} [GeV];", 50,0,500);
        m_hSvc.create1D(Form("muonPhi_%s_%s",region[r].c_str(),channel[ch].c_str()), "; muon #phi;", 32,-3.2,3.2);
        m_hSvc.create1D(Form("nocalibMuonPt_%s_%s",region[r].c_str(),channel[ch].c_str()), "; muon p_{T [GeV]};", 50,0,500);
        m_hSvc.create1D(Form("nocalibMuonPhi_%s_%s",region[r].c_str(),channel[ch].c_str()), "; muon #phi;", 32,-3.2,3.2);
      }

      m_hSvc.create1D(Form("electronPt_%s",channel[ch].c_str()), "; electron p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("electronEta_%s",channel[ch].c_str()), ";electron #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("electronPhi_%s",channel[ch].c_str()), ";electron #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("electronEtaHighPt_%s",channel[ch].c_str()), ";electron #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("electronPhiHighPt_%s",channel[ch].c_str()), ";electron #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("electronPtBoost_%s",channel[ch].c_str()), "; electron p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("electronEtaBoost_%s",channel[ch].c_str()), ";electron #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("electronPhiBoost_%s",channel[ch].c_str()), ";electron #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("nElectron_%s",channel[ch].c_str()), ";# of electron;", 4,0,4);
      m_hSvc.create1D(Form("electronMinJetDr_%s",channel[ch].c_str()), "; min #DeltaR electron jet;", 50,0,2);
      m_hSvc.create2D(Form("electronMinJetDrVsElectronPt_%s",channel[ch].c_str()), "; electron p_{T} [GeV];min #DeltaR electron jet", 50,0,500,50,0,2);
      m_hSvc.create1D(Form("electronSd0_%s",channel[ch].c_str()), "; sd0;", 50,0,5);
      m_hSvc.create1D(Form("electronDz0_%s",channel[ch].c_str()), "; dz0 [mm];", 20,0,2);
      m_hSvc.create1D(Form("electronIsoValue_%s",channel[ch].c_str()), "; ptvarcone20/pt;", 20,0,0.2);
      m_hSvc.create1D(Form("leadingElectronPt_%s",channel[ch].c_str()), "; leading electron p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("leadingElectronEta_%s",channel[ch].c_str()), "; leading electron #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("leadingElectronPhi_%s",channel[ch].c_str()), ";leading electron #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("electronPtResidual_%s",channel[ch].c_str()), "; residual electron p_{T};", 100,-0.2,0.2);

      m_hSvc.create1D(Form("nocalibElectronPt_%s",channel[ch].c_str()), "; electron p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("nocalibElectronEta_%s",channel[ch].c_str()), "; electron #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("nocalibElectronPhi_%s",channel[ch].c_str()), "; electron #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("nocalibElectronPtTightIso_%s",channel[ch].c_str()), "; electron p_{T} [GeV];", 50,0,500);
      m_hSvc.create1D(Form("nocalibElectronEtaTightIso_%s",channel[ch].c_str()), "; electron #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("nocalibElectronPhiTightIso_%s",channel[ch].c_str()), "; electron #phi;", 32,-3.2,3.2);
      m_hSvc.create1D(Form("nNocalibElectron_%s",channel[ch].c_str()), ";# of nocalibrated electron;", 4,0,4);
      m_hSvc.create1D(Form("nocalibElectronSd0_%s",channel[ch].c_str()), "; sd0;", 50,0,5);
      m_hSvc.create1D(Form("nocalibElectronDz0_%s",channel[ch].c_str()), "; dz0 [mm];", 20,0,2);
      m_hSvc.create2D(Form("electronMinJetDrVsNocalibElectronPt_%s",channel[ch].c_str()), "; nocalibrated electron p_{T} [GeV];min #DeltaR electron jet", 50,0,500,50,0,2);
      m_hSvc.create1D(Form("leadingNocalibElectronPt_%s",channel[ch].c_str()), "; leading no selection electron p_{T} [GeV];", 50,0,500);

      m_hSvc.create1D(Form("trueLeptonPt_%s",channel[ch].c_str()), "; true lepton p_{T} [GeV];", 50,0,1000);
      m_hSvc.create1D(Form("trueLeptonEta_%s",channel[ch].c_str()), "; true lepton #eta;", 25,-2.5,2.5);
      m_hSvc.create1D(Form("trueLeptonPhi_%s",channel[ch].c_str()), "; true lepton #phi;", 32,-3.2,3.2);
      /////////////////////////

      for(int c=0;c<6;c++){
        m_hSvc.create1D(Form("muonPt_cut%d_%s",c+1,channel[ch].c_str()), "; muon p_{T} [GeV];", 20,0,500);
        m_hSvc.create1D(Form("muonEta_cut%d_%s",c+1,channel[ch].c_str()), "; muon #eta;", 25,-2.5,2.5);
        m_hSvc.create1D(Form("muonPhi_cut%d_%s",c+1,channel[ch].c_str()), ";muon #phi;", 32,-3.2,3.2);
        m_hSvc.create1D(Form("electronPt_cut%d_%s",c+1,channel[ch].c_str()), "; electron p_{T} [GeV];", 50,0,1000);
        m_hSvc.create1D(Form("electronEta_cut%d_%s",c+1,channel[ch].c_str()), "; electron #eta;", 25,-2.5,2.5);
        m_hSvc.create1D(Form("electronPhi_cut%d_%s",c+1,channel[ch].c_str()), ";electron #phi;", 32,-3.2,3.2);
        m_hSvc.create1D(Form("smallJetPt_cut%d_%s",c+1,channel[ch].c_str()), "; small jet p_{T} [GeV]; Events", 50, 0, 500);
        m_hSvc.create1D(Form("smallJetEta_cut%d_%s",c+1,channel[ch].c_str()), "; small jet #eta ; Events", 30,-3,3);
        m_hSvc.create1D(Form("smallJetPhi_cut%d_%s",c+1,channel[ch].c_str()), "; small jet #phi ; Events", 32, -3.2, 3.2);
        m_hSvc.create1D(Form("smallJetMass_cut%d_%s",c+1,channel[ch].c_str()), "; small jet mass [GeV]; Events", 30, 0, 300);
        m_hSvc.create1D(Form("largeJetPt_cut%d_%s",c+1,channel[ch].c_str()), "; large jet p_{T} [GeV]; Events", 100, 0, 2000);
        m_hSvc.create1D(Form("largeJetEta_cut%d_%s",c+1,channel[ch].c_str()), "; large jet #eta ; Events", 30,-3,3);
        m_hSvc.create1D(Form("largeJetPhi_cut%d_%s",c+1,channel[ch].c_str()), "; large jet #phi ; Events", 32, -3.2, 3.2);
        m_hSvc.create1D(Form("largeJetMass_cut%d_%s",c+1,channel[ch].c_str()), "; large jet mass [GeV]; Events", 30, 0, 300);
        m_hSvc.create1D(Form("trackJetPt_cut%d_%s",c+1,channel[ch].c_str()), "; track jet p_{T} [GeV]; Events", 50, 0, 500);
        m_hSvc.create1D(Form("trackJetEta_cut%d_%s",c+1,channel[ch].c_str()), "; track jet #eta ; Events", 30,-3,3);
        m_hSvc.create1D(Form("trackJetPhi_cut%d_%s",c+1,channel[ch].c_str()), "; track jet #phi ; Events", 32, -3.2, 3.2);
        m_hSvc.create1D(Form("trackJetMv2c10_cut%d_%s",c+1,channel[ch].c_str()), "; MV2c10 value of track jet; Events", 40, -1, 1);
        m_hSvc.create1D(Form("MET_cut%d_%s",c+1,channel[ch].c_str()), "; missing E_{T} [GeV]; Events", 100,0,1000);
        m_hSvc.create1D(Form("METPhi_cut%d_%s",c+1,channel[ch].c_str()), "; missing E_{T} #phi ; Events", 32,-3.2,3.2);
        m_hSvc.create1D(Form("nMuon_cut%d_%s",c+1,channel[ch].c_str()), ";# of #mu;", 4,0,4);
        m_hSvc.create1D(Form("nElectron_cut%d_%s",c+1,channel[ch].c_str()), ";# of electron;", 4,0,4);
        m_hSvc.create1D(Form("nSmallJet_cut%d_%s",c+1,channel[ch].c_str()), ";# of small jet;", 10,0,10);
        m_hSvc.create1D(Form("nTrackJet_cut%d_%s",c+1,channel[ch].c_str()), ";# of track jet;", 20,0,20);
        m_hSvc.create1D(Form("nocalibSmallJetPt_cut%d_%s",c+1,channel[ch].c_str()), "; small jet p_{T} [GeV]; Events", 50, 0, 500);
        m_hSvc.create1D(Form("nocalibSmallJetPtNoJVT_cut%d_%s",c+1,channel[ch].c_str()), "; small jet p_{T} [GeV]; Events", 50, 0, 500);
        m_hSvc.create1D(Form("muonPtResidual_cut%d_%s",c+1,channel[ch].c_str()), "; residual muon p_{T};", 100,-0.4,0.4);
        m_hSvc.create1D(Form("electronPtResidual_cut%d_%s",c+1,channel[ch].c_str()), "; residual electron p_{T};", 100,-0.2,0.2);
      }

      //plot for analysis
      m_hSvc.create1D(Form("nGoodLargeJet_ana_%s",channel[ch].c_str()), ";# of good large jet;", 7,0,7);
      m_hSvc.create1D(Form("nLargeJet_ana_%s",channel[ch].c_str()), ";# of large jet;", 7,0,7);
      m_hSvc.create1D(Form("largeJetPt_ana_%s",channel[ch].c_str()), ";large jet p_{T} [GeV];", 100,0,2000);
      m_hSvc.create1D(Form("largeJetEta_ana_%s",channel[ch].c_str()), ";large jet #eta;", 20,-2,2);
      m_hSvc.create1D(Form("largeJetMass_ana_%s",channel[ch].c_str()), ";large jet mass [GeV];", 30,0,300);
      m_hSvc.create1D(Form("largeJetGood_ana_%s",channel[ch].c_str()), ";large jet good;", 3,0,3);
      m_hSvc.create1D(Form("largeJetTau32_ana_%s",channel[ch].c_str()), "; tau32 ; Events", 100,0,1);
      m_hSvc.create1D(Form("goodLargeJetPt_ana_%s",channel[ch].c_str()), ";large jet p_{T} [GeV];", 100,0,2000);
      m_hSvc.create1D(Form("goodLargeJetEta_ana_%s",channel[ch].c_str()), ";large jet #eta;", 20,-2,2);
      m_hSvc.create1D(Form("goodLargeJetTau32_ana_%s",channel[ch].c_str()), "; tau32 ; Events", 100,0,1);
      m_hSvc.create1D(Form("goodLargeJetMass_ana_%s",channel[ch].c_str()), ";large jet mass [GeV];", 30,0,300);
      for(int p=0;p<num_pt_range-1;p++){
        m_hSvc.create1D(Form("goodLargeJetTau32_ana_%s_ptrange%d",channel[ch].c_str(),p), "; tau32 ; Events", 100,0,1);
        m_hSvc.create1D(Form("goodLargeJetMass_ana_%s_ptrange%d",channel[ch].c_str(),p), ";large jet mass [GeV];", 30,0,300);
      }
      m_hSvc.create1D(Form("leadingLargeJetTau32_ana_%s",channel[ch].c_str()), "; leading large jet tau32 ; Events", 100,0,1);

      m_hSvc.create1D(Form("cutFlow_%s",channel[ch].c_str()), ";Number of events; step", 20,0,20);
    }

  }

AnaTtres21::~AnaTtres21() {
}

void AnaTtres21::run(const Event &evt, double weight, const std::string &s, int is2016run) {

  double mu_weight[] = {1.25318,1.25318,1.25318,1.25318,1.26584,1.26584,1.26584,1.26584,0.241111,0.249024,0.592541,0.843912,1.21718,1.46061,0.989351,0.726253,0.641527,0.629048,0.723405,0.739132,0.754858,0.770584,0.723405,0.784223,0.776426,0.889655,0.819547,0.913004,1.28034,1.94749,3.70029,13.2429,11.6849,10.5164,7.78995,5.06346,4.28446,3.11597,2.72647,2.33697,999,999,999,999,999,999,999,999,999,999};
  int n_mu_weight = sizeof(mu_weight)/sizeof(double);
  //std::cout << "n_mu_weight=" << n_mu_weight << std::endl;
  //std::cout << "AnaTtres21::run" << std::endl;
  //std::cout << "weight=" << weight << std::endl;
  HistogramService *h = &m_hSvc;
  std::string suffix = s;

  //std::cout << "Event number is " << evt.eventNumber() << std::endl;
  h->h1D("EventNumber", "", suffix)->Fill(evt.eventNumber());

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
  TLorentzVector b_from_t,b_from_tbar,t,tbar;
  b_from_t = evt.MC_b_from_t();
  b_from_tbar = evt.MC_b_from_tbar();
  t = evt.MC_t();
  tbar = evt.MC_tbar();
  float t_pt = t.Perp()*1e-3, t_eta = t.Eta(), t_phi = t.Phi();
  float tbar_pt = tbar.Perp()*1e-3, tbar_eta = tbar.Eta(), tbar_phi = tbar.Phi();
  float b_from_t_pt = b_from_t.Perp()*1e-3, b_from_t_eta = b_from_t.Eta(), b_from_t_phi = b_from_t.Phi();
  float b_from_tbar_pt = b_from_tbar.Perp()*1e-3, b_from_tbar_eta = b_from_tbar.Eta(), b_from_tbar_phi = b_from_tbar.Phi();
  //std::cout << "isLeptonic top/tbar=" << isLeptonicT << "/" << isLeptonicTbar << std::endl;

  bool isFullHad = (!isLeptonicT && !isLeptonicTbar);
  bool isFullLep = (isLeptonicT && isLeptonicTbar);
  bool isSemiLep = (!isFullHad && !isFullLep)? true : false;
  bool containTau = (Wdecay2_from_t_pdgId==-15 || Wdecay1_from_tbar_pdgId==15)? true : false;
  bool containMu = (Wdecay2_from_t_pdgId==-13 || Wdecay1_from_tbar_pdgId==13)? true : false;
  bool containEl = (Wdecay2_from_t_pdgId==-11 || Wdecay1_from_tbar_pdgId==11)? true : false;

  TLorentzVector true_lepton;
  if(isLeptonicT) true_lepton = evt.MC_Wdecay2_from_t();
  else if(isLeptonicTbar) true_lepton = evt.MC_Wdecay1_from_tbar();

  std::string channel;
  if(isFullHad) channel="FullHad";
  else if(isFullLep) channel="FullLep";
  else if(containTau) channel="SemiLepTau";
  else if(containMu) channel="SemiLepMu";
  else if(containEl) channel="SemiLepEl";
  //std::cout << "channel is " << channel << std::endl;

  //truth tjet matching
  int bt_tj_idx=-1, btbar_tj_idx=-1;
  float min_dr_tj_bt = 999., min_dr_tj_btbar = 999.;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
    float pt = evt.tjet()[bidx].mom().Pt();
    float eta = evt.tjet()[bidx].mom().Eta();
    float phi = evt.tjet()[bidx].mom().Phi();
    float drt = TuDoAtlas::calc_delta_r(eta,b_from_t_eta,phi,b_from_t_phi);
    float drtbar = TuDoAtlas::calc_delta_r(eta,b_from_tbar_eta,phi,b_from_tbar_phi);
    int isBtag = evt.tjet()[bidx].btag_mv2c10_70_trk();
    if(drt<0.2 && drt<min_dr_tj_bt) {
      min_dr_tj_bt = drt;
      bt_tj_idx=bidx;
    }
    if(drtbar<0.2 && drtbar<min_dr_tj_btbar) {
      min_dr_tj_btbar = drtbar;
      btbar_tj_idx=bidx;
    }
  }

  h->h1D(Form("trueLeptonPt_%s",channel.c_str()), "", suffix)->Fill(true_lepton.Perp()*1e-3,weight);
  h->h1D(Form("trueLeptonEta_%s",channel.c_str()), "", suffix)->Fill(true_lepton.Eta(),weight);
  h->h1D(Form("trueLeptonPhi_%s",channel.c_str()), "", suffix)->Fill(true_lepton.Phi(),weight);

  //pre event selection
  //pre trigger
  bool passtrigmu = evt.trigger("HLT_mu26_ivarmedium") || evt.trigger("HLT_mu50");

  //pre muon
  //std::cout << "Number of muon is " << evt.muon().size() << std::endl;
  int pre_numMuon=0;
  int pre_lmu_id=-1;
  float pre_lmu_pt=0.;
  for(int im=0;im<evt.muon().size();im++){
    TLorentzVector mu = evt.muon()[im].mom();
    if(mu.Perp()*1e-3<30) continue;
    if(mu.Perp()>pre_lmu_pt){
      pre_lmu_id=im;
      pre_lmu_pt=mu.Perp();
    }
  }
  bool hasMuon = (pre_lmu_id!=-1)? true : false;

  //pre small jets
  int pre_nJet=0;
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    float jet_pt = evt.jet()[jidx].mom().Perp()*1e-3;
    if(jet_pt>25) pre_nJet++;
  }
  bool hasSmallJet = (pre_nJet)? true : false;

  //pre large jet
  bool hasHighPtLargeJet=false;
  bool hasTopLikeLargeJet=false;
  for (int ilj=0; ilj < evt.largeJet().size(); ++ilj) {
    TLorentzVector lj = evt.largeJet()[ilj].mom();
    if(lj.Perp()*1e-3<100) continue;
    if(fabs(lj.Eta())>2) continue;
    if(lj.Perp()*1e-3>300 && fabs(lj.Eta())<2.0) hasHighPtLargeJet=true;
    if(lj.M()*1e-3>140 && lj.M()*1e-3<190) hasTopLikeLargeJet=true;
  }

  //pre MET
  bool hasMET = (evt.met().Perp()*1e-3>20)? true : false;
  
  //pre kinematic value
  bool hasMETplusMWT=false,hasJetCloseToLepton=false,hasGoodLargeJet=false;
  if(hasMuon){
    //MWT
    TLorentzVector l;
    l.SetPtEtaPhiM(evt.muon()[pre_lmu_id].mom().Perp(),evt.muon()[pre_lmu_id].mom().Eta(),evt.muon()[pre_lmu_id].mom().Phi(),105.6);
    float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3;
    float MET = evt.met().Perp()*1e-3;
    float MET_plus_mWt = mWt + MET;
    if(MET_plus_mWt>60) hasMETplusMWT=true;

    //Jet close to lepton
    size_t close_idx = 0;
    for (; close_idx < evt.jet().size(); ++close_idx){    
      float dphi = l.DeltaPhi(evt.jet()[close_idx].mom());
      float dy = l.Rapidity() - evt.jet()[close_idx].mom().Rapidity();
      float deltaR_tmp = std::sqrt(dy*dy + dphi*dphi);
      if(evt.jet()[close_idx].mom().Pt()*1e-3 > 25 && deltaR_tmp < 1.5) break;
    }//for     
    if(close_idx!=evt.jet().size()) hasJetCloseToLepton=true;

    //Angular cut for large jet
    if(hasJetCloseToLepton){
      TLorentzVector sj = evt.jet()[close_idx].mom();
      for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
        const TLorentzVector &lj = evt.largeJet()[ljid].mom();
        float delta_phi_large_jet_lepton = lj.DeltaPhi(l);
        float delta_r_large_jet_small_jet = lj.DeltaR(sj);
        if(lj.Perp()*1e-3 < 300 || fabs(lj.Eta())>2.0 || lj.M()*1e-3<140 || lj.M()*1e-3>190) continue;
        if(fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) hasGoodLargeJet=true; 
      }
    }
  }

  int nCut=6;
  bool cut[nCut];
  cut[0] = hasHighPtLargeJet;
  cut[1] = hasTopLikeLargeJet;
  cut[2] = hasHighPtLargeJet && hasTopLikeLargeJet;
  cut[3] = passtrigmu && hasMuon; 
  cut[4] = hasMET && hasMETplusMWT;
  cut[5] = cut[3] && cut[4] && hasJetCloseToLepton && hasGoodLargeJet;

  if(!evt.passes("no_cut")) return;
  //if(!evt.passes("bmujets_2016")) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(0);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(0);

  bool passtrig_mu1 = evt.trigger("HLT_mu26_ivarmedium");
  bool passtrig_mu2 = evt.trigger("HLT_mu50");
  //if(!passtrig_mu1 && !passtrig_mu2) return;
  bool passtrig_el1 = evt.trigger("HLT_e26_lhtight_nod0_ivarloose"); 
  bool passtrig_el2 = evt.trigger("HLT_e60_lhmedium_nod0"); 
  bool passtrig_el3 = evt.trigger("HLT_e140_lhloose_nod0"); 
  //if(!passtrig_el1 && !passtrig_el2 && !passtrig_el3) return;

  //plot for validation
  //std::cout << "mu=" << evt.mu() << std::endl;
  //std::cout << "default weight=" << weight << std::endl;
  for(int b=0;b<n_mu_weight;b++){
    if(evt.mu()>b && evt.mu()<b+1) {
      //weight = weight / mu_weight[b];
      //std::cout << "new weight=" << weight << std::endl;
    }
  }

  h->h1D("NPV_All", "", suffix)->Fill(evt.npv(),weight);
  h->h1D("VtxZ_All", "", suffix)->Fill(evt.vtxz(),weight);
  h->h1D("mu_All", "", suffix)->Fill(evt.mu(),weight);

 
  //large jet
  int numLargeJet=0,numGoodLargeJet=0;
  //std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
  int llj_id=-1;
  float llj_pt=0.;
  for (int ilj=0; ilj < evt.largeJet().size(); ++ilj) {
    TLorentzVector lj = evt.largeJet()[ilj].mom();
    if(lj.Perp()>llj_pt){
      llj_id=ilj;
      llj_pt=lj.Perp();
    }
    //std::cout << "large jet pt/eta/phi=" << lj.Perp() << "/" << lj.Eta() << "/" << lj.Phi() << std::endl;
    if(lj.Perp()*1e-3<100) continue;
    if(fabs(lj.Eta())>2) continue;
    h->h1D("largeJetPt_All", "", suffix)->Fill(lj.Perp()*1e-3,weight);
    h->h1D("largeJetEta_All", "", suffix)->Fill(lj.Eta(),weight);
    h->h1D("largeJetPhi_All", "", suffix)->Fill(lj.Phi(),weight);
    h->h1D("largeJetMass_All", "", suffix)->Fill(lj.M()*1e-3,weight);
    h->h1D(Form("largeJetPt_%s",channel.c_str()), "", suffix)->Fill(lj.Perp()*1e-3,weight);
    h->h1D(Form("largeJetEta_%s",channel.c_str()), "", suffix)->Fill(lj.Eta(),weight);
    h->h1D(Form("largeJetPhi_%s",channel.c_str()), "", suffix)->Fill(lj.Phi(),weight);
    h->h1D(Form("largeJetMass_%s",channel.c_str()), "", suffix)->Fill(lj.M()*1e-3,weight);
    //h->h1D(Form("largeJetTau32_%s",channel.c_str()), "", suffix)->Fill(evt.largeJet()[ilj].subs("tau32_wta"),weight);
    numLargeJet++;
    for(int c=0;c<nCut;c++){
      if(cut[c]){
        h->h1D(Form("largeJetPt_cut%d_All",c+1), "", suffix)->Fill(lj.Perp()*1e-3,weight);
        h->h1D(Form("largeJetEta_cut%d_All",c+1), "", suffix)->Fill(lj.Eta(),weight);
        h->h1D(Form("largeJetPhi_cut%d_All",c+1), "", suffix)->Fill(lj.Phi(),weight);
        h->h1D(Form("largeJetMass_cut%d_All",c+1), "", suffix)->Fill(lj.M()*1e-3,weight);
        h->h1D(Form("largeJetPt_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(lj.Perp()*1e-3,weight);
        h->h1D(Form("largeJetEta_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(lj.Eta(),weight);
        h->h1D(Form("largeJetPhi_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(lj.Phi(),weight);
        h->h1D(Form("largeJetMass_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(lj.M()*1e-3,weight);
      }
    }
    if(lj.Perp()*1e-3>400 && fabs(lj.Eta())<2.0) {
      h->h1D("largeJetEtaHighPt_All", "", suffix)->Fill(lj.Eta(),weight);
      h->h1D("largeJetPhiHighPt_All", "", suffix)->Fill(lj.Phi(),weight);
      h->h1D("largeJetMassHighPt_All", "", suffix)->Fill(lj.M()*1e-3,weight);
      h->h1D(Form("largeJetEtaHighPt_%s",channel.c_str()), "", suffix)->Fill(lj.Eta(),weight);
      h->h1D(Form("largeJetPhiHighPt_%s",channel.c_str()), "", suffix)->Fill(lj.Phi(),weight);
      h->h1D(Form("largeJetMassHighPt_%s",channel.c_str()), "", suffix)->Fill(lj.M()*1e-3,weight);
    }
    if(lj.Perp()*1e-3>300 && fabs(lj.Eta())<2.0 && evt.largeJet()[ilj].good()) numGoodLargeJet++;
  }
  if(llj_id!=-1){
    TLorentzVector leading_lj = evt.largeJet()[llj_id].mom();
    h->h1D("leadingLargeJetPt_All", "", suffix)->Fill(leading_lj.Perp()*1e-3,weight);
    h->h1D("leadingLargeJetEta_All", "", suffix)->Fill(leading_lj.Eta(),weight);
    h->h1D("leadingLargeJetPhi_All", "", suffix)->Fill(leading_lj.Phi(),weight);
    h->h1D("leadingLargeJetMass_All", "", suffix)->Fill(leading_lj.M()*1e-3,weight);
    h->h1D(Form("leadingLargeJetPt_%s",channel.c_str()), "", suffix)->Fill(leading_lj.Perp()*1e-3,weight);
    h->h1D(Form("leadingLargeJetEta_%s",channel.c_str()), "", suffix)->Fill(leading_lj.Eta(),weight);
    h->h1D(Form("leadingLargeJetPhi_%s",channel.c_str()), "", suffix)->Fill(leading_lj.Phi(),weight);
    h->h1D(Form("leadingLargeJetMass_%s",channel.c_str()), "", suffix)->Fill(leading_lj.M()*1e-3,weight);
  }
  h->h1D("nLargeJet_All", "", suffix)->Fill(numLargeJet,weight);
  h->h1D(Form("nLargeJet_%s",channel.c_str()), "", suffix)->Fill(numLargeJet,weight);
  h->h1D("nGoodLargeJet_All", "", suffix)->Fill(numGoodLargeJet,weight);
  h->h1D(Form("nGoodLargeJet_%s",channel.c_str()), "", suffix)->Fill(numGoodLargeJet,weight);
  for(int c=0;c<nCut;c++){
    if(cut[c]){
      h->h1D(Form("nLargeJet_cut%d_All",c+1), "", suffix)->Fill(numLargeJet,weight);
      h->h1D(Form("nLargeJet_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(numLargeJet,weight);
    }
  }
  //std::cout << "numLargeJet=" << numLargeJet << std::endl;

  //muon
  //std::cout << "Number of muon is " << evt.muon().size() << std::endl;
  int numMuon=0,numMuon25=0;
  int lmu_id=-1;
  float lmu_pt=0.;
  for(int im=0;im<evt.muon().size();im++){
    TLorentzVector mu = evt.muon()[im].mom();
    if(mu.Perp()>lmu_pt){
      lmu_id=im;
      lmu_pt=mu.Perp();
    }
    float iso  = evt.muon()[im].ptvarcone30() / mu.Perp();
    //std::cout << "iso=" << iso << std::endl;
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float deta = mu.Eta() - evt.jet()[jet_idx].mom().Eta(); 
      float dphi = acos(cos(mu.Phi() - evt.jet()[jet_idx].mom().Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    std::string region = (fabs(mu.Eta())>1.05)? "Endcap" : "Barrel";
    h->h1D("muonPt_All", "", suffix)->Fill(mu.Perp()*1e-3,weight);
    h->h1D("muonEta_All", "", suffix)->Fill(mu.Eta(),weight);
    h->h1D("muonPhi_All", "", suffix)->Fill(mu.Phi(),weight);
    h->h1D("muonIsoValue_All", "", suffix)->Fill(iso,weight);
    h->h1D("muonSd0_All", "", suffix)->Fill(evt.muon()[im].sd0(),weight);
    h->h1D("muonDz0_All", "", suffix)->Fill(evt.muon()[im].Dz0(),weight);
    h->h1D(Form("muonPt_%s",channel.c_str()), "", suffix)->Fill(mu.Perp()*1e-3,weight);
    h->h1D(Form("muonEta_%s",channel.c_str()), "", suffix)->Fill(mu.Eta(),weight);
    h->h1D(Form("muonPhi_%s",channel.c_str()), "", suffix)->Fill(mu.Phi(),weight);
    h->h1D(Form("muonIsoValue_%s",channel.c_str()), "", suffix)->Fill(iso,weight);
    h->h1D(Form("muonQuality_%s",channel.c_str()), "", suffix)->Fill(evt.muon()[im].quality(),weight);
    h->h1D(Form("muonAccept_%s",channel.c_str()), "", suffix)->Fill(evt.muon()[im].accept(),weight);
    h->h1D(Form("muonMinJetDr_%s",channel.c_str()), "", suffix)->Fill(closejl_deltaR,weight);
    h->h2D(Form("muonMinJetDrVsMuonPt_%s",channel.c_str()), "", suffix)->Fill(mu.Perp()*1e-3,closejl_deltaR,weight);
    h->h1D(Form("muonSd0_%s",channel.c_str()), "", suffix)->Fill(evt.muon()[im].sd0(),weight);
    h->h1D(Form("muonDz0_%s",channel.c_str()), "", suffix)->Fill(evt.muon()[im].Dz0(),weight);
    h->h1D(Form("muonPhi_%s_%s",region.c_str(),channel.c_str()), "", suffix)->Fill(mu.Phi(),weight);
    for(int c=0;c<nCut;c++){
      if(cut[c]){
        h->h1D(Form("muonPt_cut%d_All",c+1), "", suffix)->Fill(mu.Perp()*1e-3,weight);
        h->h1D(Form("muonEta_cut%d_All",c+1), "", suffix)->Fill(mu.Eta(),weight);
        h->h1D(Form("muonPhi_cut%d_All",c+1), "", suffix)->Fill(mu.Phi(),weight);
        h->h1D(Form("muonPt_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(mu.Perp()*1e-3,weight);
        h->h1D(Form("muonEta_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(mu.Eta(),weight);
        h->h1D(Form("muonPhi_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(mu.Phi(),weight);
      }
    }
    if(mu.Perp()*1e-3>30) numMuon++;
    if(mu.Perp()*1e-3>25) numMuon25++;
    //truth matching
    float deta = mu.Eta() - true_lepton.Eta();
    float dphi = acos(cos(mu.Phi() - true_lepton.Phi()));
    float dr = sqrt(deta*deta+dphi*dphi);
    float pt_residual = (mu.Perp() - true_lepton.Perp())/fabs(true_lepton.Perp());
    if(dr<0.2){
      h->h1D("muonPtResidual_All", "", suffix)->Fill(pt_residual,weight);
      h->h1D(Form("muonPtResidual_%s",channel.c_str()), "", suffix)->Fill(pt_residual,weight);
      for(int c=0;c<nCut;c++){
        if(cut[c]){
          h->h1D(Form("muonPtResidual_cut%d_All",c+1), "", suffix)->Fill(pt_residual,weight);
          h->h1D(Form("muonPtResidual_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(pt_residual,weight);
        }
      }
    }
  }
  if(lmu_id!=-1){
    TLorentzVector leading_mu = evt.muon()[lmu_id].mom();
    h->h1D("leadingMuonPt_All", "", suffix)->Fill(leading_mu.Perp()*1e-3,weight);
    h->h1D("leadingMuonEta_All", "", suffix)->Fill(leading_mu.Eta(),weight);
    h->h1D("leadingMuonPhi_All", "", suffix)->Fill(leading_mu.Phi(),weight);
    h->h1D(Form("leadingMuonPt_%s",channel.c_str()), "", suffix)->Fill(leading_mu.Perp()*1e-3,weight);
    h->h1D(Form("leadingMuonEta_%s",channel.c_str()), "", suffix)->Fill(leading_mu.Eta(),weight);
    h->h1D(Form("leadingMuonPhi_%s",channel.c_str()), "", suffix)->Fill(leading_mu.Phi(),weight);
  }
  h->h1D("nMuon_All", "", suffix)->Fill(numMuon,weight);
  h->h1D(Form("nMuon_%s",channel.c_str()), "", suffix)->Fill(numMuon,weight);
  for(int c=0;c<nCut;c++){
    if(cut[c]){
      h->h1D(Form("nMuon_cut%d_All",c+1), "", suffix)->Fill(numMuon,weight);
      h->h1D(Form("nMuon_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(numMuon,weight);
    }
  }


  //no selection muon
  //std::cout << "Number of no selection muon is " << evt.nocalibmu().size() << std::endl;
  int numNocalibMuon=0,muon_id_nocalib=-1;
  int leadmu_id=-1;
  float leadmu_pt=0.;
  bool hasMuonCloseToJet=false;
  for(int im=0;im<evt.nocalibmu().size();im++){
    TLorentzVector mu = evt.nocalibmu()[im].mom();
    float pt = mu.Perp();
    float eta = mu.Eta();
    float phi = mu.Phi();
    //std::cout << "pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl;
    float ptvarcone30 = evt.nocalibmu()[im].ptvarcone30();
    float iso = ptvarcone30 / pt;
    int quality = evt.nocalibmu()[im].quality();
    int accept = evt.nocalibmu()[im].accept();
    float d0sig = evt.nocalibmu()[im].sd0();
    float z0sintheta = evt.nocalibmu()[im].Dz0();
    //std::cout << "iso/quality/accept/d0sig/z0sintheta=" << iso << "/" << quality << "/" << accept << "/" << d0sig << "/" << z0sintheta << std::endl;
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float deta = eta - evt.jet()[jet_idx].mom().Eta(); 
      float dphi = acos(cos(phi - evt.jet()[jet_idx].mom().Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    if(pt*1e-3<25) continue;
    if(fabs(eta)>2.5) continue;
    h->h1D("nocalibMuonIsoValue_All", "", suffix)->Fill(iso,weight);
    h->h1D("nocalibMuonQuality_All", "", suffix)->Fill(quality,weight);
    h->h1D("nocalibMuonAccept_All", "", suffix)->Fill(accept,weight);
    h->h1D("nocalibMuonMinJetDr_All", "", suffix)->Fill(closejl_deltaR,weight);
    h->h2D("nocalibMuonMinJetDrVsnocalibMuonPt_All", "", suffix)->Fill(mu.Perp()*1e-3,closejl_deltaR,weight);
    h->h1D("nocalibMuonSd0_All", "", suffix)->Fill(d0sig,weight);
    h->h1D("nocalibMuonDz0_All", "", suffix)->Fill(z0sintheta,weight);
    h->h1D(Form("nocalibMuonIsoValue_%s",channel.c_str()), "", suffix)->Fill(iso,weight);
    h->h1D(Form("nocalibMuonQuality_%s",channel.c_str()), "", suffix)->Fill(quality,weight);
    h->h1D(Form("nocalibMuonAccept_%s",channel.c_str()), "", suffix)->Fill(accept,weight);
    h->h1D(Form("nocalibMuonMinJetDr_%s",channel.c_str()), "", suffix)->Fill(closejl_deltaR,weight);
    h->h2D(Form("nocalibMuonMinJetDrVsnocalibMuonPt_%s",channel.c_str()), "", suffix)->Fill(mu.Perp()*1e-3,closejl_deltaR,weight);
    h->h1D(Form("nocalibMuonSd0_%s",channel.c_str()), "", suffix)->Fill(d0sig,weight);
    h->h1D(Form("nocalibMuonDz0_%s",channel.c_str()), "", suffix)->Fill(z0sintheta,weight);
    if(iso>0.06) continue;
    if(quality>1) continue;
    if(!accept) continue;
    if(d0sig>3) continue;
    if(z0sintheta>0.5) continue;
    if(closejl_deltaR<1.5) hasMuonCloseToJet=true;
    if(closejl_deltaR<0.04+10./(pt*1e-3)) continue;
    if(pt>leadmu_pt) {
      leadmu_pt = pt;
      leadmu_id=im;
    }
    //if(closejl_deltaR<0.4) continue;
    h->h1D("nocalibMuonPt_All", "", suffix)->Fill(mu.Perp()*1e-3,weight);
    h->h1D("nocalibMuonEta_All", "", suffix)->Fill(mu.Eta(),weight);
    h->h1D("nocalibMuonPhi_All", "", suffix)->Fill(mu.Phi(),weight);
    h->h1D(Form("nocalibMuonPt_%s",channel.c_str()), "", suffix)->Fill(mu.Perp()*1e-3,weight);
    h->h1D(Form("nocalibMuonEta_%s",channel.c_str()), "", suffix)->Fill(mu.Eta(),weight);
    h->h1D(Form("nocalibMuonPhi_%s",channel.c_str()), "", suffix)->Fill(mu.Phi(),weight);
    if(pt*1e-3>50){
      h->h1D("nocalibMuonEtaHighPt_All", "", suffix)->Fill(mu.Eta(),weight);
      h->h1D("nocalibMuonPhiHighPt_All", "", suffix)->Fill(mu.Phi(),weight);
      h->h1D(Form("nocalibMuonEtaHighPt_%s",channel.c_str()), "", suffix)->Fill(mu.Eta(),weight);
      h->h1D(Form("nocalibMuonPhiHighPt_%s",channel.c_str()), "", suffix)->Fill(mu.Phi(),weight);
    }
    if(iso<0.03) {
      h->h1D("nocalibMuonPtTightIso_All", "", suffix)->Fill(mu.Perp()*1e-3,weight);
      h->h1D("nocalibMuonEtaTightIso_All", "", suffix)->Fill(mu.Eta(),weight);
      h->h1D("nocalibMuonPhiTightIso_All", "", suffix)->Fill(mu.Phi(),weight);
      h->h1D(Form("nocalibMuonPtTightIso_%s",channel.c_str()), "", suffix)->Fill(mu.Perp()*1e-3,weight);
      h->h1D(Form("nocalibMuonEtaTightIso_%s",channel.c_str()), "", suffix)->Fill(mu.Eta(),weight);
      h->h1D(Form("nocalibMuonPhiTightIso_%s",channel.c_str()), "", suffix)->Fill(mu.Phi(),weight);
    }
    if(numNocalibMuon==0) muon_id_nocalib=im;
    numNocalibMuon++;
  }
  if(leadmu_id!=-1){
    TLorentzVector leading_mu = evt.nocalibmu()[leadmu_id].mom();
    h->h1D("leadingNocalibMuonPt_All", "", suffix)->Fill(leading_mu.Perp()*1e-3,weight);
    h->h1D("leadingNocalibMuonEta_All", "", suffix)->Fill(leading_mu.Eta(),weight);
    h->h1D("leadingNocalibMuonPhi_All", "", suffix)->Fill(leading_mu.Phi(),weight);
    h->h1D(Form("leadingNocalibMuonPt_%s",channel.c_str()), "", suffix)->Fill(leading_mu.Perp()*1e-3,weight);
    h->h1D(Form("leadingNocalibMuonEta_%s",channel.c_str()), "", suffix)->Fill(leading_mu.Eta(),weight);
    h->h1D(Form("leadingNocalibMuonPhi_%s",channel.c_str()), "", suffix)->Fill(leading_mu.Phi(),weight);
  }
  h->h1D("nNocalibMuon_All", "", suffix)->Fill(numNocalibMuon,weight);
  h->h1D(Form("nNocalibMuon_%s",channel.c_str()), "", suffix)->Fill(numNocalibMuon,weight);
  
  bool isBoost = hasHighPtLargeJet && hasMuonCloseToJet;
  //std::cout << "isBoost=" << isBoost << std::endl;

  //electron
  //std::cout << "Number of electron is " << evt.electron().size() << std::endl;
  int numElectron=0;
  int lel_id=-1;
  float lel_pt=0.;
  for(int ie=0;ie<evt.electron().size();ie++){
    TLorentzVector el = evt.electron()[ie].mom();
    if(el.Perp()>lel_pt){
      lel_id=ie;
      lel_pt=el.Perp();
    }
    //std::cout << "electron pt/eta/phi=" << el.Perp() << "/" << el.Eta() << "/" << el.Phi() << std::endl;
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float deta = el.Eta() - evt.jet()[jet_idx].mom().Eta(); 
      float dphi = acos(cos(el.Phi() - evt.jet()[jet_idx].mom().Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    if(el.Perp()*1e-3<30) continue;
    float iso  = evt.electron()[ie].ptvarcone20() / el.Perp();
    h->h1D("electronPt_All", "", suffix)->Fill(el.Perp()*1e-3,weight);
    h->h1D("electronEta_All", "", suffix)->Fill(el.Eta(),weight);
    h->h1D("electronPhi_All", "", suffix)->Fill(el.Phi(),weight);
    h->h1D("electronSd0_All", "", suffix)->Fill(evt.electron()[ie].sd0(),weight);
    h->h1D("electronDz0_All", "", suffix)->Fill(evt.electron()[ie].Dz0(),weight);
    h->h1D("electronIsoValue_All","", suffix)->Fill(iso,weight);
    h->h1D(Form("electronPt_%s",channel.c_str()), "", suffix)->Fill(el.Perp()*1e-3,weight);
    h->h1D(Form("electronEta_%s",channel.c_str()), "", suffix)->Fill(el.Eta(),weight);
    h->h1D(Form("electronPhi_%s",channel.c_str()), "", suffix)->Fill(el.Phi(),weight);
    h->h1D(Form("electronMinJetDr_%s",channel.c_str()), "", suffix)->Fill(closejl_deltaR,weight);
    h->h2D(Form("electronMinJetDrVsElectronPt_%s",channel.c_str()), "", suffix)->Fill(el.Perp()*1e-3,closejl_deltaR,weight);
    h->h1D(Form("electronSd0_%s",channel.c_str()), "", suffix)->Fill(evt.electron()[ie].sd0(),weight);
    h->h1D(Form("electronDz0_%s",channel.c_str()), "", suffix)->Fill(evt.electron()[ie].Dz0(),weight);
    //std::cout << "electron ptvarcone20/iso=" << evt.electron()[ie].ptvarcone20() << "/" << iso << std::endl;
    h->h1D(Form("electronIsoValue_%s",channel.c_str()), "", suffix)->Fill(iso,weight);
    if(el.Perp()*1e-3>50){
      h->h1D("electronEtaHighPt_All", "", suffix)->Fill(el.Eta(),weight);
      h->h1D("electronPhiHighPt_All", "", suffix)->Fill(el.Phi(),weight);
      h->h1D(Form("electronEtaHighPt_%s",channel.c_str()), "", suffix)->Fill(el.Eta(),weight);
      h->h1D(Form("electronPhiHighPt_%s",channel.c_str()), "", suffix)->Fill(el.Phi(),weight);
    }
    for(int c=0;c<nCut;c++){
      if(cut[c]){
        h->h1D(Form("electronPt_cut%d_All",c+1), "", suffix)->Fill(el.Perp()*1e-3,weight);
        h->h1D(Form("electronEta_cut%d_All",c+1), "", suffix)->Fill(el.Eta(),weight);
        h->h1D(Form("electronPhi_cut%d_All",c+1), "", suffix)->Fill(el.Phi(),weight);
        h->h1D(Form("electronPt_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(el.Perp()*1e-3,weight);
        h->h1D(Form("electronEta_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(el.Eta(),weight);
        h->h1D(Form("electronPhi_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(el.Phi(),weight);
      }
    }
    //truth matching
    float deta = el.Eta() - true_lepton.Eta();
    float dphi = acos(cos(el.Phi() - true_lepton.Phi()));
    float dr = sqrt(deta*deta+dphi*dphi);
    float pt_residual = (el.Perp() - true_lepton.Perp())/fabs(true_lepton.Perp());
    if(dr<0.2){
      h->h1D("electronPtResidual_All", "", suffix)->Fill(pt_residual,weight);
      h->h1D(Form("electronPtResidual_%s",channel.c_str()), "", suffix)->Fill(pt_residual,weight);
      for(int c=0;c<nCut;c++){
        if(cut[c]){
          h->h1D(Form("electronPtResidual_cut%d_All",c+1), "", suffix)->Fill(pt_residual,weight);
          h->h1D(Form("electronPtResidual_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(pt_residual,weight);
        }
      }
    }
    numElectron++;
  }
  if(lel_id!=-1){
    TLorentzVector leading_el = evt.electron()[lel_id].mom();
    h->h1D("leadingElectronPt_All", "", suffix)->Fill(leading_el.Perp()*1e-3,weight);
    h->h1D("leadingElectronEta_All", "", suffix)->Fill(leading_el.Eta(),weight);
    h->h1D("leadingElectronPhi_All", "", suffix)->Fill(leading_el.Phi(),weight);
    h->h1D(Form("leadingElectronPt_%s",channel.c_str()), "", suffix)->Fill(leading_el.Perp()*1e-3,weight);
    h->h1D(Form("leadingElectronEta_%s",channel.c_str()), "", suffix)->Fill(leading_el.Eta(),weight);
    h->h1D(Form("leadingElectronPhi_%s",channel.c_str()), "", suffix)->Fill(leading_el.Phi(),weight);
  }
  h->h1D("nElectron_All", "", suffix)->Fill(numElectron,weight);
  h->h1D(Form("nElectron_%s",channel.c_str()), "", suffix)->Fill(numElectron,weight);
  for(int c=0;c<nCut;c++){
    if(cut[c]){
      h->h1D(Form("nElectron_cut%d_All",c+1), "", suffix)->Fill(numElectron,weight);
      h->h1D(Form("nElectron_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(numElectron,weight);
    }
  }

  //no selection electron
  int numNocalibElectron=0;
  //std::cout << "Number of no selection electron is " << evt.nocalibel().size() << std::endl;
  int leadel_id=-1;
  float leadel_pt=0.;
  for(int im=0;im<evt.nocalibel().size();im++){
    TLorentzVector el = evt.nocalibel()[im].mom();
    //std::cout << "pt/eta/phi=" << el.Perp() << "/" << el.Eta() << "/" << el.Phi() << std::endl;
    if(el.Perp()>leadel_pt){
      leadel_id=im;
      leadel_pt=el.Perp();
    }
    float iso  = evt.nocalibel()[im].ptvarcone20() / el.Perp();
    //std::cout << "nocalibrated electron ptvarcone20/iso=" << evt.nocalibel()[im].ptvarcone20() << "/" << iso << std::endl;
    float closejl_deltaR=999.;
    for (unsigned int jet_idx=0; jet_idx < evt.jet().size(); ++jet_idx){    
      float deta = el.Eta() - evt.jet()[jet_idx].mom().Eta(); 
      float dphi = acos(cos(el.Phi() - evt.jet()[jet_idx].mom().Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    //std::cout << "goodOQ/sd0/dz0/iso=" << evt.nocalibel()[im].goodOQ() << "/" << evt.nocalibel()[im].sd0()
      //<< "/" << evt.nocalibel()[im].Dz0() << "/" << iso << std::endl;
    //std::cout << "closejl_deltaR/passLH=" << closejl_deltaR << "/" << evt.nocalibel()[im].passLH() << std::endl;
    if(el.Perp()*1e-3<30) continue;
    h->h1D("nocalibElectronSd0_All", "", suffix)->Fill(evt.nocalibel()[im].sd0(),weight);
    h->h1D("nocalibElectronDz0_All", "", suffix)->Fill(evt.nocalibel()[im].Dz0(),weight);
    h->h1D(Form("nocalibElectronSd0_%s",channel.c_str()), "", suffix)->Fill(evt.nocalibel()[im].sd0(),weight);
    h->h1D(Form("nocalibElectronDz0_%s",channel.c_str()), "", suffix)->Fill(evt.nocalibel()[im].Dz0(),weight);
    if(!evt.nocalibel()[im].goodOQ()) continue;
    if(fabs(evt.nocalibel()[im].sd0())>5) continue;
    if(fabs(evt.nocalibel()[im].Dz0())>0.5) continue;
    if(iso>0.06) continue;
    if(closejl_deltaR<0.4) continue;
    if(!evt.nocalibel()[im].passLH()) continue;
    h->h1D("nocalibElectronPt_All", "", suffix)->Fill(el.Perp()*1e-3,weight);
    h->h1D("nocalibElectronEta_All", "", suffix)->Fill(el.Eta(),weight);
    h->h1D("nocalibElectronPhi_All", "", suffix)->Fill(el.Phi(),weight);
    h->h1D(Form("nocalibElectronPt_%s",channel.c_str()), "", suffix)->Fill(el.Perp()*1e-3,weight);
    h->h1D(Form("nocalibElectronEta_%s",channel.c_str()), "", suffix)->Fill(el.Eta(),weight);
    h->h1D(Form("nocalibElectronPhi_%s",channel.c_str()), "", suffix)->Fill(el.Phi(),weight);
    h->h2D(Form("electronMinJetDrVsNocalibElectronPt_%s",channel.c_str()), "", suffix)->Fill(el.Perp()*1e-3,closejl_deltaR,weight);
    if(iso<0.03){
      h->h1D("nocalibElectronPtTightIso_All", "", suffix)->Fill(el.Perp()*1e-3,weight);
      h->h1D("nocalibElectronEtaTightIso_All", "", suffix)->Fill(el.Eta(),weight);
      h->h1D("nocalibElectronPhiTightIso_All", "", suffix)->Fill(el.Phi(),weight);
      h->h1D(Form("nocalibElectronPtTightIso_%s",channel.c_str()), "", suffix)->Fill(el.Perp()*1e-3,weight);
      h->h1D(Form("nocalibElectronEtaTightIso_%s",channel.c_str()), "", suffix)->Fill(el.Eta(),weight);
      h->h1D(Form("nocalibElectronPhiTightIso_%s",channel.c_str()), "", suffix)->Fill(el.Phi(),weight);
    }
    numNocalibElectron++;
  }
  if(leadel_id!=-1){
    h->h1D("leadingNocalibElectronPt_All", "", suffix)->Fill(evt.nocalibel()[leadel_id].mom().Perp()*1e-3,weight);
    h->h1D(Form("leadingNocalibElectronPt_%s",channel.c_str()), "", suffix)->Fill(evt.nocalibel()[leadel_id].mom().Perp()*1e-3,weight);
  }
  h->h1D("nNocalibElectron_All", "", suffix)->Fill(numNocalibElectron,weight);
  h->h1D(Form("nNocalibElectron_%s",channel.c_str()), "", suffix)->Fill(numNocalibElectron,weight);
  
  //small jet
  int numSmallJet=0,numBtaggedSmallJet=0,numSmallJetHighPt=0;
  int lsj_id=-1;
  float lsj_pt=0.;
  //std::cout << "Number of small jet is " << evt.jet().size() << std::endl;
  for (size_t isj = 0; isj < evt.jet().size(); ++isj){
    TLorentzVector sj = evt.jet()[isj].mom();
    if(sj.Perp()*1e-3<25) continue;
    if(sj.Perp()>lsj_pt){
      lsj_id=isj;
      lsj_pt=sj.Perp();
    }
    //if(fabs(sj.Eta())>2.5) continue;
    h->h1D("smallJetPt_All", "", suffix)->Fill(sj.Perp()*1e-3,weight);
    h->h1D("smallJetEta_All", "", suffix)->Fill(sj.Eta(),weight);
    h->h1D("smallJetPhi_All", "", suffix)->Fill(sj.Phi(),weight);
    h->h1D("smallJetMass_All", "", suffix)->Fill(sj.M()*1e-3,weight);
    h->h1D(Form("smallJetPt_%s",channel.c_str()), "", suffix)->Fill(sj.Perp()*1e-3,weight);
    h->h1D(Form("smallJetEta_%s",channel.c_str()), "", suffix)->Fill(sj.Eta(),weight);
    h->h1D(Form("smallJetPhi_%s",channel.c_str()), "", suffix)->Fill(sj.Phi(),weight);
    h->h1D(Form("smallJetMass_%s",channel.c_str()), "", suffix)->Fill(sj.M()*1e-3,weight);
    if(sj.Perp()*1e-3>100){
      h->h1D("smallJetEtaHighPt_All", "", suffix)->Fill(sj.Eta(),weight);
      h->h1D("smallJetPhiHighPt_All", "", suffix)->Fill(sj.Phi(),weight);
      h->h1D("smallJetMassHighPt_All", "", suffix)->Fill(sj.M()*1e-3,weight);
      h->h1D(Form("smallJetEtaHighPt_%s",channel.c_str()), "", suffix)->Fill(sj.Eta(),weight);
      h->h1D(Form("smallJetPhiHighPt_%s",channel.c_str()), "", suffix)->Fill(sj.Phi(),weight);
      h->h1D(Form("smallJetMassHighPt_%s",channel.c_str()), "", suffix)->Fill(sj.M()*1e-3,weight);
      numSmallJetHighPt++;
    }
    if(sj.Perp()*1e-3<60 && fabs(sj.Eta())<2.4){
      h->h1D("smallJetJVTInRange_All", "", suffix)->Fill(evt.jet()[isj].jvt(),weight);
      h->h1D(Form("smallJetJVTInRange_%s",channel.c_str()), "", suffix)->Fill(evt.jet()[isj].jvt(),weight);
    }
    else {
      h->h1D("smallJetJVTOutRange_All", "", suffix)->Fill(evt.jet()[isj].jvt(),weight);
      h->h1D(Form("smallJetJVTOutRange_%s",channel.c_str()), "", suffix)->Fill(evt.jet()[isj].jvt(),weight);
    }
    if(evt.jet()[isj].btag_mv2c20_70()){
      numBtaggedSmallJet++;
    }
    numSmallJet++;
    for(int c=0;c<nCut;c++){
      if(cut[c]){
        h->h1D(Form("smallJetPt_cut%d_All",c+1), "", suffix)->Fill(sj.Perp()*1e-3,weight);
        h->h1D(Form("smallJetEta_cut%d_All",c+1), "", suffix)->Fill(sj.Eta(),weight);
        h->h1D(Form("smallJetPhi_cut%d_All",c+1), "", suffix)->Fill(sj.Phi(),weight);
        h->h1D(Form("smallJetMass_cut%d_All",c+1), "", suffix)->Fill(sj.M()*1e-3,weight);
        h->h1D(Form("smallJetPt_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(sj.Perp()*1e-3,weight);
        h->h1D(Form("smallJetEta_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(sj.Eta(),weight);
        h->h1D(Form("smallJetPhi_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(sj.Phi(),weight);
        h->h1D(Form("smallJetMass_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(sj.M()*1e-3,weight);
      }
    }
  }
  if(lsj_id!=-1){//leading
    TLorentzVector leading_sj = evt.jet()[lsj_id].mom();
    h->h1D("leadingSmallJetPt_All", "", suffix)->Fill(leading_sj.Perp()*1e-3,weight);
    h->h1D("leadingSmallJetEta_All", "", suffix)->Fill(leading_sj.Eta(),weight);
    h->h1D("leadingSmallJetPhi_All", "", suffix)->Fill(leading_sj.Phi(),weight);
    h->h1D("leadingSmallJetMass_All", "", suffix)->Fill(leading_sj.M()*1e-3,weight);
    h->h1D(Form("leadingSmallJetPt_%s",channel.c_str()), "", suffix)->Fill(leading_sj.Perp()*1e-3,weight);
    h->h1D(Form("leadingSmallJetEta_%s",channel.c_str()), "", suffix)->Fill(leading_sj.Eta(),weight);
    h->h1D(Form("leadingSmallJetPhi_%s",channel.c_str()), "", suffix)->Fill(leading_sj.Phi(),weight);
    h->h1D(Form("leadingSmallJetMass_%s",channel.c_str()), "", suffix)->Fill(leading_sj.M()*1e-3,weight);
  }
  h->h1D("nSmallJet_All", "", suffix)->Fill(numSmallJet,weight);
  h->h1D(Form("nSmallJet_%s",channel.c_str()), "", suffix)->Fill(numSmallJet,weight);
  h->h1D("nSmallJetHighPt_All", "", suffix)->Fill(numSmallJetHighPt,weight);
  h->h1D(Form("nSmallJetHighPt_%s",channel.c_str()), "", suffix)->Fill(numSmallJetHighPt,weight);
  h->h1D("nBtaggedSmallJet_All", "", suffix)->Fill(numBtaggedSmallJet,weight);
  h->h1D(Form("nBtaggedSmallJet_%s",channel.c_str()), "", suffix)->Fill(numBtaggedSmallJet,weight);
  for(int c=0;c<nCut;c++){
    if(cut[c]){
      h->h1D(Form("nSmallJet_cut%d_All",c+1), "", suffix)->Fill(numSmallJet,weight);
      h->h1D(Form("nSmallJet_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(numSmallJet,weight);
    }
  }

  //no selection small jet
  //std::cout << "Number of no selection small jet is " << evt.nocalibjet().size() << std::endl;
  int numNocalibSmallJet=0;
  int leadsj_id=-1;
  float leadsj_pt=0;
  for (size_t isj = 0; isj < evt.nocalibjet().size(); ++isj){
    TLorentzVector sj = evt.nocalibjet()[isj].mom();
    //std::cout << "pt/eta/phi=" << sj.Perp() << "/" << sj.Eta() << "/" << sj.Phi() << std::endl;
    if(sj.Perp()>leadsj_pt){
      leadsj_id=isj;
      leadsj_pt=sj.Perp();
    }
    float closejl_deltaR=999.;
    for (unsigned int el_idx=0; el_idx < evt.nocalibel().size(); ++el_idx){    
      TLorentzVector el = evt.nocalibel()[el_idx].mom();
      float iso  = evt.nocalibel()[el_idx].ptvarcone20() / el.Perp();
      if(el.Perp()*1e-3<30) continue;
      if(!evt.nocalibel()[el_idx].goodOQ()) continue;
      if(fabs(evt.nocalibel()[el_idx].sd0())>5) continue;
      if(fabs(evt.nocalibel()[el_idx].Dz0())>0.5) continue;
      if(iso>0.06) continue;
      if(!evt.nocalibel()[el_idx].passLH()) continue;
      float deta = sj.Eta() - evt.nocalibel()[el_idx].mom().Eta(); 
      float dphi = acos(cos(sj.Phi() - evt.nocalibel()[el_idx].mom().Phi())); 
      float deltaR_tmp = sqrt(deta*deta+dphi*dphi);
      if (deltaR_tmp < closejl_deltaR) closejl_deltaR = deltaR_tmp;
    }//for     
    if(sj.Perp()*1e-3<25) continue;
    if(fabs(sj.Eta())>2.5) continue;
    if(closejl_deltaR<0.2) continue;
    h->h1D("nocalibSmallJetPtNoJVT_All", "", suffix)->Fill(sj.Perp()*1e-3,weight);
    h->h1D(Form("nocalibSmallJetPtNoJVT_%s",channel.c_str()), "", suffix)->Fill(sj.Perp()*1e-3,weight);
    h->h1D("nocalibSmallJetMassNoJVT_All", "", suffix)->Fill(sj.M()*1e-3,weight);
    h->h1D(Form("nocalibSmallJetMassNoJVT_%s",channel.c_str()), "", suffix)->Fill(sj.M()*1e-3,weight);
    for(int c=0;c<nCut;c++){
      if(cut[c]){
        h->h1D(Form("nocalibSmallJetPtNoJVT_cut%d_All",c+1), "", suffix)->Fill(sj.Perp()*1e-3,weight);
        h->h1D(Form("nocalibSmallJetPtNoJVT_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(sj.Perp()*1e-3,weight);
      }
    }
    if(sj.Perp()*1e-3<60 && fabs(sj.Eta())<2.4){
      h->h1D("nocalibSmallJetJVTInRange_All", "", suffix)->Fill(evt.nocalibjet()[isj].jvt(),weight);
      h->h1D(Form("nocalibSmallJetJVTInRange_%s",channel.c_str()), "", suffix)->Fill(evt.nocalibjet()[isj].jvt(),weight);
    }
    else {
      h->h1D("nocalibSmallJetJVTOutRange_All", "", suffix)->Fill(evt.nocalibjet()[isj].jvt(),weight);
      h->h1D(Form("nocalibSmallJetJVTOutRange_%s",channel.c_str()), "", suffix)->Fill(evt.nocalibjet()[isj].jvt(),weight);
    }
    if(sj.Perp()*1e-3<60 && fabs(sj.Eta())<2.4 && evt.nocalibjet()[isj].jvt()<0.59) continue;
    h->h1D("nocalibSmallJetPt_All", "", suffix)->Fill(sj.Perp()*1e-3,weight);
    h->h1D("nocalibSmallJetEta_All", "", suffix)->Fill(sj.Eta(),weight);
    h->h1D("nocalibSmallJetPhi_All", "", suffix)->Fill(sj.Phi(),weight);
    h->h1D("nocalibSmallJetMass_All", "", suffix)->Fill(sj.M()*1e-3,weight);
    h->h1D(Form("nocalibSmallJetPt_%s",channel.c_str()), "", suffix)->Fill(sj.Perp()*1e-3,weight);
    h->h1D(Form("nocalibSmallJetEta_%s",channel.c_str()), "", suffix)->Fill(sj.Eta(),weight);
    h->h1D(Form("nocalibSmallJetPhi_%s",channel.c_str()), "", suffix)->Fill(sj.Phi(),weight);
    h->h1D(Form("nocalibSmallJetMass_%s",channel.c_str()), "", suffix)->Fill(sj.M()*1e-3,weight);
    //std::cout << "nocalibrated jet pt/eta/phi=" << sj.Perp()*1e-3 << "/" << sj.Eta() << "/" << sj.Phi() << std::endl;
    //std::cout << "no selection jet jvt=" << evt.nocalibjet()[isj].jvt() << std::endl;
    h->h1D(Form("nocalibSmallJetMinNoCalibElectronDr_%s",channel.c_str()), "", suffix)->Fill(closejl_deltaR,weight);
    for(int c=0;c<nCut;c++){
      if(cut[c]){
        h->h1D(Form("nocalibSmallJetPt_cut%d_All",c+1), "", suffix)->Fill(sj.Perp()*1e-3,weight);
        h->h1D(Form("nocalibSmallJetPt_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(sj.Perp()*1e-3,weight);
      }
    }
    numNocalibSmallJet++;
  }
  if(leadsj_id!=-1){
    h->h1D("leadingNocalibSmallJetPt_All", "", suffix)->Fill(evt.nocalibjet()[leadsj_id].mom().Perp()*1e-3,weight);
    h->h1D(Form("leadingNocalibSmallJetPt_%s",channel.c_str()), "", suffix)->Fill(evt.nocalibjet()[leadsj_id].mom().Perp()*1e-3,weight);
  }
  h->h1D("nNocalibSmallJet_All", "", suffix)->Fill(numNocalibSmallJet,weight);
  h->h1D(Form("nNocalibSmallJet_%s",channel.c_str()), "", suffix)->Fill(numNocalibSmallJet,weight);

  //no selection large jet
  int numNocalibLargeJet=0;
  //std::cout << "Number of no selection large jet is " << evt.nocalibljet().size() << std::endl;
  int leadlj_id=-1;
  float leadlj_pt=0;
  for (size_t ilj = 0; ilj < evt.nocalibljet().size(); ++ilj){
    TLorentzVector lj = evt.nocalibljet()[ilj].mom();
    if(lj.Perp()>leadlj_pt){
      leadlj_id=ilj;
      leadlj_pt=lj.Perp();
    }
    //std::cout << "pt/eta/phi=" << lj.Perp() << "/" << lj.Eta() << "/" << lj.Phi() << std::endl;
    //if(lj.Perp()*1e-3<100) continue;
    //if(fabs(lj.Eta())>2) continue;
    h->h1D("nocalibLargeJetPt_All", "", suffix)->Fill(lj.Perp()*1e-3,weight);
    h->h1D("nocalibLargeJetEta_All", "", suffix)->Fill(lj.Eta(),weight);
    h->h1D("nocalibLargeJetPhi_All", "", suffix)->Fill(lj.Phi(),weight);
    h->h1D("nocalibLargeJetMass_All", "", suffix)->Fill(lj.M()*1e-3,weight);
    h->h1D(Form("nocalibLargeJetPt_%s",channel.c_str()), "", suffix)->Fill(lj.Perp()*1e-3,weight);
    h->h1D(Form("nocalibLargeJetEta_%s",channel.c_str()), "", suffix)->Fill(lj.Eta(),weight);
    h->h1D(Form("nocalibLargeJetPhi_%s",channel.c_str()), "", suffix)->Fill(lj.Phi(),weight);
    h->h1D(Form("nocalibLargeJetMass_%s",channel.c_str()), "", suffix)->Fill(lj.M()*1e-3,weight);
    if(lj.Perp()*1e-3>300) numNocalibLargeJet++;
  }
  if(leadlj_id!=-1){
    h->h1D("leadingNocalibLargeJetPt_All", "", suffix)->Fill(evt.nocalibljet()[leadlj_id].mom().Perp()*1e-3,weight);
    h->h1D(Form("leadingNocalibLargeJetPt_%s",channel.c_str()), "", suffix)->Fill(evt.nocalibljet()[leadlj_id].mom().Perp()*1e-3,weight);
  }
  //std::cout << "numLargeJet nocalib=" << numNocalibLargeJet << std::endl;
  h->h1D("nNocalibLargeJet_All", "", suffix)->Fill(numNocalibLargeJet,weight);
  h->h1D(Form("nNocalibLargeJet_%s",channel.c_str()), "", suffix)->Fill(numNocalibLargeJet,weight);


  //Track jet
  int nTrackJets=0,nBtaggedTrackJets=0;
  //std::cout << "Number of track jets is " << evt.tjet().size() << std::endl;
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    TLorentzVector tj = evt.tjet()[jidx].mom();
    if(tj.Perp()*1e-3<10) continue;
    if(fabs(tj.Eta())>2.4) continue;
    h->h1D("trackJetPt_All", "", suffix)->Fill(tj.Perp()*1e-3,weight);
    h->h1D("trackJetEta_All", "", suffix)->Fill(tj.Eta(),weight);
    h->h1D("trackJetPhi_All", "", suffix)->Fill(tj.Phi(),weight);
    h->h1D("trackJetNumConstituents_All", "", suffix)->Fill(evt.tjet()[jidx].numConstituents(),weight);
    h->h1D("trackJetMv2c10_All", "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
    h->h2D("trackJetPtVsMv2c10_All", "", suffix)->Fill(tj.Perp()*1e-3,evt.tjet()[jidx].mv2c10(),weight);
    h->h1D(Form("trackJetPt_%s",channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,weight);
    h->h1D(Form("trackJetEta_%s",channel.c_str()), "", suffix)->Fill(tj.Eta(),weight);
    h->h1D(Form("trackJetPhi_%s",channel.c_str()), "", suffix)->Fill(tj.Phi(),weight);
    h->h1D(Form("trackJetNumConstituents_%s",channel.c_str()), "", suffix)->Fill(evt.tjet()[jidx].numConstituents(),weight);
    h->h1D(Form("trackJetMv2c10_%s",channel.c_str()), "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
    h->h2D(Form("trackJetPtVsMv2c10_%s",channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,evt.tjet()[jidx].mv2c10(),weight);
    nTrackJets++;
    for(int c=0;c<nCut;c++){
      if(cut[c]){
        h->h1D(Form("trackJetPt_cut%d_All",c+1), "", suffix)->Fill(tj.Perp()*1e-3,weight);
        h->h1D(Form("trackJetEta_cut%d_All",c+1), "", suffix)->Fill(tj.Eta(),weight);
        h->h1D(Form("trackJetPhi_cut%d_All",c+1), "", suffix)->Fill(tj.Phi(),weight);
        h->h1D(Form("trackJetMv2c10_cut%d_All",c+1), "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
        h->h1D(Form("trackJetPt_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,weight);
        h->h1D(Form("trackJetEta_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(tj.Eta(),weight);
        h->h1D(Form("trackJetPhi_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(tj.Phi(),weight);
        h->h1D(Form("trackJetMv2c10_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
      }
    }
    if(evt.tjet()[jidx].btag_mv2c10_70_trk()){
      h->h1D("btaggedTrackJetPt_All", "", suffix)->Fill(tj.Perp()*1e-3,weight);
      h->h1D("btaggedTrackJetEta_All", "", suffix)->Fill(tj.Eta(),weight);
      h->h1D("btaggedTrackJetPhi_All", "", suffix)->Fill(tj.Phi(),weight);
      h->h1D(Form("btaggedTrackJetPt_%s",channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,weight);
      h->h1D(Form("btaggedTrackJetEta_%s",channel.c_str()), "", suffix)->Fill(tj.Eta(),weight);
      h->h1D(Form("btaggedTrackJetPhi_%s",channel.c_str()), "", suffix)->Fill(tj.Phi(),weight);
      nBtaggedTrackJets++;
    }
    if(jidx==bt_tj_idx || jidx==btbar_tj_idx){//true b
      h->h1D("trackBjetPt_All", "", suffix)->Fill(tj.Perp()*1e-3,weight);
      h->h1D("trackBjetEta_All", "", suffix)->Fill(tj.Eta(),weight);
      h->h1D("trackBjetPhi_All", "", suffix)->Fill(tj.Phi(),weight);
      h->h1D("trackBjetNumConstituents_All", "", suffix)->Fill(evt.tjet()[jidx].numConstituents(),weight);
      h->h1D("trackBjetMv2c10_All", "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
      h->h1D(Form("trackBjetPt_%s",channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,weight);
      h->h1D(Form("trackBjetEta_%s",channel.c_str()), "", suffix)->Fill(tj.Eta(),weight);
      h->h1D(Form("trackBjetPhi_%s",channel.c_str()), "", suffix)->Fill(tj.Phi(),weight);
      h->h1D(Form("trackBjetNumConstituents_%s",channel.c_str()), "", suffix)->Fill(evt.tjet()[jidx].numConstituents(),weight);
      h->h1D(Form("trackBjetMv2c10_%s",channel.c_str()), "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
      if(evt.tjet()[jidx].btag_mv2c10_70_trk()){
        h->h1D("btaggedTrackBjetPt_All", "", suffix)->Fill(tj.Perp()*1e-3,weight);
        h->h1D("btaggedTrackBjetEta_All", "", suffix)->Fill(tj.Eta(),weight);
        h->h1D("btaggedTrackBjetPhi_All", "", suffix)->Fill(tj.Phi(),weight);
        h->h1D(Form("btaggedTrackBjetPt_%s",channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,weight);
        h->h1D(Form("btaggedTrackBjetEta_%s",channel.c_str()), "", suffix)->Fill(tj.Eta(),weight);
        h->h1D(Form("btaggedTrackBjetPhi_%s",channel.c_str()), "", suffix)->Fill(tj.Phi(),weight);
      }
    }
    else {
      h->h1D("trackLightJetMv2c10_All", "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
      h->h1D(Form("trackLightJetMv2c10_%s",channel.c_str()), "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
    }
  }
  h->h1D("nTrackJet_All", "", suffix)->Fill(nTrackJets,weight);
  h->h1D(Form("nTrackJet_%s",channel.c_str()), "", suffix)->Fill(nTrackJets,weight);
  h->h1D("nBtaggedTrackJet_All", "", suffix)->Fill(nBtaggedTrackJets,weight);
  h->h1D(Form("nBtaggedTrackJet_%s",channel.c_str()), "", suffix)->Fill(nBtaggedTrackJets,weight);
  for(int c=0;c<nCut;c++){
    if(cut[c]){
      h->h1D(Form("nTrackJet_cut%d_All",c+1), "", suffix)->Fill(nTrackJets,weight);
      h->h1D(Form("nTrackJet_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(nTrackJets,weight);
    }
  }

  //No selection track jet
  int numTrackJets=0,numBtaggedTrackJets=0;
  //std::cout << "Number of no selectoin tjet is " << evt.nocalibtjet().size() << std::endl;
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
    if(closejmu_deltaR_def>0.2){
      h->h1D("nocalibTrackJetPt_def_All", "", suffix)->Fill(tj.Perp()*1e-3,weight);
      h->h1D("nocalibTrackJetEta_def_All", "", suffix)->Fill(tj.Eta(),weight);
      h->h1D("nocalibTrackJetPhi_def_All", "", suffix)->Fill(tj.Phi(),weight);
      h->h1D(Form("nocalibTrackJetPt_def_%s",channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,weight);
      h->h1D(Form("nocalibTrackJetEta_def_%s",channel.c_str()), "", suffix)->Fill(tj.Eta(),weight);
      h->h1D(Form("nocalibTrackJetPhi_def_%s",channel.c_str()), "", suffix)->Fill(tj.Phi(),weight);
    }
    if(closejmu_deltaR<0.2) continue;
    h->h1D("nocalibTrackJetPt_All", "", suffix)->Fill(tj.Perp()*1e-3,weight);
    h->h1D("nocalibTrackJetEta_All", "", suffix)->Fill(tj.Eta(),weight);
    h->h1D("nocalibTrackJetPhi_All", "", suffix)->Fill(tj.Phi(),weight);
    h->h1D("nocalibTrackJetNumConstituents_All", "", suffix)->Fill(evt.tjet()[jidx].numConstituents(),weight);
    h->h1D("nocalibTrackJetMv2c10_All", "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
    h->h1D(Form("nocalibTrackJetPt_%s",channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,weight);
    h->h1D(Form("nocalibTrackJetEta_%s",channel.c_str()), "", suffix)->Fill(tj.Eta(),weight);
    h->h1D(Form("nocalibTrackJetPhi_%s",channel.c_str()), "", suffix)->Fill(tj.Phi(),weight);
    h->h1D(Form("nocalibTrackJetNumConstituents_%s",channel.c_str()), "", suffix)->Fill(evt.tjet()[jidx].numConstituents(),weight);
    h->h1D(Form("nocalibTrackJetMv2c10_%s",channel.c_str()), "", suffix)->Fill(evt.tjet()[jidx].mv2c10(),weight);
    numTrackJets++;
    if(evt.nocalibtjet()[jidx].btag_mv2c10_70_trk()){
      h->h1D("btaggedNocalibTrackJetPt_All", "", suffix)->Fill(tj.Perp()*1e-3,weight);
      h->h1D("btaggedNocalibTrackJetEta_All", "", suffix)->Fill(tj.Eta(),weight);
      h->h1D("btaggedNocalibTrackJetPhi_All", "", suffix)->Fill(tj.Phi(),weight);
      h->h1D(Form("btaggedNocalibTrackJetPt_%s",channel.c_str()), "", suffix)->Fill(tj.Perp()*1e-3,weight);
      h->h1D(Form("btaggedNocalibTrackJetEta_%s",channel.c_str()), "", suffix)->Fill(tj.Eta(),weight);
      h->h1D(Form("btaggedNocalibTrackJetPhi_%s",channel.c_str()), "", suffix)->Fill(tj.Phi(),weight);
      numBtaggedTrackJets++;
    }
  }
  h->h1D("nNocalibTrackJet_All", "", suffix)->Fill(numTrackJets,weight);
  h->h1D(Form("nNocalibTrackJet_%s",channel.c_str()), "", suffix)->Fill(numTrackJets,weight);
  h->h1D("nBtaggedNocalibTrackJet_All", "", suffix)->Fill(numBtaggedTrackJets,weight);
  h->h1D(Form("nBtaggedNocalibTrackJet_%s",channel.c_str()), "", suffix)->Fill(numBtaggedTrackJets,weight);

  //MET
  h->h1D("MET_All", "", suffix)->Fill(evt.met().Perp()*1e-3,weight);
  h->h1D("METPhi_All", "", suffix)->Fill(evt.met().Phi(),weight);
  h->h1D(Form("MET_%s",channel.c_str()), "", suffix)->Fill(evt.met().Perp()*1e-3,weight);
  h->h1D(Form("METPhi_%s",channel.c_str()), "", suffix)->Fill(evt.met().Phi(),weight);
  for(int c=0;c<nCut;c++){
    if(cut[c]){
      h->h1D(Form("MET_cut%d_All",c+1), "", suffix)->Fill(evt.met().Perp()*1e-3,weight);
      h->h1D(Form("METPhi_cut%d_All",c+1), "", suffix)->Fill(evt.met().Phi(),weight);
      h->h1D(Form("MET_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(evt.met().Perp()*1e-3,weight);
      h->h1D(Form("METPhi_cut%d_%s",c+1,channel.c_str()), "", suffix)->Fill(evt.met().Phi(),weight);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  //Event selection

  //Trigger selection
  if(!passtrig_mu1 && !passtrig_mu2) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(1);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(1);
  //std::cout << "pass trigger!!" << std::endl; 


  //Number of muon > 0
  //if(numNocalibMuon==0) return;//temporal
  if(numMuon==0) return;//temporal
  h->h1D("cutFlow_All", "", suffix)->Fill(2);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(2);
  //std::cout << "nMuon > 0" << std::endl;

  //Number of muon == 1
  //if(numNocalibMuon>1) return;//temporal
  if(numMuon25>1) return;//temporal
  h->h1D("cutFlow_All", "", suffix)->Fill(3);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(3);
  //std::cout << "nMuon==1" << std::endl;

  //Number of electron == 0
  int nElectron = 0,electron_id=-1;
  for(int ie=0;ie<evt.electron().size();ie++){
    //std::cout << "electron " << ie << std::endl;
    TLorentzVector el = evt.electron()[ie].mom();
    if(el.Perp()*1e-3 > 25) {
      if(nElectron==0) electron_id=ie;
      nElectron++;
    }
    bool isTight = evt.electron()[ie].isTightPP();
    //std::cout << "pt/isTight=" << "/" << el.Perp()*1e-3 << "/" << isTight << std::endl;
  }
  if(nElectron!=0) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(4);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(4);
  //std::cout << "nElectron==0" << std::endl;

  //Trigger matching
  bool trig1 = evt.muon()[lmu_id].HLT_mu26_ivarmedium();
  bool trig2 = evt.muon()[lmu_id].HLT_mu50();
  //if(!trig1 && !trig2) return;
  //if(!evt.nocalibmu()[muon_id_nocalib].trigger_match()) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(5);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(5);
  //std::cout << "Pass Trigger Matching" << std::endl;

  //MET & MET+MWT 
  TLorentzVector l;
  l.SetPtEtaPhiM(evt.muon()[lmu_id].mom().Perp(),evt.muon()[lmu_id].mom().Eta(),evt.muon()[lmu_id].mom().Phi(),105.6);
  //l.SetPtEtaPhiM(evt.nocalibmu()[muon_id_nocalib].mom().Perp(),evt.nocalibmu()[muon_id_nocalib].mom().Eta(),evt.nocalibmu()[muon_id_nocalib].mom().Phi(),105.6);//temporal
  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3;
  float MET = evt.met().Perp()*1e-3;
  float MET_plus_mWt = mWt + MET;
  //std::cout << "MET/MWT=" << MET << "/" << mWt << std::endl;
  if(MET<20) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(6);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(6);
  //std::cout << "Pass MET" << std::endl;
  if(MET_plus_mWt<60) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(7);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(7);
  //std::cout << "Pass MET+mWt" << std::endl;
  
  //Number of small jets
  //std::cout << "Number of small jet is " << evt.jet().size() << std::endl;
  int nJet=0;
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    float jet_pt = evt.jet()[jidx].mom().Perp()*1e-3;
    float jet_eta = evt.jet()[jidx].mom().Eta();
    float jet_phi = evt.jet()[jidx].mom().Phi();
    if(jet_pt>25) nJet++;
  }
  if(nJet==0) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(8);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(8);
  //std::cout << "Pass #Small Jet" << std::endl;

  //Jet close to lepton
  size_t close_idx = 0;
  for (; close_idx < evt.jet().size(); ++close_idx){    
    float dphi = l.DeltaPhi(evt.jet()[close_idx].mom());
    float dy = l.Rapidity() - evt.jet()[close_idx].mom().Rapidity();
    float deltaR_tmp = std::sqrt(dy*dy + dphi*dphi);
    if(evt.jet()[close_idx].mom().Pt()*1e-3 > 25 && deltaR_tmp < 1.5) break;
  }//for     
  //std::cout << "close_idx/nJet=" << close_idx << "/" << evt.jet().size() << std::endl;
  if(close_idx==evt.jet().size()) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(9);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(9);
  //std::cout << "Pass Jet Close to Lepton" << std::endl;

  float pt_range[] = { 250000.000,325000.000,375000.000,425000.000,475000.000,525000.000,575000.000,625000.000,675000.000,725000.000,775000.000,850000.000,950000.000,1100000.000,1300000.000,1680000.000,1e10} ;
  int num_pt_range = sizeof(pt_range)/sizeof(float);
  float mass_thr[] = {67888.967,72014.026,74764.066,76769.667,78354.344,79170.000,79530.000,80158.525,81195.851,82779.245,84890.965,88747.162,94262.629,102710.787,113868.253,135067.438  };
  int num_mass_thr = sizeof(mass_thr)/sizeof(float);
  float tau32_thr[] = { 0.879,0.831,0.799,0.770,0.746,0.727,0.714,0.706,0.701,0.698,0.698,0.699,0.700,0.701,0.699,0.696  };
  int num_tau32_thr = sizeof(tau32_thr)/sizeof(float);
  //std::cout << "num pt_range/mass_thr/tau32_thr=" << num_pt_range << "/" << num_mass_thr << "/" << num_tau32_thr << std::endl;
  //Number of large jet
  int nLargeJet=0,nLargeJetPt=0;
  int ljetid=-1;
  int closest_ljetid=-1;
  float mindr_ljet_lepton=999.;
  //std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
  for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
    const TLorentzVector &lj = evt.largeJet()[ljid].mom();
    int good = evt.largeJet()[ljid].good();//default
    //std::cout << "pt/eta/phi=" << lj.Perp() << "/" << lj.Eta() << "/" << lj.Phi() << std::endl;
    //std::cout << "mass/tau32=" << lj.M() << "/" << evt.largeJet()[ljid].subs("tau32_wta") << std::endl;
    bool passMassCut=false,passTau32Cut=false;
    for(int t=0;t<num_pt_range-1;t++){
      if(lj.Perp()>pt_range[t] && lj.Perp()<pt_range[t+1]){
        //std::cout << pt_range[t] << " < pt < " << pt_range[t+1] << std::endl;
        //std::cout << "thr mass/tau32=" << mass_thr[t] << "/" << tau32_thr[t] << std::endl;
        if(lj.M()>mass_thr[t]) passMassCut=true;
        //if(evt.largeJet()[ljid].subs("tau32_wta")<tau32_thr[t]) passTau32Cut=true;
      }
    }
    bool good_tmp = passMassCut && passTau32Cut;
    //std::cout << "good def/tmp=" << good << "/" << good_tmp << std::endl;
    good=good_tmp;//temporal
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0) {
      h->h1D(Form("largeJetPt_ana_%s",channel.c_str()), "", suffix)->Fill(lj.Perp()*1e-3,weight);
      h->h1D(Form("largeJetEta_ana_%s",channel.c_str()), "", suffix)->Fill(lj.Eta(),weight);
      h->h1D(Form("largeJetGood_ana_%s",channel.c_str()), "", suffix)->Fill(good,weight);
      h->h1D(Form("largeJetMass_ana_%s",channel.c_str()), "", suffix)->Fill(lj.M()*1e-3,weight);
      //h->h1D(Form("largeJetTau32_ana_%s",channel.c_str()), "", suffix)->Fill(evt.largeJet()[ljid].subs("tau32_wta"),weight);
      nLargeJetPt++; 
    }
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0 && good) {
      nLargeJet++;
      if(ljetid==-1) ljetid=ljid;
      h->h1D(Form("goodLargeJetPt_ana_%s",channel.c_str()), "", suffix)->Fill(lj.Perp()*1e-3,weight);
      h->h1D(Form("goodLargeJetMass_ana_%s",channel.c_str()), "", suffix)->Fill(lj.M()*1e-3,weight);
      //h->h1D(Form("goodLargeJetTau32_ana_%s",channel.c_str()), "", suffix)->Fill(evt.largeJet()[ljid].subs("tau32_wta"),weight);
      h->h1D(Form("goodLargeJetEta_ana_%s",channel.c_str()), "", suffix)->Fill(lj.Eta(),weight);
      //std::cout << "pt=" << lj.Perp() << std::endl;
      for(int p=0;p<num_pt_range-1;p++){
        //std::cout << pt_range[p] << " < pt < " << pt_range[p+1] << std::endl; 
        if(lj.Perp()>pt_range[p] && lj.Perp()<pt_range[p+1]){
          //std::cout << "In the range!!" << std::endl;
          h->h1D(Form("goodLargeJetMass_ana_%s_ptrange%d",channel.c_str(),p), "", suffix)->Fill(lj.M()*1e-3,weight);
          h->h1D(Form("goodLargeJetTau32_ana_%s_ptrange%d",channel.c_str(),p), "", suffix)->Fill(evt.largeJet()[ljid].subs("tau32_wta"),weight);
        }
      }
    }
  }
  h->h1D(Form("nLargeJet_ana_%s",channel.c_str()), "", suffix)->Fill(nLargeJetPt,weight);
  h->h1D(Form("nGoodLargeJet_ana_%s",channel.c_str()), "", suffix)->Fill(nLargeJet,weight);
  if(nLargeJetPt==0) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(10);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(10);
  //std::cout << "Pass #Large Jet" << std::endl;

  //Angular cut
  int nGoodLargeJets=0;
  int goodljetid=-1;
  TLorentzVector sj = evt.jet()[close_idx].mom();
  for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
    const TLorentzVector &lj = evt.largeJet()[ljid].mom();
    bool passMassCut=false,passTau32Cut=false;
    for(int t=0;t<num_pt_range-1;t++){
      if(lj.Perp()>pt_range[t] && lj.Perp()<pt_range[t+1]){
        //std::cout << pt_range[t] << " < pt < " << pt_range[t+1] << std::endl;
        //std::cout << "thr mass/tau32=" << mass_thr[t] << "/" << tau32_thr[t] << std::endl;
        if(lj.M()>mass_thr[t]) passMassCut=true;
        //if(evt.largeJet()[ljid].subs("tau32_wta")<tau32_thr[t]) passTau32Cut=true;
      }
    }
    bool good_tmp = passMassCut && passTau32Cut;
    //int good = evt.largeJet()[ljid].good();//default
    int good=good_tmp;//temporal
    float delta_phi_large_jet_lepton = lj.DeltaPhi(l);
    float delta_r_large_jet_small_jet = lj.DeltaR(sj);
    //std::cout << "large jet pt/eta/good=" << lj.Perp() << "/" << lj.Eta() << "/" << good << std::endl;
    //std::cout << "dr/dphi=" << delta_r_large_jet_small_jet << "/" << delta_phi_large_jet_lepton << std::endl;
    //if(lj.Perp()*1e-3 < 300 || fabs(lj.Eta())>2.0 || !good) continue;
    if(lj.Perp()*1e-3 < 300 || fabs(lj.Eta())>2.0) continue;//temporal
    if(fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) {
      if(goodljetid==-1) goodljetid = ljid;
      nGoodLargeJets++;
    }
  }
  //std::cout << "nGoodLargeJets=" << nGoodLargeJets << std::endl;
  if(nGoodLargeJets==0) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(11);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(11);
  //std::cout << "Pass Angular Cut" << std::endl;

  //Track B-tag
  int nTrackBtagJets=0;
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    float pt = evt.tjet()[jidx].mom().Pt()*1e-3;
    float eta = evt.tjet()[jidx].mom().Eta();
    float phi = evt.tjet()[jidx].mom().Phi();
    //std::cout << "tjet pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl; 
    //std::cout << "b-tag value = " << evt.tjet()[jidx].mv2c10() << std::endl;
    //std::cout << "numConstituents=" << evt.tjet()[jidx].numConstituents() << std::endl;
    bool is_trackBtag = evt.tjet()[jidx].btag_mv2c10_70_trk();
    if (pt > 10 && is_trackBtag) nTrackBtagJets++; 
  }
  //std::cout << "Number of track b-tagged jets is " << nTrackBtagJets << std::endl;

  if(nTrackBtagJets==0) return;
  //if(numBtaggedTrackJets==0) return;
  h->h1D("cutFlow_All", "", suffix)->Fill(12);
  h->h1D(Form("cutFlow_%s",channel.c_str()), "", suffix)->Fill(12);
  //std::cout << "Pass trk B-tag" << std::endl;
  //h->h1D(Form("leadingLargeJetTau32_ana_%s",channel.c_str()), "", suffix)->Fill(evt.largeJet()[goodljetid].subs("tau32_wta"),weight);

}
