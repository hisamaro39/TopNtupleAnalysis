/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
//#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnalysisValidation.h"
#include "TopNtupleAnalysis/AnaTtresMultiJets.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"
#include "TopNtupleAnalysis/UserFunktions.h"


AnaTtresMultiJets::AnaTtresMultiJets(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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

    m_hSvc.create1D("decayType", "; type; Events", 5,0,5);
    m_hSvc.create1D("Wdecay1_from_t_pdgId", "; pdgId; Events", 40,-20,20);
    m_hSvc.create1D("mu", "; #mu; Events", 50,0,50);

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
    m_hSvc.create1D("allSmallJetPt", "; small jet p_{T} [GeV]; Events", 50, 0, 1000);
    m_hSvc.create1D("allSmallJetEta", "; small jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("allSmallJetPhi", "; small jet #phi ; Events", 32, -3.2, 3.2);

    m_hSvc.create1D("largeJetPt", "; large jet p_{T} [GeV] ; Events", 50, 0, 2000);
    m_hSvc.create1D("largeJetEta", "; large jet #eta ; Events", 20, -2., 2.);
    m_hSvc.create1D("largeJetPhi", "; large jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("largeJetMass", "; large jet mass [GeV] ; Events", 30, 0, 300);
    m_hSvc.create1D("largeJetTau32", "; large jet tau32  ; Events", 20, 0, 1);
    m_hSvc.create1D("nLargeJet", "; # of large jets ; Events", 10, 0, 10);
    m_hSvc.create1D("topMatchedLargeJetPt", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("topMatchedLargeJetPtPassSub80", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("topMatchedLargeJetPtPassSub50", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("topMatchedLargeJetPtPassSmoothMT80", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("topMatchedLargeJetPtPassSmoothMT50", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("wMatchedLargeJetPt", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("qcdMatchedLargeJetPt", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("wMatchedLargeJetPtPassSub80", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("qcdMatchedLargeJetPtPassSub80", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("wMatchedLargeJetPtPassSub50", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("qcdMatchedLargeJetPtPassSub50", "; large jet p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("topMatchedLargeJetParticlePt", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("wMatchedLargeJetParticlePt", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("qcdMatchedLargeJetParticlePt", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("topMatchedLargeJetParticlePtPassSub80", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("topMatchedLargeJetParticlePtPassSmoothMT80", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("wMatchedLargeJetParticlePtPassSub80", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("qcdMatchedLargeJetParticlePtPassSub80", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("topMatchedLargeJetParticlePtPassSub50", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("wMatchedLargeJetParticlePtPassSub50", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("qcdMatchedLargeJetParticlePtPassSub50", "; particle p_{T} [GeV] ; Events", 11, 300, 2500);
    m_hSvc.create1D("htopMatchedLargeJetParticlePt", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassSub80", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassSub50", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassSmoothMT80", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassSmoothMT50", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassSmoothTS80", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassSmoothTS50", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassSmoothQT80", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassSmoothQT50", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassBDT80", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetParticlePtPassDNN80", "; particle p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPt", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassSub80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassSub50", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassSmoothMT80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassSmoothMT50", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassSmoothTS80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassSmoothTS50", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassSmoothQT80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassSmoothQT50", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassBDT80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("htopMatchedLargeJetPtPassDNN80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create2D("htopMatchedLargeJetPtVsMass", "; large jet p_{T} [GeV]; large jet mass [GeV] ; Events", 100, 0, 2000,100,0,300);
    m_hSvc.create2D("htopMatchedLargeJetPtVsTau32", "; large jet p_{T} [GeV]; large jet tau32 ; Events", 100, 0, 2000,100,0,1);
    m_hSvc.create1D("allLargeJetPt", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassSub80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassSub50", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassSmoothMT80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassSmoothMT50", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassSmoothTS80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassSmoothTS50", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassSmoothQT80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassSmoothQT50", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassBDT80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    m_hSvc.create1D("allLargeJetPtPassDNN80", "; large jet p_{T} [GeV] ; Events", 12, 350, 2750);
    for(int p=0;p<20;p++){
      m_hSvc.create1D(Form("htopMatchedLargeJetMass_pt%d_%d",100*p,100*(p+1)), "; large jet mass [GeV] ; Events", 50, 0, 300);
      m_hSvc.create1D(Form("htopMatchedLargeJetTau32_pt%d_%d",100*p,100*(p+1)), "; large jet tau32  ; Events", 50, 0, 1);
    }

    m_hSvc.create2D("allLargeJetPtVsMassPassSmoothMT80", "; large jet p_{T} [GeV]; large jet mass [GeV] ; Events", 100, 0, 2000,100,0,300);
    m_hSvc.create2D("allLargeJetPtVsTau32PassSmoothMT80", "; large jet p_{T} [GeV]; large jet tau32 ; Events", 100, 0, 2000,100,0,1);
    m_hSvc.create2D("allLargeJetPtVsMassPassSub80", "; large jet p_{T} [GeV]; large jet mass [GeV] ; Events", 100, 0, 2000,100,0,300);
    m_hSvc.create2D("allLargeJetPtVsTau32PassSub80", "; large jet p_{T} [GeV]; large jet tau32 ; Events", 100, 0, 2000,100,0,1);

    m_hSvc.create1D("trackJetPt", "; track jet p_{T} [GeV]; Events", 25, 0, 500);
    m_hSvc.create1D("trackJetEta", "; track jet #eta ; Events", 25,-2.5,2.5);
    m_hSvc.create1D("trackJetPhi", "; track jet #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("trackJetMv2c10", "; track jet Mv2c10 ; Events", 20, -1, 1);
    m_hSvc.create1D("nTrackJet", "; # of track jets ; Events", 15, 0, 15);  
    m_hSvc.create1D("nBtaggedTrackJet", "; # of track jets ; Events", 10, 0, 10);  

    m_hSvc.create1D("MET", "; missing E_[T] [GeV]; Events", 25, 0, 1000);
    m_hSvc.create1D("metPhi", "; missing E_{T} #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("mwt", "; M_{T,W} [GeV]; Events", 25, 0, 250);
    m_hSvc.create1D("wpt", "; p_{T,W} [GeV]; Events", 100, 0, 1000);
    m_hSvc.create1D("wptPassMet110", "; p_{T,W} [GeV]; Events", 100, 0, 1000);
    m_hSvc.create1D("wptPassPufitMet110", "; p_{T,W} [GeV]; Events", 100, 0, 1000);

    m_hSvc.create1D("deltaR_lepSmallJet", "; #DeltaR_{lep,sjet}; Events", 20, 0, 2);
    m_hSvc.create1D("deltaR_largeJetSmallJet", "; #DeltaR_{ljet,sjet}; Events", 50, 1, 6);
    m_hSvc.create1D("deltaPhi_largeJetLepton", "; #Delta#phi_{ljet,lep}; Events", 24, 2, 3.2);

    m_hSvc.create1DVar("lbMatchedTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedTrackJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedTrackJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("otherMatchedTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedTrackJetPtPassTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedTrackJetPtPassTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedTrackJetParticlePtPassTrackBtag", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedTrackJetParticlePtPassTrackBtag", "; truth jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("otherMatchedTrackJetPtPassTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("lbMatchedTrackJetMv2c10", "; track jet Mv2c10; Events", 20,-1,1);
    m_hSvc.create1D("hbMatchedTrackJetMv2c10", "; track jet Mv2c10; Events", 20,-1,1);
    m_hSvc.create1D("otherMatchedTrackJetMv2c10", "; track jet Mv2c10; Events", 20,-1,1);
    m_hSvc.create1DVar("lbMatchedSmallJetPt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedSmallJetPt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedSmallJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedSmallJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedClosestSmallJetPt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedClosestSmallJetPt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedClosestSmallJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedClosestSmallJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("otherMatchedSmallJetPt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedSmallJetPtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedSmallJetPtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedSmallJetParticlePtPassCaloBtagFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedSmallJetParticlePtPassCaloBtagFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedClosestSmallJetPtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedClosestSmallJetPtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedClosestSmallJetParticlePtPassCaloBtagFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedClosestSmallJetParticlePtPassCaloBtagFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lbMatchedSmallJetPtPassCaloBtagHyb70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("hbMatchedSmallJetPtPassCaloBtagHyb70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("otherMatchedSmallJetPtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("otherMatchedSmallJetPtPassCaloBtagHyb70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("lbMatchedSmallJetMv2c10", "; small jet Mv2c10; Events", 20,-1,1);
    m_hSvc.create1D("hbMatchedSmallJetMv2c10", "; small jet Mv2c10; Events", 20,-1,1);
    m_hSvc.create1D("otherMatchedSmallJetMv2c10", "; small jet Mv2c10; Events", 20,-1,1);

    m_hSvc.create1DVar("bMatchedTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedTrackJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedTrackJetPtPassTrackBtag", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedTrackJetParticlePtPassTrackBtag", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("bMatchedTrackJetMv2c10", "; track jet Mv2c10; Events", 100,-1,1);
    m_hSvc.create1DVar("cMatchedTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedTrackJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedTrackJetPtPassTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedTrackJetParticlePtPassTrackBtag", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("cMatchedTrackJetMv2c10", "; track jet Mv2c10; Events", 100,-1,1);
    m_hSvc.create1DVar("lightMatchedTrackJetPt", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedTrackJetParticlePt", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedTrackJetPtPassTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedTrackJetParticlePtPassTrackBtag", "; track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("lightMatchedTrackJetMv2c10", "; track jet Mv2c10; Events", 100,-1,1);
    int ptRange[] = {0,100,200,300,500,1000,3000};
    int nPtRange = sizeof(ptRange)/sizeof(int);
    std::string btagType[] = {"Mv2c10","Mv2c10mu","Mv2c10rnn"};
    int nBtagType = sizeof(btagType)/sizeof(std::string);
    std::string dl1pType[] = {"Dl1Pu","Dl1Pb","Dl1Pc","Dl1muPu","Dl1muPb","Dl1muPc","Dl1rnnPu","Dl1rnnPb","Dl1rnnPc"};
    int nDl1pType = sizeof(dl1pType)/sizeof(std::string);
    std::string dl1wType[] = {"Dl1Weight","Dl1muWeight","Dl1rnnWeight"};
    int nDl1wType = sizeof(dl1wType)/sizeof(std::string);
    for(int p=0;p<nPtRange-1;p++){
      for(int b=0;b<nBtagType;b++){
        m_hSvc.create1D(Form("bMatchedSmallJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("cMatchedSmallJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("lightMatchedSmallJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("bMatchedTrackJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("cMatchedTrackJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("lightMatchedTrackJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("bMatchedVRTrackJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; vr track jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("cMatchedVRTrackJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; vr track jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("lightMatchedVRTrackJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; vr track jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("bMatchedTrackJetInSmallJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("cMatchedTrackJetInSmallJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",btagType[b].c_str()), 100,-1,1);
        m_hSvc.create1D(Form("lightMatchedTrackJetInSmallJet%s_pt%d_%d",btagType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",btagType[b].c_str()), 100,-1,1);
      }
      for(int b=0;b<nDl1pType;b++){
        m_hSvc.create1D(Form("bMatchedSmallJet%s_pt%d_%d",dl1pType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",dl1pType[b].c_str()), 100,0,1);
        m_hSvc.create1D(Form("cMatchedSmallJet%s_pt%d_%d",dl1pType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",dl1pType[b].c_str()), 100,0,1);
        m_hSvc.create1D(Form("lightMatchedSmallJet%s_pt%d_%d",dl1pType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",dl1pType[b].c_str()), 100,0,1);
        m_hSvc.create1D(Form("bMatchedTrackJet%s_pt%d_%d",dl1pType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",dl1pType[b].c_str()), 100,0,1);
        m_hSvc.create1D(Form("cMatchedTrackJet%s_pt%d_%d",dl1pType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",dl1pType[b].c_str()), 100,0,1);
        m_hSvc.create1D(Form("lightMatchedTrackJet%s_pt%d_%d",dl1pType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",dl1pType[b].c_str()), 100,0,1);
      }
      for(int b=0;b<nDl1wType;b++){
        m_hSvc.create1D(Form("bMatchedSmallJet%s_pt%d_%d",dl1wType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",dl1wType[b].c_str()), 100,-5,10);
        m_hSvc.create1D(Form("cMatchedSmallJet%s_pt%d_%d",dl1wType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",dl1wType[b].c_str()), 100,-5,10);
        m_hSvc.create1D(Form("lightMatchedSmallJet%s_pt%d_%d",dl1wType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; small jet %s; Events",dl1wType[b].c_str()), 100,-5,10);
        m_hSvc.create1D(Form("bMatchedTrackJet%s_pt%d_%d",dl1wType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",dl1wType[b].c_str()), 100,-5,10);
        m_hSvc.create1D(Form("cMatchedTrackJet%s_pt%d_%d",dl1wType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",dl1wType[b].c_str()), 100,-5,10);
        m_hSvc.create1D(Form("lightMatchedTrackJet%s_pt%d_%d",dl1wType[b].c_str(),ptRange[p],ptRange[p+1]), Form("; track jet %s; Events",dl1wType[b].c_str()), 100,-5,10);
      }
    }

    m_hSvc.create1DVar("bMatchedVRTrackJetPt", "; VR track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedVRTrackJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedVRTrackJetPtPassVRTrackBtag", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedVRTrackJetParticlePtPassVRTrackBtag", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("bMatchedVRTrackJetMv2c10", "; VR track jet Mv2c10; Events", 100,-1,1);
    m_hSvc.create1DVar("cMatchedVRTrackJetPt", "; VR track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedVRTrackJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedVRTrackJetPtPassVRTrackBtag", "; VR track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedVRTrackJetParticlePtPassVRTrackBtag", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("cMatchedVRTrackJetMv2c10", "; VR track jet Mv2c10; Events", 100,-1,1);
    m_hSvc.create1DVar("lightMatchedVRTrackJetPt", "; VR track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedVRTrackJetParticlePt", "; VR track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedVRTrackJetPtPassVRTrackBtag", "; VR track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedVRTrackJetParticlePtPassVRTrackBtag", "; VR track jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("lightMatchedVRTrackJetMv2c10", "; VR track jet Mv2c10; Events", 100,-1,1);

    m_hSvc.create1DVar("bMatchedSmallJetPt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetPtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetPtPassCaloBtagHyb70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetPtPassCaloBtagMuFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetPtPassCaloBtagRnnFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetPtPassCaloBtagDl1Fix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetPtPassCaloBtagDl1muFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetPtPassCaloBtagDl1rnnFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetParticlePtPassCaloBtagFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetParticlePtPassCaloBtagHyb70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetParticlePtPassCaloBtagMuFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetParticlePtPassCaloBtagRnnFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetParticlePtPassCaloBtagDl1Fix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetParticlePtPassCaloBtagDl1muFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("bMatchedSmallJetParticlePtPassCaloBtagDl1rnnFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("bMatchedSmallJetMv2c10", "; small jet Mv2c10; Events", 100,-1,1);
    m_hSvc.create1DVar("cMatchedSmallJetPt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetParticlePt", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetPtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetPtPassCaloBtagHyb70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetPtPassCaloBtagMuFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetPtPassCaloBtagRnnFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetPtPassCaloBtagDl1Fix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetPtPassCaloBtagDl1muFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetPtPassCaloBtagDl1rnnFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetParticlePtPassCaloBtagFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetParticlePtPassCaloBtagHyb70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetParticlePtPassCaloBtagMuFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetParticlePtPassCaloBtagRnnFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetParticlePtPassCaloBtagDl1Fix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetParticlePtPassCaloBtagDl1muFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("cMatchedSmallJetParticlePtPassCaloBtagDl1rnnFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("cMatchedSmallJetMv2c10", "; small jet Mv2c10; Events", 100,-1,1);
    m_hSvc.create1DVar("lightMatchedSmallJetPt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetParticlePt", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetPtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetPtPassCaloBtagHyb70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetPtPassCaloBtagMuFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetPtPassCaloBtagRnnFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetPtPassCaloBtagDl1Fix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetPtPassCaloBtagDl1muFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetPtPassCaloBtagDl1rnnFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetParticlePtPassCaloBtagFix70", "; small jet p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetParticlePtPassCaloBtagHyb70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetParticlePtPassCaloBtagMuFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetParticlePtPassCaloBtagRnnFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetParticlePtPassCaloBtagDl1Fix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetParticlePtPassCaloBtagDl1muFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1DVar("lightMatchedSmallJetParticlePtPassCaloBtagDl1rnnFix70", "; truth p_{T} [GeV]; Events", varBinN8, varBin8);
    m_hSvc.create1D("lightMatchedSmallJetMv2c10", "; small jet Mv2c10; Events", 100,-1,1);

    std::string channel[6] = {"All","FullHad","FullLep","SemiLepTau","SemiLepMu","SemiLepEl"};
    for(int ch=0;ch<6;ch++) {
      m_hSvc.create1D(Form("cutFlow%s",channel[ch].c_str()), ";Number of events; step", 20,0,20);
      m_hSvc.create1D(Form("mtt%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mttSmooth%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mttDnn%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mlt%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mht%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mttReso%s",channel[ch].c_str()), "; mass of the top-antitop system [GeV]; Events", 60,0,6000);
      m_hSvc.create1D(Form("mltReso%s",channel[ch].c_str()), "; mass of leptonic top [GeV]; Events", 50,0,500);
      m_hSvc.create1D(Form("mhtReso%s",channel[ch].c_str()), "; mass of hadronic top [GeV]; Events", 50,0,500);
    }
    m_hSvc.create1DVar("mtt", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
    m_hSvc.create1D("chi2", "; #chi^{2}; Events", 50,0,10);
    m_hSvc.create1D("log10chi2", "; log#chi^{2}; Events", 40,-1,1);
    m_hSvc.create1DVar("mttReso", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
  }

AnaTtresMultiJets::~AnaTtresMultiJets() {
}

void AnaTtresMultiJets::run(const Event &evt, double weight, const std::string &s, int mode, TSpline* m_spline) {
  //std::cout << "AnaTtresMultiJets::run" << std::endl;
  //std::cout << "analysis mode = " << mode << std::endl;
  //std::cout << "m_boosted/m_electron=" << m_boosted << "/" << m_electron << std::endl;

  HistogramService *h = &m_hSvc;
  std::string suffix = s;

  //MET trigger
  bool pass_met_trigger = evt.trigger("HLT_xe110_mht_L1XE50"); 
  bool pass_pufit_met_trigger = evt.trigger("HLT_xe110_pufit_L1XE50"); 
  //std::cout << "pass_met_trigger=" << pass_met_trigger << std::endl;
  //std::cout << "pass_pufit_met_trigger=" << pass_pufit_met_trigger << std::endl;
  //Muon trigger
  bool pass_muon_trigger = evt.trigger("HLT_mu26_ivarmedium") || evt.trigger("HLT_mu50"); 
  //std::cout << "pass muon trigger=" << pass_muon_trigger << std::endl;
  //Electron trigger
  bool pass_electron_trigger = evt.trigger("HLT_e26_lhtight_nod0_ivarloose") || evt.trigger("HLT_e60_lhmedium_nod0") || evt.trigger("HLT_e140_lhloose_nod0"); 

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

  //MET
  TLorentzVector vmet = evt.met();
  float met = vmet.Perp()*1e-3;
  float mwt=0.;
  TLorentzVector l;
  if(m_electron && electron_id!=-1) l = evt.electron()[electron_id].mom();
  if(!m_electron && muon_id!=-1) l = evt.muon()[muon_id].mom();
  if(l.Perp()>1e-5) mwt = sqrt(2. * l.Perp() * vmet.Perp() * (1. - cos(vmet.DeltaPhi(l)) ))*1e-3;
  float wpt = (vmet+l).Perp()*1e-3;

  //Number of small jets
  int nSmallJet=0;
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    float jet_pt = evt.jet()[jidx].mom().Perp()*1e-3;
    if(jet_pt>25) nSmallJet++;
  }

  //Jet close to lepton
  /*size_t close_idx = 0;
  for (; close_idx < evt.jet().size(); ++close_idx){
    if (evt.jet()[close_idx].closeToLepton()) break;
  }
  bool jetCloseToLepton = (close_idx==evt.jet().size())? false : true;
  TLorentzVector sj;
  if(jetCloseToLepton) sj = evt.jet()[close_idx].mom();*/

  //Jet close to lepton
  size_t close_idx = 0;
  for (; close_idx < evt.jet().size(); ++close_idx){    
    float dphi = l.DeltaPhi(evt.jet()[close_idx].mom());
    float dy = l.Rapidity() - evt.jet()[close_idx].mom().Rapidity();
    float deltaR_tmp = std::sqrt(dy*dy + dphi*dphi);
    if(evt.jet()[close_idx].mom().Pt()*1e-3 > 25 && deltaR_tmp < 1.5) break;
  }//for     
  bool jetCloseToLepton = (close_idx==evt.jet().size())? false : true;
  TLorentzVector sj;
  if(jetCloseToLepton) sj = evt.jet()[close_idx].mom();

  //Number of large jets
  int nLargeJet=0,nLargeJetPt=0,nLargeJetSmooth=0,nLargeJetDnn=0;
  //std::cout << "# of large jet is " << evt.largeJet().size() << std::endl;
  for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
    const TLorentzVector &lj = evt.largeJet()[ljid].mom();
    bool good = evt.largeJet()[ljid].good_sub80();//default
    bool good_smooth = evt.largeJet()[ljid].good_smooth_mt80();//smooth mt 80
    bool good_dnn = evt.largeJet()[ljid].good_dnn80();//dnn 80
    //std::cout << "large jet pt/eta/phi=" << lj.Perp() << "/" << lj.Eta() << "/" << lj.Phi() << std::endl;
    //std::cout << "large jet mass/tau32=" << lj.M() << "/" << evt.largeJet()[ljid].subs("tau32_wta") << std::endl;;
    //std::cout << "good=" << good << std::endl;
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0 && good) nLargeJet++;
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0) nLargeJetPt++;
    if(lj.Perp()*1e-3 > 350 && fabs(lj.Eta())<2.0 && good_smooth) nLargeJetSmooth++;
    if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0 && good_dnn) nLargeJetDnn++;
  }

  //Angular cut
  int nGoodLargeJet=0,nGoodLargeJetSmooth=0,nGoodLargeJetDnn=0;
  int ljetid=-1,ljetid_smooth=-1,ljetid_dnn=-1;
  if(sj.Perp()>1e-5 && l.Perp()>1e-5){
    for (int ljid=0; ljid < evt.largeJet().size(); ++ljid) {
      const TLorentzVector &lj = evt.largeJet()[ljid].mom();
      bool good = evt.largeJet()[ljid].good_sub80();//default
      bool good_smooth = evt.largeJet()[ljid].good_smooth_mt80();//smooth mt 80
      bool good_dnn = evt.largeJet()[ljid].good_dnn80();//dnn 80
      float delta_phi_large_jet_lepton = lj.DeltaPhi(l);
      float delta_r_large_jet_small_jet = lj.DeltaR(sj);
      if(lj.Perp()*1e-3 > 300 && fabs(lj.Eta())<2.0 && fabs(delta_phi_large_jet_lepton)>2.3 && delta_r_large_jet_small_jet>1.5) {
        if(good){
          if(nGoodLargeJet==0) ljetid=ljid;
          nGoodLargeJet++;
        }
        if(good_smooth){
          if(nGoodLargeJetSmooth==0) ljetid_smooth=ljid;
          nGoodLargeJetSmooth++;
        }
        if(good_dnn){
          if(nGoodLargeJetDnn==0) ljetid_dnn=ljid;
          nGoodLargeJetDnn++;
        }
      }
    }
  }

  //Track btag
  int nTrackBtagJets=0;
  float temp_thr=-1;
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    float pt = evt.tjet()[jidx].mom().Pt()*1e-3;
    float eta = evt.tjet()[jidx].mom().Eta();
    float phi = evt.tjet()[jidx].mom().Phi();
    if(pt<300) temp_thr = 0.6455;
    else temp_thr = -0.3;
    //bool is_trackBtag = (isRel21)? evt.tjet()[jidx].btag_mv2c10_70_trk_rel21() : evt.tjet()[jidx].btag_mv2c10_70_trk();
    bool is_trackBtag = evt.tjet()[jidx].btag_mv2c10_70_trk();
    if (is_trackBtag && pt > 10 && fabs(eta) < 2.5 && evt.tjet()[jidx].numConstituents() >= 2) nTrackBtagJets++; 
    //if (is_trackBtag && pt > 10 && fabs(eta) < 2.5) nTrackBtagJets++; 
  }
  bool pass_fixed_btag = nTrackBtagJets > 0 ;
  //std::cout << "nTrackBtagJets=" << nTrackBtagJets << std::endl;

  //Calo B-tag
  int nCaloBtagJets=0,nCaloHybBtagJets,nCaloDl1RnnBtagJets=0,nCaloDl1BtagJets=0;
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    bool is_caloBtag = evt.jet()[jidx].mv2c10()>0.824;
    float pt = evt.jet()[jidx].mv2c10()*1e-3;
    float mv2c10 = evt.jet()[jidx].mv2c10();
    float dl1_pu = evt.jet()[jidx].dl1_pu();
    float dl1_pb = evt.jet()[jidx].dl1_pb();
    float dl1_pc = evt.jet()[jidx].dl1_pc();
    float dl1_weight = TuDoAtlas::calc_dl1_weight(dl1_pb, dl1_pc, dl1_pu, 1);
    float dl1rnn_pu = evt.jet()[jidx].dl1rnn_pu();
    float dl1rnn_pb = evt.jet()[jidx].dl1rnn_pb();
    float dl1rnn_pc = evt.jet()[jidx].dl1rnn_pc();
    float dl1rnn_weight = TuDoAtlas::calc_dl1_weight(dl1rnn_pb, dl1rnn_pc, dl1rnn_pu, 1);
    //std::cout << "dl1rnn_weight=" << dl1rnn_weight << std::endl;
    if(is_caloBtag) nCaloBtagJets++;
    if(mv2c10>m_spline->Eval(pt)) nCaloHybBtagJets++;
    if(dl1rnn_weight>2.98) nCaloDl1RnnBtagJets++;
    if(dl1_weight>2.01) nCaloDl1BtagJets++;
  }
  //std::cout << "nBtagJets calo/dl1rnn=" << nCaloBtagJets << "/" << nCaloDl1RnnBtagJets << std::endl;

  //Calo B-tag for resolved channel
  std::vector<TLorentzVector*> vjets;
  std::vector<bool> vjets_btagged;
  vjets.clear();vjets_btagged.clear();
  for (int jid=0; jid < evt.jet().size(); ++jid){    
    //std::cout << "jid=" << jid << std::endl;
    bool is_caloBtag = evt.jet()[jid].btag_mv2c10_70();
    //std::cout << "jet pt/eta/phi=" << evt.jet()[jid].mom().Perp() << "/" << evt.jet()[jid].mom().Eta() << "/" << evt.jet()[jid].mom().Phi() << std::endl;
    vjets.push_back(new TLorentzVector(0,0,0,0));
    vjets[jid]->SetPtEtaPhiE(evt.jet()[jid].mom().Perp(), evt.jet()[jid].mom().Eta(), evt.jet()[jid].mom().Phi(), evt.jet()[jid].mom().E());
    bool is_btagged(false);
    for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx)
    {
      TLorentzVector tmpTJet;
      tmpTJet.SetPtEtaPhiE(evt.tjet()[bidx].mom().Perp(),evt.tjet()[bidx].mom().Eta(),evt.tjet()[bidx].mom().Phi(),evt.tjet()[bidx].mom().E());
      if(tmpTJet.DeltaR(*vjets[jid]) <= 0.4) {
        //bool is_trackBtag = (isRel21)? evt.tjet()[bidx].btag_mv2c10_70_trk_rel21() : evt.tjet()[bidx].btag_mv2c10_70_trk();
        bool is_trackBtag = evt.tjet()[bidx].btag_mv2c10_70_trk();
        if(evt.tjet()[bidx].mom().Pt() > 10e3 && fabs(evt.tjet()[bidx].mom().Eta())<2.5 &&
            is_trackBtag && evt.tjet()[bidx].numConstituents()>=2)  {
          is_btagged = true ;
          break; 
        }
      } 
    } 
    vjets_btagged.push_back(is_btagged);
  }

  TLorentzVector tmpmet;
  tmpmet.SetPtEtaPhiM(evt.met().Perp(), 0, evt.met().Phi(), 0);
  int igj3, igj4, igb3, igb4, ign1;
  double chi2ming1, chi2ming1H, chi2ming1L;
  bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &tmpmet, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L); 

  float log10chi2 = (chi2ming1>0)? log10(chi2ming1) : 100;

  h->h1D("mu", "", suffix)->Fill(evt.mu());

  for(int t=0;t<evt.truth().size();t++){
    int id = evt.truth()[t].id();
    TLorentzVector par;
    par.SetPxPyPzE(evt.truth()[t].px(),evt.truth()[t].py(),evt.truth()[t].pz(),evt.truth()[t].e());
    if(abs(id)!=3) continue;
    //std::cout << "particle id/pt/eta/phi=" << id << "/" << par.Perp() << "/" << par.Eta() << "/" << par.Phi() << std::endl;
  }
  //std::cout << "multijet pass mu/e=" << evt.passes("bmujets_notrigger_nobtag_noang") << "/" << evt.passes("bejets_nobtag_noang") << std::endl;

  //some selections
  bool pass_bmujets = evt.passes("bmujets_2016");
  bool pass_bejets = evt.passes("bejets_2016");
  bool pass_bmujets_dl1btag = evt.passes("bmujets_dl1btag_2016");
  bool pass_bejets_dl1btag = evt.passes("bejets_dl1btag_2016");
  bool pass_bmujets_dl1rnnbtag = evt.passes("bmujets_dl1rnnbtag_2016");
  bool pass_bejets_dl1rnnbtag = evt.passes("bejets_dl1rnnbtag_2016");
  bool pass_bmujets_smoothtop = evt.passes("bmujets_smoothtop_2016");
  bool pass_bejets_smoothtop = evt.passes("bejets_smoothtop_2016");
  bool pass_bmujets_dnntop = evt.passes("bmujets_dnntop_2016");
  bool pass_bejets_dnntop = evt.passes("bejets_dnntop_2016");
  //std::cout << "muon def/dl1/dl1rnn/smoothtop/dnntop=" 
    //<< pass_bmujets << "/" << pass_bmujets_dl1btag << "/" << pass_bmujets_dl1rnnbtag << "/" << pass_bmujets_smoothtop << "/" << pass_bmujets_dnntop << std::endl;
  //std::cout << "electron def/dl1/dl1rnn/smoothtop/dnntop=" 
    //<< pass_bejets << "/" << pass_bejets_dl1btag << "/" << pass_bejets_dl1rnnbtag << "/" << pass_bejets_smoothtop << "/" << pass_bejets_dnntop << std::endl;

  //Boosted Muon 
  if(!m_electron && m_boosted){
    if(mode==0){
      if(!evt.passes("bmujets_2015") && !evt.passes("bmujets_2016")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
  }

  //Boosted Electron 
  if(m_electron && m_boosted){
    if(mode==0){
      if(!evt.passes("bejets_2015") && !evt.passes("bejets_2016")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
  }

  TLorentzVector ljt,ljt_smooth,ljt_dnn;
  if(nGoodLargeJet!=0) ljt = evt.jet()[ljetid].mom();
  if(nGoodLargeJetSmooth!=0) ljt_smooth = evt.jet()[ljetid_smooth].mom();
  if(nGoodLargeJetDnn!=0) ljt_dnn = evt.jet()[ljetid_dnn].mom();

  //Lepton
  if(m_electron) h->h1D("nLepton", "", suffix)->Fill(nElectron30, weight);
  else h->h1D("nLepton", "", suffix)->Fill(nMuon30, weight);
  if(l.Perp()>1e-5){
    h->h1D("lepPt", "", suffix)->Fill(l.Perp()*1e-3, weight);
    h->h1D("lepEta", "", suffix)->Fill(l.Eta(), weight);
    h->h1D("lepPhi", "", suffix)->Fill(l.Phi(), weight);
  }

  //MET
  h->h1D("MET", "", suffix)->Fill(met, weight);
  h->h1D("metPhi", "", suffix)->Fill(vmet.Phi(), weight);
  if(l.Perp()>1e-5){
    h->h1D("mwt", "", suffix)->Fill(mwt, weight);
    h->h1D("wpt", "", suffix)->Fill(wpt, weight);
    if(pass_met_trigger) h->h1D("wptPassMet110", "", suffix)->Fill(wpt, weight);
    if(pass_pufit_met_trigger) h->h1D("wptPassPufitMet110", "", suffix)->Fill(wpt, weight);
    h->h1D("metPlusMwt", "", suffix)->Fill(met+mwt, weight);
  }

  //small-R jet
  //std::cout << "# of small-R jet is " << evt.jet().size() << std::endl;
  h->h1D("nSmallJet", "", suffix)->Fill(evt.jet().size(), weight);
  if(sj.Perp()>1e-5){
    h->h1D("smallJetPt", "", suffix)->Fill(sj.Perp()*1e-3, weight);
    h->h1D("smallJetEta", "", suffix)->Fill(sj.Eta(), weight);
    h->h1D("smallJetPhi", "", suffix)->Fill(sj.Phi(), weight);
    h->h1D("smallJetMass", "", suffix)->Fill(sj.M()*1e-3, weight);
  }
  if(evt.jet().size()>0){
    const TLorentzVector &lead_sj = evt.jet()[0].mom();
    h->h1D("leadSmallJetPt", "", suffix)->Fill(lead_sj.Perp()*1e-3, weight);
    h->h1D("leadSmallJetEta", "", suffix)->Fill(lead_sj.Eta(), weight);
    h->h1D("leadSmallJetPhi", "", suffix)->Fill(lead_sj.Phi(), weight);
    h->h1D("leadSmallJetMass", "", suffix)->Fill(lead_sj.M()*1e-3, weight);
  }
  if(evt.jet().size()>1){
    const TLorentzVector &second_lead_sj = evt.jet()[1].mom();
    h->h1D("secondLeadSmallJetPt", "", suffix)->Fill(second_lead_sj.Perp()*1e-3, weight);
    h->h1D("secondLeadSmallJetEta", "", suffix)->Fill(second_lead_sj.Eta(), weight);
    h->h1D("secondLeadSmallJetPhi", "", suffix)->Fill(second_lead_sj.Phi(), weight);
    h->h1D("secondLeadSmallJetMass", "", suffix)->Fill(second_lead_sj.M()*1e-3, weight);
  }
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    h->h1D("allSmallJetPt", "", suffix)->Fill(evt.jet()[jidx].mom().Perp()*1e-3, weight);
    h->h1D("allSmallJetEta", "", suffix)->Fill(evt.jet()[jidx].mom().Eta(), weight);
    h->h1D("allSmallJetPhi", "", suffix)->Fill(evt.jet()[jidx].mom().Phi(), weight);
  }

  //largeR jet
  //std::cout << "# of large-R jet is " << evt.largeJet().size() << std::endl;
  h->h1D("nLargeJet", "", suffix)->Fill(nLargeJet, weight);
  if(ljt.Perp()>1e-5){
    h->h1D("largeJetPt", "", suffix)->Fill(ljt.Perp()*1e-3, weight);
    h->h1D("largeJetEta", "", suffix)->Fill(ljt.Eta(), weight);
    h->h1D("largeJetPhi", "", suffix)->Fill(ljt.Phi(), weight);
    h->h1D("largeJetMass", "", suffix)->Fill(ljt.M()*1e-3, weight);
    //h->h1D("largeJetTau32", "", suffix)->Fill(evt.largeJet()[ljetid].subs("tau32_wta"), weight);
  }

  //kinematic values
  if(l.Perp()>1e-3 && sj.Perp()>1e-3) h->h1D("deltaR_lepSmallJet", "",suffix)->Fill(l.DeltaR(sj),weight);
  if(ljt.Perp()>1e-3 && sj.Perp()>1e-3) h->h1D("deltaR_largeJetSmallJet", "",suffix)->Fill(ljt.DeltaR(sj),weight);
  if(ljt.Perp()>1e-3 && l.Perp()>1e-3) h->h1D("deltaPhi_largeJetLepton", "",suffix)->Fill(fabs(ljt.DeltaPhi(l)),weight);

  //track jet
  //std::cout << "# of track jet is " << evt.tjet().size() << std::endl;
  for(int itj=0;itj<evt.tjet().size();itj++){
    TLorentzVector tjet = evt.tjet()[itj].mom();
    h->h1D("trackJetPt", "", suffix)->Fill(tjet.Perp()*1e-3, weight);
    h->h1D("trackJetEta", "", suffix)->Fill(tjet.Eta(), weight);
    h->h1D("trackJetPhi", "", suffix)->Fill(tjet.Phi(), weight);
    h->h1D("trackJetMv2c10", "", suffix)->Fill(evt.tjet()[itj].mv2c10(), weight);
  }
  h->h1D("nTrackJet", "", suffix)->Fill(evt.tjet().size(), weight);
  h->h1D("nBtaggedTrackJet", "", suffix)->Fill(nTrackBtagJets, weight);

  //mass of ttbar (boosted)
  std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), true);   
  TLorentzVector nu(0,0,0,0);
  if (vec_nu.size() > 0) {
    nu = *(vec_nu[0]);
    for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
    vec_nu.clear();
  }
  float mtt = (ljt+sj+nu+l).M();
  float mtt_smooth = (ljt_smooth+sj+nu+l).M();
  float mtt_dnn = (ljt_dnn+sj+nu+l).M();
  float mlt = (sj+nu+l).M();
  float mht = ljt.M();
  if(ljt.Perp()>1e-3 && sj.Perp()>1e-3 && l.Perp()>1e-3){
    h->h1D("mtt", "", suffix)->Fill(mtt*1e-3, weight);
    h->h1D("mttAll", "", suffix)->Fill(mtt*1e-3, weight);
  }
  if(ljt_smooth.Perp()>1e-3 && sj.Perp()>1e-3 && l.Perp()>1e-3){
    h->h1D("mttSmoothAll", "", suffix)->Fill(mtt*1e-3, weight);
  }
  if(ljt_dnn.Perp()>1e-3 && sj.Perp()>1e-3 && l.Perp()>1e-3){
    h->h1D("mttDnnAll", "", suffix)->Fill(mtt*1e-3, weight);
  }
  if(sj.Perp()>1e-3 && l.Perp()>1e-3){
    h->h1D("mltAll", "", suffix)->Fill(mlt*1e-3, weight);
  }
  if(ljt.Perp()>1e-3){
    h->h1D("mhtAll", "", suffix)->Fill(mht*1e-3, weight);
  }

  //mass of ttbar (resolved)
  if (status){
    float mlt_reso = m_chi2.getResult_Mtl();
    float mht_reso = m_chi2.getResult_Mth();
    float mtt_reso = m_chi2.getResult_Mtt();
    h->h1D("mttReso", "", suffix)->Fill(mtt_reso*1e-3, weight);
    h->h1D("mttResoAll", "", suffix)->Fill(mtt_reso*1e-3, weight);
    h->h1D("mltResoAll", "", suffix)->Fill(mlt_reso*1e-3, weight);
    h->h1D("mhtResoAll", "", suffix)->Fill(mht_reso*1e-3, weight);
    h->h1D("chi2", "", suffix)->Fill(chi2ming1, weight);
    h->h1D("log10chi2", "", suffix)->Fill(log10chi2, weight);
    //std::cout << "resolved mlt/mht/mtt=" << mlt_reso << "/" << mht_reso << "/" << mtt_reso << std::endl;
  }

}

