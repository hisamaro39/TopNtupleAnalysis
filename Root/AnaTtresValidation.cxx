/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
//#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnalysisValidation.h"
#include "TopNtupleAnalysis/AnaTtresValidation.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"
#include "TopNtupleAnalysis/UserFunktions.h"


AnaTtresValidation::AnaTtresValidation(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
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

    m_hSvc.create1D("met", "; missing E_[T] [GeV]; Events", 25, 0, 1000);
    m_hSvc.create1D("metPhi", "; missing E_{T} #phi ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D("mwt", "; M_{T,W} [GeV]; Events", 25, 0, 250);
    m_hSvc.create1D("wpt", "; p_{T,W} [GeV]; Events", 100, 0, 1000);
    m_hSvc.create1D("wptPassMet110", "; p_{T,W} [GeV]; Events", 100, 0, 1000);
    m_hSvc.create1D("wptPassPufitMet110", "; p_{T,W} [GeV]; Events", 100, 0, 1000);

    m_hSvc.create1D("deltaR_lepSmallJet", "; #DeltaR_{lep,sjet}; Events", 20, 0, 2);
    m_hSvc.create1D("deltaR_largeJetSmallJet", "; #DeltaR_{ljet,sjet}; Events", 50, 1, 6);
    m_hSvc.create1D("deltaPhi_largeJetLepton", "; #Delta#phi_{ljet,lep}; Events", 24, 2, 3.2);

    m_hSvc.create1D("deltaR_trackJetTruthHB", "; #DeltaR_{lep,sjet}; Events", 60, 0, 6);
    m_hSvc.create1D("deltaR_trackJetTruthLB", "; #DeltaR_{lep,sjet}; Events", 60, 0, 6);

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

AnaTtresValidation::~AnaTtresValidation() {
}

void AnaTtresValidation::run(const Event &evt, double weight, const std::string &s, int mode, TSpline* m_spline) {
  //std::cout << "AnaTtresValidation::run" << std::endl;
  //std::cout << "analysis mode = " << mode << std::endl;
  //std::cout << "m_boosted/m_electron=" << m_boosted << "/" << m_electron << std::endl;

  float btag_thr_calo70 = 0.831;
  if(mode==0 || mode==3 || mode==14 || mode==17) btag_thr_calo70 = 0.824;

  //std::cout << "isRel21=" << isRel21 << std::endl;
  int Wdecay1_from_t_pdgId = evt.MC_Wdecay1_from_t_pdgId();
  int Wdecay2_from_t_pdgId = evt.MC_Wdecay2_from_t_pdgId();
  int Wdecay1_from_tbar_pdgId = evt.MC_Wdecay1_from_tbar_pdgId();
  int Wdecay2_from_tbar_pdgId = evt.MC_Wdecay2_from_tbar_pdgId();
  bool isLeptonicT = (abs(Wdecay1_from_t_pdgId)<7)? false : true;
  bool isLeptonicTbar = (abs(Wdecay1_from_tbar_pdgId)<7)? false : true;
  TLorentzVector b_from_t,b_from_tbar,t,tbar,W_from_t,W_from_tbar,Wdecay1_from_t,Wdecay2_from_t,Wdecay1_from_tbar,Wdecay2_from_tbar;
  b_from_t = evt.MC_b_from_t();
  b_from_tbar = evt.MC_b_from_tbar();
  W_from_t = evt.MC_Wdecay1_from_t() + evt.MC_Wdecay2_from_t(); 
  W_from_tbar = evt.MC_Wdecay1_from_tbar() + evt.MC_Wdecay2_from_tbar(); 
  Wdecay1_from_t = evt.MC_Wdecay1_from_t();
  Wdecay2_from_t = evt.MC_Wdecay2_from_t();
  Wdecay1_from_tbar = evt.MC_Wdecay1_from_tbar();
  Wdecay2_from_tbar = evt.MC_Wdecay2_from_tbar();
  t = evt.MC_t();
  tbar = evt.MC_tbar();
  float t_pt = t.Perp()*1e-3, t_eta = t.Eta(), t_phi = t.Phi();
  float tbar_pt = tbar.Perp()*1e-3, tbar_eta = tbar.Eta(), tbar_phi = tbar.Phi();
  float b_from_t_pt = b_from_t.Perp()*1e-3, b_from_t_eta = b_from_t.Eta(), b_from_t_phi = b_from_t.Phi();
  float b_from_tbar_pt = b_from_tbar.Perp()*1e-3, b_from_tbar_eta = b_from_tbar.Eta(), b_from_tbar_phi = b_from_tbar.Phi();
  //if(Wdecay1_from_t_pdgId==0) return;//temporal

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
  else std::cout << "Unknow channel" << std::endl;
  //std::cout << "channel is " << channel << std::endl;

  TLorentzVector b_from_htop,b_from_ltop,htop,ltop,Wdecay1_from_htop,Wdecay2_from_htop;
  if(isSemiLep){
    if(Wdecay2_from_t_pdgId==-13 || Wdecay2_from_t_pdgId==-11 || Wdecay2_from_t_pdgId==-15){//t is leptonic
      b_from_ltop = b_from_t;
      ltop = t;
      b_from_htop = b_from_tbar;
      htop = W_from_tbar + b_from_tbar;
      Wdecay1_from_htop = Wdecay1_from_tbar;
      Wdecay2_from_htop = Wdecay2_from_tbar;
    }
    if(Wdecay2_from_t_pdgId==-1 || Wdecay2_from_t_pdgId==-3){//t is hadronic
      b_from_ltop = b_from_tbar;
      ltop = tbar;
      b_from_htop = b_from_t;
      htop = W_from_t + b_from_t;
      Wdecay1_from_htop = Wdecay1_from_t;
      Wdecay2_from_htop = Wdecay2_from_t;
    }
  }
  //std::cout << "Wdecay1/Wdecay2 from t is " << Wdecay1_from_t_pdgId << "/" << Wdecay2_from_t_pdgId << std::endl;
  //std::cout << "Wdecay1/Wdecay2 from tbar is " << Wdecay1_from_tbar_pdgId << "/" << Wdecay2_from_tbar_pdgId << std::endl;
  //std::cout << "pt of ltop/htop=" << ltop.Perp() << "/" << htop.Perp() << std::endl;

  HistogramService *h = &m_hSvc;
  bool isTight = false;
  std::string suffix = s;

  if(channel=="FullHad") h->h1D("decayType", "", suffix)->Fill(0);
  if(channel=="SemiLepEl") h->h1D("decayType", "", suffix)->Fill(1);
  if(channel=="SemiLepMu") h->h1D("decayType", "", suffix)->Fill(2);
  if(channel=="SemiLepTau") h->h1D("decayType", "", suffix)->Fill(3);
  if(channel=="FullLep") h->h1D("decayType", "", suffix)->Fill(4);
  if(channel!="SemiLepMu") return;
  
  h->h1D("Wdecay1_from_t_pdgId", "", suffix)->Fill(Wdecay1_from_t_pdgId);

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
    bool is_caloBtag = evt.jet()[jidx].mv2c10()>btag_thr_calo70;
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
    if(mode==0 || mode==1 || mode==12 || mode==13 || mode==18 || mode==19 || mode==22){//apply all selection
      if(!evt.passes("no_cut")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      if(!pass_muon_trigger) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(1);
      //if(evt.randomRunNumber()<297730) return;
      //std::cout << "randomRunNumber is " << evt.randomRunNumber() << std::endl;
      h->h1D("cutFlowAll", "", suffix)->Fill(2);
      if(nMuon30==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(3);
      if(nMuon25>1) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(4);
      if(nElectron25!=0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(5);
      if(!evt.muon()[muon_id].HLT_mu26_ivarmedium() && !evt.muon()[muon_id].HLT_mu50()) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(6);
      if(!(evt.passes("jet_clean"))) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(7);
      if(met<20) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(8);
      if(met+mwt<60) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(9);
      if(nSmallJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(10);
      /*if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(mode==22 && nLargeJetSmooth==0) return;
      else if(nLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(mode==22 && nGoodLargeJetSmooth==0) return;
      else if(nGoodLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(mode==12 && nCaloBtagJets==0) return;
      else if(mode==13 && nCaloHybBtagJets==0) return;
      else if(mode==18 && nCaloDl1BtagJets==0) return;
      else if(mode==19 && nCaloDl1RnnBtagJets==0) return;
      else if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);*/
    }
    if(mode==2){//use met trigger
      if(!evt.passes("bmujets_notrigger_nobtag")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      if(wpt<220) {
        if(!pass_muon_trigger) return;
      }
      else if(!pass_met_trigger) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(1);
      if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(2);
    }
    if(mode==3){//check tagging
      //if(channel!="SemiLepMu") return;
      if(!evt.passes("bmujets_check_tagging")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      /*if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(nLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(nGoodLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);*/
    }
    if(mode==4){//check met trigger efficiency
      if(!evt.passes("no_cut")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      //if(!pass_muon_trigger) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(1);
      //if(evt.randomRunNumber()<297730) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(2);
      if(nMuon30==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(3);
      if(nMuon25>1) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(4);
      if(nElectron25!=0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(5);
      //if(!evt.muon()[muon_id].HLT_mu26_ivarmedium() && !evt.muon()[muon_id].HLT_mu50()) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(6);
      if(!(evt.passes("jet_clean"))) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(7);
      if(met<20) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(8);
      if(met+mwt<60) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(9);
      if(nSmallJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(10);
      //if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(nLargeJetPt==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      //if(nGoodLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      //if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);
    }
    if(mode==5 || mode==6){//apply all selection use met trigger
      if(!evt.passes("no_cut")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      //if(!pass_muon_trigger) return;
      if(wpt<220) {
        if(!pass_muon_trigger) return;
      }else {
        if(mode==5 && !pass_met_trigger) return;
        if(mode==6 && !pass_pufit_met_trigger) return;
      }
      h->h1D("cutFlowAll", "", suffix)->Fill(1);
      //if(evt.randomRunNumber()<297730) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(2);
      if(nMuon30==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(3);
      if(nMuon25>1) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(4);
      if(nElectron25!=0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(5);
      //if(!evt.muon()[muon_id].HLT_mu26_ivarmedium() && !evt.muon()[muon_id].HLT_mu50()) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(6);
      if(!(evt.passes("jet_clean"))) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(7);
      if(met<20) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(8);
      if(met+mwt<60) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(9);
      if(nSmallJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(10);
      if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(nLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(nGoodLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);
    }
    if(mode==9){//default
      if(!evt.passes("bmujets_2016")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
    if(mode==10){//default
      //if(!evt.passes("bmujets_notrigger_nobtag")) return;
      if(!evt.passes("bmujets_2016")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
    if(mode==14 || mode==15 || mode==16){//use new b-tag
      if(mode==14 && !pass_bmujets) return;
      if(mode==15 && !pass_bmujets_dl1rnnbtag) return;
      if(mode==16 && !pass_bmujets_dl1btag) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      /*if(!evt.passes("bmujets_check_tagging")) return;
      if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(nLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(nGoodLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(mode==14 && nCaloBtagJets==0) return;
      if(mode==15 && nCaloDl1RnnBtagJets==0) return;
      if(mode==16 && nCaloDl1BtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);*/
    }
    if(mode==20 || mode==21){//use new top-tag
      if(!evt.passes("bmujets_check_tagging")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(mode==20 && nLargeJet==0) return;
      if(mode==21 && nLargeJetSmooth==0) return;
      if(mode==24 && nLargeJetDnn==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if((mode==20) && nGoodLargeJet==0) return;
      if(mode==21 && nGoodLargeJetSmooth==0) return;
      if(mode==24 && nGoodLargeJetDnn==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);
    }
    if(mode==23 || mode==24){//use new top-tag
      if(mode==23 && !pass_bmujets_smoothtop) return;
      if(mode==24 && !pass_bmujets_dnntop) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
  }

  //Boosted Electron 
  if(m_electron && m_boosted){
    if(mode==0 || mode==1 || mode==12 || mode==13 || mode==18 || mode==19 || mode==22){//apply all selection
      if(!evt.passes("no_cut")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      if(!pass_electron_trigger) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(1);
      //if(evt.randomRunNumber()<297730) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(2);
      if(nElectron30==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(3);
      if(nElectron25>1) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(4);
      if(nMuon25!=0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(5);
      if(!evt.electron()[electron_id].HLT_e26_lhtight_nod0_ivarloose() && !evt.electron()[electron_id].HLT_e60_lhmedium_nod0() 
          && !evt.electron()[electron_id].HLT_e140_lhloose_nod0()) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(6);
      if(!(evt.passes("jet_clean"))) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(7);
      if(met<20) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(8);
      if(met+mwt<60) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(9);
      if(nSmallJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(10);
      if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(mode==22 && nLargeJetSmooth==0) return;
      else if(nLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(mode==22 && nGoodLargeJetSmooth==0) return;
      else if(nGoodLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(mode==12 && nCaloBtagJets==0) return;
      else if(mode==13 && nCaloHybBtagJets==0) return;
      else if(mode==18 && nCaloDl1BtagJets==0) return;
      else if(mode==19 && nCaloDl1RnnBtagJets==0) return;
      else if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);
    }
    if(mode==2){//use met trigger
      if(!evt.passes("bejets_2016")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
    if(mode==3){//check tagging
      //if(channel!="SemiLepEl") return;
      if(!evt.passes("bejets_check_tagging")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
    if(mode==9){//default
      if(!evt.passes("bejets_2016")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
    if(mode==10){//default
      //if(!evt.passes("bejets_nobtag_noang")) return;
      if(!evt.passes("bejets_2016")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
    }
    if(mode==14 || mode==15 || mode==16){//use new b-tag
      if(mode==14 && !pass_bejets) return;
      if(mode==15 && !pass_bejets_dl1rnnbtag) return;
      if(mode==16 && !pass_bejets_dl1btag) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      /*if(!evt.passes("bejets_check_tagging")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(nLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(nGoodLargeJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(mode==14 && nCaloBtagJets==0) return;
      if(mode==15 && nCaloDl1RnnBtagJets==0) return;
      if(mode==16 && nCaloDl1BtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);*/
    }
    if(mode==20 || mode==21){//use new top-tag
      if(!evt.passes("bejets_check_tagging")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      if(!jetCloseToLepton) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(mode==20 && nLargeJet==0) return;
      if(mode==21 && nLargeJetSmooth==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(mode==20 && nGoodLargeJet==0) return;
      if(mode==21 && nGoodLargeJetSmooth==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);
    }
  }

  //Resolved Muon 
  if(!m_electron && !m_boosted){
    if(mode==0 || mode==1){//apply all selection
      if(!evt.passes("no_cut")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      if(!pass_muon_trigger) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(1);
      //if(evt.randomRunNumber()<297730) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(2);
      if(nMuon30==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(3);
      if(nMuon25>1) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(4);
      if(nElectron25!=0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(5);
      if(!evt.muon()[muon_id].HLT_mu26_ivarmedium() && !evt.muon()[muon_id].HLT_mu50()) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(6);
      if(!(evt.passes("jet_clean"))) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(7);
      if(met<20) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(8);
      if(met+mwt<60) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(9);
      if(nSmallJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(10);
      if(nSmallJet==1) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(nSmallJet==2) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(nSmallJet==3) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);
      if(log10chi2>0.9) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(15);
      //std::cout << "log10chi2=" << log10chi2 << std::endl;
    }
  }

  //Resolved Electron 
  if(m_electron && !m_boosted){
    if(mode==0 || mode==1){//apply all selection
      if(!evt.passes("no_cut")) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(0);
      if(!pass_electron_trigger) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(1);
      //if(evt.randomRunNumber()<297730) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(2);
      if(nElectron30==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(3);
      if(nElectron25>1) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(4);
      if(nMuon25!=0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(5);
      if(!evt.electron()[electron_id].HLT_e26_lhtight_nod0_ivarloose() && !evt.electron()[electron_id].HLT_e60_lhmedium_nod0() 
          && !evt.electron()[electron_id].HLT_e140_lhloose_nod0()) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(6);
      if(!(evt.passes("jet_clean"))) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(7);
      if(met<20) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(8);
      if(met+mwt<60) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(9);
      if(nSmallJet==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(10);
      if(nSmallJet==1) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(11);
      if(nSmallJet==2) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(12);
      if(nSmallJet==3) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(13);
      if(nTrackBtagJets==0) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(14);
      if(log10chi2>0.9) return;
      h->h1D("cutFlowAll", "", suffix)->Fill(15);
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
  h->h1D("met", "", suffix)->Fill(met, weight);
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
    h->h1D(Form("mtt%s",channel.c_str()), "", suffix)->Fill(mtt*1e-3, weight);
  }
  if(ljt_smooth.Perp()>1e-3 && sj.Perp()>1e-3 && l.Perp()>1e-3){
    h->h1D("mttSmoothAll", "", suffix)->Fill(mtt*1e-3, weight);
  }
  if(ljt_dnn.Perp()>1e-3 && sj.Perp()>1e-3 && l.Perp()>1e-3){
    h->h1D("mttDnnAll", "", suffix)->Fill(mtt*1e-3, weight);
  }
  if(sj.Perp()>1e-3 && l.Perp()>1e-3){
    h->h1D("mltAll", "", suffix)->Fill(mlt*1e-3, weight);
    h->h1D(Form("mlt%s",channel.c_str()), "", suffix)->Fill(mlt*1e-3, weight);
  }
  if(ljt.Perp()>1e-3){
    h->h1D(Form("mht%s",channel.c_str()), "", suffix)->Fill(mht*1e-3, weight);
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


  //track jet and hb and lb matching
  if(isSemiLep){
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    float pt = evt.tjet()[jidx].mom().Pt()*1e-3;
    float eta = evt.tjet()[jidx].mom().Eta();
    float phi = evt.tjet()[jidx].mom().Phi();
    float mv2c10 = evt.tjet()[jidx].mv2c10();
    float dr_hb = TuDoAtlas::calc_delta_r(eta,b_from_htop.Eta(),phi,b_from_htop.Phi());
    float dr_lb = TuDoAtlas::calc_delta_r(eta,b_from_ltop.Eta(),phi,b_from_ltop.Phi());
    h->h1D("deltaR_trackJetTruthHB", "", suffix)->Fill(dr_hb, weight);
    h->h1D("deltaR_trackJetTruthLB", "", suffix)->Fill(dr_lb, weight);
    //std::cout << "track jet pt/eta/phi/dr_hb/dr_lb=" << pt << "/" << eta << "/" << phi << "/" << dr_hb << "/" << dr_lb << std::endl;
    if(dr_lb<0.3){
      h->h1D("lbMatchedTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("lbMatchedTrackJetParticlePt", "", suffix)->Fill(b_from_ltop.Perp()*1e-3, weight);
      h->h1D("lbMatchedTrackJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      if(mv2c10>0.6455) h->h1D("lbMatchedTrackJetPtPassTrackBtag", "", suffix)->Fill(pt, weight);
      if(mv2c10>0.6455) h->h1D("lbMatchedTrackJetParticlePtPassTrackBtag", "", suffix)->Fill(b_from_ltop.Perp()*1e-3, weight);
    }
    if(dr_hb<0.3){
      h->h1D("hbMatchedTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("hbMatchedTrackJetParticlePt", "", suffix)->Fill(b_from_htop.Perp()*1e-3, weight);
      h->h1D("hbMatchedTrackJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      if(mv2c10>0.6455) h->h1D("hbMatchedTrackJetPtPassTrackBtag", "", suffix)->Fill(pt, weight);
      if(mv2c10>0.6455) h->h1D("hbMatchedTrackJetParticlePtPassTrackBtag", "", suffix)->Fill(b_from_htop.Perp()*1e-3, weight);
    }
    if(dr_lb>0.3 && dr_hb>0.3){
      h->h1D("otherMatchedTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("otherMatchedTrackJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      if(evt.tjet()[jidx].mv2c10()>0.6455) h->h1D("otherMatchedTrackJetPtPassTrackBtag", "", suffix)->Fill(pt, weight);
    }
  }
  }

  //small jet and hb and lb matching
  if(isSemiLep){
    int lb_match_id=-1,hb_match_id=-1;
    float min_dr_lb=999,min_dr_hb=999;
    for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
      float pt = evt.jet()[jidx].mom().Pt()*1e-3;
      float eta = evt.jet()[jidx].mom().Eta();
      float phi = evt.jet()[jidx].mom().Phi();
      float dr_hb = TuDoAtlas::calc_delta_r(eta,b_from_htop.Eta(),phi,b_from_htop.Phi());
      float dr_lb = TuDoAtlas::calc_delta_r(eta,b_from_ltop.Eta(),phi,b_from_ltop.Phi());
      float mv2c10 = evt.jet()[jidx].mv2c10();
      if(dr_lb<min_dr_lb){
        min_dr_lb=dr_lb;
        lb_match_id=jidx;
      }
      if(dr_hb<min_dr_hb){
        min_dr_hb=dr_hb;
        hb_match_id=jidx;
      }
      if(dr_lb<0.3){
        h->h1D("lbMatchedSmallJetPt", "", suffix)->Fill(pt, weight);
        h->h1D("lbMatchedSmallJetParticlePt", "", suffix)->Fill(b_from_ltop.Perp()*1e-3, weight);
        h->h1D("lbMatchedSmallJetMv2c10", "", suffix)->Fill(mv2c10, weight);
        if(mv2c10>btag_thr_calo70) h->h1D("lbMatchedSmallJetPtPassCaloBtagFix70", "", suffix)->Fill(pt, weight);
        if(mv2c10>btag_thr_calo70) h->h1D("lbMatchedSmallJetParticlePtPassCaloBtagFix70", "", suffix)->Fill(b_from_ltop.Perp()*1e-3, weight);
        if(mv2c10>m_spline->Eval(pt)) h->h1D("lbMatchedSmallJetPtPassCaloBtagHyb70", "", suffix)->Fill(pt, weight);
        //std::cout << "lb matched small jet pt/eta/phi/dr_hb/dr_lb=" << pt << "/" << eta << "/" << phi << "/" << dr_hb << "/" << dr_lb << std::endl;
      }
      if(dr_hb<0.3){
        h->h1D("hbMatchedSmallJetPt", "", suffix)->Fill(pt, weight);
        h->h1D("hbMatchedSmallJetParticlePt", "", suffix)->Fill(b_from_htop.Perp()*1e-3, weight);
        h->h1D("hbMatchedSmallJetMv2c10", "", suffix)->Fill(mv2c10, weight);
        if(mv2c10>btag_thr_calo70) h->h1D("hbMatchedSmallJetPtPassCaloBtagFix70", "", suffix)->Fill(pt, weight);
        if(mv2c10>btag_thr_calo70) h->h1D("hbMatchedSmallJetParticlePtPassCaloBtagFix70", "", suffix)->Fill(b_from_htop.Perp()*1e-3, weight);
        if(mv2c10>m_spline->Eval(pt)) h->h1D("hbMatchedSmallJetPtPassCaloBtagHyb70", "", suffix)->Fill(pt, weight);
        //std::cout << "hb matched small jet pt/eta/phi/dr_hb/dr_lb=" << pt << "/" << eta << "/" << phi << "/" << dr_hb << "/" << dr_lb << std::endl;
      }
      if(dr_lb>0.3 && dr_hb>0.3){
        h->h1D("otherMatchedSmallJetPt", "", suffix)->Fill(pt, weight);
        h->h1D("otherMatchedSmallJetMv2c10", "", suffix)->Fill(evt.jet()[jidx].mv2c10(), weight);
        if(mv2c10>btag_thr_calo70) h->h1D("otherMatchedSmallJetPtPassCaloBtagFix70", "", suffix)->Fill(pt, weight);
        if(evt.jet()[jidx].mv2c10()>m_spline->Eval(pt)) h->h1D("otherMatchedSmallJetPtPassCaloBtagHyb70", "", suffix)->Fill(pt, weight);
      }
    }
    if(min_dr_lb<0.3){
      float pt = evt.jet()[lb_match_id].mom().Perp()*1e-3; 
      float mv2c10 = evt.jet()[lb_match_id].mv2c10();
      h->h1D("lbMatchedClosestSmallJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("lbMatchedClosestSmallJetParticlePt", "", suffix)->Fill(b_from_ltop.Perp()*1e-3, weight);
      if(mv2c10>btag_thr_calo70) h->h1D("lbMatchedClosestSmallJetPtPassCaloBtagFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10>btag_thr_calo70) h->h1D("lbMatchedClosestSmallJetParticlePtPassCaloBtagFix70", "", suffix)->Fill(b_from_ltop.Perp()*1e-3, weight);
    }
    if(min_dr_hb<0.3){
      float pt = evt.jet()[hb_match_id].mom().Perp()*1e-3; 
      float mv2c10 = evt.jet()[hb_match_id].mv2c10();
      h->h1D("hbMatchedClosestSmallJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("hbMatchedClosestSmallJetParticlePt", "", suffix)->Fill(b_from_ltop.Perp()*1e-3, weight);
      if(mv2c10>btag_thr_calo70) h->h1D("hbMatchedClosestSmallJetPtPassCaloBtagFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10>btag_thr_calo70) h->h1D("hbMatchedClosestSmallJetParticlePtPassCaloBtagFix70", "", suffix)->Fill(b_from_ltop.Perp()*1e-3, weight);
    }
  }

  //std::cout << "Number of truth jet is " << evt.truthjet().size() << std::endl;
  for (int t=0;t<evt.truthjet().size();t++){
    TLorentzVector truth_jet = evt.truthjet()[t].mom();
    //std::cout << "truth jet pt/eta/phi=" << truth_jet.Perp() << "/" << truth_jet.Eta() << "/" << truth_jet.Phi() << std::endl;
  }

  for (int t=0;t<evt.truthlargeJet().size();t++){
    TLorentzVector truth_ljet = evt.truthlargeJet()[t].mom();
    //std::cout << "truth large jet pt/eta/phi=" << truth_ljet.Perp() << "/" << truth_ljet.Eta() << "/" << truth_ljet.Phi() << std::endl;
  }

  TLorentzVector par;
  //track jet & truth matching
  std::vector<int> tjet_match_id,tjet_truth_index;
  //std::cout << "Number of track jets is " << evt.tjet().size() << std::endl;
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    float pt = evt.tjet()[jidx].mom().Perp()*1e-3;
    float eta = evt.tjet()[jidx].mom().Eta();
    float phi = evt.tjet()[jidx].mom().Phi();
    float max_pt_b=0,max_pt_c=0,max_pt_tau=0,max_pt_light=0;
    int match_id_b=-99,match_id_c=-99,match_id_tau=-99,match_id_light=-99;
    int match_index_b=-1,match_index_c=-1,match_index_tau=-1,match_index_light=-1;
    for(int t=0;t<evt.truth().size();t++){
      par.SetPxPyPzE(evt.truth()[t].px(),evt.truth()[t].py(),evt.truth()[t].pz(),evt.truth()[t].e());
      int id = evt.truth()[t].id();
      if(abs(id)==24 || abs(id)==6) continue;
      float dr = TuDoAtlas::calc_delta_r(eta,par.Eta(),phi,par.Phi());
      if(dr<0.3 && abs(id)==5 && par.Perp()>max_pt_b){
        max_pt_b = par.Perp();
        match_id_b = id;
        match_index_b = t;
      }
      if(dr<0.3 && abs(id)==4 && par.Perp()>max_pt_c){
        max_pt_c = par.Perp();
        match_id_c = id;
        match_index_c = t;
      }
      if(dr<0.3 && abs(id)==15 && par.Perp()>max_pt_tau){
        max_pt_tau = par.Perp();
        match_id_tau = id;
        match_index_tau = t;
      }
      if (dr<0.3 && par.Perp()>max_pt_light && (abs(id)<5 || abs(id)==21)){
        max_pt_light = par.Perp();
        match_id_light = id;
        match_index_light = t;
      }
    }
    if(match_index_b!=-1){
      tjet_match_id.push_back(match_id_b);
      tjet_truth_index.push_back(match_index_b);
    }
    else if(match_index_c!=-1){
      tjet_match_id.push_back(match_id_c);
      tjet_truth_index.push_back(match_index_c);
    }
    else if(match_index_tau!=-1){
      tjet_match_id.push_back(match_id_tau);
      tjet_truth_index.push_back(match_index_tau);
    }
    else if(match_index_light!=-1){
      tjet_match_id.push_back(match_id_light);
      tjet_truth_index.push_back(match_index_light);
    }
    else {
      tjet_match_id.push_back(-99);
      tjet_truth_index.push_back(-1);
    }
  }
  int ptRange[] = {0,100,200,300,500,1000,3000};
  int nPtRange = sizeof(ptRange)/sizeof(int);
  for (size_t jidx = 0; jidx < evt.tjet().size(); ++jidx){
    if(tjet_truth_index[jidx]==-1) continue;
    float pt = evt.tjet()[jidx].mom().Perp()*1e-3;
    float truth_pt = evt.truth()[tjet_truth_index[jidx]].pt()*1e-3;
    float mv2c10 = evt.tjet()[jidx].mv2c10();
    float mv2c10mu = evt.tjet()[jidx].mv2c10mu();
    float mv2c10rnn = evt.tjet()[jidx].mv2c10rnn();
    float dl1_pu = evt.tjet()[jidx].dl1_pu();
    float dl1_pb = evt.tjet()[jidx].dl1_pb();
    float dl1_pc = evt.tjet()[jidx].dl1_pc();
    float dl1_weight = TuDoAtlas::calc_dl1_weight(dl1_pb, dl1_pc, dl1_pu, 0);
    float dl1mu_pu = evt.tjet()[jidx].dl1mu_pu();
    float dl1mu_pb = evt.tjet()[jidx].dl1mu_pb();
    float dl1mu_pc = evt.tjet()[jidx].dl1mu_pc();
    float dl1mu_weight = TuDoAtlas::calc_dl1_weight(dl1mu_pb, dl1mu_pc, dl1mu_pu, 0);
    float dl1rnn_pu = evt.tjet()[jidx].dl1rnn_pu();
    float dl1rnn_pb = evt.tjet()[jidx].dl1rnn_pb();
    float dl1rnn_pc = evt.tjet()[jidx].dl1rnn_pc();
    float dl1rnn_weight = TuDoAtlas::calc_dl1_weight(dl1rnn_pb, dl1rnn_pc, dl1rnn_pu, 1);
    if(tjet_match_id[jidx]==-99) continue;
    if(abs(tjet_match_id[jidx])==5){
      h->h1D("bMatchedTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("bMatchedTrackJetParticlePt", "", suffix)->Fill(truth_pt, weight);
      h->h1D("bMatchedTrackJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      for(int p=0;p<nPtRange-1;p++){ 
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]){
          h->h1D(Form("bMatchedTrackJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10, weight);
          h->h1D(Form("bMatchedTrackJetMv2c10mu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10mu, weight);
          h->h1D(Form("bMatchedTrackJetMv2c10rnn_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10rnn, weight);
          h->h1D(Form("bMatchedTrackJetDl1Pu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pu, weight);
          h->h1D(Form("bMatchedTrackJetDl1Pb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pb, weight);
          h->h1D(Form("bMatchedTrackJetDl1Pc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pc, weight);
          h->h1D(Form("bMatchedTrackJetDl1Weight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_weight, weight);
          h->h1D(Form("bMatchedTrackJetDl1muPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pu, weight);
          h->h1D(Form("bMatchedTrackJetDl1muPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pb, weight);
          h->h1D(Form("bMatchedTrackJetDl1muPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pc, weight);
          h->h1D(Form("bMatchedTrackJetDl1muWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_weight, weight);
          h->h1D(Form("bMatchedTrackJetDl1rnnPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pu, weight);
          h->h1D(Form("bMatchedTrackJetDl1rnnPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pb, weight);
          h->h1D(Form("bMatchedTrackJetDl1rnnPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pc, weight);
          h->h1D(Form("bMatchedTrackJetDl1rnnWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_weight, weight);
        }
      }
      if(mv2c10>0.6455) {
        h->h1D("bMatchedTrackJetPtPassTrackBtag", "", suffix)->Fill(pt, weight);
        h->h1D("bMatchedTrackJetParticlePtPassTrackBtag", "", suffix)->Fill(truth_pt, weight);
      }
    }
    else if(abs(tjet_match_id[jidx])==4){
      h->h1D("cMatchedTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("cMatchedTrackJetParticlePt", "", suffix)->Fill(truth_pt, weight);
      h->h1D("cMatchedTrackJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      for(int p=0;p<nPtRange-1;p++){ 
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]){
          h->h1D(Form("cMatchedTrackJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10, weight);
          h->h1D(Form("cMatchedTrackJetMv2c10mu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10mu, weight);
          h->h1D(Form("cMatchedTrackJetMv2c10rnn_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10rnn, weight);
          h->h1D(Form("cMatchedTrackJetDl1Pu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pu, weight);
          h->h1D(Form("cMatchedTrackJetDl1Pb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pb, weight);
          h->h1D(Form("cMatchedTrackJetDl1Pc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pc, weight);
          h->h1D(Form("cMatchedTrackJetDl1Weight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_weight, weight);
          h->h1D(Form("cMatchedTrackJetDl1muPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pu, weight);
          h->h1D(Form("cMatchedTrackJetDl1muPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pb, weight);
          h->h1D(Form("cMatchedTrackJetDl1muPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pc, weight);
          h->h1D(Form("cMatchedTrackJetDl1muWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_weight, weight);
          h->h1D(Form("cMatchedTrackJetDl1rnnPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pu, weight);
          h->h1D(Form("cMatchedTrackJetDl1rnnPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pb, weight);
          h->h1D(Form("cMatchedTrackJetDl1rnnPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pc, weight);
          h->h1D(Form("cMatchedTrackJetDl1rnnWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_weight, weight);
        }
      }
      if(mv2c10>0.6455) {
        h->h1D("cMatchedTrackJetPtPassTrackBtag", "", suffix)->Fill(pt, weight);
        h->h1D("cMatchedTrackJetParticlePtPassTrackBtag", "", suffix)->Fill(truth_pt, weight);
      }
    }
    else if(abs(tjet_match_id[jidx])<5 || abs(tjet_match_id[jidx]==21)){
      //std::cout << "tjet match_id/truth_index=" << tjet_match_id[jidx] << "/" << tjet_truth_index[jidx] << std::endl;
      h->h1D("lightMatchedTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("lightMatchedTrackJetParticlePt", "", suffix)->Fill(truth_pt, weight);
      h->h1D("lightMatchedTrackJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      for(int p=0;p<nPtRange-1;p++){ 
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]){
          h->h1D(Form("lightMatchedTrackJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10, weight);
          h->h1D(Form("lightMatchedTrackJetMv2c10mu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10mu, weight);
          h->h1D(Form("lightMatchedTrackJetMv2c10rnn_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10rnn, weight);
          h->h1D(Form("lightMatchedTrackJetDl1Pu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pu, weight);
          h->h1D(Form("lightMatchedTrackJetDl1Pb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pb, weight);
          h->h1D(Form("lightMatchedTrackJetDl1Pc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pc, weight);
          h->h1D(Form("lightMatchedTrackJetDl1Weight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_weight, weight);
          h->h1D(Form("lightMatchedTrackJetDl1muPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pu, weight);
          h->h1D(Form("lightMatchedTrackJetDl1muPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pb, weight);
          h->h1D(Form("lightMatchedTrackJetDl1muPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pc, weight);
          h->h1D(Form("lightMatchedTrackJetDl1muWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_weight, weight);
          h->h1D(Form("lightMatchedTrackJetDl1rnnPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pu, weight);
          h->h1D(Form("lightMatchedTrackJetDl1rnnPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pb, weight);
          h->h1D(Form("lightMatchedTrackJetDl1rnnPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pc, weight);
          h->h1D(Form("lightMatchedTrackJetDl1rnnWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_weight, weight);
        }
      }
      if(mv2c10>0.6455) {
        h->h1D("lightMatchedTrackJetPtPassTrackBtag", "", suffix)->Fill(pt, weight);
        h->h1D("lightMatchedTrackJetParticlePtPassTrackBtag", "", suffix)->Fill(truth_pt, weight);
      }
    }
  }

  //VR track jet & truth matching
  std::vector<int> vrtjet_match_id,vrtjet_truth_index;
  for (size_t jidx = 0; jidx < evt.vrtjet().size(); ++jidx){
    float pt = evt.vrtjet()[jidx].mom().Perp()*1e-3;
    float eta = evt.vrtjet()[jidx].mom().Eta();
    float phi = evt.vrtjet()[jidx].mom().Phi();
    float max_pt=0;
    float max_pt_b=0,max_pt_c=0,max_pt_tau=0,max_pt_light=0;
    int match_id_b=-99,match_id_c=-99,match_id_tau=-99,match_id_light=-99;
    int match_index_b=-1,match_index_c=-1,match_index_tau=-1,match_index_light=-1;
    for(int t=0;t<evt.truth().size();t++){
      par.SetPxPyPzE(evt.truth()[t].px(),evt.truth()[t].py(),evt.truth()[t].pz(),evt.truth()[t].e());
      int id = evt.truth()[t].id();
      if(abs(id)==24 || abs(id)==6) continue;
      float dr = TuDoAtlas::calc_delta_r(eta,par.Eta(),phi,par.Phi());
      if(dr<0.3 && abs(id)==5 && par.Perp()>max_pt_b){
        max_pt_b = par.Perp();
        match_id_b = id;
        match_index_b = t;
      }
      if(dr<0.3 && abs(id)==4 && par.Perp()>max_pt_c){
        max_pt_c = par.Perp();
        match_id_c = id;
        match_index_c = t;
      }
      if(dr<0.3 && abs(id)==15 && par.Perp()>max_pt_tau){
        max_pt_tau = par.Perp();
        match_id_tau = id;
        match_index_tau = t;
      }
      if (dr<0.3 && par.Perp()>max_pt_light && (abs(id)<5 || abs(id)==21)){
        max_pt_light = par.Perp();
        match_id_light = id;
        match_index_light = t;
      }
    }
    if(match_index_b!=-1){
      vrtjet_match_id.push_back(match_id_b);
      vrtjet_truth_index.push_back(match_index_b);
    }
    else if(match_index_c!=-1){
      vrtjet_match_id.push_back(match_id_c);
      vrtjet_truth_index.push_back(match_index_c);
    }
    else if(match_index_tau!=-1){
      vrtjet_match_id.push_back(match_id_tau);
      vrtjet_truth_index.push_back(match_index_tau);
    }
    else if(match_index_light!=-1){
      vrtjet_match_id.push_back(match_id_light);
      vrtjet_truth_index.push_back(match_index_light);
    }
    else {
      vrtjet_match_id.push_back(-99);
      vrtjet_truth_index.push_back(-1);
    }
  }
  for (size_t jidx = 0; jidx < evt.vrtjet().size(); ++jidx){
    if(vrtjet_truth_index[jidx]==-1) continue;
    float pt = evt.vrtjet()[jidx].mom().Perp()*1e-3;
    float truth_pt = evt.truth()[vrtjet_truth_index[jidx]].pt()*1e-3;
    //std::cout << "match id is " << vrtjet_match_id[jidx] << std::endl;
    if(vrtjet_match_id[jidx]==-99) continue;
    if(abs(vrtjet_match_id[jidx])==5){
      h->h1D("bMatchedVRTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("bMatchedVRTrackJetParticlePt", "", suffix)->Fill(evt.truth()[vrtjet_truth_index[jidx]].pt()*1e-3, weight);
      h->h1D("bMatchedVRTrackJetMv2c10", "", suffix)->Fill(evt.vrtjet()[jidx].mv2c10(), weight);
      for(int p=0;p<nPtRange-1;p++){ 
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]) {
          h->h1D(Form("bMatchedVRTrackJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(evt.vrtjet()[jidx].mv2c10(), weight);
        }
      }
      if(evt.vrtjet()[jidx].mv2c10()>0.8492) {
        h->h1D("bMatchedVRTrackJetPtPassVRTrackBtag", "", suffix)->Fill(pt, weight);
        h->h1D("bMatchedVRTrackJetParticlePtPassVRTrackBtag", "", suffix)->Fill(evt.truth()[vrtjet_truth_index[jidx]].pt()*1e-3, weight);
      }
    }
    else if(abs(vrtjet_match_id[jidx])==4){
      h->h1D("cMatchedVRTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("cMatchedVRTrackJetParticlePt", "", suffix)->Fill(evt.truth()[vrtjet_truth_index[jidx]].pt()*1e-3, weight);
      h->h1D("cMatchedVRTrackJetMv2c10", "", suffix)->Fill(evt.vrtjet()[jidx].mv2c10(), weight);
      for(int p=0;p<nPtRange-1;p++){ 
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]) {
          h->h1D(Form("cMatchedVRTrackJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(evt.vrtjet()[jidx].mv2c10(), weight);
        }
      }
      if(evt.tjet()[jidx].mv2c10()>0.8492) {
        h->h1D("cMatchedVRTrackJetPtPassVRTrackBtag", "", suffix)->Fill(pt, weight);
        h->h1D("cMatchedVRTrackJetParticlePtPassVRTrackBtag", "", suffix)->Fill(evt.truth()[vrtjet_truth_index[jidx]].pt()*1e-3, weight);
      }
    }
    else if(abs(vrtjet_match_id[jidx])<5 || abs(vrtjet_match_id[jidx]==21)){
      h->h1D("lightMatchedVRTrackJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("lightMatchedVRTrackJetParticlePt", "", suffix)->Fill(evt.truth()[vrtjet_truth_index[jidx]].pt()*1e-3, weight);
      h->h1D("lightMatchedVRTrackJetMv2c10", "", suffix)->Fill(evt.vrtjet()[jidx].mv2c10(), weight);
      for(int p=0;p<nPtRange-1;p++){ 
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]) {
          h->h1D(Form("lightMatchedVRTrackJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(evt.vrtjet()[jidx].mv2c10(), weight);
        }
      }
      if(evt.vrtjet()[jidx].mv2c10()>0.8492) {
        h->h1D("lightMatchedVRTrackJetPtPassVRTrackBtag", "", suffix)->Fill(pt, weight);
        h->h1D("lightMatchedVRTrackJetParticlePtPassVRTrackBtag", "", suffix)->Fill(evt.truth()[vrtjet_truth_index[jidx]].pt()*1e-3, weight);
      }
    }
  }

  //small jet & truth matching
  std::vector<int> sjet_match_id,sjet_truth_index;
  std::vector<float> tjet_mv2c10_in_sjet;
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    float pt = evt.jet()[jidx].mom().Perp()*1e-3;
    float eta = evt.jet()[jidx].mom().Eta();
    float phi = evt.jet()[jidx].mom().Phi();
    float max_pt_b=0,max_pt_c=0,max_pt_tau=0,max_pt_light=0;
    int match_id_b=-99,match_id_c=-99,match_id_tau=-99,match_id_light=-99;
    int match_index_b=-1,match_index_c=-1,match_index_tau=-1,match_index_light=-1;
    //std::cout << "small jet pt/eta/phi=" << pt << "/" << eta << "/" << phi << std::endl;
    for(int t=0;t<evt.truth().size();t++){
      par.SetPxPyPzE(evt.truth()[t].px(),evt.truth()[t].py(),evt.truth()[t].pz(),evt.truth()[t].e());
      int id = evt.truth()[t].id();
      if(abs(id)==24 || abs(id)==6) continue;
      float dr = TuDoAtlas::calc_delta_r(eta,par.Eta(),phi,par.Phi());
      if(dr<0.3 && abs(id)==5 && par.Perp()>max_pt_b){
        max_pt_b = par.Perp();
        match_id_b = id;
        match_index_b = t;
      }
      if(dr<0.3 && abs(id)==4 && par.Perp()>max_pt_c){
        max_pt_c = par.Perp();
        match_id_c = id;
        match_index_c = t;
      }
      if(dr<0.3 && abs(id)==15 && par.Perp()>max_pt_tau){
        max_pt_tau = par.Perp();
        match_id_tau = id;
        match_index_tau = t;
      }
      if (dr<0.3 && par.Perp()>max_pt_light && (abs(id)<5 || abs(id)==21)){
        max_pt_light = par.Perp();
        match_id_light = id;
        match_index_light = t;
      }
    }
    //std::cout << "match_index b/c/tau/light=" << match_index_b << "/" << match_index_c << "/" << match_index_tau << "/" << match_index_light << std::endl;
    if(match_index_b!=-1){
      sjet_match_id.push_back(match_id_b);
      sjet_truth_index.push_back(match_index_b);
    }
    else if(match_index_c!=-1){
      sjet_match_id.push_back(match_id_c);
      sjet_truth_index.push_back(match_index_c);
    }
    else if(match_index_tau!=-1){
      sjet_match_id.push_back(match_id_tau);
      sjet_truth_index.push_back(match_index_tau);
    }
    else if(match_index_light!=-1){
      sjet_match_id.push_back(match_id_light);
      sjet_truth_index.push_back(match_index_light);
    }
    else {
      sjet_match_id.push_back(-99);
      sjet_truth_index.push_back(-1);
    }
    //track jets in smalljets
    float max_mv2c10=-99;
    for (size_t tj = 0; tj < evt.tjet().size(); ++tj){
      float tjet_pt = evt.tjet()[tj].mom().Perp()*1e-3;
      float tjet_eta = evt.tjet()[tj].mom().Eta();
      float tjet_phi = evt.tjet()[tj].mom().Phi();
      float tjet_mv2c10 = evt.tjet()[tj].mv2c10();
      float tjet_dr = TuDoAtlas::calc_delta_r(eta,tjet_eta,phi,tjet_phi);
      if(tjet_dr<0.4 && tjet_mv2c10>max_mv2c10){
        max_mv2c10=tjet_mv2c10;
        //std::cout << "tjet in small jet pt/eta/phi/dr=" << tjet_pt << "/" <<  tjet_eta << "/" << tjet_phi << "/" << tjet_dr << std::endl;
      }
    }
    tjet_mv2c10_in_sjet.push_back(max_mv2c10);
  }
  for (size_t jidx = 0; jidx < evt.jet().size(); ++jidx){
    if(sjet_truth_index[jidx]==-1) continue;
    float pt = evt.jet()[jidx].mom().Perp()*1e-3;
    float truth_pt = evt.truth()[sjet_truth_index[jidx]].pt()*1e-3;
    float mv2c10 = evt.jet()[jidx].mv2c10();
    float mv2c10mu = evt.jet()[jidx].mv2c10mu();
    float mv2c10rnn = evt.jet()[jidx].mv2c10rnn();
    float dl1_pu = evt.jet()[jidx].dl1_pu();
    float dl1_pb = evt.jet()[jidx].dl1_pb();
    float dl1_pc = evt.jet()[jidx].dl1_pc();
    float dl1_weight = TuDoAtlas::calc_dl1_weight(dl1_pb, dl1_pc, dl1_pu, 0);
    float dl1mu_pu = evt.jet()[jidx].dl1mu_pu();
    float dl1mu_pb = evt.jet()[jidx].dl1mu_pb();
    float dl1mu_pc = evt.jet()[jidx].dl1mu_pc();
    float dl1mu_weight = TuDoAtlas::calc_dl1_weight(dl1mu_pb, dl1mu_pc, dl1mu_pu, 0);
    float dl1rnn_pu = evt.jet()[jidx].dl1rnn_pu();
    float dl1rnn_pb = evt.jet()[jidx].dl1rnn_pb();
    float dl1rnn_pc = evt.jet()[jidx].dl1rnn_pc();
    float dl1rnn_weight = TuDoAtlas::calc_dl1_weight(dl1rnn_pb, dl1rnn_pc, dl1rnn_pu, 1);
    //std::cout << "dl1_weight=" << dl1_weight << std::endl;
    //std::cout << "sjet_match_id=" << sjet_match_id[jidx] << std::endl;
    if(sjet_match_id[jidx]==-99) continue;
    if(abs(sjet_match_id[jidx])==5){
      h->h1D("bMatchedSmallJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("bMatchedSmallJetParticlePt", "", suffix)->Fill(truth_pt, weight);
      h->h1D("bMatchedSmallJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      for(int p=0;p<nPtRange-1;p++){ 
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]) {
          h->h1D(Form("bMatchedSmallJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10, weight);
          h->h1D(Form("bMatchedSmallJetMv2c10mu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10mu, weight);
          h->h1D(Form("bMatchedSmallJetMv2c10rnn_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10rnn, weight);
          h->h1D(Form("bMatchedSmallJetDl1Pu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pu, weight);
          h->h1D(Form("bMatchedSmallJetDl1Pb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pb, weight);
          h->h1D(Form("bMatchedSmallJetDl1Pc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pc, weight);
          h->h1D(Form("bMatchedSmallJetDl1Weight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_weight, weight);
          h->h1D(Form("bMatchedSmallJetDl1muPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pu, weight);
          h->h1D(Form("bMatchedSmallJetDl1muPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pb, weight);
          h->h1D(Form("bMatchedSmallJetDl1muPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pc, weight);
          h->h1D(Form("bMatchedSmallJetDl1muWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_weight, weight);
          h->h1D(Form("bMatchedSmallJetDl1rnnPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pu, weight);
          h->h1D(Form("bMatchedSmallJetDl1rnnPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pb, weight);
          h->h1D(Form("bMatchedSmallJetDl1rnnPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pc, weight);
          h->h1D(Form("bMatchedSmallJetDl1rnnWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_weight, weight);
          h->h1D(Form("bMatchedTrackJetInSmallJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(tjet_mv2c10_in_sjet[jidx], weight);
        }
      }
      if(mv2c10>btag_thr_calo70) h->h1D("bMatchedSmallJetPtPassCaloBtagFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10>btag_thr_calo70) h->h1D("bMatchedSmallJetParticlePtPassCaloBtagFix70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10>m_spline->Eval(pt)) h->h1D("bMatchedSmallJetPtPassCaloBtagHyb70", "", suffix)->Fill(pt, weight);
      if(mv2c10>m_spline->Eval(pt)) h->h1D("bMatchedSmallJetParticlePtPassCaloBtagHyb70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10mu>0.869) h->h1D("bMatchedSmallJetPtPassCaloBtagMuFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10mu>0.869) h->h1D("bMatchedSmallJetParticlePtPassCaloBtagMuFix70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10rnn>0.868) h->h1D("bMatchedSmallJetPtPassCaloBtagRnnFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10rnn>0.868) h->h1D("bMatchedSmallJetParticlePtPassCaloBtagRnnFix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1_weight>2.01) h->h1D("bMatchedSmallJetPtPassCaloBtagDl1Fix70", "", suffix)->Fill(pt, weight);
      if(dl1_weight>2.01) h->h1D("bMatchedSmallJetParticlePtPassCaloBtagDl1Fix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1mu_weight>1.83) h->h1D("bMatchedSmallJetPtPassCaloBtagDl1muFix70", "", suffix)->Fill(pt, weight);
      if(dl1mu_weight>1.83) h->h1D("bMatchedSmallJetParticlePtPassCaloBtagDl1muFix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1rnn_weight>2.98) h->h1D("bMatchedSmallJetPtPassCaloBtagDl1rnnFix70", "", suffix)->Fill(pt, weight);
      if(dl1rnn_weight>2.98) h->h1D("bMatchedSmallJetParticlePtPassCaloBtagDl1rnnFix70", "", suffix)->Fill(truth_pt, weight);
    }
    else if(abs(sjet_match_id[jidx])==4){
      h->h1D("cMatchedSmallJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("cMatchedSmallJetParticlePt", "", suffix)->Fill(truth_pt, weight);
      h->h1D("cMatchedSmallJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      for(int p=0;p<nPtRange-1;p++){ 
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]) {
          h->h1D(Form("cMatchedSmallJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10, weight);
          h->h1D(Form("cMatchedSmallJetMv2c10mu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10mu, weight);
          h->h1D(Form("cMatchedSmallJetMv2c10rnn_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10rnn, weight);
          h->h1D(Form("cMatchedSmallJetDl1Pu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pu, weight);
          h->h1D(Form("cMatchedSmallJetDl1Pb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pb, weight);
          h->h1D(Form("cMatchedSmallJetDl1Pc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pc, weight);
          h->h1D(Form("cMatchedSmallJetDl1Weight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_weight, weight);
          h->h1D(Form("cMatchedSmallJetDl1muPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pu, weight);
          h->h1D(Form("cMatchedSmallJetDl1muPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pb, weight);
          h->h1D(Form("cMatchedSmallJetDl1muPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pc, weight);
          h->h1D(Form("cMatchedSmallJetDl1muWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_weight, weight);
          h->h1D(Form("cMatchedSmallJetDl1rnnPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pu, weight);
          h->h1D(Form("cMatchedSmallJetDl1rnnPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pb, weight);
          h->h1D(Form("cMatchedSmallJetDl1rnnPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pc, weight);
          h->h1D(Form("cMatchedSmallJetDl1rnnWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_weight, weight);
          h->h1D(Form("cMatchedTrackJetInSmallJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(tjet_mv2c10_in_sjet[jidx], weight);
        }
      }
      if(mv2c10>btag_thr_calo70) h->h1D("cMatchedSmallJetPtPassCaloBtagFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10>btag_thr_calo70) h->h1D("cMatchedSmallJetParticlePtPassCaloBtagFix70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10>m_spline->Eval(pt)) h->h1D("cMatchedSmallJetPtPassCaloBtagHyb70", "", suffix)->Fill(pt, weight);
      if(mv2c10>m_spline->Eval(pt)) h->h1D("cMatchedSmallJetParticlePtPassCaloBtagHyb70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10mu>0.869) h->h1D("cMatchedSmallJetPtPassCaloBtagMuFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10mu>0.869) h->h1D("cMatchedSmallJetParticlePtPassCaloBtagMuFix70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10rnn>0.868) h->h1D("cMatchedSmallJetPtPassCaloBtagRnnFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10rnn>0.868) h->h1D("cMatchedSmallJetParticlePtPassCaloBtagRnnFix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1_weight>2.01) h->h1D("cMatchedSmallJetPtPassCaloBtagDl1Fix70", "", suffix)->Fill(pt, weight);
      if(dl1_weight>2.01) h->h1D("cMatchedSmallJetParticlePtPassCaloBtagDl1Fix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1mu_weight>1.83) h->h1D("cMatchedSmallJetPtPassCaloBtagDl1muFix70", "", suffix)->Fill(pt, weight);
      if(dl1mu_weight>1.83) h->h1D("cMatchedSmallJetParticlePtPassCaloBtagDl1muFix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1rnn_weight>2.98) h->h1D("cMatchedSmallJetPtPassCaloBtagDl1rnnFix70", "", suffix)->Fill(pt, weight);
      if(dl1rnn_weight>2.98) h->h1D("cMatchedSmallJetParticlePtPassCaloBtagDl1rnnFix70", "", suffix)->Fill(truth_pt, weight);
    }
    else if(abs(sjet_match_id[jidx])<5 || abs(sjet_match_id[jidx]==21)){
      h->h1D("lightMatchedSmallJetPt", "", suffix)->Fill(pt, weight);
      h->h1D("lightMatchedSmallJetParticlePt", "", suffix)->Fill(truth_pt, weight);
      h->h1D("lightMatchedSmallJetMv2c10", "", suffix)->Fill(mv2c10, weight);
      for(int p=0;p<nPtRange-1;p++){
        if(truth_pt>ptRange[p] && truth_pt<ptRange[p+1]){
          h->h1D(Form("lightMatchedSmallJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10, weight);
          h->h1D(Form("lightMatchedSmallJetMv2c10mu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10mu, weight);
          h->h1D(Form("lightMatchedSmallJetMv2c10rnn_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(mv2c10rnn, weight);
          h->h1D(Form("lightMatchedSmallJetDl1Pu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pu, weight);
          h->h1D(Form("lightMatchedSmallJetDl1Pb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pb, weight);
          h->h1D(Form("lightMatchedSmallJetDl1Pc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_pc, weight);
          h->h1D(Form("lightMatchedSmallJetDl1Weight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1_weight, weight);
          h->h1D(Form("lightMatchedSmallJetDl1muPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pu, weight);
          h->h1D(Form("lightMatchedSmallJetDl1muPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pb, weight);
          h->h1D(Form("lightMatchedSmallJetDl1muPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_pc, weight);
          h->h1D(Form("lightMatchedSmallJetDl1muWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1mu_weight, weight);
          h->h1D(Form("lightMatchedSmallJetDl1rnnPu_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pu, weight);
          h->h1D(Form("lightMatchedSmallJetDl1rnnPb_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pb, weight);
          h->h1D(Form("lightMatchedSmallJetDl1rnnPc_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_pc, weight);
          h->h1D(Form("lightMatchedSmallJetDl1rnnWeight_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(dl1rnn_weight, weight);
          h->h1D(Form("lightMatchedTrackJetInSmallJetMv2c10_pt%d_%d",ptRange[p],ptRange[p+1]), "", suffix)->Fill(tjet_mv2c10_in_sjet[jidx], weight);
        }
      }
      if(mv2c10>btag_thr_calo70) h->h1D("lightMatchedSmallJetPtPassCaloBtagFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10>btag_thr_calo70) h->h1D("lightMatchedSmallJetParticlePtPassCaloBtagFix70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10>m_spline->Eval(pt)) h->h1D("lightMatchedSmallJetPtPassCaloBtagHyb70", "", suffix)->Fill(pt, weight);
      if(mv2c10>m_spline->Eval(pt)) h->h1D("lightMatchedSmallJetParticlePtPassCaloBtagHyb70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10mu>0.869) h->h1D("lightMatchedSmallJetPtPassCaloBtagMuFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10mu>0.869) h->h1D("lightMatchedSmallJetParticlePtPassCaloBtagMuFix70", "", suffix)->Fill(truth_pt, weight);
      if(mv2c10rnn>0.868) h->h1D("lightMatchedSmallJetPtPassCaloBtagRnnFix70", "", suffix)->Fill(pt, weight);
      if(mv2c10rnn>0.868) h->h1D("lightMatchedSmallJetParticlePtPassCaloBtagRnnFix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1_weight>2.01) h->h1D("lightMatchedSmallJetPtPassCaloBtagDl1Fix70", "", suffix)->Fill(pt, weight);
      if(dl1_weight>2.01) h->h1D("lightMatchedSmallJetParticlePtPassCaloBtagDl1Fix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1mu_weight>1.83) h->h1D("lightMatchedSmallJetPtPassCaloBtagDl1muFix70", "", suffix)->Fill(pt, weight);
      if(dl1mu_weight>1.83) h->h1D("lightMatchedSmallJetParticlePtPassCaloBtagDl1muFix70", "", suffix)->Fill(truth_pt, weight);
      if(dl1rnn_weight>2.98) h->h1D("lightMatchedSmallJetPtPassCaloBtagDl1rnnFix70", "", suffix)->Fill(pt, weight);
      if(dl1rnn_weight>2.98) h->h1D("lightMatchedSmallJetParticlePtPassCaloBtagDl1rnnFix70", "", suffix)->Fill(truth_pt, weight);
    }
  }

  //large jet & truth matching
  /*float min_dr_htop_truth_ljet=100;
    for (int t=0;t<evt.truthlargeJet().size();t++){
    TLorentzVector truth_ljet = evt.truthlargeJet()[t].mom();
    float dr_htop_truth_ljet = TuDoAtlas::calc_delta_r(truth_ljet.Eta(),htop.Eta(),truth_ljet.Phi(),htop.Phi());
    if(dr_htop_truth_ljet<min_dr_htop_truth_ljet) min_dr_htop_truth_ljet = dr_htop_truth_ljet;
    }
    bool isClose_htop_truth_ljet = min_dr_htop_truth_ljet<0.75;*/
  //std::cout << "min_dr_htop_truth_ljet=" << min_dr_htop_truth_ljet << std::endl;
  std::vector<int> ljet_match_id,ljet_truth_index;
  bool htop_pass_sub80=false,htop_pass_sub50=false;
  bool htop_pass_smooth_mt80=false,htop_pass_smooth_mt50=false;
  bool htop_pass_smooth_ts80=false,htop_pass_smooth_ts50=false;
  bool htop_pass_smooth_qt80=false,htop_pass_smooth_qt50=false;
  bool htop_pass_bdt80=false,htop_pass_dnn80=false;
  int htop_match_id=-1,min_dr_htop=100;
  //std::cout << "Number of large jet is " << evt.largeJet().size() << std::endl;
  for (size_t jidx = 0; jidx < evt.largeJet().size(); ++jidx){
    const TLorentzVector &lj = evt.largeJet()[jidx].mom();
    float pt = lj.Perp();
    float eta = lj.Eta();
    float phi = lj.Phi();
    float mass = lj.M();
    float tau32 = evt.largeJet()[jidx].tau32();
    if(pt*1e-3<300 || fabs(eta)>2) continue;
    //if(pt*1e-3>400) continue;
    float max_pt=0;
    int match_id=-99,match_index=-1;
    //std::cout << "large jet pt/eta/phi/id=" << pt << "/" << eta << "/" << phi << std::endl;
    for(int t=0;t<evt.truth().size();t++){
      //std::cout << "particle px/py/pz=" << evt.truth()[t].px() << "/" << evt.truth()[t].py() << "/" << evt.truth()[t].pz() << std::endl;
      par.SetPxPyPzE(evt.truth()[t].px(),evt.truth()[t].py(),evt.truth()[t].pz(),evt.truth()[t].e());
      int id = evt.truth()[t].id();
      float dr = TuDoAtlas::calc_delta_r(eta,par.Eta(),phi,par.Phi());
      if(dr<1. && par.Perp()>max_pt){
        max_pt = par.Perp();
        match_id = id;
        match_index = t;
        //std::cout << "matched particle id/pt/eta/phi/dr=" << id << "/" << par.Perp() << "/" << par.Eta() << "/" << par.Phi() << "/" << dr << std::endl;
      }
    }
    ljet_match_id.push_back(match_id);
    ljet_truth_index.push_back(match_index);
    //std::cout << "matched id=" << match_id << std::endl; 

    float dr_htop = TuDoAtlas::calc_delta_r(eta,htop.Eta(),phi,htop.Phi());
    float dr_b_from_htop = TuDoAtlas::calc_delta_r(eta,b_from_htop.Eta(),phi,b_from_htop.Phi());
    float dr_Wdecay1_from_htop = TuDoAtlas::calc_delta_r(eta,Wdecay1_from_htop.Eta(),phi,Wdecay1_from_htop.Phi());
    float dr_Wdecay2_from_htop = TuDoAtlas::calc_delta_r(eta,Wdecay2_from_htop.Eta(),phi,Wdecay2_from_htop.Phi());
    float min_dr_reco_truth=100,min_dr_htop_truth_ljet=100;
    for (int t=0;t<evt.truthlargeJet().size();t++){
      TLorentzVector truth_ljet = evt.truthlargeJet()[t].mom();
      float dr_reco_truth = TuDoAtlas::calc_delta_r(truth_ljet.Eta(),eta,truth_ljet.Phi(),phi);
      float dr_htop_truth_ljet = TuDoAtlas::calc_delta_r(truth_ljet.Eta(),htop.Eta(),truth_ljet.Phi(),htop.Phi());
      if(dr_reco_truth<min_dr_reco_truth) min_dr_reco_truth = dr_reco_truth;
      if(dr_htop_truth_ljet<min_dr_htop_truth_ljet) min_dr_htop_truth_ljet = dr_htop_truth_ljet;
    }
    //std::cout << "dr ljet between truth&reco=" << min_dr_reco_truth << std::endl;
    //std::cout << "dr ljet between truth ljet&htop=" << min_dr_htop_truth_ljet << std::endl;
    //std::cout << "dr t/b/Wdecay1/Wdecay2=" << dr_htop << "/" << dr_b_from_htop << "/" << dr_Wdecay1_from_htop << "/" << dr_Wdecay2_from_htop << std::endl;
    //std::cout << "good sub80/smoothmt80=" << evt.largeJet()[jidx].good_sub80() << "/" << evt.largeJet()[jidx].good_smooth_mt80() << std::endl;
    if(dr_b_from_htop<0.75 && dr_Wdecay1_from_htop<0.75 && dr_Wdecay2_from_htop<0.75 && min_dr_reco_truth<0.75 && min_dr_htop_truth_ljet<0.75){
      //if(dr_b_from_htop<0.75 && dr_Wdecay1_from_htop<0.75 && dr_Wdecay2_from_htop<0.75){
      htop_match_id = jidx;
      min_dr_htop=dr_htop;
      if(evt.largeJet()[jidx].good_sub80()) htop_pass_sub80=true; 
      if(evt.largeJet()[jidx].good_sub50()) htop_pass_sub50=true; 
      if(evt.largeJet()[jidx].good_smooth_mt80()) htop_pass_smooth_mt80=true; 
      if(evt.largeJet()[jidx].good_smooth_mt50()) htop_pass_smooth_mt50=true; 
      if(evt.largeJet()[jidx].good_smooth_ts80()) htop_pass_smooth_ts80=true; 
      if(evt.largeJet()[jidx].good_smooth_ts50()) htop_pass_smooth_ts50=true; 
      if(evt.largeJet()[jidx].good_smooth_qt80()) htop_pass_smooth_qt80=true; 
      if(evt.largeJet()[jidx].good_smooth_qt50()) htop_pass_smooth_qt50=true; 
      if(evt.largeJet()[jidx].good_bdt80()) htop_pass_bdt80=true; 
      if(evt.largeJet()[jidx].good_dnn80()) htop_pass_dnn80=true; 
    }
    //std::cout << "good sub80/bdt80/dnn80=" << htop_pass_bdt80 << "/" << htop_pass_dnn80 << std::endl;
    h->h1D("allLargeJetPt", "", suffix)->Fill(pt*1e-3, weight);
    if(evt.largeJet()[jidx].good_sub80()) {
      h->h1D("allLargeJetPtPassSub80", "", suffix)->Fill(pt*1e-3, weight);
      h->h2D("allLargeJetPtVsMassPassSub80", "", suffix)->Fill(pt*1e-3,mass*1e-3);
      h->h2D("allLargeJetPtVsTau32PassSub80", "", suffix)->Fill(pt*1e-3,tau32);
    }
    if(evt.largeJet()[jidx].good_sub50()) h->h1D("allLargeJetPtPassSub50", "", suffix)->Fill(pt*1e-3, weight);
    if(evt.largeJet()[jidx].good_smooth_mt80()) {
      h->h1D("allLargeJetPtPassSmoothMT80", "", suffix)->Fill(pt*1e-3, weight);
      h->h2D("allLargeJetPtVsMassPassSmoothMT80", "", suffix)->Fill(pt*1e-3,mass*1e-3);
      h->h2D("allLargeJetPtVsTau32PassSmoothMT80", "", suffix)->Fill(pt*1e-3,tau32);
    }
    if(evt.largeJet()[jidx].good_smooth_mt50()) h->h1D("allLargeJetPtPassSmoothMT50", "", suffix)->Fill(pt*1e-3, weight);
    if(evt.largeJet()[jidx].good_smooth_ts80()) h->h1D("allLargeJetPtPassSmoothTS80", "", suffix)->Fill(pt*1e-3, weight);
    if(evt.largeJet()[jidx].good_smooth_ts50()) h->h1D("allLargeJetPtPassSmoothTS50", "", suffix)->Fill(pt*1e-3, weight);
    if(evt.largeJet()[jidx].good_smooth_qt80()) h->h1D("allLargeJetPtPassSmoothQT80", "", suffix)->Fill(pt*1e-3, weight);
    if(evt.largeJet()[jidx].good_smooth_qt50()) h->h1D("allLargeJetPtPassSmoothQT50", "", suffix)->Fill(pt*1e-3, weight);
    if(evt.largeJet()[jidx].good_bdt80()) h->h1D("allLargeJetPtPassBDT80", "", suffix)->Fill(pt*1e-3, weight);
    if(evt.largeJet()[jidx].good_dnn80()) h->h1D("allLargeJetPtPassDNN80", "", suffix)->Fill(pt*1e-3, weight);
    }
    //std::cout << "htop pt/eta/phi=" << htop.Perp() << "/" << htop.Eta() << "/" << htop.Phi() << std::endl;
    if(htop_match_id!=-1){
      h->h1D("htopMatchedLargeJetPt", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      h->h1D("htopMatchedLargeJetParticlePt", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      h->h2D("htopMatchedLargeJetPtVsMass", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3,evt.largeJet()[htop_match_id].mom().M()*1e-3);
      h->h2D("htopMatchedLargeJetPtVsTau32", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3,evt.largeJet()[htop_match_id].tau32());
      if(htop_pass_sub80) h->h1D("htopMatchedLargeJetParticlePtPassSub80", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_sub80) h->h1D("htopMatchedLargeJetPtPassSub80", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_sub50) h->h1D("htopMatchedLargeJetParticlePtPassSub50", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_sub50) h->h1D("htopMatchedLargeJetPtPassSub50", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_smooth_mt80) h->h1D("htopMatchedLargeJetParticlePtPassSmoothMT80", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_smooth_mt80) h->h1D("htopMatchedLargeJetPtPassSmoothMT80", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_smooth_mt50) h->h1D("htopMatchedLargeJetParticlePtPassSmoothMT50", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_smooth_mt50) h->h1D("htopMatchedLargeJetPtPassSmoothMT50", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_smooth_ts80) h->h1D("htopMatchedLargeJetParticlePtPassSmoothTS80", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_smooth_ts80) h->h1D("htopMatchedLargeJetPtPassSmoothTS80", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_smooth_ts50) h->h1D("htopMatchedLargeJetParticlePtPassSmoothTS50", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_smooth_ts50) h->h1D("htopMatchedLargeJetPtPassSmoothTS50", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_smooth_qt80) h->h1D("htopMatchedLargeJetParticlePtPassSmoothQT80", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_smooth_qt80) h->h1D("htopMatchedLargeJetPtPassSmoothQT80", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_smooth_qt50) h->h1D("htopMatchedLargeJetParticlePtPassSmoothQT50", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_smooth_qt50) h->h1D("htopMatchedLargeJetPtPassSmoothQT50", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_bdt80) h->h1D("htopMatchedLargeJetParticlePtPassBDT80", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_bdt80) h->h1D("htopMatchedLargeJetPtPassBDT80", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      if(htop_pass_dnn80) h->h1D("htopMatchedLargeJetParticlePtPassDNN80", "", suffix)->Fill(htop.Perp()*1e-3, weight);
      if(htop_pass_dnn80) h->h1D("htopMatchedLargeJetPtPassDNN80", "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().Perp()*1e-3, weight);
      for(int p=0;p<20;p++){
        if(evt.largeJet()[htop_match_id].mom().Perp()*1e-3>100*p && evt.largeJet()[htop_match_id].mom().Perp()*1e-3<100*(p+1)){
          h->h1D(Form("htopMatchedLargeJetMass_pt%d_%d",100*p,100*(p+1)), "", suffix)->Fill(evt.largeJet()[htop_match_id].mom().M()*1e-3, weight);
          h->h1D(Form("htopMatchedLargeJetTau32_pt%d_%d",100*p,100*(p+1)), "", suffix)->Fill(evt.largeJet()[htop_match_id].tau32(), weight);
        }
      }
    }

    for (size_t jidx = 0; jidx < evt.largeJet().size(); ++jidx){
      const TLorentzVector &lj = evt.largeJet()[jidx].mom();
      float pt = lj.Perp();
      float eta = lj.Eta();
      float phi = lj.Phi();
      if(pt*1e-3<300 || fabs(eta)>2) continue;
      if(ljet_match_id[jidx]==-99) continue;
      if(abs(ljet_match_id[jidx])==6){
        h->h1D("topMatchedLargeJetPt", "", suffix)->Fill(pt*1e-3, weight);
        h->h1D("topMatchedLargeJetParticlePt", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        if(evt.largeJet()[jidx].good_sub80()) {
          h->h1D("topMatchedLargeJetPtPassSub80", "", suffix)->Fill(pt*1e-3, weight);
          h->h1D("topMatchedLargeJetParticlePtPassSub80", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        }
        if(evt.largeJet()[jidx].good_sub50()) {
          h->h1D("topMatchedLargeJetPtPassSub50", "", suffix)->Fill(pt*1e-3, weight);
          h->h1D("topMatchedLargeJetParticlePtPassSub50", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        }
        if(evt.largeJet()[jidx].good_smooth_mt80()) {
          h->h1D("topMatchedLargeJetPtPassSmoothMT80", "", suffix)->Fill(pt*1e-3, weight);
          h->h1D("topMatchedLargeJetParticlePtPassSmoothMT80", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        }
      }
      else if(abs(ljet_match_id[jidx])==24){
        h->h1D("wMatchedLargeJetPt", "", suffix)->Fill(pt*1e-3, weight);
        h->h1D("wMatchedLargeJetParticlePt", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        if(evt.largeJet()[jidx].good_sub80()) {
          h->h1D("wMatchedLargeJetPtPassSub80", "", suffix)->Fill(pt*1e-3, weight);
          h->h1D("wMatchedLargeJetParticlePtPassSub80", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        }
        if(evt.largeJet()[jidx].good_sub50()) {
          h->h1D("wMatchedLargeJetPtPassSub50", "", suffix)->Fill(pt*1e-3, weight);
          h->h1D("wMatchedLargeJetParticlePtPassSub50", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        }
      }
      else if(abs(ljet_match_id[jidx])<6 || abs(ljet_match_id[jidx])==21){
        h->h1D("qcdMatchedLargeJetPt", "", suffix)->Fill(pt*1e-3, weight);
        h->h1D("qcdMatchedLargeJetParticlePt", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        if(evt.largeJet()[jidx].good_sub80()) {
          h->h1D("qcdMatchedLargeJetPtPassSub80", "", suffix)->Fill(pt*1e-3, weight);
          h->h1D("qcdMatchedLargeJetParticlePtPassSub80", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        }
        if(evt.largeJet()[jidx].good_sub50()) {
          h->h1D("qcdMatchedLargeJetPtPassSub50", "", suffix)->Fill(pt*1e-3, weight);
          h->h1D("qcdMatchedLargeJetParticlePtPassSub50", "", suffix)->Fill(evt.truth()[ljet_truth_index[jidx]].pt()*1e-3, weight);
        }
      }
    }

  }
