/**
 * @brief Code to implement the Matrix method for the QCD estimation
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef MMUTILS_H
#define MMUTILS_H

#include "TopNtupleAnalysis/Event.h"

#include <string>
#include <getopt.h>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"

class MMUtils{

  public:
    MMUtils(const int isBtagged, const std::string &eff_filename2015, const std::string &fake_filename2015, const std::string &eff_filename2016, const std::string &fake_filename2016); 
    ~MMUtils();

    float getMMweights(const Event &evt, const int runMM_StatErr, const bool isElectron, const bool isBoosted, const unsigned int runNumber);
    
    void get2Drates(float &rate, float &rate_err, TH2F* rate_map, float x, float y);
    void get1Drates(float &rate, float &rate_err, TH1F* rate_map, float x);
    
    void getRatesBoostedMu(float &realRate, float &realRate_err, float &fakeRate, float &fakeRate_err, float lepPt, float closejl_DR, const unsigned int runNumber);
    void getRatesBoostedEl(float &realRate, float &realRate_err, float &fakeRate, float &fakeRate_err, float lepPt, float closejl_DR, float absEta, float cosDPhi, const unsigned int runNumber);
    void getRatesResolvedMu(float &realRate, float &realRate_err, float &fakeRate, float &fakeRate_err, float lepPt, float closejl_DR, float absEta, float cosDPhi, float met, float mwt, float DPhi, float topoetcone, const unsigned int runNumber);
    void getRatesResolvedEl(float &realRate, float &realRate_err, float &fakeRate, float &fakeRate_err, float lepPt, float closejl_DR, float closejl_pT, float cosDPhi, float mwt, float topoetcone, const unsigned int runNumber);
    
  private:

    TH2F * eff_map_resolved_e_2015;
    TH2F * eff_map_resolved_mu_2015;
    TH2F * eff_map_boosted_e_2015;
    TH2F * eff_map_boosted_mu_2015;

    TH2F * eff_map_resolved_mu_2015_lDR;
    TH2F * eff_map_resolved_mu_2015_mDR;
    TH2F * eff_map_resolved_mu_2015_hDR;
	
    TH1F * fake_map_resolved_mu_2015_tmp;
    TH2F * fake_map_resolved_mu_2015_lDR;
    TH2F * fake_map_resolved_mu_2015_mDR;
    TH2F * fake_map_resolved_mu_2015_hDR;
    
    TH2F * fake_map_resolved_e_2015;
    TH2F * fake_map2_resolved_e_2015;
    TH2F * fake_map_resolved_e_2015_lEta;
    TH2F * fake_map_resolved_e_2015_hEta;
    
//     TH1F * fake_minDeltaR_resolved_e_hEta;
//     TH1F * fake_pt_resolved_e;
//     
//     TH1F * fake_pt_hmWt_resolved_e;
//     TH1F * fake_pt_mmWt_resolved_e;
//     TH1F * fake_pt_lmWt_resolved_e;
//     
     TH1F * fake_pt_boosted_e_2015;
     TH1F * fake_dr_boosted_mu_2015;

    TH2F * eff_map_resolved_e_2016;
    TH2F * eff_map_resolved_mu_2016;
    TH2F * eff_map_resolved_mu_2016_lDR;
    TH2F * eff_map_resolved_mu_2016_mDR;
    TH2F * eff_map_resolved_mu_2016_hDR;
    TH2F * eff_map_boosted_e_2016;
    TH2F * eff_map_boosted_mu_2016;
	
    TH2F * fake_map_resolved_e_2016_lDR;
    TH2F * fake_map_resolved_e_2016_mDR;
    TH2F * fake_map_resolved_e_2016_hDR;
    
    TH1F * fake_map_resolved_mu_2016_tmp;
    TH2F * fake_map_resolved_mu_2016_lDR;
    TH2F * fake_map_resolved_mu_2016_mDR;
    TH2F * fake_map_resolved_mu_2016_hDR;
    
    TH2F * fake_map_resolved_e_2016;
    TH2F * fake_map2_resolved_e_2016;
    TH2F * fake_map_resolved_e_2016_lEta;
    TH2F * fake_map_resolved_e_2016_hEta;
    
//     TH1F * fake_minDeltaR_resolved_e_hEta;
//     TH1F * fake_pt_resolved_e;
//     
//     TH1F * fake_pt_hmWt_resolved_e;
//     TH1F * fake_pt_mmWt_resolved_e;
//     TH1F * fake_pt_lmWt_resolved_e;
//     
     TH1F * fake_pt_boosted_e_2016;
     TH1F * fake_dr_boosted_mu_2016;
    
    int m_isBtagged; 
        
};

#endif

