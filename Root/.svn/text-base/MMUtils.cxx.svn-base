/**
 * @brief Code to parse input parameters without BOOST to avoid the unnecessary overload of libraries needed.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch> 
 */
#include "TopNtupleAnalysis/MMUtils.h"

#include <getopt.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"


MMUtils::MMUtils(const int isBtagged, const std::string &eff_filename2015, const std::string &fake_filename2015, const std::string &eff_filename2016, const std::string &fake_filename2016) {
    m_isBtagged = isBtagged;
    
    std::cout << "MMUtils:: Going to use  " << std::endl;
    std::cout << "  for 2015:" << eff_filename2015 << " and " << fake_filename2015 << std::endl;
    std::cout << "  for 2016:" << eff_filename2016 << " and " << fake_filename2016 << std::endl;
    // ------------------------ 2015 files ----------------------------------------------------//
    TFile m_eff_rootfile2015(eff_filename2015.c_str(), "r");
    
    //eff_map_resolved_e_2015 = 	(TH2F*) m_eff_rootfile2015.Get("eff_pTdr_resolved_e");
    //if(eff_map_resolved_e_2015)	eff_map_resolved_e_2015->SetDirectory(0);

    eff_map_resolved_e_2015 = 	(TH2F*) m_eff_rootfile2015.Get("eff_lepPt_topoetcone_resolved_e");
    if(eff_map_resolved_e_2015)	eff_map_resolved_e_2015->SetDirectory(0);
    
    eff_map_resolved_mu_2015 = 	(TH2F*) m_eff_rootfile2015.Get("eff_LepPt_DR_resolved_mu");
    if(eff_map_resolved_mu_2015)	eff_map_resolved_mu_2015->SetDirectory(0);

    eff_map_resolved_mu_2015_lDR =  (TH2F*) m_eff_rootfile2015.Get("eff_lepPt_topoetcone_lowDR_resolved_mu");
    if(eff_map_resolved_mu_2015_lDR)        eff_map_resolved_mu_2015_lDR->SetDirectory(0);

    eff_map_resolved_mu_2015_mDR =  (TH2F*) m_eff_rootfile2015.Get("eff_lepPt_topoetcone_medDR_resolved_mu");
    if(eff_map_resolved_mu_2015_mDR)        eff_map_resolved_mu_2015_mDR->SetDirectory(0);

    eff_map_resolved_mu_2015_hDR =  (TH2F*) m_eff_rootfile2015.Get("eff_lepPt_topoetcone_highDR_resolved_mu");
    if(eff_map_resolved_mu_2015_hDR)        eff_map_resolved_mu_2015_hDR->SetDirectory(0);

          
    eff_map_boosted_e_2015 = 	(TH2F*) m_eff_rootfile2015.Get("eff_pTdr_boosted_e");
    if(eff_map_boosted_e_2015)	eff_map_boosted_e_2015->SetDirectory(0);
    
    eff_map_boosted_mu_2015 = 	(TH2F*) m_eff_rootfile2015.Get("eff_pTdr_boosted_mu");    
    if(eff_map_boosted_mu_2015)	eff_map_boosted_mu_2015->SetDirectory(0);
      
    TFile m_fake_rootfile2015(fake_filename2015.c_str(), "r");

    fake_map_resolved_e_2015 =        (TH2F*)m_fake_rootfile2015.Get("2Dfake_lepPt_topoetcone_resolved_e");
    if(fake_map_resolved_e_2015)      fake_map_resolved_e_2015->SetDirectory(0);

    /*
    fake_map_resolved_mu_2015_lDR =    (TH2F*)m_fake_rootfile2015.Get("2Dfake_mwt_met_map_lowDR_resolved_mu");
    if(fake_map_resolved_mu_2015_lDR)  fake_map_resolved_mu_2015_lDR->SetDirectory(0);
    fake_map_resolved_mu_2015_mDR =    (TH2F*)m_fake_rootfile2015.Get("2Dfake_mwt_met_map_medDR_resolved_mu");
    if(fake_map_resolved_mu_2015_mDR)  fake_map_resolved_mu_2015_mDR->SetDirectory(0);
    fake_map_resolved_mu_2015_hDR =    (TH2F*)m_fake_rootfile2015.Get("2Dfake_mwt_met_map_highDR_resolved_mu");
    if(fake_map_resolved_mu_2015_hDR)  fake_map_resolved_mu_2015_hDR->SetDirectory(0);
    */
    fake_map_resolved_mu_2015_lDR =    (TH2F*)m_fake_rootfile2015.Get("2Dfake_lepPt_topoetcone_lowDR_resolved_mu");
    if(fake_map_resolved_mu_2015_lDR)  fake_map_resolved_mu_2015_lDR->SetDirectory(0);
    fake_map_resolved_mu_2015_mDR =    (TH2F*)m_fake_rootfile2015.Get("2Dfake_lepPt_topoetcone_resolved_mu");
    if(fake_map_resolved_mu_2015_mDR)  fake_map_resolved_mu_2015_mDR->SetDirectory(0);
    fake_map_resolved_mu_2015_hDR =    (TH2F*)m_fake_rootfile2015.Get("2Dfake_lepPt_topoetcone_highDR_resolved_mu");
    if(fake_map_resolved_mu_2015_hDR)  fake_map_resolved_mu_2015_hDR->SetDirectory(0);

       
 
    fake_pt_boosted_e_2015 = 	(TH1F*)m_fake_rootfile2015.Get("fakeRate_pt_boosted_e");
    if(fake_pt_boosted_e_2015)	fake_pt_boosted_e_2015->SetDirectory(0);
    
    fake_dr_boosted_mu_2015  = 	(TH1F*)m_fake_rootfile2015.Get("fakeRate_dr_boosted_mu");
    if(fake_dr_boosted_mu_2015)	fake_dr_boosted_mu_2015->SetDirectory(0); 

    // --------------------- 2016 files ------------------------------------------- //
    TFile m_eff_rootfile2016(eff_filename2016.c_str(), "r");
    
    //eff_map_resolved_e_2016 = 	(TH2F*) m_eff_rootfile2016.Get("eff_pTdr_resolved_e");
    //if(eff_map_resolved_e_2016)	eff_map_resolved_e_2016->SetDirectory(0);

    eff_map_resolved_e_2016 = 	(TH2F*) m_eff_rootfile2016.Get("eff_lepPt_topoetcone_resolved_e");
    if(eff_map_resolved_e_2016)	eff_map_resolved_e_2016->SetDirectory(0);
    
    eff_map_resolved_mu_2016 = 	(TH2F*) m_eff_rootfile2016.Get("eff_LepPt_DR_sbin_resolved_mu");
    if(eff_map_resolved_mu_2016)	eff_map_resolved_mu_2016->SetDirectory(0);

    eff_map_resolved_mu_2016_lDR =  (TH2F*) m_eff_rootfile2016.Get("eff_lepPt_topoetcone_lowDR_resolved_mu");
    if(eff_map_resolved_mu_2016_lDR)        eff_map_resolved_mu_2016_lDR->SetDirectory(0);

    eff_map_resolved_mu_2016_mDR =  (TH2F*) m_eff_rootfile2016.Get("eff_lepPt_topoetcone_medDR_resolved_mu");
    if(eff_map_resolved_mu_2016_mDR)        eff_map_resolved_mu_2016_mDR->SetDirectory(0);
          
    eff_map_resolved_mu_2016_hDR =  (TH2F*) m_eff_rootfile2016.Get("eff_lepPt_topoetcone_highDR_resolved_mu");
    if(eff_map_resolved_mu_2016_hDR)        eff_map_resolved_mu_2016_hDR->SetDirectory(0);

    eff_map_boosted_e_2016 = 	(TH2F*) m_eff_rootfile2016.Get("eff_pTdr_boosted_e");
    if(eff_map_boosted_e_2016)	eff_map_boosted_e_2016->SetDirectory(0);
    
    eff_map_boosted_mu_2016 = 	(TH2F*) m_eff_rootfile2016.Get("eff_pTdr_boosted_mu");    
    if(eff_map_boosted_mu_2016)	eff_map_boosted_mu_2016->SetDirectory(0);
      
    TFile m_fake_rootfile2016(fake_filename2016.c_str(), "r");

    fake_map_resolved_e_2016 =        (TH2F*)m_fake_rootfile2016.Get("2Dfake_lepPt_topoetcone_resolved_e");
    if(fake_map_resolved_e_2016)      fake_map_resolved_e_2016->SetDirectory(0);

    fake_map_resolved_mu_2016_lDR =    (TH2F*)m_fake_rootfile2016.Get("2Dfake_lepPt_topoetcone_lowDR_resolved_mu");
    if(fake_map_resolved_mu_2016_lDR)  fake_map_resolved_mu_2016_lDR->SetDirectory(0);
    fake_map_resolved_mu_2016_mDR =    (TH2F*)m_fake_rootfile2016.Get("2Dfake_lepPt_topoetcone_resolved_mu");
    if(fake_map_resolved_mu_2016_mDR)  fake_map_resolved_mu_2016_mDR->SetDirectory(0);
    fake_map_resolved_mu_2016_hDR =    (TH2F*)m_fake_rootfile2016.Get("2Dfake_lepPt_topoetcone_highDR_resolved_mu");
    if(fake_map_resolved_mu_2016_hDR)  fake_map_resolved_mu_2016_hDR->SetDirectory(0);
    
    fake_pt_boosted_e_2016 = 	(TH1F*)m_fake_rootfile2016.Get("fakeRate_pt_boosted_e");
    if(fake_pt_boosted_e_2016)	fake_pt_boosted_e_2016->SetDirectory(0);
    
    fake_dr_boosted_mu_2016  = 	(TH1F*)m_fake_rootfile2016.Get("fakeRate_dr_boosted_mu");
    if(fake_dr_boosted_mu_2016)	fake_dr_boosted_mu_2016->SetDirectory(0); 
  
}
  
MMUtils::~MMUtils(){ 
  /*
  delete eff_map_resolved_e_2015  ;
  delete eff_map_resolved_mu_2015  ;
  delete eff_map_boosted_e_2015 ;
  delete eff_map_boosted_mu_2015 ;
  
  delete fake_map_resolved_mu_2015_lDR;
  delete fake_map_resolved_mu_2015_mDR;
  delete fake_map_resolved_mu_2015_hDR;

  delete eff_map_resolved_e_2016  ;
  delete eff_map_resolved_mu_2016  ;
  delete eff_map_resolved_mu_2016_lDR  ;
  delete eff_map_resolved_mu_2016_mDR  ;
  delete eff_map_resolved_mu_2016_hDR  ;
  delete eff_map_boosted_e_2016 ;
  delete eff_map_boosted_mu_2016 ;
  
  delete fake_map_resolved_e_2016_lDR;
  delete fake_map_resolved_e_2016_mDR;
  delete fake_map_resolved_e_2016_hDR;
  delete fake_map_resolved_mu_2016_lDR;
  delete fake_map_resolved_mu_2016_mDR;
  delete fake_map_resolved_mu_2016_hDR;
    
  delete fake_pt_boosted_e_2015;  
  delete fake_dr_boosted_mu_2015;
  */
} 

void MMUtils::get2Drates(float &rate, float &rate_err, TH2F* rate_map, float x, float y){
   
   int binx(0);
   int biny(0);
      
   if(x > rate_map->GetXaxis()->GetXmax()){
          binx = rate_map->GetXaxis()->FindBin(0.99*rate_map->GetXaxis()->GetXmax());	  
   }else{
          binx = rate_map->GetXaxis()->FindBin(x);
   }
   
   if(y >= rate_map->GetYaxis()->GetXmax()){	
	  biny = rate_map->GetYaxis()->FindBin(0.99*rate_map->GetYaxis()->GetXmax());
   }else{
	  biny = rate_map->GetYaxis()->FindBin(y);
   }
            
   int bin      = rate_map->GetBin(binx, biny, 0);   
   rate    	= rate_map->GetBinContent(bin);
   rate_err 	= rate_map->GetBinError(bin); 
   
   return;
}

void MMUtils::get1Drates(float &rate, float &rate_err, TH1F* rate_map, float x){
   
   int binx(0);
   
   if(x >= rate_map->GetXaxis()->GetXmax()){
          binx = rate_map->GetXaxis()->FindBin(0.99*rate_map->GetXaxis()->GetXmax());	  
   }else{
          binx = rate_map->GetXaxis()->FindBin(x);
   }
      
   rate    	= rate_map->GetBinContent(binx);
   rate_err 	= rate_map->GetBinError(binx); 
   return;
}

void MMUtils::getRatesBoostedMu(float &realRate, float &realRate_err, float &fakeRate, float &fakeRate_err, float lepPt, float closejl_DR, const unsigned int runNumber){

    if (runNumber<290000){//2015
    	get2Drates(realRate, realRate_err, eff_map_boosted_mu_2015, lepPt, closejl_DR);
    	get1Drates(fakeRate, fakeRate_err, fake_dr_boosted_mu_2015, closejl_DR);
    }
    else{
    	get2Drates(realRate, realRate_err, eff_map_boosted_mu_2016, lepPt, closejl_DR);
    	get1Drates(fakeRate, fakeRate_err, fake_dr_boosted_mu_2016, closejl_DR);
    }
    return;
}

void MMUtils::getRatesBoostedEl(float &realRate, float &realRate_err, float &fakeRate, float &fakeRate_err, float lepPt, float closejl_DR, float absEta, float cosDPhi, const unsigned int runNumber){

    if (runNumber<290000){//2015
    	get2Drates(realRate, realRate_err, eff_map_boosted_e_2015, lepPt, closejl_DR);
   
    	if(absEta > 1.8)	get2Drates(fakeRate, fakeRate_err, fake_map_resolved_e_2015_hEta, lepPt, cosDPhi);
    	else			get2Drates(fakeRate, fakeRate_err, fake_map_resolved_e_2015_lEta, lepPt, cosDPhi);
    }
    else{
    	get2Drates(realRate, realRate_err, eff_map_boosted_e_2016, lepPt, closejl_DR);
   
    	if(absEta > 1.8)	get2Drates(fakeRate, fakeRate_err, fake_map_resolved_e_2016_hEta, lepPt, cosDPhi);
    	else			get2Drates(fakeRate, fakeRate_err, fake_map_resolved_e_2016_lEta, lepPt, cosDPhi);
   }
    return;
}

void MMUtils::getRatesResolvedMu(float &realRate, float &realRate_err, float &fakeRate, float &fakeRate_err, float lepPt, float closejl_DR, float closejl_pT, float cosDPhi, float met, float mwt, float DPhi, float topoetcone, const unsigned int runNumber){

 if (runNumber<290000){//2015

    float met_min(0);
    float met_limit(95);

    float mwt_min(0);
    float mwt_limit(95);

    float lepPt_min(0.);
    float lepPt_limit(100);

    float topoet_min(-5.);
    float topoet_limit(10);

    lepPt_min = std::min(lepPt, lepPt_limit);
    topoet_min = std::min(topoetcone, topoet_limit);


    if(closejl_DR < 0.4) {
    if(lepPt_min > 30. && lepPt_min <= 35. && topoet_min > 3 && topoet_min <=10) topoet_min = 2;
    if(lepPt_min > 35. && lepPt_min <= 40. && topoet_min > 6 && topoet_min <=10) topoet_min = 5;
    }
    if(closejl_DR < 0.6 && closejl_DR > 0.4) {
    if(lepPt_min > 85. && lepPt_min <= 100. && topoet_min > 10 && topoet_min <=30) topoet_min = 7;
    }
   
    if(closejl_DR > 0.6)
    get2Drates(realRate, realRate_err, eff_map_resolved_mu_2015_hDR, lepPt_min, topoet_min);
    else if(closejl_DR < 0.6 && closejl_DR > 0.4)
    get2Drates(realRate, realRate_err, eff_map_resolved_mu_2015_mDR, lepPt_min, topoet_min);
    else
    get2Drates(realRate, realRate_err, eff_map_resolved_mu_2015_lDR, lepPt_min, topoet_min);

    // ---- fake rates 2015 ----//
    lepPt_min = 0.;
    lepPt_limit = 99.;

    topoet_min = -5.;
    topoet_limit = 30;


    if(m_isBtagged==0) {

      lepPt_min = std::min(lepPt, lepPt_limit);
      topoet_min = std::min(topoetcone, topoet_limit);
      if(topoet_min > 1 && topoet_min < 3 && lepPt_min > 50 && lepPt_min <=100)
       topoet_min = 0;

    } // if(m_isBtagged==0)
    else {
      lepPt_min = std::min(lepPt, lepPt_limit);
      topoet_min = std::min(topoetcone, topoet_limit);

    }

    if(m_isBtagged==0)
     get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2015_mDR, lepPt_min, topoet_min);
    else {
    if(closejl_DR < 0.4)
    get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2015_lDR, lepPt_min, topoet_min);
    else
    get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2015_hDR, lepPt_min, topoet_min);
    }

    
    if (fakeRate<0) fakeRate=0;


  }
  else{//2016
    float met_min(0);
    float met_limit(95);
    
    float mwt_min(0);
    float mwt_limit(95);

    float lepPt_min(0.);
    float lepPt_limit(600);

    float topoet_min(-5.);
    float topoet_limit(30);

    if(closejl_DR > 0.6)
    get2Drates(realRate, realRate_err, eff_map_resolved_mu_2016_hDR, lepPt, topoetcone); 
    else if(closejl_DR < 0.6 && closejl_DR > 0.4)
    get2Drates(realRate, realRate_err, eff_map_resolved_mu_2016_mDR, lepPt, topoetcone);
    else
    get2Drates(realRate, realRate_err, eff_map_resolved_mu_2016_lDR, lepPt, topoetcone);
    //get2Drates(realRate, realRate_err, eff_map_resolved_mu_2016, lepPt, closejl_DR);
    
  
    if(m_isBtagged==0) {

      lepPt_min = std::min(lepPt, lepPt_limit);
      topoet_min = std::min(topoetcone, topoet_limit);
      if(topoet_min < 1 && lepPt_min >=100)
       topoet_min = 2;

    } // if(m_isBtagged==0)
    else {

      //if(closejl_DR >= 0.4) {

       lepPt_min = std::min(lepPt, (float) 99.);
       topoet_min = std::min(topoetcone, topoet_limit);
   

        if(closejl_DR >= 0.4) {
        if(topoet_min > 10 && lepPt_min >=50 && lepPt_min <100)
        topoet_min = 8;
        } // if(closejl_DR >= 0.4)

    } // if(m_isBtagged==0)




    if(m_isBtagged==0)
     get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2016_mDR, lepPt_min, topoet_min);
    else {
    if(closejl_DR < 0.4)
    get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2016_lDR, lepPt_min, topoet_min);
    else
    get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2016_hDR, lepPt_min, topoet_min);
    }
    if (fakeRate<0) fakeRate=0; 

    if(fakeRate>realRate) {
    std::cout<<"REACHED HERE = "<<std::endl;
    if(closejl_DR < 0.4) {
    int lepPtBin = fake_map_resolved_mu_2016_lDR->GetXaxis()->FindBin(lepPt_min);
    int topoetBin = fake_map_resolved_mu_2016_lDR->GetYaxis()->FindBin(topoet_min);
    std::cout<<"DR < 0.4"<<std::endl;
    std::cout<<"FakeRate = "<<fakeRate<<", RealRate = "<<realRate<<", xbin = "<<lepPtBin<<", ybin = "<<topoetBin<<std::endl;
    } // if(closejl_DR < 0.4)


    if(closejl_DR > 0.4) {
    int lepPtBin = fake_map_resolved_mu_2016_hDR->GetXaxis()->FindBin(lepPt_min);
    int topoetBin = fake_map_resolved_mu_2016_hDR->GetYaxis()->FindBin(topoet_min);
    std::cout<<"DR > 0.4"<<std::endl; 
    std::cout<<"FakeRate = "<<fakeRate<<", RealRate = "<<realRate<<", xbin = "<<lepPtBin<<", ybin = "<<topoetBin<<std::endl;
    } // if(closejl_DR > 0.4)
    }

    //if (fakeRate>1) fakeRate=0.9;
    /*
    if(closejl_DR > 0.6){

      mwt_limit = 195;
      met_limit=140;
      
      if(met>140)mwt_limit=40;
      
      mwt_min = std::min(mwt, mwt_limit);
      met_min = std::min(met, met_limit);
      
      get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2016_hDR, met_min, mwt_min);
      
      if (fakeRate<0) fakeRate=0;
      
    }else if(closejl_DR > 0.4){
      if(m_isBtagged==0){
         if(met<40)	mwt_limit = 95;
         else 		mwt_limit = 195;
      } 
      else {
      		mwt_limit = 195;
      		met_limit = 195;
      }
      mwt_min = std::min(mwt, mwt_limit);
      met_min = std::min(met, met_limit);
      
      get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2016_mDR, met_min, mwt_min);
      
      if (fakeRate<0) fakeRate=0;
            
    }else{
    
      	mwt_limit = 195;
        met_limit = 140;
    
      mwt_min = std::min(mwt, mwt_limit);
      met_min = std::min(met, met_limit);
      
      get2Drates(fakeRate, fakeRate_err, fake_map_resolved_mu_2016_lDR, met_min, mwt_min);
 
      if (fakeRate<0) fakeRate=0;
            
    }//if

    */
  }
  
  return;
        
}

void MMUtils::getRatesResolvedEl(float &realRate, float &realRate_err, float &fakeRate, float &fakeRate_err, float lepPt, float closejl_DR, float absEta, float cosDPhi, float mwt, float topoetcone, const unsigned int runNumber){

 if (runNumber<290000){//2015
    get2Drates(realRate, realRate_err, eff_map_resolved_e_2015, lepPt, topoetcone);  

    float lepPt_min(0.);
    //float lepPt_limit(700);
    float lepPt_limit(695);

    float topoet_min(-5.);
    float topoet_limit(30);
 
    lepPt_min = std::min(lepPt, lepPt_limit);
    topoet_min = std::min(topoetcone, topoet_limit);

    if(m_isBtagged==0) {
    if(lepPt_min > 120 && lepPt_min < 150 && topoet_min > -8. && topoet_min < 1.)
    topoet_min = 2.;
    if(lepPt_min > 50 && lepPt_min < 60 && topoet_min > 6. && topoet_min < 10.)
    lepPt_min = 45.;
    }
    else { // if b-tagged == 1
    if(lepPt_min > 120  && topoet_min > -8. && topoet_min < 1.)
    topoet_min = 2.;
    
    } 

    get2Drates(fakeRate, fakeRate_err, fake_map_resolved_e_2015, lepPt_min, topoet_min);
 }
 else {
    get2Drates(realRate, realRate_err, eff_map_resolved_e_2016, lepPt, topoetcone);  

    float lepPt_min(0.);
    //float lepPt_limit(700);
    float lepPt_limit(695);

    float topoet_min(-5.);
    float topoet_limit(30);

    lepPt_min = std::min(lepPt, lepPt_limit);
    topoet_min = std::min(topoetcone, topoet_limit);

    /*
    if(m_isBtagged==0) {
    if(lepPt_min > 40 && lepPt_min < 50 && topoet_min > 10. && topoet_min < 30.)
    topoet_min = 8.;
    } */
 
    get2Drates(fakeRate, fakeRate_err, fake_map_resolved_e_2016, lepPt_min, topoet_min);
 } 
    return;
}

float MMUtils::getMMweights(const Event &evt, const int runMM_StatErr, const bool isElectron, const bool isBoosted, const unsigned int runNumber) {
   
   bool isTight;
   float d0sig;
   float topoetcone;

   if ( isElectron && evt.electron().size()<1 ) return 0;
   if (!isElectron && evt.muon().size()<1 )     return 0;
   
   TLorentzVector lepP4; 
    
   if (isElectron)    {
       lepP4 = evt.electron()[0].mom();
       isTight = evt.electron()[0].isTightPP();
       d0sig = evt.electron()[0].sd0();
       topoetcone = evt.electron()[0].topoetcone20()*1e-3;
       
   }else {
       lepP4 = evt.muon()[0].mom();
       isTight = evt.muon()[0].isTight();
       d0sig = evt.muon()[0].sd0();
       topoetcone = evt.muon()[0].topoetcone20()*1e-3;
   }
   
   float lepPt  = lepP4.Perp()*1e-3; 
   float absEta = fabs(lepP4.Eta());
   
   //min DR(lept, jet)
   float closejl_pT(0); 
   float closejl_DR(99);
   
   float deltaRapidity2 = 99;
   float deltaPhi2	= 99;
   float deltaR_tmp     = 99;

   size_t jet_idx = 0;
   for (; jet_idx < evt.jet().size(); ++jet_idx){  

       deltaRapidity2 = pow(evt.jet()[jet_idx].mom().Rapidity() - lepP4.Rapidity(), 2);    
       deltaPhi2 = pow(evt.jet()[jet_idx].mom().DeltaPhi(lepP4), 2);	     
       deltaR_tmp = sqrt(deltaPhi2 + deltaRapidity2);

       if (deltaR_tmp < closejl_DR){
   	  closejl_DR = deltaR_tmp;
          closejl_pT = evt.jet()[jet_idx].mom().Perp()*1e-3;
       }  
   }//for 
   
   //Cos(deltaPhi(met,lep)) 
   float deltaPhi = evt.met().DeltaPhi(lepP4);   
   float cosDPhi = std::cos(deltaPhi);   
   float MET = evt.met().Pt()*1e-3;
   float mWt = sqrt(2. * lepP4.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(lepP4)) ))*1e-3; 
        
   //Getting rates
   float fakeRate(0);
   float fakeRate_err(0);
   float realRate(0);
   float realRate_err(0);  
   
//    if (isBoosted){
//    	if (isElectron)	getRatesBoostedEl(realRate, realRate_err, fakeRate, fakeRate_err, lepPt, closejl_DR, absEta, cosDPhi, runNumber);
//    	else		getRatesBoostedMu(realRate, realRate_err, fakeRate, fakeRate_err, lepPt, closejl_DR, runNumber);
//    }
//    else{
   	if (isElectron)	getRatesResolvedEl(realRate, realRate_err, fakeRate, fakeRate_err, lepPt, closejl_DR, absEta, cosDPhi, mWt, topoetcone, runNumber);
   	else		getRatesResolvedMu(realRate, realRate_err, fakeRate, fakeRate_err, lepPt, closejl_DR, closejl_pT, cosDPhi, MET, mWt, deltaPhi, topoetcone, runNumber);
// 		
//    }//isBoosted	

   //Implementing weights
   float Weight = 1;

   //Maximazing the stat error for runMM_StatErr==1 || runMM_StatErr==2
   if (runMM_StatErr==1){
     realRate += realRate_err;
     fakeRate -= fakeRate_err;
   } else if (runMM_StatErr==2){
     realRate -= realRate_err;
     fakeRate += fakeRate_err;
   }
   	
   if (isTight){
     if(realRate>0 && fakeRate>0)     Weight = fakeRate*(realRate - 1)/(realRate - fakeRate);      
     else{
     	   //  if(realRate==0)std::cerr << "Error: realRate equal to " << realRate << " (tight) " << mWt << std::endl;	   
     	   //  if(fakeRate==0)std::cerr << "Error: fakeRate equal to " << fakeRate << " (tight) " << lepPt << " " << topoetcone << std::endl;	    
	     Weight = 0;
     }
   }
   else {	
     if(realRate>0 && fakeRate>0)  Weight = fakeRate*realRate/(realRate - fakeRate);
     else{
	  //   if(realRate==0)std::cerr << "Error: realRate equal to " << realRate << "  (anti-tight) " << mWt << std::endl;	   
     	     //if(fakeRate==0)std::cerr << "Error: fakeRate equal to " << fakeRate << "  (anti-tight) " << mWt << " " << cosDPhi << std::endl;
     	     Weight = 0;
     }       
   
   }//isTight
   
   //if(!isTight)	std::cout << "is tight? " << isTight << " - QCD weight" << Weight << std::endl;
   
   return Weight;  

}//getMMweights

