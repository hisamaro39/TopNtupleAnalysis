/** 
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresQCD.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"

AnaTtresQCD::AnaTtresQCD(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList, int dsid)
  : Analysis(filename, systList), m_electron(electron), m_boosted(boosted), m_dsid(dsid),
    m_neutrinoBuilder("MeV"),  m_chi2("MeV") {
  
  m_chi2.Init(TtresChi2::DATA2015_MC15C);

  std::string suffix(""), btag("");
  IniHistograms(suffix, btag);
  std::string suffix0(""), btag0("btag0_");
  IniHistograms(suffix0, btag0);
  std::string suffix1(""), btag1("btag1_");
  IniHistograms(suffix1, btag1);


  /*
  suffix = "corr1.0_";
  IniHistograms(suffix);
  
  suffix = "corr2.0_";
  IniHistograms(suffix);
  
  suffix = "corr2.1_";
  IniHistograms(suffix);
*/
}//AnaTtresQCD::AnaTtresQCD

AnaTtresQCD::~AnaTtresQCD() {
}

void AnaTtresQCD::run(const Event &evt, double weight, const std::string &suffix){
}


//**************************
//******* Real Rates *******
//**************************



void AnaTtresQCD::runRealRateWQCDCR_2015(const Event &evt, double weight, const std::string &suffix){
  // check channel
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted)
    if (!(evt.passes("bejetsWCR_2015") || evt.passes("bmujetsWCR_2015")))
      return;
  
  if (!m_boosted)
    if (!(evt.passes("rejetsWCR_2015") || evt.passes("rmujetsWCR_2015")))
      return;

  if (!m_boosted)	if(evt.jet().size()<2)	return; 
  
  HistogramService *h = &m_hSvc;
	  
  //Pre-selection: 
  int ttbar_type = evt.MC_ttbar_type();
 
  ///Electrons  
  if (evt.MC_w1l_pdgId()==11){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCe_pt_eta", "", suffix)->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 

  }else if(evt.MC_w2l_pdgId()==-11){  

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCe_pt_eta", "", suffix)->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if
    
  ///Muons 
  if (evt.MC_w1l_pdgId()==13){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCmu_pt_eta", "", suffix) ->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 
  		
  }else if(evt.MC_w2l_pdgId()==-13){

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCmu_pt_eta", "", suffix) ->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if
   
  ///Taus 
  if (evt.MC_w1l_pdgId()==15){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCtau_pt_eta", "", suffix) ->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 
  		
  }else if(evt.MC_w2l_pdgId()==-15){

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCtau_pt_eta", "", suffix) ->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if


  ///----------------------------------
  //Matching truth and reco lepton
  ///----------------------------------

  TLorentzVector lept;    
  int leptMa_pdgId = 0;

  float dr = 99;
  float lep_drMax = 0.2;
  float tau_drMax = 0.4;

  float d0sig(0); 
  
  bool isTight(0);

  if (m_electron) {  
    lept = evt.electron()[0].mom();	  
    d0sig = evt.electron()[0].sd0();

    if (evt.MC_w1l_pdgId()==11){
  	dr = lept.DeltaR(evt.MC_w1l());
  	if (dr<lep_drMax)	leptMa_pdgId = evt.MC_w1l_pdgId();
    
    }else if (evt.MC_w1l_pdgId()==15){
  	dr = lept.DeltaR(evt.MC_w1l());
  	if (dr<tau_drMax)	leptMa_pdgId = evt.MC_w1l_pdgId();

    }else if (evt.MC_w2l_pdgId()==-11){
  	dr = lept.DeltaR(evt.MC_w2l());
  	if (dr<lep_drMax)	leptMa_pdgId = evt.MC_w2l_pdgId();

    }else if (evt.MC_w2l_pdgId()==-15){
  	dr = lept.DeltaR(evt.MC_w2l());
  	if (dr<tau_drMax)	leptMa_pdgId = evt.MC_w2l_pdgId();


    }else if (abs(evt.MC_w1l_pdgId())==13 || abs(evt.MC_w2l_pdgId())==13)	std::cout << "reco electron and truth muon" << std::endl;
  	   
  } else {
    lept = evt.muon()[0].mom();   
    d0sig = evt.muon()[0].sd0();
    
    if (evt.MC_w1l_pdgId()==13 && lept.DeltaR(evt.MC_w1l())<lep_drMax){
  	leptMa_pdgId = evt.MC_w1l_pdgId();
    }else if (evt.MC_w2l_pdgId()==-13 && lept.DeltaR(evt.MC_w2l())<lep_drMax){
  	leptMa_pdgId = evt.MC_w2l_pdgId();
    }
    else if (evt.MC_w1l_pdgId()==15 && lept.DeltaR(evt.MC_w1l())<tau_drMax){
  	leptMa_pdgId = evt.MC_w1l_pdgId();
    }else if (evt.MC_w2l_pdgId()==-15 && lept.DeltaR(evt.MC_w2l())<tau_drMax){
  	leptMa_pdgId = evt.MC_w2l_pdgId();
    }
    else if (abs(evt.MC_w1l_pdgId())==11 || abs(evt.MC_w2l_pdgId())==11)	std::cout << "reco muon and truth electron" << std::endl;    


  }//m_electron

  h->h1D("eff_d0sig", "", suffix)  ->Fill(d0sig, weight);
  
  //nB-tagged jets 
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       	if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          	nTrkBtagged += 1;
  }
  
  //if (nTrkBtagged!=0)	return;			
			
  if (leptMa_pdgId!=0)	{

        GetRealHistograms(evt, weight, suffix, "");	
   }
}//AnaTtresQCD::runRealRateWQCDCR_2015


void AnaTtresQCD::runRealRateQCDCR_2015(const Event &evt, double weight, const std::string &suffix){
  // check channel
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted)
    if (!(evt.passes("bejets_2015") || evt.passes("bmujets_2015")))
      return;
  //std::cout<<"HERE1"<<std::endl;
  if (!m_boosted)
    if (!(evt.passes("rejets_2015") || evt.passes("rmujets_2015")))
      return;

  if (!m_boosted)	if(evt.jet().size()<4)	return; 
  
  HistogramService *h = &m_hSvc;
	  
  //Pre-selection: 
  int ttbar_type = evt.MC_ttbar_type();
 
  ///Electrons  
  if (evt.MC_w1l_pdgId()==11){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCe_pt_eta", "", suffix)->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 

  }else if(evt.MC_w2l_pdgId()==-11){  

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCe_pt_eta", "", suffix)->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if
    
  ///Muons 
  if (evt.MC_w1l_pdgId()==13){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCmu_pt_eta", "", suffix) ->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 
  		
  }else if(evt.MC_w2l_pdgId()==-13){

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCmu_pt_eta", "", suffix) ->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if
   
  ///Taus 
  if (evt.MC_w1l_pdgId()==15){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCtau_pt_eta", "", suffix) ->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 
  		
  }else if(evt.MC_w2l_pdgId()==-15){

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCtau_pt_eta", "", suffix) ->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if


  ///----------------------------------
  //Matching truth and reco lepton
  ///----------------------------------

  TLorentzVector lept;    
  int leptMa_pdgId = 0;

  float dr = 99;
  float lep_drMax = 0.2;
  float tau_drMax = 0.4;

  float d0sig(0); 
  
  bool isTight(0);

  if (m_electron) {  
    lept = evt.electron()[0].mom();	  
    d0sig = evt.electron()[0].sd0();

    if (evt.MC_w1l_pdgId()==11){
  	dr = lept.DeltaR(evt.MC_w1l());
  	if (dr<lep_drMax)	leptMa_pdgId = evt.MC_w1l_pdgId();
    
    }else if (evt.MC_w1l_pdgId()==15){
  	dr = lept.DeltaR(evt.MC_w1l());
  	if (dr<tau_drMax)	leptMa_pdgId = evt.MC_w1l_pdgId();

    }else if (evt.MC_w2l_pdgId()==-11){
  	dr = lept.DeltaR(evt.MC_w2l());
  	if (dr<lep_drMax)	leptMa_pdgId = evt.MC_w2l_pdgId();

    }else if (evt.MC_w2l_pdgId()==-15){
  	dr = lept.DeltaR(evt.MC_w2l());
  	if (dr<tau_drMax)	leptMa_pdgId = evt.MC_w2l_pdgId();


    }else if (abs(evt.MC_w1l_pdgId())==13 || abs(evt.MC_w2l_pdgId())==13)	std::cout << "reco electron and truth muon" << std::endl;
  	   
  } else {
    lept = evt.muon()[0].mom();   
    d0sig = evt.muon()[0].sd0();
    
    if (evt.MC_w1l_pdgId()==13 && lept.DeltaR(evt.MC_w1l())<lep_drMax){
  	leptMa_pdgId = evt.MC_w1l_pdgId();
    }else if (evt.MC_w2l_pdgId()==-13 && lept.DeltaR(evt.MC_w2l())<lep_drMax){
  	leptMa_pdgId = evt.MC_w2l_pdgId();
    }
    else if (evt.MC_w1l_pdgId()==15 && lept.DeltaR(evt.MC_w1l())<tau_drMax){
  	leptMa_pdgId = evt.MC_w1l_pdgId();
    }else if (evt.MC_w2l_pdgId()==-15 && lept.DeltaR(evt.MC_w2l())<tau_drMax){
  	leptMa_pdgId = evt.MC_w2l_pdgId();
    }
    else if (abs(evt.MC_w1l_pdgId())==11 || abs(evt.MC_w2l_pdgId())==11)	std::cout << "reco muon and truth electron" << std::endl;    


  }//m_electron

  h->h1D("eff_d0sig", "", suffix)  ->Fill(d0sig, weight);
  
  //nB-tagged jets 
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       	if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          	nTrkBtagged += 1;
  }
  
  //if (nTrkBtagged!=0)	return;			
		
  if (leptMa_pdgId!=0)	{

     GetRealHistograms(evt, weight, suffix, "");
     if(nTrkBtagged == 0)
     GetRealHistograms(evt, weight, suffix, "btag0_");
     if(nTrkBtagged >= 1)
     GetRealHistograms(evt, weight, suffix, "btag1_");
	
   }
	
  
}//AnaTtresQCD::runRealRateQCDCR_2015




void AnaTtresQCD::runRealRateWQCDCR_2016(const Event &evt, double weight, const std::string &suffix){
  // check channel
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;
  //std::cout<<"HERE"<<std::endl;
  if (m_boosted)
    if (!(evt.passes("bejetsWCR_2016") || evt.passes("bmujetsWCR_2016")))
      return;
  //std::cout<<"HERE1"<<std::endl;
  if (!m_boosted)
    if (!(evt.passes("rejetsWCR_2016") || evt.passes("rmujetsWCR_2016")))
      return;
  if (!m_boosted)	if(evt.jet().size()<2)	return;
 
  HistogramService *h = &m_hSvc;
  //Pre-selection: 
  int ttbar_type = evt.MC_ttbar_type();
 
  ///Electrons  
  if (evt.MC_w1l_pdgId()==11){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCe_pt_eta", "", suffix)->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 

  }else if(evt.MC_w2l_pdgId()==-11){  

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCe_pt_eta", "", suffix)->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if
    
  ///Muons 
  if (evt.MC_w1l_pdgId()==13){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCmu_pt_eta", "", suffix) ->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 
  		
  }else if(evt.MC_w2l_pdgId()==-13){

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCmu_pt_eta", "", suffix) ->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if
   
  ///Taus 
  if (evt.MC_w1l_pdgId()==15){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCtau_pt_eta", "", suffix) ->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 
  		
  }else if(evt.MC_w2l_pdgId()==-15){

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCtau_pt_eta", "", suffix) ->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if


  ///----------------------------------
  //Matching truth and reco lepton
  ///----------------------------------

  TLorentzVector lept;    
  int leptMa_pdgId = 0;

  float dr = 99;
  float lep_drMax = 0.2;
  float tau_drMax = 0.4;

  float d0sig(0); 
  
  bool isTight(0);

  if (m_electron) {  
    lept = evt.electron()[0].mom();	  
    d0sig = evt.electron()[0].sd0();

    if (evt.MC_w1l_pdgId()==11){
  	dr = lept.DeltaR(evt.MC_w1l());
  	if (dr<lep_drMax)	leptMa_pdgId = evt.MC_w1l_pdgId();
    
    }else if (evt.MC_w1l_pdgId()==15){
  	dr = lept.DeltaR(evt.MC_w1l());
  	if (dr<tau_drMax)	leptMa_pdgId = evt.MC_w1l_pdgId();

    }else if (evt.MC_w2l_pdgId()==-11){
  	dr = lept.DeltaR(evt.MC_w2l());
  	if (dr<lep_drMax)	leptMa_pdgId = evt.MC_w2l_pdgId();

    }else if (evt.MC_w2l_pdgId()==-15){
  	dr = lept.DeltaR(evt.MC_w2l());
  	if (dr<tau_drMax)	leptMa_pdgId = evt.MC_w2l_pdgId();


    }else if (abs(evt.MC_w1l_pdgId())==13 || abs(evt.MC_w2l_pdgId())==13)	std::cout << "reco electron and truth muon" << std::endl;
  	   
  } else {
    lept = evt.muon()[0].mom();   
    d0sig = evt.muon()[0].sd0();
 
    if (evt.MC_w1l_pdgId()==13 && lept.DeltaR(evt.MC_w1l())<lep_drMax){
  	leptMa_pdgId = evt.MC_w1l_pdgId();
    }else if (evt.MC_w2l_pdgId()==-13 && lept.DeltaR(evt.MC_w2l())<lep_drMax){
  	leptMa_pdgId = evt.MC_w2l_pdgId();
    }
    else if (evt.MC_w1l_pdgId()==15 && lept.DeltaR(evt.MC_w1l())<tau_drMax){
  	leptMa_pdgId = evt.MC_w1l_pdgId();
    }else if (evt.MC_w2l_pdgId()==-15 && lept.DeltaR(evt.MC_w2l())<tau_drMax){
  	leptMa_pdgId = evt.MC_w2l_pdgId();
    }
    else if (abs(evt.MC_w1l_pdgId())==11 || abs(evt.MC_w2l_pdgId())==11)	std::cout << "reco muon and truth electron" << std::endl;    


  }//m_electron

  h->h1D("eff_d0sig", "", suffix)  ->Fill(d0sig, weight);
  
  //nB-tagged jets 
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       	if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          	nTrkBtagged += 1;
  }
  
//  if (nTrkBtagged!=0)	return;	
//   if (nTrkBtagged < 1) return;		

			
  if (leptMa_pdgId!=0)	{

     GetRealHistograms(evt, weight, suffix, "");
     if(nTrkBtagged == 0)
     GetRealHistograms(evt, weight, suffix, "btag0_");
     if(nTrkBtagged >= 1)
     GetRealHistograms(evt, weight, suffix, "btag1_");
	
   }
}//AnaTtresQCD::runRealRateWQCDCR_2016


void AnaTtresQCD::runRealRateQCDCR_2016(const Event &evt, double weight, const std::string &suffix){
  // check channel
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;
  //std::cout<<"HERE"<<std::endl;
  if (m_boosted)
    if (!(evt.passes("bejets_2016") || evt.passes("bmujets_2016")))
      return;
  //std::cout<<"HERE1"<<std::endl;
  if (!m_boosted)
    if (!(evt.passes("rejets_2016") || evt.passes("rmujets_2016")))
      return;
  if (!m_boosted)	if(evt.jet().size()<4)	return;
 
  HistogramService *h = &m_hSvc;
  //Pre-selection: 
  int ttbar_type = evt.MC_ttbar_type();
 
  ///Electrons  
  if (evt.MC_w1l_pdgId()==11){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCe_pt_eta", "", suffix)->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 

  }else if(evt.MC_w2l_pdgId()==-11){  

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCe_pt_eta", "", suffix)->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if
    
  ///Muons 
  if (evt.MC_w1l_pdgId()==13){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCmu_pt_eta", "", suffix) ->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 
  		
  }else if(evt.MC_w2l_pdgId()==-13){

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCmu_pt_eta", "", suffix) ->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if
   
  ///Taus 
  if (evt.MC_w1l_pdgId()==15){

  	if (evt.MC_w1l().Perp()<15000 || fabs(evt.MC_w1l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCtau_pt_eta", "", suffix) ->Fill(evt.MC_w1l().Perp()*1e-3, evt.MC_w1l().Eta()); 
  		
  }else if(evt.MC_w2l_pdgId()==-15){

  	if (evt.MC_w2l().Perp()<15000 || fabs(evt.MC_w2l().Eta())>3.0)  return; 	
  	h->h2D("eff_MCtau_pt_eta", "", suffix) ->Fill(evt.MC_w2l().Perp()*1e-3, evt.MC_w2l().Eta()); 

  }//if


  ///----------------------------------
  //Matching truth and reco lepton
  ///----------------------------------

  TLorentzVector lept;    
  int leptMa_pdgId = 0;

  float dr = 99;
  float lep_drMax = 0.2;
  float tau_drMax = 0.4;

  float d0sig(0); 
  
  bool isTight(0);

  if (m_electron) {  
    lept = evt.electron()[0].mom();	  
    d0sig = evt.electron()[0].sd0();

    if (evt.MC_w1l_pdgId()==11){
  	dr = lept.DeltaR(evt.MC_w1l());
  	if (dr<lep_drMax)	leptMa_pdgId = evt.MC_w1l_pdgId();
    
    }else if (evt.MC_w1l_pdgId()==15){
  	dr = lept.DeltaR(evt.MC_w1l());
  	if (dr<tau_drMax)	leptMa_pdgId = evt.MC_w1l_pdgId();

    }else if (evt.MC_w2l_pdgId()==-11){
  	dr = lept.DeltaR(evt.MC_w2l());
  	if (dr<lep_drMax)	leptMa_pdgId = evt.MC_w2l_pdgId();

    }else if (evt.MC_w2l_pdgId()==-15){
  	dr = lept.DeltaR(evt.MC_w2l());
  	if (dr<tau_drMax)	leptMa_pdgId = evt.MC_w2l_pdgId();


    }else if (abs(evt.MC_w1l_pdgId())==13 || abs(evt.MC_w2l_pdgId())==13)	std::cout << "reco electron and truth muon" << std::endl;
  	   
  } else {
    lept = evt.muon()[0].mom();   
    d0sig = evt.muon()[0].sd0();
 
    if (evt.MC_w1l_pdgId()==13 && lept.DeltaR(evt.MC_w1l())<lep_drMax){
  	leptMa_pdgId = evt.MC_w1l_pdgId();
    }else if (evt.MC_w2l_pdgId()==-13 && lept.DeltaR(evt.MC_w2l())<lep_drMax){
  	leptMa_pdgId = evt.MC_w2l_pdgId();
    }
    else if (evt.MC_w1l_pdgId()==15 && lept.DeltaR(evt.MC_w1l())<tau_drMax){
  	leptMa_pdgId = evt.MC_w1l_pdgId();
    }else if (evt.MC_w2l_pdgId()==-15 && lept.DeltaR(evt.MC_w2l())<tau_drMax){
  	leptMa_pdgId = evt.MC_w2l_pdgId();
    }
    else if (abs(evt.MC_w1l_pdgId())==11 || abs(evt.MC_w2l_pdgId())==11)	std::cout << "reco muon and truth electron" << std::endl;    


  }//m_electron

  h->h1D("eff_d0sig", "", suffix)  ->Fill(d0sig, weight);
  
  //nB-tagged jets 
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       	if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          	nTrkBtagged += 1;
  }
  
//  if (nTrkBtagged!=0)	return;	
//   if (nTrkBtagged < 1) return;		

  //std::cout<<"reached HERE "<<std::endl;			
  if (leptMa_pdgId!=0)	{

     GetRealHistograms(evt, weight, suffix, "");
     if(nTrkBtagged == 0)
     GetRealHistograms(evt, weight, suffix, "btag0_");
     if(nTrkBtagged >= 1)
     GetRealHistograms(evt, weight, suffix, "btag1_");
	
   }
}//AnaTtresQCD::runRealRateQCDCR_2016


void AnaTtresQCD::GetRealHistograms(const Event &evt, double weight, const std::string &suffix, const std::string &btag){

  HistogramService *h = &m_hSvc;
  
  //lepton
  TLorentzVector lept; 
  float Dz0(0); 
  float d0sig(0); 
  float ptvarcone(-1);
  float topoetcone(-1);

  if (m_electron) {
     lept  	= evt.electron()[0].mom();
     Dz0    	= evt.electron()[0].Dz0();
     d0sig 	= evt.electron()[0].sd0();
     ptvarcone	= evt.electron()[0].ptvarcone20();
     topoetcone	= evt.electron()[0].topoetcone20();

  } else{			
     lept  	= evt.muon()[0].mom();
     Dz0    	= evt.muon()[0].Dz0();
     d0sig 	= evt.muon()[0].sd0();
     ptvarcone	= evt.muon()[0].ptvarcone30();
     topoetcone	= evt.muon()[0].topoetcone20();
  }//m_electron 
  
  h->h2D(btag + "eff_MaLep_pt_eta", "", suffix)->Fill(lept.Perp()*1e-3, lept.Eta());  

  //deltaR between lepton and the closest narrow jet
  float closejl_deltaR  = 99;
  float deltaR_tmp	= 99;
  int closejl_idx	= -1;
  float deltaRapidity2  = 99;
  float deltaPhi2	= 99;

  size_t jet_idx = 0;
  for (; jet_idx < evt.jet().size(); ++jet_idx){  

    deltaRapidity2 = pow(evt.jet()[jet_idx].mom().Rapidity() - lept.Rapidity(), 2);	     
    deltaPhi2 = pow(evt.jet()[jet_idx].mom().DeltaPhi(lept), 2);	   
    deltaR_tmp = sqrt(deltaPhi2 + deltaRapidity2);
    
    if (deltaR_tmp < closejl_deltaR){
  	closejl_deltaR = deltaR_tmp;
  	closejl_idx = jet_idx;
    }	
  }//for    
   
  if (closejl_deltaR<99)  {
  	h->h2D(btag + "eff_LepPt_DR", 	"", suffix)->Fill(lept.Perp()*1e-3, closejl_deltaR, weight);
        h->h2D(btag + "eff_LepPt_DR_mbin",   "", suffix)->Fill(lept.Perp()*1e-3, closejl_deltaR, weight);
        h->h2D(btag + "eff_LepPt_DR_sbin",   "", suffix)->Fill(lept.Perp()*1e-3, closejl_deltaR, weight);
	h->h1D(btag + "eff_lepPt", 	"", suffix)->Fill(lept.Perp()*1e-3, weight);
        h->h1D(btag + "eff_minDeltaR", "", suffix)->Fill(closejl_deltaR, weight);
        h->h1D(btag + "eff_Dz0sin", 	"", suffix)->Fill(Dz0, weight);
        h->h1D(btag + "eff_d0sig",  	"", suffix)->Fill(d0sig, weight);
	
	float cosDPhi = cos(evt.met().DeltaPhi(lept));
	float mWt = sqrt(2. * lept.Perp() * evt.met().Perp() * (1. - cosDPhi))*1e-3; 
  	float MET = evt.met().Perp()*1e-3;
	
	h->h1D(btag + "eff_MET", 	"", suffix)->Fill(MET, weight);
  	h->h1D(btag + "eff_mwt", 	"", suffix)->Fill(mWt, weight);	
	h->h1D(btag + "eff_mwt_met", 	"", suffix)->Fill(mWt+MET, weight);	
        h->h1D(btag + "eff_cosDPhi",   "", suffix)->Fill(cosDPhi, weight);
	h->h1D(btag + "eff_nJets", 	"", suffix)->Fill(evt.jet().size(), weight);
	
	h->h2D(btag + "eff_mwt_met_map", 	"", suffix)->Fill(MET, mWt, weight);
	h->h1D(btag + "eff_ptvarcone", 	"", suffix)->Fill(ptvarcone*1e-3, weight);
        h->h1D(btag + "eff_topoetcone", 	"", suffix)->Fill(topoetcone*1e-3, weight);
     
        h->h1D(btag + "eff_topoetcone_effBins","", suffix)->Fill(topoetcone*1e-3, weight);
        h->h2D(btag + "eff_lepPt_topoetcone",  "", suffix)->Fill(lept.Pt()*1e-3 ,topoetcone*1e-3 , weight);
     
	//nB-tagged jets 
  	int nTrkBtagged = 0; 
  	for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       		if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          		nTrkBtagged += 1;
  
  	}//for     
  	h->h1D(btag + "eff_nTrkBtagJets", "", suffix)->Fill(nTrkBtagged, weight); 
	
	h->h2D(btag + "eff_lepPt_cosDPhi",    "", suffix)->Fill(lept.Pt()*1e-3 ,cosDPhi , weight); 
	 
	// For muon channel
	if(closejl_deltaR>0.6){
	  h->h2D(btag + "eff_lepPt_met_highDR",	"", suffix)->Fill(lept.Pt()*1e-3, MET, weight);
	  h->h2D(btag + "eff_lepPt_cosDPhi_highDR",	"", suffix)->Fill(lept.Pt()*1e-3, cosDPhi , weight);
	  h->h2D(btag + "eff_mwt_met_map_highDR", 	"", suffix)->Fill(MET, mWt, weight);
	  h->h2D(btag + "eff_lepPt_topoetcone_highDR",  "", suffix)->Fill(lept.Pt()*1e-3 ,topoetcone*1e-3 , weight); 	  
        }else if(closejl_deltaR>0.4){
          h->h2D(btag + "eff_lepPt_met_medDR",		"", suffix)->Fill(lept.Pt()*1e-3, MET, weight);
	  h->h2D(btag + "eff_lepPt_cosDPhi_medDR",	"", suffix)->Fill(lept.Pt()*1e-3, cosDPhi , weight);
	  h->h2D(btag + "eff_mwt_met_map_medDR", 	"", suffix)->Fill(MET, mWt, weight);
	  h->h2D(btag + "eff_lepPt_topoetcone_medDR",  "", suffix)->Fill(lept.Pt()*1e-3 ,topoetcone*1e-3 , weight);	  
        }else{
          h->h2D(btag + "eff_lepPt_met_lowDR",		"", suffix)->Fill(lept.Pt()*1e-3, MET, weight);
	  h->h2D(btag + "eff_lepPt_cosDPhi_lowDR",	"", suffix)->Fill(lept.Pt()*1e-3, cosDPhi , weight);
	  h->h2D(btag + "eff_mwt_met_map_lowDR", 	"", suffix)->Fill(MET, mWt, weight);
	  h->h2D(btag + "eff_lepPt_topoetcone_lowDR",  "", suffix)->Fill(lept.Pt()*1e-3 ,topoetcone*1e-3 , weight);	  
        }//if	
			
   }//if (closejl_deltaR<99) 


}//AnaTtresQCD::GetRealHistograms 


//**************************
//******* Fake Rates *******
//**************************


void AnaTtresQCD::runFakeRateQCDCR_2015(const Event &evt, double weight, const std::string &suffix){
    
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted)
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;
  
  if (!m_boosted)	if(evt.jet().size()<4)	return;
      
    
  bool isTight;
  TLorentzVector lept; 
  float sd0(99);

  if (m_electron) {	
	isTight = evt.electron()[0].isTightPP();
	lept  = evt.electron()[0].mom();
	sd0    = evt.electron()[0].sd0();

  } else{
	isTight = evt.muon()[0].isTight();
	lept  = evt.muon()[0].mom();
        sd0    = evt.muon()[0].sd0();
		
  }//m_electron

  // ---------- GEN LEPTON MATCHING ------------- //
  if(m_dsid == 410000) {

  TLorentzVector lgen;

  if(fabs(evt.MC_Wdecay1_from_t_pdgId()) == 11 || fabs(evt.MC_Wdecay1_from_t_pdgId()) == 13 || fabs(evt.MC_Wdecay1_from_t_pdgId()) == 15)
   lgen = evt.MC_Wdecay1_from_t();

  else if(fabs(evt.MC_Wdecay2_from_t_pdgId()) == 11 || fabs(evt.MC_Wdecay2_from_t_pdgId()) == 13 || fabs(evt.MC_Wdecay2_from_t_pdgId()) == 15)
   lgen = evt.MC_Wdecay2_from_t();

  else if(fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 11 || fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 13  || fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 15)
   lgen = evt.MC_Wdecay1_from_tbar();

  else if(fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 11 || fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 13 || fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 15)
   lgen = evt.MC_Wdecay2_from_tbar();

  if (lept.DeltaR(lgen) > 0.2) return;

  } // if(m_dsid == 410000)
  // --------------- GEN LEPTON MATCHING --------- //
  //



  HistogramService *h = &m_hSvc;
  
  float mWt = sqrt(2. * lept.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(lept)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
      
  if(m_electron){
      if( (MET>20) || (MET+mWt)>60)	return;
  }else{
      if( (MET>20) || (MET+mWt)>60)	return;
      if(fabs(sd0)>3)			return;
  }//if
 
  int nTrkBtagged = 0;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
          nTrkBtagged += 1;
  }//for 

    if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }
  //if(log10(chi2Value) < 6)
  //  //std::cout<<"FILLING HISTO, logChi2 = "<<log10(chi2Value)<<std::endl;
  GetFakeHistograms(evt, weight, "", suffix, "");
  if(nTrkBtagged == 0)
  GetFakeHistograms(evt, weight, "", suffix, "btag0_");
  if(nTrkBtagged >= 1)
  GetFakeHistograms(evt, weight, "", suffix, "btag1_");
 
  
} // runFakeRateQCDCR_2015


void AnaTtresQCD::runFakeRateQCDCR_2016(const Event &evt, double weight, const std::string &suffix){
  
 
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted)
    if (!(evt.passes("bejetsIncluR_2016") || evt.passes("bmujetsQCDCR_2016")))
      return;

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2016") || evt.passes("rmujetsQCDCR_2016")))
      return;
  
  if (!m_boosted)	if(evt.jet().size()<4)	return;
         
 
  bool isTight;
  TLorentzVector lept; 
  float sd0(99);

  //---- extra addition for log10(chi2) ---- //
  int  igj3, igj4; // index for the Whad
  int igb3, igb4; // index for the b's
  int  ign1;  // index for the neutrino (because chi2 can test both pz solution)
  double chi2ming1, chi2ming1H, chi2ming1L;

  std::vector<TLorentzVector *> vjets;
  std::vector<bool> vjets_btagged;

  for (size_t z = 0; z < evt.jet().size(); ++z) {
      vjets.push_back(new TLorentzVector(0,0,0,0));
      vjets[z]->SetPtEtaPhiE(evt.jet()[z].mom().Perp(), evt.jet()[z].mom().Eta(), evt.jet()[z].mom().Phi(), evt.jet()[z].mom().E());

      TLorentzVector tmpAk4Jet = evt.jet()[z].mom();
      bool is_btagged(false), is_bmatched(false);
 
      for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx)
       {
         TLorentzVector tmpTJet = evt.tjet()[bidx].mom();
         if(tmpAk4Jet.DeltaR(tmpTJet) <= 0.4){
         if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
          is_btagged = true;
          // break;
          }
        } // for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx) 

      vjets_btagged.push_back(is_btagged);
    } // for (size_t z = 0; z < evt.jet().size(); ++z



  if (m_electron) {	
	isTight = evt.electron()[0].isTightPP();
	lept  = evt.electron()[0].mom();
	sd0    = evt.electron()[0].sd0();
  } else{
	isTight = evt.muon()[0].isTight();
	lept  = evt.muon()[0].mom();
        sd0    = evt.muon()[0].sd0();
  }//m_electron


  // -------------- GEN LEPTON MATCHING --------- //
  if(m_dsid == 410000) {

  TLorentzVector lgen;

  if(fabs(evt.MC_Wdecay1_from_t_pdgId()) == 11 || fabs(evt.MC_Wdecay1_from_t_pdgId()) == 13 || fabs(evt.MC_Wdecay1_from_t_pdgId()) == 15)
   lgen = evt.MC_Wdecay1_from_t();

  else if(fabs(evt.MC_Wdecay2_from_t_pdgId()) == 11 || fabs(evt.MC_Wdecay2_from_t_pdgId()) == 13 || fabs(evt.MC_Wdecay2_from_t_pdgId()) == 15)
   lgen = evt.MC_Wdecay2_from_t();

  else if(fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 11 || fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 13  || fabs(evt.MC_Wdecay1_from_tbar_pdgId()) == 15)
   lgen = evt.MC_Wdecay1_from_tbar();

  else if(fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 11 || fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 13 || fabs(evt.MC_Wdecay2_from_tbar_pdgId()) == 15)
   lgen = evt.MC_Wdecay2_from_tbar();

  if (lept.DeltaR(lgen) > 0.2) return;

  } // if(m_dsid == 410000)
  // --------------- GEN LEPTON MATCHING --------- //

  HistogramService *h = &m_hSvc;
 
  TLorentzVector met = evt.met();
  bool status = m_chi2.findMinChiSquare(&lept, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L);

  float chi2Value = 1000000; // log10(1000000) = 6
  if(status) 
  chi2Value = chi2ming1;
 
  float mWt = sqrt(2. * lept.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(lept)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  //std::cout<<"REACHED HERE"<<std::endl;    
  if(m_electron){
      if( (MET>20) || (MET+mWt)>60)	return;
  }else{
      if( (MET>20) || (MET+mWt)>60)	return;
      if(fabs(sd0)>3)			return;
      // if(log10(chi2Value) < 1.5 || log10(chi2Value) > 6)        return;
  }//if
      
  int nTrkBtagged = 0;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
          nTrkBtagged += 1;
  }//for 
 
    if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }
  //if(log10(chi2Value) < 6)
  //std::cout<<"FILLING HISTO, logChi2 = "<<log10(chi2Value)<<std::endl;
  GetFakeHistograms(evt, weight, "", suffix, "");
  if(nTrkBtagged == 0)
  GetFakeHistograms(evt, weight, "", suffix, "btag0_");
  if(nTrkBtagged >= 1)
  GetFakeHistograms(evt, weight, "", suffix, "btag1_");

  
} // runFakeRateQCDCR_2016


void AnaTtresQCD::runFakeRateWQCDCR_2015(const Event &evt, double weight, const std::string &suffix){
    
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted)
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;

  if (!m_boosted)	if(evt.jet().size()<2)	return;
    
  bool isTight;
  TLorentzVector lept; 
  float sd0(99);
    
  if (m_electron) {	
	isTight 	= evt.electron()[0].isTightPP();
	lept  		= evt.electron()[0].mom();
	sd0    		= evt.electron()[0].sd0();

  } else{
	isTight 	= evt.muon()[0].isTight();
	lept  		= evt.muon()[0].mom();
        sd0    		= evt.muon()[0].sd0();
			
  }//m_electron
  HistogramService *h = &m_hSvc;
  
  float mWt = sqrt(2. * lept.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(lept)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
      
  if(m_electron){
      if( (MET>20) || (MET+mWt)>60)	return;
  }else{
      if( (MET<20) || (MET+mWt)<60)	return;
      if(fabs(sd0)<5)			return;
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;		
    
  GetFakeHistograms(evt, weight, "", suffix, "");   
 
} //AnaTtresQCD::runFakeRateWQCDCR_2015


void AnaTtresQCD::runFakeRateWQCDCR_2016(const Event &evt, double weight, const std::string &suffix){
    
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted)
    if (!(evt.passes("bejetsIncluR_2016") || evt.passes("bmujetsQCDCR_2016")))
      return;

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2016") || evt.passes("rmujetsQCDCR_2016")))
      return;

  if (!m_boosted)	if(evt.jet().size()<2)	return;
        
  bool isTight;
  TLorentzVector lept; 
  float sd0(99);

  if (m_electron) {	
	isTight = evt.electron()[0].isTightPP();
	lept  = evt.electron()[0].mom();
  } else{
	isTight = evt.muon()[0].isTight();
	lept  = evt.muon()[0].mom();
        sd0    = evt.muon()[0].sd0();	
		
  }//m_electron
  HistogramService *h = &m_hSvc;
    
  float mWt = sqrt(2. * lept.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(lept)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
      
  if(m_electron){
      if( (MET>20) || (MET+mWt)>60)	return;
  }else{
      //if( (MET<20) || (MET+mWt)<60)	return;
      //if(fabs(sd0)<5)			return;
      if( (MET<20) || (MET+mWt)>60)     return;
      if(fabs(sd0)<5)                   return;
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;	
  //if (nTrkBtagged < 1) return;	


  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }


  GetFakeHistograms(evt, weight, "", suffix, "");
  if(nTrkBtagged == 0)
  GetFakeHistograms(evt, weight, "", suffix, "btag0_");
  if(nTrkBtagged >= 1)
  GetFakeHistograms(evt, weight, "", suffix, "btag1_");

      
} //AnaTtresQCD::runFakeRateWQCDCR_2016



void AnaTtresQCD::GetFakeHistograms(const Event &evt, const double weight, const std::string &prefix, const std::string &suffix, const std::string &btag){
  // check channel

  HistogramService *h = &m_hSvc;
    
  h->h1D(btag +prefix+"eventN", "", suffix)->Fill(evt.eventNumber());
  h->h1D(btag +prefix+"runN", "", suffix)->Fill(evt.runNumber());

  TLorentzVector lept; 
  float Dz0(0); 
  float d0sig(0);
  float ptvarcone(-1);
  float topoetcone(-1);
  
   //lepton
  if (m_electron) {
     lept  = evt.electron()[0].mom();
     Dz0    = evt.electron()[0].Dz0();
     d0sig = evt.electron()[0].sd0();
     ptvarcone	= evt.electron()[0].ptvarcone20();
     topoetcone	= evt.electron()[0].topoetcone20();	  
  } else{			
     lept  	= evt.muon()[0].mom();
     Dz0    	= evt.muon()[0].Dz0();
     d0sig 	= evt.muon()[0].sd0();
     ptvarcone	= evt.muon()[0].ptvarcone30();
     topoetcone	= evt.muon()[0].topoetcone20();     
  }//m_electron 
    
  float theta = 2*atan(exp(-lept.Eta()));
  h->h1D(btag +prefix+"fake_Dz0sin", "",   suffix)->Fill(Dz0, weight);
  h->h1D(btag +prefix+"fake_d0sig",  "",   suffix)->Fill(d0sig  	, weight);
  
  // MET
  float MET = evt.met().Perp()*1e-3;

  h->h1D(btag +prefix+"fake_MET", 	    "", suffix)->Fill(MET, weight);
  h->h1D(btag +prefix+"fake_MET_effBins",    "", suffix)->Fill(MET, weight);

  h->h1D(btag +prefix+"fake_MET_phi",     "", suffix)->Fill(evt.met().Phi(), weight);
  
  // Transverse W mass
  float dPhi = evt.met().DeltaPhi(lept);
  float cosDPhi = std::cos(dPhi);
  float mWt = std::sqrt(2. * lept.Perp() * evt.met().Perp() * (1. - cosDPhi ))*1e-3; 

  //if(MET < 20)   
  //std::cout<<"MET = "<<MET<<", MWT = "<<mWt<<", MET+MWT = "<<(MET+mWt)<<std::endl;
 
  h->h1D(btag +prefix+"fake_mwt", 	    	"",     suffix)->Fill(mWt, weight);
  h->h1D(btag +prefix+"fake_mwt_effBins",	"",     suffix)->Fill(mWt, weight);

  h->h1D(btag +prefix+"fake_mwt_met", 	 	"",     suffix)->Fill(mWt+MET, weight);
  h->h1D(btag +prefix+"fake_mwt_met_effBins", 	"",     suffix)->Fill(mWt+MET, weight);
  
  h->h2D(btag +prefix+"fake_mwt_met_map", 	"", 	suffix)->Fill(MET, mWt, weight);
  
  h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi", 	    	"", suffix)->Fill(cosDPhi, weight);
  h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_effBins",  "", suffix)->Fill(cosDPhi, weight);
  
  h->h1D(btag +prefix+"fake_DeltaPhi",  			"", suffix)->Fill(dPhi, weight);
  
  // n Jets
  h->h1D(btag +prefix+"fake_nJets", "", suffix)->Fill(evt.jet().size(), weight);
  
  //nB-tagged jets 
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for     
  h->h1D(btag +prefix+"fake_nTrkBtagJets", "", suffix)->Fill(nTrkBtagged, weight); 
       
  //min DR(lept, jet) 
  float deltaRapidity2   = 99;
  float deltaPhi2	 = 99;
  float closejl_deltaR   = 99;
  float deltaR_tmp	 = 99;
  int   closejl_idx	 = -1;
  float closejl_pT	 = -1;
  
  size_t jet_idx = 0;
  for (; jet_idx < evt.jet().size(); ++jet_idx){     

    deltaRapidity2 = pow(evt.jet()[jet_idx].mom().Rapidity() - lept.Rapidity(), 2);	       
    deltaPhi2 = pow(evt.jet()[jet_idx].mom().DeltaPhi(lept), 2);      
    deltaR_tmp = sqrt(deltaPhi2 + deltaRapidity2);

     if (deltaR_tmp < closejl_deltaR){
  	closejl_deltaR = deltaR_tmp;
  	closejl_idx = jet_idx;  		   
     }   
  }//for
  
  if (closejl_deltaR<99){
  
     closejl_pT  = evt.jet().at(closejl_idx).mom().Pt();
     
     h->h1D(btag +prefix+"fake_lepPt_effBins", 	"", suffix)->Fill(lept.Pt()*1e-3 , weight);
     h->h1D(btag +prefix+"fake_lepPt",	  	"", suffix)->Fill(lept.Pt()*1e-3 , weight); 
     h->h1D(btag +prefix+"fake_lepEta",  	"", suffix)->Fill(lept.Eta()     , weight);
     h->h1D(btag +prefix+"fake_lepPhi",  	"", suffix)->Fill(lept.Phi()     , weight);

     h->h1D(btag +prefix+"fake_minDeltaR_effBins", 	"", suffix)->Fill(closejl_deltaR , weight); 
     h->h1D(btag +prefix+"fake_minDeltaR",	  	"", suffix)->Fill(closejl_deltaR , weight);
     
     h->h1D(btag +prefix+"fake_closJetPt_effBins", 	"", suffix)->Fill(closejl_pT*1e-3, weight);
     h->h1D(btag +prefix+"fake_closJetPt",	  	"", suffix)->Fill(closejl_pT*1e-3, weight);
     
     h->h1D(btag +prefix+"fake_closJetJVT", 		"", suffix)->Fill(evt.jet()[closejl_idx].jvt(), weight);
     
     h->h2D(btag +prefix+"fake_lepPt_lepEta",     "", suffix)->Fill(lept.Pt()*1e-3 ,fabs(lept.Eta()), weight);
     h->h2D(btag +prefix+"fake_lepPt_closJetPt",  "", suffix)->Fill(lept.Pt()*1e-3 ,closejl_pT*1e-3  , weight);
     h->h2D(btag +prefix+"fake_lepPt_minDeltaR",  "", suffix)->Fill(lept.Pt()*1e-3 ,closejl_deltaR, weight);
     h->h2D(btag +prefix+"fake_closJetPt_cosDPhi","", suffix)->Fill(closejl_pT*1e-3 ,cosDPhi  , weight);
     h->h2D(btag +prefix+"fake_minDeltaR_cosDPhi","", suffix)->Fill(closejl_deltaR  ,cosDPhi  , weight);
     h->h2D(btag +prefix+"fake_lepPt_met",        "", suffix)->Fill(lept.Pt()*1e-3 , MET, weight);
     h->h1D(btag +prefix+"fake_ptvarcone", 	    "", suffix)->Fill(ptvarcone*1e-3, weight);
     h->h1D(btag +prefix+"fake_topoetcone", 	    "", suffix)->Fill(topoetcone*1e-3, weight);
     
     h->h1D(btag +prefix+"fake_topoetcone_effBins", 	"", suffix)->Fill(topoetcone*1e-3, weight);
     h->h2D(btag +prefix+"fake_lepPt_topoetcone",    	"", suffix)->Fill(lept.Pt()*1e-3 ,topoetcone*1e-3 , weight);
     if (closejl_deltaR<0.4)
     h->h2D(btag +prefix+"fake_lepPt_topoetcone_lowDR",       "", suffix)->Fill(lept.Pt()*1e-3 ,topoetcone*1e-3 , weight);   
     else
     h->h2D(btag +prefix+"fake_lepPt_topoetcone_highDR",       "", suffix)->Fill(lept.Pt()*1e-3 ,topoetcone*1e-3 , weight); 
 
     h->h2D(btag +prefix+"fake_lepPt_cosDPhi",    "", suffix)->Fill(lept.Pt()*1e-3 ,cosDPhi , weight); 
     
     if(abs(lept.Eta())> 1.0){
  	h->h2D(btag +prefix+"fake_lepPt_closJetPt_highEta",   "", suffix)->Fill(lept.Pt()*1e-3, closejl_pT*1e-3 , weight);
	h->h2D(btag +prefix+"fake_lepPt_cosDPhi_highEta",     "", suffix)->Fill(lept.Pt()*1e-3, cosDPhi , weight);
        h->h1D(btag +prefix+"fake_minDeltaR_highEta",	     	"", suffix)->Fill(closejl_deltaR , weight);
	h->h1D(btag +prefix+"fake_lepPt_highEta",	  	"", suffix)->Fill(lept.Pt()*1e-3 , weight); 
	h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_highEta",	"", suffix)->Fill(cosDPhi, weight);

     }else{
  	h->h2D(btag +prefix+"fake_lepPt_closJetPt_lowEta",  	"", suffix)->Fill(lept.Pt()*1e-3, closejl_pT*1e-3 , weight);
	h->h2D(btag +prefix+"fake_lepPt_cosDPhi_lowEta",    	"", suffix)->Fill(lept.Pt()*1e-3, cosDPhi , weight);
	h->h1D(btag +prefix+"fake_minDeltaR_lowEta",	   	"", suffix)->Fill(closejl_deltaR , weight);
	h->h1D(btag +prefix+"fake_lepPt_lowEta",	  	"", suffix)->Fill(lept.Pt()*1e-3 , weight); 
	h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_lowEta",	"", suffix)->Fill(cosDPhi, weight);
	
     }

     if (mWt>20){
         h->h1D(btag +prefix+"fake_lepPt_highmWt",    	"", suffix)->Fill(lept.Pt()*1e-3, weight);
	 h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_highmWt",	"", suffix)->Fill(cosDPhi, weight);
         h->h2D(btag +prefix+"fake_lepPt_cosDPhi_highmWt",     "", suffix)->Fill(lept.Pt()*1e-3, cosDPhi, weight);
     
     } else if (mWt>10){
         h->h1D(btag +prefix+"fake_lepPt_medmWt",    		"", suffix)->Fill(lept.Pt()*1e-3, weight);
	 h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_medmWt",	"", suffix)->Fill(cosDPhi, weight);
         h->h2D(btag +prefix+"fake_lepPt_cosDPhi_medmWt",     "", suffix)->Fill(lept.Pt()*1e-3, cosDPhi, weight);
     
     } else{
         h->h1D(btag +prefix+"fake_lepPt_lowmWt",    		"", suffix)->Fill(lept.Pt()*1e-3, weight);
	 h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_lowmWt",	"", suffix)->Fill(cosDPhi, weight);
         h->h2D(btag +prefix+"fake_lepPt_cosDPhi_lowmWt",     "", suffix)->Fill(lept.Pt()*1e-3, cosDPhi, weight);
	 
     }//mWt

     // For electron channel
     if(closejl_deltaR>2.6){
  	h->h2D(btag +prefix+"fake_lepPt_closJetPt_highDR",   "", suffix)->Fill(lept.Pt()*1e-3, closejl_pT*1e-3 , weight);
	h->h2D(btag +prefix+"fake_lepPt_cosDPhi_highDR",     "", suffix)->Fill(lept.Pt()*1e-3, cosDPhi , weight);
     }else if(closejl_deltaR >1){
  	h->h2D(btag +prefix+"fake_lepPt_closJetPt_medDR",  "", suffix)->Fill(lept.Pt()*1e-3, closejl_pT*1e-3 , weight);
	h->h2D(btag +prefix+"fake_lepPt_cosDPhi_medDR",    "", suffix)->Fill(lept.Pt()*1e-3, cosDPhi , weight);
     }else{
  	h->h2D(btag +prefix+"fake_lepPt_closJetPt_lowDR",  "", suffix)->Fill(lept.Pt()*1e-3, closejl_pT*1e-3 , weight);
	h->h2D(btag +prefix+"fake_lepPt_cosDPhi_lowDR",    "", suffix)->Fill(lept.Pt()*1e-3, cosDPhi , weight);
     }
     
     // For muon channel
     if(closejl_deltaR>0.6){
	h->h2D(btag +prefix+"fake_lepPt_closJetPt_highDR",  	"", suffix)->Fill(lept.Pt()*1e-3 ,closejl_pT*1e-3  , weight);
	h->h2D(btag +prefix+"fake_lepPt_met_highDR",  		"", suffix)->Fill(lept.Pt()*1e-3 ,MET  , weight);
	h->h1D(btag +prefix+"fake_minDeltaR_highDR",	  	"", suffix)->Fill(closejl_deltaR , weight);
	h->h1D(btag +prefix+"fake_lepPt_highDR",	  		"", suffix)->Fill(lept.Pt()*1e-3 , weight); 
	h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_highDR",	"", suffix)->Fill(cosDPhi, weight);
	h->h2D(btag +prefix+"fake_mwt_met_map_highDR", 		"", suffix)->Fill(MET, mWt, weight);
	
     }else if(closejl_deltaR>0.4){
	h->h2D(btag +prefix+"fake_lepPt_closJetPt_medDR",  	"", suffix)->Fill(lept.Pt()*1e-3 ,closejl_pT*1e-3  , weight);
	h->h2D(btag +prefix+"fake_lepPt_met_medDR",  		"", suffix)->Fill(lept.Pt()*1e-3 ,MET  , weight);
	h->h1D(btag +prefix+"fake_minDeltaR_medDR",	  	"", suffix)->Fill(closejl_deltaR , weight);
	h->h1D(btag +prefix+"fake_lepPt_medDR",	  		"", suffix)->Fill(lept.Pt()*1e-3 , weight); 
	h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_medDR",	"", suffix)->Fill(cosDPhi, weight);
	h->h2D(btag +prefix+"fake_mwt_met_map_medDR", 		"", suffix)->Fill(MET, mWt, weight);
	
     } else {
	h->h2D(btag +prefix+"fake_lepPt_closJetPt_lowDR",  	"", suffix)->Fill(lept.Pt()*1e-3 ,closejl_pT*1e-3  , weight);
	h->h2D(btag +prefix+"fake_lepPt_met_lowDR",  		"", suffix)->Fill(lept.Pt()*1e-3 ,MET  , weight);
	h->h1D(btag +prefix+"fake_minDeltaR_lowDR",	  	"", suffix)->Fill(closejl_deltaR , weight);
	h->h1D(btag +prefix+"fake_lepPt_lowDR",	  		"", suffix)->Fill(lept.Pt()*1e-3 , weight); 
	h->h1D(btag +prefix+"fake_cos_metPhi_lepPhi_lowDR",	"", suffix)->Fill(cosDPhi, weight);
	h->h2D(btag +prefix+"fake_mwt_met_map_lowDR", 		"", suffix)->Fill(MET, mWt, weight);
	
     }
     
	
  }//if(closejl_deltaR<99)	
  
  //min DR(lept, tjet) 
  deltaRapidity2   = 99;
  deltaPhi2	   = 99;
  deltaR_tmp	   = 99;
  
  float closetjl_deltaR   = 99;  
  int   closetjl_idx	 = -1;
  float closetjl_pT	 = -1;
  
  size_t tjet_idx = 0;
  for (; tjet_idx < evt.tjet().size(); ++tjet_idx){     

    deltaRapidity2 = pow(evt.tjet()[tjet_idx].mom().Rapidity() - lept.Rapidity(), 2);	       
    deltaPhi2 = pow(evt.tjet()[tjet_idx].mom().DeltaPhi(lept), 2);      
    deltaR_tmp = sqrt(deltaPhi2 + deltaRapidity2);

     if (deltaR_tmp < closetjl_deltaR){
  	closetjl_deltaR = deltaR_tmp;
  	closetjl_idx = tjet_idx;  		   
     }   
  }//for
  
  if (closetjl_deltaR<99){
  
     h->h1D(btag +prefix+"fake_minDeltaR_tjet_effBins", 	"", suffix)->Fill(closetjl_deltaR , weight); 
     h->h1D(btag +prefix+"fake_minDeltaR_tjet",	  	"", suffix)->Fill(closetjl_deltaR , weight);
     
     h->h1D(btag +prefix+"fake_closJetPt_tjet_effBins", 	"", suffix)->Fill(closetjl_pT*1e-3, weight);
     h->h1D(btag +prefix+"fake_closJetPt_tjet",	  	"", suffix)->Fill(closetjl_pT*1e-3, weight);
  
  }//if(closetjl_deltaR<99)
  
}//AnaTtresQCD::GetFakeHistograms

void AnaTtresQCD::get1Drates(float &rate, float &rate_err, TH1F* rate_map, float x){
   
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


void AnaTtresQCD::get2Drates(float &rate, float &rate_err, TH2F* rate_map, float x, float y){
   
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



void AnaTtresQCD::IniHistograms(std::string &suffix,  std::string &btag){

  //****Efficiency studies
  
  suffix = btag + suffix;
   
  Double_t eff_pT_bins_re[]     = {30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 120, 150, 300, 700};
  Double_t eff_pT_bins_rmu[]    = {25, 30, 35, 40, 50, 60, 70, 85, 100, 150, 300, 700};
  Double_t eff_pT_bins_be[]     = {30, 40, 60, 120, 700};
  Double_t eff_pT_bins_bmu[]    = {25, 35, 50, 100, 700};

  Double_t DR_bins_re[]  	= {0., 0.4, 0.6, 1.0, 1.5, 7.0};
  
  //Double_t DR_bins_rmu[] 	= {0., 0.4, 7.0};
  //Double_t DR_bins_rmu[]        = {0., 0.2, 0.35, 0.4, 0.5, 0.55, 0.6, 1.0, 1.5, 2.5, 7.0};
  Double_t DR_bins_rmu[]        = {0., 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
  Double_t DR_bins_rmu_mbin[]   = {0., 0.2, 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
  Double_t DR_bins_rmu_sbin[]   = {0., 0.1, 0.2, 0.3, 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
  
  Double_t DR_bins_be[]  	= {0., 0.4, 0.6, 0.8, 1.0, 1.5};
  Double_t DR_bins_bmu[] 	= {0., 0.2, 0.3, 0.4, 0.5, 1.5};
    
  //MC variables:  
  m_hSvc.create2D(suffix+"eff_MCe_pt_eta",   "; pt of electron(truth) [GeV]; #eta of electron(truth)", 40, 0, 500, 24, -3.5, 3.5);
  m_hSvc.create2D(suffix+"eff_MCmu_pt_eta",  "; pt of muon(truth) [GeV]; #eta of muon(truth)", 40, 0, 500, 24, -3.5, 3.5);
  m_hSvc.create2D(suffix+"eff_MCtau_pt_eta", "; pt of tau(truth) [GeV]; #eta of tau(truth)", 40, 0, 500, 24, -3.5, 3.5);
   
  m_hSvc.create1D(suffix+"eff_d0sig", 		"; d0sig; Events", 30, -15, 15);  
  m_hSvc.create1D(suffix+"eff_Dz0sin",          "; #Delta z0*sin(#theta); Events", 20, 0, 1);
  m_hSvc.create1D(suffix+"eff_MET", 	        "; Missing E_{T} [GeV]; Events", 100, 0, 500);
  m_hSvc.create1D(suffix+"eff_mwt",	        "; W transverse mass [GeV]; Events", 50, 0, 500);
  m_hSvc.create1D(suffix+"eff_mwt_met",	        "; MET + MWT [GeV]; Events", 50, 0, 500);
  m_hSvc.create1D(suffix+"eff_cosDPhi",  	"; Cos( #Delta #phi(met, lept) ); Events", 50, -1, 1);  
  m_hSvc.create1D(suffix+"eff_nJets",           "; number of jets ; Events", 10, -0.5, 9.5);
  m_hSvc.create1D(suffix+"eff_nTrkBtagJets",    "; number of b-tagged track jets ; Events", 10, 0, 10);
    
  //Matched reco leptons  
  m_hSvc.create2D(suffix+"eff_MaLep_pt_eta", "; pt of Matched lepton [GeV]; #eta of Matched lepton", 40, 0, 500, 24, -3., 3.);
  m_hSvc.create1D(suffix+"eff_ptvarcone",  "; ptvarcone; Events", 50, -10, 40); 
  m_hSvc.create1D(suffix+"eff_topoetcone", "; topoetcone; Events", 40, -10, 30);


  //****Fake studies  
   
  m_hSvc.create1D(suffix+"fake_Ma_pdgID_wl",       	"; pdg ID; Events", 40, -20, 20);   
  m_hSvc.create1D(suffix+"fake_ratioPt_wl",		"; lept Pt/true lept Pt [GeV]; Events", 60, 0, 3);
  m_hSvc.create1D(suffix+"fake_MC_ttbar_type_wl", 	"; MC ttbar type; Events", 2002, -2, 2000); 
  m_hSvc.create1D(suffix+"fake_mtt_wl", 		"; (MC) mass of the top-antitop system [GeV]; Events", 50, 0, 2000);
  
  m_hSvc.create1D(suffix+"fake_Ma_pdgID_wh",       	"; pdg ID; Events", 40, -20, 20);   
  m_hSvc.create1D(suffix+"fake_ratioPt_wh",		"; lept Pt/true lept Pt [GeV]; Events", 60, 0, 3);
  m_hSvc.create1D(suffix+"fake_MC_ttbar_type_wh", 	"; MC ttbar type; Events", 2002, -2, 2000); 
  m_hSvc.create1D(suffix+"fake_mtt_wh", 		"; (MC) mass of the top-antitop system [GeV]; Events", 50, 0, 2000);
  
  m_hSvc.create1D(suffix+"fake_Ma_pdgID_strange",       "; pdg ID; Events", 40, -20, 20);   
  m_hSvc.create1D(suffix+"fake_ratioPt_strange",	"; lept Pt/true lept Pt [GeV]; Events", 60, 0, 3);
  m_hSvc.create1D(suffix+"fake_MC_ttbar_type_strange", 	"; MC ttbar type; Events", 2002, -2, 2000); 
  m_hSvc.create1D(suffix+"fake_mtt_strange", 		"; (MC) mass of the top-antitop system [GeV]; Events", 50, 0, 2000);
  
  m_hSvc.create1D(suffix+"fake_mtt_bef", 		"; (MC) mass of the top-antitop system [GeV]; Events", 50, 0, 2000);
  m_hSvc.create1D(suffix+"fake_MC_ttbar_type_bef", 	"; MC ttbar type; Events", 2002, -2, 2000); 
  
  m_hSvc.create1D(suffix+"fake_mtt_aft", 		"; (MC) mass of the top-antitop system [GeV]; Events", 50, 0, 2000);
  m_hSvc.create1D(suffix+"fake_MC_ttbar_type_aft", 	"; MC ttbar type; Events", 2002, -2, 2000); 

  m_hSvc.create1D(suffix+"eventN",              "; N event; Events", 10000, 0, 1e10);
  m_hSvc.create1D(suffix+"runN",   	        "; N run; Events", 1000, 0, 1e10); 
  m_hSvc.create1D(suffix+"fake_Dz0sin",         "; #Delta z0*sin(#theta); Events", 20, 0, 1);
  m_hSvc.create1D(suffix+"fake_MET", 	        "; Missing E_{T} [GeV]; Events", 100, 0, 500);
  m_hSvc.create1D(suffix+"fake_MET_phi",        "; Missing E_{T} #phi; Events", 20, -4.0, 4.0);	
  m_hSvc.create1D(suffix+"fake_mwt", 	        "; W transverse mass [GeV]; Events", 50, 0, 500);
  m_hSvc.create1D(suffix+"fake_mwt_met",        "; MET + MWT [GeV]; Events", 50, 0, 500);
  m_hSvc.create1D(suffix+"fake_closJetPt",      "; Pt of closest jet to lepton [GeV]; Events", 100, 25, 525);
  m_hSvc.create1D(suffix+"fake_nJets",          "; number of jets ; Events", 10, -0.5, 9.5);
  m_hSvc.create1D(suffix+"fake_nTrkBtagJets",   "; number of b-tagged track jets ; Events", 10, 0, 10); 
   
  double varBin1[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 340, 380, 450, 700};
  int varBinN1 = sizeof(varBin1)/sizeof(double) - 1;
   
  double CdPhiBin[] = {-1.0, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 0.95, 1.0};  
  int CdPhiBinN = sizeof(CdPhiBin)/sizeof(double) - 1; 
  
  double CdPhi_lDR_Bin[] = {-1.0,  -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.0};
  int CdPhi_lDR_BinN = sizeof(CdPhi_lDR_Bin)/sizeof(double) - 1;
  
  double CdPhi_mDR_Bin[] = {-1.0,  -0.6, -0.3, 0.0, 0.3, 0.6, 0.8, 0.9, 0.95, 1.0};
  int CdPhi_mDR_BinN = sizeof(CdPhi_mDR_Bin)/sizeof(double) - 1;
  
  double CdPhi_hDR_Bin[] = {-1.0,  -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.0};
  int CdPhi_hDR_BinN = sizeof(CdPhi_hDR_Bin)/sizeof(double) - 1;
  
  double d0SigBin[] = {-100, -40, -30, -20, -15, -10, -6, -5, -3, -2, -1, -0.5, -0.2, 0.0, 0.2, 0.5, 1, 2, 3 , 5, 6, 10, 15, 20, 30, 40, 100};
  int d0SigBinN = sizeof(d0SigBin)/sizeof(double) - 1; 
  
  Double_t fake_pT_bins_re_lowDR[]  	 	= {30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 700};
  Double_t fake_closeLJpT_bins_re_lowDR[]  	= {20, 40, 60, 80, 100, 700};

  Double_t fake_pT_bins_re_medDR[]  	 	= {30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 700};
  Double_t fake_closeLJpT_bins_re_medDR[]  	= {20, 40, 60, 80, 100, 700};

  Double_t fake_pT_bins_re_highDR[]  	 	= {30, 35, 40, 50, 60, 70, 100, 120, 700};
  Double_t fake_closeLJpT_bins_re_highDR[]  	= {20, 40, 60, 80, 100, 700};
  
  m_hSvc.create1DVar(suffix+"fake_d0sig",          		"; d0sig; Events", d0SigBinN, d0SigBin);
  m_hSvc.create1DVar(suffix+"fake_lepPt",  	        	"; lepton p_{T} [GeV]; Events", varBinN1, varBin1);
  m_hSvc.create1DVar(suffix+"fake_lepPt_lowEta",  	        "; lepton p_{T} [GeV]; Events", varBinN1, varBin1);
  m_hSvc.create1DVar(suffix+"fake_lepPt_highEta",  	        "; lepton p_{T} [GeV]; Events", varBinN1, varBin1);
    
  m_hSvc.create1DVar(suffix+"fake_cos_metPhi_lepPhi_highmWt",  	        "; lepton p_{T} [GeV]; Events", CdPhi_hDR_BinN,  CdPhi_hDR_Bin);
  m_hSvc.create1DVar(suffix+"fake_cos_metPhi_lepPhi_medmWt",  	        "; lepton p_{T} [GeV]; Events", CdPhi_mDR_BinN,  CdPhi_mDR_Bin);
  m_hSvc.create1DVar(suffix+"fake_cos_metPhi_lepPhi_lowmWt",  	        "; lepton p_{T} [GeV]; Events", CdPhi_lDR_BinN,  CdPhi_lDR_Bin);
  
  
  m_hSvc.create1DVar(suffix+"fake_cos_metPhi_lepPhi_effBins",  	"; Cos( #Delta #phi(met, lept) ); Events", CdPhiBinN, CdPhiBin);
  
  m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_lowEta",  	"; Cos( #Delta #phi(met, lept) ); Events", 50, -1, 1);
  m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_highEta",  	"; Cos( #Delta #phi(met, lept) ); Events", 50, -1, 1);
  m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi",  		"; Cos( #Delta #phi(met, lept) ); Events", 50, -1, 1);
  m_hSvc.create1D(suffix+"fake_DeltaPhi",  			"; #Delta #Phi(met, lept); Events", 50, -10, 10);
  m_hSvc.create1D(suffix+"fake_closJetJVT", 			"; JVT of closest Jet to lepton ; Events", 20, -1, 1); 
  m_hSvc.create1D(suffix+"fake_minDeltaR_highEta",      	"; min #Delta R(lep, jet); Events", 70, 0, 7); 
  m_hSvc.create1D(suffix+"fake_minDeltaR_lowEta",      		"; min #Delta R(lep, jet); Events", 70, 0, 7); 
  m_hSvc.create1D(suffix+"fake_minDeltaR",      		"; min #Delta R(lep, jet); Events", 70, 0, 7); 
  m_hSvc.create1D(suffix+"fake_minDeltaR_tjet",      		"; min #Delta R(lep, tjet); Events", 70, 0, 7); 
  m_hSvc.create1D(suffix+"fake_lepPhi",         		"; lepton #phi; Events", 20, -4.0, 4.0);  
  m_hSvc.create1D(suffix+"fake_lepEta",         		"; lepton #eta; Events", 20, -4.0, 4.0); 
  m_hSvc.create1D(suffix+"fake_ptvarcone",         		"; ptvarcone; Events", 50, -10, 40); 
  m_hSvc.create1D(suffix+"fake_topoetcone",         		"; topoetcone; Events", 40, -10, 30); 
   
  
  double metBin[] = {0, 20, 30, 40, 50, 70, 100, 500};
  int metBinN = sizeof(metBin)/sizeof(double) - 1; 
  m_hSvc.create1DVar(suffix+"fake_MET_effBins", 	      "; Missing E_{T} [GeV]; Events", metBinN, metBin);

  double mwtBin[] = {0, 10, 20, 45, 60, 80, 100, 200, 600};
  int mwtBinN = sizeof(mwtBin)/sizeof(double) - 1;
  m_hSvc.create1DVar(suffix+"fake_mwt_effBins",               "; W transverse mass [GeV]; Events", mwtBinN, mwtBin);

  double mwtmetBin[] = {0, 60, 65, 70, 75, 80, 85, 90, 100, 120, 150, 400};
  int mwtmetBinN = sizeof(mwtmetBin)/sizeof(double) - 1; 
  m_hSvc.create1DVar(suffix+"fake_mwt_met_effBins", 	      "; MET + MWT [GeV]; Events", mwtmetBinN, mwtmetBin);
  
  m_hSvc.create2DVar(suffix+"eff_mwt_met_map",	       "; MET [GeV]; MWT [GeV]", metBinN, metBin, mwtBinN, mwtBin);
  m_hSvc.create2DVar(suffix+"fake_mwt_met_map",        "; MET [GeV]; MWT [GeV]", metBinN, metBin, mwtBinN, mwtBin);
    
  double metBin_lowDR[] = {0, 20, 30, 40, 60, 80, 100, 150, 500};
  int metBinN_lowDR = sizeof(metBin_lowDR)/sizeof(double) - 1; 
  double mwtBin_lowDR[] = {0, 15, 60, 100, 200, 600};
  int mwtBinN_lowDR = sizeof(mwtBin_lowDR)/sizeof(double) - 1;
    
  m_hSvc.create2DVar(suffix+"eff_mwt_met_map_lowDR",   "; MET [GeV]; MWT [GeV]", metBinN_lowDR, metBin_lowDR, mwtBinN_lowDR, mwtBin_lowDR);
  m_hSvc.create2DVar(suffix+"fake_mwt_met_map_lowDR",  "; MET [GeV]; MWT [GeV]", metBinN_lowDR, metBin_lowDR, mwtBinN_lowDR, mwtBin_lowDR);
  
  double metBin_medDR[] = {0, 20, 30, 40, 70, 100, 500};
  int metBinN_medDR = sizeof(metBin_medDR)/sizeof(double) - 1; 
  double mwtBin_medDR[] = {0, 15, 60, 100, 200, 600};
  int mwtBinN_medDR = sizeof(mwtBin_medDR)/sizeof(double) - 1;
  
  m_hSvc.create2DVar(suffix+"eff_mwt_met_map_medDR",   "; MET [GeV]; MWT [GeV]", metBinN_medDR, metBin_medDR, mwtBinN_medDR, mwtBin_medDR);
  m_hSvc.create2DVar(suffix+"fake_mwt_met_map_medDR",  "; MET [GeV]; MWT [GeV]", metBinN_medDR, metBin_medDR, mwtBinN_medDR, mwtBin_medDR);
  
  double metBin_highDR[] = {0, 20, 30, 40, 50, 60, 100, 150, 500};
  int metBinN_highDR = sizeof(metBin_highDR)/sizeof(double) - 1; 
  double mwtBin_highDR[] = {0, 15, 30, 45, 60, 75, 100, 200, 600};
  int mwtBinN_highDR = sizeof(mwtBin_highDR)/sizeof(double) - 1;
  
  m_hSvc.create2DVar(suffix+"eff_mwt_met_map_highDR",  "; MET [GeV]; MWT [GeV]", metBinN_highDR, metBin_highDR, mwtBinN_highDR, mwtBin_highDR);
  m_hSvc.create2DVar(suffix+"fake_mwt_met_map_highDR", "; MET [GeV]; MWT [GeV]", metBinN_highDR, metBin_highDR, mwtBinN_highDR, mwtBin_highDR);
    
  Double_t fake_pT_bins_re[]  	 	= {30, 35, 40, 50, 60, 80, 100, 120, 150, 700};
  Double_t fake_closeLJpT_bins_re[]  	= {20, 40, 60, 80, 100, 700};
  Double_t fake_minDeltaR_bins_re[]    	= {0.0, 0.4, 1.0, 7};
  Double_t fake_leptEta_bins_re[]    	= {0.0, 1.8, 2.0, 2.5};
 
  Double_t fake_pT_bins_re_lowEta[]  	 	= {30, 35, 40, 60, 80, 100, 120, 700};
  int pT_binsN_re_lowEta = sizeof(fake_pT_bins_re_lowEta)/sizeof(double) - 1;
  
  Double_t fake_pT_bins_re_highEta[]  	 	= {30, 35, 40, 60, 100, 120, 700};
  int pT_binsN_re_highEta = sizeof(fake_pT_bins_re_highEta)/sizeof(double) - 1;

  
  Double_t fake_pT_bins_rmu[] 	 	= {25, 30, 35, 40, 50, 100, 700};  
  Double_t fake_pT_bins_lDR_rmu[] 	= {25, 30, 40, 50, 100, 700}; 
  Double_t fake_pT_bins_mDR_rmu[] 	= {25, 30, 35, 40, 50, 100, 700}; 
  Double_t fake_pT_bins_hDR_rmu[] 	= {25, 30, 35, 40, 50, 100, 700}; 
  Double_t fake_leptEta_bins_rmu[]    	= {0.0, 1.6, 2.5};
  Double_t fake_closeLJpT_bins_rmu[] 	= {25, 35, 50, 100, 700};
  
  Double_t fake_pT_bins_be[]  	 	= {30, 60, 120, 700};
  Double_t fake_leptEta_bins_be[]    	= {0.0, 1.6, 2.5};
  Double_t fake_closeLJpT_bins_be[]  	= {20, 40, 60, 80, 100, 700};
  
  Double_t fake_pT_bins_bmu[] 	 	= {25, 50, 100, 700};
  Double_t fake_leptEta_bins_bmu[]    	= {0.0, 1.6 ,2.5};
  Double_t fake_closeLJpT_bins_bmu[] 	= {20, 40, 60, 80, 100, 700}; 

  //double topoetconeBin[] = {-8, -5, -2, 1, 3, 6, 10, 30}; 
  //double topoetconeBin[] = {-8, 1, 6, 10, 30}; // --- for 2015 
  double topoetconeBin[] = {-8, 1, 3, 6, 10, 30}; // --- for 2016
  int topoetconeBinN = sizeof(topoetconeBin)/sizeof(double) - 1;
    
  int eff_N_pT_bins(0);
  int N_DR_bins(0);
  int N_DR_bins_mbin(0);
  int N_DR_bins_sbin(0);
  int fake_N_pT_bins(0);
  int fake_N_pT_bins_lDR(0);
  int fake_N_pT_bins_mDR(0);
  int fake_N_pT_bins_hDR(0);
  int fake_N_leptEta_bins(0);
  int fake_N_closeLJpT_bins(0);
  int fake_N_pT_bins_re_lowDR(0); 
  int fake_N_closeLJpT_bins_re_lowDR(0);
  int fake_N_pT_bins_re_medDR(0); 
  int fake_N_closeLJpT_bins_re_medDR(0);
  int fake_N_pT_bins_re_highDR(0); 
  int fake_N_closeLJpT_bins_re_highDR(0);
 
  if (m_boosted){
    
     if (m_electron)  	{
     
        eff_N_pT_bins 		= sizeof(eff_pT_bins_be)/sizeof(double) -1;
	N_DR_bins 		= sizeof(DR_bins_be)/sizeof(double) -1;	
	fake_N_pT_bins 		= sizeof(fake_pT_bins_be)/sizeof(double) -1;
	fake_N_leptEta_bins 	= sizeof(fake_leptEta_bins_be)/sizeof(double) -1;
        fake_N_closeLJpT_bins 	= sizeof(fake_closeLJpT_bins_be)/sizeof(double) -1;

     	//Eff
  	m_hSvc.create2DVar(suffix+"eff_LepPt_DR", 	"; Pt of Matched lepton [GeV]; min #Delta R(lep, jet)", eff_N_pT_bins, eff_pT_bins_be, N_DR_bins, DR_bins_be);
	m_hSvc.create1DVar(suffix+"eff_lepPt",  	"; lepton p_{T} [GeV]; Events", eff_N_pT_bins, eff_pT_bins_be);
        m_hSvc.create1DVar(suffix+"eff_minDeltaR",  	"; min #Delta R(lep, jet); Events", N_DR_bins, DR_bins_be);
	
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi",	        "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi_lowDR",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi_highDR",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_be, CdPhiBinN, CdPhiBin);
        m_hSvc.create2DVar(suffix+"eff_lepPt_topoetcone",       "; Pt of lepton [GeV]; topoetcone20 [GeV]", fake_N_pT_bins, fake_pT_bins_re, topoetconeBinN, topoetconeBin); 
  	//Fake 1D     
  	m_hSvc.create1DVar(suffix+"fake_lepPt_effBins", 	  "; Pt of lepton [GeV]; Eff", fake_N_pT_bins, fake_pT_bins_be);
	
	//Fake 1D
	m_hSvc.create1DVar(suffix+"fake_minDeltaR_effBins",   	  "; min #Delta R(lep, jet); Eff", N_DR_bins, DR_bins_be);
	m_hSvc.create1DVar(suffix+"fake_minDeltaR_tjet_effBins",  "; min #Delta R(lep, tjet); Eff", N_DR_bins, DR_bins_be);
	//Fake 1D
	m_hSvc.create1DVar(suffix+"fake_closJetPt_effBins",      "; Pt of closest jet to lepton [GeV]; Eff", fake_N_closeLJpT_bins, fake_closeLJpT_bins_be);
	
	//Fake 2D
	m_hSvc.create2DVar(suffix+"fake_lepPt_lepEta",   	  "; Pt of lepton [GeV]; |#eta| of lepton", fake_N_pT_bins, fake_pT_bins_be, fake_N_leptEta_bins, fake_leptEta_bins_be);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt",   	  "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]", fake_N_pT_bins, fake_pT_bins_be, fake_N_closeLJpT_bins, fake_closeLJpT_bins_be);
	m_hSvc.create2DVar(suffix+"fake_closJetPt_cosDPhi",       "; Pt of closest jet to lepton [GeV]; Cos( #Delta #phi(met, lept) )",   fake_N_closeLJpT_bins, fake_closeLJpT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi", 	   "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi_highLepPt", "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi_lowLepPt",  "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_be, CdPhiBinN, CdPhiBin);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi",   	  "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowEta",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highEta",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowDR",     "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highDR",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_medDR",     "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_be, CdPhiBinN, CdPhiBin);
	m_hSvc.create1DVar(suffix+"fake_lepPt_highDR",            "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_be);
	m_hSvc.create1D(suffix+"fake_minDeltaR_highDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_highDR",	  "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	m_hSvc.create1DVar(suffix+"fake_lepPt_medDR",             "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_be);
	m_hSvc.create1D(suffix+"fake_minDeltaR_medDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_medDR",    "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	m_hSvc.create1DVar(suffix+"fake_lepPt_lowDR",             "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_be);
	m_hSvc.create1D(suffix+"fake_minDeltaR_lowDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_lowDR",    "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_lowDR",   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]", fake_N_pT_bins, fake_pT_bins_be, fake_N_closeLJpT_bins, fake_closeLJpT_bins_be);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_medDR",   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]", fake_N_pT_bins, fake_pT_bins_be, fake_N_closeLJpT_bins, fake_closeLJpT_bins_be);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_highDR",  "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]", fake_N_pT_bins, fake_pT_bins_be, fake_N_closeLJpT_bins, fake_closeLJpT_bins_be);

	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR",   	  "; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_be, N_DR_bins, DR_bins_be);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_lowDPhi", "; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_be, N_DR_bins, DR_bins_be);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_highDPhi","; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_be, N_DR_bins, DR_bins_be);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_met",	   	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_be, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_lowDR", 	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_be, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_medDR",	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_be, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_highDR",	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_be, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_met_highLepPt", "; min #Delta R(lep, jet); MET [GeV]", N_DR_bins, DR_bins_be, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_met_lowLepPt",  "; min #Delta R(lep, jet); MET [GeV]", N_DR_bins, DR_bins_be, metBinN, metBin);
	
	//ljet
	m_hSvc.create1D(suffix+"fake_ljet_pt",    	"; Pt of largeR jet [GeV]; Events", 100, 0, 1000);
	m_hSvc.create1D(suffix+"fake_ljet_m",    	"; mass of largeR jet [GeV]; Events", 50, 0, 200);
	m_hSvc.create1D(suffix+"fake_ljet_tau32", 	"; #tau_{32}; Events", 20, 0, 1);
	m_hSvc.create1D(suffix+"fake_ljet_tau32_wta", 	"; #tau_{32} wta; Events", 20, 0, 1);

     }else{
     
        eff_N_pT_bins 		= sizeof(eff_pT_bins_bmu)/sizeof(double) -1;
	N_DR_bins 		= sizeof(DR_bins_bmu)/sizeof(double) -1;	
	fake_N_pT_bins 		= sizeof(fake_pT_bins_bmu)/sizeof(double) -1;
        fake_N_leptEta_bins 	= sizeof(fake_leptEta_bins_bmu)/sizeof(double) -1;
	fake_N_closeLJpT_bins 	= sizeof(fake_closeLJpT_bins_bmu)/sizeof(double) -1;
	
     	//Eff
  	m_hSvc.create2DVar(suffix+"eff_LepPt_DR", 	"; Pt of Matched lepton [GeV]; min #Delta R(lep, jet)", eff_N_pT_bins, eff_pT_bins_bmu, N_DR_bins, DR_bins_bmu);
	m_hSvc.create1DVar(suffix+"eff_lepPt",  	"; lepton p_{T} [GeV]; Events", eff_N_pT_bins, eff_pT_bins_bmu);
        m_hSvc.create1DVar(suffix+"eff_minDeltaR",  	"; min #Delta R(lep, jet); Events", N_DR_bins, DR_bins_bmu);
	
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi",	        "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi_lowDR",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi_highDR",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_bmu, CdPhiBinN, CdPhiBin);
  
  	//Fake 1D    
  	m_hSvc.create1DVar(suffix+"fake_lepPt_effBins", 	  	"; Pt of lepton [GeV]; Eff", fake_N_pT_bins, fake_pT_bins_bmu);
	m_hSvc.create1DVar(suffix+"fake_minDeltaR_tjet_effBins",   	"; min #Delta R(lep, tjet); Eff", N_DR_bins, DR_bins_bmu);
	m_hSvc.create1DVar(suffix+"fake_minDeltaR_effBins",   		"; min #Delta R(lep, jet); Eff", N_DR_bins, DR_bins_bmu);
	m_hSvc.create1DVar(suffix+"fake_closJetPt_effBins",   		"; Pt of closest jet to lepton [GeV]; Eff", fake_N_closeLJpT_bins, fake_closeLJpT_bins_bmu);
	
	//Fake 2D
	m_hSvc.create2DVar(suffix+"fake_lepPt_lepEta",      "; Pt of lepton [GeV]; |#eta| of lepton",		   fake_N_pT_bins, fake_pT_bins_bmu, fake_N_leptEta_bins, fake_leptEta_bins_bmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt",   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins, fake_pT_bins_bmu, fake_N_closeLJpT_bins, fake_closeLJpT_bins_bmu);		
	m_hSvc.create2DVar(suffix+"fake_closJetPt_cosDPhi", "; Pt of closest jet to lepton [GeV]; Cos( #Delta #phi(met, lept) )",   fake_N_closeLJpT_bins, fake_closeLJpT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi", "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi_highLepPt", "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi_lowLepPt",  "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_bmu, CdPhiBinN, CdPhiBin);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi",   	  "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowEta",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highEta",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowDR",     "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highDR",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_medDR",     "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_bmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create1DVar(suffix+"fake_lepPt_highDR",            "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_bmu);
	m_hSvc.create1D(suffix+"fake_minDeltaR_highDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_highDR",	  "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	m_hSvc.create1DVar(suffix+"fake_lepPt_medDR",             "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_bmu);
	m_hSvc.create1D(suffix+"fake_minDeltaR_medDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_medDR",    "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	m_hSvc.create1DVar(suffix+"fake_lepPt_lowDR",             "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_bmu);
	m_hSvc.create1D(suffix+"fake_minDeltaR_lowDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_lowDR",    "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_lowDR",   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins, fake_pT_bins_bmu, fake_N_closeLJpT_bins, fake_closeLJpT_bins_bmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_medDR",   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins, fake_pT_bins_bmu, fake_N_closeLJpT_bins, fake_closeLJpT_bins_bmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_highDR",  "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins, fake_pT_bins_bmu, fake_N_closeLJpT_bins, fake_closeLJpT_bins_bmu);

	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR",   	   "; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_bmu, N_DR_bins, DR_bins_bmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_lowDPhi",  "; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_bmu, N_DR_bins, DR_bins_bmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_highDPhi", "; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_bmu, N_DR_bins, DR_bins_bmu);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_met",	   	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_bmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_lowDR", 	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_bmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_medDR",	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_bmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_highDR",	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_bmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_met_highLepPt", "; min #Delta R(lep, jet); MET [GeV]", N_DR_bins, DR_bins_bmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_met_lowLepPt",  "; min #Delta R(lep, jet); MET [GeV]", N_DR_bins, DR_bins_bmu, metBinN, metBin);
	
	//ljet
	m_hSvc.create1D(suffix+"fake_ljet_pt",    	"; Pt of largeR jet [GeV]; Events", 100, 0, 1000);
	m_hSvc.create1D(suffix+"fake_ljet_m",    	"; mass of largeR jet [GeV]; Events", 50, 0, 200);
	m_hSvc.create1D(suffix+"fake_ljet_tau32", 	"; #tau_{32}; Events", 20, 0, 1);
	m_hSvc.create1D(suffix+"fake_ljet_tau32_wta", 	"; #tau_{32} wta; Events", 20, 0, 1);
	
     }
  }else{
  
    if (m_electron){
        eff_N_pT_bins 			= sizeof(eff_pT_bins_re)/sizeof(double) -1;
	N_DR_bins 			= sizeof(DR_bins_re)/sizeof(double) -1;	
        fake_N_pT_bins 			= sizeof(fake_pT_bins_re)/sizeof(double) -1;
	fake_N_closeLJpT_bins	 	= sizeof(fake_closeLJpT_bins_re)/sizeof(double) -1;
	fake_N_leptEta_bins 		= sizeof(fake_leptEta_bins_re)/sizeof(double) -1;

	fake_N_pT_bins_re_lowDR 		= sizeof(fake_pT_bins_re_lowDR)/sizeof(double) -1;
	fake_N_closeLJpT_bins_re_lowDR 		= sizeof(fake_closeLJpT_bins_re_lowDR)/sizeof(double) -1;
        fake_N_pT_bins_re_medDR 		= sizeof(fake_pT_bins_re_medDR)/sizeof(double) -1;
	fake_N_closeLJpT_bins_re_medDR 		= sizeof(fake_closeLJpT_bins_re_medDR)/sizeof(double) -1;
	fake_N_pT_bins_re_highDR 		= sizeof(fake_pT_bins_re_highDR)/sizeof(double) -1;
	fake_N_closeLJpT_bins_re_highDR 	= sizeof(fake_closeLJpT_bins_re_highDR)/sizeof(double) -1;

       	//Eff
 	m_hSvc.create2DVar(suffix+"eff_LepPt_DR", 	"; Pt of Matched lepton [GeV]; min #Delta R(lep, jet)", eff_N_pT_bins, eff_pT_bins_re, N_DR_bins, DR_bins_re);
	m_hSvc.create1DVar(suffix+"eff_lepPt",  	"; lepton p_{T} [GeV]; Events", eff_N_pT_bins, eff_pT_bins_re);
        m_hSvc.create1DVar(suffix+"eff_minDeltaR",  	"; min #Delta R(lep, jet); Events", N_DR_bins, DR_bins_re);
	
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi",	        "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_re, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi_lowDR",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_re, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi_highDR",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_re, CdPhiBinN, CdPhiBin);
        m_hSvc.create2DVar(suffix+"eff_lepPt_topoetcone",   	"; Pt of lepton [GeV]; topoetcone20 [GeV]", fake_N_pT_bins, fake_pT_bins_re, topoetconeBinN, topoetconeBin);
	m_hSvc.create1DVar(suffix+"eff_topoetcone_effBins",    "; topoetcone [GeV]; Events", topoetconeBinN, topoetconeBin); 

  	//Fake 1D     
  	m_hSvc.create1DVar(suffix+"fake_lepPt_effBins",   		"; Pt of lepton [GeV]; Eff", fake_N_pT_bins, fake_pT_bins_re);
	m_hSvc.create1DVar(suffix+"fake_lepEta_effBins",  		"; |#eta| of lepton; Eff", fake_N_leptEta_bins, fake_leptEta_bins_re);
	m_hSvc.create1DVar(suffix+"fake_minDeltaR_effBins", 		"; min #Delta R(lep, jet); Eff", N_DR_bins, DR_bins_re);
	m_hSvc.create1DVar(suffix+"fake_minDeltaR_tjet_effBins", 	"; min #Delta R(lep, tjet); Eff", N_DR_bins, DR_bins_re);
	m_hSvc.create1DVar(suffix+"fake_closJetPt_effBins", 		"; Pt of closest jet to lepton [GeV]; Eff", fake_N_closeLJpT_bins, fake_closeLJpT_bins_re);
	
	//Fake 2D
	m_hSvc.create2DVar(suffix+"fake_lepPt_lepEta",   	   "; Pt of lepton [GeV]; |#eta| of lepton", fake_N_pT_bins, fake_pT_bins_re, fake_N_leptEta_bins, fake_leptEta_bins_re);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt",   	   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]", fake_N_pT_bins, fake_pT_bins_re, fake_N_closeLJpT_bins, fake_closeLJpT_bins_re);
	m_hSvc.create2DVar(suffix+"fake_closJetPt_cosDPhi", 	   "; Pt of closest jet to lepton [GeV]; Cos( #Delta #phi(met, lept) )",   fake_N_closeLJpT_bins, fake_closeLJpT_bins_re, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi",       "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_re, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi_highLepPt", "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_re, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi_lowLepPt",  "; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_re, CdPhiBinN, CdPhiBin);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi",   	  "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_re, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowEta",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", pT_binsN_re_lowEta, fake_pT_bins_re_lowEta, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highEta",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", pT_binsN_re_highEta, fake_pT_bins_re_highEta, CdPhiBinN, CdPhiBin);

	m_hSvc.create2DVar(suffix+"fake_lepPt_topoetcone",   	"; Pt of lepton [GeV]; topoetcone20 [GeV]", fake_N_pT_bins, fake_pT_bins_re, topoetconeBinN, topoetconeBin);
        m_hSvc.create2DVar(suffix+"fake_lepPt_topoetcone_lowDR",      "; Pt of lepton [GeV]; topoetcone20 [GeV]", fake_N_pT_bins, fake_pT_bins_re, topoetconeBinN, topoetconeBin);
        m_hSvc.create2DVar(suffix+"fake_lepPt_topoetcone_highDR",      "; Pt of lepton [GeV]; topoetcone20 [GeV]", fake_N_pT_bins, fake_pT_bins_re, topoetconeBinN, topoetconeBin);
	m_hSvc.create1DVar(suffix+"fake_topoetcone_effBins",    "; topoetcone [GeV]; Events", topoetconeBinN, topoetconeBin); 
	
	m_hSvc.create1DVar(suffix+"fake_lepPt_highDR",            "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_re);
	m_hSvc.create1D(suffix+"fake_minDeltaR_highDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_highDR",	  "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	m_hSvc.create1DVar(suffix+"fake_lepPt_medDR",             "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_re);
	m_hSvc.create1D(suffix+"fake_minDeltaR_medDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_medDR",    "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	m_hSvc.create1DVar(suffix+"fake_lepPt_lowDR",             "; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_re);
	m_hSvc.create1D(suffix+"fake_minDeltaR_lowDR",      	  "; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_lowDR",    "; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowDR",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )",  fake_N_pT_bins_re_lowDR,  fake_pT_bins_re_lowDR,  CdPhi_lDR_BinN, CdPhi_lDR_Bin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_medDR",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )",  fake_N_pT_bins_re_medDR,  fake_pT_bins_re_medDR,  CdPhi_mDR_BinN, CdPhi_mDR_Bin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highDR",  "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )",  fake_N_pT_bins_re_highDR, fake_pT_bins_re_highDR, CdPhi_hDR_BinN, CdPhi_hDR_Bin);

	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowmWt",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )",  fake_N_pT_bins_re_lowDR,  fake_pT_bins_re_lowDR,  CdPhi_lDR_BinN, CdPhi_lDR_Bin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_medmWt",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )",  fake_N_pT_bins_re_highDR, fake_pT_bins_re_highDR, CdPhi_hDR_BinN, CdPhi_hDR_Bin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highmWt",  "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )",  fake_N_pT_bins_re_highDR, fake_pT_bins_re_highDR, CdPhi_hDR_BinN, CdPhi_hDR_Bin);

        m_hSvc.create1DVar(suffix+"fake_lepPt_highmWt",  	        "; lepton p_{T} [GeV]; Events", fake_N_pT_bins_re_lowDR,  fake_pT_bins_re_lowDR);
        m_hSvc.create1DVar(suffix+"fake_lepPt_medmWt",  	        "; lepton p_{T} [GeV]; Events", fake_N_pT_bins_re_medDR,  fake_pT_bins_re_medDR);
        m_hSvc.create1DVar(suffix+"fake_lepPt_lowmWt",  	        "; lepton p_{T} [GeV]; Events", fake_N_pT_bins_re_highDR, fake_pT_bins_re_highDR);

        m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_lowDR",   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins_re_lowDR,  fake_pT_bins_re_lowDR,  fake_N_closeLJpT_bins_re_lowDR, fake_closeLJpT_bins_re_lowDR);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_medDR",   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins_re_medDR,  fake_pT_bins_re_medDR,  fake_N_closeLJpT_bins_re_medDR, fake_closeLJpT_bins_re_medDR);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_highDR",  "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins_re_highDR, fake_pT_bins_re_highDR, fake_N_closeLJpT_bins_re_highDR, fake_closeLJpT_bins_re_highDR);

	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_lowEta",   "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins_re_lowDR, fake_pT_bins_re_lowDR, fake_N_closeLJpT_bins_re_lowDR, fake_closeLJpT_bins_re_lowDR);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_highEta",  "; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins_re_highDR, fake_pT_bins_re_highDR, fake_N_closeLJpT_bins_re_highDR, fake_closeLJpT_bins_re_highDR);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR",   	  "; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_re, N_DR_bins, DR_bins_re);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_lowDPhi", "; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_re, N_DR_bins, DR_bins_re);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_highDPhi","; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_re, N_DR_bins, DR_bins_re);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_met",	   	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_re, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_lowDR", 	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_re, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_medDR",	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_re, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_highDR",	   "; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_re, metBinN, metBin);

	//chi2
	m_hSvc.create1D(suffix+"fake_chi2", "; log(#chi^{2}) ; Events", 50, -3, 7);

     }else{
     
        eff_N_pT_bins 		= sizeof(eff_pT_bins_rmu)/sizeof(double) -1;
	N_DR_bins 		= sizeof(DR_bins_rmu)/sizeof(double) -1;
        N_DR_bins_mbin          = sizeof(DR_bins_rmu_mbin)/sizeof(double) -1;	
        N_DR_bins_sbin          = sizeof(DR_bins_rmu_sbin)/sizeof(double) -1;
	fake_N_pT_bins 		= sizeof(fake_pT_bins_rmu)/sizeof(double) -1;
	fake_N_pT_bins_lDR 	= sizeof(fake_pT_bins_lDR_rmu)/sizeof(double) -1;
	fake_N_pT_bins_mDR 	= sizeof(fake_pT_bins_mDR_rmu)/sizeof(double) -1;
	fake_N_pT_bins_hDR 	= sizeof(fake_pT_bins_hDR_rmu)/sizeof(double) -1;
	fake_N_leptEta_bins 	= sizeof(fake_leptEta_bins_rmu)/sizeof(double) -1;
	fake_N_closeLJpT_bins 	= sizeof(fake_closeLJpT_bins_rmu)/sizeof(double) -1;
     
     	//Eff
  	m_hSvc.create2DVar(suffix+"eff_LepPt_DR", 	"; Pt of Matched lepton [GeV]; min #Delta R(lep, jet)", eff_N_pT_bins, eff_pT_bins_rmu, N_DR_bins, DR_bins_rmu);
        m_hSvc.create2DVar(suffix+"eff_LepPt_DR_mbin",       "; Pt of Matched lepton [GeV]; min #Delta R(lep, jet)", eff_N_pT_bins, eff_pT_bins_rmu, N_DR_bins_mbin, DR_bins_rmu_mbin);
        m_hSvc.create2DVar(suffix+"eff_LepPt_DR_sbin",       "; Pt of Matched lepton [GeV]; min #Delta R(lep, jet)", eff_N_pT_bins, eff_pT_bins_rmu, N_DR_bins_sbin, DR_bins_rmu_sbin);
        m_hSvc.create1DVar(suffix+"eff_lepPt",  	"; lepton p_{T} [GeV]; Events", eff_N_pT_bins, eff_pT_bins_rmu);
        m_hSvc.create1DVar(suffix+"eff_minDeltaR",  	"; min #Delta R(lep, jet); Events", N_DR_bins, DR_bins_rmu);
	
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi",	        "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_rmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi_lowDR",    "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_rmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"eff_lepPt_cosDPhi_highDR",   "; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", eff_N_pT_bins, eff_pT_bins_rmu, CdPhiBinN, CdPhiBin);
	
        m_hSvc.create2DVar(suffix+"eff_lepPt_topoetcone",      "; Pt of lepton [GeV]; topoetcone20 [GeV]", eff_N_pT_bins, eff_pT_bins_rmu, topoetconeBinN, topoetconeBin);
        m_hSvc.create2DVar(suffix+"eff_lepPt_topoetcone_lowDR",      "; Pt of lepton [GeV]; topoetcone20 [GeV]", eff_N_pT_bins, eff_pT_bins_rmu, topoetconeBinN, topoetconeBin);
        m_hSvc.create2DVar(suffix+"eff_lepPt_topoetcone_medDR",      "; Pt of lepton [GeV]; topoetcone20 [GeV]", eff_N_pT_bins, eff_pT_bins_rmu, topoetconeBinN, topoetconeBin);
        m_hSvc.create2DVar(suffix+"eff_lepPt_topoetcone_highDR",      "; Pt of lepton [GeV]; topoetcone20 [GeV]", eff_N_pT_bins, eff_pT_bins_rmu, topoetconeBinN, topoetconeBin);
        m_hSvc.create1DVar(suffix+"eff_topoetcone_effBins",    "; topoetcone [GeV]; Events", topoetconeBinN, topoetconeBin);


	
	
  	//Fake 1D     
  	m_hSvc.create1DVar(suffix+"fake_lepPt_effBins",  		"; Pt of lepton [GeV]; Eff", fake_N_pT_bins, fake_pT_bins_rmu);
	m_hSvc.create1DVar(suffix+"fake_minDeltaR_effBins",    		"; min #Delta R(lep, jet); Eff", N_DR_bins, DR_bins_rmu);
	m_hSvc.create1DVar(suffix+"fake_minDeltaR_tjet_effBins",    	"; min #Delta R(lep, tjet); Eff", N_DR_bins, DR_bins_rmu);
	m_hSvc.create1DVar(suffix+"fake_lepEta_effBins",  		"; |#eta| of lepton; Eff", fake_N_leptEta_bins, fake_leptEta_bins_rmu);
	m_hSvc.create1DVar(suffix+"fake_closJetPt_effBins",    		"; Pt of closest jet to lepton [GeV]; Eff", fake_N_closeLJpT_bins, fake_closeLJpT_bins_rmu);

	m_hSvc.create1D(suffix+"fake_minDeltaR_highDR",      	  	"; min #Delta R(lep, jet); Events", 70, 0, 7);	
	m_hSvc.create1DVar(suffix+"fake_lepPt_highDR",            	"; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_rmu);
	m_hSvc.create1DVar(suffix+"fake_cos_metPhi_lepPhi_highDR",	"; Cos( #Delta #phi(met, lept) ); Events", CdPhiBinN, CdPhiBin);
	m_hSvc.create1DVar(suffix+"fake_mwt_met_highDR", 	  	"; MET + MWT [GeV]; Events", mwtmetBinN, mwtmetBin);
	m_hSvc.create1DVar(suffix+"fake_met_highDR", 	  	  	"; MET [GeV]; Events", metBinN, metBin);

	m_hSvc.create1D(suffix+"fake_minDeltaR_lowDR",      	  	"; min #Delta R(lep, jet); Events", 70, 0, 7);
	m_hSvc.create1DVar(suffix+"fake_lepPt_lowDR",             	"; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_rmu);
	m_hSvc.create1DVar(suffix+"fake_cos_metPhi_lepPhi_lowDR", 	"; Cos( #Delta #phi(met, lept) ); Events", CdPhiBinN, CdPhiBin);
	m_hSvc.create1DVar(suffix+"fake_mwt_met_lowDR", 	  	"; MET + MWT [GeV]; Events", mwtmetBinN, mwtmetBin);
	m_hSvc.create1DVar(suffix+"fake_met_lowDR", 	  	  	"; MET [GeV]; Events", metBinN, metBin);

	//m_hSvc.create1DVar(suffix+"fake_lepPt_medDR",           	"; lepton p_{T} [GeV]; Events", fake_N_pT_bins, fake_pT_bins_rmu);
	//m_hSvc.create1D(suffix+"fake_minDeltaR_medDR",      	  	"; min #Delta R(lep, jet); Events", 70, 0, 7);
	//m_hSvc.create1D(suffix+"fake_cos_metPhi_lepPhi_medDR",    	"; Cos( #Delta #phi(met, lept) ); Events", 100, -1, 1);
		
	//Fake 2D
	m_hSvc.create2DVar(suffix+"fake_lepPt_lepEta",   	  	"; Pt of lepton [GeV]; |#eta| of lepton", fake_N_pT_bins, fake_pT_bins_rmu, fake_N_leptEta_bins, fake_leptEta_bins_rmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt",   	  	"; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]", fake_N_pT_bins, fake_pT_bins_rmu, fake_N_closeLJpT_bins, fake_closeLJpT_bins_rmu);
	m_hSvc.create2DVar(suffix+"fake_closJetPt_cosDPhi", 	  	"; Pt of closest jet to lepton [GeV]; Cos( #Delta #phi(met, lept) )",   fake_N_closeLJpT_bins, fake_closeLJpT_bins_rmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi",      	"; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_rmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi_highLepPt", 	"; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_rmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_cosDPhi_lowLepPt",  	"; min #Delta R(lep, jet); Cos( #Delta #phi(met, lept) )", N_DR_bins, DR_bins_rmu, CdPhiBinN, CdPhiBin);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi",   	  	"; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_rmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowEta",    	"; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_rmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highEta",   	"; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_rmu, CdPhiBinN, CdPhiBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_lowDR",     	"; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins_lDR, fake_pT_bins_lDR_rmu, CdPhi_lDR_BinN, CdPhi_lDR_Bin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_highDR",    	"; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_rmu, CdPhi_hDR_BinN, CdPhi_hDR_Bin);
	//m_hSvc.create2DVar(suffix+"fake_lepPt_cosDPhi_medDR",     	"; Pt of lepton [GeV]; Cos( #Delta #phi(met, lept) )", fake_N_pT_bins, fake_pT_bins_rmu, CdPhi_mDR_BinN, CdPhi_mDR_Bin);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_lowDR",   	"; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins, fake_pT_bins_rmu, fake_N_closeLJpT_bins, fake_closeLJpT_bins_rmu);
	//m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_medDR",   	"; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins, fake_pT_bins_rmu, fake_N_closeLJpT_bins, fake_closeLJpT_bins_rmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_closJetPt_highDR",  	"; Pt of lepton [GeV]; Pt of closest jet to lepton [GeV]",  fake_N_pT_bins, fake_pT_bins_rmu, fake_N_closeLJpT_bins, fake_closeLJpT_bins_rmu);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_lowDR", 	"; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_rmu, N_DR_bins, DR_bins_rmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_highDR",	"; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_rmu, N_DR_bins, DR_bins_rmu);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR",   	  	"; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_rmu, N_DR_bins, DR_bins_rmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_lowDPhi", 	"; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_rmu, N_DR_bins, DR_bins_rmu);
	m_hSvc.create2DVar(suffix+"fake_lepPt_minDeltaR_highDPhi",	"; Pt of lepton [GeV]; min #Delta R(lep, jet)", fake_N_pT_bins, fake_pT_bins_rmu, N_DR_bins, DR_bins_rmu);
	
	m_hSvc.create2DVar(suffix+"fake_lepPt_met",	   	   	"; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_rmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_lowDR", 	   	"; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_rmu, metBinN, metBin);
	//m_hSvc.create2DVar(suffix+"fake_lepPt_met_medDR",	   	"; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_rmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_lepPt_met_highDR",	   	"; Pt of lepton [GeV]; MET [GeV] ", fake_N_pT_bins, fake_pT_bins_rmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_met_highLepPt", 	"; min #Delta R(lep, jet); MET [GeV]", N_DR_bins, DR_bins_rmu, metBinN, metBin);
	m_hSvc.create2DVar(suffix+"fake_minDeltaR_met_lowLepPt",  	"; min #Delta R(lep, jet); MET [GeV]", N_DR_bins, DR_bins_rmu, metBinN, metBin);
	
	//chi2
	m_hSvc.create1D(suffix+"fake_chi2", "; log(#chi^{2}) ; Events", 50, -3, 7);
	m_hSvc.create2DVar(suffix+"fake_lepPt_topoetcone",   	"; Pt of lepton [GeV]; topoetcone20 [GeV]", fake_N_pT_bins, fake_pT_bins_rmu, topoetconeBinN, topoetconeBin);
        m_hSvc.create2DVar(suffix+"fake_lepPt_topoetcone_lowDR",      "; Pt of lepton [GeV]; topoetcone20 [GeV]", fake_N_pT_bins_lDR, fake_pT_bins_lDR_rmu, topoetconeBinN, topoetconeBin);
        m_hSvc.create2DVar(suffix+"fake_lepPt_topoetcone_highDR",      "; Pt of lepton [GeV]; topoetcone20 [GeV]", fake_N_pT_bins, fake_pT_bins_rmu, topoetconeBinN, topoetconeBin);
	m_hSvc.create1DVar(suffix+"fake_topoetcone_effBins",    "; topoetcone [GeV]; Events", topoetconeBinN, topoetconeBin);


     }	
  }//if (m_boosted) 







}//AnaTtresQCD::IniHistograms


