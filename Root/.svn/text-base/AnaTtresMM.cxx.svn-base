/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include <sstream>
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresMM.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"

AnaTtresMM::AnaTtresMM(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList, int dsid)
  : Analysis(filename, systList), m_electron(electron), m_boosted(boosted), m_dsid(dsid),
    m_neutrinoBuilder("MeV"), m_chi2("MeV") {

  m_chi2.Init(TtresChi2::DATA2015_MC15C);

  std::string suffix = "";
  std::string btag0 = "btag0_";
  std::string btag1 = "btag1_";
  IniHistograms(suffix);
  IniHistograms(btag0);
  IniHistograms(btag1);
  
  
  
  /*
  std::string name_hDR = "hDR_";
  IniHistograms(name_hDR);
  
  std::string name_mDR = "mDR_";
  IniHistograms(name_mDR);
  
  std::string name_lDR = "lDR_";
  IniHistograms(name_lDR);
  
  std::string name_hEta = "hEta_";
  IniHistograms(name_hEta);
  
  std::string name_lEta = "lEta_";
  IniHistograms(name_lEta);
  */
}//AnaTtresMM

AnaTtresMM::~AnaTtresMM() {
}

void AnaTtresMM::run(const Event &evt, double weight, const std::string &suffix){
}

//***********************************
//**Multijet background generation **
//***********************************

// ----- begin functions with 2015 configuration ---------- //

void AnaTtresMM::runMatrixMethod_QCDCR2j_2015(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;

  if (m_boosted)	return;
  else			if(evt.jet().size()<2)	return;

  HistogramService *h = &m_hSvc;

  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
               
  }//m_electron


 // ----------------- GEN LEPTON MATCHING -------------------- //
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

  if (l.DeltaR(lgen) > 0.2) return;

  } // if(m_dsid == 410000)
 // ------------- GEN LEPTON MATCHING --------- //


  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET>20) || (MET+mWt)>60)	return;
      //if( (MET<20) || (MET+mWt)>60)	return;
  }else{
      if( (MET>20) || (MET+mWt)>60)	return;
      if(fabs(d0sig)>3)			return;
  }//if

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

  int nTrkBtagged = 0;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
          nTrkBtagged += 1;

  }//for

  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);
    
    
  
}//runMatrixMethod_QCDCR2j_2015


void AnaTtresMM::runMatrixMethod_QCDCR4j_2015(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;

  if (m_boosted)       if(evt.largeJet().size()<1)  return;
  if (!m_boosted)      if(evt.jet().size()<4)	return;

  HistogramService *h = &m_hSvc;

  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
               
  }//m_electron


  // ----------------- GEN LEPTON MATCHING -------------------- //
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

  if (l.DeltaR(lgen) > 0.2) return;

  } // if(m_dsid == 410000)
 // ------------- GEN LEPTON MATCHING --------- //


  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET>20) || (MET+mWt)>60)	return;
      //if( (MET<20) || (MET+mWt)>60)	return;
  }else{
      if( (MET>20) || (MET+mWt)>60)	return;
      if(fabs(d0sig)>3)			return;
  }//if

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

  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }

  int nTrkBtagged = 0;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
          nTrkBtagged += 1;

  }//for

  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);
    
    
  
}//runMatrixMethod_QCDCR4j_2015




void AnaTtresMM::runMatrixMethod_QCDVR1_2j_2015(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;

  if (m_boosted)	return;
  else			if(evt.jet().size()<2)	return;

  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
           
  }//m_electron

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET>20) || (MET+mWt)<60)		return; //VR1
      //if( (MET<20) || (MET+mWt)>60)		return; //VR2
  }else{
      if( (MET>20) || (MET+mWt)<60)		return;
      //if( (MET<20) || (MET+mWt)>60)		return; //VR2
      if(fabs(d0sig)>3)      	return;
  }//if
  
  //nB-tagged jets 
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       	if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          	nTrkBtagged += 1;
  }
  
  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);

}//AnaTtresMM::runMatrixMethod_QCDVR1_2j_2015


void AnaTtresMM::runMatrixMethod_QCDVR2_2j_2015(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;

  if (m_boosted)	return;
  else			if(evt.jet().size()<2)	return;

  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
           
  }//m_electron

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      //if( (MET>20) || (MET+mWt)<60)		return; //VR1
      if( (MET<20) || (MET+mWt)>60)		return; //VR2
  }else{
      //if( (MET>20) || (MET+mWt)<60)		return;
      if( (MET<20) || (MET+mWt)>60)		return; //VR2
      if(fabs(d0sig)>3)      	return;
  }//if
  
  //nB-tagged jets 
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       	if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          	nTrkBtagged += 1;
  }
  
  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);

}//AnaTtresMM::runMatrixMethod_QCDVR2_2j_2015


void AnaTtresMM::runMatrixMethod_WjetsCR2j_2015(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsWCR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;

  if (m_boosted)	return;  
  else			if(evt.jet().size()<2)	return;

  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
               
  }//m_electron

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET<20) || (MET+mWt)<60)		return;
  }else{
      if( (MET<20) || (MET+mWt)<60)		return;
      if(fabs(d0sig)>3)			      	return;
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
 //if (nTrkBtagged!=0)	return;	
 
  GetHistograms(evt, weight, "", suffix);

}//runMatrixMethod_WjetsCR2j_2015


void AnaTtresMM::runMatrixMethod_SR4j_2015(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsWCR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;

  int nLargeRjets(0);
  if (m_boosted){
  
    size_t goodljet_idx = 0;
    for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx) {
      if (evt.largeJet()[goodljet_idx].good())
        if (evt.largeJet()[goodljet_idx].mom().Pt()*1e-3 > 300)
          nLargeRjets += 1;
    }
    
    if (nLargeRjets==0)	return;
  
  } else	if(evt.jet().size()<4)	return;
  
  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
               
  }//m_electron

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET<20) || (MET+mWt)<60)		return;
  }else{
      if( (MET<20) || (MET+mWt)<60)		return;
      if(fabs(d0sig)>3)			      	return;
  }//if
  
  GetHistograms(evt, weight, "", suffix);

}//runMatrixMethod_SR4j_2015

void AnaTtresMM::runMatrixMethod_QCDVR1_4j_2015(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;
/*
  int nLargeRjets(0);
  if (m_boosted){
  
    size_t goodljet_idx = 0;
    for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx) {
      if (evt.largeJet()[goodljet_idx].good())
        if (evt.largeJet()[goodljet_idx].mom().Pt()*1e-3 > 300)
          nLargeRjets += 1;
    }
    
    if (nLargeRjets==0)	return;
  
  } else	if(evt.jet().size()<4)	return;
*/
  if (m_boosted)       if(evt.largeJet().size()<1)  return;
  if(!m_boosted)       if(evt.jet().size()<4)  return;

  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
               
  }//m_electron

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET>20) || (MET+mWt)<60)		return; //Validation region 1
      //if( (MET<20) || (MET+mWt)>60)		return; //Validation region 2
  }else{
      if( (MET > 20) || (MET+mWt)<60)             return; // VR1
      //if( (MET < 20) || (MET+mWt)>60)             return; // VR2
      if(fabs(d0sig)>3)                          return;
      //
      
  }//if

  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }

  int nTrkBtagged = 0;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
          nTrkBtagged += 1;

  }//for
  //std::cout<<"REACHED HERE, wt = "<<weight<<std::endl;
  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);


}//AnaTtresMM::runMatrixMethod_QCDVR1_4j_2015


void AnaTtresMM::runMatrixMethod_QCDVR2_4j_2015(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bmujetsQCDCR_2015")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rmujetsQCDCR_2015")))
      return;
/*
  int nLargeRjets(0);
  if (m_boosted){
  
    size_t goodljet_idx = 0;
    for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx) {
      if (evt.largeJet()[goodljet_idx].good())
        if (evt.largeJet()[goodljet_idx].mom().Pt()*1e-3 > 300)
          nLargeRjets += 1;
    }
    
    if (nLargeRjets==0)	return;
  
  } else	if(evt.jet().size()<4)	return;
*/
  if (m_boosted)       if(evt.largeJet().size()<1)  return;
  if(!m_boosted)       if(evt.jet().size()<4)  return;


  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
               
  }//m_electron

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      //if( (MET>20) || (MET+mWt)<60)		return; //Validation region 1
      if( (MET<20) || (MET+mWt)>60)		return; //Validation region 2
  }else{
      //if( (MET > 20) || (MET+mWt)<60)             return; // VR1
      if( (MET < 20) || (MET+mWt)>60)             return; // VR2
      if(fabs(d0sig)>3)                          return;
      //
      
  }//if

  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }

  int nTrkBtagged = 0;
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())
          nTrkBtagged += 1;

  }//for
  //std::cout<<"REACHED HERE, wt = "<<weight<<std::endl;
  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);


}//AnaTtresMM::runMatrixMethod_QCDVR2_4j_2015


// ----- end functions with 2015 configuration ---------- //

// ----- functions with 2016 version ---------- //

void AnaTtresMM::runMatrixMethod_QCDCR2j_2016(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2016") || evt.passes("bmujetsQCDCR_2016")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2016") || evt.passes("rmujetsQCDCR_2016")))
      return;

  if (!m_boosted)	if(evt.jet().size()<2)	return;

  HistogramService *h = &m_hSvc;

  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
    
           
  }//m_electron


  // ----------------- GEN LEPTON MATCHING -------------------- //
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

  if (l.DeltaR(lgen) > 0.2) return;

  } // if(m_dsid == 410000)
 // ------------- GEN LEPTON MATCHING --------- //


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
       
  TLorentzVector met = evt.met();
  bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L);

  float chi2Value = 1000000; // log10(1000000) = 6
  if(status)
  chi2Value = chi2ming1;



  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET>20) || (MET+mWt)>60)	return;
  }else{
      if( (MET>20) || (MET+mWt)>60)	return;
      if(fabs(d0sig)>3)			return;
  }//if

  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;	

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
   
  //if (closejl_deltaR>0.6)	GetHistograms(evt, weight, "hDR_", suffix);
  //else if (closejl_deltaR>0.4)	GetHistograms(evt, weight, "mDR_", suffix);
  //else				GetHistograms(evt, weight, "lDR_", suffix);

  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }
 

  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix); 
 
} // runMatrixMethod_QCDCR2j_2016



void AnaTtresMM::runMatrixMethod_QCDVR1_2j_2016(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2016") || evt.passes("bmujetsQCDCR_2016")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2016") || evt.passes("rmujetsQCDCR_2016")))
      return;

  if (!m_boosted)	if(evt.jet().size()<2)	return;

  HistogramService *h = &m_hSvc;

  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
           
  }//m_electron


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
       
  TLorentzVector met = evt.met();
  bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L);

  float chi2Value = 1000000; // log10(1000000) = 6
  if(status)
  chi2Value = chi2ming1;

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
 
  if(m_electron){
      if( (MET>20) && (MET+mWt)<60)		return; // VR1
      //if( (MET < 20) || (MET+mWt)>60)             return; // VR2
  }else{
      if( (MET > 20) || (MET+mWt)<60)             return; // VR1
      //if( (MET < 20) || (MET+mWt)>60)             return; // VR2
      if(fabs(d0sig)>3)                        return;
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;	

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

  //if (closejl_deltaR>0.6)	GetHistograms(evt, weight, "hDR_", suffix);
  //else if (closejl_deltaR>0.4)	GetHistograms(evt, weight, "mDR_", suffix);
  //else				GetHistograms(evt, weight, "lDR_", suffix);
  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }

  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0) 
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);
  
}//AnaTtresMM::runMatrixMethod_QCDVR1_2j_2016



void AnaTtresMM::runMatrixMethod_QCDVR2_2j_2016(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2016") || evt.passes("bmujetsQCDCR_2016")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2016") || evt.passes("rmujetsQCDCR_2016")))
      return;

  if (!m_boosted)	if(evt.jet().size()<2)	return;

  HistogramService *h = &m_hSvc;

  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
           
  }//m_electron


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
       
  TLorentzVector met = evt.met();
  bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L);

  float chi2Value = 1000000; // log10(1000000) = 6
  if(status)
  chi2Value = chi2ming1;

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
 
  if(m_electron){
      //if( (MET>20) && (MET+mWt)<60)		return;
      if( (MET < 20) || (MET+mWt)>60)             return; // VR2
  }else{
      //if( (MET > 20) || (MET+mWt)<60)             return; // VR1
      if( (MET < 20) || (MET+mWt)>60)             return; // VR2
      if(fabs(d0sig)>3)                        return;
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;	

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

  //if (closejl_deltaR>0.6)	GetHistograms(evt, weight, "hDR_", suffix);
  //else if (closejl_deltaR>0.4)	GetHistograms(evt, weight, "mDR_", suffix);
  //else				GetHistograms(evt, weight, "lDR_", suffix);
  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }

  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0) 
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);
  
}//AnaTtresMM::runMatrixMethod_QCDVR2_2j_2016



void AnaTtresMM::runMatrixMethod_QCDSR2j_2016(const Event &evt, double weight, const std::string &suffix) {

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2015") || evt.passes("bejetsWCR_2016")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2015") || evt.passes("rejetsWCR_2016")))
      return;
  
  if (!m_boosted)	if(evt.jet().size()<2)	return;

  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
    /*
    //Muon trigers
    trig1 = evt.muon()[0].HLT_mu20_L1MU15(); //prescaled
    trig2 = evt.muon()[0].HLT_mu50();
    trig3 = evt.muon()[0].HLT_mu20_iloose_L1MU15();
    
    bool trig_prescaled   = trig1;
    bool trig_unprescaled = trig2 || trig3;
    
    //if (isTight)
       if (trig_prescaled && !trig_unprescaled)	return;
     */      
  }//m_electron

  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET<20) || (MET+mWt)<60)		return;
  }else{
      if( (MET<20) || (MET+mWt)<60)		return;
      if(fabs(d0sig)>3)			      	return;
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;	

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

  //if (closejl_deltaR>0.6)	GetHistograms(evt, weight, "hDR_", suffix);
  //else if (closejl_deltaR>0.4)	GetHistograms(evt, weight, "mDR_", suffix);
  //else				GetHistograms(evt, weight, "lDR_", suffix);
  GetHistograms(evt, weight, "", suffix);
}//runMatrixMethod_QCDSR2j_2016


void AnaTtresMM::runMatrixMethod_QCDCR4j_2016(const Event &evt, double weight, const std::string &suffix) {

  //if(evt.runNumber_or_RandomRunNumber() > 302872) return;

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2016") || evt.passes("bmujetsQCDCR_2016")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2016") || evt.passes("rmujetsQCDCR_2016")))
      return;

  if (m_boosted)       if(evt.largeJet().size()<1)  return;
  if (!m_boosted)	if(evt.jet().size()<4)	return;

  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
    
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
           
  }//m_electron


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
       
  TLorentzVector met = evt.met();
  bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L);

  float chi2Value = 1000000; // log10(1000000) = 6
  if(status)
  chi2Value = chi2ming1;

  // ----------------- GEN LEPTON MATCHING -------------------- //
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

  if (l.DeltaR(lgen) > 0.2) return;

  } // if(m_dsid == 410000)
 // ------------- GEN LEPTON MATCHING --------- //


  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET>20) || (MET+mWt)>60)	return;
  }else{
      if( (MET>20) || (MET+mWt)>60)	return;
      if(fabs(d0sig)>3)			return;
      //if(log10(chi2Value) < 1.5 || log10(chi2Value) > 6)        return;
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;	

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

  //if (closejl_deltaR>0.6)	GetHistograms(evt, weight, "hDR_", suffix);
  //else if (closejl_deltaR>0.4)	GetHistograms(evt, weight, "mDR_", suffix);
  //else				GetHistograms(evt, weight, "lDR_", suffix);
  //std::cout<<"Real weight = "<<weight<<std::endl; 
  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }
  //std::cout<<"FILLING HISTO, weight = "<<weight<<std::endl;
  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);
  


}//AnaTtresMM::runMatrixMethod_QCDCR4j_2016


void AnaTtresMM::runMatrixMethod_QCDVR1_4j_2016(const Event &evt, double weight, const std::string &suffix) {

  //if(evt.runNumber_or_RandomRunNumber() > 302872) return; 

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2016") || evt.passes("bmujetsQCDCR_2016")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2016") || evt.passes("rmujetsQCDCR_2016")))
      return;

  if (m_boosted)       if(evt.largeJet().size()<1)  return;
  if (!m_boosted)	if(evt.jet().size()<4)	return;

  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
  
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
    
  }//m_electron


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
       
  TLorentzVector met = evt.met();
  bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L);

  float chi2Value = 1000000; // log10(1000000) = 6
  if(status)
  chi2Value = chi2ming1;


  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      if( (MET>20) || (MET+mWt)<60)		return; // VR1
      //if( (MET < 20) || (MET+mWt)>60)             return; // VR2
  }else{
      if( (MET > 20) || (MET+mWt)<60)             return; // VR1
      //if( (MET < 20) || (MET+mWt)>60)             return; // VR2
      if(fabs(d0sig)>3)                          return; 
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;	

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

  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }


  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);

}//runMatrixMethod_QCDVR1_4j_2016


void AnaTtresMM::runMatrixMethod_QCDVR2_4j_2016(const Event &evt, double weight, const std::string &suffix) {

  //if(evt.runNumber_or_RandomRunNumber() > 302872) return;

  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted) {
    if (!(evt.passes("bejetsIncluR_2016") || evt.passes("bmujetsQCDCR_2016")))
      return;
  }

  if (!m_boosted)
    if (!(evt.passes("rejetsIncluR_2016") || evt.passes("rmujetsQCDCR_2016")))
      return;

  if (m_boosted)       if(evt.largeJet().size()<1)  return;
  if (!m_boosted)	if(evt.jet().size()<4)	return;

  HistogramService *h = &m_hSvc;
  
  bool isTight = false;
  float d0sig(99);
  
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
    isTight = evt.electron()[0].isTightPP();
    d0sig = evt.electron()[0].sd0();

  } else {
    l = evt.muon()[0].mom();
    isTight = evt.muon()[0].isTight();
    d0sig = evt.muon()[0].sd0();
    
  }//m_electron


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
       
  TLorentzVector met = evt.met();
  bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L);

  float chi2Value = 1000000; // log10(1000000) = 6
  if(status)
  chi2Value = chi2ming1;


  float mWt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l)) ))*1e-3; 
  float MET = evt.met().Perp()*1e-3;
  
  if(m_electron){
      //if( (MET>20) || (MET+mWt)<60)		return; // VR1
      if( (MET < 20) || (MET+mWt)>60)             return; // VR2
  }else{
      //if( (MET > 20) || (MET+mWt)<60)             return; // VR1
      if( (MET < 20) || (MET+mWt)>60)             return; // VR2
      if(fabs(d0sig)>3)                          return; 
  }//if
  
  int nTrkBtagged = 0; 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx){
       if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk())	     
          nTrkBtagged += 1;
  
  }//for 
  
  //if (nTrkBtagged!=0)	return;	

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

  if(std::isnan(weight)) {
  std::cout<<"Weight is Nan, making it  0"<<std::endl;
  weight = 0;
  }

  if(std::isinf(weight)) {
  std::cout<<"Weight is Inf, making it  0"<<std::endl;
  weight = 0;
  }


  GetHistograms(evt, weight, "", suffix);
  if(nTrkBtagged == 0)
  GetHistograms(evt, weight, "btag0_",suffix );
  if(nTrkBtagged >= 1)
  GetHistograms(evt, weight, "btag1_", suffix);

}//runMatrixMethod_QCDVR2_4j_2016




// ===== end of v-2016 functions -------------- //




void AnaTtresMM::GetHistograms(const Event &evt, double weight, const std::string &prefix, const std::string &suffix){

  HistogramService *h = &m_hSvc;
  
  //lepton
  TLorentzVector l;   
  float d0sig(0);
  float ptvarcone(-1);
  float topoetcone(-1);

  if (m_electron) {
     l  = evt.electron()[0].mom();
     d0sig = evt.electron()[0].sd0();
     ptvarcone	= evt.electron()[0].ptvarcone20();
     topoetcone	= evt.electron()[0].topoetcone20();
  } else{			
     l  = evt.muon()[0].mom();
     d0sig = evt.muon()[0].sd0();
     ptvarcone	= evt.muon()[0].ptvarcone30();
     topoetcone	= evt.muon()[0].topoetcone20();

  }//m_electron 

  h->h1D(prefix+"lepPt", "", 		suffix)->Fill(l.Perp()*1e-3, weight);
  h->h1D(prefix+"lepPt_effBins", "", 	suffix)->Fill(l.Perp()*1e-3, weight);
  h->h1D(prefix+"lepEta", "", 		suffix)->Fill(l.Eta(), weight);
  h->h1D(prefix+"lepEta_effBins", "", 	suffix)->Fill(fabs(l.Eta()), weight);
  h->h1D(prefix+"lepPhi", "", 		suffix)->Fill(l.Phi(), weight);

  h->h1D(prefix+"MET_phi", "", 	suffix)->Fill(evt.met().Phi(), weight);
  h->h1D(prefix+"ptvarcone", 	    "", suffix)->Fill(ptvarcone*1e-3, weight);
  h->h1D(prefix+"topoetcone", 	    "", suffix)->Fill(topoetcone*1e-3, weight);


  const TLorentzVector &j = evt.jet()[0].mom();
  h->h1D(prefix+"leadJetPt", "", suffix)->Fill(j.Perp()*1e-3, weight);
    
  // for now
  int nJets = evt.jet().size(); //njets 
  h->h1D(prefix+"nJets", "", suffix)->Fill(nJets, weight);
  
  int nBtagged = 0; //nB-tagged jets 
  for (size_t bidx = 0; bidx < evt.jet().size(); ++bidx)
    if (evt.jet()[bidx].btag_mv2c20_70()){
      if(nBtagged==0)h->h1D(prefix+"leadbJetPt", "", suffix)->Fill(evt.jet()[bidx].mom().Perp()*1e-3, weight);
      nBtagged += 1;
    }
  h->h1D(prefix+"nBtagJets", "", suffix)->Fill(nBtagged, weight);

  std::vector<float> tjetPt_vector;
  int nTrkBtagged = 0; //nB-tagged jets 
  for (size_t bidx = 0; bidx < evt.tjet().size(); ++bidx) {
    tjetPt_vector.push_back(evt.tjet()[bidx].mom().Perp());
    if (evt.tjet()[bidx].btag_mv2c10_70_trk() && evt.tjet()[bidx].pass_trk()){
      if(nTrkBtagged==0)h->h1D(prefix+"leadTrkbJetPt", "", suffix)->Fill(evt.tjet()[bidx].mom().Perp()*1e-3, weight);
      nTrkBtagged += 1;
    }
  }
  h->h1D(prefix+"nTrkBtagJets", "", suffix)->Fill(nTrkBtagged, weight);

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
     std::string nameJet_m = prefix+"jet" + ss.str()+"_m";
     h->h1D(nameJet_m, "", suffix)->Fill(jetMass_vector[i]*1e-3, weight);
     
     std::string nameJet_pt = prefix+"jet" + ss.str()+"_pt";
     h->h1D(nameJet_pt, "", suffix)->Fill(jetPt_vector[i]*1e-3, weight);
     
     std::string nameJet_eta = prefix+"jet" + ss.str()+"_eta";
     h->h1D(nameJet_eta, "", suffix)->Fill(jetEta_vector[i], weight);
     
     std::string nameJet_phi = prefix+"jet" + ss.str()+"_phi";
     h->h1D(nameJet_phi, "", suffix)->Fill(jetPhi_vector[i], weight);    
  
  }//for
  
  float mtt = -1;
  
  float met = evt.met().Perp()*1e-3;
  
  //missing et
  h->h1D(prefix+"MET", 	 	"", suffix)->Fill(met, weight);
  h->h1D(prefix+"MET_effBins",  "", suffix)->Fill(met, weight);

  //transverse W mass
  float mwt = sqrt(2. * l.Perp() * evt.met().Perp() * (1. - cos(evt.met().DeltaPhi(l))))*1e-3;
  
  h->h1D(prefix+"mwt", 	 	"", suffix)->Fill(mwt, weight); 
  h->h1D(prefix+"mwt_effBins",  "", suffix)->Fill(mwt, weight);

  h->h1D(prefix+"mwt_met", "", suffix)->Fill(mwt+met, weight);
  
  h->h1D(prefix+"yields", "", suffix)->Fill(1, weight);
  
  //mu
  h->h1D(prefix+"mu", "", suffix)->Fill(evt.mu()*1.16, weight); 
  h->h1D(prefix+"mu_original", "", suffix)->Fill(evt.mu_original(), weight); 
    
  //npv
  h->h1D(prefix+"npv", "", suffix)->Fill(evt.npv(), weight);
  
  //z prosition of primary vertex
  h->h1D(prefix+"vtxz", "", suffix)->Fill(evt.vtxz(), weight);
  
  //d0sig
  h->h1D(prefix+"d0sig", "", suffix)->Fill(d0sig, weight);
  
  // Cos(metPhi, lepPhi)
  h->h1D(prefix+"CosDPhi_effBins","", suffix)->Fill(cos(evt.met().DeltaPhi(l)), weight);
  h->h1D(prefix+"CosDPhi", 	   "", suffix)->Fill(cos(evt.met().DeltaPhi(l)), weight);
  
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
  h->h1D(prefix+"closejl_minDeltaR", 		"", suffix)->Fill(closejl_deltaR, weight);
  h->h1D(prefix+"closejl_minDeltaR_effBins", 	"", suffix)->Fill(closejl_deltaR, weight);
  h->h1D(prefix+"closejl_pt", 			"", suffix)->Fill(closejl_pt*1e-3, weight);
  h->h1D(prefix+"closejl_jvt", 		"", suffix)->Fill(evt.jet()[closejl_idx].jvt(), weight);

  h->h1D(prefix+"trueMtt", "", suffix)->Fill(evt.MC_ttbar_beforeFSR().M()*1e-3, weight);
  _tree_truemtt = evt.MC_ttbar_beforeFSR().M()*1e-3;
  if (m_boosted) {
    
    size_t close_idx = 0;
    for (; close_idx < evt.jet().size(); ++close_idx){
      if (evt.jet()[close_idx].closeToLepton())
        break;
    }
    const TLorentzVector &sj = evt.jet()[close_idx].mom();
    h->h1D(prefix+"closeJetPt", "", suffix)->Fill(sj.Perp()*1e-3, weight);
    
    size_t goodljet_idx = 0;
    for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx) {
      //std::cout << "jet " << goodljet_idx << " pt = " << evt.largeJet()[goodljet_idx].mom().Perp() << ", good = " << evt.largeJet()[goodljet_idx].good() << std::endl;
      if (evt.largeJet()[goodljet_idx].good())
        break;
    }
    
    //std::cout << "idx = " << goodljet_idx << "/" << evt.largeJet().size() << std::endl;
    const TLorentzVector &lj = evt.largeJet()[goodljet_idx].mom();
    h->h1D(prefix+"largeJetPt", "", suffix)->Fill(lj.Perp()*1e-3, weight);
    h->h1D(prefix+"largeJetM", "", suffix)->Fill(lj.M()*1e-3, weight);
    h->h1D(prefix+"largeJetEta", "", suffix)->Fill(lj.Eta(), weight);
    h->h1D(prefix+"largeJetPhi", "", suffix)->Fill(lj.Phi(), weight);
    h->h1D(prefix+"largeJetSd12", "", suffix)->Fill(evt.largeJet()[goodljet_idx].split12()*1e-3, weight);
    if(evt.largeJet().size() > 0) {  
    h->h1D(prefix+"largeJet_tau32", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau32"), weight);
    h->h1D(prefix+"largeJet_tau32_wta", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau32_wta"), weight);

    h->h1D("largeJet_tau21", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau21"), weight);
    h->h1D("largeJet_tau21_wta", "", suffix)->Fill(evt.largeJet()[goodljet_idx].subs("tau21_wta"), weight);
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
    h->h1D(prefix+"mtt", "", suffix)->Fill(mtt*1e-3, weight);
    h->h1D(prefix+"mtlep_boo", "", suffix)->Fill((sj+nu+l).M()*1e-3, weight);

    _tree_mtt = mtt*1e-3;
    _tree_weight = weight;
    _tree_cat = -1;
    if (evt.passes("bejets") && m_boosted && m_electron) _tree_cat = 0;
    if (evt.passes("bmujets") && m_boosted && !m_electron) _tree_cat = 1;
    _tree_syst = suffix;
    h->m_tree->Fill();
    
  }else{
    
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

    std::vector<TLorentzVector *> vjets;
    std::vector<bool> vjets_btagged;
    for (size_t z = 0; z < evt.jet().size(); ++z) {
      vjets.push_back(new TLorentzVector(0,0,0,0));
      vjets[z]->SetPtEtaPhiE(evt.jet()[z].mom().Perp(), evt.jet()[z].mom().Eta(), evt.jet()[z].mom().Phi(), evt.jet()[z].mom().E());
      // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BTagingxAODEDM
      // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks
      vjets_btagged.push_back(evt.jet()[z].btag_mv2c20_60());
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
    h->h1D(prefix+"mtt", "", suffix)->Fill(mtt*1e-3, weight);
    h->h1D(prefix+"mtlep_res", "", suffix)->Fill(mtl*1e-3, weight);
    h->h1D(prefix+"mthad_res", "", suffix)->Fill(mth*1e-3, weight);
    h->h1D(prefix+"mwhad_res", "", suffix)->Fill(mwh*1e-3, weight);
    h->h1D(prefix+"chi2", "", suffix)->Fill(log10(chi2Value), weight);

    _tree_mtt = mtt*1e-3;
    _tree_weight = weight;
    _tree_syst = suffix;
    _tree_cat = -1;
    if (evt.passes("rejets") && !m_boosted && m_electron) _tree_cat = 2;
    if (evt.passes("rmujets") && !m_boosted && !m_electron) _tree_cat = 3;
    h->m_tree->Fill();

  }
 
}//GetHistograms


void AnaTtresMM::IniHistograms(std::string &suffix){

  m_hSvc.m_tree->Branch("truemtt",    &_tree_truemtt);
  m_hSvc.m_tree->Branch("mtt",    &_tree_mtt);
  m_hSvc.m_tree->Branch("weight", &_tree_weight);
  m_hSvc.m_tree->Branch("cat",    &_tree_cat);
  m_hSvc.m_tree->Branch("syst",   &_tree_syst);
  
  //double varBin1[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 340, 380, 450, 500};
  double varBin1[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 190, 500};
  int varBinN1 = sizeof(varBin1)/sizeof(double) - 1;
  double varBin2[] = {250, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 540, 580, 620, 660, 700, 800, 1e3, 1.2e3, 1.5e3};
  int varBinN2 = sizeof(varBin2)/sizeof(double) - 1;
  double varBin3[] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500};
  int varBinN3 = sizeof(varBin3)/sizeof(double) - 1;
  double varBin4[] = {80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500};
  int varBinN4 = sizeof(varBin4)/sizeof(double) - 1;

  double varBin5[] = {0, 80, 160, 240, 320, 400, 480, 560,640,720,800,920,1040,1160,1280,1400,1550,1700,2000,2300,2600,2900,3200,3600,4100,4600,5100,6000};
  int varBinN5 = sizeof(varBin5)/sizeof(double) - 1;

  //double CdPhiBin[] = {-1.0, -0.98, -0.8, -0.20, 0.20, 0.40, 0.60, 0.70, 0.80, 0.85, 0.90, 0.92, 0.94, 0.96, 0.98, 1.0};
  double CdPhiBin[] = {0.0, 0.5, 0.7, 0.9, 1.0};
  int CdPhiBinN = sizeof(CdPhiBin)/sizeof(double) - 1; 

  //Plots for QCD validation
  if (m_electron){
     if(m_boosted){
        double varBin6[]  = {30, 40, 60, 120, 700};
	double varBin7[]  = {0., 0.4, 0.6, 0.8, 1.0, 1.5}; 
	double varBin8[]  = {0.0, 1.6, 2.5};
	int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
	int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;
	int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;	
	m_hSvc.create1DVar(suffix+"lepPt_effBins", 		"; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
	m_hSvc.create1DVar(suffix+"closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);
	m_hSvc.create1DVar(suffix+"lepEta_effBins",  		"; |#eta| of lepton [GeV]; Eff", varBinN8, varBin8);
	m_hSvc.create1DVar(suffix+"CosDPhi_effBins",  		"; Cos( #Delta #phi(met, lept) ); Events", CdPhiBinN, CdPhiBin);
	
     }else{
        double varBin6[]  = {30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 700};
	double varBin7[]  = {0.4, 1.0, 2.6, 7.0};
	double varBin8[]  = {0.0, 1.8, 2.5};
	int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
	int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;
	int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;	
	m_hSvc.create1DVar(suffix+"lepPt_effBins", 		"; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
	m_hSvc.create1DVar(suffix+"closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);
	m_hSvc.create1DVar(suffix+"lepEta_effBins",  		"; |#eta| of lepton [GeV]; Eff", varBinN8, varBin8);
	m_hSvc.create1DVar(suffix+"CosDPhi_effBins",  		"; Cos( #Delta #phi(met, lept) ); Events", CdPhiBinN, CdPhiBin);
	
     }//m_boosted

  }else{
     if(m_boosted){
        double varBin6[] = {25, 35, 40, 50, 100, 700};
	double varBin7[] = {0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5};
	double varBin8[]  = {0.0, 1.6, 2.5};
	int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
	int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;
	int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;
	m_hSvc.create1DVar(suffix+"lepPt_effBins", 		"; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
	m_hSvc.create1DVar(suffix+"closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);
	m_hSvc.create1DVar(suffix+"lepEta_effBins",  		"; |#eta| of lepton [GeV]; Eff", varBinN8, varBin8);
	m_hSvc.create1DVar(suffix+"CosDPhi_effBins",  		"; Cos( #Delta #phi(met, lept) ); Events", CdPhiBinN, CdPhiBin);

     }else{
        double varBin6[] = {25, 30, 35, 40, 50, 100, 700};
	//double varBin7[] = {0., 0.4, 0.6, 1.0, 1.5, 2.5, 7.0};
	//double varBin7[] = {0., 0.2, 0.35, 0.4, 0.5, 0.55, 0.6, 1.0, 1.5, 2.5, 7.0};
	double varBin7[] = {0., 0.4, 0.6, 7.0};
        double varBin8[]  = {0.0, 1.6, 2.5};
	int varBinN6 = sizeof(varBin6)/sizeof(double) - 1;
	int varBinN7 = sizeof(varBin7)/sizeof(double) - 1;
	int varBinN8 = sizeof(varBin8)/sizeof(double) - 1;
	m_hSvc.create1DVar(suffix+"lepPt_effBins", 		"; lepton p_{T} [GeV]; Events", varBinN6, varBin6);
	m_hSvc.create1DVar(suffix+"closejl_minDeltaR_effBins", "; min #Delta R(lep, jet); Events", varBinN7, varBin7);
	m_hSvc.create1DVar(suffix+"lepEta_effBins",  		"; |#eta| of lepton [GeV]; Eff", varBinN8, varBin8);
	m_hSvc.create1DVar(suffix+"CosDPhi_effBins",  		"; Cos( #Delta #phi(met, lept) ); Events", CdPhiBinN, CdPhiBin);
	
     }//m_boosted
  }//m_electron

  m_hSvc.create1D(suffix+"yields", "; One ; Events", 1, 0.5, 1.5);

  m_hSvc.create1DVar(suffix+"lepPt", 		"; lepton p_{T} [GeV]; Events", varBinN1, varBin1);  
  //m_hSvc.create1D("lepPt",    	"; Pt of lept [GeV]; Events", 100, 25, 525);
  m_hSvc.create1D(suffix+"lepEta", 		"; lepton #eta ; Events", 20, -2.5, 2.5);
  m_hSvc.create1D(suffix+"lepPhi", 		"; lepton #phi [rd] ; Events", 32, -3.2, 3.2);
  m_hSvc.create1DVar(suffix+"leadJetPt", 	"; leading Jet p_{T} [GeV]; Events", varBinN1, varBin1); 
  m_hSvc.create1D(suffix+"nJets", 		"; number of jets ; Events", 10, -0.5, 9.5);
  
  m_hSvc.create1DVar(suffix+"leadbJetPt", 	"; leading b-jet p_{T} [GeV]; Events", varBinN1, varBin1);
  m_hSvc.create1DVar(suffix+"leadTrkbJetPt", 	"; leading b-jet p_{T} [GeV]; Events", varBinN1, varBin1);
  m_hSvc.create1D(suffix+"nBtagJets", 		"; number of b-tagged jets ; Events", 10, 0, 10);  
  m_hSvc.create1D(suffix+"nTrkBtagJets", 	"; number of b-tagged track jets ; Events", 10, 0, 10);  
  
  m_hSvc.create1DVar(suffix+"MET", 	"; missing E_{T} [GeV]; Events", varBinN1, varBin1);
  m_hSvc.create1D(suffix+"MET_phi", 	"; missing E_{T} #phi [rd] ; Events", 32, -3.2, 3.2);

  m_hSvc.create1D(suffix+"ptvarcone",         		"; ptvarcone [GeV]; Events", 50, -10, 40); 
  m_hSvc.create1D(suffix+"topoetcone",         		"; topoetcone [GeV]; Events", 40, -10, 30); 

  //double metBin[] = {20, 30, 40, 50, 60, 100, 500};
  double metBin[] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68};
  int metBinN = sizeof(metBin)/sizeof(double) - 1; 
  m_hSvc.create1DVar(suffix+"MET_effBins",               "; Missing E_{T} [GeV]; Events", metBinN, metBin);

  double mwtBin[] = {0, 15, 40, 60, 100, 600};
  int mwtBinN = sizeof(mwtBin)/sizeof(double) - 1;
  m_hSvc.create1DVar(suffix+"mwt_effBins",               "; W transverse mass [GeV]; Events", mwtBinN, mwtBin);

  double mwtmetBin[] = {0, 10, 20, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 150, 400};
  int mwtmetBinN = sizeof(mwtmetBin)/sizeof(double) - 1; 
  m_hSvc.create1DVar(suffix+"mwt_met", 	      "; MET + MWT [GeV]; Events", mwtmetBinN, mwtmetBin);
  
  double d0SigBin[] = {-100, -50, -30, -20, -10, -8, -6, -5, -4, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.7, -0.3, 0.0, 0.3, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4, 5, 6, 8, 10, 20, 30, 50, 100};
  int d0SigBinN = sizeof(d0SigBin)/sizeof(double) - 1;  
  m_hSvc.create1DVar(suffix+"d0sig", 		"; d0sig; Events", d0SigBinN, d0SigBin);

  m_hSvc.create1D(suffix+"CosDPhi",  		"; Cos( #Delta #phi(met, lept) ); Events", 50, -1, 1);  
  m_hSvc.create1D(suffix+"mwt", 		"; W transverse mass [GeV]; Events", 30, 0, 300);
  m_hSvc.create1D(suffix+"mu", 			"; <mu>; Events", 100, 0, 100);
  m_hSvc.create1D(suffix+"mu_original", 	"; <mu_origianl>; Events", 100, 0, 100);
  m_hSvc.create1D(suffix+"vtxz", 		";Z position of truth primary vertex; Events", 40, -400, 400);
  m_hSvc.create1D(suffix+"npv", 		"; npv; Events", 50, 0, 50);
  
  m_hSvc.create1D(suffix+"closejl_minDeltaR", 	"; min #Delta R(lep, jet); Events", 50, 0, 5);
  m_hSvc.create1DVar(suffix+"closejl_pt", 	"; Pt of closest jet to lep [GeV]; Events", varBinN1, varBin1);
  m_hSvc.create1D(suffix+"closejl_jvt", 	"; JVT of closest Jet to lepton ; Events", 20, -1, 1);  

  m_hSvc.create1D(suffix+"weight_leptSF", "; QCD weights; Events", 200, 0, 2);
  m_hSvc.create1D(suffix+"weight", "; QCD weights; Events", 2000, -100, 100);
  m_hSvc.create2D(suffix+"weight_leptPt", ";lept Pt (GeV); QCD weights",100, 25, 525, 2000, -100, 100);
  m_hSvc.create2D(suffix+"weight_mwt", ";MWT (GeV); QCD weights",30, 25, 625, 200, -10, 10); 
  m_hSvc.create2D(suffix+"weight_met", ";MET (GeV); QCD weights",30, 25, 625, 200, -10, 10);

  m_hSvc.create1D(suffix+"jet0_m", "; mass of the leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
  m_hSvc.create1D(suffix+"jet1_m", "; mass of the sub-leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
  m_hSvc.create1D(suffix+"jet2_m", "; mass of the third leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
  m_hSvc.create1D(suffix+"jet3_m", "; mass of the fourth leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
  m_hSvc.create1D(suffix+"jet4_m", "; mass of the fifth leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
  m_hSvc.create1D(suffix+"jet5_m", "; mass of the sixth leading R=0.4 calo jet [GeV]; Events", 20, 0, 100); 
  
  m_hSvc.create1DVar(suffix+"jet0_pt", "; p_{T} of the leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3);
  m_hSvc.create1DVar(suffix+"jet1_pt", "; p_{T} of the sub-leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
  m_hSvc.create1DVar(suffix+"jet2_pt", "; p_{T} of the third leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
  m_hSvc.create1DVar(suffix+"jet3_pt", "; p_{T} of the fourth leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
  m_hSvc.create1DVar(suffix+"jet4_pt", "; p_{T} of the fifth leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
  m_hSvc.create1DVar(suffix+"jet5_pt", "; p_{T} of the sixth leading R=0.4 calo jet [GeV]; Events", varBinN3, varBin3); 
  
  m_hSvc.create1D(suffix+"jet0_eta", "; #eta of the leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
  m_hSvc.create1D(suffix+"jet1_eta", "; #eta of the sub-leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
  m_hSvc.create1D(suffix+"jet2_eta", "; #eta of the third leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
  m_hSvc.create1D(suffix+"jet3_eta", "; #eta of the fourth leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
  m_hSvc.create1D(suffix+"jet4_eta", "; #eta of the fifth leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
  m_hSvc.create1D(suffix+"jet5_eta", "; #eta of the sixth leading R=0.4 calo jet ; Events", 25, -2.5, 2.5); 
  
  m_hSvc.create1D(suffix+"jet0_phi", "; #phi of the leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
  m_hSvc.create1D(suffix+"jet1_phi", "; #phi of the sub-leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
  m_hSvc.create1D(suffix+"jet2_phi", "; #phi of the third leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
  m_hSvc.create1D(suffix+"jet3_phi", "; #phi of the fourth leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
  m_hSvc.create1D(suffix+"jet4_phi", "; #phi of the fifth leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 
  m_hSvc.create1D(suffix+"jet5_phi", "; #phi of the sixth leading R=0.4 calo jet ; Events", 32, -3.2, 3.2); 

  m_hSvc.create1DVar(suffix+"tjet0_pt", "; p_{T} of the leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3);
  m_hSvc.create1DVar(suffix+"tjet1_pt", "; p_{T} of the sub-leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
  m_hSvc.create1DVar(suffix+"tjet2_pt", "; p_{T} of the third leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
  m_hSvc.create1DVar(suffix+"tjet3_pt", "; p_{T} of the fourth leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
  m_hSvc.create1DVar(suffix+"tjet4_pt", "; p_{T} of the fifth leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
  m_hSvc.create1DVar(suffix+"tjet5_pt", "; p_{T} of the sixth leading R=0.2 track jet [GeV]; Events", varBinN3, varBin3); 
  
  if (m_boosted) {
    m_hSvc.create1DVar(suffix+"closeJetPt", "; selected Jet p_{T} [GeV] ; Events", varBinN1, varBin1);
    m_hSvc.create1DVar(suffix+"largeJetPt", "; large jet p_{T} [GeV] ; Events", varBinN2, varBin2);
    m_hSvc.create1D(suffix+"largeJetM", "; large jet mass [GeV] ; Events", 30, 0, 300);
    m_hSvc.create1D(suffix+"largeJetEta", "; large jet #eta ; Events", 20, -2., 2.);
    m_hSvc.create1D(suffix+"largeJetPhi", "; large jet #phi [rd] ; Events", 32, -3.2, 3.2);
    m_hSvc.create1D(suffix+"largeJetSd12", "; large jet #sqrt{d_{12}} [GeV] ; Events", 20, 0, 200);
    m_hSvc.create1DVar(suffix+"mtlep_boo", "; leptonic top mass [GeV] ; Events", varBinN4, varBin4);

  } else {
    m_hSvc.create1DVar(suffix+"mtlep_res", "; leptonic top mass [GeV]; Events", varBinN4, varBin4);
    m_hSvc.create1DVar(suffix+"mthad_res", "; hadronic top mass [GeV]; Events", varBinN4, varBin4);
    m_hSvc.create1DVar(suffix+"mwhad_res", "; hadronic W boson mass [GeV]; Events", 40, 0, 400);
    m_hSvc.create1D(suffix+"chi2", "; log(#chi^{2}) ; Events", 50, -3, 7);
  }

  m_hSvc.create1DVar(suffix+"mtt", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
  m_hSvc.create1DVar(suffix+"trueMtt", "; mass of the top-antitop system [GeV]; Events", varBinN5, varBin5);
  
  m_hSvc.create1D(suffix+"largeJet_tau32", "; #tau_{32}; Events", 20, 0, 1);
  m_hSvc.create1D(suffix+"largeJet_tau32_wta", "; #tau_{32} wta; Events", 20, 0, 1);

  m_hSvc.create1D(suffix+"largeJet_tau21", "; #tau_{21}; Events", 20, 0, 1);
  m_hSvc.create1D(suffix+"largeJet_tau21_wta", "; #tau_{21} wta; Events", 20, 0, 1);
  

}//IniHistograms

