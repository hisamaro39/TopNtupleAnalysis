/**
 * @brief Event class, with information read off the input file read using MiniTree.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef EVENT_H
#define EVENT_H

#include <vector>
#include "TopNtupleAnalysis/Electron.h"
#include "TopNtupleAnalysis/Muon.h"
#include "TopNtupleAnalysis/Jet.h"
#include "TopNtupleAnalysis/LargeJet.h"

class Event {
  public:
    Event();
    virtual ~Event();

    void clear();

    std::vector<Electron> &electron();
    std::vector<Muon> &muon();
    std::vector<Jet> &jet();
    std::vector<LargeJet> &largeJet();

    std::vector<Jet> &tjet();
    std::vector<Jet> &akt4truthjet();
    
    const bool trigger(const std::string &t) const;
    void setTrigger(const std::string &t, bool passes);

    const std::vector<Electron> &electron() const;
    const std::vector<Muon> &muon() const;
    const std::vector<Jet> &jet() const;
    const std::vector<LargeJet> &largeJet() const;

    const std::vector<Jet> &tjet() const;
    const std::vector<Jet> &akt4truthjet() const;

    void met(float met_x, float met_y);
    TLorentzVector met() const;

    unsigned int &runNumber();
    const unsigned int runNumber() const;
    unsigned int &randomRunNumber();
    const unsigned int randomRunNumber() const;
    const unsigned int runNumber_or_RandomRunNumber() const;

    ULong64_t &eventNumber();
    const ULong64_t eventNumber() const;

    bool &isData();
    const bool isData() const;

    int &channelNumber();
    const int channelNumber() const;

    float &weight_pileup();
    const float weight_pileup() const;

    float &mu();
    const float mu() const;
    
    float &mu_original();
    const float mu_original() const;

    float &weight_mc();
    const float weight_mc() const;
    
    float &weight_bTagSF();
    const float weight_bTagSF() const;
    
    float &weight_leptonSF();
    const float weight_leptonSF() const;
    
    float &weight_Sherpa22_corr();
    const float weight_Sherpa22_corr() const;

    int &Wfilter_Sherpa_nT();
    const int Wfilter_Sherpa_nT() const;

    std::vector<std::string> &passes();
    const bool passes(const std::string &selection) const;

    // not implemented yet:
    unsigned int &lbn();
    const unsigned int lbn() const;
 
    const int hfor() const;
    int &hfor();

    int &npv();
    const int npv() const;

    float &vtxz();
    const float vtxz() const;
    
    //MC
    TLorentzVector &MC_Wdecay1_from_t();
    const TLorentzVector &MC_Wdecay1_from_t() const;
    int &MC_Wdecay1_from_t_pdgId();
    const int MC_Wdecay1_from_t_pdgId() const;
    
    TLorentzVector &MC_Wdecay2_from_t();
    const TLorentzVector &MC_Wdecay2_from_t() const;
    int &MC_Wdecay2_from_t_pdgId();
    const int MC_Wdecay2_from_t_pdgId() const;
    
    TLorentzVector &MC_b_from_t();
    const TLorentzVector &MC_b_from_t() const;
    
    TLorentzVector &MC_Wdecay1_from_tbar();
    const TLorentzVector &MC_Wdecay1_from_tbar() const;
    int &MC_Wdecay1_from_tbar_pdgId();
    const int MC_Wdecay1_from_tbar_pdgId() const;
    
    TLorentzVector &MC_Wdecay2_from_tbar();
    const TLorentzVector &MC_Wdecay2_from_tbar() const;
    int &MC_Wdecay2_from_tbar_pdgId();
    const int MC_Wdecay2_from_tbar_pdgId() const;
    
    TLorentzVector &MC_b_from_tbar();
    const TLorentzVector &MC_b_from_tbar() const;
    
    //MC
    TLorentzVector &MC_w1h();
    const TLorentzVector &MC_w1h() const;
    int &MC_w1h_pdgId();
    const int MC_w1h_pdgId() const;
    
    TLorentzVector &MC_w2h();
    const TLorentzVector &MC_w2h() const;
    int &MC_w2h_pdgId();
    const int MC_w2h_pdgId() const;
    
    TLorentzVector &MC_bh();
    const TLorentzVector &MC_bh() const;
    
    TLorentzVector &MC_w1l();
    const TLorentzVector &MC_w1l() const;
    int &MC_w1l_pdgId();
    const int MC_w1l_pdgId() const;
    
    TLorentzVector &MC_w2l();
    const TLorentzVector &MC_w2l() const;
    int &MC_w2l_pdgId();
    const int MC_w2l_pdgId() const;
    
    TLorentzVector &MC_bl();
    const TLorentzVector &MC_bl() const;

    TLorentzVector &MC_ttbar_beforeFSR();
    const TLorentzVector &MC_ttbar_beforeFSR() const;
    
    int &MC_ttbar_type();
    const int MC_ttbar_type() const;
    
    //MA
    TLorentzVector &MA_w1h();
    const TLorentzVector &MA_w1h() const;
    int &MA_w1h_pdgId();
    const int MA_w1h_pdgId() const;
    
    TLorentzVector &MA_w2h();
    const TLorentzVector &MA_w2h() const;
    int &MA_w2h_pdgId();
    const int MA_w2h_pdgId() const;
    
    TLorentzVector &MA_bh();
    const TLorentzVector &MA_bh() const;
    
    TLorentzVector &MA_w1l();
    const TLorentzVector &MA_w1l() const;
    
    int &MA_w1l_pdgId();
    const int MA_w1l_pdgId() const;
    
    TLorentzVector &MA_w2l();
    const TLorentzVector &MA_w2l() const;
    
    int &MA_w2l_pdgId();
    const int MA_w2l_pdgId() const;
        
    TLorentzVector &MA_bl();
    const TLorentzVector &MA_bl() const;

    TLorentzVector &MC_t();
    const TLorentzVector &MC_t() const;

    TLorentzVector &MC_tbar();
    const TLorentzVector &MC_tbar() const;

  protected:

    std::vector<std::string> m_trigger;
    std::vector<std::string> m_passes;

    int m_hfor;

    std::vector<Electron> m_electron;
    std::vector<Muon> m_muon;
    std::vector<Jet> m_jet;
    std::vector<LargeJet> m_largeJet;

    TLorentzVector m_met;

    std::vector<Jet> m_tjet;
    
    std::vector<Jet> m_akt4truthjet;
    
    TLorentzVector m_MC_w1h;
    int m_MC_w1h_pdgId;
    TLorentzVector m_MC_w2h;
    int m_MC_w2h_pdgId;   
    TLorentzVector m_MC_bh;
    TLorentzVector m_MC_w1l;
    int m_MC_w1l_pdgId;
    TLorentzVector m_MC_w2l;
    int m_MC_w2l_pdgId;   
    TLorentzVector m_MC_bl;    
    
    TLorentzVector m_MC_ttbar_beforeFSR;

    TLorentzVector m_MC_t;
    TLorentzVector m_MC_tbar;
    int m_MC_ttbar_type;
    
    TLorentzVector m_MC_Wdecay1_from_t;
    int m_MC_Wdecay1_from_t_pdgId;
    TLorentzVector m_MC_Wdecay2_from_t;
    int m_MC_Wdecay2_from_t_pdgId;   
    TLorentzVector m_MC_b_from_t;
    TLorentzVector m_MC_Wdecay1_from_tbar;
    int m_MC_Wdecay1_from_tbar_pdgId;
    TLorentzVector m_MC_Wdecay2_from_tbar;
    int m_MC_Wdecay2_from_tbar_pdgId;   
    TLorentzVector m_MC_b_from_tbar;
    
    TLorentzVector m_MA_w1h;
    int m_MA_w1h_pdgId;
    TLorentzVector m_MA_w2h;
    int m_MA_w2h_pdgId;   
    TLorentzVector m_MA_bh;
    TLorentzVector m_MA_w1l;
    int m_MA_w1l_pdgId;
    TLorentzVector m_MA_w2l;   
    int m_MA_w2l_pdgId;   
    TLorentzVector m_MA_bl; 

    unsigned int m_runNumber;
    unsigned int m_randomRunNumber;
    ULong64_t m_eventNumber;
    int m_channelNumber;
    bool m_isData;
    float m_mu;
    float m_mu_original;

    int m_npv;
    float m_vtxz;
    float m_rho;

    float m_weight_mc;
    float m_weight_pileup;
    float m_weight_bTagSF;
    float m_weight_leptonSF;
    float m_weight_Sherpa22_corr;
    int m_Wfilter_Sherpa_nT;
    
    unsigned int m_lbn;
};

#endif

