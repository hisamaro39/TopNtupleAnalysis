/**
 * @brief Event class, with information read off the input file read using MiniTree.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/Event.h"
#include <vector>
#include <cmath>
#include <algorithm>

Event::Event() {
}

Event::~Event() {
}

std::vector<std::string> &Event::passes() {
  return m_passes;
}

const bool Event::trigger(const std::string &t) const {
  return (std::find(m_trigger.begin(), m_trigger.end(), t) != m_trigger.end());
}

void Event::setTrigger(const std::string &t, bool passes) {
  std::vector<std::string>::iterator it = std::find(m_trigger.begin(), m_trigger.end(), t);
  if (it != m_trigger.end() && !passes) {
    m_trigger.erase(it);
  } else if (it == m_trigger.end() && passes) {
    m_trigger.push_back(t);
  }
}

const bool Event::passes(const std::string &selection) const {
  return std::find(m_passes.begin(), m_passes.end(), selection) != m_passes.end();
}

const int Event::hfor() const {
  return m_hfor;
}

int &Event::hfor() {
  return m_hfor;
}

void Event::clear() {
  m_trigger.clear();
  m_electron.clear();
  m_muon.clear();
  m_jet.clear();
  m_truthjet.clear();
  m_truthlargeJet.clear();
  m_tjet.clear();
  m_nocalibtjet.clear();
  m_vrtjet.clear();
  m_nocalibjet.clear();
  m_nocalibljet.clear();
  m_nocalibmu.clear();
  m_nocalibel.clear();
  m_largeJet.clear();
  m_akt4truthjet.clear();
  m_truth.clear();
  m_inimu.clear();
  m_met.SetPxPyPzE(0,0,0,0);
  m_passes.clear();
}

float &Event::weight_pileup() {
  return m_weight_pileup;
}
const float Event::weight_pileup() const {
  return m_weight_pileup;
}

float &Event::weight_mc() {
  return m_weight_mc;
}

const float Event::weight_mc() const {
  return m_weight_mc;
}

float &Event::MC_met() {
  return m_MC_met;
}

const float Event::MC_met() const {
  return m_MC_met;
}

float &Event::chi2_all() {
  return m_chi2_all;
}

const float Event::chi2_all() const {
  return m_chi2_all;
}

float &Event::MC_met_nomuon() {
  return m_MC_met_nomuon;
}

const float Event::MC_met_nomuon() const {
  return m_MC_met_nomuon;
}

float &Event::MC_met_muon() {
  return m_MC_met_muon;
}

const float Event::MC_met_muon() const {
  return m_MC_met_muon;
}


float &Event::weight_bTagSF() {
  return m_weight_bTagSF;
}

const float Event::weight_bTagSF() const {
  return m_weight_bTagSF;
}

float &Event::weight_leptonSF() {
  return m_weight_leptonSF;
}

const float Event::weight_leptonSF() const {
  return m_weight_leptonSF;
}

float &Event::weight_Sherpa22_corr() {
  return m_weight_Sherpa22_corr;
}

const float Event::weight_Sherpa22_corr() const {
  return m_weight_Sherpa22_corr;
}

int &Event::Wfilter_Sherpa_nT() {
  return m_Wfilter_Sherpa_nT;
}

const int Event::Wfilter_Sherpa_nT() const {
  return m_Wfilter_Sherpa_nT;
}

std::vector<Electron> &Event::electron() {
  return m_electron;
}

std::vector<Muon> &Event::muon() {
  return m_muon;
}

std::vector<Truth> &Event::truth() {
  return m_truth;
}

std::vector<IniMuon> &Event::inimu() {
  return m_inimu;
}

std::vector<Jet> &Event::jet() {
  return m_jet;
}

std::vector<Jet> &Event::truthjet() {
  return m_truthjet;
}

std::vector<LargeJet> &Event::truthlargeJet() {
  return m_truthlargeJet;
}


std::vector<Jet> &Event::tjet() {
  return m_tjet;
}

std::vector<Jet> &Event::nocalibtjet() {
  return m_nocalibtjet;
}

std::vector<Jet> &Event::vrtjet() {
  return m_vrtjet;
}

std::vector<Jet> &Event::nocalibjet() {
  return m_nocalibjet;
}

std::vector<Jet> &Event::nocalibljet() {
  return m_nocalibljet;
}

std::vector<Muon> &Event::nocalibmu() {
  return m_nocalibmu;
}

std::vector<Electron> &Event::nocalibel() {
  return m_nocalibel;
}

std::vector<LargeJet> &Event::largeJet() {
  return m_largeJet;
}

//truth jets
std::vector<Jet> &Event::akt4truthjet() {
  return m_akt4truthjet;
}

void Event::met(float met_x, float met_y) {
  m_met.SetPxPyPzE(met_x, met_y, 0, std::sqrt(std::pow(met_x, 2) + std::pow(met_y, 2)));
}

TLorentzVector Event::met() const {
  return m_met;
}

const std::vector<Electron> &Event::electron() const {
  return m_electron;
}

const std::vector<Muon> &Event::muon() const {
  return m_muon;
}

const std::vector<Truth> &Event::truth() const {
  return m_truth;
}

const std::vector<IniMuon> &Event::inimu() const {
  return m_inimu;
}

const std::vector<Jet> &Event::jet() const {
  return m_jet;
}

const std::vector<Jet> &Event::truthjet() const {
  return m_truthjet;
}

const std::vector<LargeJet> &Event::truthlargeJet() const {
  return m_truthlargeJet;
}

const std::vector<Jet> &Event::tjet() const {
  return m_tjet;
}

const std::vector<Jet> &Event::nocalibtjet() const {
  return m_nocalibtjet;
}

const std::vector<Jet> &Event::vrtjet() const {
  return m_vrtjet;
}

const std::vector<Jet> &Event::nocalibjet() const {
  return m_nocalibjet;
}

const std::vector<Jet> &Event::nocalibljet() const {
  return m_nocalibljet;
}

const std::vector<Muon> &Event::nocalibmu() const {
  return m_nocalibmu;
}

const std::vector<Electron> &Event::nocalibel() const {
  return m_nocalibel;
}

const std::vector<LargeJet> &Event::largeJet() const {
  return m_largeJet;
}

unsigned int &Event::runNumber() {
  return m_runNumber;
}

const unsigned int Event::runNumber() const {
  return m_runNumber;
}

unsigned int &Event::randomRunNumber() {
  return  m_randomRunNumber;
}

const unsigned int Event::randomRunNumber() const {
  return  m_randomRunNumber;
}

const unsigned int Event::runNumber_or_RandomRunNumber() const {
  return m_isData ? m_runNumber : m_randomRunNumber;
}

ULong64_t &Event::eventNumber() {
  return m_eventNumber;
}

const ULong64_t Event::eventNumber() const {
  return m_eventNumber;
}

bool &Event::isData() {
  return m_isData;
}
const bool Event::isData() const {
  return m_isData;
}

int &Event::channelNumber() {
  return m_channelNumber;
}
const int Event::channelNumber() const {
  return m_channelNumber;
}

float &Event::mu() {
  return m_mu;
}
const float Event::mu() const {
  return m_mu;
}

float &Event::mu_original() {
  return m_mu_original;
}
const float Event::mu_original() const {
  return m_mu_original;
}

int &Event::npv() {
  return m_npv;
}
const int Event::npv() const {
  return m_npv;
}

float &Event::vtxz(){
  return m_vtxz;
}

const float Event::vtxz() const{
  return m_vtxz;
}

unsigned int &Event::lbn() {
  return m_lbn;
}

const unsigned int Event::lbn() const {
  return m_lbn;
}

//Variables from truth

TLorentzVector &Event::MC_Wdecay1_from_t() {
  return m_MC_Wdecay1_from_t;
}
const TLorentzVector &Event::MC_Wdecay1_from_t() const {
  return m_MC_Wdecay1_from_t;
}

int &Event::MC_Wdecay1_from_t_pdgId() {
  return m_MC_Wdecay1_from_t_pdgId;
}
const int Event::MC_Wdecay1_from_t_pdgId() const {
  return m_MC_Wdecay1_from_t_pdgId;
}

TLorentzVector &Event::MC_Wdecay2_from_t() {
  return m_MC_Wdecay2_from_t;
}
const TLorentzVector &Event::MC_Wdecay2_from_t() const {
  return m_MC_Wdecay2_from_t;
}

int &Event::MC_Wdecay2_from_t_pdgId() {
  return m_MC_Wdecay2_from_t_pdgId;
}
const int Event::MC_Wdecay2_from_t_pdgId() const {
  return m_MC_Wdecay2_from_t_pdgId;
}

TLorentzVector &Event::MC_b_from_t() {
  return m_MC_b_from_t;
}
const TLorentzVector &Event::MC_b_from_t() const {
  return m_MC_b_from_t;
}

TLorentzVector &Event::MC_Wdecay1_from_tbar() {
  return m_MC_Wdecay1_from_tbar;
}
const TLorentzVector &Event::MC_Wdecay1_from_tbar() const {
  return m_MC_Wdecay1_from_tbar;
}

int &Event::MC_Wdecay1_from_tbar_pdgId() {
  return m_MC_Wdecay1_from_tbar_pdgId;
}
const int Event::MC_Wdecay1_from_tbar_pdgId() const {
  return m_MC_Wdecay1_from_tbar_pdgId;
}

TLorentzVector &Event::MC_Wdecay2_from_tbar() {
  return m_MC_Wdecay2_from_tbar;
}
const TLorentzVector &Event::MC_Wdecay2_from_tbar() const {
  return m_MC_Wdecay2_from_tbar;
}

int &Event::MC_Wdecay2_from_tbar_pdgId() {
  return m_MC_Wdecay2_from_tbar_pdgId;
}
const int Event::MC_Wdecay2_from_tbar_pdgId() const {
  return m_MC_Wdecay2_from_tbar_pdgId;
}

TLorentzVector &Event::MC_b_from_tbar() {
  return m_MC_b_from_tbar;
}
const TLorentzVector &Event::MC_b_from_tbar() const {
  return m_MC_b_from_tbar;
}

//Variables from truth: lepton+jet channel

//hadronic top decay
TLorentzVector &Event::MC_w1h() {
  return m_MC_w1h;
}
const TLorentzVector &Event::MC_w1h() const {
  return m_MC_w1h;
}

int &Event::MC_w1h_pdgId() {
  return m_MC_w1h_pdgId;
}
const int Event::MC_w1h_pdgId() const {
  return m_MC_w1h_pdgId;
}

TLorentzVector &Event::MC_w2h() {
  return m_MC_w2h;
}
const TLorentzVector &Event::MC_w2h() const {
  return m_MC_w2h;
}

int &Event::MC_w2h_pdgId() {
  return m_MC_w2h_pdgId;
}
const int Event::MC_w2h_pdgId() const {
  return m_MC_w2h_pdgId;
}

TLorentzVector &Event::MC_bh() {
  return m_MC_bh;
}
const TLorentzVector &Event::MC_bh() const {
  return m_MC_bh;
}

//leptonic top decay
TLorentzVector &Event::MC_w1l() {
  return m_MC_w1l;
}
const TLorentzVector &Event::MC_w1l() const {
  return m_MC_w1l;
}

signed int &Event::MC_w1l_pdgId() {
  return m_MC_w1l_pdgId;
}
const signed int Event::MC_w1l_pdgId() const {
  return m_MC_w1l_pdgId;
}

TLorentzVector &Event::MC_w2l() {
  return m_MC_w2l;
}
const TLorentzVector &Event::MC_w2l() const {
  return m_MC_w2l;
}

int &Event::MC_w2l_pdgId() {
  return m_MC_w2l_pdgId;
}
const int Event::MC_w2l_pdgId() const {
  return m_MC_w2l_pdgId;
}

TLorentzVector &Event::MC_bl() {
  return m_MC_bl;
}
const TLorentzVector &Event::MC_bl() const {
  return m_MC_bl;
}

TLorentzVector &Event::MC_ttbar_beforeFSR() {
  return m_MC_ttbar_beforeFSR;
}
const TLorentzVector &Event::MC_ttbar_beforeFSR() const {
  return m_MC_ttbar_beforeFSR;
}

TLorentzVector &Event::MC_t() {
  return m_MC_t;
}
const TLorentzVector &Event::MC_t() const {
  return m_MC_t;
}

TLorentzVector &Event::MC_th() {
  return m_MC_th;
}
const TLorentzVector &Event::MC_th() const {
  return m_MC_th;
}

TLorentzVector &Event::MC_tl() {
  return m_MC_tl;
}
const TLorentzVector &Event::MC_tl() const {
  return m_MC_tl;
}

TLorentzVector &Event::MC_tbar() {
  return m_MC_tbar;
}
const TLorentzVector &Event::MC_tbar() const {
  return m_MC_tbar;
}

int &Event::MC_ttbar_type() {
  return m_MC_ttbar_type;
}
const int Event::MC_ttbar_type() const {
  return m_MC_ttbar_type;
}

//MA
//hadronic top decay

TLorentzVector &Event::MA_w1h() {
  return m_MA_w1h;
}
const TLorentzVector &Event::MA_w1h() const {
  return m_MA_w1h;
}

int &Event::MA_w1h_pdgId() {
  return m_MA_w1h_pdgId;
}
const int Event::MA_w1h_pdgId() const {
  return m_MA_w1h_pdgId;
}

TLorentzVector &Event::MA_w2h() {
  return m_MA_w2h;
}
const TLorentzVector &Event::MA_w2h() const {
  return m_MA_w2h;
}

int &Event::MA_w2h_pdgId() {
  return m_MA_w2h_pdgId;
}
const int Event::MA_w2h_pdgId() const {
  return m_MA_w2h_pdgId;
}

TLorentzVector &Event::MA_bh() {
  return m_MA_bh;
}
const TLorentzVector &Event::MA_bh() const {
  return m_MA_bh;
}

//leptonic top decay
TLorentzVector &Event::MA_w1l() {
  return m_MA_w1l;
}
const TLorentzVector &Event::MA_w1l() const {
  return m_MA_w1l;
}

signed int &Event::MA_w1l_pdgId() {
  return m_MA_w1l_pdgId;
}
const signed int Event::MA_w1l_pdgId() const {
  return m_MA_w1l_pdgId;
}

TLorentzVector &Event::MA_w2l() {
  return m_MA_w2l;
}
const TLorentzVector &Event::MA_w2l() const {
  return m_MA_w2l;
}

int &Event::MA_w2l_pdgId() {
  return m_MA_w2l_pdgId;
}
const int Event::MA_w2l_pdgId() const {
  return m_MA_w2l_pdgId;
}

TLorentzVector &Event::MA_bl() {
  return m_MA_bl;
}
const TLorentzVector &Event::MA_bl() const {
  return m_MA_bl;
}

float &Event::BDT_TOPtag() {
  return m_BDT_TOPtag;
}

const float Event::BDT_TOPtag() const {
  return m_BDT_TOPtag;
}

