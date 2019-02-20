/**
 * @brief Large-R jet representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MObject.h"
#include "TopNtupleAnalysis/LargeJet.h"
#include <cmath>
#include "TLorentzVector.h"

LargeJet::LargeJet()
  : MObject() {
  m_type = MObject::largejet;
}

LargeJet::LargeJet(const TLorentzVector &v)
  : MObject(v, MObject::largejet) {
}

LargeJet::LargeJet(const LargeJet &l)
  : MObject(l.mom(), MObject::largejet) {
  m_split12 = l.m_split12;
  m_tau32 = l.m_tau32;
  m_BDT_TOPtag = l.m_BDT_TOPtag;
  m_trueFlavour = l.m_trueFlavour;
  m_good = l.m_good;
  m_good50 = l.m_good50;
  m_good_sub50 = l.m_good_sub50;
  m_good_sub80 = l.m_good_sub80;
  m_good_smooth_mt80 = l.m_good_smooth_mt80;
  m_good_smooth_mt50 = l.m_good_smooth_mt50;
  m_good_smooth_ts80 = l.m_good_smooth_ts80;
  m_good_smooth_ts50 = l.m_good_smooth_ts50;
  m_good_smooth_qt80 = l.m_good_smooth_qt80;
  m_good_smooth_qt50 = l.m_good_smooth_qt50;
  m_good_bdt80 = l.m_good_bdt80;
  m_good_dnn80 = l.m_good_dnn80;
  m_subs = l.m_subs;
}

LargeJet::~LargeJet() {
}

void LargeJet::setGood(bool b) {
  m_good = b;
}

void LargeJet::setGood50(bool b) {
  m_good50 = b;
}

bool LargeJet::good() const {
  return m_good;
}

bool &LargeJet::good() {
  return m_good;
}

bool LargeJet::good50() const {
  return m_good50;
}

bool &LargeJet::good50() {
  return m_good50;
}

bool LargeJet::good_sub80() const {return m_good_sub80;}

bool &LargeJet::good_sub80() {return m_good_sub80;}

bool LargeJet::good_sub50() const {return m_good_sub50;}

bool &LargeJet::good_sub50() {return m_good_sub50;}

bool LargeJet::good_smooth_mt80() const {return m_good_smooth_mt80;}

bool &LargeJet::good_smooth_mt80() {return m_good_smooth_mt80;}

bool LargeJet::good_smooth_mt50() const {return m_good_smooth_mt50;}

bool &LargeJet::good_smooth_mt50() {return m_good_smooth_mt50;}

bool LargeJet::good_smooth_ts80() const {return m_good_smooth_ts80;}

bool &LargeJet::good_smooth_ts80() {return m_good_smooth_ts80;}

bool LargeJet::good_smooth_ts50() const {return m_good_smooth_ts50;}

bool &LargeJet::good_smooth_ts50() {return m_good_smooth_ts50;}

bool LargeJet::good_smooth_qt80() const {return m_good_smooth_qt80;}

bool &LargeJet::good_smooth_qt80() {return m_good_smooth_qt80;}

bool LargeJet::good_smooth_qt50() const {return m_good_smooth_qt50;}

bool &LargeJet::good_smooth_qt50() {return m_good_smooth_qt50;}

bool LargeJet::good_bdt80() const {return m_good_bdt80;}

bool &LargeJet::good_bdt80() {return m_good_bdt80;}

bool LargeJet::good_dnn80() const {return m_good_dnn80;}

bool &LargeJet::good_dnn80() {return m_good_dnn80;}

const int LargeJet::trueFlavour() const {
  return m_trueFlavour;
}

int &LargeJet::trueFlavour() {
  return m_trueFlavour;
}

bool LargeJet::pass() const {
  if (std::fabs(mom().Eta()) > 1.2) return false;
  if (mom().Perp() < 200e3) return false;
  return true;
}

bool LargeJet::passLoose() const {
  if (std::fabs(mom().Eta()) > 2.0) return false;
  if (mom().Perp() < 200e3) return false;
  return true;
}

bool LargeJet::passFakeSelection(const TLorentzVector &lept, const TLorentzVector &selJet) const {
  if (std::fabs(mom().Eta()) > 2.0) 			return false;
  if (mom().Perp() < 200e3) 				return false;
  if (mom().M() > 70e3) 				return false;  
  if(std::fabs(mom().DeltaPhi(lept)) < 2.3)		return false;  
  if(mom().DeltaR(selJet) < 1.5)			return false;
        
  return true;
}

double &LargeJet::split12() {
  return m_split12;
}
const double LargeJet::split12() const {
  return m_split12;
}

double &LargeJet::tau32() {
  return m_tau32;
}
const double LargeJet::tau32() const {
  return m_tau32;
}

double &LargeJet::BDT_TOPtag() {
  return m_BDT_TOPtag;
}
const double LargeJet::BDT_TOPtag() const {
  return m_BDT_TOPtag;
}

float &LargeJet::subs(const std::string &s) {
  return m_subs[s];
}

const float LargeJet::subs(const std::string &s) const {
  return m_subs.at(s);
}

