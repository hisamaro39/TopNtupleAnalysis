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
  m_trueFlavour = l.m_trueFlavour;
  m_good = l.m_good;
  m_subs = l.m_subs;
}

LargeJet::~LargeJet() {
}

void LargeJet::setGood(bool b) {
  m_good = b;
}

bool LargeJet::good() const {
  return m_good;
}

bool &LargeJet::good() {
  return m_good;
}

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

float &LargeJet::subs(const std::string &s) {
  return m_subs[s];
}

const float LargeJet::subs(const std::string &s) const {
  return m_subs.at(s);
}

