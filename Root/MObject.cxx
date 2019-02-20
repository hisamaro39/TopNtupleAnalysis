/**
 * @brief Base class for all Analysis objects (electrons, muons, jets, large-R jets, etc).
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MObject.h"

MObject::MObject()
  : m_mom(0,0,0,0), m_type(MObject::jet) {
}

MObject::MObject(const TLorentzVector &v, const Type t)
  :m_mom(v), m_type(t) {
}

MObject::~MObject() {
}

MObject::Type &MObject::type() {
  return m_type;
}

const MObject::Type &MObject::type() const {
  return m_type;
}

const TLorentzVector &MObject::mom() const {
  return m_mom;
}

TLorentzVector &MObject::mom() {
  return m_mom;
}

// must be in the header file
/*
template <class C>
float MObject::minDeltaR(const std::vector<C> &o) const {
  float dR = 99;
  for (int k = 0; k < o.size(); ++k) {
    float tdR = mom().DeltaR(o[k].mom());
    if (tdR < dR) {
      dR = tdR;
    }
  }
  return dR;
}
*/

