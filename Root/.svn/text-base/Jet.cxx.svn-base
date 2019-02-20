/**
 * @brief Jet representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MObject.h"
#include "TopNtupleAnalysis/Jet.h"
#include <cmath>

Jet::Jet()
  : MObject() {
  m_type = MObject::jet;
  m_trueflavour = 1;
  m_mv1 = 0;
  m_mv2c20 = 0;
  m_ip3dsv1 = 0;
  m_jvt = 1;
  m_closeToLepton = 0;
  m_numConstituents = 0;
}

Jet::Jet(const TLorentzVector &v)
  : MObject(v, MObject::jet) {
  m_trueflavour = 1;
  m_mv1 = 0;
  m_mv2c20 = 0;
  m_ip3dsv1 = 0;
  m_jvt = 1;
  m_closeToLepton = 0;
  m_numConstituents = 0;
}

Jet::~Jet() {
}

int &Jet::trueFlavour() {
  return m_trueflavour;
}
const int Jet::trueFlavour() const {
  return m_trueflavour;
}

int &Jet::numConstituents() {
  return m_numConstituents;
}
const int Jet::numConstituents() const {
  return m_numConstituents;
}

bool Jet::pass() const {
  if (std::fabs(mom().Eta()) > 2.5) return false;
  if (mom().Perp() < 25e3) return false;
  if ( (mom().Perp() <= 50e3) && (std::fabs(mom().Eta()) < 2.4) )//&& \
       (std::fabs(jvt()) <= 0.5) )
    return false;
  return true;
}


bool Jet::pass_trk() const {
  if (std::fabs(mom().Eta()) > 2.5) return false;
  if (mom().Perp() < 10e3) return false;
  if (numConstituents() < 2) return false;
  return true;
}
bool Jet::btag() const {
  if (ip3dsv1() < 1.85) return false;
  return true;
}

bool Jet::btag_mv2c20_60() const {
  if (mv2c20() < 0.5102) return false;
  return true;
}

bool Jet::btag_mv2c20_70() const {
  if (mv2c20() < -0.0436) return false;
  return true;
}

bool Jet::btag_mv2c20_70_trk() const {
  if (mv2c20() < -0.3098) return false;
  return true;
}

bool Jet::btag_mv2c10_70_trk() const {
  if (mv2c10() < 0.6455) return false; // June 23 2016
  return true;
}


const float Jet::mv1() const {
  return m_mv1;
}
float &Jet::mv1() {
  return m_mv1;
}

const float Jet::mv2c20() const {
  return m_mv2c20;
}
float &Jet::mv2c20() {
  return m_mv2c20;
}

const float Jet::mv2c10() const {
  return m_mv2c10;
}
float &Jet::mv2c10() {
  return m_mv2c10;
}

const float Jet::ip3dsv1() const {
  return m_ip3dsv1;
}
float &Jet::ip3dsv1() {
  return m_ip3dsv1;
}

const int Jet::closeToLepton() const {
  return m_closeToLepton;
}
int &Jet::closeToLepton() {
  return m_closeToLepton;
}

const float Jet::jvt() const {
  return m_jvt;
}
float &Jet::jvt() {
  return m_jvt;
}


