/**
 * @brief Truth representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MObject.h"
#include "TopNtupleAnalysis/Truth.h"
#include <cmath>

Truth::Truth()
  : MObject() {
  m_type = MObject::truth;
}

Truth::Truth(const TLorentzVector &v)
  : MObject(v, MObject::truth) {
}

Truth::~Truth() {
}

const float Truth::px() const {
  return m_px;
}

float &Truth::px() {
  return m_px;
}

const float Truth::py() const {
  return m_py;
}

float &Truth::py() {
  return m_py;
}

const float Truth::pz() const {
  return m_pz;
}

float &Truth::pz() {
  return m_pz;
}

const float Truth::pt() const {
  float pt = sqrt(m_px*m_px+m_py*m_py);
  return pt;
}

float &Truth::pt() {
  float pt = sqrt(m_px*m_px+m_py*m_py);
  return pt;
}

const float Truth::e() const {
  return m_e;
}

float &Truth::e() {
  return m_e;
}

const float Truth::m() const {
  return m_m;
}

float &Truth::m() {
  return m_m;
}

const int Truth::id() const {
  return m_id;
}

int &Truth::id() {
  return m_id;
}

const int Truth::status() const {
  return m_status;
}

int &Truth::status() {
  return m_status;
}

const int Truth::muon_type() const {
  return m_muon_type;
}

int &Truth::muon_type() {
  return m_muon_type;
}

