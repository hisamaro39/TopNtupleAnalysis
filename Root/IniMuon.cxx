/**
 * @brief IniMuon representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MObject.h"
#include "TopNtupleAnalysis/IniMuon.h"
#include <cmath>

IniMuon::IniMuon()
  : MObject() {
  m_type = MObject::inimu;
}

IniMuon::IniMuon(const TLorentzVector &v)
  : MObject(v, MObject::inimu) {
}

IniMuon::~IniMuon() {
}

const float IniMuon::pt() const {
  return m_pt;
}

float &IniMuon::pt() {
  return m_pt;
}

const float IniMuon::eta() const {
  return m_eta;
}

float &IniMuon::eta() {
  return m_eta;
}

const float IniMuon::phi() const {
  return m_phi;
}

float &IniMuon::phi() {
  return m_phi;
}

const float IniMuon::charge() const {
  return m_charge;
}

float &IniMuon::charge() {
  return m_charge;
}

const float IniMuon::ptvarcone30() const {
  return m_ptvarcone30;
}

float &IniMuon::ptvarcone30() {
  return m_ptvarcone30;
}

const float IniMuon::d0sig() const {
  return m_d0sig;
}

float &IniMuon::d0sig() {
  return m_d0sig;
}

const float IniMuon::z0sintheta() const {
  return m_z0sintheta;
}

float &IniMuon::z0sintheta() {
  return m_z0sintheta;
}

const int IniMuon::quality() const {
  return m_quality;
}

int &IniMuon::quality() {
  return m_quality;
}

const int IniMuon::accept() const {
  return m_accept;
}

int &IniMuon::accept() {
  return m_accept;
}

const int IniMuon::trigger_match() const {
  return m_trigger_match;
}

int &IniMuon::trigger_match() {
  return m_trigger_match;
}
