/**
 * @brief Muon representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MObject.h"
#include "TopNtupleAnalysis/Muon.h"
#include <cmath>

Muon::Muon()
  : MObject(), m_mi(-1), m_tight(false), m_Dz0(0), m_z0_exPV(0), m_d0(0), m_sd0(1), m_author(0), m_passTrkCuts(false) {
  m_type = MObject::mu;
}

Muon::Muon(const TLorentzVector &v)
  : MObject(v, MObject::mu), m_mi(-1), m_tight(false), m_Dz0(0), m_z0_exPV(0), m_d0(0), m_sd0(1), m_author(0), m_passTrkCuts(false) {
}

Muon::~Muon() {
}

void Muon::setMI(float iso) {
  m_mi = iso;
}

float Muon::mi() const {
  return m_mi;
}

void Muon::setTight(char t) {
  if (t)	m_tight = true;
  else		m_tight = false;
}

bool Muon::isTight() const {
  return m_tight;
}

const float Muon::Dz0() const {
  return m_Dz0;
}

float &Muon::Dz0() {
  return m_Dz0;
}

const float Muon::z0_exPV() const {
  return m_z0_exPV;
}

float &Muon::z0_exPV() {
  return m_z0_exPV;
}

const float Muon::d0() const {
  return m_d0;
}

float &Muon::d0() {
  return m_d0;
}

const float Muon::sd0() const {
  return m_sd0;
}
float &Muon::sd0() {
  return m_sd0;
}

const float Muon::ptvarcone30() const {
  return m_ptvarcone30;
}
float &Muon::ptvarcone30() {
  return m_ptvarcone30;
}

const float Muon::topoetcone20() const {
  return m_topoetcone20;
}
float &Muon::topoetcone20() {
  return m_topoetcone20;
}

//Triggers

const bool Muon::HLT_mu50() const {
  return m_HLT_mu50;
}
bool &Muon::HLT_mu50() {
  return m_HLT_mu50;
}

const bool Muon::HLT_mu20_iloose_L1MU15() const {
  return m_HLT_mu20_iloose_L1MU15;
}
bool &Muon::HLT_mu20_iloose_L1MU15() {
  return m_HLT_mu20_iloose_L1MU15;
}

const bool Muon::HLT_mu20_L1MU15() const {
  return m_HLT_mu20_L1MU15;
}
bool &Muon::HLT_mu20_L1MU15() {
  return m_HLT_mu20_L1MU15;
}


const int Muon::author() const {
  return m_author;
}
int &Muon::author() {
  return m_author;
}

const bool Muon::passTrkCuts() const {
  return m_passTrkCuts;
}

bool &Muon::passTrkCuts() {
  return m_passTrkCuts;
}


bool Muon::pass() const {
  if (mom().Perp() < 25e3) return false;
  float eta = std::fabs(mom().Eta());
  if (eta > 2.5) return false;
  if (!isTight()) return false;

  if (author() != 12) return false;

  if (!passTrkCuts()) return false;

  if (mi()/mom().Perp() > 0.05) return false;

  if (std::fabs(d0()/sd0()) > 3) return false;
  if (std::fabs(z0_exPV()) > 2) return false;

  return true;
}

bool Muon::passLoose() const {
  if (mom().Perp() < 25e3) return false;
  float eta = std::fabs(mom().Eta());
  if (eta > 2.5) return false;
  if (!isTight()) return false;
  if (author() != 12) return false;
  if (!passTrkCuts()) return false;
  if (std::fabs(d0()/sd0()) > 3) return false;
  if (std::fabs(z0_exPV()) > 2) return false;
  //if (mi()/mom().Perp() > 0.05) return false;
  return true;
}

TLorentzVector &Muon::momME() {
  return m_momME;
}

const TLorentzVector &Muon::momME() const {
  return m_momME;
}

TLorentzVector &Muon::momMECorr() {
  return m_momMECorr;
}
const TLorentzVector &Muon::momMECorr() const {
  return m_momMECorr;
}

TLorentzVector &Muon::momMS() {
  return m_momMS;
}
const TLorentzVector &Muon::momMS() const {
  return m_momMS;
}

TLorentzVector &Muon::momID() {
  return m_momID;
}
const TLorentzVector &Muon::momID() const {
  return m_momID;
}

TLorentzVector &Muon::momTrk() {
  return m_momTrk;
}
const TLorentzVector &Muon::momTrk() const {
  return m_momTrk;
}

int Muon::charge() const {
  return m_charge;
}
int &Muon::charge() {
  return m_charge;
}

void Muon::setST(bool b) {
  m_st = b;
}
bool Muon::st() const {
  return m_st;
}

void Muon::setSA(bool b) {
  m_sa = b;
}
bool Muon::sa() const {
  return m_sa;
}

void Muon::setCB(bool b) {
  m_cb = b;
}

bool Muon::cb() const {
  return m_cb;
}
