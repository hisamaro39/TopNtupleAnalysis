/**
 * @brief Electron representation for information read off input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MObject.h"
#include "TopNtupleAnalysis/Electron.h"
#include <cmath>
#include "TopNtupleAnalysis/Jet.h"
#include <iostream>

Electron::Electron()
  : MObject(), m_mi(-1), m_isTightPP(false), m_mom_calo(0,0,0,0), m_mom_trk(0,0,0,0) {
  m_type = MObject::e;
}

Electron::Electron(const TLorentzVector &v)
  : MObject(v, MObject::e), m_mi(-1), m_isTightPP(false), m_mom_calo(v), m_mom_trk(v) {
}

Electron::~Electron() {
}

void Electron::setMI(float iso) {
  m_mi = iso;
}

float Electron::mi() const {
  return m_mi;
}

void Electron::setTightPP(char t) {  
  if (t)	m_isTightPP = true;
  else		m_isTightPP = false;
}

bool Electron::isTightPP() const {
  return m_isTightPP;
}

void Electron::setMediumPP(bool t) {
  m_isMediumPP = t;
}

bool Electron::isMediumPP() const {
  return m_isMediumPP;
}

const TLorentzVector &Electron::caloMom() const {
  return m_mom_calo;
}

TLorentzVector &Electron::caloMom() {
  return m_mom_calo;
}

const TLorentzVector &Electron::trkMom() const {
  return m_mom_trk;
}

TLorentzVector &Electron::trkMom() {
  return m_mom_trk;
}

const float Electron::Dz0() const {
  return m_Dz0;
}

float &Electron::Dz0() {
  return m_Dz0;
}

const float Electron::d0() const {
  return m_d0;
}

float &Electron::d0() {
  return m_d0;
}

const float Electron::sd0() const {
  return m_sd0;
}
float &Electron::sd0() {
  return m_sd0;
}

const int Electron::author() const {
  return m_author;
}

int &Electron::author() {
  return m_author;
}

const int Electron::nSiHits() const {
  return m_nSiHits;
}

int &Electron::nSiHits() {
  return m_nSiHits;
}

const int Electron::oq() const {
  return m_oq;
}

int &Electron::oq() {
  return m_oq;
}

const int Electron::isEM() const {
  return m_isEM;
}

int &Electron::isEM() {
  return m_isEM;
}

const float Electron::ptvarcone20() const {
  return m_ptvarcone20;
}
float &Electron::ptvarcone20() {
  return m_ptvarcone20;
}

const float Electron::topoetcone20() const {
  return m_topoetcone20;
}
float &Electron::topoetcone20() {
  return m_topoetcone20;
}

//Triggers

const bool Electron::HLT_e24_lhmedium_iloose_L1EM20VH() const {
  return m_HLT_e24_lhmedium_iloose_L1EM20VH;
}

bool &Electron::HLT_e24_lhmedium_iloose_L1EM20VH(){
  return m_HLT_e24_lhmedium_iloose_L1EM20VH;
}


const bool Electron::HLT_e24_lhmedium_L1EM18VH() const {
  return m_HLT_e24_lhmedium_L1EM18VH;
}

bool &Electron::HLT_e24_lhmedium_L1EM18VH() {
  return m_HLT_e24_lhmedium_L1EM18VH;
}


const bool Electron::HLT_e24_lhmedium_L1EM20VH() const {
  return m_HLT_e24_lhmedium_L1EM20VH;
}

bool &Electron::HLT_e24_lhmedium_L1EM20VH() {
  return m_HLT_e24_lhmedium_L1EM20VH;
}


const bool Electron::HLT_e60_lhmedium() const {
  return m_HLT_e60_lhmedium;
}

bool &Electron::HLT_e60_lhmedium() {
  return m_HLT_e60_lhmedium;
}


const bool Electron::HLT_e120_lhloose() const {
  return m_HLT_e120_lhloose;
}

bool &Electron::HLT_e120_lhloose() {
  return m_HLT_e120_lhloose;
}


bool Electron::pass() const {
  if (mom().Perp() < 25e3) return false;
  if ( (author() != 1) && (author() != 3) ) return false;
  float eta = std::fabs(caloMom().Eta());
  if (eta > 1.37 && eta < 1.52) return false;
  if (eta > 2.47) return false;
  if (!isTightPP()) return false;
  if (std::fabs(Dz0()) > 0.5) return false;
  if (mi()/mom().Perp() > 0.05) return false;
  if ( (oq() & 1446) != 0 ) return false;
  return true;
}

bool Electron::passLoose() const {
  if (mom().Perp() < 25e3) return false;
  if ( (author() != 1) && (author() != 3) ) return false;
  float eta = std::fabs(caloMom().Eta());
  if (eta > 1.37 && eta < 1.52) return false;
  if (eta > 2.47) return false;
  if (std::fabs(Dz0()) > 0.5) return false;
  if ( (oq() & 1446) != 0 ) return false;

  if (!isMediumPP()) return false;
  if ((isEM() & (1<<0x1)) != 0) return false; 
  return true;
}

