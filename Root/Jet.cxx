/**
 * @brief Jet representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MObject.h"
#include "TopNtupleAnalysis/Jet.h"
#include <cmath>
#include <iostream>

Jet::Jet()
  : MObject() {
  m_type = MObject::jet;
  m_trueflavour = 1;
  m_mv1 = 0;
  m_mv2c20 = 0;
  m_ip3dsv1 = 0;
  m_jvt = -1;
  m_closeToLepton = 0;
  m_numConstituents = 0;
  //m_inf = TFile::Open("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/xAODBTaggingEfficiency/13TeV/2017-21-13TeV-MC16-CDI-2017-07-02_v1.root", "read");
  //m_spline = (TSpline3*) m_inf->Get("MV2c10/AntiKt4EMTopoJets/HybBEff_70/cutprofile");
  //m_inf = TFile::Open("./share/cutprofile_MV2c10_AntiKet4EMTopoJets_HybBEff_70.root", "read");//made by myself
  //m_spline=nullptr;
  //m_spline = (TSpline3*) m_inf->Get("cutprofile");
}

Jet::Jet(const TLorentzVector &v)
  : MObject(v, MObject::jet) {
  m_trueflavour = 1;
  m_mv1 = 0;
  m_mv2c20 = 0;
  m_ip3dsv1 = 0;
  m_jvt = -1;
  m_closeToLepton = 0;
  m_numConstituents = 0;
}

Jet::~Jet() {
  //m_inf->Close();
  //delete m_spline;
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
  //if (mv2c20() < -0.0436) return false;//old
  if (mv2c20() < 0.7110) return false;
  return true;
}

bool Jet::btag_mv2c10_70() const {
  //if (mv2c20() < -0.0436) return false;//old
  if (mv2c10() < 0.831) return false;//rel21 
  return true;
}

bool Jet::btag_mv2c20_70_trk() const {
  if (mv2c20() < -0.3098) return false;
  return true;
}

bool Jet::btag_mv2c10_70_trk() const {
  //std::cout << "btag for rel20" << std::endl;
  if (mv2c10() < 0.6455) return false; // June 23 2016
  return true;
}

bool Jet::btag_mv2c10_70_trk_rel21() const {
  //std::cout << "btag for rel21" << std::endl;
  if (mv2c10() < 0.831) return false; // June 23 2016
  return true;
}

bool Jet::btag_mv2c10_70_hyb(float pt) const {
  //float cutvalue = m_spline->Eval(pt); 
  //std::cout << "cutvalue=" << cutvalue << std::endl;
  //if (mv2c10() < cutvalue) return false;
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

const float Jet::mv2c10mu() const {
  return m_mv2c10mu;
}
float &Jet::mv2c10mu() {
  return m_mv2c10mu;
}

const float Jet::mv2c10rnn() const {
  return m_mv2c10rnn;
}
float &Jet::mv2c10rnn() {
  return m_mv2c10rnn;
}


const float Jet::dl1_pu() const {
  return m_dl1_pu;
}
float &Jet::dl1_pu() {
  return m_dl1_pu;
}

const float Jet::dl1_pb() const {
  return m_dl1_pb;
}
float &Jet::dl1_pb() {
  return m_dl1_pb;
}

const float Jet::dl1_pc() const {
  return m_dl1_pc;
}
float &Jet::dl1_pc() {
  return m_dl1_pc;
}

const float Jet::dl1mu_pu() const {
  return m_dl1mu_pu;
}
float &Jet::dl1mu_pu() {
  return m_dl1mu_pu;
}

const float Jet::dl1mu_pb() const {
  return m_dl1mu_pb;
}
float &Jet::dl1mu_pb() {
  return m_dl1mu_pb;
}

const float Jet::dl1mu_pc() const {
  return m_dl1mu_pc;
}
float &Jet::dl1mu_pc() {
  return m_dl1mu_pc;
}

const float Jet::dl1rnn_pu() const {
  return m_dl1rnn_pu;
}
float &Jet::dl1rnn_pu() {
  return m_dl1rnn_pu;
}

const float Jet::dl1rnn_pb() const {
  return m_dl1rnn_pb;
}
float &Jet::dl1rnn_pb() {
  return m_dl1rnn_pb;
}

const float Jet::dl1rnn_pc() const {
  return m_dl1rnn_pc;
}
float &Jet::dl1rnn_pc() {
  return m_dl1rnn_pc;
}

const int Jet::is_flatbtag() const {
  return m_is_flatbtag;
}
int &Jet::is_flatbtag() {
  return m_is_flatbtag;
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


