/**
 * @brief Truth representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef TRUTH_H
#define TRUTH_H

#include "TopNtupleAnalysis/MObject.h"
#include "TLorentzVector.h"

class Truth : public MObject {
  public:
    Truth();
    Truth(const TLorentzVector &v);
    virtual ~Truth();

    const float px() const;
    float &px();
    const float py() const;
    float &py();
    const float pz() const;
    float &pz();
    const float pt() const;
    float &pt();
    const float e() const;
    float &e();
    const float m() const;
    float &m();
    const int id() const;
    int &id();
    const int status() const;
    int &status();
    const int muon_type() const;
    int &muon_type();

  protected:
    float m_px;
    float m_py;
    float m_pz;
    float m_e;
    float m_m;
    int m_id;
    int m_status;
    int m_muon_type;

};

#endif
