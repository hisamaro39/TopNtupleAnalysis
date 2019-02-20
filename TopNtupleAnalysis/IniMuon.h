/**
 * @brief IniMuon representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef INIMUON_H
#define INIMUON_H

#include "TopNtupleAnalysis/MObject.h"
#include "TLorentzVector.h"

class IniMuon : public MObject {
  public:
    IniMuon();
    IniMuon(const TLorentzVector &v);
    virtual ~IniMuon();

    const float pt() const;
    float &pt();

    const float eta() const;
    float &eta();

    const float phi() const;
    float &phi();

    const float charge() const;
    float &charge();

    const float ptvarcone30() const;
    float &ptvarcone30();

    const float d0sig() const;
    float &d0sig();
    
    const float z0sintheta() const;
    float &z0sintheta();

    const int quality() const;
    int &quality();

    const int accept() const;
    int &accept();

    const int trigger_match() const;
    int &trigger_match();

  protected:
    float m_pt;
    float m_eta;
    float m_phi;
    float m_charge;
    float m_ptvarcone30;
    float m_d0sig;
    float m_z0sintheta;
    int m_quality;
    int m_accept;
    int m_trigger_match;

};

#endif
