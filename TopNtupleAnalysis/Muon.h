/**
 * @brief Muon representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef MUON_H
#define MUON_H

#include "TopNtupleAnalysis/MObject.h"
#include "TLorentzVector.h"

class Muon : public MObject {
  public:
    Muon();
    Muon(const TLorentzVector &v);
    virtual ~Muon();

    void setMI(float iso);
    float mi() const;

    bool isTight() const;
    void setTight(char t);

    const float Dz0() const;
    float &Dz0();

    const float z0_exPV() const;
    float &z0_exPV();

    const float d0() const;
    float &d0();

    const float sd0() const;
    float &sd0();
    
    const bool HLT_mu50() const;
    bool &HLT_mu50();

    const bool HLT_mu26_ivarmedium() const;
    bool &HLT_mu26_ivarmedium();
    
    const bool HLT_mu20_iloose_L1MU15() const;
    bool &HLT_mu20_iloose_L1MU15();
    
    const bool HLT_mu20_L1MU15() const;
    bool &HLT_mu20_L1MU15();

    const int author() const;
    int &author();

    const int quality() const;
    int &quality();

    const int accept() const;
    int &accept();

    const int trigger_match() const;
    int &trigger_match();

    const bool passTrkCuts() const;
    bool &passTrkCuts();
    
    const float ptvarcone30() const;
    float &ptvarcone30();
    
    const float topoetcone20() const;
    float &topoetcone20();
    
    bool pass() const;
    bool passLoose() const;

    TLorentzVector &momME();
    const TLorentzVector &momME() const;

    TLorentzVector &momMECorr();
    const TLorentzVector &momMECorr() const;

    TLorentzVector &momMS();
    const TLorentzVector &momMS() const;

    TLorentzVector &momID();
    const TLorentzVector &momID() const;

    TLorentzVector &momTrk();
    const TLorentzVector &momTrk() const;

    int charge() const;
    int &charge();

    void setST(bool b);
    bool st() const;

    void setSA(bool b);
    bool sa() const;

    void setCB(bool b);
    bool cb() const;

  protected:
    float m_mi;
    bool m_tight;
    float m_Dz0;
    float m_z0_exPV;
    float m_d0;
    float m_sd0;
    float m_ptvarcone30;
    float m_topoetcone20;
    
    bool m_HLT_mu50;
    bool m_HLT_mu20_iloose_L1MU15;
    bool m_HLT_mu20_L1MU15;
    bool m_HLT_mu26_ivarmedium;
    
    int m_author;
    int m_quality;
    int m_accept;
    int m_trigger_match;
    bool m_passTrkCuts;

    TLorentzVector m_momME;
    TLorentzVector m_momMECorr;
    TLorentzVector m_momMS;
    TLorentzVector m_momID;
    TLorentzVector m_momTrk;

    int m_charge;
    bool m_st;
    bool m_sa;
    bool m_cb;

};

#endif
