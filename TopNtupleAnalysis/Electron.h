/**
 * @brief Electron representation for information read off input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef ELECTRON_H
#define ELECTRON_H

#include "TopNtupleAnalysis/MObject.h"
#include "TLorentzVector.h"

class Electron : public MObject {
  public:
    Electron();
    Electron(const TLorentzVector &v);
    virtual ~Electron();

    void setMI(float iso);
    float mi() const;
    
    void setTightPP(char isTightPP);
    bool isTightPP() const;

    void setMediumPP(bool isMediumPP);
    bool isMediumPP() const;

    const TLorentzVector &caloMom() const;
    TLorentzVector &caloMom();
    const TLorentzVector &trkMom() const;
    TLorentzVector &trkMom();
    const float Dz0() const;
    float &Dz0();
    
    const float d0() const;
    float &d0();

    const float sd0() const;
    float &sd0();
    
    const int author() const;
    int &author();

    int &nSiHits();
    const int nSiHits() const;

    int &oq();
    const int oq() const;

    int &isEM();
    const int isEM() const;

    int &goodOQ();
    const int goodOQ() const;

    int &passLH();
    const int passLH() const;
    
    const float ptvarcone20() const;
    float &ptvarcone20();
    
    const float topoetcone20() const;
    float &topoetcone20();
    
    const bool HLT_e24_lhmedium_iloose_L1EM20VH() const;
    bool &HLT_e24_lhmedium_iloose_L1EM20VH();
    
    const bool HLT_e24_lhmedium_L1EM18VH() const;
    bool &HLT_e24_lhmedium_L1EM18VH();
    
    const bool HLT_e24_lhmedium_L1EM20VH() const;
    bool &HLT_e24_lhmedium_L1EM20VH();
    
    const bool HLT_e60_lhmedium() const;
    bool &HLT_e60_lhmedium();
    
    const bool HLT_e120_lhloose() const;  
    bool &HLT_e120_lhloose();

    const bool HLT_e26_lhtight_nod0_ivarloose() const;
    bool &HLT_e26_lhtight_nod0_ivarloose();

    const bool HLT_e60_lhmedium_nod0() const;  
    bool &HLT_e60_lhmedium_nod0();

    const bool HLT_e140_lhloose_nod0() const;  
    bool &HLT_e140_lhloose_nod0();
    
    bool pass() const;
    bool passLoose() const;

  protected:
    float m_mi;
    bool m_isTightPP;
    
    bool m_isMediumPP;
    TLorentzVector m_mom_calo;
    TLorentzVector m_mom_trk;
    float m_Dz0;
    float m_d0;
    float m_sd0;
    float m_ptvarcone20;
    float m_topoetcone20;
    int m_author;

    int m_nSiHits;
    int m_oq;
    int m_isEM;
    int m_goodOQ;
    int m_passLH;
    
    bool m_HLT_e24_lhmedium_iloose_L1EM20VH;    
    bool m_HLT_e24_lhmedium_L1EM18VH;
    bool m_HLT_e24_lhmedium_L1EM20VH;
    bool m_HLT_e60_lhmedium;
    bool m_HLT_e120_lhloose;
    bool m_HLT_e26_lhtight_nod0_ivarloose;
    bool m_HLT_e60_lhmedium_nod0;
    bool m_HLT_e140_lhloose_nod0;
    
    int m_GSF_trk_index;
};

#endif
