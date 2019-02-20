/**
 * @brief Jet representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef JET_H
#define JET_H

#include "TopNtupleAnalysis/MObject.h"
#include "TLorentzVector.h"

class Jet : public MObject {
  public:
    Jet();
    Jet(const TLorentzVector &v);
    virtual ~Jet();

    int &trueFlavour();
    const int trueFlavour() const;

    int &numConstituents();
    const int numConstituents() const;

    const float mv1() const;
    float &mv1();
    
    const float mv2c10() const;
    float &mv2c10();
    
    const float mv2c20() const;
    float &mv2c20();
    
    const float ip3dsv1() const;
    float &ip3dsv1();

    const float jvt() const;
    float &jvt();

    bool pass() const;
    bool pass_trk() const;
    bool btag() const;
    bool btag_mv2c20_60() const;
    bool btag_mv2c20_70() const;
    bool btag_mv2c10_70_trk() const;
    bool btag_mv2c20_70_trk() const;

    int &closeToLepton();
    const int closeToLepton() const;

  protected:
    int m_trueflavour;
    float m_mv1;
    float m_mv2c10;
    float m_mv2c20;
    float m_ip3dsv1;
    float m_jvt;
    int m_closeToLepton;
    int m_numConstituents;

};

#endif
