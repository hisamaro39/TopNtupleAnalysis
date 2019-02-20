/**
 * @brief Jet representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef JET_H
#define JET_H

#include "TopNtupleAnalysis/MObject.h"
#include "TLorentzVector.h"
#include "TSpline.h"
#include "TFile.h"

class Jet : public MObject {
  public:
    Jet();
    Jet(const TLorentzVector &v);
    virtual ~Jet();

    //TFile *m_inf;
    //TSpline *m_spline;

    int &trueFlavour();
    const int trueFlavour() const;

    int &numConstituents();
    const int numConstituents() const;

    const float mv1() const;
    float &mv1();
    
    const float mv2c10() const;
    float &mv2c10();

    const float mv2c10mu() const;
    float &mv2c10mu();

    const float mv2c10rnn() const;
    float &mv2c10rnn();

    const float dl1_pu() const;
    float &dl1_pu();

    const float dl1_pb() const;
    float &dl1_pb();

    const float dl1_pc() const;
    float &dl1_pc();

    const float dl1mu_pu() const;
    float &dl1mu_pu();

    const float dl1mu_pb() const;
    float &dl1mu_pb();

    const float dl1mu_pc() const;
    float &dl1mu_pc();

    const float dl1rnn_pu() const;
    float &dl1rnn_pu();

    const float dl1rnn_pb() const;
    float &dl1rnn_pb();

    const float dl1rnn_pc() const;
    float &dl1rnn_pc();

    const int is_flatbtag() const;
    int &is_flatbtag();
    
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
    bool btag_mv2c10_70() const;
    bool btag_mv2c10_70_trk() const;
    bool btag_mv2c10_70_trk_rel21() const;
    bool btag_mv2c20_70_trk() const;
    bool btag_mv2c10_70_hyb(float pt) const;

    int &closeToLepton();
    const int closeToLepton() const;

  protected:
    int m_trueflavour;
    float m_mv1;
    float m_mv2c10;
    float m_mv2c10mu;
    float m_mv2c10rnn;
    float m_dl1_pu;
    float m_dl1_pb;
    float m_dl1_pc;
    float m_dl1mu_pu;
    float m_dl1mu_pb;
    float m_dl1mu_pc;
    float m_dl1rnn_pu;
    float m_dl1rnn_pb;
    float m_dl1rnn_pc;
    int m_is_flatbtag;
    float m_mv2c20;
    float m_ip3dsv1;
    float m_jvt;
    int m_closeToLepton;
    int m_numConstituents;

};

#endif
