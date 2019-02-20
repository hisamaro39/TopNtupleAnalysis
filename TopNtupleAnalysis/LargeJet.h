/**
 * @brief Large-R jet representation to be read off the input file.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef LARGEJET_H
#define LARGEJET_H

#include "TopNtupleAnalysis/MObject.h"
#include "TLorentzVector.h"

class LargeJet : public MObject {
  public:
    LargeJet();
    LargeJet(const TLorentzVector &v);
    LargeJet(const LargeJet &l);
    virtual ~LargeJet();

    int &trueFlavour();
    const int trueFlavour() const;

    bool pass() const;
    bool passLoose() const;
    bool passFakeSelection(const TLorentzVector &lept, const TLorentzVector &selJet) const;
    
    double &split12();
    const double split12() const;

    double &tau32();
    const double tau32() const;

    double &BDT_TOPtag();
    const double BDT_TOPtag() const;

    bool good() const;
    bool &good();
    void setGood(bool);

    bool good50() const;
    bool &good50();
    void setGood50(bool);

    bool good_sub80() const;
    bool &good_sub80();

    bool good_sub50() const;
    bool &good_sub50();

    bool good_smooth_mt80() const;
    bool &good_smooth_mt80();

    bool good_smooth_mt50() const;
    bool &good_smooth_mt50();

    bool good_smooth_ts80() const;
    bool &good_smooth_ts80();

    bool good_smooth_ts50() const;
    bool &good_smooth_ts50();

    bool good_smooth_qt80() const;
    bool &good_smooth_qt80();

    bool good_smooth_qt50() const;
    bool &good_smooth_qt50();
    
    bool good_bdt80() const;
    bool &good_bdt80();

    bool good_dnn80() const;
    bool &good_dnn80();

    float &subs(const std::string &s);
    const float subs(const std::string &s) const;

  protected:
    double m_split12;
    double m_tau32;
    double m_BDT_TOPtag;

    int m_trueFlavour;
    bool m_good;
    bool m_good50;
    bool m_good_sub80;
    bool m_good_sub50;
    bool m_good_smooth_mt80;
    bool m_good_smooth_mt50;
    bool m_good_smooth_ts80;
    bool m_good_smooth_ts50;
    bool m_good_smooth_qt80;
    bool m_good_smooth_qt50;
    bool m_good_bdt80;
    bool m_good_dnn80;

    std::map<std::string, float> m_subs;
};

#endif
