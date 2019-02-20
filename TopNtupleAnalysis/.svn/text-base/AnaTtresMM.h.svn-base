/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef ANATTRESMM_H
#define ANATTRESMM_H

#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TopNtupleAnalysis/Event.h"
#include "TopNtupleAnalysis/Analysis.h"

#include "TopNtupleAnalysis/TtresNeutrinoBuilder.h"
#include "TopNtupleAnalysis/TtresChi2.h"

class AnaTtresMM : public Analysis {
  public:
    AnaTtresMM(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList, int dsid);
    virtual ~AnaTtresMM();

    void run(const Event &e, double weight, const std::string &syst);
    void terminate() {};
    void setIsData(bool isData) {};
    
    virtual void runMatrixMethod_QCDVR1_2j_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDVR2_2j_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDCR2j_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDCR4j_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_WjetsCR2j_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_SR4j_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDVR1_4j_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDVR2_4j_2015(const Event &e, double weight, const std::string &suffix);

    virtual void runMatrixMethod_QCDVR1_2j_2016(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDVR2_2j_2016(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDCR2j_2016(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDSR2j_2016(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDVR1_4j_2016(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDVR2_4j_2016(const Event &e, double weight, const std::string &suffix);
    virtual void runMatrixMethod_QCDCR4j_2016(const Event &e, double weight, const std::string &suffix);

    virtual void GetHistograms(const Event &evt, const double weight,const std::string &prefix, const std::string &suffix);
    virtual void IniHistograms(std::string &suffix);
    
    inline bool isElectron(){return m_electron;};
    inline bool isBoosted(){return m_boosted;};
    
  protected:
    bool m_electron;
    bool m_boosted;

    TtresNeutrinoBuilder m_neutrinoBuilder;
    TtresChi2 m_chi2;

    double _tree_truemtt;
    double _tree_mtt;
    double _tree_weight;
    int     _tree_cat;
    std::string _tree_syst;

    int m_dsid;
};

#endif

