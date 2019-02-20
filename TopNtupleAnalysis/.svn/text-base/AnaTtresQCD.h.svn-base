/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef ANATTRESQCD_H
#define ANATTRESQCD_H

#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TopNtupleAnalysis/Event.h"
#include "TopNtupleAnalysis/Analysis.h"

#include "TopNtupleAnalysis/TtresNeutrinoBuilder.h"
#include "TopNtupleAnalysis/TtresChi2.h"

class AnaTtresQCD : public Analysis {
  public:
    AnaTtresQCD(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList, int dsid);
    virtual ~AnaTtresQCD();

    void run(const Event &e, double weight, const std::string &syst);
        
    virtual void runRealRateWQCDCR_2015(const Event &e, double weight, const std::string &suffix);    
    virtual void runRealRateWQCDCR_2016(const Event &e, double weight, const std::string &suffix);
    virtual void runRealRateQCDCR_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runRealRateQCDCR_2016(const Event &e, double weight, const std::string &suffix);
    virtual void GetRealHistograms(const Event &evt, const double weight, const std::string &suffix, const std::string &btag);
 
    virtual void runFakeRateQCDCR_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runFakeRateWQCDCR_2015(const Event &e, double weight, const std::string &suffix);
    virtual void runFakeRateQCDCR_2016(const Event &e, double weight, const std::string &suffix);
    virtual void runFakeRateWQCDCR_2016(const Event &e, double weight, const std::string &suffix);
    virtual void get1Drates(float &rate, float &rate_err, TH1F* rate_map, float x);
    virtual void get2Drates(float &rate, float &rate_err, TH2F* rate_map, float x, float y);
    virtual void GetFakeHistograms(const Event &e, double weight, const std::string &suffix_corr, const std::string &suffix, const std::string &btag);

    virtual void IniHistograms(std::string &suffix,  std::string &btag);
 
    void terminate() {};
    void setIsData(bool isData) {};

  protected:
    bool m_electron;
    bool m_boosted;

    TtresNeutrinoBuilder m_neutrinoBuilder;
    TtresChi2 m_chi2;

    int m_dsid;
};

#endif

