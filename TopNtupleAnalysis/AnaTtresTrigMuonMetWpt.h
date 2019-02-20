/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef ANATTRESTRIGMUONMETWPT_H
#define ANATTRESTRIGMUONMETWPT_H

#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TopNtupleAnalysis/Event.h"
//#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnalysisTrig.h"

#include "TopNtupleAnalysis/TtresNeutrinoBuilder.h"
#include "TopNtupleAnalysis/TtresChi2.h"


class AnaTtresTrigMuonMetWpt : public AnalysisTrig {
  public:
    AnaTtresTrigMuonMetWpt(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList);
    virtual ~AnaTtresTrigMuonMetWpt();

    void run(const Event &e, double weight, const std::string &syst, int is2016run, int trthr);
    void terminate() {};
    void setIsData(bool isData) {};
   
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
};

#endif

