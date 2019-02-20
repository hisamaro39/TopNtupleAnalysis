/**
 * @brief Analysis class for TU Dortmund to run analysis with 2 jets and 1 lepton
 * @author Hendrik Esch <hendrik.esch@tu-dortmund.de>
 */
#ifndef ANATUDOTTRESBOOSTED_H
#define ANATUDOTTRESBOOSTED_H

#include <string>
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TopNtupleAnalysis/Event.h"
#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/MiniTree.h"
#include "TopNtupleAnalysis/KinematicUtils.h"
#include "TopNtupleAnalysis/UserFunktions.h"

#include "TopNtupleAnalysis/TtresNeutrinoBuilder.h"
#include "TopNtupleAnalysis/TtresChi2.h"
// #include "TuDoBase/HistoListSvc.h"

class AnaTuDoTtresBoosted : public Analysis {
  public:
    AnaTuDoTtresBoosted(const std::string &filename, bool electron, MiniTree *mini);
    virtual ~AnaTuDoTtresBoosted();

    void run(const Event &e, double weight, const std::string &syst);
    void terminate();
    void setIsData(bool isData);

  protected:
    bool m_electron;
    bool m_isData;

    TtresNeutrinoBuilder m_neutrinoBuilder;
    TtresChi2 m_chi2;
//     HistoListSvc *m_histoSvc;
    MiniTree *m_mini;
    
    TH1F *h_TuDoCutFlow;
    
//     HistoList1D *m_leppt;
//     HistoList1D *m_lepeta;
//     HistoList1D *m_lepphi;
//     HistoList1D *m_lepcharge;
//     
//     HistoList1D *m_leadJetPt;
//     HistoList1D *m_leadbJetPt;
//     HistoList1D *m_met;
//     HistoList1D *m_met_phi;
//     
//     HistoList1D *m_closeJetPt;
//     HistoList1D *m_largeJetPt;
//     HistoList1D *m_largeJetM;
//     HistoList1D *m_largeJetSd12;
//     
//     HistoList1D *m_mtlep_boo;
//     
//     HistoList1D *m_mtt;
};

#endif
