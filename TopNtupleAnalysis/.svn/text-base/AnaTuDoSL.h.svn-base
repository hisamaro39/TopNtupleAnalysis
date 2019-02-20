/**
 * @brief Analysis class for TU Dortmund to run analysis with 2 jets and 1 lepton
 * @author Hendrik Esch <hendrik.esch@tu-dortmund.de>
 */
#ifndef ANATUDOSL_H
#define ANATUDOSL_H

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

class AnaTuDoSL : public Analysis {
  public:
    AnaTuDoSL(const std::string &filename, bool electron, MiniTree *mini);
    virtual ~AnaTuDoSL();

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
    
//     HistoList2D *m_lep_DeltaPhi;
//     HistoList2D *m_lep_DeltaPhiCut;
//     
//     HistoList1D *m_leppt;
//     HistoList1D *m_lepeta;
//     HistoList1D *m_lepphi;
//     HistoList1D *m_lepcharge;
//     HistoList1D *m_lepmetphi;
//     
//     HistoList1D *m_nupt;
//     HistoList1D *m_nueta;
//     HistoList1D *m_nuphi;
//     
//     HistoList1D *m_jpt_1;
//     HistoList1D *m_jeta_1;
//     HistoList1D *m_jphi_1;
//     HistoList1D *m_jpt_2;
//     HistoList1D *m_jeta_2;
//     HistoList1D *m_jphi_2;
//     
//     HistoList1D *m_bpt;
//     HistoList1D *m_beta;
//     HistoList1D *m_bphi;
//     HistoList1D *m_lightpt;
//     HistoList1D *m_lighteta;
//     HistoList1D *m_lightphi;
//        
//     HistoList1D *m_met;
//     HistoList1D *m_metx;
//     HistoList1D *m_mety;
//     HistoList1D *m_metphi;
//     HistoList1D *m_mtw;
//     
//     HistoList1D *m_deltaR_j1_lnu;
//     HistoList1D *m_ht;
//     HistoList1D *m_deltaR_j_b;
//     HistoList1D *m_mlnub;
//     HistoList1D *m_mb;
//     HistoList1D *m_mjb;
//     HistoList1D *m_etalnu;
//     HistoList1D *m_etalnub;
//     HistoList1D *m_mlb;
//     HistoList1D *m_deltaR_l_lnub;
//     HistoList1D *m_pt_lnu;
//     HistoList1D *m_pt_lnub;
//     HistoList1D *m_mnu;
//     HistoList1D *m_costheta_lj_top;
//     HistoList1D *m_WHelicity;
//     HistoList1D *m_aplanarity;
//     HistoList1D *m_sphericity;
//     HistoList1D *m_mlnubj;
//     HistoList1D *m_etaj;
//     HistoList1D *m_deltaEta_l_b;
//     HistoList1D *m_deltaPhi_l_b;
//     HistoList1D *m_deltaR_l_b;
};

#endif
