//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_TopEffTh_H
#define Parameters_TopEffTh_H

#include <complex> 

#include "TopNtupleAnalysis/read_slha.h"
using namespace std; 

class Parameters_TopEffTh
{
  public:

    static Parameters_TopEffTh * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_ymb, aS,
        mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB, mdl_C1qt,
        mdl_C1qd, mdl_C1qu, mdl_C8dt, mdl_C8ut, mdl_C83qq, mdl_C81qq,
        mdl_C13qq, mdl_CphiG, mdl_CG, mdl_ICtG, mdl_RCtG, mdl_ICtW, mdl_RCtW,
        mdl_IC3phiq, mdl_RC3phiq, mdl_Lambda, mdl_MZ__exp__2, mdl_MZ__exp__4,
        mdl_sqrt__2, mdl_MH__exp__2, mdl_Lambda__exp__2, mdl_aEW, mdl_MW,
        mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2,
        mdl_sw, mdl_g1, mdl_gw, mdl_vev, mdl_vev__exp__2, mdl_lam, mdl_yb,
        mdl_yt, mdl_ytau, mdl_muH, mdl_ee__exp__2, mdl_sw__exp__2,
        mdl_cw__exp__2;
    std::complex<double> mdl_complexi, mdl_C3phiq, mdl_CtG, mdl_CtW,
        mdl_conjg__C3phiq, mdl_conjg__CtG, mdl_conjg__CtW;
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G, mdl_G__exp__2, mdl_G__exp__3; 
    // Model couplings independent of aS
    std::complex<double> GC_1, GC_2, GC_4, GC_57, GC_58, GC_59, GC_60, GC_69,
        GC_78, GC_79, GC_80, GC_81, GC_87, GC_89, GC_91, GC_92, GC_93, GC_95,
        GC_108, GC_128, GC_138, GC_145, GC_165, GC_166, GC_168, GC_170, GC_10,
        GC_11, GC_12, GC_13, GC_14, GC_15, GC_16, GC_17, GC_23, GC_24, GC_25,
        GC_26, GC_27, GC_28;
    // Model couplings dependent on aS
    std::complex<double> GC_146, GC_7, GC_6, GC_8, GC_50, GC_49, GC_44, GC_84,
        GC_85, GC_86;

    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_TopEffTh * instance; 
}; 

#endif  // Parameters_TopEffTh_H

