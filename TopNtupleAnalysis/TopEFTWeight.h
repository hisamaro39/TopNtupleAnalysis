#ifndef TOPEFTWEIGHT_H
#define TOPEFTWEIGHT_H
// Using EFT for the top production at LO
// Based on the following papers
// Please cite them if used.
//
// %\cite{Degrande:2010kt}
// \bibitem{Degrande:2010kt}
//   C.~Degrande, J.~M.~Gerard, C.~Grojean, F.~Maltoni and G.~Servant,
//   %``Non-resonant New Physics in Top Pair Production at Hadron Colliders,''
//   JHEP {\bf 1103} (2011) 125
//   doi:10.1007/JHEP03(2011)125
//   [arXiv:1010.6304 [hep-ph]].
//   %%CITATION = doi:10.1007/JHEP03(2011)125;%%
//
// %\cite{Zhang:2010dr}
// \bibitem{Zhang:2010dr}
//   C.~Zhang and S.~Willenbrock,
//   %``Effective-Field-Theory Approach to Top-Quark Production and Decay,''
//   Phys.\ Rev.\ D {\bf 83} (2011) 034006
//   doi:10.1103/PhysRevD.83.034006
//   [arXiv:1008.3869 [hep-ph]].
//   %%CITATION = doi:10.1103/PhysRevD.83.034006;%%


#include "TLorentzVector.h"
#include <map>

class TopEFTWeight { 
  
public:
    TopEFTWeight(double cvv = 1.0, double lambda = 5000); // lambda in GeV    
    ~TopEFTWeight();
    std::pair<double, double> getScaleFactor(
        float  MC_t_pt,
        float  MC_t_eta,
        float  MC_t_phi,
        float  MC_t_m,
        float  MC_tbar_pt,
        float  MC_tbar_eta,
        float  MC_tbar_phi,
        float  MC_tbar_m,
        float  MC_i1_px,
        float  MC_i1_py,
        float  MC_i1_pz,
        float  MC_i1_m,
        float  MC_i2_px,
        float  MC_i2_py,
        float  MC_i2_pz,
        float  MC_i2_m,
	int    MC_i1_pid,
	int    MC_i2_pid);
private:
    double m_cvv;
    double m_caa;
    double m_cvvp;
    double m_caap;
    double m_lambda;

    double m2(TLorentzVector a, TLorentzVector b);
    double m2diff(TLorentzVector a, TLorentzVector b);
    double alphas(double scale);
};

#endif //TOPEFTWEIGHT_H
