#ifndef USERFUNCTIONS_H
#define USERFUNCTIONS_H

//stl
//#include <utility>

//root
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TLorentzVector.h"

//wuto
// #include "WuToObjects/WuToObjectsDict.h"
// #include "WuToObjects/Def.h"

namespace TuDoAtlas{

  //common functions____________________________________________________________________

  //wrappers of root genvector functions________________________________________________

  double angle(TLorentzVector v1, TLorentzVector v2);
  double cos_theta(TLorentzVector v1, TLorentzVector v2);
  double delta_phi(TLorentzVector v1, TLorentzVector v2);
  double invariant_mass(TLorentzVector v1, TLorentzVector v2);
  double delta_r(TLorentzVector v1, TLorentzVector v2);
  double delta_eta(TLorentzVector v1, TLorentzVector v2);

  //some reco_____________________
  double WHelicity(TLorentzVector top4, TLorentzVector W4, TLorentzVector l4);
  
  double Polarization(TLorentzVector lqjet, TLorentzVector lepton, TLorentzVector top);

  double mT_W(TLorentzVector l4, TLorentzVector nu4);

  double m_nu(TLorentzVector l4, TLorentzVector nu4);

  //triangular cut
  bool PassSlideMeT(double dPhiAnchor, double MeTAnchor, double eldPhiMet, double met);
  
  //aplanarity + spherisity
  std::pair<double, double> calc_apl_shere(std::vector<TLorentzVector> jets = std::vector<TLorentzVector>(),
                                           std::vector<TLorentzVector> leps = std::vector<TLorentzVector>(),
                                           std::vector<TLorentzVector> mets = std::vector<TLorentzVector>());

  double aplanarity(std::vector<TLorentzVector> jets = std::vector<TLorentzVector>(),
                    std::vector<TLorentzVector> leps = std::vector<TLorentzVector>(),
                    std::vector<TLorentzVector> mets = std::vector<TLorentzVector>());

  double aplanarity2(std::vector<TLorentzVector> jets = std::vector<TLorentzVector>(),
                    TLorentzVector lep = TLorentzVector(),
                    TLorentzVector met = TLorentzVector());
  
  double aplanarity3(TLorentzVector bjet = TLorentzVector(),
                     TLorentzVector nonbjet = TLorentzVector(),
                     TLorentzVector lep = TLorentzVector(),
                     TLorentzVector met = TLorentzVector());

  double spherisity(std::vector<TLorentzVector> jets = std::vector<TLorentzVector>(),
                    std::vector<TLorentzVector> leps = std::vector<TLorentzVector>(),
                    std::vector<TLorentzVector> mets = std::vector<TLorentzVector>());

  double spherisity2(std::vector<TLorentzVector> jets = std::vector<TLorentzVector>(),
                    TLorentzVector lep = TLorentzVector(),
                    TLorentzVector met = TLorentzVector());
  
  double spherisity3(TLorentzVector bjet = TLorentzVector(),
                     TLorentzVector nonbjet = TLorentzVector(),
                     TLorentzVector lep = TLorentzVector(),
                     TLorentzVector met = TLorentzVector());
  
  double jetprobRND();
  
  double ht(std::vector<TLorentzVector>& leps,
        std::vector<TLorentzVector>& jets,
        std::vector<TLorentzVector>& met);
  
  double ht2(TLorentzVector bjet = TLorentzVector(),
            TLorentzVector nonbjet = TLorentzVector(),
            TLorentzVector lep = TLorentzVector(),
            TLorentzVector met = TLorentzVector());

}
#endif //USERFUNCTIONS_H
