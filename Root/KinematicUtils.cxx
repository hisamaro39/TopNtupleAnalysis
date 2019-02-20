#include "TopNtupleAnalysis/KinematicUtils.h"
#include "TMath.h"
#include <cmath>

namespace KinematicUtils {

double deltaR(double eta1, double eta2, double phi1, double phi2) {
  double dPhi=deltaPhi(phi1,phi2);
  double dEta=std::fabs(eta1-eta2);
  double dR=std::sqrt(std::pow(dEta,2)+std::pow(dPhi,2));
  return dR;
}
  
//-------------------------------------------------------------------------------

double deltaPhi(double phi1, double phi2) {
  double dPhi=std::fabs(phi1-phi2);
  if (dPhi>M_PI) dPhi=2*M_PI-dPhi;
  return dPhi;
}

//-------------------------------------------------------------------------------

double transMass(double ptLep, double phiLep, double met, double phiMet) {
  return std::sqrt(2.0*ptLep*met*( 1 - std::cos( phiLep-phiMet ) ) );
}

//-------------------------------------------------------------------------------

void rotateXaxis(float &PX, float &PY, float &PZ, float px, float py, float pz, float theta, float phi){
  theta = TMath::Pi()-fabs(theta);
  float xx = cos(theta)*cos(phi);
  float xy = cos(theta)*sin(phi);
  float xz = sin(theta);
  float yx = -1*sin(phi);
  float yy = cos(phi);
  float yz = 0;
  float zx = -1*sin(theta)*cos(phi);
  float zy = -1*sin(theta)*sin(theta);
  float zz = cos(theta);
  PX = xx*px + xy*py+ xz*pz;
  PY = yx*px + yy*py + yz*pz;
  PZ = -1*(zx*px + zy*py + zz*pz);
}

}
