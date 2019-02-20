#ifndef KINEMATICUTILS_H
#define KINEMATICUTILS_H

/** 
 ** @author Alison Lister <Alison.Lister@cern.ch>
 ** @author W. H. Bell <W.Bell@cern.ch>
 **
 ** @brief A collection of useful methods to calculate kinematic quantities.
 **
 ** The functions collected in the KinematicUtils namespace can be used directly,
 ** they do not require any initialization or depend on any external settings.
 */
namespace KinematicUtils {
  
  /** A function to calculate delta R between two objects.
   *  @param eta1 eta of the first object.
   *  @param eta2 eta of the second object.
   *  @param phi1 phi of the first object.
   *  @param phi2 phi of the second object.
   *  @return deltaR = sqrt(deltaEta**2 + deltaPhi**2) between the two objects.
   */
  double deltaR(double eta1, double eta2, double phi1, double phi2);
  
  /** A function to calculate delta phi between two objects.
   *  The input angles must use the same interval ([0..2 PI] or [-PI..PI]).
   *  @param phi1 angle phi of the first object.
   *  @param phi2 angle phi of the second object.
   *  @return angle between phi1 and phi2, limited to the interval [0..PI].
   */
  double deltaPhi(double phi1, double phi2);

  /** A function to calculate the transverse mass.
   * @param ptLep lepton pT
   * @param phiLep lepton phi
   * @param met Missing ET
   * @param phiMet the azimuthal angle of the MET vector
   * @return the transverse mass. 
   */
  double transMass(double ptLep, double phiLep, double met, double phiMet);
}

#endif 
