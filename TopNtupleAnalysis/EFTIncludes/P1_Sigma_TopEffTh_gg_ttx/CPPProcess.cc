//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

//#include "CPPProcess.h"
////#include "HelAmps_TopEffTh.h"

using namespace MG5_TopEffTh; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ NP<=2 @1

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_TopEffTh::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->mdl_MT); 
  jamp2[0] = new double[3]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
    pars->printDependentParameters(); 
    pars->printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 3; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 16; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  //std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1,
      -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1},
      {-1, 1, 1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1,
      1, -1}, {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1,
      1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {256}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_1_gg_ttx(); 

        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_1_gg_ttx(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 21 && id2 == 21)
  {
    // Add matrix elements for processes with beams (21, 21)
    return matrix_element[0]; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  //int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  VVS2_3(w[0], w[1], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[4]); 
  VVV2P0_1(w[0], w[1], pars->GC_28, pars->ZERO, pars->ZERO, w[5]); 
  VVV1P0_1(w[0], w[1], pars->GC_6, pars->ZERO, pars->ZERO, w[6]); 
  FFV5_1(w[2], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[7]); 
  FFV3_8_1(w[2], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[8]);
  FFV5_2(w[3], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[9]); 
  FFV3_8_2(w[3], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[10]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFVV1_2_0(w[3], w[2], w[0], w[1], pars->GC_146, pars->GC_85, amp[0]); 
  FFS2_0(w[3], w[2], w[4], pars->GC_95, amp[1]); 
  FFV5_0(w[3], w[2], w[5], pars->GC_7, amp[2]); 
  FFV5_0(w[3], w[2], w[6], pars->GC_7, amp[3]); 
  FFV3_8_0(w[3], w[2], w[6], pars->GC_145, pars->GC_79, amp[4]); 
  FFV5_0(w[3], w[7], w[1], pars->GC_7, amp[5]); 
  FFV3_8_0(w[3], w[7], w[1], pars->GC_145, pars->GC_79, amp[6]); 
  FFV5_0(w[3], w[8], w[1], pars->GC_7, amp[7]); 
  FFV5_0(w[9], w[2], w[1], pars->GC_7, amp[8]); 
  FFV3_8_0(w[9], w[2], w[1], pars->GC_145, pars->GC_79, amp[9]); 
  FFV5_0(w[10], w[2], w[1], pars->GC_7, amp[10]); 

}
double CPPProcess::matrix_1_gg_ttx() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 11; 
  const int ncolor = 3; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {3, 3, 1}; 
  static const double cf[ncolor][ncolor] = {{16, -2, 6}, {-2, 16, 6}, {2, 2,
      6}};

  // Calculate color flows
  jamp[0] = +std::complex<double> (0, 1) * amp[0] + std::complex<double> (0, 1)
      * amp[2] + std::complex<double> (0, 1) * amp[3] + std::complex<double>
      (0, 1) * amp[4] - amp[5] - amp[6] - amp[7];
  jamp[1] = -std::complex<double> (0, 1) * amp[0] - std::complex<double> (0, 1)
      * amp[2] - std::complex<double> (0, 1) * amp[3] - std::complex<double>
      (0, 1) * amp[4] - amp[8] - amp[9] - amp[10];
  jamp[2] = +2. * (-amp[1]); 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



