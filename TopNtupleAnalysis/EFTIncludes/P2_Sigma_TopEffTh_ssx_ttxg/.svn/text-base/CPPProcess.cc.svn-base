//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

//#include "CPPProcess.h"
//#include "HelAmps_TopEffTh.h"

using namespace MG5_TopEffTh; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: s s~ > t t~ g NP<=2 @2

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
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[4]; 
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
  for(int i = 0; i < 4; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 32; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  //std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1}, {-1,
      -1, -1, -1, 1}, {-1, -1, -1, 1, -1}, {-1, -1, -1, 1, 1}, {-1, -1, 1, -1,
      -1}, {-1, -1, 1, -1, 1}, {-1, -1, 1, 1, -1}, {-1, -1, 1, 1, 1}, {-1, 1,
      -1, -1, -1}, {-1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1},
      {-1, 1, 1, -1, -1}, {-1, 1, 1, -1, 1}, {-1, 1, 1, 1, -1}, {-1, 1, 1, 1,
      1}, {1, -1, -1, -1, -1}, {1, -1, -1, -1, 1}, {1, -1, -1, 1, -1}, {1, -1,
      -1, 1, 1}, {1, -1, 1, -1, -1}, {1, -1, 1, -1, 1}, {1, -1, 1, 1, -1}, {1,
      -1, 1, 1, 1}, {1, 1, -1, -1, -1}, {1, 1, -1, -1, 1}, {1, 1, -1, 1, -1},
      {1, 1, -1, 1, 1}, {1, 1, 1, -1, -1}, {1, 1, 1, -1, 1}, {1, 1, 1, 1, -1},
      {1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {36, 36}; 

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
        t[0] = matrix_2_ssx_ttxg(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_2_ssx_ttxg(); 
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
      t[0] = matrix_2_ssx_ttxg(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_2_ssx_ttxg(); 
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
  if(id1 == -3 && id2 == 3)
  {
    // Add matrix elements for processes with beams (-3, 3)
    return matrix_element[1]; 
  }
  else if(id1 == 3 && id2 == -3)
  {
    // Add matrix elements for processes with beams (3, -3)
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
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  FFV1_2(w[0], w[4], pars->GC_7, pars->ZERO, pars->ZERO, w[5]); 
  FFV5P0_3(w[3], w[2], pars->GC_2, pars->ZERO, pars->ZERO, w[6]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_170, pars->GC_91, pars->ZERO, pars->ZERO,
      w[7]);
  FFV5P0_3(w[3], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[8]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[9]);
  FFV3_8_3(w[3], w[2], pars->GC_166, pars->GC_81, pars->mdl_MZ, pars->mdl_WZ,
      w[10]);
  FFV2_7_3(w[3], w[2], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[11]);
  FFV2_3(w[3], w[2], pars->GC_138, pars->mdl_MZ, pars->mdl_WZ, w[12]); 
  FFV5_1(w[2], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[13]); 
  FFV1P0_3(w[0], w[1], pars->GC_1, pars->ZERO, pars->ZERO, w[14]); 
  FFV1P0_3(w[0], w[1], pars->GC_7, pars->ZERO, pars->ZERO, w[15]); 
  FFV2_4_3(w[0], w[1], pars->GC_58, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[16]);
  FFV3_8_1(w[2], w[4], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[17]);
  FFV5_2(w[3], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[18]); 
  FFV3_8_2(w[3], w[4], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[19]);
  FFV1_1(w[1], w[4], pars->GC_7, pars->ZERO, pars->ZERO, w[20]); 
  FFS2_3(w[3], w[2], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[21]); 
  FFVV1_2P0_3(w[3], w[2], w[4], pars->GC_146, pars->GC_85, pars->ZERO,
      pars->ZERO, w[22]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[5], w[1], w[6], pars->GC_1, amp[0]); 
  FFV1_0(w[5], w[1], w[7], pars->GC_1, amp[1]); 
  FFV1_0(w[5], w[1], w[8], pars->GC_7, amp[2]); 
  FFV1_0(w[5], w[1], w[9], pars->GC_7, amp[3]); 
  FFV2_4_0(w[5], w[1], w[10], pars->GC_58, pars->GC_69, amp[4]); 
  FFV2_4_0(w[5], w[1], w[11], pars->GC_58, pars->GC_69, amp[5]); 
  FFV2_4_0(w[5], w[1], w[12], pars->GC_58, pars->GC_69, amp[6]); 
  FFFF4_5_0(w[5], w[2], w[3], w[1], pars->GC_15, pars->GC_16, amp[7]); 
  FFFF1_0(w[5], w[2], w[3], w[1], -pars->GC_12, amp[8]); 
  FFFF1_6_0(w[5], w[2], w[3], w[1], -pars->GC_10, -pars->GC_26, amp[9]); 
  FFV5_0(w[3], w[13], w[14], pars->GC_2, amp[10]); 
  FFV3_8_0(w[3], w[13], w[14], pars->GC_170, pars->GC_91, amp[11]); 
  FFV5_0(w[3], w[13], w[15], pars->GC_7, amp[12]); 
  FFV3_8_0(w[3], w[13], w[15], pars->GC_145, pars->GC_79, amp[13]); 
  FFV3_8_0(w[3], w[13], w[16], pars->GC_166, pars->GC_81, amp[14]); 
  FFV2_7_0(w[3], w[13], w[16], pars->GC_59, pars->GC_69, amp[15]); 
  FFV2_0(w[3], w[13], w[16], pars->GC_138, amp[16]); 
  FFV5_0(w[3], w[17], w[14], pars->GC_2, amp[17]); 
  FFV5_0(w[3], w[17], w[15], pars->GC_7, amp[18]); 
  FFV2_7_0(w[3], w[17], w[16], pars->GC_59, pars->GC_69, amp[19]); 
  FFFF4_5_0(w[0], w[13], w[3], w[1], pars->GC_15, pars->GC_16, amp[20]); 
  FFFF1_0(w[0], w[13], w[3], w[1], -pars->GC_12, amp[21]); 
  FFFF1_6_0(w[0], w[13], w[3], w[1], -pars->GC_10, -pars->GC_26, amp[22]); 
  FFV5_0(w[18], w[2], w[14], pars->GC_2, amp[23]); 
  FFV3_8_0(w[18], w[2], w[14], pars->GC_170, pars->GC_91, amp[24]); 
  FFV5_0(w[18], w[2], w[15], pars->GC_7, amp[25]); 
  FFV3_8_0(w[18], w[2], w[15], pars->GC_145, pars->GC_79, amp[26]); 
  FFV3_8_0(w[18], w[2], w[16], pars->GC_166, pars->GC_81, amp[27]); 
  FFV2_7_0(w[18], w[2], w[16], pars->GC_59, pars->GC_69, amp[28]); 
  FFV2_0(w[18], w[2], w[16], pars->GC_138, amp[29]); 
  FFV5_0(w[19], w[2], w[14], pars->GC_2, amp[30]); 
  FFV5_0(w[19], w[2], w[15], pars->GC_7, amp[31]); 
  FFV2_7_0(w[19], w[2], w[16], pars->GC_59, pars->GC_69, amp[32]); 
  FFFF4_5_0(w[0], w[2], w[18], w[1], pars->GC_15, pars->GC_16, amp[33]); 
  FFFF1_0(w[0], w[2], w[18], w[1], -pars->GC_12, amp[34]); 
  FFFF1_6_0(w[0], w[2], w[18], w[1], -pars->GC_10, -pars->GC_26, amp[35]); 
  FFV1_0(w[0], w[20], w[6], pars->GC_1, amp[36]); 
  FFV1_0(w[0], w[20], w[7], pars->GC_1, amp[37]); 
  FFV1_0(w[0], w[20], w[8], pars->GC_7, amp[38]); 
  FFV1_0(w[0], w[20], w[9], pars->GC_7, amp[39]); 
  FFV2_4_0(w[0], w[20], w[10], pars->GC_58, pars->GC_69, amp[40]); 
  FFV2_4_0(w[0], w[20], w[11], pars->GC_58, pars->GC_69, amp[41]); 
  FFV2_4_0(w[0], w[20], w[12], pars->GC_58, pars->GC_69, amp[42]); 
  FFFF4_5_0(w[0], w[2], w[3], w[20], pars->GC_15, pars->GC_16, amp[43]); 
  FFFF1_0(w[0], w[2], w[3], w[20], -pars->GC_12, amp[44]); 
  FFFF1_6_0(w[0], w[2], w[3], w[20], -pars->GC_10, -pars->GC_26, amp[45]); 
  VVS2_0(w[4], w[15], w[21], pars->GC_78, amp[46]); 
  VVV2_0(w[4], w[15], w[8], pars->GC_28, amp[47]); 
  VVV1_0(w[4], w[15], w[8], pars->GC_6, amp[48]); 
  VVV1_0(w[4], w[15], w[9], pars->GC_6, amp[49]); 
  FFV1_0(w[0], w[1], w[22], pars->GC_7, amp[50]); 

}
double CPPProcess::matrix_2_ssx_ttxg() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 51; 
  const int ncolor = 4; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1, 1, 1, 1}; 
  static const double cf[ncolor][ncolor] = {{12, 4, 4, 0}, {4, 12, 0, 4}, {4,
      0, 12, 4}, {0, 4, 4, 12}};

  // Calculate color flows
  jamp[0] = -amp[10] - amp[11] + 1./6. * amp[12] + 1./6. * amp[13] - amp[14] -
      amp[15] - amp[16] - amp[17] + 1./6. * amp[18] - amp[19] + amp[21] - 1./6.
      * amp[22] - amp[23] - amp[24] + 1./6. * amp[25] + 1./6. * amp[26] -
      amp[27] - amp[28] - amp[29] - amp[30] + 1./6. * amp[31] - amp[32] +
      amp[34] - 1./6. * amp[35];
  jamp[1] = -1./2. * amp[2] - 1./2. * amp[3] + amp[7] + 1./2. * amp[9] - 1./2.
      * amp[12] - 1./2. * amp[13] - 1./2. * amp[18] + amp[20] + 1./2. * amp[22]
      + 1./2. * std::complex<double> (0, 1) * amp[47] + 1./2. *
      std::complex<double> (0, 1) * amp[48] + 1./2. * std::complex<double> (0,
      1) * amp[49] - 1./2. * std::complex<double> (0, 1) * amp[50];
  jamp[2] = -1./2. * amp[25] - 1./2. * amp[26] - 1./2. * amp[31] + amp[33] +
      1./2. * amp[35] - 1./2. * amp[38] - 1./2. * amp[39] + amp[43] + 1./2. *
      amp[45] - 1./2. * std::complex<double> (0, 1) * amp[47] - 1./2. *
      std::complex<double> (0, 1) * amp[48] - 1./2. * std::complex<double> (0,
      1) * amp[49] + 1./2. * std::complex<double> (0, 1) * amp[50];
  jamp[3] = -amp[0] - amp[1] + 1./6. * amp[2] + 1./6. * amp[3] - amp[4] -
      amp[5] - amp[6] + amp[8] - 1./6. * amp[9] - amp[36] - amp[37] + 1./6. *
      amp[38] + 1./6. * amp[39] - amp[40] - amp[41] - amp[42] + amp[44] - 1./6.
      * amp[45] - amp[46];

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



