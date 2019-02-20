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
// Process: c~ d~ > t t~ u~ s~ NP<=2 @3

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
  mME.push_back(pars->ZERO); 
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
  const int ncomb = 64; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  //std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
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
        t[0] = matrix_3_cxdx_ttxuxsx(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_3_cxdx_ttxuxsx(); 
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
      t[0] = matrix_3_cxdx_ttxuxsx(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_3_cxdx_ttxuxsx(); 
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
  if(id1 == -1 && id2 == -4)
  {
    // Add matrix elements for processes with beams (-1, -4)
    return matrix_element[1]; 
  }
  else if(id1 == -4 && id2 == -1)
  {
    // Add matrix elements for processes with beams (-4, -1)
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
  oxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[4]); 
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]); 
  FFV2_3(w[4], w[1], pars->GC_57, pars->mdl_MW, pars->mdl_WW, w[6]); 
  FFV2_3(w[5], w[0], pars->GC_57, pars->mdl_MW, pars->mdl_WW, w[7]); 
  FFV3_1(w[2], w[6], pars->GC_165, pars->mdl_MB, pars->ZERO, w[8]); 
  FFV2_1(w[2], w[6], pars->GC_57, pars->mdl_MB, pars->ZERO, w[9]); 
  FFV2_1(w[2], w[6], pars->GC_92, pars->mdl_MB, pars->ZERO, w[10]); 
  FFS2_3(w[3], w[2], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[11]); 
  FFS1_3(w[3], w[2], pars->GC_128, pars->mdl_MH, pars->mdl_WH, w[12]); 
  FFV5P0_3(w[3], w[2], pars->GC_2, pars->ZERO, pars->ZERO, w[13]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_170, pars->GC_91, pars->ZERO, pars->ZERO,
      w[14]);
  FFV3_8_3(w[3], w[2], pars->GC_166, pars->GC_81, pars->mdl_MZ, pars->mdl_WZ,
      w[15]);
  FFV2_7_3(w[3], w[2], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[16]);
  FFV2_3(w[3], w[2], pars->GC_138, pars->mdl_MZ, pars->mdl_WZ, w[17]); 
  FFV2_2(w[5], w[6], pars->GC_57, pars->ZERO, pars->ZERO, w[18]); 
  FFV2_1(w[0], w[6], pars->GC_57, pars->ZERO, pars->ZERO, w[19]); 
  FFV5P0_3(w[3], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[20]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[21]);
  FFFF4_5_4(w[5], w[2], w[3], pars->GC_15, pars->GC_16, pars->ZERO, pars->ZERO,
      w[22]);
  FFFF1_4(w[5], w[2], w[3], -pars->GC_12, pars->ZERO, pars->ZERO, w[23]); 
  FFFF1_6_4(w[5], w[2], w[3], -pars->GC_10, -pars->GC_26, pars->ZERO,
      pars->ZERO, w[24]);
  FFFF3_2(w[3], w[5], w[0], pars->GC_14, pars->mdl_MB, pars->ZERO, w[25]); 
  FFFF3_2(w[3], w[5], w[0], pars->GC_25, pars->mdl_MB, pars->ZERO, w[26]); 
  FFFF4_5_1(w[2], w[3], w[0], pars->GC_17, pars->GC_16, pars->ZERO, pars->ZERO,
      w[27]);
  FFFF3_1(w[2], w[3], w[0], pars->GC_13, pars->ZERO, pars->ZERO, w[28]); 
  FFFF3_7_1(w[2], w[3], w[0], pars->GC_11, pars->GC_27, pars->ZERO, pars->ZERO,
      w[29]);
  FFV2_2(w[4], w[7], pars->GC_57, pars->ZERO, pars->ZERO, w[30]); 
  FFV1_2(w[4], w[13], pars->GC_2, pars->ZERO, pars->ZERO, w[31]); 
  FFV1_2(w[4], w[14], pars->GC_2, pars->ZERO, pars->ZERO, w[32]); 
  FFV1_2(w[4], w[20], pars->GC_7, pars->ZERO, pars->ZERO, w[33]); 
  FFV1_2(w[4], w[21], pars->GC_7, pars->ZERO, pars->ZERO, w[34]); 
  FFV2_7_2(w[4], w[15], pars->GC_59, pars->GC_69, pars->ZERO, pars->ZERO,
      w[35]);
  FFV2_7_2(w[4], w[16], pars->GC_59, pars->GC_69, pars->ZERO, pars->ZERO,
      w[36]);
  FFV2_7_2(w[4], w[17], pars->GC_59, pars->GC_69, pars->ZERO, pars->ZERO,
      w[37]);
  FFFF4_5_1(w[2], w[3], w[1], pars->GC_15, pars->GC_16, pars->ZERO, pars->ZERO,
      w[38]);
  FFFF1_1(w[2], w[3], w[1], -pars->GC_12, pars->ZERO, pars->ZERO, w[39]); 
  FFFF1_2_6_1(w[2], w[3], w[1], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->ZERO, pars->ZERO, w[40]);
  FFFF4_5_2(w[3], w[4], w[2], pars->GC_16, pars->GC_17, pars->ZERO, pars->ZERO,
      w[41]);
  FFFF1_2(w[3], w[4], w[2], -pars->GC_13, pars->ZERO, pars->ZERO, w[42]); 
  FFFF1_6_2(w[3], w[4], w[2], -pars->GC_11, -pars->GC_27, pars->ZERO,
      pars->ZERO, w[43]);
  FFFF3_3(w[4], w[1], w[2], pars->GC_14, pars->mdl_MB, pars->ZERO, w[44]); 
  FFFF3_3(w[4], w[1], w[2], pars->GC_25, pars->mdl_MB, pars->ZERO, w[45]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFVV1_2_0(w[3], w[2], w[6], w[7], pars->GC_168, pars->GC_89, amp[0]); 
  FFV2_0(w[3], w[8], w[7], pars->GC_57, amp[1]); 
  FFV8_0(w[3], w[9], w[7], pars->GC_80, amp[2]); 
  FFV2_0(w[3], w[9], w[7], pars->GC_57, amp[3]); 
  FFV2_0(w[3], w[9], w[7], pars->GC_108, amp[4]); 
  FFV2_0(w[3], w[10], w[7], pars->GC_57, amp[5]); 
  VVS1_0(w[6], w[7], w[11], pars->GC_87, amp[6]); 
  VVS1_0(w[6], w[7], w[12], pars->GC_87, amp[7]); 
  VVV1_0(w[13], w[6], w[7], pars->GC_4, amp[8]); 
  VVV1_0(w[14], w[6], w[7], pars->GC_4, amp[9]); 
  VVV1_0(w[6], w[7], w[15], pars->GC_60, amp[10]); 
  VVV1_0(w[6], w[7], w[16], pars->GC_60, amp[11]); 
  VVV1_0(w[6], w[7], w[17], pars->GC_60, amp[12]); 
  FFV1_0(w[18], w[0], w[13], pars->GC_2, amp[13]); 
  FFV1_0(w[5], w[19], w[13], pars->GC_1, amp[14]); 
  FFV1_0(w[18], w[0], w[14], pars->GC_2, amp[15]); 
  FFV1_0(w[5], w[19], w[14], pars->GC_1, amp[16]); 
  FFV1_0(w[18], w[0], w[20], pars->GC_7, amp[17]); 
  FFV1_0(w[5], w[19], w[20], pars->GC_7, amp[18]); 
  FFV1_0(w[18], w[0], w[21], pars->GC_7, amp[19]); 
  FFV1_0(w[5], w[19], w[21], pars->GC_7, amp[20]); 
  FFV2_7_0(w[18], w[0], w[15], pars->GC_59, pars->GC_69, amp[21]); 
  FFV2_4_0(w[5], w[19], w[15], pars->GC_58, pars->GC_69, amp[22]); 
  FFV2_7_0(w[18], w[0], w[16], pars->GC_59, pars->GC_69, amp[23]); 
  FFV2_4_0(w[5], w[19], w[16], pars->GC_58, pars->GC_69, amp[24]); 
  FFV2_7_0(w[18], w[0], w[17], pars->GC_59, pars->GC_69, amp[25]); 
  FFV2_4_0(w[5], w[19], w[17], pars->GC_58, pars->GC_69, amp[26]); 
  FFV2_0(w[22], w[0], w[6], pars->GC_57, amp[27]); 
  FFV2_0(w[23], w[0], w[6], pars->GC_57, amp[28]); 
  FFV2_0(w[24], w[0], w[6], pars->GC_57, amp[29]); 
  FFV2_0(w[25], w[2], w[6], pars->GC_57, amp[30]); 
  FFV2_0(w[26], w[2], w[6], pars->GC_57, amp[31]); 
  FFV2_0(w[5], w[27], w[6], pars->GC_57, amp[32]); 
  FFV2_0(w[5], w[28], w[6], pars->GC_57, amp[33]); 
  FFV2_0(w[5], w[29], w[6], pars->GC_57, amp[34]); 
  FFV1_0(w[30], w[1], w[13], pars->GC_1, amp[35]); 
  FFV2_0(w[31], w[1], w[7], pars->GC_57, amp[36]); 
  FFV1_0(w[30], w[1], w[14], pars->GC_1, amp[37]); 
  FFV2_0(w[32], w[1], w[7], pars->GC_57, amp[38]); 
  FFV1_0(w[30], w[1], w[20], pars->GC_7, amp[39]); 
  FFV2_0(w[33], w[1], w[7], pars->GC_57, amp[40]); 
  FFV1_0(w[30], w[1], w[21], pars->GC_7, amp[41]); 
  FFV2_0(w[34], w[1], w[7], pars->GC_57, amp[42]); 
  FFV2_4_0(w[30], w[1], w[15], pars->GC_58, pars->GC_69, amp[43]); 
  FFV2_0(w[35], w[1], w[7], pars->GC_57, amp[44]); 
  FFV2_4_0(w[30], w[1], w[16], pars->GC_58, pars->GC_69, amp[45]); 
  FFV2_0(w[36], w[1], w[7], pars->GC_57, amp[46]); 
  FFV2_4_0(w[30], w[1], w[17], pars->GC_58, pars->GC_69, amp[47]); 
  FFV2_0(w[37], w[1], w[7], pars->GC_57, amp[48]); 
  FFV2_0(w[4], w[38], w[7], pars->GC_57, amp[49]); 
  FFV2_0(w[4], w[39], w[7], pars->GC_57, amp[50]); 
  FFV2_0(w[4], w[40], w[7], pars->GC_57, amp[51]); 
  FFV2_0(w[41], w[1], w[7], pars->GC_57, amp[52]); 
  FFV2_0(w[42], w[1], w[7], pars->GC_57, amp[53]); 
  FFV2_0(w[43], w[1], w[7], pars->GC_57, amp[54]); 
  FFV2_0(w[3], w[44], w[7], pars->GC_57, amp[55]); 
  FFV2_0(w[3], w[45], w[7], pars->GC_57, amp[56]); 

}
double CPPProcess::matrix_3_cxdx_ttxuxsx() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 57; 
  const int ncolor = 3; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1, 1, 1}; 
  static const double cf[ncolor][ncolor] = {{27, 3, 9}, {3, 27, 9}, {9, 9,
      27}};

  // Calculate color flows
  jamp[0] = -1./2. * amp[17] - 1./2. * amp[18] - 1./2. * amp[19] - 1./2. *
      amp[20] + amp[27] + 1./2. * amp[29] - 1./2. * amp[31] + amp[32] + 1./2. *
      amp[34];
  jamp[1] = -1./2. * amp[39] - 1./2. * amp[40] - 1./2. * amp[41] - 1./2. *
      amp[42] + amp[49] + 1./2. * amp[51] + amp[52] + 1./2. * amp[54] - 1./2. *
      amp[56];
  jamp[2] = -amp[0] - amp[1] - amp[2] - amp[3] - amp[4] - amp[5] - amp[6] -
      amp[7] - amp[8] - amp[9] - amp[10] - amp[11] - amp[12] - amp[13] -
      amp[14] - amp[15] - amp[16] + 1./6. * amp[17] + 1./6. * amp[18] + 1./6. *
      amp[19] + 1./6. * amp[20] - amp[21] - amp[22] - amp[23] - amp[24] -
      amp[25] - amp[26] + amp[28] - 1./6. * amp[29] + 1./6. * amp[31] - amp[30]
      + amp[33] - 1./6. * amp[34] - amp[35] - amp[36] - amp[37] - amp[38] +
      1./6. * amp[39] + 1./6. * amp[40] + 1./6. * amp[41] + 1./6. * amp[42] -
      amp[43] - amp[44] - amp[45] - amp[46] - amp[47] - amp[48] + amp[50] -
      1./6. * amp[51] + amp[53] - 1./6. * amp[54] + 1./6. * amp[56] - amp[55];

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



