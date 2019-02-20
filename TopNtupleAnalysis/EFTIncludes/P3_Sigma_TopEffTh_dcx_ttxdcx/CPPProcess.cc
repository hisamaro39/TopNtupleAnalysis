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
// Process: d c~ > t t~ d c~ NP<=2 @3

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
  jamp2[0] = new double[6]; 
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
  for(int i = 0; i < 6; i++ )
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
        t[0] = matrix_3_dcx_ttxdcx(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_3_dcx_ttxdcx(); 
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
      t[0] = matrix_3_dcx_ttxdcx(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_3_dcx_ttxdcx(); 
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
  if(id1 == -4 && id2 == 1)
  {
    // Add matrix elements for processes with beams (-4, 1)
    return matrix_element[1]; 
  }
  else if(id1 == 1 && id2 == -4)
  {
    // Add matrix elements for processes with beams (1, -4)
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
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]); 
  FFV1P0_3(w[5], w[1], pars->GC_2, pars->ZERO, pars->ZERO, w[6]); 
  FFV1P0_3(w[0], w[4], pars->GC_1, pars->ZERO, pars->ZERO, w[7]); 
  FFV5_1(w[2], w[6], pars->GC_2, pars->mdl_MT, pars->mdl_WT, w[8]); 
  FFV3_8_1(w[2], w[6], pars->GC_170, pars->GC_91, pars->mdl_MT, pars->mdl_WT,
      w[9]);
  FFV5_2(w[3], w[6], pars->GC_2, pars->mdl_MT, pars->mdl_WT, w[10]); 
  FFV3_8_2(w[3], w[6], pars->GC_170, pars->GC_91, pars->mdl_MT, pars->mdl_WT,
      w[11]);
  FFV1P0_3(w[0], w[4], pars->GC_7, pars->ZERO, pars->ZERO, w[12]); 
  FFV2_4_3(w[0], w[4], pars->GC_58, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[13]);
  FFV1P0_3(w[5], w[1], pars->GC_7, pars->ZERO, pars->ZERO, w[14]); 
  FFV5_1(w[2], w[14], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[15]); 
  FFV3_8_1(w[2], w[14], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[16]);
  FFV5_2(w[3], w[14], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[17]); 
  FFV3_8_2(w[3], w[14], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[18]);
  FFV2_7_3(w[5], w[1], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[19]);
  FFV3_8_1(w[2], w[19], pars->GC_166, pars->GC_81, pars->mdl_MT, pars->mdl_WT,
      w[20]);
  FFV2_7_1(w[2], w[19], pars->GC_59, pars->GC_69, pars->mdl_MT, pars->mdl_WT,
      w[21]);
  FFV2_1(w[2], w[19], pars->GC_138, pars->mdl_MT, pars->mdl_WT, w[22]); 
  FFV3_8_2(w[3], w[19], pars->GC_166, pars->GC_81, pars->mdl_MT, pars->mdl_WT,
      w[23]);
  FFV2_7_2(w[3], w[19], pars->GC_59, pars->GC_69, pars->mdl_MT, pars->mdl_WT,
      w[24]);
  FFV2_2(w[3], w[19], pars->GC_138, pars->mdl_MT, pars->mdl_WT, w[25]); 
  FFS2_3(w[3], w[2], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[26]); 
  FFV5P0_3(w[3], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[27]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[28]);
  FFS1_3(w[3], w[2], pars->GC_128, pars->mdl_MH, pars->mdl_WH, w[29]); 
  FFV5P0_3(w[3], w[2], pars->GC_2, pars->ZERO, pars->ZERO, w[30]); 
  FFV1_2(w[0], w[6], pars->GC_1, pars->ZERO, pars->ZERO, w[31]); 
  FFV1_1(w[4], w[6], pars->GC_1, pars->ZERO, pars->ZERO, w[32]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_170, pars->GC_91, pars->ZERO, pars->ZERO,
      w[33]);
  FFV3_8_3(w[3], w[2], pars->GC_166, pars->GC_81, pars->mdl_MZ, pars->mdl_WZ,
      w[34]);
  FFV2_7_3(w[3], w[2], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[35]);
  FFV2_3(w[3], w[2], pars->GC_138, pars->mdl_MZ, pars->mdl_WZ, w[36]); 
  FFV1_2(w[0], w[14], pars->GC_7, pars->ZERO, pars->ZERO, w[37]); 
  FFV1_1(w[4], w[14], pars->GC_7, pars->ZERO, pars->ZERO, w[38]); 
  FFV2_4_2(w[0], w[19], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[39]);
  FFV2_4_1(w[4], w[19], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[40]);
  FFFF4_5_4(w[0], w[2], w[3], pars->GC_15, pars->GC_16, pars->ZERO, pars->ZERO,
      w[41]);
  FFFF1_4(w[0], w[2], w[3], -pars->GC_12, pars->ZERO, pars->ZERO, w[42]); 
  FFFF1_2_6_4(w[0], w[2], w[3], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->ZERO, pars->ZERO, w[43]);
  FFFF4_5_3(w[0], w[2], w[4], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[44]);
  FFFF1_3(w[0], w[2], w[4], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[45]); 
  FFFF1_2_6_3(w[0], w[2], w[4], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[46]);
  FFFF4_5_2(w[0], w[3], w[4], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[47]);
  FFFF1_2(w[0], w[3], w[4], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[48]); 
  FFFF1_2_6_2(w[0], w[3], w[4], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[49]);
  FFFF4_5_1(w[2], w[3], w[4], pars->GC_15, pars->GC_16, pars->ZERO, pars->ZERO,
      w[50]);
  FFFF1_1(w[2], w[3], w[4], -pars->GC_12, pars->ZERO, pars->ZERO, w[51]); 
  FFFF1_2_6_1(w[2], w[3], w[4], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->ZERO, pars->ZERO, w[52]);
  FFV1_2(w[5], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[53]); 
  FFV1_2(w[5], w[30], pars->GC_2, pars->ZERO, pars->ZERO, w[54]); 
  FFV1_2(w[5], w[33], pars->GC_2, pars->ZERO, pars->ZERO, w[55]); 
  FFV1_2(w[5], w[27], pars->GC_7, pars->ZERO, pars->ZERO, w[56]); 
  FFV1_2(w[5], w[28], pars->GC_7, pars->ZERO, pars->ZERO, w[57]); 
  FFV2_7_2(w[5], w[34], pars->GC_59, pars->GC_69, pars->ZERO, pars->ZERO,
      w[58]);
  FFV2_7_2(w[5], w[35], pars->GC_59, pars->GC_69, pars->ZERO, pars->ZERO,
      w[59]);
  FFV2_7_2(w[5], w[36], pars->GC_59, pars->GC_69, pars->ZERO, pars->ZERO,
      w[60]);
  FFV1_2(w[5], w[12], pars->GC_7, pars->ZERO, pars->ZERO, w[61]); 
  FFV2_7_2(w[5], w[13], pars->GC_59, pars->GC_69, pars->ZERO, pars->ZERO,
      w[62]);
  FFFF4_5_1(w[2], w[3], w[1], pars->GC_17, pars->GC_16, pars->ZERO, pars->ZERO,
      w[63]);
  FFFF3_1(w[2], w[3], w[1], pars->GC_13, pars->ZERO, pars->ZERO, w[64]); 
  FFFF3_7_1(w[2], w[3], w[1], pars->GC_11, pars->GC_27, pars->ZERO, pars->ZERO,
      w[65]);
  FFFF4_5_4(w[5], w[2], w[3], pars->GC_17, pars->GC_16, pars->ZERO, pars->ZERO,
      w[66]);
  FFFF3_4(w[5], w[2], w[3], pars->GC_13, pars->ZERO, pars->ZERO, w[67]); 
  FFFF3_7_4(w[5], w[2], w[3], pars->GC_11, pars->GC_27, pars->ZERO, pars->ZERO,
      w[68]);
  FFFF4_5_3(w[5], w[2], w[1], pars->GC_17, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[69]);
  FFFF3_3(w[5], w[2], w[1], pars->GC_13, pars->mdl_MT, pars->mdl_WT, w[70]); 
  FFFF3_7_3(w[5], w[2], w[1], pars->GC_11, pars->GC_27, pars->mdl_MT,
      pars->mdl_WT, w[71]);
  FFFF4_5_2(w[5], w[3], w[1], pars->GC_17, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[72]);
  FFFF3_2(w[5], w[3], w[1], pars->GC_13, pars->mdl_MT, pars->mdl_WT, w[73]); 
  FFFF3_7_2(w[5], w[3], w[1], pars->GC_11, pars->GC_27, pars->mdl_MT,
      pars->mdl_WT, w[74]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV5_0(w[3], w[8], w[7], pars->GC_2, amp[0]); 
  FFV3_8_0(w[3], w[8], w[7], pars->GC_170, pars->GC_91, amp[1]); 
  FFV5_0(w[3], w[9], w[7], pars->GC_2, amp[2]); 
  FFV5_0(w[10], w[2], w[7], pars->GC_2, amp[3]); 
  FFV3_8_0(w[10], w[2], w[7], pars->GC_170, pars->GC_91, amp[4]); 
  FFV5_0(w[11], w[2], w[7], pars->GC_2, amp[5]); 
  FFV5_0(w[3], w[8], w[12], pars->GC_7, amp[6]); 
  FFV3_8_0(w[3], w[8], w[12], pars->GC_145, pars->GC_79, amp[7]); 
  FFV5_0(w[3], w[9], w[12], pars->GC_7, amp[8]); 
  FFV5_0(w[10], w[2], w[12], pars->GC_7, amp[9]); 
  FFV3_8_0(w[10], w[2], w[12], pars->GC_145, pars->GC_79, amp[10]); 
  FFV5_0(w[11], w[2], w[12], pars->GC_7, amp[11]); 
  FFV3_8_0(w[3], w[8], w[13], pars->GC_166, pars->GC_81, amp[12]); 
  FFV2_7_0(w[3], w[8], w[13], pars->GC_59, pars->GC_69, amp[13]); 
  FFV2_0(w[3], w[8], w[13], pars->GC_138, amp[14]); 
  FFV2_7_0(w[3], w[9], w[13], pars->GC_59, pars->GC_69, amp[15]); 
  FFV3_8_0(w[10], w[2], w[13], pars->GC_166, pars->GC_81, amp[16]); 
  FFV2_7_0(w[10], w[2], w[13], pars->GC_59, pars->GC_69, amp[17]); 
  FFV2_0(w[10], w[2], w[13], pars->GC_138, amp[18]); 
  FFV2_7_0(w[11], w[2], w[13], pars->GC_59, pars->GC_69, amp[19]); 
  FFV5_0(w[3], w[15], w[7], pars->GC_2, amp[20]); 
  FFV3_8_0(w[3], w[15], w[7], pars->GC_170, pars->GC_91, amp[21]); 
  FFV5_0(w[3], w[16], w[7], pars->GC_2, amp[22]); 
  FFV5_0(w[17], w[2], w[7], pars->GC_2, amp[23]); 
  FFV3_8_0(w[17], w[2], w[7], pars->GC_170, pars->GC_91, amp[24]); 
  FFV5_0(w[18], w[2], w[7], pars->GC_2, amp[25]); 
  FFVV1_2_0(w[3], w[2], w[14], w[12], pars->GC_146, pars->GC_85, amp[26]); 
  FFV5_0(w[3], w[15], w[12], pars->GC_7, amp[27]); 
  FFV3_8_0(w[3], w[15], w[12], pars->GC_145, pars->GC_79, amp[28]); 
  FFV5_0(w[3], w[16], w[12], pars->GC_7, amp[29]); 
  FFV5_0(w[17], w[2], w[12], pars->GC_7, amp[30]); 
  FFV3_8_0(w[17], w[2], w[12], pars->GC_145, pars->GC_79, amp[31]); 
  FFV5_0(w[18], w[2], w[12], pars->GC_7, amp[32]); 
  FFV3_8_0(w[3], w[15], w[13], pars->GC_166, pars->GC_81, amp[33]); 
  FFV2_7_0(w[3], w[15], w[13], pars->GC_59, pars->GC_69, amp[34]); 
  FFV2_0(w[3], w[15], w[13], pars->GC_138, amp[35]); 
  FFV2_7_0(w[3], w[16], w[13], pars->GC_59, pars->GC_69, amp[36]); 
  FFV3_8_0(w[17], w[2], w[13], pars->GC_166, pars->GC_81, amp[37]); 
  FFV2_7_0(w[17], w[2], w[13], pars->GC_59, pars->GC_69, amp[38]); 
  FFV2_0(w[17], w[2], w[13], pars->GC_138, amp[39]); 
  FFV2_7_0(w[18], w[2], w[13], pars->GC_59, pars->GC_69, amp[40]); 
  FFV5_0(w[3], w[20], w[7], pars->GC_2, amp[41]); 
  FFV5_0(w[3], w[21], w[7], pars->GC_2, amp[42]); 
  FFV3_8_0(w[3], w[21], w[7], pars->GC_170, pars->GC_91, amp[43]); 
  FFV5_0(w[3], w[22], w[7], pars->GC_2, amp[44]); 
  FFV5_0(w[23], w[2], w[7], pars->GC_2, amp[45]); 
  FFV5_0(w[24], w[2], w[7], pars->GC_2, amp[46]); 
  FFV3_8_0(w[24], w[2], w[7], pars->GC_170, pars->GC_91, amp[47]); 
  FFV5_0(w[25], w[2], w[7], pars->GC_2, amp[48]); 
  FFV5_0(w[3], w[20], w[12], pars->GC_7, amp[49]); 
  FFV5_0(w[3], w[21], w[12], pars->GC_7, amp[50]); 
  FFV3_8_0(w[3], w[21], w[12], pars->GC_145, pars->GC_79, amp[51]); 
  FFV5_0(w[3], w[22], w[12], pars->GC_7, amp[52]); 
  FFV5_0(w[23], w[2], w[12], pars->GC_7, amp[53]); 
  FFV5_0(w[24], w[2], w[12], pars->GC_7, amp[54]); 
  FFV3_8_0(w[24], w[2], w[12], pars->GC_145, pars->GC_79, amp[55]); 
  FFV5_0(w[25], w[2], w[12], pars->GC_7, amp[56]); 
  FFV2_7_0(w[3], w[20], w[13], pars->GC_59, pars->GC_69, amp[57]); 
  FFV3_8_0(w[3], w[21], w[13], pars->GC_166, pars->GC_81, amp[58]); 
  FFV2_7_0(w[3], w[21], w[13], pars->GC_59, pars->GC_69, amp[59]); 
  FFV2_0(w[3], w[21], w[13], pars->GC_138, amp[60]); 
  FFV2_7_0(w[3], w[22], w[13], pars->GC_59, pars->GC_69, amp[61]); 
  FFV2_7_0(w[23], w[2], w[13], pars->GC_59, pars->GC_69, amp[62]); 
  FFV3_8_0(w[24], w[2], w[13], pars->GC_166, pars->GC_81, amp[63]); 
  FFV2_7_0(w[24], w[2], w[13], pars->GC_59, pars->GC_69, amp[64]); 
  FFV2_0(w[24], w[2], w[13], pars->GC_138, amp[65]); 
  FFV2_7_0(w[25], w[2], w[13], pars->GC_59, pars->GC_69, amp[66]); 
  VVS2_0(w[14], w[12], w[26], pars->GC_78, amp[67]); 
  VVV2_0(w[14], w[12], w[27], pars->GC_28, amp[68]); 
  VVV1_0(w[14], w[12], w[27], pars->GC_6, amp[69]); 
  VVV1_0(w[14], w[12], w[28], pars->GC_6, amp[70]); 
  VVS1_0(w[19], w[13], w[26], pars->GC_93, amp[71]); 
  VVS1_0(w[19], w[13], w[29], pars->GC_93, amp[72]); 
  FFV1_0(w[31], w[4], w[30], pars->GC_1, amp[73]); 
  FFV1_0(w[0], w[32], w[30], pars->GC_1, amp[74]); 
  FFV1_0(w[31], w[4], w[33], pars->GC_1, amp[75]); 
  FFV1_0(w[0], w[32], w[33], pars->GC_1, amp[76]); 
  FFV1_0(w[31], w[4], w[27], pars->GC_7, amp[77]); 
  FFV1_0(w[0], w[32], w[27], pars->GC_7, amp[78]); 
  FFV1_0(w[31], w[4], w[28], pars->GC_7, amp[79]); 
  FFV1_0(w[0], w[32], w[28], pars->GC_7, amp[80]); 
  FFV2_4_0(w[31], w[4], w[34], pars->GC_58, pars->GC_69, amp[81]); 
  FFV2_4_0(w[0], w[32], w[34], pars->GC_58, pars->GC_69, amp[82]); 
  FFV2_4_0(w[31], w[4], w[35], pars->GC_58, pars->GC_69, amp[83]); 
  FFV2_4_0(w[0], w[32], w[35], pars->GC_58, pars->GC_69, amp[84]); 
  FFV2_4_0(w[31], w[4], w[36], pars->GC_58, pars->GC_69, amp[85]); 
  FFV2_4_0(w[0], w[32], w[36], pars->GC_58, pars->GC_69, amp[86]); 
  FFV1_0(w[37], w[4], w[30], pars->GC_1, amp[87]); 
  FFV1_0(w[0], w[38], w[30], pars->GC_1, amp[88]); 
  FFV1_0(w[37], w[4], w[33], pars->GC_1, amp[89]); 
  FFV1_0(w[0], w[38], w[33], pars->GC_1, amp[90]); 
  FFV1_0(w[37], w[4], w[27], pars->GC_7, amp[91]); 
  FFV1_0(w[0], w[38], w[27], pars->GC_7, amp[92]); 
  FFV1_0(w[37], w[4], w[28], pars->GC_7, amp[93]); 
  FFV1_0(w[0], w[38], w[28], pars->GC_7, amp[94]); 
  FFV2_4_0(w[37], w[4], w[34], pars->GC_58, pars->GC_69, amp[95]); 
  FFV2_4_0(w[0], w[38], w[34], pars->GC_58, pars->GC_69, amp[96]); 
  FFV2_4_0(w[37], w[4], w[35], pars->GC_58, pars->GC_69, amp[97]); 
  FFV2_4_0(w[0], w[38], w[35], pars->GC_58, pars->GC_69, amp[98]); 
  FFV2_4_0(w[37], w[4], w[36], pars->GC_58, pars->GC_69, amp[99]); 
  FFV2_4_0(w[0], w[38], w[36], pars->GC_58, pars->GC_69, amp[100]); 
  FFV1_0(w[39], w[4], w[30], pars->GC_1, amp[101]); 
  FFV1_0(w[0], w[40], w[30], pars->GC_1, amp[102]); 
  FFV1_0(w[39], w[4], w[33], pars->GC_1, amp[103]); 
  FFV1_0(w[0], w[40], w[33], pars->GC_1, amp[104]); 
  FFV1_0(w[39], w[4], w[27], pars->GC_7, amp[105]); 
  FFV1_0(w[0], w[40], w[27], pars->GC_7, amp[106]); 
  FFV1_0(w[39], w[4], w[28], pars->GC_7, amp[107]); 
  FFV1_0(w[0], w[40], w[28], pars->GC_7, amp[108]); 
  FFV2_4_0(w[39], w[4], w[34], pars->GC_58, pars->GC_69, amp[109]); 
  FFV2_4_0(w[0], w[40], w[34], pars->GC_58, pars->GC_69, amp[110]); 
  FFV2_4_0(w[39], w[4], w[35], pars->GC_58, pars->GC_69, amp[111]); 
  FFV2_4_0(w[0], w[40], w[35], pars->GC_58, pars->GC_69, amp[112]); 
  FFV2_4_0(w[39], w[4], w[36], pars->GC_58, pars->GC_69, amp[113]); 
  FFV2_4_0(w[0], w[40], w[36], pars->GC_58, pars->GC_69, amp[114]); 
  FFV1_0(w[41], w[4], w[6], pars->GC_1, amp[115]); 
  FFV1_0(w[42], w[4], w[6], pars->GC_1, amp[116]); 
  FFV1_0(w[43], w[4], w[6], pars->GC_1, amp[117]); 
  FFV1_0(w[41], w[4], w[14], pars->GC_7, amp[118]); 
  FFV1_0(w[42], w[4], w[14], pars->GC_7, amp[119]); 
  FFV1_0(w[43], w[4], w[14], pars->GC_7, amp[120]); 
  FFV2_4_0(w[41], w[4], w[19], pars->GC_58, pars->GC_69, amp[121]); 
  FFV2_4_0(w[42], w[4], w[19], pars->GC_58, pars->GC_69, amp[122]); 
  FFV2_4_0(w[43], w[4], w[19], pars->GC_58, pars->GC_69, amp[123]); 
  FFV5_0(w[3], w[44], w[6], pars->GC_2, amp[124]); 
  FFV5_0(w[3], w[45], w[6], pars->GC_2, amp[125]); 
  FFV5_0(w[3], w[46], w[6], pars->GC_2, amp[126]); 
  FFV5_0(w[3], w[44], w[14], pars->GC_7, amp[127]); 
  FFV5_0(w[3], w[45], w[14], pars->GC_7, amp[128]); 
  FFV5_0(w[3], w[46], w[14], pars->GC_7, amp[129]); 
  FFV2_7_0(w[3], w[44], w[19], pars->GC_59, pars->GC_69, amp[130]); 
  FFV2_7_0(w[3], w[45], w[19], pars->GC_59, pars->GC_69, amp[131]); 
  FFV2_7_0(w[3], w[46], w[19], pars->GC_59, pars->GC_69, amp[132]); 
  FFV5_0(w[47], w[2], w[6], pars->GC_2, amp[133]); 
  FFV5_0(w[48], w[2], w[6], pars->GC_2, amp[134]); 
  FFV5_0(w[49], w[2], w[6], pars->GC_2, amp[135]); 
  FFV5_0(w[47], w[2], w[14], pars->GC_7, amp[136]); 
  FFV5_0(w[48], w[2], w[14], pars->GC_7, amp[137]); 
  FFV5_0(w[49], w[2], w[14], pars->GC_7, amp[138]); 
  FFV2_7_0(w[47], w[2], w[19], pars->GC_59, pars->GC_69, amp[139]); 
  FFV2_7_0(w[48], w[2], w[19], pars->GC_59, pars->GC_69, amp[140]); 
  FFV2_7_0(w[49], w[2], w[19], pars->GC_59, pars->GC_69, amp[141]); 
  FFV1_0(w[0], w[50], w[6], pars->GC_1, amp[142]); 
  FFV1_0(w[0], w[51], w[6], pars->GC_1, amp[143]); 
  FFV1_0(w[0], w[52], w[6], pars->GC_1, amp[144]); 
  FFV1_0(w[0], w[50], w[14], pars->GC_7, amp[145]); 
  FFV1_0(w[0], w[51], w[14], pars->GC_7, amp[146]); 
  FFV1_0(w[0], w[52], w[14], pars->GC_7, amp[147]); 
  FFV2_4_0(w[0], w[50], w[19], pars->GC_58, pars->GC_69, amp[148]); 
  FFV2_4_0(w[0], w[51], w[19], pars->GC_58, pars->GC_69, amp[149]); 
  FFV2_4_0(w[0], w[52], w[19], pars->GC_58, pars->GC_69, amp[150]); 
  FFV1_0(w[53], w[1], w[30], pars->GC_2, amp[151]); 
  FFV1_0(w[54], w[1], w[7], pars->GC_2, amp[152]); 
  FFV1_0(w[53], w[1], w[33], pars->GC_2, amp[153]); 
  FFV1_0(w[55], w[1], w[7], pars->GC_2, amp[154]); 
  FFV1_0(w[53], w[1], w[27], pars->GC_7, amp[155]); 
  FFV1_0(w[56], w[1], w[7], pars->GC_2, amp[156]); 
  FFV1_0(w[53], w[1], w[28], pars->GC_7, amp[157]); 
  FFV1_0(w[57], w[1], w[7], pars->GC_2, amp[158]); 
  FFV2_7_0(w[53], w[1], w[34], pars->GC_59, pars->GC_69, amp[159]); 
  FFV1_0(w[58], w[1], w[7], pars->GC_2, amp[160]); 
  FFV2_7_0(w[53], w[1], w[35], pars->GC_59, pars->GC_69, amp[161]); 
  FFV1_0(w[59], w[1], w[7], pars->GC_2, amp[162]); 
  FFV2_7_0(w[53], w[1], w[36], pars->GC_59, pars->GC_69, amp[163]); 
  FFV1_0(w[60], w[1], w[7], pars->GC_2, amp[164]); 
  FFV1_0(w[61], w[1], w[30], pars->GC_2, amp[165]); 
  FFV1_0(w[54], w[1], w[12], pars->GC_7, amp[166]); 
  FFV1_0(w[61], w[1], w[33], pars->GC_2, amp[167]); 
  FFV1_0(w[55], w[1], w[12], pars->GC_7, amp[168]); 
  FFV1_0(w[61], w[1], w[27], pars->GC_7, amp[169]); 
  FFV1_0(w[56], w[1], w[12], pars->GC_7, amp[170]); 
  FFV1_0(w[61], w[1], w[28], pars->GC_7, amp[171]); 
  FFV1_0(w[57], w[1], w[12], pars->GC_7, amp[172]); 
  FFV2_7_0(w[61], w[1], w[34], pars->GC_59, pars->GC_69, amp[173]); 
  FFV1_0(w[58], w[1], w[12], pars->GC_7, amp[174]); 
  FFV2_7_0(w[61], w[1], w[35], pars->GC_59, pars->GC_69, amp[175]); 
  FFV1_0(w[59], w[1], w[12], pars->GC_7, amp[176]); 
  FFV2_7_0(w[61], w[1], w[36], pars->GC_59, pars->GC_69, amp[177]); 
  FFV1_0(w[60], w[1], w[12], pars->GC_7, amp[178]); 
  FFV1_0(w[62], w[1], w[30], pars->GC_2, amp[179]); 
  FFV2_7_0(w[54], w[1], w[13], pars->GC_59, pars->GC_69, amp[180]); 
  FFV1_0(w[62], w[1], w[33], pars->GC_2, amp[181]); 
  FFV2_7_0(w[55], w[1], w[13], pars->GC_59, pars->GC_69, amp[182]); 
  FFV1_0(w[62], w[1], w[27], pars->GC_7, amp[183]); 
  FFV2_7_0(w[56], w[1], w[13], pars->GC_59, pars->GC_69, amp[184]); 
  FFV1_0(w[62], w[1], w[28], pars->GC_7, amp[185]); 
  FFV2_7_0(w[57], w[1], w[13], pars->GC_59, pars->GC_69, amp[186]); 
  FFV2_7_0(w[62], w[1], w[34], pars->GC_59, pars->GC_69, amp[187]); 
  FFV2_7_0(w[58], w[1], w[13], pars->GC_59, pars->GC_69, amp[188]); 
  FFV2_7_0(w[62], w[1], w[35], pars->GC_59, pars->GC_69, amp[189]); 
  FFV2_7_0(w[59], w[1], w[13], pars->GC_59, pars->GC_69, amp[190]); 
  FFV2_7_0(w[62], w[1], w[36], pars->GC_59, pars->GC_69, amp[191]); 
  FFV2_7_0(w[60], w[1], w[13], pars->GC_59, pars->GC_69, amp[192]); 
  FFV1_0(w[5], w[63], w[7], pars->GC_2, amp[193]); 
  FFV1_0(w[5], w[64], w[7], pars->GC_2, amp[194]); 
  FFV1_0(w[5], w[65], w[7], pars->GC_2, amp[195]); 
  FFV1_0(w[5], w[63], w[12], pars->GC_7, amp[196]); 
  FFV1_0(w[5], w[64], w[12], pars->GC_7, amp[197]); 
  FFV1_0(w[5], w[65], w[12], pars->GC_7, amp[198]); 
  FFV2_7_0(w[5], w[63], w[13], pars->GC_59, pars->GC_69, amp[199]); 
  FFV2_7_0(w[5], w[64], w[13], pars->GC_59, pars->GC_69, amp[200]); 
  FFV2_7_0(w[5], w[65], w[13], pars->GC_59, pars->GC_69, amp[201]); 
  FFV1_0(w[66], w[1], w[7], pars->GC_2, amp[202]); 
  FFV1_0(w[67], w[1], w[7], pars->GC_2, amp[203]); 
  FFV1_0(w[68], w[1], w[7], pars->GC_2, amp[204]); 
  FFV1_0(w[66], w[1], w[12], pars->GC_7, amp[205]); 
  FFV1_0(w[67], w[1], w[12], pars->GC_7, amp[206]); 
  FFV1_0(w[68], w[1], w[12], pars->GC_7, amp[207]); 
  FFV2_7_0(w[66], w[1], w[13], pars->GC_59, pars->GC_69, amp[208]); 
  FFV2_7_0(w[67], w[1], w[13], pars->GC_59, pars->GC_69, amp[209]); 
  FFV2_7_0(w[68], w[1], w[13], pars->GC_59, pars->GC_69, amp[210]); 
  FFV5_0(w[3], w[69], w[7], pars->GC_2, amp[211]); 
  FFV5_0(w[3], w[70], w[7], pars->GC_2, amp[212]); 
  FFV5_0(w[3], w[71], w[7], pars->GC_2, amp[213]); 
  FFV5_0(w[3], w[69], w[12], pars->GC_7, amp[214]); 
  FFV5_0(w[3], w[70], w[12], pars->GC_7, amp[215]); 
  FFV5_0(w[3], w[71], w[12], pars->GC_7, amp[216]); 
  FFV2_7_0(w[3], w[69], w[13], pars->GC_59, pars->GC_69, amp[217]); 
  FFV2_7_0(w[3], w[70], w[13], pars->GC_59, pars->GC_69, amp[218]); 
  FFV2_7_0(w[3], w[71], w[13], pars->GC_59, pars->GC_69, amp[219]); 
  FFV5_0(w[72], w[2], w[7], pars->GC_2, amp[220]); 
  FFV5_0(w[73], w[2], w[7], pars->GC_2, amp[221]); 
  FFV5_0(w[74], w[2], w[7], pars->GC_2, amp[222]); 
  FFV5_0(w[72], w[2], w[12], pars->GC_7, amp[223]); 
  FFV5_0(w[73], w[2], w[12], pars->GC_7, amp[224]); 
  FFV5_0(w[74], w[2], w[12], pars->GC_7, amp[225]); 
  FFV2_7_0(w[72], w[2], w[13], pars->GC_59, pars->GC_69, amp[226]); 
  FFV2_7_0(w[73], w[2], w[13], pars->GC_59, pars->GC_69, amp[227]); 
  FFV2_7_0(w[74], w[2], w[13], pars->GC_59, pars->GC_69, amp[228]); 

}
double CPPProcess::matrix_3_dcx_ttxdcx() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 229; 
  const int ncolor = 6; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1, 1, 1, 1, 1, 1}; 
  static const double cf[ncolor][ncolor] = {{27, 9, 9, 3, 3, 9}, {9, 27, 3, 9,
      9, 3}, {9, 3, 27, 9, 9, 3}, {3, 9, 9, 27, 3, 9}, {3, 9, 9, 3, 27, 9}, {9,
      3, 3, 9, 9, 27}};

  // Calculate color flows
  jamp[0] = -1./2. * amp[67] - 1./2. * amp[87] - 1./2. * amp[88] - 1./2. *
      amp[89] - 1./2. * amp[90] + 1./12. * amp[91] + 1./12. * amp[92] + 1./12.
      * amp[93] + 1./12. * amp[94] - 1./2. * amp[95] - 1./2. * amp[96] - 1./2.
      * amp[97] - 1./2. * amp[98] - 1./2. * amp[99] - 1./2. * amp[100] + 1./2.
      * amp[119] - 1./12. * amp[120] + 1./2. * amp[146] - 1./12. * amp[147] -
      1./2. * amp[165] - 1./2. * amp[166] - 1./2. * amp[167] - 1./2. * amp[168]
      + 1./12. * amp[169] + 1./12. * amp[170] + 1./12. * amp[171] + 1./12. *
      amp[172] - 1./2. * amp[173] - 1./2. * amp[174] - 1./2. * amp[175] - 1./2.
      * amp[176] - 1./2. * amp[177] - 1./2. * amp[178] + 1./2. * amp[197] -
      1./12. * amp[198] + 1./2. * amp[206] - 1./12. * amp[207];
  jamp[1] = +1./4. * std::complex<double> (0, 1) * amp[26] - 1./4. * amp[27] -
      1./4. * amp[28] - 1./4. * amp[29] + 1./4. * std::complex<double> (0, 1) *
      amp[68] + 1./4. * std::complex<double> (0, 1) * amp[69] + 1./4. *
      std::complex<double> (0, 1) * amp[70] - 1./4. * amp[91] - 1./4. * amp[93]
      + 1./2. * amp[136] + 1./4. * amp[138] + 1./2. * amp[145] + 1./4. *
      amp[147] - 1./4. * amp[170] - 1./4. * amp[172] + 1./2. * amp[205] + 1./4.
      * amp[207] + 1./2. * amp[214] + 1./4. * amp[216];
  jamp[2] = -1./4. * std::complex<double> (0, 1) * amp[26] - 1./4. * amp[30] -
      1./4. * amp[31] - 1./4. * amp[32] - 1./4. * std::complex<double> (0, 1) *
      amp[68] - 1./4. * std::complex<double> (0, 1) * amp[69] - 1./4. *
      std::complex<double> (0, 1) * amp[70] - 1./4. * amp[92] - 1./4. * amp[94]
      + 1./2. * amp[118] + 1./4. * amp[120] + 1./2. * amp[127] + 1./4. *
      amp[129] - 1./4. * amp[169] - 1./4. * amp[171] + 1./2. * amp[196] + 1./4.
      * amp[198] + 1./2. * amp[223] + 1./4. * amp[225];
  jamp[3] = -1./2. * amp[20] - 1./2. * amp[21] - 1./2. * amp[22] - 1./2. *
      amp[23] - 1./2. * amp[24] - 1./2. * amp[25] + 1./12. * amp[27] + 1./12. *
      amp[28] + 1./12. * amp[29] + 1./12. * amp[30] + 1./12. * amp[31] + 1./12.
      * amp[32] - 1./2. * amp[33] - 1./2. * amp[34] - 1./2. * amp[35] - 1./2. *
      amp[36] - 1./2. * amp[37] - 1./2. * amp[38] - 1./2. * amp[39] - 1./2. *
      amp[40] + 1./2. * amp[128] - 1./12. * amp[129] + 1./2. * amp[137] -
      1./12. * amp[138] - 1./2. * amp[155] - 1./2. * amp[156] - 1./2. *
      amp[157] - 1./2. * amp[158] + 1./12. * amp[169] + 1./12. * amp[170] +
      1./12. * amp[171] + 1./12. * amp[172] - 1./2. * amp[183] - 1./2. *
      amp[184] - 1./2. * amp[185] - 1./2. * amp[186] + amp[193] + 1./2. *
      amp[195] - 1./6. * amp[196] - 1./12. * amp[198] + amp[199] + 1./2. *
      amp[201] + amp[202] + 1./2. * amp[204] - 1./6. * amp[205] - 1./12. *
      amp[207] + amp[208] + 1./2. * amp[210] + amp[211] + 1./2. * amp[213] -
      1./6. * amp[214] - 1./12. * amp[216] + amp[217] + 1./2. * amp[219] +
      amp[220] + 1./2. * amp[222] - 1./6. * amp[223] - 1./12. * amp[225] +
      amp[226] + 1./2. * amp[228];
  jamp[4] = -1./2. * amp[6] - 1./2. * amp[7] - 1./2. * amp[8] - 1./2. * amp[9]
      - 1./2. * amp[10] - 1./2. * amp[11] + 1./12. * amp[27] + 1./12. * amp[28]
      + 1./12. * amp[29] + 1./12. * amp[30] + 1./12. * amp[31] + 1./12. *
      amp[32] - 1./2. * amp[49] - 1./2. * amp[50] - 1./2. * amp[51] - 1./2. *
      amp[52] - 1./2. * amp[53] - 1./2. * amp[54] - 1./2. * amp[55] - 1./2. *
      amp[56] - 1./2. * amp[77] - 1./2. * amp[78] - 1./2. * amp[79] - 1./2. *
      amp[80] + 1./12. * amp[91] + 1./12. * amp[92] + 1./12. * amp[93] + 1./12.
      * amp[94] - 1./2. * amp[105] - 1./2. * amp[106] - 1./2. * amp[107] -
      1./2. * amp[108] + amp[115] + 1./2. * amp[117] - 1./6. * amp[118] -
      1./12. * amp[120] + amp[121] + 1./2. * amp[123] + amp[124] + 1./2. *
      amp[126] - 1./6. * amp[127] - 1./12. * amp[129] + amp[130] + 1./2. *
      amp[132] + amp[133] + 1./2. * amp[135] - 1./6. * amp[136] - 1./12. *
      amp[138] + amp[139] + 1./2. * amp[141] + amp[142] + 1./2. * amp[144] -
      1./6. * amp[145] - 1./12. * amp[147] + amp[148] + 1./2. * amp[150] +
      1./2. * amp[215] - 1./12. * amp[216] + 1./2. * amp[224] - 1./12. *
      amp[225];
  jamp[5] = -amp[0] - amp[1] - amp[2] - amp[3] - amp[4] - amp[5] + 1./6. *
      amp[6] + 1./6. * amp[7] + 1./6. * amp[8] + 1./6. * amp[9] + 1./6. *
      amp[10] + 1./6. * amp[11] - amp[12] - amp[13] - amp[14] - amp[15] -
      amp[16] - amp[17] - amp[18] - amp[19] + 1./6. * amp[20] + 1./6. * amp[21]
      + 1./6. * amp[22] + 1./6. * amp[23] + 1./6. * amp[24] + 1./6. * amp[25] -
      1./36. * amp[27] - 1./36. * amp[28] - 1./36. * amp[29] - 1./36. * amp[30]
      - 1./36. * amp[31] - 1./36. * amp[32] + 1./6. * amp[33] + 1./6. * amp[34]
      + 1./6. * amp[35] + 1./6. * amp[36] + 1./6. * amp[37] + 1./6. * amp[38] +
      1./6. * amp[39] + 1./6. * amp[40] - amp[41] - amp[42] - amp[43] - amp[44]
      - amp[45] - amp[46] - amp[47] - amp[48] + 1./6. * amp[49] + 1./6. *
      amp[50] + 1./6. * amp[51] + 1./6. * amp[52] + 1./6. * amp[53] + 1./6. *
      amp[54] + 1./6. * amp[55] + 1./6. * amp[56] - amp[57] - amp[58] - amp[59]
      - amp[60] - amp[61] - amp[62] - amp[63] - amp[64] - amp[65] - amp[66] +
      1./6. * amp[67] - amp[71] - amp[72] - amp[73] - amp[74] - amp[75] -
      amp[76] + 1./6. * amp[77] + 1./6. * amp[78] + 1./6. * amp[79] + 1./6. *
      amp[80] - amp[81] - amp[82] - amp[83] - amp[84] - amp[85] - amp[86] +
      1./6. * amp[87] + 1./6. * amp[88] + 1./6. * amp[89] + 1./6. * amp[90] -
      1./36. * amp[91] - 1./36. * amp[92] - 1./36. * amp[93] - 1./36. * amp[94]
      + 1./6. * amp[95] + 1./6. * amp[96] + 1./6. * amp[97] + 1./6. * amp[98] +
      1./6. * amp[99] + 1./6. * amp[100] - amp[101] - amp[102] - amp[103] -
      amp[104] + 1./6. * amp[105] + 1./6. * amp[106] + 1./6. * amp[107] + 1./6.
      * amp[108] - amp[109] - amp[110] - amp[111] - amp[112] - amp[113] -
      amp[114] + amp[116] - 1./6. * amp[117] - 1./6. * amp[119] + 1./36. *
      amp[120] + amp[122] - 1./6. * amp[123] + amp[125] - 1./6. * amp[126] -
      1./6. * amp[128] + 1./36. * amp[129] + amp[131] - 1./6. * amp[132] +
      amp[134] - 1./6. * amp[135] - 1./6. * amp[137] + 1./36. * amp[138] +
      amp[140] - 1./6. * amp[141] + amp[143] - 1./6. * amp[144] - 1./6. *
      amp[146] + 1./36. * amp[147] + amp[149] - 1./6. * amp[150] - amp[151] -
      amp[152] - amp[153] - amp[154] + 1./6. * amp[155] + 1./6. * amp[156] +
      1./6. * amp[157] + 1./6. * amp[158] - amp[159] - amp[160] - amp[161] -
      amp[162] - amp[163] - amp[164] + 1./6. * amp[165] + 1./6. * amp[166] +
      1./6. * amp[167] + 1./6. * amp[168] - 1./36. * amp[169] - 1./36. *
      amp[170] - 1./36. * amp[171] - 1./36. * amp[172] + 1./6. * amp[173] +
      1./6. * amp[174] + 1./6. * amp[175] + 1./6. * amp[176] + 1./6. * amp[177]
      + 1./6. * amp[178] - amp[179] - amp[180] - amp[181] - amp[182] + 1./6. *
      amp[183] + 1./6. * amp[184] + 1./6. * amp[185] + 1./6. * amp[186] -
      amp[187] - amp[188] - amp[189] - amp[190] - amp[191] - amp[192] +
      amp[194] - 1./6. * amp[195] - 1./6. * amp[197] + 1./36. * amp[198] +
      amp[200] - 1./6. * amp[201] + amp[203] - 1./6. * amp[204] - 1./6. *
      amp[206] + 1./36. * amp[207] + amp[209] - 1./6. * amp[210] + amp[212] -
      1./6. * amp[213] - 1./6. * amp[215] + 1./36. * amp[216] + amp[218] -
      1./6. * amp[219] + amp[221] - 1./6. * amp[222] - 1./6. * amp[224] +
      1./36. * amp[225] + amp[227] - 1./6. * amp[228];

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



