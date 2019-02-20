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
// Process: d~ d~ > t t~ d~ d~ NP<=2 @3

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
  const int denominators[nprocesses] = {72}; 

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
        t[0] = matrix_3_dxdx_ttxdxdx(); 

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
      t[0] = matrix_3_dxdx_ttxdxdx(); 

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
  if(id1 == -1 && id2 == -1)
  {
    // Add matrix elements for processes with beams (-1, -1)
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
  FFV1P0_3(w[4], w[0], pars->GC_1, pars->ZERO, pars->ZERO, w[6]); 
  FFV1P0_3(w[5], w[1], pars->GC_1, pars->ZERO, pars->ZERO, w[7]); 
  FFV5_1(w[2], w[6], pars->GC_2, pars->mdl_MT, pars->mdl_WT, w[8]); 
  FFV3_8_1(w[2], w[6], pars->GC_170, pars->GC_91, pars->mdl_MT, pars->mdl_WT,
      w[9]);
  FFV5_2(w[3], w[6], pars->GC_2, pars->mdl_MT, pars->mdl_WT, w[10]); 
  FFV3_8_2(w[3], w[6], pars->GC_170, pars->GC_91, pars->mdl_MT, pars->mdl_WT,
      w[11]);
  FFV1P0_3(w[5], w[1], pars->GC_7, pars->ZERO, pars->ZERO, w[12]); 
  FFV2_4_3(w[5], w[1], pars->GC_58, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[13]);
  FFV1P0_3(w[4], w[0], pars->GC_7, pars->ZERO, pars->ZERO, w[14]); 
  FFV5_1(w[2], w[14], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[15]); 
  FFV3_8_1(w[2], w[14], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[16]);
  FFV5_2(w[3], w[14], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[17]); 
  FFV3_8_2(w[3], w[14], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[18]);
  FFV2_4_3(w[4], w[0], pars->GC_58, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
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
  FFV1_2(w[5], w[6], pars->GC_1, pars->ZERO, pars->ZERO, w[31]); 
  FFV1_1(w[1], w[6], pars->GC_1, pars->ZERO, pars->ZERO, w[32]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_170, pars->GC_91, pars->ZERO, pars->ZERO,
      w[33]);
  FFV3_8_3(w[3], w[2], pars->GC_166, pars->GC_81, pars->mdl_MZ, pars->mdl_WZ,
      w[34]);
  FFV2_7_3(w[3], w[2], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[35]);
  FFV2_3(w[3], w[2], pars->GC_138, pars->mdl_MZ, pars->mdl_WZ, w[36]); 
  FFV1_2(w[5], w[14], pars->GC_7, pars->ZERO, pars->ZERO, w[37]); 
  FFV1_1(w[1], w[14], pars->GC_7, pars->ZERO, pars->ZERO, w[38]); 
  FFV2_4_2(w[5], w[19], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[39]);
  FFV2_4_1(w[1], w[19], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[40]);
  FFFF4_5_4(w[5], w[2], w[3], pars->GC_15, pars->GC_16, pars->ZERO, pars->ZERO,
      w[41]);
  FFFF1_4(w[5], w[2], w[3], -pars->GC_12, pars->ZERO, pars->ZERO, w[42]); 
  FFFF1_2_6_4(w[5], w[2], w[3], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->ZERO, pars->ZERO, w[43]);
  FFFF4_5_3(w[5], w[2], w[1], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[44]);
  FFFF1_3(w[5], w[2], w[1], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[45]); 
  FFFF1_2_6_3(w[5], w[2], w[1], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[46]);
  FFFF4_5_2(w[5], w[3], w[1], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[47]);
  FFFF1_2(w[5], w[3], w[1], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[48]); 
  FFFF1_2_6_2(w[5], w[3], w[1], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[49]);
  FFFF4_5_1(w[2], w[3], w[1], pars->GC_15, pars->GC_16, pars->ZERO, pars->ZERO,
      w[50]);
  FFFF1_1(w[2], w[3], w[1], -pars->GC_12, pars->ZERO, pars->ZERO, w[51]); 
  FFFF1_2_6_1(w[2], w[3], w[1], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->ZERO, pars->ZERO, w[52]);
  FFV1P0_3(w[4], w[1], pars->GC_1, pars->ZERO, pars->ZERO, w[53]); 
  FFV1P0_3(w[5], w[0], pars->GC_1, pars->ZERO, pars->ZERO, w[54]); 
  FFV5_1(w[2], w[53], pars->GC_2, pars->mdl_MT, pars->mdl_WT, w[55]); 
  FFV3_8_1(w[2], w[53], pars->GC_170, pars->GC_91, pars->mdl_MT, pars->mdl_WT,
      w[56]);
  FFV5_2(w[3], w[53], pars->GC_2, pars->mdl_MT, pars->mdl_WT, w[57]); 
  FFV3_8_2(w[3], w[53], pars->GC_170, pars->GC_91, pars->mdl_MT, pars->mdl_WT,
      w[58]);
  FFV1P0_3(w[5], w[0], pars->GC_7, pars->ZERO, pars->ZERO, w[59]); 
  FFV2_4_3(w[5], w[0], pars->GC_58, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[60]);
  FFV1P0_3(w[4], w[1], pars->GC_7, pars->ZERO, pars->ZERO, w[61]); 
  FFV5_1(w[2], w[61], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[62]); 
  FFV3_8_1(w[2], w[61], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[63]);
  FFV5_2(w[3], w[61], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[64]); 
  FFV3_8_2(w[3], w[61], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[65]);
  FFV2_4_3(w[4], w[1], pars->GC_58, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[66]);
  FFV3_8_1(w[2], w[66], pars->GC_166, pars->GC_81, pars->mdl_MT, pars->mdl_WT,
      w[67]);
  FFV2_7_1(w[2], w[66], pars->GC_59, pars->GC_69, pars->mdl_MT, pars->mdl_WT,
      w[68]);
  FFV2_1(w[2], w[66], pars->GC_138, pars->mdl_MT, pars->mdl_WT, w[69]); 
  FFV3_8_2(w[3], w[66], pars->GC_166, pars->GC_81, pars->mdl_MT, pars->mdl_WT,
      w[70]);
  FFV2_7_2(w[3], w[66], pars->GC_59, pars->GC_69, pars->mdl_MT, pars->mdl_WT,
      w[71]);
  FFV2_2(w[3], w[66], pars->GC_138, pars->mdl_MT, pars->mdl_WT, w[72]); 
  FFV1_2(w[5], w[53], pars->GC_1, pars->ZERO, pars->ZERO, w[73]); 
  FFV1_1(w[0], w[53], pars->GC_1, pars->ZERO, pars->ZERO, w[74]); 
  FFV1_2(w[5], w[61], pars->GC_7, pars->ZERO, pars->ZERO, w[75]); 
  FFV1_1(w[0], w[61], pars->GC_7, pars->ZERO, pars->ZERO, w[76]); 
  FFV2_4_2(w[5], w[66], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[77]);
  FFV2_4_1(w[0], w[66], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[78]);
  FFFF4_5_3(w[5], w[2], w[0], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[79]);
  FFFF1_3(w[5], w[2], w[0], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[80]); 
  FFFF1_2_6_3(w[5], w[2], w[0], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[81]);
  FFFF4_5_2(w[5], w[3], w[0], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[82]);
  FFFF1_2(w[5], w[3], w[0], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[83]); 
  FFFF1_2_6_2(w[5], w[3], w[0], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[84]);
  FFFF4_5_1(w[2], w[3], w[0], pars->GC_15, pars->GC_16, pars->ZERO, pars->ZERO,
      w[85]);
  FFFF1_1(w[2], w[3], w[0], -pars->GC_12, pars->ZERO, pars->ZERO, w[86]); 
  FFFF1_2_6_1(w[2], w[3], w[0], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->ZERO, pars->ZERO, w[87]);
  FFV1_2(w[4], w[54], pars->GC_1, pars->ZERO, pars->ZERO, w[88]); 
  FFV1_2(w[4], w[30], pars->GC_1, pars->ZERO, pars->ZERO, w[89]); 
  FFV1_2(w[4], w[33], pars->GC_1, pars->ZERO, pars->ZERO, w[90]); 
  FFV1_2(w[4], w[27], pars->GC_7, pars->ZERO, pars->ZERO, w[91]); 
  FFV1_2(w[4], w[28], pars->GC_7, pars->ZERO, pars->ZERO, w[92]); 
  FFV2_4_2(w[4], w[34], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[93]);
  FFV2_4_2(w[4], w[35], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[94]);
  FFV2_4_2(w[4], w[36], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[95]);
  FFV1_2(w[4], w[59], pars->GC_7, pars->ZERO, pars->ZERO, w[96]); 
  FFV2_4_2(w[4], w[60], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[97]);
  FFV1_2(w[4], w[7], pars->GC_1, pars->ZERO, pars->ZERO, w[98]); 
  FFV1_2(w[4], w[12], pars->GC_7, pars->ZERO, pars->ZERO, w[99]); 
  FFV2_4_2(w[4], w[13], pars->GC_58, pars->GC_69, pars->ZERO, pars->ZERO,
      w[100]);
  FFFF4_5_4(w[4], w[2], w[3], pars->GC_15, pars->GC_16, pars->ZERO, pars->ZERO,
      w[101]);
  FFFF1_4(w[4], w[2], w[3], -pars->GC_12, pars->ZERO, pars->ZERO, w[102]); 
  FFFF1_2_6_4(w[4], w[2], w[3], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->ZERO, pars->ZERO, w[103]);
  FFFF4_5_3(w[4], w[2], w[0], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[104]);
  FFFF1_3(w[4], w[2], w[0], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[105]); 
  FFFF1_2_6_3(w[4], w[2], w[0], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[106]);
  FFFF4_5_3(w[4], w[2], w[1], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[107]);
  FFFF1_3(w[4], w[2], w[1], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[108]); 
  FFFF1_2_6_3(w[4], w[2], w[1], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[109]);
  FFFF4_5_2(w[4], w[3], w[0], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[110]);
  FFFF1_2(w[4], w[3], w[0], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[111]); 
  FFFF1_2_6_2(w[4], w[3], w[0], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[112]);
  FFFF4_5_2(w[4], w[3], w[1], pars->GC_15, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[113]);
  FFFF1_2(w[4], w[3], w[1], -pars->GC_12, pars->mdl_MT, pars->mdl_WT, w[114]); 
  FFFF1_2_6_2(w[4], w[3], w[1], -pars->GC_24, -pars->GC_23, -pars->GC_26,
      pars->mdl_MT, pars->mdl_WT, w[115]);

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
  FFV1_0(w[31], w[1], w[30], pars->GC_1, amp[73]); 
  FFV1_0(w[5], w[32], w[30], pars->GC_1, amp[74]); 
  FFV1_0(w[31], w[1], w[33], pars->GC_1, amp[75]); 
  FFV1_0(w[5], w[32], w[33], pars->GC_1, amp[76]); 
  FFV1_0(w[31], w[1], w[27], pars->GC_7, amp[77]); 
  FFV1_0(w[5], w[32], w[27], pars->GC_7, amp[78]); 
  FFV1_0(w[31], w[1], w[28], pars->GC_7, amp[79]); 
  FFV1_0(w[5], w[32], w[28], pars->GC_7, amp[80]); 
  FFV2_4_0(w[31], w[1], w[34], pars->GC_58, pars->GC_69, amp[81]); 
  FFV2_4_0(w[5], w[32], w[34], pars->GC_58, pars->GC_69, amp[82]); 
  FFV2_4_0(w[31], w[1], w[35], pars->GC_58, pars->GC_69, amp[83]); 
  FFV2_4_0(w[5], w[32], w[35], pars->GC_58, pars->GC_69, amp[84]); 
  FFV2_4_0(w[31], w[1], w[36], pars->GC_58, pars->GC_69, amp[85]); 
  FFV2_4_0(w[5], w[32], w[36], pars->GC_58, pars->GC_69, amp[86]); 
  FFV1_0(w[37], w[1], w[30], pars->GC_1, amp[87]); 
  FFV1_0(w[5], w[38], w[30], pars->GC_1, amp[88]); 
  FFV1_0(w[37], w[1], w[33], pars->GC_1, amp[89]); 
  FFV1_0(w[5], w[38], w[33], pars->GC_1, amp[90]); 
  FFV1_0(w[37], w[1], w[27], pars->GC_7, amp[91]); 
  FFV1_0(w[5], w[38], w[27], pars->GC_7, amp[92]); 
  FFV1_0(w[37], w[1], w[28], pars->GC_7, amp[93]); 
  FFV1_0(w[5], w[38], w[28], pars->GC_7, amp[94]); 
  FFV2_4_0(w[37], w[1], w[34], pars->GC_58, pars->GC_69, amp[95]); 
  FFV2_4_0(w[5], w[38], w[34], pars->GC_58, pars->GC_69, amp[96]); 
  FFV2_4_0(w[37], w[1], w[35], pars->GC_58, pars->GC_69, amp[97]); 
  FFV2_4_0(w[5], w[38], w[35], pars->GC_58, pars->GC_69, amp[98]); 
  FFV2_4_0(w[37], w[1], w[36], pars->GC_58, pars->GC_69, amp[99]); 
  FFV2_4_0(w[5], w[38], w[36], pars->GC_58, pars->GC_69, amp[100]); 
  FFV1_0(w[39], w[1], w[30], pars->GC_1, amp[101]); 
  FFV1_0(w[5], w[40], w[30], pars->GC_1, amp[102]); 
  FFV1_0(w[39], w[1], w[33], pars->GC_1, amp[103]); 
  FFV1_0(w[5], w[40], w[33], pars->GC_1, amp[104]); 
  FFV1_0(w[39], w[1], w[27], pars->GC_7, amp[105]); 
  FFV1_0(w[5], w[40], w[27], pars->GC_7, amp[106]); 
  FFV1_0(w[39], w[1], w[28], pars->GC_7, amp[107]); 
  FFV1_0(w[5], w[40], w[28], pars->GC_7, amp[108]); 
  FFV2_4_0(w[39], w[1], w[34], pars->GC_58, pars->GC_69, amp[109]); 
  FFV2_4_0(w[5], w[40], w[34], pars->GC_58, pars->GC_69, amp[110]); 
  FFV2_4_0(w[39], w[1], w[35], pars->GC_58, pars->GC_69, amp[111]); 
  FFV2_4_0(w[5], w[40], w[35], pars->GC_58, pars->GC_69, amp[112]); 
  FFV2_4_0(w[39], w[1], w[36], pars->GC_58, pars->GC_69, amp[113]); 
  FFV2_4_0(w[5], w[40], w[36], pars->GC_58, pars->GC_69, amp[114]); 
  FFV1_0(w[41], w[1], w[6], pars->GC_1, amp[115]); 
  FFV1_0(w[42], w[1], w[6], pars->GC_1, amp[116]); 
  FFV1_0(w[43], w[1], w[6], pars->GC_1, amp[117]); 
  FFV1_0(w[41], w[1], w[14], pars->GC_7, amp[118]); 
  FFV1_0(w[42], w[1], w[14], pars->GC_7, amp[119]); 
  FFV1_0(w[43], w[1], w[14], pars->GC_7, amp[120]); 
  FFV2_4_0(w[41], w[1], w[19], pars->GC_58, pars->GC_69, amp[121]); 
  FFV2_4_0(w[42], w[1], w[19], pars->GC_58, pars->GC_69, amp[122]); 
  FFV2_4_0(w[43], w[1], w[19], pars->GC_58, pars->GC_69, amp[123]); 
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
  FFV1_0(w[5], w[50], w[6], pars->GC_1, amp[142]); 
  FFV1_0(w[5], w[51], w[6], pars->GC_1, amp[143]); 
  FFV1_0(w[5], w[52], w[6], pars->GC_1, amp[144]); 
  FFV1_0(w[5], w[50], w[14], pars->GC_7, amp[145]); 
  FFV1_0(w[5], w[51], w[14], pars->GC_7, amp[146]); 
  FFV1_0(w[5], w[52], w[14], pars->GC_7, amp[147]); 
  FFV2_4_0(w[5], w[50], w[19], pars->GC_58, pars->GC_69, amp[148]); 
  FFV2_4_0(w[5], w[51], w[19], pars->GC_58, pars->GC_69, amp[149]); 
  FFV2_4_0(w[5], w[52], w[19], pars->GC_58, pars->GC_69, amp[150]); 
  FFV5_0(w[3], w[55], w[54], pars->GC_2, amp[151]); 
  FFV3_8_0(w[3], w[55], w[54], pars->GC_170, pars->GC_91, amp[152]); 
  FFV5_0(w[3], w[56], w[54], pars->GC_2, amp[153]); 
  FFV5_0(w[57], w[2], w[54], pars->GC_2, amp[154]); 
  FFV3_8_0(w[57], w[2], w[54], pars->GC_170, pars->GC_91, amp[155]); 
  FFV5_0(w[58], w[2], w[54], pars->GC_2, amp[156]); 
  FFV5_0(w[3], w[55], w[59], pars->GC_7, amp[157]); 
  FFV3_8_0(w[3], w[55], w[59], pars->GC_145, pars->GC_79, amp[158]); 
  FFV5_0(w[3], w[56], w[59], pars->GC_7, amp[159]); 
  FFV5_0(w[57], w[2], w[59], pars->GC_7, amp[160]); 
  FFV3_8_0(w[57], w[2], w[59], pars->GC_145, pars->GC_79, amp[161]); 
  FFV5_0(w[58], w[2], w[59], pars->GC_7, amp[162]); 
  FFV3_8_0(w[3], w[55], w[60], pars->GC_166, pars->GC_81, amp[163]); 
  FFV2_7_0(w[3], w[55], w[60], pars->GC_59, pars->GC_69, amp[164]); 
  FFV2_0(w[3], w[55], w[60], pars->GC_138, amp[165]); 
  FFV2_7_0(w[3], w[56], w[60], pars->GC_59, pars->GC_69, amp[166]); 
  FFV3_8_0(w[57], w[2], w[60], pars->GC_166, pars->GC_81, amp[167]); 
  FFV2_7_0(w[57], w[2], w[60], pars->GC_59, pars->GC_69, amp[168]); 
  FFV2_0(w[57], w[2], w[60], pars->GC_138, amp[169]); 
  FFV2_7_0(w[58], w[2], w[60], pars->GC_59, pars->GC_69, amp[170]); 
  FFV5_0(w[3], w[62], w[54], pars->GC_2, amp[171]); 
  FFV3_8_0(w[3], w[62], w[54], pars->GC_170, pars->GC_91, amp[172]); 
  FFV5_0(w[3], w[63], w[54], pars->GC_2, amp[173]); 
  FFV5_0(w[64], w[2], w[54], pars->GC_2, amp[174]); 
  FFV3_8_0(w[64], w[2], w[54], pars->GC_170, pars->GC_91, amp[175]); 
  FFV5_0(w[65], w[2], w[54], pars->GC_2, amp[176]); 
  FFVV1_2_0(w[3], w[2], w[61], w[59], pars->GC_146, pars->GC_85, amp[177]); 
  FFV5_0(w[3], w[62], w[59], pars->GC_7, amp[178]); 
  FFV3_8_0(w[3], w[62], w[59], pars->GC_145, pars->GC_79, amp[179]); 
  FFV5_0(w[3], w[63], w[59], pars->GC_7, amp[180]); 
  FFV5_0(w[64], w[2], w[59], pars->GC_7, amp[181]); 
  FFV3_8_0(w[64], w[2], w[59], pars->GC_145, pars->GC_79, amp[182]); 
  FFV5_0(w[65], w[2], w[59], pars->GC_7, amp[183]); 
  FFV3_8_0(w[3], w[62], w[60], pars->GC_166, pars->GC_81, amp[184]); 
  FFV2_7_0(w[3], w[62], w[60], pars->GC_59, pars->GC_69, amp[185]); 
  FFV2_0(w[3], w[62], w[60], pars->GC_138, amp[186]); 
  FFV2_7_0(w[3], w[63], w[60], pars->GC_59, pars->GC_69, amp[187]); 
  FFV3_8_0(w[64], w[2], w[60], pars->GC_166, pars->GC_81, amp[188]); 
  FFV2_7_0(w[64], w[2], w[60], pars->GC_59, pars->GC_69, amp[189]); 
  FFV2_0(w[64], w[2], w[60], pars->GC_138, amp[190]); 
  FFV2_7_0(w[65], w[2], w[60], pars->GC_59, pars->GC_69, amp[191]); 
  FFV5_0(w[3], w[67], w[54], pars->GC_2, amp[192]); 
  FFV5_0(w[3], w[68], w[54], pars->GC_2, amp[193]); 
  FFV3_8_0(w[3], w[68], w[54], pars->GC_170, pars->GC_91, amp[194]); 
  FFV5_0(w[3], w[69], w[54], pars->GC_2, amp[195]); 
  FFV5_0(w[70], w[2], w[54], pars->GC_2, amp[196]); 
  FFV5_0(w[71], w[2], w[54], pars->GC_2, amp[197]); 
  FFV3_8_0(w[71], w[2], w[54], pars->GC_170, pars->GC_91, amp[198]); 
  FFV5_0(w[72], w[2], w[54], pars->GC_2, amp[199]); 
  FFV5_0(w[3], w[67], w[59], pars->GC_7, amp[200]); 
  FFV5_0(w[3], w[68], w[59], pars->GC_7, amp[201]); 
  FFV3_8_0(w[3], w[68], w[59], pars->GC_145, pars->GC_79, amp[202]); 
  FFV5_0(w[3], w[69], w[59], pars->GC_7, amp[203]); 
  FFV5_0(w[70], w[2], w[59], pars->GC_7, amp[204]); 
  FFV5_0(w[71], w[2], w[59], pars->GC_7, amp[205]); 
  FFV3_8_0(w[71], w[2], w[59], pars->GC_145, pars->GC_79, amp[206]); 
  FFV5_0(w[72], w[2], w[59], pars->GC_7, amp[207]); 
  FFV2_7_0(w[3], w[67], w[60], pars->GC_59, pars->GC_69, amp[208]); 
  FFV3_8_0(w[3], w[68], w[60], pars->GC_166, pars->GC_81, amp[209]); 
  FFV2_7_0(w[3], w[68], w[60], pars->GC_59, pars->GC_69, amp[210]); 
  FFV2_0(w[3], w[68], w[60], pars->GC_138, amp[211]); 
  FFV2_7_0(w[3], w[69], w[60], pars->GC_59, pars->GC_69, amp[212]); 
  FFV2_7_0(w[70], w[2], w[60], pars->GC_59, pars->GC_69, amp[213]); 
  FFV3_8_0(w[71], w[2], w[60], pars->GC_166, pars->GC_81, amp[214]); 
  FFV2_7_0(w[71], w[2], w[60], pars->GC_59, pars->GC_69, amp[215]); 
  FFV2_0(w[71], w[2], w[60], pars->GC_138, amp[216]); 
  FFV2_7_0(w[72], w[2], w[60], pars->GC_59, pars->GC_69, amp[217]); 
  VVS2_0(w[61], w[59], w[26], pars->GC_78, amp[218]); 
  VVV2_0(w[61], w[59], w[27], pars->GC_28, amp[219]); 
  VVV1_0(w[61], w[59], w[27], pars->GC_6, amp[220]); 
  VVV1_0(w[61], w[59], w[28], pars->GC_6, amp[221]); 
  VVS1_0(w[66], w[60], w[26], pars->GC_93, amp[222]); 
  VVS1_0(w[66], w[60], w[29], pars->GC_93, amp[223]); 
  FFV1_0(w[73], w[0], w[30], pars->GC_1, amp[224]); 
  FFV1_0(w[5], w[74], w[30], pars->GC_1, amp[225]); 
  FFV1_0(w[73], w[0], w[33], pars->GC_1, amp[226]); 
  FFV1_0(w[5], w[74], w[33], pars->GC_1, amp[227]); 
  FFV1_0(w[73], w[0], w[27], pars->GC_7, amp[228]); 
  FFV1_0(w[5], w[74], w[27], pars->GC_7, amp[229]); 
  FFV1_0(w[73], w[0], w[28], pars->GC_7, amp[230]); 
  FFV1_0(w[5], w[74], w[28], pars->GC_7, amp[231]); 
  FFV2_4_0(w[73], w[0], w[34], pars->GC_58, pars->GC_69, amp[232]); 
  FFV2_4_0(w[5], w[74], w[34], pars->GC_58, pars->GC_69, amp[233]); 
  FFV2_4_0(w[73], w[0], w[35], pars->GC_58, pars->GC_69, amp[234]); 
  FFV2_4_0(w[5], w[74], w[35], pars->GC_58, pars->GC_69, amp[235]); 
  FFV2_4_0(w[73], w[0], w[36], pars->GC_58, pars->GC_69, amp[236]); 
  FFV2_4_0(w[5], w[74], w[36], pars->GC_58, pars->GC_69, amp[237]); 
  FFV1_0(w[75], w[0], w[30], pars->GC_1, amp[238]); 
  FFV1_0(w[5], w[76], w[30], pars->GC_1, amp[239]); 
  FFV1_0(w[75], w[0], w[33], pars->GC_1, amp[240]); 
  FFV1_0(w[5], w[76], w[33], pars->GC_1, amp[241]); 
  FFV1_0(w[75], w[0], w[27], pars->GC_7, amp[242]); 
  FFV1_0(w[5], w[76], w[27], pars->GC_7, amp[243]); 
  FFV1_0(w[75], w[0], w[28], pars->GC_7, amp[244]); 
  FFV1_0(w[5], w[76], w[28], pars->GC_7, amp[245]); 
  FFV2_4_0(w[75], w[0], w[34], pars->GC_58, pars->GC_69, amp[246]); 
  FFV2_4_0(w[5], w[76], w[34], pars->GC_58, pars->GC_69, amp[247]); 
  FFV2_4_0(w[75], w[0], w[35], pars->GC_58, pars->GC_69, amp[248]); 
  FFV2_4_0(w[5], w[76], w[35], pars->GC_58, pars->GC_69, amp[249]); 
  FFV2_4_0(w[75], w[0], w[36], pars->GC_58, pars->GC_69, amp[250]); 
  FFV2_4_0(w[5], w[76], w[36], pars->GC_58, pars->GC_69, amp[251]); 
  FFV1_0(w[77], w[0], w[30], pars->GC_1, amp[252]); 
  FFV1_0(w[5], w[78], w[30], pars->GC_1, amp[253]); 
  FFV1_0(w[77], w[0], w[33], pars->GC_1, amp[254]); 
  FFV1_0(w[5], w[78], w[33], pars->GC_1, amp[255]); 
  FFV1_0(w[77], w[0], w[27], pars->GC_7, amp[256]); 
  FFV1_0(w[5], w[78], w[27], pars->GC_7, amp[257]); 
  FFV1_0(w[77], w[0], w[28], pars->GC_7, amp[258]); 
  FFV1_0(w[5], w[78], w[28], pars->GC_7, amp[259]); 
  FFV2_4_0(w[77], w[0], w[34], pars->GC_58, pars->GC_69, amp[260]); 
  FFV2_4_0(w[5], w[78], w[34], pars->GC_58, pars->GC_69, amp[261]); 
  FFV2_4_0(w[77], w[0], w[35], pars->GC_58, pars->GC_69, amp[262]); 
  FFV2_4_0(w[5], w[78], w[35], pars->GC_58, pars->GC_69, amp[263]); 
  FFV2_4_0(w[77], w[0], w[36], pars->GC_58, pars->GC_69, amp[264]); 
  FFV2_4_0(w[5], w[78], w[36], pars->GC_58, pars->GC_69, amp[265]); 
  FFV1_0(w[41], w[0], w[53], pars->GC_1, amp[266]); 
  FFV1_0(w[42], w[0], w[53], pars->GC_1, amp[267]); 
  FFV1_0(w[43], w[0], w[53], pars->GC_1, amp[268]); 
  FFV1_0(w[41], w[0], w[61], pars->GC_7, amp[269]); 
  FFV1_0(w[42], w[0], w[61], pars->GC_7, amp[270]); 
  FFV1_0(w[43], w[0], w[61], pars->GC_7, amp[271]); 
  FFV2_4_0(w[41], w[0], w[66], pars->GC_58, pars->GC_69, amp[272]); 
  FFV2_4_0(w[42], w[0], w[66], pars->GC_58, pars->GC_69, amp[273]); 
  FFV2_4_0(w[43], w[0], w[66], pars->GC_58, pars->GC_69, amp[274]); 
  FFV5_0(w[3], w[79], w[53], pars->GC_2, amp[275]); 
  FFV5_0(w[3], w[80], w[53], pars->GC_2, amp[276]); 
  FFV5_0(w[3], w[81], w[53], pars->GC_2, amp[277]); 
  FFV5_0(w[3], w[79], w[61], pars->GC_7, amp[278]); 
  FFV5_0(w[3], w[80], w[61], pars->GC_7, amp[279]); 
  FFV5_0(w[3], w[81], w[61], pars->GC_7, amp[280]); 
  FFV2_7_0(w[3], w[79], w[66], pars->GC_59, pars->GC_69, amp[281]); 
  FFV2_7_0(w[3], w[80], w[66], pars->GC_59, pars->GC_69, amp[282]); 
  FFV2_7_0(w[3], w[81], w[66], pars->GC_59, pars->GC_69, amp[283]); 
  FFV5_0(w[82], w[2], w[53], pars->GC_2, amp[284]); 
  FFV5_0(w[83], w[2], w[53], pars->GC_2, amp[285]); 
  FFV5_0(w[84], w[2], w[53], pars->GC_2, amp[286]); 
  FFV5_0(w[82], w[2], w[61], pars->GC_7, amp[287]); 
  FFV5_0(w[83], w[2], w[61], pars->GC_7, amp[288]); 
  FFV5_0(w[84], w[2], w[61], pars->GC_7, amp[289]); 
  FFV2_7_0(w[82], w[2], w[66], pars->GC_59, pars->GC_69, amp[290]); 
  FFV2_7_0(w[83], w[2], w[66], pars->GC_59, pars->GC_69, amp[291]); 
  FFV2_7_0(w[84], w[2], w[66], pars->GC_59, pars->GC_69, amp[292]); 
  FFV1_0(w[5], w[85], w[53], pars->GC_1, amp[293]); 
  FFV1_0(w[5], w[86], w[53], pars->GC_1, amp[294]); 
  FFV1_0(w[5], w[87], w[53], pars->GC_1, amp[295]); 
  FFV1_0(w[5], w[85], w[61], pars->GC_7, amp[296]); 
  FFV1_0(w[5], w[86], w[61], pars->GC_7, amp[297]); 
  FFV1_0(w[5], w[87], w[61], pars->GC_7, amp[298]); 
  FFV2_4_0(w[5], w[85], w[66], pars->GC_58, pars->GC_69, amp[299]); 
  FFV2_4_0(w[5], w[86], w[66], pars->GC_58, pars->GC_69, amp[300]); 
  FFV2_4_0(w[5], w[87], w[66], pars->GC_58, pars->GC_69, amp[301]); 
  FFV1_0(w[88], w[1], w[30], pars->GC_1, amp[302]); 
  FFV1_0(w[89], w[1], w[54], pars->GC_1, amp[303]); 
  FFV1_0(w[88], w[1], w[33], pars->GC_1, amp[304]); 
  FFV1_0(w[90], w[1], w[54], pars->GC_1, amp[305]); 
  FFV1_0(w[88], w[1], w[27], pars->GC_7, amp[306]); 
  FFV1_0(w[91], w[1], w[54], pars->GC_1, amp[307]); 
  FFV1_0(w[88], w[1], w[28], pars->GC_7, amp[308]); 
  FFV1_0(w[92], w[1], w[54], pars->GC_1, amp[309]); 
  FFV2_4_0(w[88], w[1], w[34], pars->GC_58, pars->GC_69, amp[310]); 
  FFV1_0(w[93], w[1], w[54], pars->GC_1, amp[311]); 
  FFV2_4_0(w[88], w[1], w[35], pars->GC_58, pars->GC_69, amp[312]); 
  FFV1_0(w[94], w[1], w[54], pars->GC_1, amp[313]); 
  FFV2_4_0(w[88], w[1], w[36], pars->GC_58, pars->GC_69, amp[314]); 
  FFV1_0(w[95], w[1], w[54], pars->GC_1, amp[315]); 
  FFV1_0(w[96], w[1], w[30], pars->GC_1, amp[316]); 
  FFV1_0(w[89], w[1], w[59], pars->GC_7, amp[317]); 
  FFV1_0(w[96], w[1], w[33], pars->GC_1, amp[318]); 
  FFV1_0(w[90], w[1], w[59], pars->GC_7, amp[319]); 
  FFV1_0(w[96], w[1], w[27], pars->GC_7, amp[320]); 
  FFV1_0(w[91], w[1], w[59], pars->GC_7, amp[321]); 
  FFV1_0(w[96], w[1], w[28], pars->GC_7, amp[322]); 
  FFV1_0(w[92], w[1], w[59], pars->GC_7, amp[323]); 
  FFV2_4_0(w[96], w[1], w[34], pars->GC_58, pars->GC_69, amp[324]); 
  FFV1_0(w[93], w[1], w[59], pars->GC_7, amp[325]); 
  FFV2_4_0(w[96], w[1], w[35], pars->GC_58, pars->GC_69, amp[326]); 
  FFV1_0(w[94], w[1], w[59], pars->GC_7, amp[327]); 
  FFV2_4_0(w[96], w[1], w[36], pars->GC_58, pars->GC_69, amp[328]); 
  FFV1_0(w[95], w[1], w[59], pars->GC_7, amp[329]); 
  FFV1_0(w[97], w[1], w[30], pars->GC_1, amp[330]); 
  FFV2_4_0(w[89], w[1], w[60], pars->GC_58, pars->GC_69, amp[331]); 
  FFV1_0(w[97], w[1], w[33], pars->GC_1, amp[332]); 
  FFV2_4_0(w[90], w[1], w[60], pars->GC_58, pars->GC_69, amp[333]); 
  FFV1_0(w[97], w[1], w[27], pars->GC_7, amp[334]); 
  FFV2_4_0(w[91], w[1], w[60], pars->GC_58, pars->GC_69, amp[335]); 
  FFV1_0(w[97], w[1], w[28], pars->GC_7, amp[336]); 
  FFV2_4_0(w[92], w[1], w[60], pars->GC_58, pars->GC_69, amp[337]); 
  FFV2_4_0(w[97], w[1], w[34], pars->GC_58, pars->GC_69, amp[338]); 
  FFV2_4_0(w[93], w[1], w[60], pars->GC_58, pars->GC_69, amp[339]); 
  FFV2_4_0(w[97], w[1], w[35], pars->GC_58, pars->GC_69, amp[340]); 
  FFV2_4_0(w[94], w[1], w[60], pars->GC_58, pars->GC_69, amp[341]); 
  FFV2_4_0(w[97], w[1], w[36], pars->GC_58, pars->GC_69, amp[342]); 
  FFV2_4_0(w[95], w[1], w[60], pars->GC_58, pars->GC_69, amp[343]); 
  FFV1_0(w[4], w[50], w[54], pars->GC_1, amp[344]); 
  FFV1_0(w[4], w[51], w[54], pars->GC_1, amp[345]); 
  FFV1_0(w[4], w[52], w[54], pars->GC_1, amp[346]); 
  FFV1_0(w[4], w[50], w[59], pars->GC_7, amp[347]); 
  FFV1_0(w[4], w[51], w[59], pars->GC_7, amp[348]); 
  FFV1_0(w[4], w[52], w[59], pars->GC_7, amp[349]); 
  FFV2_4_0(w[4], w[50], w[60], pars->GC_58, pars->GC_69, amp[350]); 
  FFV2_4_0(w[4], w[51], w[60], pars->GC_58, pars->GC_69, amp[351]); 
  FFV2_4_0(w[4], w[52], w[60], pars->GC_58, pars->GC_69, amp[352]); 
  FFV1_0(w[98], w[0], w[30], pars->GC_1, amp[353]); 
  FFV1_0(w[89], w[0], w[7], pars->GC_1, amp[354]); 
  FFV1_0(w[98], w[0], w[33], pars->GC_1, amp[355]); 
  FFV1_0(w[90], w[0], w[7], pars->GC_1, amp[356]); 
  FFV1_0(w[98], w[0], w[27], pars->GC_7, amp[357]); 
  FFV1_0(w[91], w[0], w[7], pars->GC_1, amp[358]); 
  FFV1_0(w[98], w[0], w[28], pars->GC_7, amp[359]); 
  FFV1_0(w[92], w[0], w[7], pars->GC_1, amp[360]); 
  FFV2_4_0(w[98], w[0], w[34], pars->GC_58, pars->GC_69, amp[361]); 
  FFV1_0(w[93], w[0], w[7], pars->GC_1, amp[362]); 
  FFV2_4_0(w[98], w[0], w[35], pars->GC_58, pars->GC_69, amp[363]); 
  FFV1_0(w[94], w[0], w[7], pars->GC_1, amp[364]); 
  FFV2_4_0(w[98], w[0], w[36], pars->GC_58, pars->GC_69, amp[365]); 
  FFV1_0(w[95], w[0], w[7], pars->GC_1, amp[366]); 
  FFV1_0(w[99], w[0], w[30], pars->GC_1, amp[367]); 
  FFV1_0(w[89], w[0], w[12], pars->GC_7, amp[368]); 
  FFV1_0(w[99], w[0], w[33], pars->GC_1, amp[369]); 
  FFV1_0(w[90], w[0], w[12], pars->GC_7, amp[370]); 
  FFV1_0(w[99], w[0], w[27], pars->GC_7, amp[371]); 
  FFV1_0(w[91], w[0], w[12], pars->GC_7, amp[372]); 
  FFV1_0(w[99], w[0], w[28], pars->GC_7, amp[373]); 
  FFV1_0(w[92], w[0], w[12], pars->GC_7, amp[374]); 
  FFV2_4_0(w[99], w[0], w[34], pars->GC_58, pars->GC_69, amp[375]); 
  FFV1_0(w[93], w[0], w[12], pars->GC_7, amp[376]); 
  FFV2_4_0(w[99], w[0], w[35], pars->GC_58, pars->GC_69, amp[377]); 
  FFV1_0(w[94], w[0], w[12], pars->GC_7, amp[378]); 
  FFV2_4_0(w[99], w[0], w[36], pars->GC_58, pars->GC_69, amp[379]); 
  FFV1_0(w[95], w[0], w[12], pars->GC_7, amp[380]); 
  FFV1_0(w[100], w[0], w[30], pars->GC_1, amp[381]); 
  FFV2_4_0(w[89], w[0], w[13], pars->GC_58, pars->GC_69, amp[382]); 
  FFV1_0(w[100], w[0], w[33], pars->GC_1, amp[383]); 
  FFV2_4_0(w[90], w[0], w[13], pars->GC_58, pars->GC_69, amp[384]); 
  FFV1_0(w[100], w[0], w[27], pars->GC_7, amp[385]); 
  FFV2_4_0(w[91], w[0], w[13], pars->GC_58, pars->GC_69, amp[386]); 
  FFV1_0(w[100], w[0], w[28], pars->GC_7, amp[387]); 
  FFV2_4_0(w[92], w[0], w[13], pars->GC_58, pars->GC_69, amp[388]); 
  FFV2_4_0(w[100], w[0], w[34], pars->GC_58, pars->GC_69, amp[389]); 
  FFV2_4_0(w[93], w[0], w[13], pars->GC_58, pars->GC_69, amp[390]); 
  FFV2_4_0(w[100], w[0], w[35], pars->GC_58, pars->GC_69, amp[391]); 
  FFV2_4_0(w[94], w[0], w[13], pars->GC_58, pars->GC_69, amp[392]); 
  FFV2_4_0(w[100], w[0], w[36], pars->GC_58, pars->GC_69, amp[393]); 
  FFV2_4_0(w[95], w[0], w[13], pars->GC_58, pars->GC_69, amp[394]); 
  FFV1_0(w[4], w[85], w[7], pars->GC_1, amp[395]); 
  FFV1_0(w[4], w[86], w[7], pars->GC_1, amp[396]); 
  FFV1_0(w[4], w[87], w[7], pars->GC_1, amp[397]); 
  FFV1_0(w[4], w[85], w[12], pars->GC_7, amp[398]); 
  FFV1_0(w[4], w[86], w[12], pars->GC_7, amp[399]); 
  FFV1_0(w[4], w[87], w[12], pars->GC_7, amp[400]); 
  FFV2_4_0(w[4], w[85], w[13], pars->GC_58, pars->GC_69, amp[401]); 
  FFV2_4_0(w[4], w[86], w[13], pars->GC_58, pars->GC_69, amp[402]); 
  FFV2_4_0(w[4], w[87], w[13], pars->GC_58, pars->GC_69, amp[403]); 
  FFV1_0(w[101], w[1], w[54], pars->GC_1, amp[404]); 
  FFV1_0(w[102], w[1], w[54], pars->GC_1, amp[405]); 
  FFV1_0(w[103], w[1], w[54], pars->GC_1, amp[406]); 
  FFV1_0(w[101], w[1], w[59], pars->GC_7, amp[407]); 
  FFV1_0(w[102], w[1], w[59], pars->GC_7, amp[408]); 
  FFV1_0(w[103], w[1], w[59], pars->GC_7, amp[409]); 
  FFV2_4_0(w[101], w[1], w[60], pars->GC_58, pars->GC_69, amp[410]); 
  FFV2_4_0(w[102], w[1], w[60], pars->GC_58, pars->GC_69, amp[411]); 
  FFV2_4_0(w[103], w[1], w[60], pars->GC_58, pars->GC_69, amp[412]); 
  FFV1_0(w[101], w[0], w[7], pars->GC_1, amp[413]); 
  FFV1_0(w[102], w[0], w[7], pars->GC_1, amp[414]); 
  FFV1_0(w[103], w[0], w[7], pars->GC_1, amp[415]); 
  FFV1_0(w[101], w[0], w[12], pars->GC_7, amp[416]); 
  FFV1_0(w[102], w[0], w[12], pars->GC_7, amp[417]); 
  FFV1_0(w[103], w[0], w[12], pars->GC_7, amp[418]); 
  FFV2_4_0(w[101], w[0], w[13], pars->GC_58, pars->GC_69, amp[419]); 
  FFV2_4_0(w[102], w[0], w[13], pars->GC_58, pars->GC_69, amp[420]); 
  FFV2_4_0(w[103], w[0], w[13], pars->GC_58, pars->GC_69, amp[421]); 
  FFV5_0(w[3], w[104], w[7], pars->GC_2, amp[422]); 
  FFV5_0(w[3], w[105], w[7], pars->GC_2, amp[423]); 
  FFV5_0(w[3], w[106], w[7], pars->GC_2, amp[424]); 
  FFV5_0(w[3], w[104], w[12], pars->GC_7, amp[425]); 
  FFV5_0(w[3], w[105], w[12], pars->GC_7, amp[426]); 
  FFV5_0(w[3], w[106], w[12], pars->GC_7, amp[427]); 
  FFV2_7_0(w[3], w[104], w[13], pars->GC_59, pars->GC_69, amp[428]); 
  FFV2_7_0(w[3], w[105], w[13], pars->GC_59, pars->GC_69, amp[429]); 
  FFV2_7_0(w[3], w[106], w[13], pars->GC_59, pars->GC_69, amp[430]); 
  FFV5_0(w[3], w[107], w[54], pars->GC_2, amp[431]); 
  FFV5_0(w[3], w[108], w[54], pars->GC_2, amp[432]); 
  FFV5_0(w[3], w[109], w[54], pars->GC_2, amp[433]); 
  FFV5_0(w[3], w[107], w[59], pars->GC_7, amp[434]); 
  FFV5_0(w[3], w[108], w[59], pars->GC_7, amp[435]); 
  FFV5_0(w[3], w[109], w[59], pars->GC_7, amp[436]); 
  FFV2_7_0(w[3], w[107], w[60], pars->GC_59, pars->GC_69, amp[437]); 
  FFV2_7_0(w[3], w[108], w[60], pars->GC_59, pars->GC_69, amp[438]); 
  FFV2_7_0(w[3], w[109], w[60], pars->GC_59, pars->GC_69, amp[439]); 
  FFV5_0(w[110], w[2], w[7], pars->GC_2, amp[440]); 
  FFV5_0(w[111], w[2], w[7], pars->GC_2, amp[441]); 
  FFV5_0(w[112], w[2], w[7], pars->GC_2, amp[442]); 
  FFV5_0(w[110], w[2], w[12], pars->GC_7, amp[443]); 
  FFV5_0(w[111], w[2], w[12], pars->GC_7, amp[444]); 
  FFV5_0(w[112], w[2], w[12], pars->GC_7, amp[445]); 
  FFV2_7_0(w[110], w[2], w[13], pars->GC_59, pars->GC_69, amp[446]); 
  FFV2_7_0(w[111], w[2], w[13], pars->GC_59, pars->GC_69, amp[447]); 
  FFV2_7_0(w[112], w[2], w[13], pars->GC_59, pars->GC_69, amp[448]); 
  FFV5_0(w[113], w[2], w[54], pars->GC_2, amp[449]); 
  FFV5_0(w[114], w[2], w[54], pars->GC_2, amp[450]); 
  FFV5_0(w[115], w[2], w[54], pars->GC_2, amp[451]); 
  FFV5_0(w[113], w[2], w[59], pars->GC_7, amp[452]); 
  FFV5_0(w[114], w[2], w[59], pars->GC_7, amp[453]); 
  FFV5_0(w[115], w[2], w[59], pars->GC_7, amp[454]); 
  FFV2_7_0(w[113], w[2], w[60], pars->GC_59, pars->GC_69, amp[455]); 
  FFV2_7_0(w[114], w[2], w[60], pars->GC_59, pars->GC_69, amp[456]); 
  FFV2_7_0(w[115], w[2], w[60], pars->GC_59, pars->GC_69, amp[457]); 

}
double CPPProcess::matrix_3_dxdx_ttxdxdx() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 458; 
  const int ncolor = 6; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1, 1, 1, 1, 1, 1}; 
  static const double cf[ncolor][ncolor] = {{27, 9, 9, 3, 3, 9}, {9, 27, 3, 9,
      9, 3}, {9, 3, 27, 9, 9, 3}, {3, 9, 9, 27, 3, 9}, {3, 9, 9, 3, 27, 9}, {9,
      3, 3, 9, 9, 27}};

  // Calculate color flows
  jamp[0] = +1./4. * std::complex<double> (0, 1) * amp[26] + 1./4. * amp[30] +
      1./4. * amp[31] + 1./4. * amp[32] + 1./4. * std::complex<double> (0, 1) *
      amp[68] + 1./4. * std::complex<double> (0, 1) * amp[69] + 1./4. *
      std::complex<double> (0, 1) * amp[70] + 1./4. * amp[92] + 1./4. * amp[94]
      - 1./2. * amp[118] - 1./4. * amp[120] - 1./2. * amp[127] - 1./4. *
      amp[129] - 1./2. * amp[157] - 1./2. * amp[158] - 1./2. * amp[159] - 1./2.
      * amp[160] - 1./2. * amp[161] - 1./2. * amp[162] + 1./12. * amp[178] +
      1./12. * amp[179] + 1./12. * amp[180] + 1./12. * amp[181] + 1./12. *
      amp[182] + 1./12. * amp[183] - 1./2. * amp[200] - 1./2. * amp[201] -
      1./2. * amp[202] - 1./2. * amp[203] - 1./2. * amp[204] - 1./2. * amp[205]
      - 1./2. * amp[206] - 1./2. * amp[207] - 1./2. * amp[228] - 1./2. *
      amp[229] - 1./2. * amp[230] - 1./2. * amp[231] + 1./12. * amp[242] +
      1./12. * amp[243] + 1./12. * amp[244] + 1./12. * amp[245] - 1./2. *
      amp[256] - 1./2. * amp[257] - 1./2. * amp[258] - 1./2. * amp[259] +
      amp[266] + 1./2. * amp[268] - 1./6. * amp[269] - 1./12. * amp[271] +
      amp[272] + 1./2. * amp[274] + amp[275] + 1./2. * amp[277] - 1./6. *
      amp[278] - 1./12. * amp[280] + amp[281] + 1./2. * amp[283] + amp[284] +
      1./2. * amp[286] - 1./6. * amp[287] - 1./12. * amp[289] + amp[290] +
      1./2. * amp[292] + amp[293] + 1./2. * amp[295] - 1./6. * amp[296] -
      1./12. * amp[298] + amp[299] + 1./2. * amp[301] + 1./4. * amp[371] +
      1./4. * amp[373] - 1./2. * amp[398] - 1./4. * amp[400] + 1./2. * amp[435]
      - 1./12. * amp[436] - 1./2. * amp[443] - 1./4. * amp[445] + 1./2. *
      amp[453] - 1./12. * amp[454];
  jamp[1] = +1./2. * amp[20] + 1./2. * amp[21] + 1./2. * amp[22] + 1./2. *
      amp[23] + 1./2. * amp[24] + 1./2. * amp[25] - 1./12. * amp[27] - 1./12. *
      amp[28] - 1./12. * amp[29] - 1./12. * amp[30] - 1./12. * amp[31] - 1./12.
      * amp[32] + 1./2. * amp[33] + 1./2. * amp[34] + 1./2. * amp[35] + 1./2. *
      amp[36] + 1./2. * amp[37] + 1./2. * amp[38] + 1./2. * amp[39] + 1./2. *
      amp[40] - 1./2. * amp[128] + 1./12. * amp[129] - 1./2. * amp[137] +
      1./12. * amp[138] + 1./4. * std::complex<double> (0, 1) * amp[177] -
      1./4. * amp[178] - 1./4. * amp[179] - 1./4. * amp[180] + 1./4. *
      std::complex<double> (0, 1) * amp[219] + 1./4. * std::complex<double> (0,
      1) * amp[220] + 1./4. * std::complex<double> (0, 1) * amp[221] - 1./4. *
      amp[242] - 1./4. * amp[244] + 1./2. * amp[287] + 1./4. * amp[289] + 1./2.
      * amp[296] + 1./4. * amp[298] - 1./4. * amp[321] - 1./4. * amp[323] +
      1./2. * amp[357] + 1./2. * amp[358] + 1./2. * amp[359] + 1./2. * amp[360]
      - 1./12. * amp[371] - 1./12. * amp[372] - 1./12. * amp[373] - 1./12. *
      amp[374] + 1./2. * amp[385] + 1./2. * amp[386] + 1./2. * amp[387] + 1./2.
      * amp[388] - amp[395] - 1./2. * amp[397] + 1./6. * amp[398] + 1./12. *
      amp[400] - amp[401] - 1./2. * amp[403] + 1./2. * amp[407] + 1./4. *
      amp[409] - amp[413] - 1./2. * amp[415] + 1./6. * amp[416] + 1./12. *
      amp[418] - amp[419] - 1./2. * amp[421] - amp[422] - 1./2. * amp[424] +
      1./6. * amp[425] + 1./12. * amp[427] - amp[428] - 1./2. * amp[430] +
      1./2. * amp[434] + 1./4. * amp[436] - amp[440] - 1./2. * amp[442] + 1./6.
      * amp[443] + 1./12. * amp[445] - amp[446] - 1./2. * amp[448];
  jamp[2] = +1./2. * amp[6] + 1./2. * amp[7] + 1./2. * amp[8] + 1./2. * amp[9]
      + 1./2. * amp[10] + 1./2. * amp[11] - 1./12. * amp[27] - 1./12. * amp[28]
      - 1./12. * amp[29] - 1./12. * amp[30] - 1./12. * amp[31] - 1./12. *
      amp[32] + 1./2. * amp[49] + 1./2. * amp[50] + 1./2. * amp[51] + 1./2. *
      amp[52] + 1./2. * amp[53] + 1./2. * amp[54] + 1./2. * amp[55] + 1./2. *
      amp[56] + 1./2. * amp[77] + 1./2. * amp[78] + 1./2. * amp[79] + 1./2. *
      amp[80] - 1./12. * amp[91] - 1./12. * amp[92] - 1./12. * amp[93] - 1./12.
      * amp[94] + 1./2. * amp[105] + 1./2. * amp[106] + 1./2. * amp[107] +
      1./2. * amp[108] - amp[115] - 1./2. * amp[117] + 1./6. * amp[118] +
      1./12. * amp[120] - amp[121] - 1./2. * amp[123] - amp[124] - 1./2. *
      amp[126] + 1./6. * amp[127] + 1./12. * amp[129] - amp[130] - 1./2. *
      amp[132] - amp[133] - 1./2. * amp[135] + 1./6. * amp[136] + 1./12. *
      amp[138] - amp[139] - 1./2. * amp[141] - amp[142] - 1./2. * amp[144] +
      1./6. * amp[145] + 1./12. * amp[147] - amp[148] - 1./2. * amp[150] -
      1./4. * std::complex<double> (0, 1) * amp[177] - 1./4. * amp[181] - 1./4.
      * amp[182] - 1./4. * amp[183] - 1./4. * std::complex<double> (0, 1) *
      amp[219] - 1./4. * std::complex<double> (0, 1) * amp[220] - 1./4. *
      std::complex<double> (0, 1) * amp[221] - 1./4. * amp[243] - 1./4. *
      amp[245] + 1./2. * amp[269] + 1./4. * amp[271] + 1./2. * amp[278] + 1./4.
      * amp[280] - 1./4. * amp[320] - 1./4. * amp[322] + 1./2. * amp[347] +
      1./4. * amp[349] - 1./2. * amp[426] + 1./12. * amp[427] - 1./2. *
      amp[444] + 1./12. * amp[445] + 1./2. * amp[452] + 1./4. * amp[454];
  jamp[3] = +amp[0] + amp[1] + amp[2] + amp[3] + amp[4] + amp[5] - 1./6. *
      amp[6] - 1./6. * amp[7] - 1./6. * amp[8] - 1./6. * amp[9] - 1./6. *
      amp[10] - 1./6. * amp[11] + amp[12] + amp[13] + amp[14] + amp[15] +
      amp[16] + amp[17] + amp[18] + amp[19] - 1./6. * amp[20] - 1./6. * amp[21]
      - 1./6. * amp[22] - 1./6. * amp[23] - 1./6. * amp[24] - 1./6. * amp[25] +
      1./36. * amp[27] + 1./36. * amp[28] + 1./36. * amp[29] + 1./36. * amp[30]
      + 1./36. * amp[31] + 1./36. * amp[32] - 1./6. * amp[33] - 1./6. * amp[34]
      - 1./6. * amp[35] - 1./6. * amp[36] - 1./6. * amp[37] - 1./6. * amp[38] -
      1./6. * amp[39] - 1./6. * amp[40] + amp[41] + amp[42] + amp[43] + amp[44]
      + amp[45] + amp[46] + amp[47] + amp[48] - 1./6. * amp[49] - 1./6. *
      amp[50] - 1./6. * amp[51] - 1./6. * amp[52] - 1./6. * amp[53] - 1./6. *
      amp[54] - 1./6. * amp[55] - 1./6. * amp[56] + amp[57] + amp[58] + amp[59]
      + amp[60] + amp[61] + amp[62] + amp[63] + amp[64] + amp[65] + amp[66] -
      1./6. * amp[67] + amp[71] + amp[72] + amp[73] + amp[74] + amp[75] +
      amp[76] - 1./6. * amp[77] - 1./6. * amp[78] - 1./6. * amp[79] - 1./6. *
      amp[80] + amp[81] + amp[82] + amp[83] + amp[84] + amp[85] + amp[86] -
      1./6. * amp[87] - 1./6. * amp[88] - 1./6. * amp[89] - 1./6. * amp[90] +
      1./36. * amp[91] + 1./36. * amp[92] + 1./36. * amp[93] + 1./36. * amp[94]
      - 1./6. * amp[95] - 1./6. * amp[96] - 1./6. * amp[97] - 1./6. * amp[98] -
      1./6. * amp[99] - 1./6. * amp[100] + amp[101] + amp[102] + amp[103] +
      amp[104] - 1./6. * amp[105] - 1./6. * amp[106] - 1./6. * amp[107] - 1./6.
      * amp[108] + amp[109] + amp[110] + amp[111] + amp[112] + amp[113] +
      amp[114] - amp[116] + 1./6. * amp[117] + 1./6. * amp[119] - 1./36. *
      amp[120] - amp[122] + 1./6. * amp[123] - amp[125] + 1./6. * amp[126] +
      1./6. * amp[128] - 1./36. * amp[129] - amp[131] + 1./6. * amp[132] -
      amp[134] + 1./6. * amp[135] + 1./6. * amp[137] - 1./36. * amp[138] -
      amp[140] + 1./6. * amp[141] - amp[143] + 1./6. * amp[144] + 1./6. *
      amp[146] - 1./36. * amp[147] - amp[149] + 1./6. * amp[150] - 1./2. *
      amp[218] - 1./2. * amp[238] - 1./2. * amp[239] - 1./2. * amp[240] - 1./2.
      * amp[241] + 1./12. * amp[242] + 1./12. * amp[243] + 1./12. * amp[244] +
      1./12. * amp[245] - 1./2. * amp[246] - 1./2. * amp[247] - 1./2. *
      amp[248] - 1./2. * amp[249] - 1./2. * amp[250] - 1./2. * amp[251] + 1./2.
      * amp[270] - 1./12. * amp[271] + 1./2. * amp[297] - 1./12. * amp[298] -
      1./2. * amp[316] - 1./2. * amp[317] - 1./2. * amp[318] - 1./2. * amp[319]
      + 1./12. * amp[320] + 1./12. * amp[321] + 1./12. * amp[322] + 1./12. *
      amp[323] - 1./2. * amp[324] - 1./2. * amp[325] - 1./2. * amp[326] - 1./2.
      * amp[327] - 1./2. * amp[328] - 1./2. * amp[329] + 1./2. * amp[348] -
      1./12. * amp[349] + amp[353] + amp[354] + amp[355] + amp[356] - 1./6. *
      amp[357] - 1./6. * amp[358] - 1./6. * amp[359] - 1./6. * amp[360] +
      amp[361] + amp[362] + amp[363] + amp[364] + amp[365] + amp[366] - 1./6. *
      amp[367] - 1./6. * amp[368] - 1./6. * amp[369] - 1./6. * amp[370] +
      1./36. * amp[371] + 1./36. * amp[372] + 1./36. * amp[373] + 1./36. *
      amp[374] - 1./6. * amp[375] - 1./6. * amp[376] - 1./6. * amp[377] - 1./6.
      * amp[378] - 1./6. * amp[379] - 1./6. * amp[380] + amp[381] + amp[382] +
      amp[383] + amp[384] - 1./6. * amp[385] - 1./6. * amp[386] - 1./6. *
      amp[387] - 1./6. * amp[388] + amp[389] + amp[390] + amp[391] + amp[392] +
      amp[393] + amp[394] - amp[396] + 1./6. * amp[397] + 1./6. * amp[399] -
      1./36. * amp[400] - amp[402] + 1./6. * amp[403] + 1./2. * amp[408] -
      1./12. * amp[409] - amp[414] + 1./6. * amp[415] + 1./6. * amp[417] -
      1./36. * amp[418] - amp[420] + 1./6. * amp[421] - amp[423] + 1./6. *
      amp[424] + 1./6. * amp[426] - 1./36. * amp[427] - amp[429] + 1./6. *
      amp[430] - amp[441] + 1./6. * amp[442] + 1./6. * amp[444] - 1./36. *
      amp[445] - amp[447] + 1./6. * amp[448];
  jamp[4] = -1./4. * std::complex<double> (0, 1) * amp[26] + 1./4. * amp[27] +
      1./4. * amp[28] + 1./4. * amp[29] - 1./4. * std::complex<double> (0, 1) *
      amp[68] - 1./4. * std::complex<double> (0, 1) * amp[69] - 1./4. *
      std::complex<double> (0, 1) * amp[70] + 1./4. * amp[91] + 1./4. * amp[93]
      - 1./2. * amp[136] - 1./4. * amp[138] - 1./2. * amp[145] - 1./4. *
      amp[147] - 1./2. * amp[171] - 1./2. * amp[172] - 1./2. * amp[173] - 1./2.
      * amp[174] - 1./2. * amp[175] - 1./2. * amp[176] + 1./12. * amp[178] +
      1./12. * amp[179] + 1./12. * amp[180] + 1./12. * amp[181] + 1./12. *
      amp[182] + 1./12. * amp[183] - 1./2. * amp[184] - 1./2. * amp[185] -
      1./2. * amp[186] - 1./2. * amp[187] - 1./2. * amp[188] - 1./2. * amp[189]
      - 1./2. * amp[190] - 1./2. * amp[191] + 1./2. * amp[279] - 1./12. *
      amp[280] + 1./2. * amp[288] - 1./12. * amp[289] - 1./2. * amp[306] -
      1./2. * amp[307] - 1./2. * amp[308] - 1./2. * amp[309] + 1./12. *
      amp[320] + 1./12. * amp[321] + 1./12. * amp[322] + 1./12. * amp[323] -
      1./2. * amp[334] - 1./2. * amp[335] - 1./2. * amp[336] - 1./2. * amp[337]
      + amp[344] + 1./2. * amp[346] - 1./6. * amp[347] - 1./12. * amp[349] +
      amp[350] + 1./2. * amp[352] + 1./4. * amp[372] + 1./4. * amp[374] +
      amp[404] + 1./2. * amp[406] - 1./6. * amp[407] - 1./12. * amp[409] +
      amp[410] + 1./2. * amp[412] - 1./2. * amp[416] - 1./4. * amp[418] - 1./2.
      * amp[425] - 1./4. * amp[427] + amp[431] + 1./2. * amp[433] - 1./6. *
      amp[434] - 1./12. * amp[436] + amp[437] + 1./2. * amp[439] + amp[449] +
      1./2. * amp[451] - 1./6. * amp[452] - 1./12. * amp[454] + amp[455] +
      1./2. * amp[457];
  jamp[5] = +1./2. * amp[67] + 1./2. * amp[87] + 1./2. * amp[88] + 1./2. *
      amp[89] + 1./2. * amp[90] - 1./12. * amp[91] - 1./12. * amp[92] - 1./12.
      * amp[93] - 1./12. * amp[94] + 1./2. * amp[95] + 1./2. * amp[96] + 1./2.
      * amp[97] + 1./2. * amp[98] + 1./2. * amp[99] + 1./2. * amp[100] - 1./2.
      * amp[119] + 1./12. * amp[120] - 1./2. * amp[146] + 1./12. * amp[147] -
      amp[151] - amp[152] - amp[153] - amp[154] - amp[155] - amp[156] + 1./6. *
      amp[157] + 1./6. * amp[158] + 1./6. * amp[159] + 1./6. * amp[160] + 1./6.
      * amp[161] + 1./6. * amp[162] - amp[163] - amp[164] - amp[165] - amp[166]
      - amp[167] - amp[168] - amp[169] - amp[170] + 1./6. * amp[171] + 1./6. *
      amp[172] + 1./6. * amp[173] + 1./6. * amp[174] + 1./6. * amp[175] + 1./6.
      * amp[176] - 1./36. * amp[178] - 1./36. * amp[179] - 1./36. * amp[180] -
      1./36. * amp[181] - 1./36. * amp[182] - 1./36. * amp[183] + 1./6. *
      amp[184] + 1./6. * amp[185] + 1./6. * amp[186] + 1./6. * amp[187] + 1./6.
      * amp[188] + 1./6. * amp[189] + 1./6. * amp[190] + 1./6. * amp[191] -
      amp[192] - amp[193] - amp[194] - amp[195] - amp[196] - amp[197] -
      amp[198] - amp[199] + 1./6. * amp[200] + 1./6. * amp[201] + 1./6. *
      amp[202] + 1./6. * amp[203] + 1./6. * amp[204] + 1./6. * amp[205] + 1./6.
      * amp[206] + 1./6. * amp[207] - amp[208] - amp[209] - amp[210] - amp[211]
      - amp[212] - amp[213] - amp[214] - amp[215] - amp[216] - amp[217] + 1./6.
      * amp[218] - amp[222] - amp[223] - amp[224] - amp[225] - amp[226] -
      amp[227] + 1./6. * amp[228] + 1./6. * amp[229] + 1./6. * amp[230] + 1./6.
      * amp[231] - amp[232] - amp[233] - amp[234] - amp[235] - amp[236] -
      amp[237] + 1./6. * amp[238] + 1./6. * amp[239] + 1./6. * amp[240] + 1./6.
      * amp[241] - 1./36. * amp[242] - 1./36. * amp[243] - 1./36. * amp[244] -
      1./36. * amp[245] + 1./6. * amp[246] + 1./6. * amp[247] + 1./6. *
      amp[248] + 1./6. * amp[249] + 1./6. * amp[250] + 1./6. * amp[251] -
      amp[252] - amp[253] - amp[254] - amp[255] + 1./6. * amp[256] + 1./6. *
      amp[257] + 1./6. * amp[258] + 1./6. * amp[259] - amp[260] - amp[261] -
      amp[262] - amp[263] - amp[264] - amp[265] + amp[267] - 1./6. * amp[268] -
      1./6. * amp[270] + 1./36. * amp[271] + amp[273] - 1./6. * amp[274] +
      amp[276] - 1./6. * amp[277] - 1./6. * amp[279] + 1./36. * amp[280] +
      amp[282] - 1./6. * amp[283] + amp[285] - 1./6. * amp[286] - 1./6. *
      amp[288] + 1./36. * amp[289] + amp[291] - 1./6. * amp[292] + amp[294] -
      1./6. * amp[295] - 1./6. * amp[297] + 1./36. * amp[298] + amp[300] -
      1./6. * amp[301] - amp[302] - amp[303] - amp[304] - amp[305] + 1./6. *
      amp[306] + 1./6. * amp[307] + 1./6. * amp[308] + 1./6. * amp[309] -
      amp[310] - amp[311] - amp[312] - amp[313] - amp[314] - amp[315] + 1./6. *
      amp[316] + 1./6. * amp[317] + 1./6. * amp[318] + 1./6. * amp[319] -
      1./36. * amp[320] - 1./36. * amp[321] - 1./36. * amp[322] - 1./36. *
      amp[323] + 1./6. * amp[324] + 1./6. * amp[325] + 1./6. * amp[326] + 1./6.
      * amp[327] + 1./6. * amp[328] + 1./6. * amp[329] - amp[330] - amp[331] -
      amp[332] - amp[333] + 1./6. * amp[334] + 1./6. * amp[335] + 1./6. *
      amp[336] + 1./6. * amp[337] - amp[338] - amp[339] - amp[340] - amp[341] -
      amp[342] - amp[343] + amp[345] - 1./6. * amp[346] - 1./6. * amp[348] +
      1./36. * amp[349] + amp[351] - 1./6. * amp[352] + 1./2. * amp[367] +
      1./2. * amp[368] + 1./2. * amp[369] + 1./2. * amp[370] - 1./12. *
      amp[371] - 1./12. * amp[372] - 1./12. * amp[373] - 1./12. * amp[374] +
      1./2. * amp[375] + 1./2. * amp[376] + 1./2. * amp[377] + 1./2. * amp[378]
      + 1./2. * amp[379] + 1./2. * amp[380] - 1./2. * amp[399] + 1./12. *
      amp[400] + amp[405] - 1./6. * amp[406] - 1./6. * amp[408] + 1./36. *
      amp[409] + amp[411] - 1./6. * amp[412] - 1./2. * amp[417] + 1./12. *
      amp[418] + amp[432] - 1./6. * amp[433] - 1./6. * amp[435] + 1./36. *
      amp[436] + amp[438] - 1./6. * amp[439] + amp[450] - 1./6. * amp[451] -
      1./6. * amp[453] + 1./36. * amp[454] + amp[456] - 1./6. * amp[457];

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



