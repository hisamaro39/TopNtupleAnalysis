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
// Process: g g > t t~ g g NP<=2 @3

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
  jamp2[0] = new double[50]; 
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
  for(int i = 0; i < 50; i++ )
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
  const int denominators[nprocesses] = {512}; 

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
        t[0] = matrix_3_gg_ttxgg(); 

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
      t[0] = matrix_3_gg_ttxgg(); 

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
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  VVV2P0_1(w[0], w[1], pars->GC_28, pars->ZERO, pars->ZERO, w[6]); 
  FFV5P0_3(w[3], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[7]); 
  VVV1P0_1(w[6], w[4], pars->GC_6, pars->ZERO, pars->ZERO, w[8]); 
  VVV1P0_1(w[6], w[5], pars->GC_6, pars->ZERO, pars->ZERO, w[9]); 
  VVV1P0_1(w[0], w[1], pars->GC_6, pars->ZERO, pars->ZERO, w[10]); 
  FFS2_3(w[3], w[2], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[11]); 
  VVV1P0_1(w[10], w[4], pars->GC_6, pars->ZERO, pars->ZERO, w[12]); 
  VVV1P0_1(w[10], w[5], pars->GC_6, pars->ZERO, pars->ZERO, w[13]); 
  VVV2P0_1(w[10], w[4], pars->GC_28, pars->ZERO, pars->ZERO, w[14]); 
  VVV2P0_1(w[10], w[5], pars->GC_28, pars->ZERO, pars->ZERO, w[15]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[16]);
  VVV1P0_1(w[4], w[5], pars->GC_6, pars->ZERO, pars->ZERO, w[17]); 
  VVV2P0_1(w[4], w[5], pars->GC_28, pars->ZERO, pars->ZERO, w[18]); 
  VVS2_3(w[0], w[1], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[19]); 
  FFV5_1(w[2], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[20]); 
  FFS2_2(w[3], w[19], pars->GC_95, pars->mdl_MT, pars->mdl_WT, w[21]); 
  FFV5_2(w[3], w[6], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[22]); 
  FFV5_2(w[3], w[10], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[23]); 
  FFV3_8_2(w[3], w[10], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[24]);
  VVS2_3(w[10], w[5], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[25]); 
  FFV3_8_1(w[2], w[4], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[26]);
  FFV5_2(w[3], w[5], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[27]); 
  FFV3_8_2(w[3], w[5], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[28]);
  FFV5_1(w[2], w[5], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[29]); 
  VVS2_3(w[10], w[4], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[30]); 
  FFV3_8_1(w[2], w[5], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[31]);
  FFV5_2(w[3], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[32]); 
  FFV3_8_2(w[3], w[4], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[33]);
  FFS2_1(w[2], w[19], pars->GC_95, pars->mdl_MT, pars->mdl_WT, w[34]); 
  FFV5_1(w[2], w[6], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[35]); 
  FFV5_1(w[2], w[10], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[36]); 
  FFV3_8_1(w[2], w[10], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[37]);
  VVS2_3(w[4], w[5], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[38]); 
  FFVV1_2P0_3(w[3], w[2], w[4], pars->GC_146, pars->GC_85, pars->ZERO,
      pars->ZERO, w[39]);
  FFVV1_2P0_3(w[3], w[2], w[5], pars->GC_146, pars->GC_85, pars->ZERO,
      pars->ZERO, w[40]);
  FFVV1_2_1(w[2], w[4], w[5], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[41]);
  FFVV1_2_2(w[3], w[4], w[5], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[42]);
  FFV5_1(w[2], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[43]); 
  FFV5_2(w[3], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[44]); 
  FFV5_1(w[43], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[45]); 
  FFV3_8_1(w[43], w[4], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[46]);
  FFV5_1(w[43], w[5], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[47]); 
  FFV3_8_1(w[43], w[5], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[48]);
  FFV3_8_2(w[3], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[49]);
  FFV3_8_1(w[2], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[50]);
  FFV5_1(w[50], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[51]); 
  FFV5_1(w[50], w[5], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[52]); 
  VVS2_3(w[1], w[4], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[53]); 
  VVV2P0_1(w[1], w[4], pars->GC_28, pars->ZERO, pars->ZERO, w[54]); 
  FFV5P0_3(w[3], w[43], pars->GC_7, pars->ZERO, pars->ZERO, w[55]); 
  VVV1P0_1(w[1], w[4], pars->GC_6, pars->ZERO, pars->ZERO, w[56]); 
  FFS2_3(w[3], w[43], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[57]); 
  FFV3_8P0_3(w[3], w[43], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[58]);
  FFV5P0_3(w[3], w[50], pars->GC_7, pars->ZERO, pars->ZERO, w[59]); 
  VVS2_3(w[1], w[5], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[60]); 
  VVV2P0_1(w[1], w[5], pars->GC_28, pars->ZERO, pars->ZERO, w[61]); 
  VVV1P0_1(w[1], w[5], pars->GC_6, pars->ZERO, pars->ZERO, w[62]); 
  FFV5_1(w[43], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[63]); 
  FFV3_8_1(w[43], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[64]);
  FFV5_1(w[50], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[65]); 
  FFVV1_2_2(w[3], w[1], w[4], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[66]);
  FFVV1_2_2(w[3], w[1], w[5], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[67]);
  VVVV2P0_1(w[1], w[4], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[68]); 
  VVVV8P0_1(w[1], w[4], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[69]); 
  VVVV7P0_1(w[1], w[4], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[70]); 
  VVVV1P0_1(w[1], w[4], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[71]); 
  VVVV4P0_1(w[1], w[4], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[72]); 
  VVVV5P0_1(w[1], w[4], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[73]); 
  VVVS1_4(w[1], w[4], w[5], pars->GC_84, pars->mdl_MH, pars->mdl_WH, w[74]); 
  FFV5_2(w[3], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[75]); 
  FFV5_1(w[2], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[76]); 
  FFV5_2(w[75], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[77]); 
  FFV3_8_2(w[75], w[4], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[78]);
  FFV5_2(w[75], w[5], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[79]); 
  FFV3_8_2(w[75], w[5], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[80]);
  FFV3_8_1(w[2], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[81]);
  FFV3_8_2(w[3], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[82]);
  FFV5_2(w[82], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[83]); 
  FFV5_2(w[82], w[5], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[84]); 
  FFV5P0_3(w[75], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[85]); 
  FFS2_3(w[75], w[2], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[86]); 
  FFV3_8P0_3(w[75], w[2], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[87]);
  FFV5P0_3(w[82], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[88]); 
  FFV5_2(w[75], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[89]); 
  FFV3_8_2(w[75], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[90]);
  FFV5_2(w[82], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[91]); 
  FFVV1_2_1(w[2], w[1], w[4], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[92]);
  FFVV1_2_1(w[2], w[1], w[5], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[93]);
  VVS2_3(w[0], w[4], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[94]); 
  FFS2_2(w[3], w[94], pars->GC_95, pars->mdl_MT, pars->mdl_WT, w[95]); 
  VVV2P0_1(w[0], w[4], pars->GC_28, pars->ZERO, pars->ZERO, w[96]); 
  FFV5_2(w[3], w[96], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[97]); 
  VVV1P0_1(w[96], w[5], pars->GC_6, pars->ZERO, pars->ZERO, w[98]); 
  VVV1P0_1(w[0], w[4], pars->GC_6, pars->ZERO, pars->ZERO, w[99]); 
  FFV5_2(w[3], w[99], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[100]); 
  FFV3_8_2(w[3], w[99], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[101]);
  VVS2_3(w[99], w[5], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[102]); 
  VVV2P0_1(w[99], w[5], pars->GC_28, pars->ZERO, pars->ZERO, w[103]); 
  VVV1P0_1(w[99], w[5], pars->GC_6, pars->ZERO, pars->ZERO, w[104]); 
  FFS2_1(w[2], w[94], pars->GC_95, pars->mdl_MT, pars->mdl_WT, w[105]); 
  FFV5_1(w[2], w[96], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[106]); 
  FFV5_1(w[2], w[99], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[107]); 
  FFV3_8_1(w[2], w[99], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[108]);
  VVV1P0_1(w[96], w[1], pars->GC_6, pars->ZERO, pars->ZERO, w[109]); 
  VVV1P0_1(w[99], w[1], pars->GC_6, pars->ZERO, pars->ZERO, w[110]); 
  VVV2P0_1(w[99], w[1], pars->GC_28, pars->ZERO, pars->ZERO, w[111]); 
  VVS2_3(w[99], w[1], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[112]); 
  FFVV1_2P0_3(w[3], w[2], w[1], pars->GC_146, pars->GC_85, pars->ZERO,
      pars->ZERO, w[113]);
  VVS2_3(w[0], w[5], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[114]); 
  FFS2_2(w[3], w[114], pars->GC_95, pars->mdl_MT, pars->mdl_WT, w[115]); 
  VVV2P0_1(w[0], w[5], pars->GC_28, pars->ZERO, pars->ZERO, w[116]); 
  FFV5_2(w[3], w[116], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[117]); 
  VVV1P0_1(w[116], w[4], pars->GC_6, pars->ZERO, pars->ZERO, w[118]); 
  VVV1P0_1(w[0], w[5], pars->GC_6, pars->ZERO, pars->ZERO, w[119]); 
  FFV5_2(w[3], w[119], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[120]); 
  FFV3_8_2(w[3], w[119], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[121]);
  VVS2_3(w[119], w[4], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[122]); 
  VVV2P0_1(w[119], w[4], pars->GC_28, pars->ZERO, pars->ZERO, w[123]); 
  VVV1P0_1(w[119], w[4], pars->GC_6, pars->ZERO, pars->ZERO, w[124]); 
  FFS2_1(w[2], w[114], pars->GC_95, pars->mdl_MT, pars->mdl_WT, w[125]); 
  FFV5_1(w[2], w[116], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[126]); 
  FFV5_1(w[2], w[119], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[127]); 
  FFV3_8_1(w[2], w[119], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[128]);
  VVV1P0_1(w[116], w[1], pars->GC_6, pars->ZERO, pars->ZERO, w[129]); 
  VVV1P0_1(w[119], w[1], pars->GC_6, pars->ZERO, pars->ZERO, w[130]); 
  VVV2P0_1(w[119], w[1], pars->GC_28, pars->ZERO, pars->ZERO, w[131]); 
  VVS2_3(w[119], w[1], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[132]); 
  FFV5_1(w[76], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[133]); 
  FFV3_8_1(w[76], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[134]);
  FFV5_2(w[32], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[135]); 
  FFV3_8_2(w[32], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[136]);
  FFV5_2(w[33], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[137]); 
  FFV5_1(w[81], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[138]); 
  FFV5_2(w[27], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[139]); 
  FFV3_8_2(w[27], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[140]);
  FFV5_2(w[28], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[141]); 
  VVV1P0_1(w[0], w[18], pars->GC_6, pars->ZERO, pars->ZERO, w[142]); 
  VVS2_3(w[0], w[17], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[143]); 
  VVV2P0_1(w[0], w[17], pars->GC_28, pars->ZERO, pars->ZERO, w[144]); 
  VVV1P0_1(w[0], w[17], pars->GC_6, pars->ZERO, pars->ZERO, w[145]); 
  FFV5_2(w[44], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[146]); 
  FFV3_8_2(w[44], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[147]);
  FFV5_1(w[20], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[148]); 
  FFV3_8_1(w[20], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[149]);
  FFV5_1(w[26], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[150]); 
  FFV5_2(w[49], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[151]); 
  FFV5_1(w[29], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[152]); 
  FFV3_8_1(w[29], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[153]);
  FFV5_1(w[31], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[154]); 
  VVV1P0_1(w[0], w[54], pars->GC_6, pars->ZERO, pars->ZERO, w[155]); 
  VVV1P0_1(w[0], w[7], pars->GC_6, pars->ZERO, pars->ZERO, w[156]); 
  VVV1P0_1(w[0], w[56], pars->GC_6, pars->ZERO, pars->ZERO, w[157]); 
  VVS2P0_1(w[0], w[11], pars->GC_78, pars->ZERO, pars->ZERO, w[158]); 
  VVV2P0_1(w[0], w[56], pars->GC_28, pars->ZERO, pars->ZERO, w[159]); 
  VVV2P0_1(w[0], w[7], pars->GC_28, pars->ZERO, pars->ZERO, w[160]); 
  VVV1P0_1(w[0], w[16], pars->GC_6, pars->ZERO, pars->ZERO, w[161]); 
  VVS2_3(w[0], w[56], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[162]); 
  VVV1P0_1(w[0], w[61], pars->GC_6, pars->ZERO, pars->ZERO, w[163]); 
  VVV1P0_1(w[0], w[62], pars->GC_6, pars->ZERO, pars->ZERO, w[164]); 
  VVV2P0_1(w[0], w[62], pars->GC_28, pars->ZERO, pars->ZERO, w[165]); 
  VVS2_3(w[0], w[62], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[166]); 
  FFVV1_2_1(w[2], w[0], w[1], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[167]);
  FFVV1_2_2(w[3], w[0], w[1], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[168]);
  VVVV2P0_1(w[0], w[1], w[4], pars->GC_44, pars->ZERO, pars->ZERO, w[169]); 
  VVVV8P0_1(w[0], w[1], w[4], pars->GC_44, pars->ZERO, pars->ZERO, w[170]); 
  VVVV7P0_1(w[0], w[1], w[4], pars->GC_44, pars->ZERO, pars->ZERO, w[171]); 
  VVVV1P0_1(w[0], w[1], w[4], pars->GC_8, pars->ZERO, pars->ZERO, w[172]); 
  VVVV4P0_1(w[0], w[1], w[4], pars->GC_8, pars->ZERO, pars->ZERO, w[173]); 
  VVVV5P0_1(w[0], w[1], w[4], pars->GC_8, pars->ZERO, pars->ZERO, w[174]); 
  VVVS1_4(w[0], w[1], w[4], pars->GC_84, pars->mdl_MH, pars->mdl_WH, w[175]); 
  VVVV2P0_1(w[0], w[1], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[176]); 
  VVVV8P0_1(w[0], w[1], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[177]); 
  VVVV7P0_1(w[0], w[1], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[178]); 
  VVVV1P0_1(w[0], w[1], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[179]); 
  VVVV4P0_1(w[0], w[1], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[180]); 
  VVVV5P0_1(w[0], w[1], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[181]); 
  VVVS1_4(w[0], w[1], w[5], pars->GC_84, pars->mdl_MH, pars->mdl_WH, w[182]); 
  FFVV1_2P0_3(w[3], w[2], w[0], pars->GC_146, pars->GC_85, pars->ZERO,
      pars->ZERO, w[183]);
  FFVV1_2_1(w[2], w[0], w[4], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[184]);
  FFVV1_2_1(w[2], w[0], w[5], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[185]);
  FFVV1_2_2(w[3], w[0], w[4], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[186]);
  FFVV1_2_2(w[3], w[0], w[5], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[187]);
  VVVV2P0_1(w[0], w[4], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[188]); 
  VVVV8P0_1(w[0], w[4], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[189]); 
  VVVV7P0_1(w[0], w[4], w[5], pars->GC_44, pars->ZERO, pars->ZERO, w[190]); 
  VVVV1P0_1(w[0], w[4], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[191]); 
  VVVV4P0_1(w[0], w[4], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[192]); 
  VVVV5P0_1(w[0], w[4], w[5], pars->GC_8, pars->ZERO, pars->ZERO, w[193]); 
  VVVS1_4(w[0], w[4], w[5], pars->GC_84, pars->mdl_MH, pars->mdl_WH, w[194]); 
  VVVVS1_5(w[0], w[1], w[4], w[5], pars->GC_86, pars->mdl_MH, pars->mdl_WH,
      w[195]);
  VVVVS2_5(w[0], w[1], w[4], w[5], pars->GC_86, pars->mdl_MH, pars->mdl_WH,
      w[196]);
  VVVVS3_5(w[0], w[1], w[4], w[5], pars->GC_86, pars->mdl_MH, pars->mdl_WH,
      w[197]);
  VVVVV3P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[198]);
  VVVVV13P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[199]);
  VVVVV15P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[200]);
  VVVVV2P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[201]);
  VVVVV12P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[202]);
  VVVVV14P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[203]);
  VVVVV8P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[204]);
  VVVVV11P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[205]);
  VVVVV7P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[206]);
  VVVVV9P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[207]);
  VVVVV10P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[208]);
  VVVVV6P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[209]);
  VVVVV1P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[210]);
  VVVVV6P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[211]);
  VVVVV7P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[212]);
  VVVVV4P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[213]);
  VVVVV10P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[214]);
  VVVVV14P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[215]);
  VVVVV5P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[216]);
  VVVVV11P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[217]);
  VVVVV12P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[218]);
  VVVVV5P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[219]);
  VVVVV9P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[220]);
  VVVVV15P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[221]);
  VVVVV4P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[222]);
  VVVVV8P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[223]);
  VVVVV13P0_1(w[0], w[1], w[4], w[5], pars->GC_49, pars->ZERO, pars->ZERO,
      w[224]);
  VVVVV1P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[225]);
  VVVVV2P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[226]);
  VVVVV3P0_1(w[0], w[1], w[4], w[5], pars->GC_50, pars->ZERO, pars->ZERO,
      w[227]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVVV1_0(w[6], w[7], w[4], w[5], pars->GC_8, amp[0]); 
  VVVV4_0(w[6], w[7], w[4], w[5], pars->GC_8, amp[1]); 
  VVVV5_0(w[6], w[7], w[4], w[5], pars->GC_8, amp[2]); 
  VVV1_0(w[7], w[5], w[8], pars->GC_6, amp[3]); 
  VVV1_0(w[7], w[4], w[9], pars->GC_6, amp[4]); 
  VVVS1_0(w[10], w[4], w[5], w[11], pars->GC_84, amp[5]); 
  VVS2_0(w[5], w[12], w[11], pars->GC_78, amp[6]); 
  VVS2_0(w[4], w[13], w[11], pars->GC_78, amp[7]); 
  VVVV2_0(w[10], w[7], w[4], w[5], pars->GC_44, amp[8]); 
  VVVV8_0(w[10], w[7], w[4], w[5], pars->GC_44, amp[9]); 
  VVVV7_0(w[10], w[7], w[4], w[5], pars->GC_44, amp[10]); 
  VVVV1_0(w[10], w[7], w[4], w[5], pars->GC_8, amp[11]); 
  VVVV4_0(w[10], w[7], w[4], w[5], pars->GC_8, amp[12]); 
  VVVV5_0(w[10], w[7], w[4], w[5], pars->GC_8, amp[13]); 
  VVV1_0(w[7], w[5], w[14], pars->GC_6, amp[14]); 
  VVV2_0(w[7], w[5], w[12], pars->GC_28, amp[15]); 
  VVV1_0(w[7], w[5], w[12], pars->GC_6, amp[16]); 
  VVV1_0(w[7], w[4], w[15], pars->GC_6, amp[17]); 
  VVV2_0(w[7], w[4], w[13], pars->GC_28, amp[18]); 
  VVV1_0(w[7], w[4], w[13], pars->GC_6, amp[19]); 
  VVVV1_0(w[10], w[16], w[4], w[5], pars->GC_8, amp[20]); 
  VVVV4_0(w[10], w[16], w[4], w[5], pars->GC_8, amp[21]); 
  VVVV5_0(w[10], w[16], w[4], w[5], pars->GC_8, amp[22]); 
  VVV1_0(w[16], w[5], w[12], pars->GC_6, amp[23]); 
  VVV1_0(w[16], w[4], w[13], pars->GC_6, amp[24]); 
  VVV1_0(w[6], w[7], w[17], pars->GC_6, amp[25]); 
  VVS2_0(w[10], w[17], w[11], pars->GC_78, amp[26]); 
  VVV1_0(w[10], w[7], w[18], pars->GC_6, amp[27]); 
  VVV2_0(w[10], w[7], w[17], pars->GC_28, amp[28]); 
  VVV1_0(w[10], w[7], w[17], pars->GC_6, amp[29]); 
  VVV1_0(w[10], w[16], w[17], pars->GC_6, amp[30]); 
  FFV5_0(w[21], w[20], w[5], pars->GC_7, amp[31]); 
  FFV5_0(w[22], w[20], w[5], pars->GC_7, amp[32]); 
  FFV5_0(w[3], w[20], w[9], pars->GC_7, amp[33]); 
  FFVV1_2_0(w[3], w[20], w[10], w[5], pars->GC_146, pars->GC_85, amp[34]); 
  FFV5_0(w[23], w[20], w[5], pars->GC_7, amp[35]); 
  FFV3_8_0(w[23], w[20], w[5], pars->GC_145, pars->GC_79, amp[36]); 
  FFV5_0(w[24], w[20], w[5], pars->GC_7, amp[37]); 
  FFS2_0(w[3], w[20], w[25], pars->GC_95, amp[38]); 
  FFV5_0(w[3], w[20], w[15], pars->GC_7, amp[39]); 
  FFV5_0(w[3], w[20], w[13], pars->GC_7, amp[40]); 
  FFV3_8_0(w[3], w[20], w[13], pars->GC_145, pars->GC_79, amp[41]); 
  FFV5_0(w[23], w[26], w[5], pars->GC_7, amp[42]); 
  FFV5_0(w[3], w[26], w[13], pars->GC_7, amp[43]); 
  FFS2_0(w[27], w[20], w[19], pars->GC_95, amp[44]); 
  FFV5_0(w[27], w[20], w[6], pars->GC_7, amp[45]); 
  FFV5_0(w[27], w[20], w[10], pars->GC_7, amp[46]); 
  FFV3_8_0(w[27], w[20], w[10], pars->GC_145, pars->GC_79, amp[47]); 
  FFV5_0(w[28], w[20], w[10], pars->GC_7, amp[48]); 
  FFV5_0(w[27], w[26], w[10], pars->GC_7, amp[49]); 
  FFV5_0(w[21], w[29], w[4], pars->GC_7, amp[50]); 
  FFV5_0(w[22], w[29], w[4], pars->GC_7, amp[51]); 
  FFV5_0(w[3], w[29], w[8], pars->GC_7, amp[52]); 
  FFVV1_2_0(w[3], w[29], w[10], w[4], pars->GC_146, pars->GC_85, amp[53]); 
  FFV5_0(w[23], w[29], w[4], pars->GC_7, amp[54]); 
  FFV3_8_0(w[23], w[29], w[4], pars->GC_145, pars->GC_79, amp[55]); 
  FFV5_0(w[24], w[29], w[4], pars->GC_7, amp[56]); 
  FFS2_0(w[3], w[29], w[30], pars->GC_95, amp[57]); 
  FFV5_0(w[3], w[29], w[14], pars->GC_7, amp[58]); 
  FFV5_0(w[3], w[29], w[12], pars->GC_7, amp[59]); 
  FFV3_8_0(w[3], w[29], w[12], pars->GC_145, pars->GC_79, amp[60]); 
  FFV5_0(w[23], w[31], w[4], pars->GC_7, amp[61]); 
  FFV5_0(w[3], w[31], w[12], pars->GC_7, amp[62]); 
  FFS2_0(w[32], w[29], w[19], pars->GC_95, amp[63]); 
  FFV5_0(w[32], w[29], w[6], pars->GC_7, amp[64]); 
  FFV5_0(w[32], w[29], w[10], pars->GC_7, amp[65]); 
  FFV3_8_0(w[32], w[29], w[10], pars->GC_145, pars->GC_79, amp[66]); 
  FFV5_0(w[33], w[29], w[10], pars->GC_7, amp[67]); 
  FFV5_0(w[32], w[31], w[10], pars->GC_7, amp[68]); 
  FFV5_0(w[32], w[34], w[5], pars->GC_7, amp[69]); 
  FFV5_0(w[32], w[35], w[5], pars->GC_7, amp[70]); 
  FFV5_0(w[32], w[2], w[9], pars->GC_7, amp[71]); 
  FFVV1_2_0(w[32], w[2], w[10], w[5], pars->GC_146, pars->GC_85, amp[72]); 
  FFV5_0(w[32], w[36], w[5], pars->GC_7, amp[73]); 
  FFV3_8_0(w[32], w[36], w[5], pars->GC_145, pars->GC_79, amp[74]); 
  FFV5_0(w[32], w[37], w[5], pars->GC_7, amp[75]); 
  FFS2_0(w[32], w[2], w[25], pars->GC_95, amp[76]); 
  FFV5_0(w[32], w[2], w[15], pars->GC_7, amp[77]); 
  FFV5_0(w[32], w[2], w[13], pars->GC_7, amp[78]); 
  FFV3_8_0(w[32], w[2], w[13], pars->GC_145, pars->GC_79, amp[79]); 
  FFV5_0(w[33], w[36], w[5], pars->GC_7, amp[80]); 
  FFV5_0(w[33], w[2], w[13], pars->GC_7, amp[81]); 
  FFV5_0(w[27], w[34], w[4], pars->GC_7, amp[82]); 
  FFV5_0(w[27], w[35], w[4], pars->GC_7, amp[83]); 
  FFV5_0(w[27], w[2], w[8], pars->GC_7, amp[84]); 
  FFVV1_2_0(w[27], w[2], w[10], w[4], pars->GC_146, pars->GC_85, amp[85]); 
  FFV5_0(w[27], w[36], w[4], pars->GC_7, amp[86]); 
  FFV3_8_0(w[27], w[36], w[4], pars->GC_145, pars->GC_79, amp[87]); 
  FFV5_0(w[27], w[37], w[4], pars->GC_7, amp[88]); 
  FFS2_0(w[27], w[2], w[30], pars->GC_95, amp[89]); 
  FFV5_0(w[27], w[2], w[14], pars->GC_7, amp[90]); 
  FFV5_0(w[27], w[2], w[12], pars->GC_7, amp[91]); 
  FFV3_8_0(w[27], w[2], w[12], pars->GC_145, pars->GC_79, amp[92]); 
  FFV5_0(w[28], w[36], w[4], pars->GC_7, amp[93]); 
  FFV5_0(w[28], w[2], w[12], pars->GC_7, amp[94]); 
  FFV5_0(w[3], w[34], w[17], pars->GC_7, amp[95]); 
  FFV5_0(w[21], w[2], w[17], pars->GC_7, amp[96]); 
  FFV5_0(w[3], w[35], w[17], pars->GC_7, amp[97]); 
  FFV5_0(w[22], w[2], w[17], pars->GC_7, amp[98]); 
  FFS2_0(w[3], w[36], w[38], pars->GC_95, amp[99]); 
  FFS2_0(w[23], w[2], w[38], pars->GC_95, amp[100]); 
  FFV5_0(w[3], w[36], w[18], pars->GC_7, amp[101]); 
  FFV5_0(w[23], w[2], w[18], pars->GC_7, amp[102]); 
  FFVV1_2_0(w[3], w[2], w[10], w[17], pars->GC_146, pars->GC_85, amp[103]); 
  FFV5_0(w[3], w[36], w[17], pars->GC_7, amp[104]); 
  FFV3_8_0(w[3], w[36], w[17], pars->GC_145, pars->GC_79, amp[105]); 
  FFV5_0(w[3], w[37], w[17], pars->GC_7, amp[106]); 
  FFV5_0(w[23], w[2], w[17], pars->GC_7, amp[107]); 
  FFV3_8_0(w[23], w[2], w[17], pars->GC_145, pars->GC_79, amp[108]); 
  FFV5_0(w[24], w[2], w[17], pars->GC_7, amp[109]); 
  VVV1_0(w[10], w[39], w[5], pars->GC_6, amp[110]); 
  VVV1_0(w[10], w[40], w[4], pars->GC_6, amp[111]); 
  FFV5_0(w[3], w[41], w[10], pars->GC_7, amp[112]); 
  FFV5_0(w[42], w[2], w[10], pars->GC_7, amp[113]); 
  FFVV1_2_0(w[44], w[43], w[4], w[5], pars->GC_146, pars->GC_85, amp[114]); 
  FFV5_0(w[44], w[45], w[5], pars->GC_7, amp[115]); 
  FFV3_8_0(w[44], w[45], w[5], pars->GC_145, pars->GC_79, amp[116]); 
  FFV5_0(w[44], w[46], w[5], pars->GC_7, amp[117]); 
  FFV5_0(w[44], w[47], w[4], pars->GC_7, amp[118]); 
  FFV3_8_0(w[44], w[47], w[4], pars->GC_145, pars->GC_79, amp[119]); 
  FFV5_0(w[44], w[48], w[4], pars->GC_7, amp[120]); 
  FFV5_0(w[49], w[45], w[5], pars->GC_7, amp[121]); 
  FFV5_0(w[49], w[47], w[4], pars->GC_7, amp[122]); 
  FFV5_0(w[44], w[51], w[5], pars->GC_7, amp[123]); 
  FFV5_0(w[44], w[52], w[4], pars->GC_7, amp[124]); 
  FFS2_0(w[44], w[43], w[38], pars->GC_95, amp[125]); 
  FFV5_0(w[44], w[43], w[18], pars->GC_7, amp[126]); 
  FFV5_0(w[44], w[43], w[17], pars->GC_7, amp[127]); 
  FFV3_8_0(w[44], w[43], w[17], pars->GC_145, pars->GC_79, amp[128]); 
  FFV5_0(w[49], w[43], w[17], pars->GC_7, amp[129]); 
  FFV5_0(w[44], w[50], w[17], pars->GC_7, amp[130]); 
  FFS2_0(w[3], w[47], w[53], pars->GC_95, amp[131]); 
  VVV1_0(w[54], w[5], w[55], pars->GC_6, amp[132]); 
  FFV5_0(w[3], w[47], w[54], pars->GC_7, amp[133]); 
  FFVV1_2_0(w[3], w[43], w[56], w[5], pars->GC_146, pars->GC_85, amp[134]); 
  VVS2_0(w[56], w[5], w[57], pars->GC_78, amp[135]); 
  VVV2_0(w[56], w[5], w[55], pars->GC_28, amp[136]); 
  VVV1_0(w[56], w[5], w[55], pars->GC_6, amp[137]); 
  VVV1_0(w[56], w[5], w[58], pars->GC_6, amp[138]); 
  FFV5_0(w[3], w[47], w[56], pars->GC_7, amp[139]); 
  FFV3_8_0(w[3], w[47], w[56], pars->GC_145, pars->GC_79, amp[140]); 
  FFV5_0(w[3], w[48], w[56], pars->GC_7, amp[141]); 
  VVV1_0(w[56], w[5], w[59], pars->GC_6, amp[142]); 
  FFV5_0(w[3], w[52], w[56], pars->GC_7, amp[143]); 
  FFS2_0(w[27], w[43], w[53], pars->GC_95, amp[144]); 
  FFV5_0(w[27], w[43], w[54], pars->GC_7, amp[145]); 
  FFV5_0(w[27], w[43], w[56], pars->GC_7, amp[146]); 
  FFV3_8_0(w[27], w[43], w[56], pars->GC_145, pars->GC_79, amp[147]); 
  FFV5_0(w[28], w[43], w[56], pars->GC_7, amp[148]); 
  FFV5_0(w[27], w[50], w[56], pars->GC_7, amp[149]); 
  FFS2_0(w[3], w[45], w[60], pars->GC_95, amp[150]); 
  VVV1_0(w[61], w[4], w[55], pars->GC_6, amp[151]); 
  FFV5_0(w[3], w[45], w[61], pars->GC_7, amp[152]); 
  FFVV1_2_0(w[3], w[43], w[62], w[4], pars->GC_146, pars->GC_85, amp[153]); 
  VVS2_0(w[62], w[4], w[57], pars->GC_78, amp[154]); 
  VVV2_0(w[62], w[4], w[55], pars->GC_28, amp[155]); 
  VVV1_0(w[62], w[4], w[55], pars->GC_6, amp[156]); 
  VVV1_0(w[62], w[4], w[58], pars->GC_6, amp[157]); 
  FFV5_0(w[3], w[45], w[62], pars->GC_7, amp[158]); 
  FFV3_8_0(w[3], w[45], w[62], pars->GC_145, pars->GC_79, amp[159]); 
  FFV5_0(w[3], w[46], w[62], pars->GC_7, amp[160]); 
  VVV1_0(w[62], w[4], w[59], pars->GC_6, amp[161]); 
  FFV5_0(w[3], w[51], w[62], pars->GC_7, amp[162]); 
  FFS2_0(w[32], w[43], w[60], pars->GC_95, amp[163]); 
  FFV5_0(w[32], w[43], w[61], pars->GC_7, amp[164]); 
  FFV5_0(w[32], w[43], w[62], pars->GC_7, amp[165]); 
  FFV3_8_0(w[32], w[43], w[62], pars->GC_145, pars->GC_79, amp[166]); 
  FFV5_0(w[33], w[43], w[62], pars->GC_7, amp[167]); 
  FFV5_0(w[32], w[50], w[62], pars->GC_7, amp[168]); 
  FFVV1_2_0(w[32], w[43], w[1], w[5], pars->GC_146, pars->GC_85, amp[169]); 
  FFV5_0(w[32], w[63], w[5], pars->GC_7, amp[170]); 
  FFV3_8_0(w[32], w[63], w[5], pars->GC_145, pars->GC_79, amp[171]); 
  FFV5_0(w[32], w[64], w[5], pars->GC_7, amp[172]); 
  FFV5_0(w[32], w[47], w[1], pars->GC_7, amp[173]); 
  FFV3_8_0(w[32], w[47], w[1], pars->GC_145, pars->GC_79, amp[174]); 
  FFV5_0(w[32], w[48], w[1], pars->GC_7, amp[175]); 
  FFV5_0(w[33], w[63], w[5], pars->GC_7, amp[176]); 
  FFV5_0(w[33], w[47], w[1], pars->GC_7, amp[177]); 
  FFV5_0(w[32], w[65], w[5], pars->GC_7, amp[178]); 
  FFV5_0(w[32], w[52], w[1], pars->GC_7, amp[179]); 
  FFVV1_2_0(w[27], w[43], w[1], w[4], pars->GC_146, pars->GC_85, amp[180]); 
  FFV5_0(w[27], w[63], w[4], pars->GC_7, amp[181]); 
  FFV3_8_0(w[27], w[63], w[4], pars->GC_145, pars->GC_79, amp[182]); 
  FFV5_0(w[27], w[64], w[4], pars->GC_7, amp[183]); 
  FFV5_0(w[27], w[45], w[1], pars->GC_7, amp[184]); 
  FFV3_8_0(w[27], w[45], w[1], pars->GC_145, pars->GC_79, amp[185]); 
  FFV5_0(w[27], w[46], w[1], pars->GC_7, amp[186]); 
  FFV5_0(w[28], w[63], w[4], pars->GC_7, amp[187]); 
  FFV5_0(w[28], w[45], w[1], pars->GC_7, amp[188]); 
  FFV5_0(w[27], w[65], w[4], pars->GC_7, amp[189]); 
  FFV5_0(w[27], w[51], w[1], pars->GC_7, amp[190]); 
  FFS2_0(w[3], w[63], w[38], pars->GC_95, amp[191]); 
  FFV5_0(w[3], w[63], w[18], pars->GC_7, amp[192]); 
  VVV1_0(w[1], w[18], w[55], pars->GC_6, amp[193]); 
  FFVV1_2_0(w[3], w[43], w[1], w[17], pars->GC_146, pars->GC_85, amp[194]); 
  FFV5_0(w[3], w[63], w[17], pars->GC_7, amp[195]); 
  FFV3_8_0(w[3], w[63], w[17], pars->GC_145, pars->GC_79, amp[196]); 
  FFV5_0(w[3], w[64], w[17], pars->GC_7, amp[197]); 
  VVS2_0(w[1], w[17], w[57], pars->GC_78, amp[198]); 
  VVV2_0(w[1], w[17], w[55], pars->GC_28, amp[199]); 
  VVV1_0(w[1], w[17], w[55], pars->GC_6, amp[200]); 
  VVV1_0(w[1], w[17], w[58], pars->GC_6, amp[201]); 
  FFV5_0(w[3], w[65], w[17], pars->GC_7, amp[202]); 
  VVV1_0(w[1], w[17], w[59], pars->GC_6, amp[203]); 
  FFV5_0(w[66], w[43], w[5], pars->GC_7, amp[204]); 
  FFV5_0(w[67], w[43], w[4], pars->GC_7, amp[205]); 
  FFV5_0(w[3], w[43], w[68], pars->GC_7, amp[206]); 
  FFV5_0(w[3], w[43], w[69], pars->GC_7, amp[207]); 
  FFV5_0(w[3], w[43], w[70], pars->GC_7, amp[208]); 
  FFV5_0(w[3], w[43], w[71], pars->GC_7, amp[209]); 
  FFV5_0(w[3], w[43], w[72], pars->GC_7, amp[210]); 
  FFV5_0(w[3], w[43], w[73], pars->GC_7, amp[211]); 
  FFV3_8_0(w[3], w[43], w[71], pars->GC_145, pars->GC_79, amp[212]); 
  FFV3_8_0(w[3], w[43], w[72], pars->GC_145, pars->GC_79, amp[213]); 
  FFV3_8_0(w[3], w[43], w[73], pars->GC_145, pars->GC_79, amp[214]); 
  FFS2_0(w[3], w[43], w[74], pars->GC_95, amp[215]); 
  FFV5_0(w[3], w[50], w[71], pars->GC_7, amp[216]); 
  FFV5_0(w[3], w[50], w[72], pars->GC_7, amp[217]); 
  FFV5_0(w[3], w[50], w[73], pars->GC_7, amp[218]); 
  FFV5_0(w[42], w[43], w[1], pars->GC_7, amp[219]); 
  FFVV1_2_0(w[75], w[76], w[4], w[5], pars->GC_146, pars->GC_85, amp[220]); 
  FFV5_0(w[77], w[76], w[5], pars->GC_7, amp[221]); 
  FFV3_8_0(w[77], w[76], w[5], pars->GC_145, pars->GC_79, amp[222]); 
  FFV5_0(w[78], w[76], w[5], pars->GC_7, amp[223]); 
  FFV5_0(w[79], w[76], w[4], pars->GC_7, amp[224]); 
  FFV3_8_0(w[79], w[76], w[4], pars->GC_145, pars->GC_79, amp[225]); 
  FFV5_0(w[80], w[76], w[4], pars->GC_7, amp[226]); 
  FFV5_0(w[77], w[81], w[5], pars->GC_7, amp[227]); 
  FFV5_0(w[79], w[81], w[4], pars->GC_7, amp[228]); 
  FFV5_0(w[83], w[76], w[5], pars->GC_7, amp[229]); 
  FFV5_0(w[84], w[76], w[4], pars->GC_7, amp[230]); 
  FFS2_0(w[75], w[76], w[38], pars->GC_95, amp[231]); 
  FFV5_0(w[75], w[76], w[18], pars->GC_7, amp[232]); 
  FFV5_0(w[75], w[76], w[17], pars->GC_7, amp[233]); 
  FFV3_8_0(w[75], w[76], w[17], pars->GC_145, pars->GC_79, amp[234]); 
  FFV5_0(w[75], w[81], w[17], pars->GC_7, amp[235]); 
  FFV5_0(w[82], w[76], w[17], pars->GC_7, amp[236]); 
  FFS2_0(w[79], w[2], w[53], pars->GC_95, amp[237]); 
  VVV1_0(w[54], w[5], w[85], pars->GC_6, amp[238]); 
  FFV5_0(w[79], w[2], w[54], pars->GC_7, amp[239]); 
  FFVV1_2_0(w[75], w[2], w[56], w[5], pars->GC_146, pars->GC_85, amp[240]); 
  VVS2_0(w[56], w[5], w[86], pars->GC_78, amp[241]); 
  VVV2_0(w[56], w[5], w[85], pars->GC_28, amp[242]); 
  VVV1_0(w[56], w[5], w[85], pars->GC_6, amp[243]); 
  VVV1_0(w[56], w[5], w[87], pars->GC_6, amp[244]); 
  FFV5_0(w[79], w[2], w[56], pars->GC_7, amp[245]); 
  FFV3_8_0(w[79], w[2], w[56], pars->GC_145, pars->GC_79, amp[246]); 
  FFV5_0(w[80], w[2], w[56], pars->GC_7, amp[247]); 
  VVV1_0(w[56], w[5], w[88], pars->GC_6, amp[248]); 
  FFV5_0(w[84], w[2], w[56], pars->GC_7, amp[249]); 
  FFS2_0(w[75], w[29], w[53], pars->GC_95, amp[250]); 
  FFV5_0(w[75], w[29], w[54], pars->GC_7, amp[251]); 
  FFV5_0(w[75], w[29], w[56], pars->GC_7, amp[252]); 
  FFV3_8_0(w[75], w[29], w[56], pars->GC_145, pars->GC_79, amp[253]); 
  FFV5_0(w[75], w[31], w[56], pars->GC_7, amp[254]); 
  FFV5_0(w[82], w[29], w[56], pars->GC_7, amp[255]); 
  FFS2_0(w[77], w[2], w[60], pars->GC_95, amp[256]); 
  VVV1_0(w[61], w[4], w[85], pars->GC_6, amp[257]); 
  FFV5_0(w[77], w[2], w[61], pars->GC_7, amp[258]); 
  FFVV1_2_0(w[75], w[2], w[62], w[4], pars->GC_146, pars->GC_85, amp[259]); 
  VVS2_0(w[62], w[4], w[86], pars->GC_78, amp[260]); 
  VVV2_0(w[62], w[4], w[85], pars->GC_28, amp[261]); 
  VVV1_0(w[62], w[4], w[85], pars->GC_6, amp[262]); 
  VVV1_0(w[62], w[4], w[87], pars->GC_6, amp[263]); 
  FFV5_0(w[77], w[2], w[62], pars->GC_7, amp[264]); 
  FFV3_8_0(w[77], w[2], w[62], pars->GC_145, pars->GC_79, amp[265]); 
  FFV5_0(w[78], w[2], w[62], pars->GC_7, amp[266]); 
  VVV1_0(w[62], w[4], w[88], pars->GC_6, amp[267]); 
  FFV5_0(w[83], w[2], w[62], pars->GC_7, amp[268]); 
  FFS2_0(w[75], w[20], w[60], pars->GC_95, amp[269]); 
  FFV5_0(w[75], w[20], w[61], pars->GC_7, amp[270]); 
  FFV5_0(w[75], w[20], w[62], pars->GC_7, amp[271]); 
  FFV3_8_0(w[75], w[20], w[62], pars->GC_145, pars->GC_79, amp[272]); 
  FFV5_0(w[75], w[26], w[62], pars->GC_7, amp[273]); 
  FFV5_0(w[82], w[20], w[62], pars->GC_7, amp[274]); 
  FFVV1_2_0(w[75], w[20], w[1], w[5], pars->GC_146, pars->GC_85, amp[275]); 
  FFV5_0(w[89], w[20], w[5], pars->GC_7, amp[276]); 
  FFV3_8_0(w[89], w[20], w[5], pars->GC_145, pars->GC_79, amp[277]); 
  FFV5_0(w[90], w[20], w[5], pars->GC_7, amp[278]); 
  FFV5_0(w[79], w[20], w[1], pars->GC_7, amp[279]); 
  FFV3_8_0(w[79], w[20], w[1], pars->GC_145, pars->GC_79, amp[280]); 
  FFV5_0(w[80], w[20], w[1], pars->GC_7, amp[281]); 
  FFV5_0(w[89], w[26], w[5], pars->GC_7, amp[282]); 
  FFV5_0(w[79], w[26], w[1], pars->GC_7, amp[283]); 
  FFV5_0(w[91], w[20], w[5], pars->GC_7, amp[284]); 
  FFV5_0(w[84], w[20], w[1], pars->GC_7, amp[285]); 
  FFVV1_2_0(w[75], w[29], w[1], w[4], pars->GC_146, pars->GC_85, amp[286]); 
  FFV5_0(w[89], w[29], w[4], pars->GC_7, amp[287]); 
  FFV3_8_0(w[89], w[29], w[4], pars->GC_145, pars->GC_79, amp[288]); 
  FFV5_0(w[90], w[29], w[4], pars->GC_7, amp[289]); 
  FFV5_0(w[77], w[29], w[1], pars->GC_7, amp[290]); 
  FFV3_8_0(w[77], w[29], w[1], pars->GC_145, pars->GC_79, amp[291]); 
  FFV5_0(w[78], w[29], w[1], pars->GC_7, amp[292]); 
  FFV5_0(w[89], w[31], w[4], pars->GC_7, amp[293]); 
  FFV5_0(w[77], w[31], w[1], pars->GC_7, amp[294]); 
  FFV5_0(w[91], w[29], w[4], pars->GC_7, amp[295]); 
  FFV5_0(w[83], w[29], w[1], pars->GC_7, amp[296]); 
  FFS2_0(w[89], w[2], w[38], pars->GC_95, amp[297]); 
  FFV5_0(w[89], w[2], w[18], pars->GC_7, amp[298]); 
  VVV1_0(w[1], w[18], w[85], pars->GC_6, amp[299]); 
  FFVV1_2_0(w[75], w[2], w[1], w[17], pars->GC_146, pars->GC_85, amp[300]); 
  FFV5_0(w[89], w[2], w[17], pars->GC_7, amp[301]); 
  FFV3_8_0(w[89], w[2], w[17], pars->GC_145, pars->GC_79, amp[302]); 
  FFV5_0(w[90], w[2], w[17], pars->GC_7, amp[303]); 
  VVS2_0(w[1], w[17], w[86], pars->GC_78, amp[304]); 
  VVV2_0(w[1], w[17], w[85], pars->GC_28, amp[305]); 
  VVV1_0(w[1], w[17], w[85], pars->GC_6, amp[306]); 
  VVV1_0(w[1], w[17], w[87], pars->GC_6, amp[307]); 
  FFV5_0(w[91], w[2], w[17], pars->GC_7, amp[308]); 
  VVV1_0(w[1], w[17], w[88], pars->GC_6, amp[309]); 
  FFV5_0(w[75], w[92], w[5], pars->GC_7, amp[310]); 
  FFV5_0(w[75], w[93], w[4], pars->GC_7, amp[311]); 
  FFV5_0(w[75], w[2], w[68], pars->GC_7, amp[312]); 
  FFV5_0(w[75], w[2], w[69], pars->GC_7, amp[313]); 
  FFV5_0(w[75], w[2], w[70], pars->GC_7, amp[314]); 
  FFV5_0(w[75], w[2], w[71], pars->GC_7, amp[315]); 
  FFV5_0(w[75], w[2], w[72], pars->GC_7, amp[316]); 
  FFV5_0(w[75], w[2], w[73], pars->GC_7, amp[317]); 
  FFV3_8_0(w[75], w[2], w[71], pars->GC_145, pars->GC_79, amp[318]); 
  FFV3_8_0(w[75], w[2], w[72], pars->GC_145, pars->GC_79, amp[319]); 
  FFV3_8_0(w[75], w[2], w[73], pars->GC_145, pars->GC_79, amp[320]); 
  FFS2_0(w[75], w[2], w[74], pars->GC_95, amp[321]); 
  FFV5_0(w[82], w[2], w[71], pars->GC_7, amp[322]); 
  FFV5_0(w[82], w[2], w[72], pars->GC_7, amp[323]); 
  FFV5_0(w[82], w[2], w[73], pars->GC_7, amp[324]); 
  FFV5_0(w[75], w[41], w[1], pars->GC_7, amp[325]); 
  FFV5_0(w[95], w[76], w[5], pars->GC_7, amp[326]); 
  FFV5_0(w[97], w[76], w[5], pars->GC_7, amp[327]); 
  FFV5_0(w[3], w[76], w[98], pars->GC_7, amp[328]); 
  FFVV1_2_0(w[3], w[76], w[99], w[5], pars->GC_146, pars->GC_85, amp[329]); 
  FFV5_0(w[100], w[76], w[5], pars->GC_7, amp[330]); 
  FFV3_8_0(w[100], w[76], w[5], pars->GC_145, pars->GC_79, amp[331]); 
  FFV5_0(w[101], w[76], w[5], pars->GC_7, amp[332]); 
  FFS2_0(w[3], w[76], w[102], pars->GC_95, amp[333]); 
  FFV5_0(w[3], w[76], w[103], pars->GC_7, amp[334]); 
  FFV5_0(w[3], w[76], w[104], pars->GC_7, amp[335]); 
  FFV3_8_0(w[3], w[76], w[104], pars->GC_145, pars->GC_79, amp[336]); 
  FFV5_0(w[100], w[81], w[5], pars->GC_7, amp[337]); 
  FFV5_0(w[3], w[81], w[104], pars->GC_7, amp[338]); 
  FFS2_0(w[27], w[76], w[94], pars->GC_95, amp[339]); 
  FFV5_0(w[27], w[76], w[96], pars->GC_7, amp[340]); 
  FFV5_0(w[27], w[76], w[99], pars->GC_7, amp[341]); 
  FFV3_8_0(w[27], w[76], w[99], pars->GC_145, pars->GC_79, amp[342]); 
  FFV5_0(w[28], w[76], w[99], pars->GC_7, amp[343]); 
  FFV5_0(w[27], w[81], w[99], pars->GC_7, amp[344]); 
  FFV5_0(w[44], w[105], w[5], pars->GC_7, amp[345]); 
  FFV5_0(w[44], w[106], w[5], pars->GC_7, amp[346]); 
  FFV5_0(w[44], w[2], w[98], pars->GC_7, amp[347]); 
  FFVV1_2_0(w[44], w[2], w[99], w[5], pars->GC_146, pars->GC_85, amp[348]); 
  FFV5_0(w[44], w[107], w[5], pars->GC_7, amp[349]); 
  FFV3_8_0(w[44], w[107], w[5], pars->GC_145, pars->GC_79, amp[350]); 
  FFV5_0(w[44], w[108], w[5], pars->GC_7, amp[351]); 
  FFS2_0(w[44], w[2], w[102], pars->GC_95, amp[352]); 
  FFV5_0(w[44], w[2], w[103], pars->GC_7, amp[353]); 
  FFV5_0(w[44], w[2], w[104], pars->GC_7, amp[354]); 
  FFV3_8_0(w[44], w[2], w[104], pars->GC_145, pars->GC_79, amp[355]); 
  FFV5_0(w[49], w[107], w[5], pars->GC_7, amp[356]); 
  FFV5_0(w[49], w[2], w[104], pars->GC_7, amp[357]); 
  FFS2_0(w[44], w[29], w[94], pars->GC_95, amp[358]); 
  FFV5_0(w[44], w[29], w[96], pars->GC_7, amp[359]); 
  FFV5_0(w[44], w[29], w[99], pars->GC_7, amp[360]); 
  FFV3_8_0(w[44], w[29], w[99], pars->GC_145, pars->GC_79, amp[361]); 
  FFV5_0(w[44], w[31], w[99], pars->GC_7, amp[362]); 
  FFV5_0(w[49], w[29], w[99], pars->GC_7, amp[363]); 
  FFV5_0(w[3], w[105], w[62], pars->GC_7, amp[364]); 
  FFV5_0(w[95], w[2], w[62], pars->GC_7, amp[365]); 
  FFV5_0(w[3], w[106], w[62], pars->GC_7, amp[366]); 
  FFV5_0(w[97], w[2], w[62], pars->GC_7, amp[367]); 
  FFS2_0(w[3], w[107], w[60], pars->GC_95, amp[368]); 
  FFS2_0(w[100], w[2], w[60], pars->GC_95, amp[369]); 
  FFV5_0(w[3], w[107], w[61], pars->GC_7, amp[370]); 
  FFV5_0(w[100], w[2], w[61], pars->GC_7, amp[371]); 
  FFVV1_2_0(w[3], w[2], w[99], w[62], pars->GC_146, pars->GC_85, amp[372]); 
  FFV5_0(w[3], w[107], w[62], pars->GC_7, amp[373]); 
  FFV3_8_0(w[3], w[107], w[62], pars->GC_145, pars->GC_79, amp[374]); 
  FFV5_0(w[3], w[108], w[62], pars->GC_7, amp[375]); 
  FFV5_0(w[100], w[2], w[62], pars->GC_7, amp[376]); 
  FFV3_8_0(w[100], w[2], w[62], pars->GC_145, pars->GC_79, amp[377]); 
  FFV5_0(w[101], w[2], w[62], pars->GC_7, amp[378]); 
  VVV1_0(w[96], w[62], w[7], pars->GC_6, amp[379]); 
  VVV1_0(w[99], w[61], w[7], pars->GC_6, amp[380]); 
  VVS2_0(w[99], w[62], w[11], pars->GC_78, amp[381]); 
  VVV2_0(w[99], w[62], w[7], pars->GC_28, amp[382]); 
  VVV1_0(w[99], w[62], w[7], pars->GC_6, amp[383]); 
  VVV1_0(w[99], w[62], w[16], pars->GC_6, amp[384]); 
  VVVV1_0(w[96], w[1], w[7], w[5], pars->GC_8, amp[385]); 
  VVVV4_0(w[96], w[1], w[7], w[5], pars->GC_8, amp[386]); 
  VVVV5_0(w[96], w[1], w[7], w[5], pars->GC_8, amp[387]); 
  VVV1_0(w[7], w[5], w[109], pars->GC_6, amp[388]); 
  VVV1_0(w[1], w[7], w[98], pars->GC_6, amp[389]); 
  VVVS1_0(w[99], w[1], w[5], w[11], pars->GC_84, amp[390]); 
  VVS2_0(w[5], w[110], w[11], pars->GC_78, amp[391]); 
  VVS2_0(w[1], w[104], w[11], pars->GC_78, amp[392]); 
  VVVV2_0(w[99], w[1], w[7], w[5], pars->GC_44, amp[393]); 
  VVVV8_0(w[99], w[1], w[7], w[5], pars->GC_44, amp[394]); 
  VVVV7_0(w[99], w[1], w[7], w[5], pars->GC_44, amp[395]); 
  VVVV1_0(w[99], w[1], w[7], w[5], pars->GC_8, amp[396]); 
  VVVV4_0(w[99], w[1], w[7], w[5], pars->GC_8, amp[397]); 
  VVVV5_0(w[99], w[1], w[7], w[5], pars->GC_8, amp[398]); 
  VVV1_0(w[7], w[5], w[111], pars->GC_6, amp[399]); 
  VVV2_0(w[7], w[5], w[110], pars->GC_28, amp[400]); 
  VVV1_0(w[7], w[5], w[110], pars->GC_6, amp[401]); 
  VVV1_0(w[1], w[7], w[103], pars->GC_6, amp[402]); 
  VVV2_0(w[1], w[7], w[104], pars->GC_28, amp[403]); 
  VVV1_0(w[1], w[7], w[104], pars->GC_6, amp[404]); 
  VVVV1_0(w[99], w[1], w[16], w[5], pars->GC_8, amp[405]); 
  VVVV4_0(w[99], w[1], w[16], w[5], pars->GC_8, amp[406]); 
  VVVV5_0(w[99], w[1], w[16], w[5], pars->GC_8, amp[407]); 
  VVV1_0(w[16], w[5], w[110], pars->GC_6, amp[408]); 
  VVV1_0(w[1], w[16], w[104], pars->GC_6, amp[409]); 
  FFV5_0(w[95], w[29], w[1], pars->GC_7, amp[410]); 
  FFV5_0(w[3], w[29], w[109], pars->GC_7, amp[411]); 
  FFV5_0(w[97], w[29], w[1], pars->GC_7, amp[412]); 
  FFVV1_2_0(w[3], w[29], w[99], w[1], pars->GC_146, pars->GC_85, amp[413]); 
  FFS2_0(w[3], w[29], w[112], pars->GC_95, amp[414]); 
  FFV5_0(w[3], w[29], w[111], pars->GC_7, amp[415]); 
  FFV5_0(w[3], w[29], w[110], pars->GC_7, amp[416]); 
  FFV3_8_0(w[3], w[29], w[110], pars->GC_145, pars->GC_79, amp[417]); 
  FFV5_0(w[100], w[29], w[1], pars->GC_7, amp[418]); 
  FFV3_8_0(w[100], w[29], w[1], pars->GC_145, pars->GC_79, amp[419]); 
  FFV5_0(w[101], w[29], w[1], pars->GC_7, amp[420]); 
  FFV5_0(w[3], w[31], w[110], pars->GC_7, amp[421]); 
  FFV5_0(w[100], w[31], w[1], pars->GC_7, amp[422]); 
  FFV5_0(w[27], w[105], w[1], pars->GC_7, amp[423]); 
  FFV5_0(w[27], w[2], w[109], pars->GC_7, amp[424]); 
  FFV5_0(w[27], w[106], w[1], pars->GC_7, amp[425]); 
  FFVV1_2_0(w[27], w[2], w[99], w[1], pars->GC_146, pars->GC_85, amp[426]); 
  FFS2_0(w[27], w[2], w[112], pars->GC_95, amp[427]); 
  FFV5_0(w[27], w[2], w[111], pars->GC_7, amp[428]); 
  FFV5_0(w[27], w[2], w[110], pars->GC_7, amp[429]); 
  FFV3_8_0(w[27], w[2], w[110], pars->GC_145, pars->GC_79, amp[430]); 
  FFV5_0(w[27], w[107], w[1], pars->GC_7, amp[431]); 
  FFV3_8_0(w[27], w[107], w[1], pars->GC_145, pars->GC_79, amp[432]); 
  FFV5_0(w[27], w[108], w[1], pars->GC_7, amp[433]); 
  FFV5_0(w[28], w[2], w[110], pars->GC_7, amp[434]); 
  FFV5_0(w[28], w[107], w[1], pars->GC_7, amp[435]); 
  VVV1_0(w[99], w[113], w[5], pars->GC_6, amp[436]); 
  FFV5_0(w[3], w[93], w[99], pars->GC_7, amp[437]); 
  FFV5_0(w[67], w[2], w[99], pars->GC_7, amp[438]); 
  VVV1_0(w[99], w[1], w[40], pars->GC_6, amp[439]); 
  FFV5_0(w[115], w[76], w[4], pars->GC_7, amp[440]); 
  FFV5_0(w[117], w[76], w[4], pars->GC_7, amp[441]); 
  FFV5_0(w[3], w[76], w[118], pars->GC_7, amp[442]); 
  FFVV1_2_0(w[3], w[76], w[119], w[4], pars->GC_146, pars->GC_85, amp[443]); 
  FFV5_0(w[120], w[76], w[4], pars->GC_7, amp[444]); 
  FFV3_8_0(w[120], w[76], w[4], pars->GC_145, pars->GC_79, amp[445]); 
  FFV5_0(w[121], w[76], w[4], pars->GC_7, amp[446]); 
  FFS2_0(w[3], w[76], w[122], pars->GC_95, amp[447]); 
  FFV5_0(w[3], w[76], w[123], pars->GC_7, amp[448]); 
  FFV5_0(w[3], w[76], w[124], pars->GC_7, amp[449]); 
  FFV3_8_0(w[3], w[76], w[124], pars->GC_145, pars->GC_79, amp[450]); 
  FFV5_0(w[120], w[81], w[4], pars->GC_7, amp[451]); 
  FFV5_0(w[3], w[81], w[124], pars->GC_7, amp[452]); 
  FFS2_0(w[32], w[76], w[114], pars->GC_95, amp[453]); 
  FFV5_0(w[32], w[76], w[116], pars->GC_7, amp[454]); 
  FFV5_0(w[32], w[76], w[119], pars->GC_7, amp[455]); 
  FFV3_8_0(w[32], w[76], w[119], pars->GC_145, pars->GC_79, amp[456]); 
  FFV5_0(w[33], w[76], w[119], pars->GC_7, amp[457]); 
  FFV5_0(w[32], w[81], w[119], pars->GC_7, amp[458]); 
  FFV5_0(w[44], w[125], w[4], pars->GC_7, amp[459]); 
  FFV5_0(w[44], w[126], w[4], pars->GC_7, amp[460]); 
  FFV5_0(w[44], w[2], w[118], pars->GC_7, amp[461]); 
  FFVV1_2_0(w[44], w[2], w[119], w[4], pars->GC_146, pars->GC_85, amp[462]); 
  FFV5_0(w[44], w[127], w[4], pars->GC_7, amp[463]); 
  FFV3_8_0(w[44], w[127], w[4], pars->GC_145, pars->GC_79, amp[464]); 
  FFV5_0(w[44], w[128], w[4], pars->GC_7, amp[465]); 
  FFS2_0(w[44], w[2], w[122], pars->GC_95, amp[466]); 
  FFV5_0(w[44], w[2], w[123], pars->GC_7, amp[467]); 
  FFV5_0(w[44], w[2], w[124], pars->GC_7, amp[468]); 
  FFV3_8_0(w[44], w[2], w[124], pars->GC_145, pars->GC_79, amp[469]); 
  FFV5_0(w[49], w[127], w[4], pars->GC_7, amp[470]); 
  FFV5_0(w[49], w[2], w[124], pars->GC_7, amp[471]); 
  FFS2_0(w[44], w[20], w[114], pars->GC_95, amp[472]); 
  FFV5_0(w[44], w[20], w[116], pars->GC_7, amp[473]); 
  FFV5_0(w[44], w[20], w[119], pars->GC_7, amp[474]); 
  FFV3_8_0(w[44], w[20], w[119], pars->GC_145, pars->GC_79, amp[475]); 
  FFV5_0(w[44], w[26], w[119], pars->GC_7, amp[476]); 
  FFV5_0(w[49], w[20], w[119], pars->GC_7, amp[477]); 
  FFV5_0(w[3], w[125], w[56], pars->GC_7, amp[478]); 
  FFV5_0(w[115], w[2], w[56], pars->GC_7, amp[479]); 
  FFV5_0(w[3], w[126], w[56], pars->GC_7, amp[480]); 
  FFV5_0(w[117], w[2], w[56], pars->GC_7, amp[481]); 
  FFS2_0(w[3], w[127], w[53], pars->GC_95, amp[482]); 
  FFS2_0(w[120], w[2], w[53], pars->GC_95, amp[483]); 
  FFV5_0(w[3], w[127], w[54], pars->GC_7, amp[484]); 
  FFV5_0(w[120], w[2], w[54], pars->GC_7, amp[485]); 
  FFVV1_2_0(w[3], w[2], w[119], w[56], pars->GC_146, pars->GC_85, amp[486]); 
  FFV5_0(w[3], w[127], w[56], pars->GC_7, amp[487]); 
  FFV3_8_0(w[3], w[127], w[56], pars->GC_145, pars->GC_79, amp[488]); 
  FFV5_0(w[3], w[128], w[56], pars->GC_7, amp[489]); 
  FFV5_0(w[120], w[2], w[56], pars->GC_7, amp[490]); 
  FFV3_8_0(w[120], w[2], w[56], pars->GC_145, pars->GC_79, amp[491]); 
  FFV5_0(w[121], w[2], w[56], pars->GC_7, amp[492]); 
  VVV1_0(w[116], w[56], w[7], pars->GC_6, amp[493]); 
  VVV1_0(w[119], w[54], w[7], pars->GC_6, amp[494]); 
  VVS2_0(w[119], w[56], w[11], pars->GC_78, amp[495]); 
  VVV2_0(w[119], w[56], w[7], pars->GC_28, amp[496]); 
  VVV1_0(w[119], w[56], w[7], pars->GC_6, amp[497]); 
  VVV1_0(w[119], w[56], w[16], pars->GC_6, amp[498]); 
  VVVV1_0(w[116], w[1], w[7], w[4], pars->GC_8, amp[499]); 
  VVVV4_0(w[116], w[1], w[7], w[4], pars->GC_8, amp[500]); 
  VVVV5_0(w[116], w[1], w[7], w[4], pars->GC_8, amp[501]); 
  VVV1_0(w[7], w[4], w[129], pars->GC_6, amp[502]); 
  VVV1_0(w[1], w[7], w[118], pars->GC_6, amp[503]); 
  VVVS1_0(w[119], w[1], w[4], w[11], pars->GC_84, amp[504]); 
  VVS2_0(w[4], w[130], w[11], pars->GC_78, amp[505]); 
  VVS2_0(w[1], w[124], w[11], pars->GC_78, amp[506]); 
  VVVV2_0(w[119], w[1], w[7], w[4], pars->GC_44, amp[507]); 
  VVVV8_0(w[119], w[1], w[7], w[4], pars->GC_44, amp[508]); 
  VVVV7_0(w[119], w[1], w[7], w[4], pars->GC_44, amp[509]); 
  VVVV1_0(w[119], w[1], w[7], w[4], pars->GC_8, amp[510]); 
  VVVV4_0(w[119], w[1], w[7], w[4], pars->GC_8, amp[511]); 
  VVVV5_0(w[119], w[1], w[7], w[4], pars->GC_8, amp[512]); 
  VVV1_0(w[7], w[4], w[131], pars->GC_6, amp[513]); 
  VVV2_0(w[7], w[4], w[130], pars->GC_28, amp[514]); 
  VVV1_0(w[7], w[4], w[130], pars->GC_6, amp[515]); 
  VVV1_0(w[1], w[7], w[123], pars->GC_6, amp[516]); 
  VVV2_0(w[1], w[7], w[124], pars->GC_28, amp[517]); 
  VVV1_0(w[1], w[7], w[124], pars->GC_6, amp[518]); 
  VVVV1_0(w[119], w[1], w[16], w[4], pars->GC_8, amp[519]); 
  VVVV4_0(w[119], w[1], w[16], w[4], pars->GC_8, amp[520]); 
  VVVV5_0(w[119], w[1], w[16], w[4], pars->GC_8, amp[521]); 
  VVV1_0(w[16], w[4], w[130], pars->GC_6, amp[522]); 
  VVV1_0(w[1], w[16], w[124], pars->GC_6, amp[523]); 
  FFV5_0(w[115], w[20], w[1], pars->GC_7, amp[524]); 
  FFV5_0(w[3], w[20], w[129], pars->GC_7, amp[525]); 
  FFV5_0(w[117], w[20], w[1], pars->GC_7, amp[526]); 
  FFVV1_2_0(w[3], w[20], w[119], w[1], pars->GC_146, pars->GC_85, amp[527]); 
  FFS2_0(w[3], w[20], w[132], pars->GC_95, amp[528]); 
  FFV5_0(w[3], w[20], w[131], pars->GC_7, amp[529]); 
  FFV5_0(w[3], w[20], w[130], pars->GC_7, amp[530]); 
  FFV3_8_0(w[3], w[20], w[130], pars->GC_145, pars->GC_79, amp[531]); 
  FFV5_0(w[120], w[20], w[1], pars->GC_7, amp[532]); 
  FFV3_8_0(w[120], w[20], w[1], pars->GC_145, pars->GC_79, amp[533]); 
  FFV5_0(w[121], w[20], w[1], pars->GC_7, amp[534]); 
  FFV5_0(w[3], w[26], w[130], pars->GC_7, amp[535]); 
  FFV5_0(w[120], w[26], w[1], pars->GC_7, amp[536]); 
  FFV5_0(w[32], w[125], w[1], pars->GC_7, amp[537]); 
  FFV5_0(w[32], w[2], w[129], pars->GC_7, amp[538]); 
  FFV5_0(w[32], w[126], w[1], pars->GC_7, amp[539]); 
  FFVV1_2_0(w[32], w[2], w[119], w[1], pars->GC_146, pars->GC_85, amp[540]); 
  FFS2_0(w[32], w[2], w[132], pars->GC_95, amp[541]); 
  FFV5_0(w[32], w[2], w[131], pars->GC_7, amp[542]); 
  FFV5_0(w[32], w[2], w[130], pars->GC_7, amp[543]); 
  FFV3_8_0(w[32], w[2], w[130], pars->GC_145, pars->GC_79, amp[544]); 
  FFV5_0(w[32], w[127], w[1], pars->GC_7, amp[545]); 
  FFV3_8_0(w[32], w[127], w[1], pars->GC_145, pars->GC_79, amp[546]); 
  FFV5_0(w[32], w[128], w[1], pars->GC_7, amp[547]); 
  FFV5_0(w[33], w[2], w[130], pars->GC_7, amp[548]); 
  FFV5_0(w[33], w[127], w[1], pars->GC_7, amp[549]); 
  VVV1_0(w[119], w[113], w[4], pars->GC_6, amp[550]); 
  FFV5_0(w[3], w[92], w[119], pars->GC_7, amp[551]); 
  FFV5_0(w[66], w[2], w[119], pars->GC_7, amp[552]); 
  VVV1_0(w[119], w[1], w[39], pars->GC_6, amp[553]); 
  FFVV1_2_0(w[32], w[76], w[0], w[5], pars->GC_146, pars->GC_85, amp[554]); 
  FFV5_0(w[32], w[133], w[5], pars->GC_7, amp[555]); 
  FFV3_8_0(w[32], w[133], w[5], pars->GC_145, pars->GC_79, amp[556]); 
  FFV5_0(w[32], w[134], w[5], pars->GC_7, amp[557]); 
  FFV5_0(w[135], w[76], w[5], pars->GC_7, amp[558]); 
  FFV3_8_0(w[135], w[76], w[5], pars->GC_145, pars->GC_79, amp[559]); 
  FFV5_0(w[136], w[76], w[5], pars->GC_7, amp[560]); 
  FFV5_0(w[33], w[133], w[5], pars->GC_7, amp[561]); 
  FFV5_0(w[137], w[76], w[5], pars->GC_7, amp[562]); 
  FFV5_0(w[32], w[138], w[5], pars->GC_7, amp[563]); 
  FFV5_0(w[135], w[81], w[5], pars->GC_7, amp[564]); 
  FFVV1_2_0(w[27], w[76], w[0], w[4], pars->GC_146, pars->GC_85, amp[565]); 
  FFV5_0(w[27], w[133], w[4], pars->GC_7, amp[566]); 
  FFV3_8_0(w[27], w[133], w[4], pars->GC_145, pars->GC_79, amp[567]); 
  FFV5_0(w[27], w[134], w[4], pars->GC_7, amp[568]); 
  FFV5_0(w[139], w[76], w[4], pars->GC_7, amp[569]); 
  FFV3_8_0(w[139], w[76], w[4], pars->GC_145, pars->GC_79, amp[570]); 
  FFV5_0(w[140], w[76], w[4], pars->GC_7, amp[571]); 
  FFV5_0(w[28], w[133], w[4], pars->GC_7, amp[572]); 
  FFV5_0(w[141], w[76], w[4], pars->GC_7, amp[573]); 
  FFV5_0(w[27], w[138], w[4], pars->GC_7, amp[574]); 
  FFV5_0(w[139], w[81], w[4], pars->GC_7, amp[575]); 
  FFS2_0(w[3], w[133], w[38], pars->GC_95, amp[576]); 
  FFV5_0(w[3], w[133], w[18], pars->GC_7, amp[577]); 
  FFV5_0(w[3], w[76], w[142], pars->GC_7, amp[578]); 
  FFVV1_2_0(w[3], w[76], w[0], w[17], pars->GC_146, pars->GC_85, amp[579]); 
  FFV5_0(w[3], w[133], w[17], pars->GC_7, amp[580]); 
  FFV3_8_0(w[3], w[133], w[17], pars->GC_145, pars->GC_79, amp[581]); 
  FFV5_0(w[3], w[134], w[17], pars->GC_7, amp[582]); 
  FFS2_0(w[3], w[76], w[143], pars->GC_95, amp[583]); 
  FFV5_0(w[3], w[76], w[144], pars->GC_7, amp[584]); 
  FFV5_0(w[3], w[76], w[145], pars->GC_7, amp[585]); 
  FFV3_8_0(w[3], w[76], w[145], pars->GC_145, pars->GC_79, amp[586]); 
  FFV5_0(w[3], w[138], w[17], pars->GC_7, amp[587]); 
  FFV5_0(w[3], w[81], w[145], pars->GC_7, amp[588]); 
  FFV5_0(w[42], w[76], w[0], pars->GC_7, amp[589]); 
  FFVV1_2_0(w[44], w[20], w[0], w[5], pars->GC_146, pars->GC_85, amp[590]); 
  FFV5_0(w[146], w[20], w[5], pars->GC_7, amp[591]); 
  FFV3_8_0(w[146], w[20], w[5], pars->GC_145, pars->GC_79, amp[592]); 
  FFV5_0(w[147], w[20], w[5], pars->GC_7, amp[593]); 
  FFV5_0(w[44], w[148], w[5], pars->GC_7, amp[594]); 
  FFV3_8_0(w[44], w[148], w[5], pars->GC_145, pars->GC_79, amp[595]); 
  FFV5_0(w[44], w[149], w[5], pars->GC_7, amp[596]); 
  FFV5_0(w[146], w[26], w[5], pars->GC_7, amp[597]); 
  FFV5_0(w[44], w[150], w[5], pars->GC_7, amp[598]); 
  FFV5_0(w[151], w[20], w[5], pars->GC_7, amp[599]); 
  FFV5_0(w[49], w[148], w[5], pars->GC_7, amp[600]); 
  FFVV1_2_0(w[44], w[29], w[0], w[4], pars->GC_146, pars->GC_85, amp[601]); 
  FFV5_0(w[146], w[29], w[4], pars->GC_7, amp[602]); 
  FFV3_8_0(w[146], w[29], w[4], pars->GC_145, pars->GC_79, amp[603]); 
  FFV5_0(w[147], w[29], w[4], pars->GC_7, amp[604]); 
  FFV5_0(w[44], w[152], w[4], pars->GC_7, amp[605]); 
  FFV3_8_0(w[44], w[152], w[4], pars->GC_145, pars->GC_79, amp[606]); 
  FFV5_0(w[44], w[153], w[4], pars->GC_7, amp[607]); 
  FFV5_0(w[146], w[31], w[4], pars->GC_7, amp[608]); 
  FFV5_0(w[44], w[154], w[4], pars->GC_7, amp[609]); 
  FFV5_0(w[151], w[29], w[4], pars->GC_7, amp[610]); 
  FFV5_0(w[49], w[152], w[4], pars->GC_7, amp[611]); 
  FFS2_0(w[146], w[2], w[38], pars->GC_95, amp[612]); 
  FFV5_0(w[146], w[2], w[18], pars->GC_7, amp[613]); 
  FFV5_0(w[44], w[2], w[142], pars->GC_7, amp[614]); 
  FFVV1_2_0(w[44], w[2], w[0], w[17], pars->GC_146, pars->GC_85, amp[615]); 
  FFV5_0(w[146], w[2], w[17], pars->GC_7, amp[616]); 
  FFV3_8_0(w[146], w[2], w[17], pars->GC_145, pars->GC_79, amp[617]); 
  FFV5_0(w[147], w[2], w[17], pars->GC_7, amp[618]); 
  FFS2_0(w[44], w[2], w[143], pars->GC_95, amp[619]); 
  FFV5_0(w[44], w[2], w[144], pars->GC_7, amp[620]); 
  FFV5_0(w[44], w[2], w[145], pars->GC_7, amp[621]); 
  FFV3_8_0(w[44], w[2], w[145], pars->GC_145, pars->GC_79, amp[622]); 
  FFV5_0(w[151], w[2], w[17], pars->GC_7, amp[623]); 
  FFV5_0(w[49], w[2], w[145], pars->GC_7, amp[624]); 
  FFV5_0(w[44], w[41], w[0], pars->GC_7, amp[625]); 
  VVVV1_0(w[0], w[54], w[7], w[5], pars->GC_8, amp[626]); 
  VVVV4_0(w[0], w[54], w[7], w[5], pars->GC_8, amp[627]); 
  VVVV5_0(w[0], w[54], w[7], w[5], pars->GC_8, amp[628]); 
  VVV1_0(w[7], w[5], w[155], pars->GC_6, amp[629]); 
  VVV1_0(w[54], w[5], w[156], pars->GC_6, amp[630]); 
  VVVS1_0(w[0], w[56], w[5], w[11], pars->GC_84, amp[631]); 
  VVS2_0(w[5], w[157], w[11], pars->GC_78, amp[632]); 
  VVV1_0(w[56], w[5], w[158], pars->GC_6, amp[633]); 
  VVVV2_0(w[0], w[56], w[7], w[5], pars->GC_44, amp[634]); 
  VVVV8_0(w[0], w[56], w[7], w[5], pars->GC_44, amp[635]); 
  VVVV7_0(w[0], w[56], w[7], w[5], pars->GC_44, amp[636]); 
  VVVV1_0(w[0], w[56], w[7], w[5], pars->GC_8, amp[637]); 
  VVVV4_0(w[0], w[56], w[7], w[5], pars->GC_8, amp[638]); 
  VVVV5_0(w[0], w[56], w[7], w[5], pars->GC_8, amp[639]); 
  VVV1_0(w[7], w[5], w[159], pars->GC_6, amp[640]); 
  VVV2_0(w[7], w[5], w[157], pars->GC_28, amp[641]); 
  VVV1_0(w[7], w[5], w[157], pars->GC_6, amp[642]); 
  VVV1_0(w[56], w[5], w[160], pars->GC_6, amp[643]); 
  VVV2_0(w[56], w[5], w[156], pars->GC_28, amp[644]); 
  VVV1_0(w[56], w[5], w[156], pars->GC_6, amp[645]); 
  VVVV1_0(w[0], w[56], w[16], w[5], pars->GC_8, amp[646]); 
  VVVV4_0(w[0], w[56], w[16], w[5], pars->GC_8, amp[647]); 
  VVVV5_0(w[0], w[56], w[16], w[5], pars->GC_8, amp[648]); 
  VVV1_0(w[16], w[5], w[157], pars->GC_6, amp[649]); 
  VVV1_0(w[56], w[5], w[161], pars->GC_6, amp[650]); 
  FFS2_0(w[3], w[152], w[53], pars->GC_95, amp[651]); 
  FFV5_0(w[3], w[29], w[155], pars->GC_7, amp[652]); 
  FFV5_0(w[3], w[152], w[54], pars->GC_7, amp[653]); 
  FFVV1_2_0(w[3], w[29], w[0], w[56], pars->GC_146, pars->GC_85, amp[654]); 
  FFS2_0(w[3], w[29], w[162], pars->GC_95, amp[655]); 
  FFV5_0(w[3], w[29], w[159], pars->GC_7, amp[656]); 
  FFV5_0(w[3], w[29], w[157], pars->GC_7, amp[657]); 
  FFV3_8_0(w[3], w[29], w[157], pars->GC_145, pars->GC_79, amp[658]); 
  FFV5_0(w[3], w[152], w[56], pars->GC_7, amp[659]); 
  FFV3_8_0(w[3], w[152], w[56], pars->GC_145, pars->GC_79, amp[660]); 
  FFV5_0(w[3], w[153], w[56], pars->GC_7, amp[661]); 
  FFV5_0(w[3], w[31], w[157], pars->GC_7, amp[662]); 
  FFV5_0(w[3], w[154], w[56], pars->GC_7, amp[663]); 
  FFS2_0(w[139], w[2], w[53], pars->GC_95, amp[664]); 
  FFV5_0(w[27], w[2], w[155], pars->GC_7, amp[665]); 
  FFV5_0(w[139], w[2], w[54], pars->GC_7, amp[666]); 
  FFVV1_2_0(w[27], w[2], w[0], w[56], pars->GC_146, pars->GC_85, amp[667]); 
  FFS2_0(w[27], w[2], w[162], pars->GC_95, amp[668]); 
  FFV5_0(w[27], w[2], w[159], pars->GC_7, amp[669]); 
  FFV5_0(w[27], w[2], w[157], pars->GC_7, amp[670]); 
  FFV3_8_0(w[27], w[2], w[157], pars->GC_145, pars->GC_79, amp[671]); 
  FFV5_0(w[139], w[2], w[56], pars->GC_7, amp[672]); 
  FFV3_8_0(w[139], w[2], w[56], pars->GC_145, pars->GC_79, amp[673]); 
  FFV5_0(w[140], w[2], w[56], pars->GC_7, amp[674]); 
  FFV5_0(w[28], w[2], w[157], pars->GC_7, amp[675]); 
  FFV5_0(w[141], w[2], w[56], pars->GC_7, amp[676]); 
  VVV1_0(w[0], w[56], w[40], pars->GC_6, amp[677]); 
  VVVV1_0(w[0], w[61], w[7], w[4], pars->GC_8, amp[678]); 
  VVVV4_0(w[0], w[61], w[7], w[4], pars->GC_8, amp[679]); 
  VVVV5_0(w[0], w[61], w[7], w[4], pars->GC_8, amp[680]); 
  VVV1_0(w[7], w[4], w[163], pars->GC_6, amp[681]); 
  VVV1_0(w[61], w[4], w[156], pars->GC_6, amp[682]); 
  VVVS1_0(w[0], w[62], w[4], w[11], pars->GC_84, amp[683]); 
  VVS2_0(w[4], w[164], w[11], pars->GC_78, amp[684]); 
  VVV1_0(w[62], w[4], w[158], pars->GC_6, amp[685]); 
  VVVV2_0(w[0], w[62], w[7], w[4], pars->GC_44, amp[686]); 
  VVVV8_0(w[0], w[62], w[7], w[4], pars->GC_44, amp[687]); 
  VVVV7_0(w[0], w[62], w[7], w[4], pars->GC_44, amp[688]); 
  VVVV1_0(w[0], w[62], w[7], w[4], pars->GC_8, amp[689]); 
  VVVV4_0(w[0], w[62], w[7], w[4], pars->GC_8, amp[690]); 
  VVVV5_0(w[0], w[62], w[7], w[4], pars->GC_8, amp[691]); 
  VVV1_0(w[7], w[4], w[165], pars->GC_6, amp[692]); 
  VVV2_0(w[7], w[4], w[164], pars->GC_28, amp[693]); 
  VVV1_0(w[7], w[4], w[164], pars->GC_6, amp[694]); 
  VVV1_0(w[62], w[4], w[160], pars->GC_6, amp[695]); 
  VVV2_0(w[62], w[4], w[156], pars->GC_28, amp[696]); 
  VVV1_0(w[62], w[4], w[156], pars->GC_6, amp[697]); 
  VVVV1_0(w[0], w[62], w[16], w[4], pars->GC_8, amp[698]); 
  VVVV4_0(w[0], w[62], w[16], w[4], pars->GC_8, amp[699]); 
  VVVV5_0(w[0], w[62], w[16], w[4], pars->GC_8, amp[700]); 
  VVV1_0(w[16], w[4], w[164], pars->GC_6, amp[701]); 
  VVV1_0(w[62], w[4], w[161], pars->GC_6, amp[702]); 
  FFS2_0(w[3], w[148], w[60], pars->GC_95, amp[703]); 
  FFV5_0(w[3], w[20], w[163], pars->GC_7, amp[704]); 
  FFV5_0(w[3], w[148], w[61], pars->GC_7, amp[705]); 
  FFVV1_2_0(w[3], w[20], w[0], w[62], pars->GC_146, pars->GC_85, amp[706]); 
  FFS2_0(w[3], w[20], w[166], pars->GC_95, amp[707]); 
  FFV5_0(w[3], w[20], w[165], pars->GC_7, amp[708]); 
  FFV5_0(w[3], w[20], w[164], pars->GC_7, amp[709]); 
  FFV3_8_0(w[3], w[20], w[164], pars->GC_145, pars->GC_79, amp[710]); 
  FFV5_0(w[3], w[148], w[62], pars->GC_7, amp[711]); 
  FFV3_8_0(w[3], w[148], w[62], pars->GC_145, pars->GC_79, amp[712]); 
  FFV5_0(w[3], w[149], w[62], pars->GC_7, amp[713]); 
  FFV5_0(w[3], w[26], w[164], pars->GC_7, amp[714]); 
  FFV5_0(w[3], w[150], w[62], pars->GC_7, amp[715]); 
  FFS2_0(w[135], w[2], w[60], pars->GC_95, amp[716]); 
  FFV5_0(w[32], w[2], w[163], pars->GC_7, amp[717]); 
  FFV5_0(w[135], w[2], w[61], pars->GC_7, amp[718]); 
  FFVV1_2_0(w[32], w[2], w[0], w[62], pars->GC_146, pars->GC_85, amp[719]); 
  FFS2_0(w[32], w[2], w[166], pars->GC_95, amp[720]); 
  FFV5_0(w[32], w[2], w[165], pars->GC_7, amp[721]); 
  FFV5_0(w[32], w[2], w[164], pars->GC_7, amp[722]); 
  FFV3_8_0(w[32], w[2], w[164], pars->GC_145, pars->GC_79, amp[723]); 
  FFV5_0(w[135], w[2], w[62], pars->GC_7, amp[724]); 
  FFV3_8_0(w[135], w[2], w[62], pars->GC_145, pars->GC_79, amp[725]); 
  FFV5_0(w[136], w[2], w[62], pars->GC_7, amp[726]); 
  FFV5_0(w[33], w[2], w[164], pars->GC_7, amp[727]); 
  FFV5_0(w[137], w[2], w[62], pars->GC_7, amp[728]); 
  VVV1_0(w[0], w[62], w[39], pars->GC_6, amp[729]); 
  VVVS1_0(w[0], w[1], w[17], w[11], pars->GC_84, amp[730]); 
  VVV1_0(w[1], w[17], w[158], pars->GC_6, amp[731]); 
  VVS2_0(w[1], w[145], w[11], pars->GC_78, amp[732]); 
  VVVV1_0(w[0], w[1], w[7], w[18], pars->GC_8, amp[733]); 
  VVVV4_0(w[0], w[1], w[7], w[18], pars->GC_8, amp[734]); 
  VVVV5_0(w[0], w[1], w[7], w[18], pars->GC_8, amp[735]); 
  VVV1_0(w[1], w[18], w[156], pars->GC_6, amp[736]); 
  VVV1_0(w[1], w[7], w[142], pars->GC_6, amp[737]); 
  VVVV2_0(w[0], w[1], w[7], w[17], pars->GC_44, amp[738]); 
  VVVV8_0(w[0], w[1], w[7], w[17], pars->GC_44, amp[739]); 
  VVVV7_0(w[0], w[1], w[7], w[17], pars->GC_44, amp[740]); 
  VVVV1_0(w[0], w[1], w[7], w[17], pars->GC_8, amp[741]); 
  VVVV4_0(w[0], w[1], w[7], w[17], pars->GC_8, amp[742]); 
  VVVV5_0(w[0], w[1], w[7], w[17], pars->GC_8, amp[743]); 
  VVV1_0(w[1], w[17], w[160], pars->GC_6, amp[744]); 
  VVV2_0(w[1], w[17], w[156], pars->GC_28, amp[745]); 
  VVV1_0(w[1], w[17], w[156], pars->GC_6, amp[746]); 
  VVV1_0(w[1], w[7], w[144], pars->GC_6, amp[747]); 
  VVV2_0(w[1], w[7], w[145], pars->GC_28, amp[748]); 
  VVV1_0(w[1], w[7], w[145], pars->GC_6, amp[749]); 
  VVVV1_0(w[0], w[1], w[16], w[17], pars->GC_8, amp[750]); 
  VVVV4_0(w[0], w[1], w[16], w[17], pars->GC_8, amp[751]); 
  VVVV5_0(w[0], w[1], w[16], w[17], pars->GC_8, amp[752]); 
  VVV1_0(w[1], w[17], w[161], pars->GC_6, amp[753]); 
  VVV1_0(w[1], w[16], w[145], pars->GC_6, amp[754]); 
  FFVV1_2_0(w[27], w[20], w[0], w[1], pars->GC_146, pars->GC_85, amp[755]); 
  FFV5_0(w[27], w[148], w[1], pars->GC_7, amp[756]); 
  FFV3_8_0(w[27], w[148], w[1], pars->GC_145, pars->GC_79, amp[757]); 
  FFV5_0(w[27], w[149], w[1], pars->GC_7, amp[758]); 
  FFV5_0(w[139], w[20], w[1], pars->GC_7, amp[759]); 
  FFV3_8_0(w[139], w[20], w[1], pars->GC_145, pars->GC_79, amp[760]); 
  FFV5_0(w[140], w[20], w[1], pars->GC_7, amp[761]); 
  FFV5_0(w[28], w[148], w[1], pars->GC_7, amp[762]); 
  FFV5_0(w[141], w[20], w[1], pars->GC_7, amp[763]); 
  FFV5_0(w[27], w[150], w[1], pars->GC_7, amp[764]); 
  FFV5_0(w[139], w[26], w[1], pars->GC_7, amp[765]); 
  FFVV1_2_0(w[32], w[29], w[0], w[1], pars->GC_146, pars->GC_85, amp[766]); 
  FFV5_0(w[32], w[152], w[1], pars->GC_7, amp[767]); 
  FFV3_8_0(w[32], w[152], w[1], pars->GC_145, pars->GC_79, amp[768]); 
  FFV5_0(w[32], w[153], w[1], pars->GC_7, amp[769]); 
  FFV5_0(w[135], w[29], w[1], pars->GC_7, amp[770]); 
  FFV3_8_0(w[135], w[29], w[1], pars->GC_145, pars->GC_79, amp[771]); 
  FFV5_0(w[136], w[29], w[1], pars->GC_7, amp[772]); 
  FFV5_0(w[33], w[152], w[1], pars->GC_7, amp[773]); 
  FFV5_0(w[137], w[29], w[1], pars->GC_7, amp[774]); 
  FFV5_0(w[32], w[154], w[1], pars->GC_7, amp[775]); 
  FFV5_0(w[135], w[31], w[1], pars->GC_7, amp[776]); 
  FFV5_0(w[32], w[167], w[5], pars->GC_7, amp[777]); 
  FFV5_0(w[27], w[167], w[4], pars->GC_7, amp[778]); 
  FFV5_0(w[3], w[167], w[17], pars->GC_7, amp[779]); 
  FFV5_0(w[168], w[20], w[5], pars->GC_7, amp[780]); 
  FFV5_0(w[168], w[29], w[4], pars->GC_7, amp[781]); 
  FFV5_0(w[168], w[2], w[17], pars->GC_7, amp[782]); 
  VVV1_0(w[169], w[7], w[5], pars->GC_6, amp[783]); 
  VVV1_0(w[170], w[7], w[5], pars->GC_6, amp[784]); 
  VVV1_0(w[171], w[7], w[5], pars->GC_6, amp[785]); 
  VVS2_0(w[172], w[5], w[11], pars->GC_78, amp[786]); 
  VVS2_0(w[173], w[5], w[11], pars->GC_78, amp[787]); 
  VVS2_0(w[174], w[5], w[11], pars->GC_78, amp[788]); 
  VVV2_0(w[172], w[7], w[5], pars->GC_28, amp[789]); 
  VVV2_0(w[173], w[7], w[5], pars->GC_28, amp[790]); 
  VVV2_0(w[174], w[7], w[5], pars->GC_28, amp[791]); 
  VVV1_0(w[172], w[7], w[5], pars->GC_6, amp[792]); 
  VVV1_0(w[173], w[7], w[5], pars->GC_6, amp[793]); 
  VVV1_0(w[174], w[7], w[5], pars->GC_6, amp[794]); 
  VVV1_0(w[172], w[16], w[5], pars->GC_6, amp[795]); 
  VVV1_0(w[173], w[16], w[5], pars->GC_6, amp[796]); 
  VVV1_0(w[174], w[16], w[5], pars->GC_6, amp[797]); 
  FFV5_0(w[3], w[29], w[169], pars->GC_7, amp[798]); 
  FFV5_0(w[3], w[29], w[170], pars->GC_7, amp[799]); 
  FFV5_0(w[3], w[29], w[171], pars->GC_7, amp[800]); 
  FFV5_0(w[3], w[29], w[172], pars->GC_7, amp[801]); 
  FFV5_0(w[3], w[29], w[173], pars->GC_7, amp[802]); 
  FFV5_0(w[3], w[29], w[174], pars->GC_7, amp[803]); 
  FFV3_8_0(w[3], w[29], w[172], pars->GC_145, pars->GC_79, amp[804]); 
  FFV3_8_0(w[3], w[29], w[173], pars->GC_145, pars->GC_79, amp[805]); 
  FFV3_8_0(w[3], w[29], w[174], pars->GC_145, pars->GC_79, amp[806]); 
  FFV5_0(w[3], w[31], w[172], pars->GC_7, amp[807]); 
  FFV5_0(w[3], w[31], w[173], pars->GC_7, amp[808]); 
  FFV5_0(w[3], w[31], w[174], pars->GC_7, amp[809]); 
  FFS2_0(w[3], w[29], w[175], pars->GC_95, amp[810]); 
  FFV5_0(w[27], w[2], w[169], pars->GC_7, amp[811]); 
  FFV5_0(w[27], w[2], w[170], pars->GC_7, amp[812]); 
  FFV5_0(w[27], w[2], w[171], pars->GC_7, amp[813]); 
  FFV5_0(w[27], w[2], w[172], pars->GC_7, amp[814]); 
  FFV5_0(w[27], w[2], w[173], pars->GC_7, amp[815]); 
  FFV5_0(w[27], w[2], w[174], pars->GC_7, amp[816]); 
  FFV3_8_0(w[27], w[2], w[172], pars->GC_145, pars->GC_79, amp[817]); 
  FFV3_8_0(w[27], w[2], w[173], pars->GC_145, pars->GC_79, amp[818]); 
  FFV3_8_0(w[27], w[2], w[174], pars->GC_145, pars->GC_79, amp[819]); 
  FFV5_0(w[28], w[2], w[172], pars->GC_7, amp[820]); 
  FFV5_0(w[28], w[2], w[173], pars->GC_7, amp[821]); 
  FFV5_0(w[28], w[2], w[174], pars->GC_7, amp[822]); 
  FFS2_0(w[27], w[2], w[175], pars->GC_95, amp[823]); 
  FFVV1_2_0(w[3], w[2], w[5], w[172], pars->GC_146, pars->GC_85, amp[824]); 
  FFVV1_2_0(w[3], w[2], w[5], w[173], pars->GC_146, pars->GC_85, amp[825]); 
  FFVV1_2_0(w[3], w[2], w[5], w[174], pars->GC_146, pars->GC_85, amp[826]); 
  VVV1_0(w[176], w[7], w[4], pars->GC_6, amp[827]); 
  VVV1_0(w[177], w[7], w[4], pars->GC_6, amp[828]); 
  VVV1_0(w[178], w[7], w[4], pars->GC_6, amp[829]); 
  VVS2_0(w[179], w[4], w[11], pars->GC_78, amp[830]); 
  VVS2_0(w[180], w[4], w[11], pars->GC_78, amp[831]); 
  VVS2_0(w[181], w[4], w[11], pars->GC_78, amp[832]); 
  VVV2_0(w[179], w[7], w[4], pars->GC_28, amp[833]); 
  VVV2_0(w[180], w[7], w[4], pars->GC_28, amp[834]); 
  VVV2_0(w[181], w[7], w[4], pars->GC_28, amp[835]); 
  VVV1_0(w[179], w[7], w[4], pars->GC_6, amp[836]); 
  VVV1_0(w[180], w[7], w[4], pars->GC_6, amp[837]); 
  VVV1_0(w[181], w[7], w[4], pars->GC_6, amp[838]); 
  VVV1_0(w[179], w[16], w[4], pars->GC_6, amp[839]); 
  VVV1_0(w[180], w[16], w[4], pars->GC_6, amp[840]); 
  VVV1_0(w[181], w[16], w[4], pars->GC_6, amp[841]); 
  FFV5_0(w[3], w[20], w[176], pars->GC_7, amp[842]); 
  FFV5_0(w[3], w[20], w[177], pars->GC_7, amp[843]); 
  FFV5_0(w[3], w[20], w[178], pars->GC_7, amp[844]); 
  FFV5_0(w[3], w[20], w[179], pars->GC_7, amp[845]); 
  FFV5_0(w[3], w[20], w[180], pars->GC_7, amp[846]); 
  FFV5_0(w[3], w[20], w[181], pars->GC_7, amp[847]); 
  FFV3_8_0(w[3], w[20], w[179], pars->GC_145, pars->GC_79, amp[848]); 
  FFV3_8_0(w[3], w[20], w[180], pars->GC_145, pars->GC_79, amp[849]); 
  FFV3_8_0(w[3], w[20], w[181], pars->GC_145, pars->GC_79, amp[850]); 
  FFV5_0(w[3], w[26], w[179], pars->GC_7, amp[851]); 
  FFV5_0(w[3], w[26], w[180], pars->GC_7, amp[852]); 
  FFV5_0(w[3], w[26], w[181], pars->GC_7, amp[853]); 
  FFS2_0(w[3], w[20], w[182], pars->GC_95, amp[854]); 
  FFV5_0(w[32], w[2], w[176], pars->GC_7, amp[855]); 
  FFV5_0(w[32], w[2], w[177], pars->GC_7, amp[856]); 
  FFV5_0(w[32], w[2], w[178], pars->GC_7, amp[857]); 
  FFV5_0(w[32], w[2], w[179], pars->GC_7, amp[858]); 
  FFV5_0(w[32], w[2], w[180], pars->GC_7, amp[859]); 
  FFV5_0(w[32], w[2], w[181], pars->GC_7, amp[860]); 
  FFV3_8_0(w[32], w[2], w[179], pars->GC_145, pars->GC_79, amp[861]); 
  FFV3_8_0(w[32], w[2], w[180], pars->GC_145, pars->GC_79, amp[862]); 
  FFV3_8_0(w[32], w[2], w[181], pars->GC_145, pars->GC_79, amp[863]); 
  FFV5_0(w[33], w[2], w[179], pars->GC_7, amp[864]); 
  FFV5_0(w[33], w[2], w[180], pars->GC_7, amp[865]); 
  FFV5_0(w[33], w[2], w[181], pars->GC_7, amp[866]); 
  FFS2_0(w[32], w[2], w[182], pars->GC_95, amp[867]); 
  FFVV1_2_0(w[3], w[2], w[4], w[179], pars->GC_146, pars->GC_85, amp[868]); 
  FFVV1_2_0(w[3], w[2], w[4], w[180], pars->GC_146, pars->GC_85, amp[869]); 
  FFVV1_2_0(w[3], w[2], w[4], w[181], pars->GC_146, pars->GC_85, amp[870]); 
  VVV1_0(w[183], w[56], w[5], pars->GC_6, amp[871]); 
  VVV1_0(w[183], w[62], w[4], pars->GC_6, amp[872]); 
  VVV1_0(w[183], w[1], w[17], pars->GC_6, amp[873]); 
  VVVV1_0(w[1], w[4], w[5], w[183], pars->GC_8, amp[874]); 
  VVVV4_0(w[1], w[4], w[5], w[183], pars->GC_8, amp[875]); 
  VVVV5_0(w[1], w[4], w[5], w[183], pars->GC_8, amp[876]); 
  FFV5_0(w[44], w[184], w[5], pars->GC_7, amp[877]); 
  FFV5_0(w[3], w[184], w[62], pars->GC_7, amp[878]); 
  FFV5_0(w[27], w[184], w[1], pars->GC_7, amp[879]); 
  FFV5_0(w[44], w[185], w[4], pars->GC_7, amp[880]); 
  FFV5_0(w[3], w[185], w[56], pars->GC_7, amp[881]); 
  FFV5_0(w[32], w[185], w[1], pars->GC_7, amp[882]); 
  FFV5_0(w[186], w[76], w[5], pars->GC_7, amp[883]); 
  FFV5_0(w[186], w[2], w[62], pars->GC_7, amp[884]); 
  FFV5_0(w[186], w[29], w[1], pars->GC_7, amp[885]); 
  FFV5_0(w[187], w[76], w[4], pars->GC_7, amp[886]); 
  FFV5_0(w[187], w[2], w[56], pars->GC_7, amp[887]); 
  FFV5_0(w[187], w[20], w[1], pars->GC_7, amp[888]); 
  FFV5_0(w[3], w[76], w[188], pars->GC_7, amp[889]); 
  FFV5_0(w[3], w[76], w[189], pars->GC_7, amp[890]); 
  FFV5_0(w[3], w[76], w[190], pars->GC_7, amp[891]); 
  FFV5_0(w[3], w[76], w[191], pars->GC_7, amp[892]); 
  FFV5_0(w[3], w[76], w[192], pars->GC_7, amp[893]); 
  FFV5_0(w[3], w[76], w[193], pars->GC_7, amp[894]); 
  FFV3_8_0(w[3], w[76], w[191], pars->GC_145, pars->GC_79, amp[895]); 
  FFV3_8_0(w[3], w[76], w[192], pars->GC_145, pars->GC_79, amp[896]); 
  FFV3_8_0(w[3], w[76], w[193], pars->GC_145, pars->GC_79, amp[897]); 
  FFV5_0(w[3], w[81], w[191], pars->GC_7, amp[898]); 
  FFV5_0(w[3], w[81], w[192], pars->GC_7, amp[899]); 
  FFV5_0(w[3], w[81], w[193], pars->GC_7, amp[900]); 
  FFS2_0(w[3], w[76], w[194], pars->GC_95, amp[901]); 
  FFV5_0(w[44], w[2], w[188], pars->GC_7, amp[902]); 
  FFV5_0(w[44], w[2], w[189], pars->GC_7, amp[903]); 
  FFV5_0(w[44], w[2], w[190], pars->GC_7, amp[904]); 
  FFV5_0(w[44], w[2], w[191], pars->GC_7, amp[905]); 
  FFV5_0(w[44], w[2], w[192], pars->GC_7, amp[906]); 
  FFV5_0(w[44], w[2], w[193], pars->GC_7, amp[907]); 
  FFV3_8_0(w[44], w[2], w[191], pars->GC_145, pars->GC_79, amp[908]); 
  FFV3_8_0(w[44], w[2], w[192], pars->GC_145, pars->GC_79, amp[909]); 
  FFV3_8_0(w[44], w[2], w[193], pars->GC_145, pars->GC_79, amp[910]); 
  FFV5_0(w[49], w[2], w[191], pars->GC_7, amp[911]); 
  FFV5_0(w[49], w[2], w[192], pars->GC_7, amp[912]); 
  FFV5_0(w[49], w[2], w[193], pars->GC_7, amp[913]); 
  FFS2_0(w[44], w[2], w[194], pars->GC_95, amp[914]); 
  VVV1_0(w[188], w[1], w[7], pars->GC_6, amp[915]); 
  VVV1_0(w[189], w[1], w[7], pars->GC_6, amp[916]); 
  VVV1_0(w[190], w[1], w[7], pars->GC_6, amp[917]); 
  VVS2_0(w[191], w[1], w[11], pars->GC_78, amp[918]); 
  VVS2_0(w[192], w[1], w[11], pars->GC_78, amp[919]); 
  VVS2_0(w[193], w[1], w[11], pars->GC_78, amp[920]); 
  VVV2_0(w[191], w[1], w[7], pars->GC_28, amp[921]); 
  VVV2_0(w[192], w[1], w[7], pars->GC_28, amp[922]); 
  VVV2_0(w[193], w[1], w[7], pars->GC_28, amp[923]); 
  VVV1_0(w[191], w[1], w[7], pars->GC_6, amp[924]); 
  VVV1_0(w[192], w[1], w[7], pars->GC_6, amp[925]); 
  VVV1_0(w[193], w[1], w[7], pars->GC_6, amp[926]); 
  VVV1_0(w[191], w[1], w[16], pars->GC_6, amp[927]); 
  VVV1_0(w[192], w[1], w[16], pars->GC_6, amp[928]); 
  VVV1_0(w[193], w[1], w[16], pars->GC_6, amp[929]); 
  FFVV1_2_0(w[3], w[2], w[1], w[191], pars->GC_146, pars->GC_85, amp[930]); 
  FFVV1_2_0(w[3], w[2], w[1], w[192], pars->GC_146, pars->GC_85, amp[931]); 
  FFVV1_2_0(w[3], w[2], w[1], w[193], pars->GC_146, pars->GC_85, amp[932]); 
  VVV1_0(w[0], w[113], w[17], pars->GC_6, amp[933]); 
  FFV5_0(w[27], w[92], w[0], pars->GC_7, amp[934]); 
  FFV5_0(w[32], w[93], w[0], pars->GC_7, amp[935]); 
  FFV5_0(w[66], w[29], w[0], pars->GC_7, amp[936]); 
  FFV5_0(w[67], w[20], w[0], pars->GC_7, amp[937]); 
  VVV1_0(w[0], w[68], w[7], pars->GC_6, amp[938]); 
  VVV1_0(w[0], w[69], w[7], pars->GC_6, amp[939]); 
  VVV1_0(w[0], w[70], w[7], pars->GC_6, amp[940]); 
  VVS2_0(w[0], w[71], w[11], pars->GC_78, amp[941]); 
  VVS2_0(w[0], w[72], w[11], pars->GC_78, amp[942]); 
  VVS2_0(w[0], w[73], w[11], pars->GC_78, amp[943]); 
  VVV2_0(w[0], w[71], w[7], pars->GC_28, amp[944]); 
  VVV2_0(w[0], w[72], w[7], pars->GC_28, amp[945]); 
  VVV2_0(w[0], w[73], w[7], pars->GC_28, amp[946]); 
  VVV1_0(w[0], w[71], w[7], pars->GC_6, amp[947]); 
  VVV1_0(w[0], w[72], w[7], pars->GC_6, amp[948]); 
  VVV1_0(w[0], w[73], w[7], pars->GC_6, amp[949]); 
  VVV1_0(w[0], w[71], w[16], pars->GC_6, amp[950]); 
  VVV1_0(w[0], w[72], w[16], pars->GC_6, amp[951]); 
  VVV1_0(w[0], w[73], w[16], pars->GC_6, amp[952]); 
  FFS2_0(w[3], w[2], w[195], pars->GC_95, amp[953]); 
  FFS2_0(w[3], w[2], w[196], pars->GC_95, amp[954]); 
  FFS2_0(w[3], w[2], w[197], pars->GC_95, amp[955]); 
  FFV5_0(w[3], w[2], w[198], pars->GC_7, amp[956]); 
  FFV5_0(w[3], w[2], w[199], pars->GC_7, amp[957]); 
  FFV5_0(w[3], w[2], w[200], pars->GC_7, amp[958]); 
  FFV5_0(w[3], w[2], w[201], pars->GC_7, amp[959]); 
  FFV5_0(w[3], w[2], w[202], pars->GC_7, amp[960]); 
  FFV5_0(w[3], w[2], w[203], pars->GC_7, amp[961]); 
  FFV5_0(w[3], w[2], w[204], pars->GC_7, amp[962]); 
  FFV5_0(w[3], w[2], w[205], pars->GC_7, amp[963]); 
  FFV5_0(w[3], w[2], w[206], pars->GC_7, amp[964]); 
  FFV5_0(w[3], w[2], w[207], pars->GC_7, amp[965]); 
  FFV5_0(w[3], w[2], w[208], pars->GC_7, amp[966]); 
  FFV5_0(w[3], w[2], w[209], pars->GC_7, amp[967]); 
  FFV5_0(w[3], w[2], w[210], pars->GC_7, amp[968]); 
  FFV5_0(w[3], w[2], w[211], pars->GC_7, amp[969]); 
  FFV5_0(w[3], w[2], w[212], pars->GC_7, amp[970]); 
  FFV5_0(w[3], w[2], w[213], pars->GC_7, amp[971]); 
  FFV5_0(w[3], w[2], w[214], pars->GC_7, amp[972]); 
  FFV5_0(w[3], w[2], w[215], pars->GC_7, amp[973]); 
  FFV5_0(w[3], w[2], w[216], pars->GC_7, amp[974]); 
  FFV5_0(w[3], w[2], w[217], pars->GC_7, amp[975]); 
  FFV5_0(w[3], w[2], w[218], pars->GC_7, amp[976]); 
  FFV5_0(w[3], w[2], w[219], pars->GC_7, amp[977]); 
  FFV5_0(w[3], w[2], w[220], pars->GC_7, amp[978]); 
  FFV5_0(w[3], w[2], w[221], pars->GC_7, amp[979]); 
  FFV5_0(w[3], w[2], w[222], pars->GC_7, amp[980]); 
  FFV5_0(w[3], w[2], w[223], pars->GC_7, amp[981]); 
  FFV5_0(w[3], w[2], w[224], pars->GC_7, amp[982]); 
  FFV5_0(w[3], w[2], w[225], pars->GC_7, amp[983]); 
  FFV5_0(w[3], w[2], w[226], pars->GC_7, amp[984]); 
  FFV5_0(w[3], w[2], w[227], pars->GC_7, amp[985]); 

}
double CPPProcess::matrix_3_gg_ttxgg() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 986; 
  const int ncolor = 50; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {18, 54, 54, 18, 18, 54, 18, 54, 54, 18,
      54, 18, 54, 54, 18, 18, 54, 18, 54, 54, 18, 54, 6, 6, 6, 6, 6, 6, 54, 18,
      54, 54, 18, 54, 18, 18, 54, 54, 18, 54, 18, 54, 54, 18, 54, 18, 18, 54,
      54, 18};
  static const double cf[ncolor][ncolor] = {{192, 64, 64, 0, 0, -8, 24, 64, -8,
      24, 64, -24, -8, -8, 0, 0, 1, -3, -8, 1, -3, -8, 24, 24, -3, 24, -3, 24,
      1, -3, -8, 10, 24, 1, 21, -6, 64, -8, 24, 1, -3, -8, 10, 24, 1, 21, -6,
      64, -8, 24}, {192, 512, -64, 168, -48, -64, -24, 8, 8, 192, 80, -24, -64,
      8, -21, 6, 8, 3, -1, -1, -24, -10, 171, -18, -18, -18, -18, 36, 8, 3, -1,
      80, 30, -10, -21, 6, 71, 62, 192, -1, -24, -10, -10, 3, 62, 168, -48, 62,
      -28, -24}, {192, -64, 512, -48, 168, 8, 192, 80, -64, -24, 8, -24, 8,
      -64, 6, -21, -1, -24, -10, 8, 3, -1, -18, 171, -18, 36, -18, -18, -1,
      -24, -10, -10, 3, 62, 168, -48, 62, -28, -24, 8, 3, -1, 80, 30, -10, -21,
      6, 71, 62, 192}, {0, 56, -16, 168, -48, -16, 0, 56, 56, 0, -16, 0, -7, 2,
      21, -6, -7, 21, 56, 2, -6, -16, 0, 0, 0, 0, 0, 0, 2, 0, -7, 2, -6, -16,
      -6, 21, -7, 56, 21, -7, 0, 2, -7, 21, 56, 21, -6, 2, -16, -6}, {0, -16,
      56, -48, 168, 56, 0, -16, -16, 0, 56, 0, 2, -7, -6, 21, 2, -6, -16, -7,
      21, 56, 0, 0, 0, 0, 0, 0, -7, 0, 2, -7, 21, 56, 21, -6, 2, -16, -6, 2, 0,
      -7, 2, -6, -16, -6, 21, -7, 56, 21}, {-24, -64, 8, -48, 168, 512, 192,
      -64, 80, 192, 8, 3, 8, -1, -21, 6, 80, 30, -10, 71, 192, 62, -18, -18,
      171, -18, 36, -18, -64, -24, 8, 8, 3, -1, -21, 6, -1, -10, -24, -10, -24,
      -1, 62, -24, -28, -48, 168, -10, 62, 3}, {24, -8, 64, 0, 0, 64, 192, 64,
      64, 24, -8, -3, 1, -8, 21, -6, 10, 24, 1, 64, 24, -8, -3, 24, 24, 24, 24,
      -3, -8, -24, -8, 1, -3, -8, 0, 0, 1, -8, -3, -8, -3, 1, 64, 24, -8, -6,
      21, 10, 1, 24}, {192, 8, 80, 168, -48, -64, 192, 512, 8, -24, -64, -24,
      -1, -10, 168, -48, -10, 3, 62, 62, -24, -28, -18, 36, -18, 171, -18, -18,
      8, -24, -64, -1, -24, -10, 6, -21, 8, -1, 3, -1, 3, 8, 71, 192, 62, 6,
      -21, 80, -10, 30}, {-24, 8, -64, 168, -48, 80, 192, 8, 512, 192, -64, 3,
      -1, 8, 6, -21, 71, 192, 62, 80, 30, -10, -18, -18, 36, -18, 171, -18,
      -10, -24, -1, 62, -24, -28, -48, 168, -10, 62, 3, -64, -24, 8, 8, 3, -1,
      -21, 6, -1, -10, -24}, {24, 64, -8, 0, 0, 64, 24, -8, 64, 192, 64, -3,
      -8, 1, -6, 21, 64, 24, -8, 10, 24, 1, 24, -3, 24, -3, 24, 24, -8, -3, 1,
      64, 24, -8, -6, 21, 10, 1, 24, -8, -24, -8, 1, -3, -8, 0, 0, 1, -8, -3},
      {192, 80, 8, -48, 168, 8, -24, -64, -64, 192, 512, -24, -10, -1, -48,
      168, 62, -24, -28, -10, 3, 62, 36, -18, -18, -18, -18, 171, -1, 3, 8, 71,
      192, 62, 6, -21, 80, -10, 30, 8, -24, -64, -1, -24, -10, 6, -21, 8, -1,
      3}, {-24, -8, -8, 0, 0, 1, -3, -8, 1, -3, -8, 192, 64, 64, 0, 0, -8, 24,
      64, -8, 24, 64, 24, 24, -3, 24, -3, 24, 10, 24, 1, 1, -3, -8, -6, 21, -8,
      64, 24, 10, 24, 1, 1, -3, -8, -6, 21, -8, 64, 24}, {-24, -64, 8, -21, 6,
      8, 3, -1, -1, -24, -10, 192, 512, -64, 168, -48, -64, -24, 8, 8, 192, 80,
      -18, 36, -18, 171, -18, -18, 80, 30, -10, 8, 3, -1, 6, -21, 62, 71, 192,
      -10, 3, 62, -1, -24, -10, -48, 168, -28, 62, -24}, {-24, 8, -64, 6, -21,
      -1, -24, -10, 8, 3, -1, 192, -64, 512, -48, 168, 8, 192, 80, -64, -24, 8,
      36, -18, -18, -18, -18, 171, -10, 3, 62, -1, -24, -10, -48, 168, -28, 62,
      -24, 80, 30, -10, 8, 3, -1, 6, -21, 62, 71, 192}, {0, -7, 2, 21, -6, -7,
      21, 56, 2, -6, -16, 0, 56, -16, 168, -48, -16, 0, 56, 56, 0, -16, 0, 0,
      0, 0, 0, 0, 2, -6, -16, 2, 0, -7, 21, -6, 56, -7, 21, -7, 21, 56, -7, 0,
      2, -6, 21, -16, 2, -6}, {0, 2, -7, -6, 21, 2, -6, -16, -7, 21, 56, 0,
      -16, 56, -48, 168, 56, 0, -16, -16, 0, 56, 0, 0, 0, 0, 0, 0, -7, 21, 56,
      -7, 0, 2, -6, 21, -16, 2, -6, 2, -6, -16, 2, 0, -7, 21, -6, 56, -7, 21},
      {3, 8, -1, -21, 6, 80, 30, -10, 71, 192, 62, -24, -64, 8, -48, 168, 512,
      192, -64, 80, 192, 8, -18, -18, 36, -18, 171, -18, 8, 3, -1, -64, -24, 8,
      6, -21, -10, -1, -24, 62, -24, -28, -10, -24, -1, 168, -48, 62, -10, 3},
      {-3, 1, -8, 21, -6, 10, 24, 1, 64, 24, -8, 24, -8, 64, 0, 0, 64, 192, 64,
      64, 24, -8, 24, -3, 24, -3, 24, 24, 1, -3, -8, -8, -24, -8, 0, 0, -8, 1,
      -3, 64, 24, -8, -8, -3, 1, 21, -6, 1, 10, 24}, {-24, -1, -10, 168, -48,
      -10, 3, 62, 62, -24, -28, 192, 8, 80, 168, -48, -64, 192, 512, 8, -24,
      -64, 171, -18, -18, -18, -18, 36, -1, -24, -10, 8, -24, -64, -21, 6, -1,
      8, 3, 71, 192, 62, -1, 3, 8, -21, 6, -10, 80, 30}, {3, -1, 8, 6, -21, 71,
      192, 62, 80, 30, -10, -24, 8, -64, 168, -48, 80, 192, 8, 512, 192, -64,
      -18, -18, 171, -18, 36, -18, 62, -24, -28, -10, -24, -1, 168, -48, 62,
      -10, 3, 8, 3, -1, -64, -24, 8, 6, -21, -10, -1, -24}, {-3, -8, 1, -6, 21,
      64, 24, -8, 10, 24, 1, 24, 64, -8, 0, 0, 64, 24, -8, 64, 192, 64, -3, 24,
      24, 24, 24, -3, 64, 24, -8, -8, -3, 1, 21, -6, 1, 10, 24, 1, -3, -8, -8,
      -24, -8, 0, 0, -8, 1, -3}, {-24, -10, -1, -48, 168, 62, -24, -28, -10, 3,
      62, 192, 80, 8, -48, 168, 8, -24, -64, -64, 192, 512, -18, 171, -18, 36,
      -18, -18, 71, 192, 62, -1, 3, 8, -21, 6, -10, 80, 30, -1, -24, -10, 8,
      -24, -64, -21, 6, -1, 8, 3}, {8, 19, -2, 0, 0, -2, -1, -2, -2, 8, 4, 8,
      -2, 4, 0, 0, -2, 8, 19, -2, -1, -2, 57, -6, -6, -6, -6, 12, -2, -1, -2,
      4, 8, -2, 0, 0, 19, -2, 8, 19, 8, -2, -2, -1, -2, 0, 0, -2, 4, 8}, {8,
      -2, 19, 0, 0, -2, 8, 4, -2, -1, -2, 8, 4, -2, 0, 0, -2, -1, -2, -2, 8,
      19, -6, 57, -6, 12, -6, -6, 19, 8, -2, -2, -1, -2, 0, 0, -2, 4, 8, -2,
      -1, -2, 4, 8, -2, 0, 0, 19, -2, 8}, {-1, -2, -2, 0, 0, 19, 8, -2, 4, 8,
      -2, -1, -2, -2, 0, 0, 4, 8, -2, 19, 8, -2, -6, -6, 57, -6, 12, -6, -2, 8,
      4, -2, 8, 19, 0, 0, -2, -2, -1, -2, 8, 19, -2, 8, 4, 0, 0, -2, -2, -1},
      {8, -2, 4, 0, 0, -2, 8, 19, -2, -1, -2, 8, 19, -2, 0, 0, -2, -1, -2, -2,
      8, 4, -6, 12, -6, 57, -6, -6, 4, 8, -2, -2, -1, -2, 0, 0, -2, 19, 8, -2,
      -1, -2, 19, 8, -2, 0, 0, 4, -2, 8}, {-1, -2, -2, 0, 0, 4, 8, -2, 19, 8,
      -2, -1, -2, -2, 0, 0, 19, 8, -2, 4, 8, -2, -6, -6, 12, -6, 57, -6, -2, 8,
      19, -2, 8, 4, 0, 0, -2, -2, -1, -2, 8, 4, -2, 8, 19, 0, 0, -2, -2, -1},
      {8, 4, -2, 0, 0, -2, -1, -2, -2, 8, 19, 8, -2, 19, 0, 0, -2, 8, 4, -2,
      -1, -2, 12, -6, -6, -6, -6, 57, -2, -1, -2, 19, 8, -2, 0, 0, 4, -2, 8, 4,
      8, -2, -2, -1, -2, 0, 0, -2, 19, 8}, {3, 8, -1, 6, -21, -64, -24, 8, -10,
      -24, -1, 30, 80, -10, 6, -21, 8, 3, -1, 62, 192, 71, -18, 171, -18, 36,
      -18, -18, 512, 192, -64, -64, -24, 8, 168, -48, 8, 80, 192, 62, 3, -10,
      -28, -24, 62, 168, -48, -1, -10, -24}, {-3, 1, -8, 0, 0, -8, -24, -8, -8,
      -3, 1, 24, 10, 1, -6, 21, 1, -3, -8, -8, 24, 64, -3, 24, 24, 24, 24, -3,
      64, 192, 64, -8, 24, 64, 0, 0, -8, 64, 24, 1, 24, 10, -8, 24, 64, 21, -6,
      1, -8, -3}, {-24, -1, -10, -21, 6, 8, -24, -64, -1, 3, 8, 3, -10, 62,
      -48, 168, -1, -24, -10, -28, -24, 62, -18, -18, 36, -18, 171, -18, -64,
      192, 512, 8, 192, 80, -48, 168, -64, 8, -24, -10, 30, 80, 62, 192, 71,
      -21, 6, 8, -1, 3}, {30, 80, -10, 6, -21, 8, 3, -1, 62, 192, 71, 3, 8, -1,
      6, -21, -64, -24, 8, -10, -24, -1, 36, -18, -18, -18, -18, 171, -64, -24,
      8, 512, 192, -64, -48, 168, 80, 8, 192, -28, -24, 62, 62, 3, -10, -48,
      168, -10, -1, -24}, {24, 10, 1, -6, 21, 1, -3, -8, -8, 24, 64, -3, 1, -8,
      0, 0, -8, -24, -8, -8, -3, 1, 24, -3, 24, -3, 24, 24, -8, 24, 64, 64,
      192, 64, 0, 0, 64, -8, 24, -8, 24, 64, 1, 24, 10, -6, 21, -8, 1, -3}, {3,
      -10, 62, -48, 168, -1, -24, -10, -28, -24, 62, -24, -1, -10, -21, 6, 8,
      -24, -64, -1, 3, 8, -18, -18, 171, -18, 36, -18, 8, 192, 80, -64, 192,
      512, 168, -48, 8, -64, -24, 62, 192, 71, -10, 30, 80, 6, -21, -1, 8, 3},
      {21, -7, 56, -6, 21, -7, 0, 2, -16, -6, 2, -6, 2, -16, 21, -6, 2, 0, -7,
      56, 21, -7, 0, 0, 0, 0, 0, 0, 56, 0, -16, -16, 0, 56, 168, -48, 56, -16,
      0, 56, 21, -7, -16, -6, 2, 21, -6, -7, 2, 0}, {-6, 2, -16, 21, -6, 2, 0,
      -7, 56, 21, -7, 21, -7, 56, -6, 21, -7, 0, 2, -16, -6, 2, 0, 0, 0, 0, 0,
      0, -16, 0, 56, 56, 0, -16, -48, 168, -16, 56, 0, -16, -6, 2, 56, 21, -7,
      -6, 21, 2, -7, 0}, {192, 71, 62, -21, 6, -1, 3, 8, -10, 30, 80, -24, 62,
      -28, 168, -48, -10, -24, -1, 62, 3, -10, 171, -18, -18, -18, -18, 36, 8,
      -24, -64, 80, 192, 8, 168, -48, 512, -64, 192, -1, 3, 8, -10, -24, -1,
      -21, 6, -64, 8, -24}, {-24, 62, -28, 168, -48, -10, -24, -1, 62, 3, -10,
      192, 71, 62, -21, 6, -1, 3, 8, -10, 30, 80, -18, 36, -18, 171, -18, -18,
      80, 192, 8, 8, -24, -64, -48, 168, -64, 512, 192, -10, -24, -1, -1, 3, 8,
      6, -21, 8, -64, -24}, {24, 64, -8, 21, -6, -8, -3, 1, 1, 24, 10, 24, 64,
      -8, 21, -6, -8, -3, 1, 1, 24, 10, 24, 24, -3, 24, -3, 24, 64, 24, -8, 64,
      24, -8, 0, 0, 64, 64, 192, -8, -3, 1, -8, -3, 1, 0, 0, -8, -8, -24}, {3,
      -1, 8, -21, 6, -10, -24, -1, -64, -24, 8, 30, -10, 80, -21, 6, 62, 192,
      71, 8, 3, -1, 171, -18, -18, -18, -18, 36, 62, 3, -10, -28, -24, 62, 168,
      -48, -1, -10, -24, 512, 192, -64, -64, -24, 8, 168, -48, 8, 80, 192},
      {-3, -8, 1, 0, 0, -8, -3, 1, -8, -24, -8, 24, 1, 10, 21, -6, -8, 24, 64,
      1, -3, -8, 24, -3, 24, -3, 24, 24, 1, 24, 10, -8, 24, 64, 21, -6, 1, -8,
      -3, 64, 192, 64, -8, 24, 64, 0, 0, -8, 64, 24}, {-24, -10, -1, 6, -21,
      -1, 3, 8, 8, -24, -64, 3, 62, -10, 168, -48, -28, -24, 62, -1, -24, -10,
      -18, -18, 171, -18, 36, -18, -10, 30, 80, 62, 192, 71, -21, 6, 8, -1, 3,
      -64, 192, 512, 8, 192, 80, -48, 168, -64, 8, -24}, {30, -10, 80, -21, 6,
      62, 192, 71, 8, 3, -1, 3, -1, 8, -21, 6, -10, -24, -1, -64, -24, 8, -18,
      36, -18, 171, -18, -18, -28, -24, 62, 62, 3, -10, -48, 168, -10, -1, -24,
      -64, -24, 8, 512, 192, -64, -48, 168, 80, 8, 192}, {24, 1, 10, 21, -6,
      -8, 24, 64, 1, -3, -8, -3, -8, 1, 0, 0, -8, -3, 1, -8, -24, -8, -3, 24,
      24, 24, 24, -3, -8, 24, 64, 1, 24, 10, -6, 21, -8, 1, -3, -8, 24, 64, 64,
      192, 64, 0, 0, 64, -8, 24}, {3, 62, -10, 168, -48, -28, -24, 62, -1, -24,
      -10, -24, -10, -1, 6, -21, -1, 3, 8, 8, -24, -64, -18, -18, 36, -18, 171,
      -18, 62, 192, 71, -10, 30, 80, 6, -21, -1, 8, 3, 8, 192, 80, -64, 192,
      512, 168, -48, 8, -64, -24}, {21, 56, -7, 21, -6, -16, -6, 2, -7, 0, 2,
      -6, -16, 2, -6, 21, 56, 21, -7, 2, 0, -7, 0, 0, 0, 0, 0, 0, 56, 21, -7,
      -16, -6, 2, 21, -6, -7, 2, 0, 56, 0, -16, -16, 0, 56, 168, -48, 56, -16,
      0}, {-6, -16, 2, -6, 21, 56, 21, -7, 2, 0, -7, 21, 56, -7, 21, -6, -16,
      -6, 2, -7, 0, 2, 0, 0, 0, 0, 0, 0, -16, -6, 2, 56, 21, -7, -6, 21, 2, -7,
      0, -16, 0, 56, 56, 0, -16, -48, 168, -16, 56, 0}, {192, 62, 71, 6, -21,
      -10, 30, 80, -1, 3, 8, -24, -28, 62, -48, 168, 62, 3, -10, -10, -24, -1,
      -18, 171, -18, 36, -18, -18, -1, 3, 8, -10, -24, -1, -21, 6, -64, 8, -24,
      8, -24, -64, 80, 192, 8, 168, -48, 512, -64, 192}, {-24, -28, 62, -48,
      168, 62, 3, -10, -10, -24, -1, 192, 62, 71, 6, -21, -10, 30, 80, -1, 3,
      8, 36, -18, -18, -18, -18, 171, -10, -24, -1, -1, 3, 8, 6, -21, 8, -64,
      -24, 80, 192, 8, 8, -24, -64, -48, 168, -64, 512, 192}, {24, -8, 64, -6,
      21, 1, 24, 10, -8, -3, 1, 24, -8, 64, -6, 21, 1, 24, 10, -8, -3, 1, 24,
      24, -3, 24, -3, 24, -8, -3, 1, -8, -3, 1, 0, 0, -8, -8, -24, 64, 24, -8,
      64, 24, -8, 0, 0, 64, 64, 192}};

  // Calculate color flows
  jamp[0] = +2. * (+std::complex<double> (0, 1) * amp[99] +
      std::complex<double> (0, 1) * amp[100] - amp[125] - amp[191] - amp[612]);
  jamp[1] = +std::complex<double> (0, 1) * amp[0] + std::complex<double> (0, 1)
      * amp[1] + std::complex<double> (0, 1) * amp[3] + std::complex<double>
      (0, 1) * amp[8] + std::complex<double> (0, 1) * amp[9] +
      std::complex<double> (0, 1) * amp[11] + std::complex<double> (0, 1) *
      amp[12] + std::complex<double> (0, 1) * amp[14] + std::complex<double>
      (0, 1) * amp[15] + std::complex<double> (0, 1) * amp[16] +
      std::complex<double> (0, 1) * amp[20] + std::complex<double> (0, 1) *
      amp[21] + std::complex<double> (0, 1) * amp[23] + std::complex<double>
      (0, 1) * amp[25] + std::complex<double> (0, 1) * amp[27] +
      std::complex<double> (0, 1) * amp[28] + std::complex<double> (0, 1) *
      amp[29] + std::complex<double> (0, 1) * amp[30] + std::complex<double>
      (0, 1) * amp[83] + amp[84] + amp[85] + std::complex<double> (0, 1) *
      amp[86] + std::complex<double> (0, 1) * amp[87] + std::complex<double>
      (0, 1) * amp[88] + amp[90] + amp[91] + amp[92] + std::complex<double> (0,
      1) * amp[93] + amp[94] + amp[97] + amp[101] - std::complex<double> (0, 1)
      * amp[103] + amp[104] + amp[105] + amp[106] + std::complex<double> (0, 1)
      * amp[111] + amp[113] + amp[132] + amp[134] + amp[136] + amp[137] +
      amp[138] + amp[142] + std::complex<double> (0, 1) * amp[145] +
      std::complex<double> (0, 1) * amp[146] + std::complex<double> (0, 1) *
      amp[147] + std::complex<double> (0, 1) * amp[148] + std::complex<double>
      (0, 1) * amp[149] + std::complex<double> (0, 1) * amp[180] - amp[181] -
      amp[182] - amp[183] - amp[187] - amp[189] + std::complex<double> (0, 1) *
      amp[192] + amp[193] + amp[194] + std::complex<double> (0, 1) * amp[195] +
      std::complex<double> (0, 1) * amp[196] + std::complex<double> (0, 1) *
      amp[197] + amp[199] + amp[200] + amp[201] + std::complex<double> (0, 1) *
      amp[202] + amp[203] - amp[208] + amp[206] - amp[211] + amp[209] -
      amp[214] + amp[212] - amp[218] + amp[216] + std::complex<double> (0, 1) *
      amp[219] + std::complex<double> (0, 1) * amp[626] + std::complex<double>
      (0, 1) * amp[627] + std::complex<double> (0, 1) * amp[629] +
      std::complex<double> (0, 1) * amp[630] + std::complex<double> (0, 1) *
      amp[634] + std::complex<double> (0, 1) * amp[635] + std::complex<double>
      (0, 1) * amp[637] + std::complex<double> (0, 1) * amp[638] +
      std::complex<double> (0, 1) * amp[640] + std::complex<double> (0, 1) *
      amp[641] + std::complex<double> (0, 1) * amp[642] + std::complex<double>
      (0, 1) * amp[643] + std::complex<double> (0, 1) * amp[644] +
      std::complex<double> (0, 1) * amp[645] + std::complex<double> (0, 1) *
      amp[646] + std::complex<double> (0, 1) * amp[647] + std::complex<double>
      (0, 1) * amp[649] + std::complex<double> (0, 1) * amp[650] + amp[665] +
      amp[667] + amp[669] + amp[670] + amp[671] + amp[675] -
      std::complex<double> (0, 1) * amp[677] + std::complex<double> (0, 1) *
      amp[733] + std::complex<double> (0, 1) * amp[734] + std::complex<double>
      (0, 1) * amp[736] + std::complex<double> (0, 1) * amp[738] +
      std::complex<double> (0, 1) * amp[739] + std::complex<double> (0, 1) *
      amp[741] + std::complex<double> (0, 1) * amp[742] + std::complex<double>
      (0, 1) * amp[744] + std::complex<double> (0, 1) * amp[745] +
      std::complex<double> (0, 1) * amp[746] + std::complex<double> (0, 1) *
      amp[750] + std::complex<double> (0, 1) * amp[751] + std::complex<double>
      (0, 1) * amp[753] + std::complex<double> (0, 1) * amp[778] + amp[779] +
      std::complex<double> (0, 1) * amp[783] - std::complex<double> (0, 1) *
      amp[785] + std::complex<double> (0, 1) * amp[789] - std::complex<double>
      (0, 1) * amp[791] + std::complex<double> (0, 1) * amp[792] -
      std::complex<double> (0, 1) * amp[794] + std::complex<double> (0, 1) *
      amp[795] - std::complex<double> (0, 1) * amp[797] - amp[813] + amp[811] -
      amp[816] + amp[814] - amp[819] + amp[817] - amp[822] + amp[820] -
      std::complex<double> (0, 1) * amp[826] + std::complex<double> (0, 1) *
      amp[824] + std::complex<double> (0, 1) * amp[871] + std::complex<double>
      (0, 1) * amp[873] + std::complex<double> (0, 1) * amp[874] -
      std::complex<double> (0, 1) * amp[876] - std::complex<double> (0, 1) *
      amp[938] + std::complex<double> (0, 1) * amp[940] - std::complex<double>
      (0, 1) * amp[944] + std::complex<double> (0, 1) * amp[946] -
      std::complex<double> (0, 1) * amp[947] + std::complex<double> (0, 1) *
      amp[949] - std::complex<double> (0, 1) * amp[950] + std::complex<double>
      (0, 1) * amp[952] + std::complex<double> (0, 1) * amp[985] +
      std::complex<double> (0, 1) * amp[965] - std::complex<double> (0, 1) *
      amp[956] - std::complex<double> (0, 1) * amp[967] - std::complex<double>
      (0, 1) * amp[978] - std::complex<double> (0, 1) * amp[983] +
      std::complex<double> (0, 1) * amp[968] + std::complex<double> (0, 1) *
      amp[958] - std::complex<double> (0, 1) * amp[979] + std::complex<double>
      (0, 1) * amp[969];
  jamp[2] = -std::complex<double> (0, 1) * amp[0] + std::complex<double> (0, 1)
      * amp[2] + std::complex<double> (0, 1) * amp[4] - std::complex<double>
      (0, 1) * amp[8] + std::complex<double> (0, 1) * amp[10] -
      std::complex<double> (0, 1) * amp[11] + std::complex<double> (0, 1) *
      amp[13] + std::complex<double> (0, 1) * amp[17] + std::complex<double>
      (0, 1) * amp[18] + std::complex<double> (0, 1) * amp[19] -
      std::complex<double> (0, 1) * amp[20] + std::complex<double> (0, 1) *
      amp[22] + std::complex<double> (0, 1) * amp[24] - std::complex<double>
      (0, 1) * amp[25] - std::complex<double> (0, 1) * amp[27] -
      std::complex<double> (0, 1) * amp[28] - std::complex<double> (0, 1) *
      amp[29] - std::complex<double> (0, 1) * amp[30] + std::complex<double>
      (0, 1) * amp[70] + amp[71] + amp[72] + std::complex<double> (0, 1) *
      amp[73] + std::complex<double> (0, 1) * amp[74] + std::complex<double>
      (0, 1) * amp[75] + amp[77] + amp[78] + amp[79] + std::complex<double> (0,
      1) * amp[80] + amp[81] - amp[97] - amp[101] + std::complex<double> (0, 1)
      * amp[103] - amp[104] - amp[105] - amp[106] + std::complex<double> (0, 1)
      * amp[110] - amp[113] + amp[151] + amp[153] + amp[155] + amp[156] +
      amp[157] + amp[161] + std::complex<double> (0, 1) * amp[164] +
      std::complex<double> (0, 1) * amp[165] + std::complex<double> (0, 1) *
      amp[166] + std::complex<double> (0, 1) * amp[167] + std::complex<double>
      (0, 1) * amp[168] + std::complex<double> (0, 1) * amp[169] - amp[170] -
      amp[171] - amp[172] - amp[176] - amp[178] - std::complex<double> (0, 1) *
      amp[192] - amp[193] - amp[194] - std::complex<double> (0, 1) * amp[195] -
      std::complex<double> (0, 1) * amp[196] - std::complex<double> (0, 1) *
      amp[197] - amp[199] - amp[200] - amp[201] - std::complex<double> (0, 1) *
      amp[202] - amp[203] - amp[207] - amp[206] - amp[210] - amp[209] -
      amp[213] - amp[212] - amp[217] - amp[216] - std::complex<double> (0, 1) *
      amp[219] + std::complex<double> (0, 1) * amp[678] + std::complex<double>
      (0, 1) * amp[679] + std::complex<double> (0, 1) * amp[681] +
      std::complex<double> (0, 1) * amp[682] + std::complex<double> (0, 1) *
      amp[686] + std::complex<double> (0, 1) * amp[687] + std::complex<double>
      (0, 1) * amp[689] + std::complex<double> (0, 1) * amp[690] +
      std::complex<double> (0, 1) * amp[692] + std::complex<double> (0, 1) *
      amp[693] + std::complex<double> (0, 1) * amp[694] + std::complex<double>
      (0, 1) * amp[695] + std::complex<double> (0, 1) * amp[696] +
      std::complex<double> (0, 1) * amp[697] + std::complex<double> (0, 1) *
      amp[698] + std::complex<double> (0, 1) * amp[699] + std::complex<double>
      (0, 1) * amp[701] + std::complex<double> (0, 1) * amp[702] + amp[717] +
      amp[719] + amp[721] + amp[722] + amp[723] + amp[727] -
      std::complex<double> (0, 1) * amp[729] - std::complex<double> (0, 1) *
      amp[733] - std::complex<double> (0, 1) * amp[734] - std::complex<double>
      (0, 1) * amp[736] - std::complex<double> (0, 1) * amp[738] -
      std::complex<double> (0, 1) * amp[739] - std::complex<double> (0, 1) *
      amp[741] - std::complex<double> (0, 1) * amp[742] - std::complex<double>
      (0, 1) * amp[744] - std::complex<double> (0, 1) * amp[745] -
      std::complex<double> (0, 1) * amp[746] - std::complex<double> (0, 1) *
      amp[750] - std::complex<double> (0, 1) * amp[751] - std::complex<double>
      (0, 1) * amp[753] + std::complex<double> (0, 1) * amp[777] - amp[779] +
      std::complex<double> (0, 1) * amp[827] - std::complex<double> (0, 1) *
      amp[829] + std::complex<double> (0, 1) * amp[833] - std::complex<double>
      (0, 1) * amp[835] + std::complex<double> (0, 1) * amp[836] -
      std::complex<double> (0, 1) * amp[838] + std::complex<double> (0, 1) *
      amp[839] - std::complex<double> (0, 1) * amp[841] - amp[857] + amp[855] -
      amp[860] + amp[858] - amp[863] + amp[861] - amp[866] + amp[864] -
      std::complex<double> (0, 1) * amp[870] + std::complex<double> (0, 1) *
      amp[868] + std::complex<double> (0, 1) * amp[872] - std::complex<double>
      (0, 1) * amp[873] + std::complex<double> (0, 1) * amp[875] +
      std::complex<double> (0, 1) * amp[876] + std::complex<double> (0, 1) *
      amp[939] + std::complex<double> (0, 1) * amp[938] + std::complex<double>
      (0, 1) * amp[945] + std::complex<double> (0, 1) * amp[944] +
      std::complex<double> (0, 1) * amp[948] + std::complex<double> (0, 1) *
      amp[947] + std::complex<double> (0, 1) * amp[951] + std::complex<double>
      (0, 1) * amp[950] - std::complex<double> (0, 1) * amp[964] -
      std::complex<double> (0, 1) * amp[985] + std::complex<double> (0, 1) *
      amp[970] - std::complex<double> (0, 1) * amp[981] + std::complex<double>
      (0, 1) * amp[956] - std::complex<double> (0, 1) * amp[982] +
      std::complex<double> (0, 1) * amp[962] + std::complex<double> (0, 1) *
      amp[957] + std::complex<double> (0, 1) * amp[983] - std::complex<double>
      (0, 1) * amp[968];
  jamp[3] = +2. * (+std::complex<double> (0, 1) * amp[135] -
      std::complex<double> (0, 1) * amp[154] + std::complex<double> (0, 1) *
      amp[198] + std::complex<double> (0, 1) * amp[215] + std::complex<double>
      (0, 1) * amp[241] - std::complex<double> (0, 1) * amp[260] +
      std::complex<double> (0, 1) * amp[304] + std::complex<double> (0, 1) *
      amp[321]);
  jamp[4] = +2. * (-std::complex<double> (0, 1) * amp[135] +
      std::complex<double> (0, 1) * amp[154] - std::complex<double> (0, 1) *
      amp[198] - std::complex<double> (0, 1) * amp[215] - std::complex<double>
      (0, 1) * amp[241] + std::complex<double> (0, 1) * amp[260] -
      std::complex<double> (0, 1) * amp[304] - std::complex<double> (0, 1) *
      amp[321]);
  jamp[5] = -amp[132] - amp[134] - amp[136] - amp[137] - amp[138] - amp[142] -
      std::complex<double> (0, 1) * amp[145] - std::complex<double> (0, 1) *
      amp[146] - std::complex<double> (0, 1) * amp[147] - std::complex<double>
      (0, 1) * amp[148] - std::complex<double> (0, 1) * amp[149] - amp[151] +
      std::complex<double> (0, 1) * amp[152] - amp[153] - amp[155] - amp[156] -
      amp[157] + std::complex<double> (0, 1) * amp[158] + std::complex<double>
      (0, 1) * amp[159] + std::complex<double> (0, 1) * amp[160] - amp[161] +
      std::complex<double> (0, 1) * amp[162] - std::complex<double> (0, 1) *
      amp[180] - amp[184] - amp[185] - amp[186] - amp[188] - amp[190] +
      std::complex<double> (0, 1) * amp[205] + amp[208] + amp[207] + amp[211] +
      amp[210] + amp[214] + amp[213] + amp[218] + amp[217] + amp[366] +
      amp[370] - std::complex<double> (0, 1) * amp[372] + amp[373] + amp[374] +
      amp[375] - std::complex<double> (0, 1) * amp[379] - std::complex<double>
      (0, 1) * amp[380] - std::complex<double> (0, 1) * amp[382] -
      std::complex<double> (0, 1) * amp[383] - std::complex<double> (0, 1) *
      amp[384] + std::complex<double> (0, 1) * amp[385] + std::complex<double>
      (0, 1) * amp[386] + std::complex<double> (0, 1) * amp[388] +
      std::complex<double> (0, 1) * amp[393] + std::complex<double> (0, 1) *
      amp[394] + std::complex<double> (0, 1) * amp[396] + std::complex<double>
      (0, 1) * amp[397] + std::complex<double> (0, 1) * amp[399] +
      std::complex<double> (0, 1) * amp[400] + std::complex<double> (0, 1) *
      amp[401] + std::complex<double> (0, 1) * amp[405] + std::complex<double>
      (0, 1) * amp[406] + std::complex<double> (0, 1) * amp[408] + amp[424] +
      std::complex<double> (0, 1) * amp[425] + amp[426] + amp[428] + amp[429] +
      amp[430] + std::complex<double> (0, 1) * amp[431] + std::complex<double>
      (0, 1) * amp[432] + std::complex<double> (0, 1) * amp[433] + amp[434] +
      std::complex<double> (0, 1) * amp[435] + amp[438] - std::complex<double>
      (0, 1) * amp[439] - std::complex<double> (0, 1) * amp[626] -
      std::complex<double> (0, 1) * amp[627] - std::complex<double> (0, 1) *
      amp[629] - std::complex<double> (0, 1) * amp[630] - std::complex<double>
      (0, 1) * amp[634] - std::complex<double> (0, 1) * amp[635] -
      std::complex<double> (0, 1) * amp[637] - std::complex<double> (0, 1) *
      amp[638] - std::complex<double> (0, 1) * amp[640] - std::complex<double>
      (0, 1) * amp[641] - std::complex<double> (0, 1) * amp[642] -
      std::complex<double> (0, 1) * amp[643] - std::complex<double> (0, 1) *
      amp[644] - std::complex<double> (0, 1) * amp[645] - std::complex<double>
      (0, 1) * amp[646] - std::complex<double> (0, 1) * amp[647] -
      std::complex<double> (0, 1) * amp[649] - std::complex<double> (0, 1) *
      amp[650] - amp[665] - amp[667] - amp[669] - amp[670] - amp[671] -
      amp[675] + std::complex<double> (0, 1) * amp[677] - std::complex<double>
      (0, 1) * amp[679] - std::complex<double> (0, 1) * amp[680] -
      std::complex<double> (0, 1) * amp[682] - std::complex<double> (0, 1) *
      amp[687] - std::complex<double> (0, 1) * amp[688] - std::complex<double>
      (0, 1) * amp[690] - std::complex<double> (0, 1) * amp[691] -
      std::complex<double> (0, 1) * amp[695] - std::complex<double> (0, 1) *
      amp[696] - std::complex<double> (0, 1) * amp[697] - std::complex<double>
      (0, 1) * amp[699] - std::complex<double> (0, 1) * amp[700] -
      std::complex<double> (0, 1) * amp[702] - std::complex<double> (0, 1) *
      amp[784] - std::complex<double> (0, 1) * amp[783] - std::complex<double>
      (0, 1) * amp[790] - std::complex<double> (0, 1) * amp[789] -
      std::complex<double> (0, 1) * amp[793] - std::complex<double> (0, 1) *
      amp[792] - std::complex<double> (0, 1) * amp[796] - std::complex<double>
      (0, 1) * amp[795] - amp[812] - amp[811] - amp[815] - amp[814] - amp[818]
      - amp[817] - amp[821] - amp[820] - std::complex<double> (0, 1) * amp[825]
      - std::complex<double> (0, 1) * amp[824] - std::complex<double> (0, 1) *
      amp[871] - std::complex<double> (0, 1) * amp[872] - std::complex<double>
      (0, 1) * amp[875] - std::complex<double> (0, 1) * amp[874] + amp[878] +
      std::complex<double> (0, 1) * amp[879] - std::complex<double> (0, 1) *
      amp[939] - std::complex<double> (0, 1) * amp[940] - std::complex<double>
      (0, 1) * amp[945] - std::complex<double> (0, 1) * amp[946] -
      std::complex<double> (0, 1) * amp[948] - std::complex<double> (0, 1) *
      amp[949] - std::complex<double> (0, 1) * amp[951] - std::complex<double>
      (0, 1) * amp[952] - std::complex<double> (0, 1) * amp[980] -
      std::complex<double> (0, 1) * amp[965] + std::complex<double> (0, 1) *
      amp[971] - std::complex<double> (0, 1) * amp[966] + std::complex<double>
      (0, 1) * amp[982] + std::complex<double> (0, 1) * amp[972] -
      std::complex<double> (0, 1) * amp[957] + std::complex<double> (0, 1) *
      amp[978] - std::complex<double> (0, 1) * amp[958] + std::complex<double>
      (0, 1) * amp[979];
  jamp[6] = +2. * (-amp[150] - amp[163] + std::complex<double> (0, 1) *
      amp[368] + std::complex<double> (0, 1) * amp[369] - amp[716]);
  jamp[7] = +std::complex<double> (0, 1) * amp[114] - amp[115] - amp[116] -
      amp[117] - amp[121] - amp[123] + std::complex<double> (0, 1) * amp[126] +
      std::complex<double> (0, 1) * amp[127] + std::complex<double> (0, 1) *
      amp[128] + std::complex<double> (0, 1) * amp[129] + std::complex<double>
      (0, 1) * amp[130] + amp[151] - std::complex<double> (0, 1) * amp[152] +
      amp[153] + amp[155] + amp[156] + amp[157] - std::complex<double> (0, 1) *
      amp[158] - std::complex<double> (0, 1) * amp[159] - std::complex<double>
      (0, 1) * amp[160] + amp[161] - std::complex<double> (0, 1) * amp[162] -
      amp[193] - amp[194] - amp[199] - amp[200] - amp[201] - amp[203] -
      std::complex<double> (0, 1) * amp[205] - amp[207] - amp[206] - amp[210] -
      amp[209] - amp[213] - amp[212] - amp[217] - amp[216] +
      std::complex<double> (0, 1) * amp[346] + amp[347] + amp[348] +
      std::complex<double> (0, 1) * amp[349] + std::complex<double> (0, 1) *
      amp[350] + std::complex<double> (0, 1) * amp[351] + amp[353] + amp[354] +
      amp[355] + std::complex<double> (0, 1) * amp[356] + amp[357] - amp[366] -
      amp[370] + std::complex<double> (0, 1) * amp[372] - amp[373] - amp[374] -
      amp[375] + std::complex<double> (0, 1) * amp[379] + std::complex<double>
      (0, 1) * amp[380] + std::complex<double> (0, 1) * amp[382] +
      std::complex<double> (0, 1) * amp[383] + std::complex<double> (0, 1) *
      amp[384] - std::complex<double> (0, 1) * amp[386] - std::complex<double>
      (0, 1) * amp[387] - std::complex<double> (0, 1) * amp[389] -
      std::complex<double> (0, 1) * amp[394] - std::complex<double> (0, 1) *
      amp[395] - std::complex<double> (0, 1) * amp[397] - std::complex<double>
      (0, 1) * amp[398] - std::complex<double> (0, 1) * amp[402] -
      std::complex<double> (0, 1) * amp[403] - std::complex<double> (0, 1) *
      amp[404] - std::complex<double> (0, 1) * amp[406] - std::complex<double>
      (0, 1) * amp[407] - std::complex<double> (0, 1) * amp[409] +
      std::complex<double> (0, 1) * amp[436] - amp[438] + amp[614] + amp[615] +
      amp[620] + amp[621] + amp[622] + amp[624] + std::complex<double> (0, 1) *
      amp[679] + std::complex<double> (0, 1) * amp[680] + std::complex<double>
      (0, 1) * amp[682] + std::complex<double> (0, 1) * amp[687] +
      std::complex<double> (0, 1) * amp[688] + std::complex<double> (0, 1) *
      amp[690] + std::complex<double> (0, 1) * amp[691] + std::complex<double>
      (0, 1) * amp[695] + std::complex<double> (0, 1) * amp[696] +
      std::complex<double> (0, 1) * amp[697] + std::complex<double> (0, 1) *
      amp[699] + std::complex<double> (0, 1) * amp[700] + std::complex<double>
      (0, 1) * amp[702] - std::complex<double> (0, 1) * amp[734] -
      std::complex<double> (0, 1) * amp[735] - std::complex<double> (0, 1) *
      amp[736] - std::complex<double> (0, 1) * amp[737] - std::complex<double>
      (0, 1) * amp[739] - std::complex<double> (0, 1) * amp[740] -
      std::complex<double> (0, 1) * amp[742] - std::complex<double> (0, 1) *
      amp[743] - std::complex<double> (0, 1) * amp[744] - std::complex<double>
      (0, 1) * amp[745] - std::complex<double> (0, 1) * amp[746] -
      std::complex<double> (0, 1) * amp[747] - std::complex<double> (0, 1) *
      amp[748] - std::complex<double> (0, 1) * amp[749] - std::complex<double>
      (0, 1) * amp[751] - std::complex<double> (0, 1) * amp[752] -
      std::complex<double> (0, 1) * amp[753] - std::complex<double> (0, 1) *
      amp[754] + std::complex<double> (0, 1) * amp[872] - std::complex<double>
      (0, 1) * amp[873] + std::complex<double> (0, 1) * amp[875] +
      std::complex<double> (0, 1) * amp[876] + std::complex<double> (0, 1) *
      amp[877] - amp[878] - amp[904] + amp[902] - amp[907] + amp[905] -
      amp[910] + amp[908] - amp[913] + amp[911] - std::complex<double> (0, 1) *
      amp[915] + std::complex<double> (0, 1) * amp[917] - std::complex<double>
      (0, 1) * amp[921] + std::complex<double> (0, 1) * amp[923] -
      std::complex<double> (0, 1) * amp[924] + std::complex<double> (0, 1) *
      amp[926] - std::complex<double> (0, 1) * amp[927] + std::complex<double>
      (0, 1) * amp[929] - std::complex<double> (0, 1) * amp[932] +
      std::complex<double> (0, 1) * amp[930] + std::complex<double> (0, 1) *
      amp[933] + std::complex<double> (0, 1) * amp[939] + std::complex<double>
      (0, 1) * amp[938] + std::complex<double> (0, 1) * amp[945] +
      std::complex<double> (0, 1) * amp[944] + std::complex<double> (0, 1) *
      amp[948] + std::complex<double> (0, 1) * amp[947] + std::complex<double>
      (0, 1) * amp[951] + std::complex<double> (0, 1) * amp[950] +
      std::complex<double> (0, 1) * amp[959] - std::complex<double> (0, 1) *
      amp[985] + std::complex<double> (0, 1) * amp[980] - std::complex<double>
      (0, 1) * amp[971] - std::complex<double> (0, 1) * amp[961] +
      std::complex<double> (0, 1) * amp[956] - std::complex<double> (0, 1) *
      amp[982] + std::complex<double> (0, 1) * amp[957] + std::complex<double>
      (0, 1) * amp[973] - std::complex<double> (0, 1) * amp[984];
  jamp[8] = -amp[132] + std::complex<double> (0, 1) * amp[133] - amp[134] -
      amp[136] - amp[137] - amp[138] + std::complex<double> (0, 1) * amp[139] +
      std::complex<double> (0, 1) * amp[140] + std::complex<double> (0, 1) *
      amp[141] - amp[142] + std::complex<double> (0, 1) * amp[143] - amp[151] -
      amp[153] - amp[155] - amp[156] - amp[157] - amp[161] -
      std::complex<double> (0, 1) * amp[164] - std::complex<double> (0, 1) *
      amp[165] - std::complex<double> (0, 1) * amp[166] - std::complex<double>
      (0, 1) * amp[167] - std::complex<double> (0, 1) * amp[168] -
      std::complex<double> (0, 1) * amp[169] - amp[173] - amp[174] - amp[175] -
      amp[177] - amp[179] + std::complex<double> (0, 1) * amp[204] + amp[208] +
      amp[207] + amp[211] + amp[210] + amp[214] + amp[213] + amp[218] +
      amp[217] + amp[480] + amp[484] - std::complex<double> (0, 1) * amp[486] +
      amp[487] + amp[488] + amp[489] - std::complex<double> (0, 1) * amp[493] -
      std::complex<double> (0, 1) * amp[494] - std::complex<double> (0, 1) *
      amp[496] - std::complex<double> (0, 1) * amp[497] - std::complex<double>
      (0, 1) * amp[498] + std::complex<double> (0, 1) * amp[499] +
      std::complex<double> (0, 1) * amp[500] + std::complex<double> (0, 1) *
      amp[502] + std::complex<double> (0, 1) * amp[507] + std::complex<double>
      (0, 1) * amp[508] + std::complex<double> (0, 1) * amp[510] +
      std::complex<double> (0, 1) * amp[511] + std::complex<double> (0, 1) *
      amp[513] + std::complex<double> (0, 1) * amp[514] + std::complex<double>
      (0, 1) * amp[515] + std::complex<double> (0, 1) * amp[519] +
      std::complex<double> (0, 1) * amp[520] + std::complex<double> (0, 1) *
      amp[522] + amp[538] + std::complex<double> (0, 1) * amp[539] + amp[540] +
      amp[542] + amp[543] + amp[544] + std::complex<double> (0, 1) * amp[545] +
      std::complex<double> (0, 1) * amp[546] + std::complex<double> (0, 1) *
      amp[547] + amp[548] + std::complex<double> (0, 1) * amp[549] + amp[552] -
      std::complex<double> (0, 1) * amp[553] - std::complex<double> (0, 1) *
      amp[627] - std::complex<double> (0, 1) * amp[628] - std::complex<double>
      (0, 1) * amp[630] - std::complex<double> (0, 1) * amp[635] -
      std::complex<double> (0, 1) * amp[636] - std::complex<double> (0, 1) *
      amp[638] - std::complex<double> (0, 1) * amp[639] - std::complex<double>
      (0, 1) * amp[643] - std::complex<double> (0, 1) * amp[644] -
      std::complex<double> (0, 1) * amp[645] - std::complex<double> (0, 1) *
      amp[647] - std::complex<double> (0, 1) * amp[648] - std::complex<double>
      (0, 1) * amp[650] - std::complex<double> (0, 1) * amp[678] -
      std::complex<double> (0, 1) * amp[679] - std::complex<double> (0, 1) *
      amp[681] - std::complex<double> (0, 1) * amp[682] - std::complex<double>
      (0, 1) * amp[686] - std::complex<double> (0, 1) * amp[687] -
      std::complex<double> (0, 1) * amp[689] - std::complex<double> (0, 1) *
      amp[690] - std::complex<double> (0, 1) * amp[692] - std::complex<double>
      (0, 1) * amp[693] - std::complex<double> (0, 1) * amp[694] -
      std::complex<double> (0, 1) * amp[695] - std::complex<double> (0, 1) *
      amp[696] - std::complex<double> (0, 1) * amp[697] - std::complex<double>
      (0, 1) * amp[698] - std::complex<double> (0, 1) * amp[699] -
      std::complex<double> (0, 1) * amp[701] - std::complex<double> (0, 1) *
      amp[702] - amp[717] - amp[719] - amp[721] - amp[722] - amp[723] -
      amp[727] + std::complex<double> (0, 1) * amp[729] - std::complex<double>
      (0, 1) * amp[828] - std::complex<double> (0, 1) * amp[827] -
      std::complex<double> (0, 1) * amp[834] - std::complex<double> (0, 1) *
      amp[833] - std::complex<double> (0, 1) * amp[837] - std::complex<double>
      (0, 1) * amp[836] - std::complex<double> (0, 1) * amp[840] -
      std::complex<double> (0, 1) * amp[839] - amp[856] - amp[855] - amp[859] -
      amp[858] - amp[862] - amp[861] - amp[865] - amp[864] -
      std::complex<double> (0, 1) * amp[869] - std::complex<double> (0, 1) *
      amp[868] - std::complex<double> (0, 1) * amp[871] - std::complex<double>
      (0, 1) * amp[872] - std::complex<double> (0, 1) * amp[875] -
      std::complex<double> (0, 1) * amp[874] + amp[881] + std::complex<double>
      (0, 1) * amp[882] - std::complex<double> (0, 1) * amp[939] -
      std::complex<double> (0, 1) * amp[940] - std::complex<double> (0, 1) *
      amp[945] - std::complex<double> (0, 1) * amp[946] - std::complex<double>
      (0, 1) * amp[948] - std::complex<double> (0, 1) * amp[949] -
      std::complex<double> (0, 1) * amp[951] - std::complex<double> (0, 1) *
      amp[952] + std::complex<double> (0, 1) * amp[974] + std::complex<double>
      (0, 1) * amp[975] + std::complex<double> (0, 1) * amp[981] -
      std::complex<double> (0, 1) * amp[977] + std::complex<double> (0, 1) *
      amp[982] - std::complex<double> (0, 1) * amp[962] - std::complex<double>
      (0, 1) * amp[957] - std::complex<double> (0, 1) * amp[963] -
      std::complex<double> (0, 1) * amp[958] + std::complex<double> (0, 1) *
      amp[979];
  jamp[9] = +2. * (-amp[131] - amp[144] + std::complex<double> (0, 1) *
      amp[482] + std::complex<double> (0, 1) * amp[483] - amp[664]);
  jamp[10] = -std::complex<double> (0, 1) * amp[114] - amp[118] - amp[119] -
      amp[120] - amp[122] - amp[124] - std::complex<double> (0, 1) * amp[126] -
      std::complex<double> (0, 1) * amp[127] - std::complex<double> (0, 1) *
      amp[128] - std::complex<double> (0, 1) * amp[129] - std::complex<double>
      (0, 1) * amp[130] + amp[132] - std::complex<double> (0, 1) * amp[133] +
      amp[134] + amp[136] + amp[137] + amp[138] - std::complex<double> (0, 1) *
      amp[139] - std::complex<double> (0, 1) * amp[140] - std::complex<double>
      (0, 1) * amp[141] + amp[142] - std::complex<double> (0, 1) * amp[143] +
      amp[193] + amp[194] + amp[199] + amp[200] + amp[201] + amp[203] -
      std::complex<double> (0, 1) * amp[204] - amp[208] + amp[206] - amp[211] +
      amp[209] - amp[214] + amp[212] - amp[218] + amp[216] +
      std::complex<double> (0, 1) * amp[460] + amp[461] + amp[462] +
      std::complex<double> (0, 1) * amp[463] + std::complex<double> (0, 1) *
      amp[464] + std::complex<double> (0, 1) * amp[465] + amp[467] + amp[468] +
      amp[469] + std::complex<double> (0, 1) * amp[470] + amp[471] - amp[480] -
      amp[484] + std::complex<double> (0, 1) * amp[486] - amp[487] - amp[488] -
      amp[489] + std::complex<double> (0, 1) * amp[493] + std::complex<double>
      (0, 1) * amp[494] + std::complex<double> (0, 1) * amp[496] +
      std::complex<double> (0, 1) * amp[497] + std::complex<double> (0, 1) *
      amp[498] - std::complex<double> (0, 1) * amp[500] - std::complex<double>
      (0, 1) * amp[501] - std::complex<double> (0, 1) * amp[503] -
      std::complex<double> (0, 1) * amp[508] - std::complex<double> (0, 1) *
      amp[509] - std::complex<double> (0, 1) * amp[511] - std::complex<double>
      (0, 1) * amp[512] - std::complex<double> (0, 1) * amp[516] -
      std::complex<double> (0, 1) * amp[517] - std::complex<double> (0, 1) *
      amp[518] - std::complex<double> (0, 1) * amp[520] - std::complex<double>
      (0, 1) * amp[521] - std::complex<double> (0, 1) * amp[523] +
      std::complex<double> (0, 1) * amp[550] - amp[552] - amp[614] - amp[615] -
      amp[620] - amp[621] - amp[622] - amp[624] + std::complex<double> (0, 1) *
      amp[627] + std::complex<double> (0, 1) * amp[628] + std::complex<double>
      (0, 1) * amp[630] + std::complex<double> (0, 1) * amp[635] +
      std::complex<double> (0, 1) * amp[636] + std::complex<double> (0, 1) *
      amp[638] + std::complex<double> (0, 1) * amp[639] + std::complex<double>
      (0, 1) * amp[643] + std::complex<double> (0, 1) * amp[644] +
      std::complex<double> (0, 1) * amp[645] + std::complex<double> (0, 1) *
      amp[647] + std::complex<double> (0, 1) * amp[648] + std::complex<double>
      (0, 1) * amp[650] + std::complex<double> (0, 1) * amp[734] +
      std::complex<double> (0, 1) * amp[735] + std::complex<double> (0, 1) *
      amp[736] + std::complex<double> (0, 1) * amp[737] + std::complex<double>
      (0, 1) * amp[739] + std::complex<double> (0, 1) * amp[740] +
      std::complex<double> (0, 1) * amp[742] + std::complex<double> (0, 1) *
      amp[743] + std::complex<double> (0, 1) * amp[744] + std::complex<double>
      (0, 1) * amp[745] + std::complex<double> (0, 1) * amp[746] +
      std::complex<double> (0, 1) * amp[747] + std::complex<double> (0, 1) *
      amp[748] + std::complex<double> (0, 1) * amp[749] + std::complex<double>
      (0, 1) * amp[751] + std::complex<double> (0, 1) * amp[752] +
      std::complex<double> (0, 1) * amp[753] + std::complex<double> (0, 1) *
      amp[754] + std::complex<double> (0, 1) * amp[871] + std::complex<double>
      (0, 1) * amp[873] + std::complex<double> (0, 1) * amp[874] -
      std::complex<double> (0, 1) * amp[876] + std::complex<double> (0, 1) *
      amp[880] - amp[881] - amp[903] - amp[902] - amp[906] - amp[905] -
      amp[909] - amp[908] - amp[912] - amp[911] + std::complex<double> (0, 1) *
      amp[916] + std::complex<double> (0, 1) * amp[915] + std::complex<double>
      (0, 1) * amp[922] + std::complex<double> (0, 1) * amp[921] +
      std::complex<double> (0, 1) * amp[925] + std::complex<double> (0, 1) *
      amp[924] + std::complex<double> (0, 1) * amp[928] + std::complex<double>
      (0, 1) * amp[927] - std::complex<double> (0, 1) * amp[931] -
      std::complex<double> (0, 1) * amp[930] - std::complex<double> (0, 1) *
      amp[933] - std::complex<double> (0, 1) * amp[938] + std::complex<double>
      (0, 1) * amp[940] - std::complex<double> (0, 1) * amp[944] +
      std::complex<double> (0, 1) * amp[946] - std::complex<double> (0, 1) *
      amp[947] + std::complex<double> (0, 1) * amp[949] - std::complex<double>
      (0, 1) * amp[950] + std::complex<double> (0, 1) * amp[952] -
      std::complex<double> (0, 1) * amp[974] - std::complex<double> (0, 1) *
      amp[959] + std::complex<double> (0, 1) * amp[985] - std::complex<double>
      (0, 1) * amp[960] + std::complex<double> (0, 1) * amp[976] -
      std::complex<double> (0, 1) * amp[956] + std::complex<double> (0, 1) *
      amp[977] + std::complex<double> (0, 1) * amp[958] - std::complex<double>
      (0, 1) * amp[979] + std::complex<double> (0, 1) * amp[984];
  jamp[11] = +2. * (-std::complex<double> (0, 1) * amp[99] -
      std::complex<double> (0, 1) * amp[100] - amp[231] - amp[297] - amp[576]);
  jamp[12] = -std::complex<double> (0, 1) * amp[0] - std::complex<double> (0,
      1) * amp[1] - std::complex<double> (0, 1) * amp[3] - std::complex<double>
      (0, 1) * amp[8] - std::complex<double> (0, 1) * amp[9] -
      std::complex<double> (0, 1) * amp[11] - std::complex<double> (0, 1) *
      amp[12] - std::complex<double> (0, 1) * amp[14] - std::complex<double>
      (0, 1) * amp[15] - std::complex<double> (0, 1) * amp[16] -
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[21] - std::complex<double> (0, 1) * amp[23] - std::complex<double>
      (0, 1) * amp[25] - std::complex<double> (0, 1) * amp[27] -
      std::complex<double> (0, 1) * amp[28] - std::complex<double> (0, 1) *
      amp[29] - std::complex<double> (0, 1) * amp[30] - std::complex<double>
      (0, 1) * amp[83] - amp[84] - amp[85] - std::complex<double> (0, 1) *
      amp[86] - std::complex<double> (0, 1) * amp[87] - std::complex<double>
      (0, 1) * amp[88] - amp[90] - amp[91] - amp[92] - std::complex<double> (0,
      1) * amp[93] - amp[94] - amp[97] - amp[101] + std::complex<double> (0, 1)
      * amp[103] - amp[104] - amp[105] - amp[106] - std::complex<double> (0, 1)
      * amp[111] - amp[113] + amp[328] + amp[329] + amp[334] + amp[335] +
      amp[336] + amp[338] + std::complex<double> (0, 1) * amp[340] +
      std::complex<double> (0, 1) * amp[341] + std::complex<double> (0, 1) *
      amp[342] + std::complex<double> (0, 1) * amp[343] + std::complex<double>
      (0, 1) * amp[344] - std::complex<double> (0, 1) * amp[385] +
      std::complex<double> (0, 1) * amp[387] - std::complex<double> (0, 1) *
      amp[388] + std::complex<double> (0, 1) * amp[389] - std::complex<double>
      (0, 1) * amp[393] + std::complex<double> (0, 1) * amp[395] -
      std::complex<double> (0, 1) * amp[396] + std::complex<double> (0, 1) *
      amp[398] - std::complex<double> (0, 1) * amp[399] - std::complex<double>
      (0, 1) * amp[400] - std::complex<double> (0, 1) * amp[401] +
      std::complex<double> (0, 1) * amp[402] + std::complex<double> (0, 1) *
      amp[403] + std::complex<double> (0, 1) * amp[404] - std::complex<double>
      (0, 1) * amp[405] + std::complex<double> (0, 1) * amp[407] -
      std::complex<double> (0, 1) * amp[408] + std::complex<double> (0, 1) *
      amp[409] - amp[424] - amp[426] - amp[428] - amp[429] - amp[430] -
      amp[434] - std::complex<double> (0, 1) * amp[436] + std::complex<double>
      (0, 1) * amp[439] + std::complex<double> (0, 1) * amp[565] - amp[566] -
      amp[567] - amp[568] - amp[572] - amp[574] + std::complex<double> (0, 1) *
      amp[577] + amp[578] + amp[579] + std::complex<double> (0, 1) * amp[580] +
      std::complex<double> (0, 1) * amp[581] + std::complex<double> (0, 1) *
      amp[582] + amp[584] + amp[585] + amp[586] + std::complex<double> (0, 1) *
      amp[587] + amp[588] + std::complex<double> (0, 1) * amp[589] -
      std::complex<double> (0, 1) * amp[733] + std::complex<double> (0, 1) *
      amp[735] + std::complex<double> (0, 1) * amp[737] - std::complex<double>
      (0, 1) * amp[738] + std::complex<double> (0, 1) * amp[740] -
      std::complex<double> (0, 1) * amp[741] + std::complex<double> (0, 1) *
      amp[743] + std::complex<double> (0, 1) * amp[747] + std::complex<double>
      (0, 1) * amp[748] + std::complex<double> (0, 1) * amp[749] -
      std::complex<double> (0, 1) * amp[750] + std::complex<double> (0, 1) *
      amp[752] + std::complex<double> (0, 1) * amp[754] - std::complex<double>
      (0, 1) * amp[778] - amp[779] + std::complex<double> (0, 1) * amp[784] +
      std::complex<double> (0, 1) * amp[785] + std::complex<double> (0, 1) *
      amp[790] + std::complex<double> (0, 1) * amp[791] + std::complex<double>
      (0, 1) * amp[793] + std::complex<double> (0, 1) * amp[794] +
      std::complex<double> (0, 1) * amp[796] + std::complex<double> (0, 1) *
      amp[797] + amp[813] + amp[812] + amp[816] + amp[815] + amp[819] +
      amp[818] + amp[822] + amp[821] + std::complex<double> (0, 1) * amp[826] +
      std::complex<double> (0, 1) * amp[825] - amp[891] + amp[889] - amp[894] +
      amp[892] - amp[897] + amp[895] - amp[900] + amp[898] +
      std::complex<double> (0, 1) * amp[915] - std::complex<double> (0, 1) *
      amp[917] + std::complex<double> (0, 1) * amp[921] - std::complex<double>
      (0, 1) * amp[923] + std::complex<double> (0, 1) * amp[924] -
      std::complex<double> (0, 1) * amp[926] + std::complex<double> (0, 1) *
      amp[927] - std::complex<double> (0, 1) * amp[929] + std::complex<double>
      (0, 1) * amp[932] - std::complex<double> (0, 1) * amp[930] -
      std::complex<double> (0, 1) * amp[933] - std::complex<double> (0, 1) *
      amp[959] + std::complex<double> (0, 1) * amp[961] + std::complex<double>
      (0, 1) * amp[966] + std::complex<double> (0, 1) * amp[967] -
      std::complex<double> (0, 1) * amp[972] + std::complex<double> (0, 1) *
      amp[983] - std::complex<double> (0, 1) * amp[968] - std::complex<double>
      (0, 1) * amp[973] + std::complex<double> (0, 1) * amp[984] -
      std::complex<double> (0, 1) * amp[969];
  jamp[13] = +std::complex<double> (0, 1) * amp[0] - std::complex<double> (0,
      1) * amp[2] - std::complex<double> (0, 1) * amp[4] + std::complex<double>
      (0, 1) * amp[8] - std::complex<double> (0, 1) * amp[10] +
      std::complex<double> (0, 1) * amp[11] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[17] - std::complex<double>
      (0, 1) * amp[18] - std::complex<double> (0, 1) * amp[19] +
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[22] - std::complex<double> (0, 1) * amp[24] + std::complex<double>
      (0, 1) * amp[25] + std::complex<double> (0, 1) * amp[27] +
      std::complex<double> (0, 1) * amp[28] + std::complex<double> (0, 1) *
      amp[29] + std::complex<double> (0, 1) * amp[30] - std::complex<double>
      (0, 1) * amp[70] - amp[71] - amp[72] - std::complex<double> (0, 1) *
      amp[73] - std::complex<double> (0, 1) * amp[74] - std::complex<double>
      (0, 1) * amp[75] - amp[77] - amp[78] - amp[79] - std::complex<double> (0,
      1) * amp[80] - amp[81] + amp[97] + amp[101] - std::complex<double> (0, 1)
      * amp[103] + amp[104] + amp[105] + amp[106] - std::complex<double> (0, 1)
      * amp[110] + amp[113] + amp[442] + amp[443] + amp[448] + amp[449] +
      amp[450] + amp[452] + std::complex<double> (0, 1) * amp[454] +
      std::complex<double> (0, 1) * amp[455] + std::complex<double> (0, 1) *
      amp[456] + std::complex<double> (0, 1) * amp[457] + std::complex<double>
      (0, 1) * amp[458] - std::complex<double> (0, 1) * amp[499] +
      std::complex<double> (0, 1) * amp[501] - std::complex<double> (0, 1) *
      amp[502] + std::complex<double> (0, 1) * amp[503] - std::complex<double>
      (0, 1) * amp[507] + std::complex<double> (0, 1) * amp[509] -
      std::complex<double> (0, 1) * amp[510] + std::complex<double> (0, 1) *
      amp[512] - std::complex<double> (0, 1) * amp[513] - std::complex<double>
      (0, 1) * amp[514] - std::complex<double> (0, 1) * amp[515] +
      std::complex<double> (0, 1) * amp[516] + std::complex<double> (0, 1) *
      amp[517] + std::complex<double> (0, 1) * amp[518] - std::complex<double>
      (0, 1) * amp[519] + std::complex<double> (0, 1) * amp[521] -
      std::complex<double> (0, 1) * amp[522] + std::complex<double> (0, 1) *
      amp[523] - amp[538] - amp[540] - amp[542] - amp[543] - amp[544] -
      amp[548] - std::complex<double> (0, 1) * amp[550] + std::complex<double>
      (0, 1) * amp[553] + std::complex<double> (0, 1) * amp[554] - amp[555] -
      amp[556] - amp[557] - amp[561] - amp[563] - std::complex<double> (0, 1) *
      amp[577] - amp[578] - amp[579] - std::complex<double> (0, 1) * amp[580] -
      std::complex<double> (0, 1) * amp[581] - std::complex<double> (0, 1) *
      amp[582] - amp[584] - amp[585] - amp[586] - std::complex<double> (0, 1) *
      amp[587] - amp[588] - std::complex<double> (0, 1) * amp[589] +
      std::complex<double> (0, 1) * amp[733] - std::complex<double> (0, 1) *
      amp[735] - std::complex<double> (0, 1) * amp[737] + std::complex<double>
      (0, 1) * amp[738] - std::complex<double> (0, 1) * amp[740] +
      std::complex<double> (0, 1) * amp[741] - std::complex<double> (0, 1) *
      amp[743] - std::complex<double> (0, 1) * amp[747] - std::complex<double>
      (0, 1) * amp[748] - std::complex<double> (0, 1) * amp[749] +
      std::complex<double> (0, 1) * amp[750] - std::complex<double> (0, 1) *
      amp[752] - std::complex<double> (0, 1) * amp[754] - std::complex<double>
      (0, 1) * amp[777] + amp[779] + std::complex<double> (0, 1) * amp[828] +
      std::complex<double> (0, 1) * amp[829] + std::complex<double> (0, 1) *
      amp[834] + std::complex<double> (0, 1) * amp[835] + std::complex<double>
      (0, 1) * amp[837] + std::complex<double> (0, 1) * amp[838] +
      std::complex<double> (0, 1) * amp[840] + std::complex<double> (0, 1) *
      amp[841] + amp[857] + amp[856] + amp[860] + amp[859] + amp[863] +
      amp[862] + amp[866] + amp[865] + std::complex<double> (0, 1) * amp[870] +
      std::complex<double> (0, 1) * amp[869] - amp[890] - amp[889] - amp[893] -
      amp[892] - amp[896] - amp[895] - amp[899] - amp[898] -
      std::complex<double> (0, 1) * amp[916] - std::complex<double> (0, 1) *
      amp[915] - std::complex<double> (0, 1) * amp[922] - std::complex<double>
      (0, 1) * amp[921] - std::complex<double> (0, 1) * amp[925] -
      std::complex<double> (0, 1) * amp[924] - std::complex<double> (0, 1) *
      amp[928] - std::complex<double> (0, 1) * amp[927] + std::complex<double>
      (0, 1) * amp[931] + std::complex<double> (0, 1) * amp[930] +
      std::complex<double> (0, 1) * amp[933] + std::complex<double> (0, 1) *
      amp[959] + std::complex<double> (0, 1) * amp[964] - std::complex<double>
      (0, 1) * amp[970] - std::complex<double> (0, 1) * amp[975] +
      std::complex<double> (0, 1) * amp[960] - std::complex<double> (0, 1) *
      amp[976] - std::complex<double> (0, 1) * amp[983] + std::complex<double>
      (0, 1) * amp[963] + std::complex<double> (0, 1) * amp[968] -
      std::complex<double> (0, 1) * amp[984];
  jamp[14] = +2. * (+std::complex<double> (0, 1) * amp[333] +
      std::complex<double> (0, 1) * amp[352] - std::complex<double> (0, 1) *
      amp[447] - std::complex<double> (0, 1) * amp[466] + std::complex<double>
      (0, 1) * amp[583] + std::complex<double> (0, 1) * amp[619] +
      std::complex<double> (0, 1) * amp[901] + std::complex<double> (0, 1) *
      amp[914]);
  jamp[15] = +2. * (-std::complex<double> (0, 1) * amp[333] -
      std::complex<double> (0, 1) * amp[352] + std::complex<double> (0, 1) *
      amp[447] + std::complex<double> (0, 1) * amp[466] - std::complex<double>
      (0, 1) * amp[583] - std::complex<double> (0, 1) * amp[619] -
      std::complex<double> (0, 1) * amp[901] - std::complex<double> (0, 1) *
      amp[914]);
  jamp[16] = -amp[328] - amp[329] - amp[334] - amp[335] - amp[336] - amp[338] -
      std::complex<double> (0, 1) * amp[340] - std::complex<double> (0, 1) *
      amp[341] - std::complex<double> (0, 1) * amp[342] - std::complex<double>
      (0, 1) * amp[343] - std::complex<double> (0, 1) * amp[344] +
      std::complex<double> (0, 1) * amp[385] - std::complex<double> (0, 1) *
      amp[387] + std::complex<double> (0, 1) * amp[388] - std::complex<double>
      (0, 1) * amp[389] + std::complex<double> (0, 1) * amp[393] -
      std::complex<double> (0, 1) * amp[395] + std::complex<double> (0, 1) *
      amp[396] - std::complex<double> (0, 1) * amp[398] + std::complex<double>
      (0, 1) * amp[399] + std::complex<double> (0, 1) * amp[400] +
      std::complex<double> (0, 1) * amp[401] - std::complex<double> (0, 1) *
      amp[402] - std::complex<double> (0, 1) * amp[403] - std::complex<double>
      (0, 1) * amp[404] + std::complex<double> (0, 1) * amp[405] -
      std::complex<double> (0, 1) * amp[407] + std::complex<double> (0, 1) *
      amp[408] - std::complex<double> (0, 1) * amp[409] + amp[424] + amp[426] +
      amp[428] + amp[429] + amp[430] + amp[434] + std::complex<double> (0, 1) *
      amp[436] - std::complex<double> (0, 1) * amp[439] + std::complex<double>
      (0, 1) * amp[441] - amp[442] - amp[443] + std::complex<double> (0, 1) *
      amp[444] + std::complex<double> (0, 1) * amp[445] + std::complex<double>
      (0, 1) * amp[446] - amp[448] - amp[449] - amp[450] + std::complex<double>
      (0, 1) * amp[451] - amp[452] + amp[481] + amp[485] + std::complex<double>
      (0, 1) * amp[486] + amp[490] + amp[491] + amp[492] + std::complex<double>
      (0, 1) * amp[493] + std::complex<double> (0, 1) * amp[494] +
      std::complex<double> (0, 1) * amp[496] + std::complex<double> (0, 1) *
      amp[497] + std::complex<double> (0, 1) * amp[498] - std::complex<double>
      (0, 1) * amp[500] - std::complex<double> (0, 1) * amp[501] -
      std::complex<double> (0, 1) * amp[503] - std::complex<double> (0, 1) *
      amp[508] - std::complex<double> (0, 1) * amp[509] - std::complex<double>
      (0, 1) * amp[511] - std::complex<double> (0, 1) * amp[512] -
      std::complex<double> (0, 1) * amp[516] - std::complex<double> (0, 1) *
      amp[517] - std::complex<double> (0, 1) * amp[518] - std::complex<double>
      (0, 1) * amp[520] - std::complex<double> (0, 1) * amp[521] -
      std::complex<double> (0, 1) * amp[523] + std::complex<double> (0, 1) *
      amp[550] + amp[551] - std::complex<double> (0, 1) * amp[565] - amp[569] -
      amp[570] - amp[571] - amp[573] - amp[575] - std::complex<double> (0, 1) *
      amp[626] + std::complex<double> (0, 1) * amp[628] - std::complex<double>
      (0, 1) * amp[629] - std::complex<double> (0, 1) * amp[634] +
      std::complex<double> (0, 1) * amp[636] - std::complex<double> (0, 1) *
      amp[637] + std::complex<double> (0, 1) * amp[639] - std::complex<double>
      (0, 1) * amp[640] - std::complex<double> (0, 1) * amp[641] -
      std::complex<double> (0, 1) * amp[642] - std::complex<double> (0, 1) *
      amp[646] + std::complex<double> (0, 1) * amp[648] - std::complex<double>
      (0, 1) * amp[649] - amp[665] + std::complex<double> (0, 1) * amp[666] -
      amp[667] - amp[669] - amp[670] - amp[671] + std::complex<double> (0, 1) *
      amp[672] + std::complex<double> (0, 1) * amp[673] + std::complex<double>
      (0, 1) * amp[674] - amp[675] + std::complex<double> (0, 1) * amp[676] +
      std::complex<double> (0, 1) * amp[677] - std::complex<double> (0, 1) *
      amp[784] - std::complex<double> (0, 1) * amp[783] - std::complex<double>
      (0, 1) * amp[790] - std::complex<double> (0, 1) * amp[789] -
      std::complex<double> (0, 1) * amp[793] - std::complex<double> (0, 1) *
      amp[792] - std::complex<double> (0, 1) * amp[796] - std::complex<double>
      (0, 1) * amp[795] - amp[812] - amp[811] - amp[815] - amp[814] - amp[818]
      - amp[817] - amp[821] - amp[820] - std::complex<double> (0, 1) * amp[825]
      - std::complex<double> (0, 1) * amp[824] + std::complex<double> (0, 1) *
      amp[886] + amp[887] + amp[891] + amp[890] + amp[894] + amp[893] +
      amp[897] + amp[896] + amp[900] + amp[899] + std::complex<double> (0, 1) *
      amp[916] + std::complex<double> (0, 1) * amp[917] + std::complex<double>
      (0, 1) * amp[922] + std::complex<double> (0, 1) * amp[923] +
      std::complex<double> (0, 1) * amp[925] + std::complex<double> (0, 1) *
      amp[926] + std::complex<double> (0, 1) * amp[928] + std::complex<double>
      (0, 1) * amp[929] - std::complex<double> (0, 1) * amp[932] -
      std::complex<double> (0, 1) * amp[931] + std::complex<double> (0, 1) *
      amp[934] - std::complex<double> (0, 1) * amp[974] - std::complex<double>
      (0, 1) * amp[960] - std::complex<double> (0, 1) * amp[965] +
      std::complex<double> (0, 1) * amp[976] - std::complex<double> (0, 1) *
      amp[961] - std::complex<double> (0, 1) * amp[966] + std::complex<double>
      (0, 1) * amp[977] + std::complex<double> (0, 1) * amp[972] +
      std::complex<double> (0, 1) * amp[978] + std::complex<double> (0, 1) *
      amp[973];
  jamp[17] = +2. * (-amp[440] - amp[453] + std::complex<double> (0, 1) *
      amp[478] + std::complex<double> (0, 1) * amp[479] - amp[537]);
  jamp[18] = +std::complex<double> (0, 1) * amp[220] - amp[224] - amp[225] -
      amp[226] - amp[228] - amp[230] + std::complex<double> (0, 1) * amp[232] +
      std::complex<double> (0, 1) * amp[233] + std::complex<double> (0, 1) *
      amp[234] + std::complex<double> (0, 1) * amp[235] + std::complex<double>
      (0, 1) * amp[236] + amp[238] + std::complex<double> (0, 1) * amp[239] +
      amp[240] + amp[242] + amp[243] + amp[244] + std::complex<double> (0, 1) *
      amp[245] + std::complex<double> (0, 1) * amp[246] + std::complex<double>
      (0, 1) * amp[247] + amp[248] + std::complex<double> (0, 1) * amp[249] +
      amp[299] + amp[300] + amp[305] + amp[306] + amp[307] + amp[309] +
      std::complex<double> (0, 1) * amp[310] - amp[314] + amp[312] - amp[317] +
      amp[315] - amp[320] + amp[318] - amp[324] + amp[322] -
      std::complex<double> (0, 1) * amp[441] + amp[442] + amp[443] -
      std::complex<double> (0, 1) * amp[444] - std::complex<double> (0, 1) *
      amp[445] - std::complex<double> (0, 1) * amp[446] + amp[448] + amp[449] +
      amp[450] - std::complex<double> (0, 1) * amp[451] + amp[452] - amp[481] -
      amp[485] - std::complex<double> (0, 1) * amp[486] - amp[490] - amp[491] -
      amp[492] - std::complex<double> (0, 1) * amp[493] - std::complex<double>
      (0, 1) * amp[494] - std::complex<double> (0, 1) * amp[496] -
      std::complex<double> (0, 1) * amp[497] - std::complex<double> (0, 1) *
      amp[498] + std::complex<double> (0, 1) * amp[500] + std::complex<double>
      (0, 1) * amp[501] + std::complex<double> (0, 1) * amp[503] +
      std::complex<double> (0, 1) * amp[508] + std::complex<double> (0, 1) *
      amp[509] + std::complex<double> (0, 1) * amp[511] + std::complex<double>
      (0, 1) * amp[512] + std::complex<double> (0, 1) * amp[516] +
      std::complex<double> (0, 1) * amp[517] + std::complex<double> (0, 1) *
      amp[518] + std::complex<double> (0, 1) * amp[520] + std::complex<double>
      (0, 1) * amp[521] + std::complex<double> (0, 1) * amp[523] -
      std::complex<double> (0, 1) * amp[550] - amp[551] - amp[578] - amp[579] -
      amp[584] - amp[585] - amp[586] - amp[588] - std::complex<double> (0, 1) *
      amp[627] - std::complex<double> (0, 1) * amp[628] - std::complex<double>
      (0, 1) * amp[630] - std::complex<double> (0, 1) * amp[635] -
      std::complex<double> (0, 1) * amp[636] - std::complex<double> (0, 1) *
      amp[638] - std::complex<double> (0, 1) * amp[639] - std::complex<double>
      (0, 1) * amp[643] - std::complex<double> (0, 1) * amp[644] -
      std::complex<double> (0, 1) * amp[645] - std::complex<double> (0, 1) *
      amp[647] - std::complex<double> (0, 1) * amp[648] - std::complex<double>
      (0, 1) * amp[650] - std::complex<double> (0, 1) * amp[734] -
      std::complex<double> (0, 1) * amp[735] - std::complex<double> (0, 1) *
      amp[736] - std::complex<double> (0, 1) * amp[737] - std::complex<double>
      (0, 1) * amp[739] - std::complex<double> (0, 1) * amp[740] -
      std::complex<double> (0, 1) * amp[742] - std::complex<double> (0, 1) *
      amp[743] - std::complex<double> (0, 1) * amp[744] - std::complex<double>
      (0, 1) * amp[745] - std::complex<double> (0, 1) * amp[746] -
      std::complex<double> (0, 1) * amp[747] - std::complex<double> (0, 1) *
      amp[748] - std::complex<double> (0, 1) * amp[749] - std::complex<double>
      (0, 1) * amp[751] - std::complex<double> (0, 1) * amp[752] -
      std::complex<double> (0, 1) * amp[753] - std::complex<double> (0, 1) *
      amp[754] - std::complex<double> (0, 1) * amp[871] - std::complex<double>
      (0, 1) * amp[873] - std::complex<double> (0, 1) * amp[874] +
      std::complex<double> (0, 1) * amp[876] - std::complex<double> (0, 1) *
      amp[886] - amp[887] - amp[890] - amp[889] - amp[893] - amp[892] -
      amp[896] - amp[895] - amp[899] - amp[898] - std::complex<double> (0, 1) *
      amp[916] - std::complex<double> (0, 1) * amp[915] - std::complex<double>
      (0, 1) * amp[922] - std::complex<double> (0, 1) * amp[921] -
      std::complex<double> (0, 1) * amp[925] - std::complex<double> (0, 1) *
      amp[924] - std::complex<double> (0, 1) * amp[928] - std::complex<double>
      (0, 1) * amp[927] + std::complex<double> (0, 1) * amp[931] +
      std::complex<double> (0, 1) * amp[930] + std::complex<double> (0, 1) *
      amp[933] + std::complex<double> (0, 1) * amp[938] - std::complex<double>
      (0, 1) * amp[940] + std::complex<double> (0, 1) * amp[944] -
      std::complex<double> (0, 1) * amp[946] + std::complex<double> (0, 1) *
      amp[947] - std::complex<double> (0, 1) * amp[949] + std::complex<double>
      (0, 1) * amp[950] - std::complex<double> (0, 1) * amp[952] +
      std::complex<double> (0, 1) * amp[974] + std::complex<double> (0, 1) *
      amp[959] - std::complex<double> (0, 1) * amp[985] + std::complex<double>
      (0, 1) * amp[960] - std::complex<double> (0, 1) * amp[976] +
      std::complex<double> (0, 1) * amp[956] - std::complex<double> (0, 1) *
      amp[977] - std::complex<double> (0, 1) * amp[958] + std::complex<double>
      (0, 1) * amp[979] - std::complex<double> (0, 1) * amp[984];
  jamp[19] = +std::complex<double> (0, 1) * amp[327] - amp[328] - amp[329] +
      std::complex<double> (0, 1) * amp[330] + std::complex<double> (0, 1) *
      amp[331] + std::complex<double> (0, 1) * amp[332] - amp[334] - amp[335] -
      amp[336] + std::complex<double> (0, 1) * amp[337] - amp[338] + amp[367] +
      amp[371] + std::complex<double> (0, 1) * amp[372] + amp[376] + amp[377] +
      amp[378] + std::complex<double> (0, 1) * amp[379] + std::complex<double>
      (0, 1) * amp[380] + std::complex<double> (0, 1) * amp[382] +
      std::complex<double> (0, 1) * amp[383] + std::complex<double> (0, 1) *
      amp[384] - std::complex<double> (0, 1) * amp[386] - std::complex<double>
      (0, 1) * amp[387] - std::complex<double> (0, 1) * amp[389] -
      std::complex<double> (0, 1) * amp[394] - std::complex<double> (0, 1) *
      amp[395] - std::complex<double> (0, 1) * amp[397] - std::complex<double>
      (0, 1) * amp[398] - std::complex<double> (0, 1) * amp[402] -
      std::complex<double> (0, 1) * amp[403] - std::complex<double> (0, 1) *
      amp[404] - std::complex<double> (0, 1) * amp[406] - std::complex<double>
      (0, 1) * amp[407] - std::complex<double> (0, 1) * amp[409] +
      std::complex<double> (0, 1) * amp[436] + amp[437] - amp[442] - amp[443] -
      amp[448] - amp[449] - amp[450] - amp[452] - std::complex<double> (0, 1) *
      amp[454] - std::complex<double> (0, 1) * amp[455] - std::complex<double>
      (0, 1) * amp[456] - std::complex<double> (0, 1) * amp[457] -
      std::complex<double> (0, 1) * amp[458] + std::complex<double> (0, 1) *
      amp[499] - std::complex<double> (0, 1) * amp[501] + std::complex<double>
      (0, 1) * amp[502] - std::complex<double> (0, 1) * amp[503] +
      std::complex<double> (0, 1) * amp[507] - std::complex<double> (0, 1) *
      amp[509] + std::complex<double> (0, 1) * amp[510] - std::complex<double>
      (0, 1) * amp[512] + std::complex<double> (0, 1) * amp[513] +
      std::complex<double> (0, 1) * amp[514] + std::complex<double> (0, 1) *
      amp[515] - std::complex<double> (0, 1) * amp[516] - std::complex<double>
      (0, 1) * amp[517] - std::complex<double> (0, 1) * amp[518] +
      std::complex<double> (0, 1) * amp[519] - std::complex<double> (0, 1) *
      amp[521] + std::complex<double> (0, 1) * amp[522] - std::complex<double>
      (0, 1) * amp[523] + amp[538] + amp[540] + amp[542] + amp[543] + amp[544]
      + amp[548] + std::complex<double> (0, 1) * amp[550] -
      std::complex<double> (0, 1) * amp[553] - std::complex<double> (0, 1) *
      amp[554] - amp[558] - amp[559] - amp[560] - amp[562] - amp[564] -
      std::complex<double> (0, 1) * amp[678] + std::complex<double> (0, 1) *
      amp[680] - std::complex<double> (0, 1) * amp[681] - std::complex<double>
      (0, 1) * amp[686] + std::complex<double> (0, 1) * amp[688] -
      std::complex<double> (0, 1) * amp[689] + std::complex<double> (0, 1) *
      amp[691] - std::complex<double> (0, 1) * amp[692] - std::complex<double>
      (0, 1) * amp[693] - std::complex<double> (0, 1) * amp[694] -
      std::complex<double> (0, 1) * amp[698] + std::complex<double> (0, 1) *
      amp[700] - std::complex<double> (0, 1) * amp[701] - amp[717] +
      std::complex<double> (0, 1) * amp[718] - amp[719] - amp[721] - amp[722] -
      amp[723] + std::complex<double> (0, 1) * amp[724] + std::complex<double>
      (0, 1) * amp[725] + std::complex<double> (0, 1) * amp[726] - amp[727] +
      std::complex<double> (0, 1) * amp[728] + std::complex<double> (0, 1) *
      amp[729] - std::complex<double> (0, 1) * amp[828] - std::complex<double>
      (0, 1) * amp[827] - std::complex<double> (0, 1) * amp[834] -
      std::complex<double> (0, 1) * amp[833] - std::complex<double> (0, 1) *
      amp[837] - std::complex<double> (0, 1) * amp[836] - std::complex<double>
      (0, 1) * amp[840] - std::complex<double> (0, 1) * amp[839] - amp[856] -
      amp[855] - amp[859] - amp[858] - amp[862] - amp[861] - amp[865] -
      amp[864] - std::complex<double> (0, 1) * amp[869] - std::complex<double>
      (0, 1) * amp[868] + std::complex<double> (0, 1) * amp[883] + amp[884] +
      amp[891] + amp[890] + amp[894] + amp[893] + amp[897] + amp[896] +
      amp[900] + amp[899] + std::complex<double> (0, 1) * amp[916] +
      std::complex<double> (0, 1) * amp[917] + std::complex<double> (0, 1) *
      amp[922] + std::complex<double> (0, 1) * amp[923] + std::complex<double>
      (0, 1) * amp[925] + std::complex<double> (0, 1) * amp[926] +
      std::complex<double> (0, 1) * amp[928] + std::complex<double> (0, 1) *
      amp[929] - std::complex<double> (0, 1) * amp[932] - std::complex<double>
      (0, 1) * amp[931] + std::complex<double> (0, 1) * amp[935] +
      std::complex<double> (0, 1) * amp[975] + std::complex<double> (0, 1) *
      amp[980] - std::complex<double> (0, 1) * amp[960] - std::complex<double>
      (0, 1) * amp[971] + std::complex<double> (0, 1) * amp[976] +
      std::complex<double> (0, 1) * amp[981] - std::complex<double> (0, 1) *
      amp[961] - std::complex<double> (0, 1) * amp[962] - std::complex<double>
      (0, 1) * amp[963] + std::complex<double> (0, 1) * amp[973];
  jamp[20] = +2. * (-amp[326] - amp[339] + std::complex<double> (0, 1) *
      amp[364] + std::complex<double> (0, 1) * amp[365] - amp[423]);
  jamp[21] = -std::complex<double> (0, 1) * amp[220] - amp[221] - amp[222] -
      amp[223] - amp[227] - amp[229] - std::complex<double> (0, 1) * amp[232] -
      std::complex<double> (0, 1) * amp[233] - std::complex<double> (0, 1) *
      amp[234] - std::complex<double> (0, 1) * amp[235] - std::complex<double>
      (0, 1) * amp[236] + amp[257] + std::complex<double> (0, 1) * amp[258] +
      amp[259] + amp[261] + amp[262] + amp[263] + std::complex<double> (0, 1) *
      amp[264] + std::complex<double> (0, 1) * amp[265] + std::complex<double>
      (0, 1) * amp[266] + amp[267] + std::complex<double> (0, 1) * amp[268] -
      amp[299] - amp[300] - amp[305] - amp[306] - amp[307] - amp[309] +
      std::complex<double> (0, 1) * amp[311] - amp[313] - amp[312] - amp[316] -
      amp[315] - amp[319] - amp[318] - amp[323] - amp[322] -
      std::complex<double> (0, 1) * amp[327] + amp[328] + amp[329] -
      std::complex<double> (0, 1) * amp[330] - std::complex<double> (0, 1) *
      amp[331] - std::complex<double> (0, 1) * amp[332] + amp[334] + amp[335] +
      amp[336] - std::complex<double> (0, 1) * amp[337] + amp[338] - amp[367] -
      amp[371] - std::complex<double> (0, 1) * amp[372] - amp[376] - amp[377] -
      amp[378] - std::complex<double> (0, 1) * amp[379] - std::complex<double>
      (0, 1) * amp[380] - std::complex<double> (0, 1) * amp[382] -
      std::complex<double> (0, 1) * amp[383] - std::complex<double> (0, 1) *
      amp[384] + std::complex<double> (0, 1) * amp[386] + std::complex<double>
      (0, 1) * amp[387] + std::complex<double> (0, 1) * amp[389] +
      std::complex<double> (0, 1) * amp[394] + std::complex<double> (0, 1) *
      amp[395] + std::complex<double> (0, 1) * amp[397] + std::complex<double>
      (0, 1) * amp[398] + std::complex<double> (0, 1) * amp[402] +
      std::complex<double> (0, 1) * amp[403] + std::complex<double> (0, 1) *
      amp[404] + std::complex<double> (0, 1) * amp[406] + std::complex<double>
      (0, 1) * amp[407] + std::complex<double> (0, 1) * amp[409] -
      std::complex<double> (0, 1) * amp[436] - amp[437] + amp[578] + amp[579] +
      amp[584] + amp[585] + amp[586] + amp[588] - std::complex<double> (0, 1) *
      amp[679] - std::complex<double> (0, 1) * amp[680] - std::complex<double>
      (0, 1) * amp[682] - std::complex<double> (0, 1) * amp[687] -
      std::complex<double> (0, 1) * amp[688] - std::complex<double> (0, 1) *
      amp[690] - std::complex<double> (0, 1) * amp[691] - std::complex<double>
      (0, 1) * amp[695] - std::complex<double> (0, 1) * amp[696] -
      std::complex<double> (0, 1) * amp[697] - std::complex<double> (0, 1) *
      amp[699] - std::complex<double> (0, 1) * amp[700] - std::complex<double>
      (0, 1) * amp[702] + std::complex<double> (0, 1) * amp[734] +
      std::complex<double> (0, 1) * amp[735] + std::complex<double> (0, 1) *
      amp[736] + std::complex<double> (0, 1) * amp[737] + std::complex<double>
      (0, 1) * amp[739] + std::complex<double> (0, 1) * amp[740] +
      std::complex<double> (0, 1) * amp[742] + std::complex<double> (0, 1) *
      amp[743] + std::complex<double> (0, 1) * amp[744] + std::complex<double>
      (0, 1) * amp[745] + std::complex<double> (0, 1) * amp[746] +
      std::complex<double> (0, 1) * amp[747] + std::complex<double> (0, 1) *
      amp[748] + std::complex<double> (0, 1) * amp[749] + std::complex<double>
      (0, 1) * amp[751] + std::complex<double> (0, 1) * amp[752] +
      std::complex<double> (0, 1) * amp[753] + std::complex<double> (0, 1) *
      amp[754] - std::complex<double> (0, 1) * amp[872] + std::complex<double>
      (0, 1) * amp[873] - std::complex<double> (0, 1) * amp[875] -
      std::complex<double> (0, 1) * amp[876] - std::complex<double> (0, 1) *
      amp[883] - amp[884] - amp[891] + amp[889] - amp[894] + amp[892] -
      amp[897] + amp[895] - amp[900] + amp[898] + std::complex<double> (0, 1) *
      amp[915] - std::complex<double> (0, 1) * amp[917] + std::complex<double>
      (0, 1) * amp[921] - std::complex<double> (0, 1) * amp[923] +
      std::complex<double> (0, 1) * amp[924] - std::complex<double> (0, 1) *
      amp[926] + std::complex<double> (0, 1) * amp[927] - std::complex<double>
      (0, 1) * amp[929] + std::complex<double> (0, 1) * amp[932] -
      std::complex<double> (0, 1) * amp[930] - std::complex<double> (0, 1) *
      amp[933] - std::complex<double> (0, 1) * amp[939] - std::complex<double>
      (0, 1) * amp[938] - std::complex<double> (0, 1) * amp[945] -
      std::complex<double> (0, 1) * amp[944] - std::complex<double> (0, 1) *
      amp[948] - std::complex<double> (0, 1) * amp[947] - std::complex<double>
      (0, 1) * amp[951] - std::complex<double> (0, 1) * amp[950] -
      std::complex<double> (0, 1) * amp[959] + std::complex<double> (0, 1) *
      amp[985] - std::complex<double> (0, 1) * amp[980] + std::complex<double>
      (0, 1) * amp[971] + std::complex<double> (0, 1) * amp[961] -
      std::complex<double> (0, 1) * amp[956] + std::complex<double> (0, 1) *
      amp[982] - std::complex<double> (0, 1) * amp[957] - std::complex<double>
      (0, 1) * amp[973] + std::complex<double> (0, 1) * amp[984];
  jamp[22] = +2. * (+amp[5] + amp[6] - amp[7] + amp[26] - amp[495] - amp[504] -
      amp[505] + amp[506] + amp[631] + amp[632] + amp[633] + amp[730] +
      amp[731] - amp[732] + amp[786] - amp[788] + amp[831] + amp[832] -
      amp[919] - amp[918] + amp[941] - amp[943] - amp[955] + amp[953]);
  jamp[23] = +2. * (-amp[5] - amp[6] + amp[7] - amp[26] - amp[381] - amp[390] -
      amp[391] + amp[392] + amp[683] + amp[684] + amp[685] - amp[730] -
      amp[731] + amp[732] + amp[787] + amp[788] + amp[830] - amp[832] +
      amp[918] - amp[920] - amp[942] - amp[941] - amp[954] - amp[953]);
  jamp[24] = +2. * (+amp[381] + amp[390] + amp[391] - amp[392] + amp[495] +
      amp[504] + amp[505] - amp[506] - amp[631] - amp[632] - amp[633] -
      amp[683] - amp[684] - amp[685] - amp[787] - amp[786] - amp[831] -
      amp[830] + amp[919] + amp[920] + amp[942] + amp[943] + amp[955] +
      amp[954]);
  jamp[25] = +2. * (-amp[5] - amp[6] + amp[7] - amp[26] - amp[381] - amp[390] -
      amp[391] + amp[392] + amp[683] + amp[684] + amp[685] - amp[730] -
      amp[731] + amp[732] + amp[787] + amp[788] + amp[830] - amp[832] +
      amp[918] - amp[920] - amp[942] - amp[941] - amp[954] - amp[953]);
  jamp[26] = +2. * (+amp[381] + amp[390] + amp[391] - amp[392] + amp[495] +
      amp[504] + amp[505] - amp[506] - amp[631] - amp[632] - amp[633] -
      amp[683] - amp[684] - amp[685] - amp[787] - amp[786] - amp[831] -
      amp[830] + amp[919] + amp[920] + amp[942] + amp[943] + amp[955] +
      amp[954]);
  jamp[27] = +2. * (+amp[5] + amp[6] - amp[7] + amp[26] - amp[495] - amp[504] -
      amp[505] + amp[506] + amp[631] + amp[632] + amp[633] + amp[730] +
      amp[731] - amp[732] + amp[786] - amp[788] + amp[831] + amp[832] -
      amp[919] - amp[918] + amp[941] - amp[943] - amp[955] + amp[953]);
  jamp[28] = -std::complex<double> (0, 1) * amp[1] - std::complex<double> (0,
      1) * amp[2] - std::complex<double> (0, 1) * amp[3] - std::complex<double>
      (0, 1) * amp[4] - std::complex<double> (0, 1) * amp[9] -
      std::complex<double> (0, 1) * amp[10] - std::complex<double> (0, 1) *
      amp[12] - std::complex<double> (0, 1) * amp[13] - std::complex<double>
      (0, 1) * amp[14] - std::complex<double> (0, 1) * amp[15] -
      std::complex<double> (0, 1) * amp[16] - std::complex<double> (0, 1) *
      amp[17] - std::complex<double> (0, 1) * amp[18] - std::complex<double>
      (0, 1) * amp[19] - std::complex<double> (0, 1) * amp[21] -
      std::complex<double> (0, 1) * amp[22] - std::complex<double> (0, 1) *
      amp[23] - std::complex<double> (0, 1) * amp[24] + amp[33] + amp[34] +
      amp[39] + amp[40] + amp[41] + amp[43] + std::complex<double> (0, 1) *
      amp[45] + std::complex<double> (0, 1) * amp[46] + std::complex<double>
      (0, 1) * amp[47] + std::complex<double> (0, 1) * amp[48] +
      std::complex<double> (0, 1) * amp[49] - amp[84] - amp[85] - amp[90] -
      amp[91] - amp[92] - amp[94] - std::complex<double> (0, 1) * amp[110] -
      std::complex<double> (0, 1) * amp[111] - amp[366] - amp[370] +
      std::complex<double> (0, 1) * amp[372] - amp[373] - amp[374] - amp[375] +
      std::complex<double> (0, 1) * amp[379] + std::complex<double> (0, 1) *
      amp[380] + std::complex<double> (0, 1) * amp[382] + std::complex<double>
      (0, 1) * amp[383] + std::complex<double> (0, 1) * amp[384] -
      std::complex<double> (0, 1) * amp[385] - std::complex<double> (0, 1) *
      amp[386] - std::complex<double> (0, 1) * amp[388] - std::complex<double>
      (0, 1) * amp[393] - std::complex<double> (0, 1) * amp[394] -
      std::complex<double> (0, 1) * amp[396] - std::complex<double> (0, 1) *
      amp[397] - std::complex<double> (0, 1) * amp[399] - std::complex<double>
      (0, 1) * amp[400] - std::complex<double> (0, 1) * amp[401] -
      std::complex<double> (0, 1) * amp[405] - std::complex<double> (0, 1) *
      amp[406] - std::complex<double> (0, 1) * amp[408] - amp[424] -
      std::complex<double> (0, 1) * amp[425] - amp[426] - amp[428] - amp[429] -
      amp[430] - std::complex<double> (0, 1) * amp[431] - std::complex<double>
      (0, 1) * amp[432] - std::complex<double> (0, 1) * amp[433] - amp[434] -
      std::complex<double> (0, 1) * amp[435] - amp[438] + std::complex<double>
      (0, 1) * amp[439] - std::complex<double> (0, 1) * amp[678] +
      std::complex<double> (0, 1) * amp[680] - std::complex<double> (0, 1) *
      amp[681] - std::complex<double> (0, 1) * amp[686] + std::complex<double>
      (0, 1) * amp[688] - std::complex<double> (0, 1) * amp[689] +
      std::complex<double> (0, 1) * amp[691] - std::complex<double> (0, 1) *
      amp[692] - std::complex<double> (0, 1) * amp[693] - std::complex<double>
      (0, 1) * amp[694] - std::complex<double> (0, 1) * amp[698] +
      std::complex<double> (0, 1) * amp[700] - std::complex<double> (0, 1) *
      amp[701] + amp[704] + std::complex<double> (0, 1) * amp[705] + amp[706] +
      amp[708] + amp[709] + amp[710] + std::complex<double> (0, 1) * amp[711] +
      std::complex<double> (0, 1) * amp[712] + std::complex<double> (0, 1) *
      amp[713] + amp[714] + std::complex<double> (0, 1) * amp[715] +
      std::complex<double> (0, 1) * amp[729] + std::complex<double> (0, 1) *
      amp[755] - amp[756] - amp[757] - amp[758] - amp[762] - amp[764] +
      std::complex<double> (0, 1) * amp[784] + std::complex<double> (0, 1) *
      amp[785] + std::complex<double> (0, 1) * amp[790] + std::complex<double>
      (0, 1) * amp[791] + std::complex<double> (0, 1) * amp[793] +
      std::complex<double> (0, 1) * amp[794] + std::complex<double> (0, 1) *
      amp[796] + std::complex<double> (0, 1) * amp[797] + amp[813] + amp[812] +
      amp[816] + amp[815] + amp[819] + amp[818] + amp[822] + amp[821] +
      std::complex<double> (0, 1) * amp[826] + std::complex<double> (0, 1) *
      amp[825] - std::complex<double> (0, 1) * amp[827] + std::complex<double>
      (0, 1) * amp[829] - std::complex<double> (0, 1) * amp[833] +
      std::complex<double> (0, 1) * amp[835] - std::complex<double> (0, 1) *
      amp[836] + std::complex<double> (0, 1) * amp[838] - std::complex<double>
      (0, 1) * amp[839] + std::complex<double> (0, 1) * amp[841] - amp[844] +
      amp[842] - amp[847] + amp[845] - amp[850] + amp[848] - amp[853] +
      amp[851] + std::complex<double> (0, 1) * amp[870] - std::complex<double>
      (0, 1) * amp[868] - amp[878] - std::complex<double> (0, 1) * amp[879] +
      std::complex<double> (0, 1) * amp[937] + std::complex<double> (0, 1) *
      amp[964] - std::complex<double> (0, 1) * amp[970] + std::complex<double>
      (0, 1) * amp[980] - std::complex<double> (0, 1) * amp[971] +
      std::complex<double> (0, 1) * amp[981] + std::complex<double> (0, 1) *
      amp[966] - std::complex<double> (0, 1) * amp[962] + std::complex<double>
      (0, 1) * amp[967] - std::complex<double> (0, 1) * amp[972] -
      std::complex<double> (0, 1) * amp[969];
  jamp[29] = +2. * (-amp[256] - amp[269] - std::complex<double> (0, 1) *
      amp[368] - std::complex<double> (0, 1) * amp[369] - amp[703]);
  jamp[30] = -std::complex<double> (0, 1) * amp[346] - amp[347] - amp[348] -
      std::complex<double> (0, 1) * amp[349] - std::complex<double> (0, 1) *
      amp[350] - std::complex<double> (0, 1) * amp[351] - amp[353] - amp[354] -
      amp[355] - std::complex<double> (0, 1) * amp[356] - amp[357] + amp[366] +
      amp[370] - std::complex<double> (0, 1) * amp[372] + amp[373] + amp[374] +
      amp[375] - std::complex<double> (0, 1) * amp[379] - std::complex<double>
      (0, 1) * amp[380] - std::complex<double> (0, 1) * amp[382] -
      std::complex<double> (0, 1) * amp[383] - std::complex<double> (0, 1) *
      amp[384] + std::complex<double> (0, 1) * amp[386] + std::complex<double>
      (0, 1) * amp[387] + std::complex<double> (0, 1) * amp[389] +
      std::complex<double> (0, 1) * amp[394] + std::complex<double> (0, 1) *
      amp[395] + std::complex<double> (0, 1) * amp[397] + std::complex<double>
      (0, 1) * amp[398] + std::complex<double> (0, 1) * amp[402] +
      std::complex<double> (0, 1) * amp[403] + std::complex<double> (0, 1) *
      amp[404] + std::complex<double> (0, 1) * amp[406] + std::complex<double>
      (0, 1) * amp[407] + std::complex<double> (0, 1) * amp[409] -
      std::complex<double> (0, 1) * amp[436] + amp[438] - amp[461] - amp[462] -
      amp[467] - amp[468] - amp[469] - amp[471] + std::complex<double> (0, 1) *
      amp[473] + std::complex<double> (0, 1) * amp[474] + std::complex<double>
      (0, 1) * amp[475] + std::complex<double> (0, 1) * amp[476] +
      std::complex<double> (0, 1) * amp[477] - std::complex<double> (0, 1) *
      amp[499] + std::complex<double> (0, 1) * amp[501] - std::complex<double>
      (0, 1) * amp[502] + std::complex<double> (0, 1) * amp[503] -
      std::complex<double> (0, 1) * amp[507] + std::complex<double> (0, 1) *
      amp[509] - std::complex<double> (0, 1) * amp[510] + std::complex<double>
      (0, 1) * amp[512] - std::complex<double> (0, 1) * amp[513] -
      std::complex<double> (0, 1) * amp[514] - std::complex<double> (0, 1) *
      amp[515] + std::complex<double> (0, 1) * amp[516] + std::complex<double>
      (0, 1) * amp[517] + std::complex<double> (0, 1) * amp[518] -
      std::complex<double> (0, 1) * amp[519] + std::complex<double> (0, 1) *
      amp[521] - std::complex<double> (0, 1) * amp[522] + std::complex<double>
      (0, 1) * amp[523] + amp[525] + amp[527] + amp[529] + amp[530] + amp[531]
      + amp[535] - std::complex<double> (0, 1) * amp[550] +
      std::complex<double> (0, 1) * amp[553] + std::complex<double> (0, 1) *
      amp[590] - amp[594] - amp[595] - amp[596] - amp[598] - amp[600] +
      std::complex<double> (0, 1) * amp[678] - std::complex<double> (0, 1) *
      amp[680] + std::complex<double> (0, 1) * amp[681] + std::complex<double>
      (0, 1) * amp[686] - std::complex<double> (0, 1) * amp[688] +
      std::complex<double> (0, 1) * amp[689] - std::complex<double> (0, 1) *
      amp[691] + std::complex<double> (0, 1) * amp[692] + std::complex<double>
      (0, 1) * amp[693] + std::complex<double> (0, 1) * amp[694] +
      std::complex<double> (0, 1) * amp[698] - std::complex<double> (0, 1) *
      amp[700] + std::complex<double> (0, 1) * amp[701] - amp[704] -
      std::complex<double> (0, 1) * amp[705] - amp[706] - amp[708] - amp[709] -
      amp[710] - std::complex<double> (0, 1) * amp[711] - std::complex<double>
      (0, 1) * amp[712] - std::complex<double> (0, 1) * amp[713] - amp[714] -
      std::complex<double> (0, 1) * amp[715] - std::complex<double> (0, 1) *
      amp[729] + std::complex<double> (0, 1) * amp[828] + std::complex<double>
      (0, 1) * amp[827] + std::complex<double> (0, 1) * amp[834] +
      std::complex<double> (0, 1) * amp[833] + std::complex<double> (0, 1) *
      amp[837] + std::complex<double> (0, 1) * amp[836] + std::complex<double>
      (0, 1) * amp[840] + std::complex<double> (0, 1) * amp[839] - amp[843] -
      amp[842] - amp[846] - amp[845] - amp[849] - amp[848] - amp[852] -
      amp[851] + std::complex<double> (0, 1) * amp[869] + std::complex<double>
      (0, 1) * amp[868] - std::complex<double> (0, 1) * amp[877] + amp[878] +
      amp[904] + amp[903] + amp[907] + amp[906] + amp[910] + amp[909] +
      amp[913] + amp[912] - std::complex<double> (0, 1) * amp[916] -
      std::complex<double> (0, 1) * amp[917] - std::complex<double> (0, 1) *
      amp[922] - std::complex<double> (0, 1) * amp[923] - std::complex<double>
      (0, 1) * amp[925] - std::complex<double> (0, 1) * amp[926] -
      std::complex<double> (0, 1) * amp[928] - std::complex<double> (0, 1) *
      amp[929] + std::complex<double> (0, 1) * amp[932] + std::complex<double>
      (0, 1) * amp[931] - std::complex<double> (0, 1) * amp[937] -
      std::complex<double> (0, 1) * amp[975] - std::complex<double> (0, 1) *
      amp[980] + std::complex<double> (0, 1) * amp[960] + std::complex<double>
      (0, 1) * amp[971] - std::complex<double> (0, 1) * amp[976] -
      std::complex<double> (0, 1) * amp[981] + std::complex<double> (0, 1) *
      amp[961] + std::complex<double> (0, 1) * amp[962] + std::complex<double>
      (0, 1) * amp[963] - std::complex<double> (0, 1) * amp[973];
  jamp[31] = +std::complex<double> (0, 1) * amp[1] + std::complex<double> (0,
      1) * amp[2] + std::complex<double> (0, 1) * amp[3] + std::complex<double>
      (0, 1) * amp[4] + std::complex<double> (0, 1) * amp[9] +
      std::complex<double> (0, 1) * amp[10] + std::complex<double> (0, 1) *
      amp[12] + std::complex<double> (0, 1) * amp[13] + std::complex<double>
      (0, 1) * amp[14] + std::complex<double> (0, 1) * amp[15] +
      std::complex<double> (0, 1) * amp[16] + std::complex<double> (0, 1) *
      amp[17] + std::complex<double> (0, 1) * amp[18] + std::complex<double>
      (0, 1) * amp[19] + std::complex<double> (0, 1) * amp[21] +
      std::complex<double> (0, 1) * amp[22] + std::complex<double> (0, 1) *
      amp[23] + std::complex<double> (0, 1) * amp[24] - amp[33] - amp[34] -
      amp[39] - amp[40] - amp[41] - amp[43] - std::complex<double> (0, 1) *
      amp[45] - std::complex<double> (0, 1) * amp[46] - std::complex<double>
      (0, 1) * amp[47] - std::complex<double> (0, 1) * amp[48] -
      std::complex<double> (0, 1) * amp[49] + amp[84] + amp[85] + amp[90] +
      amp[91] + amp[92] + amp[94] + std::complex<double> (0, 1) * amp[110] +
      std::complex<double> (0, 1) * amp[111] - amp[481] - amp[485] -
      std::complex<double> (0, 1) * amp[486] - amp[490] - amp[491] - amp[492] -
      std::complex<double> (0, 1) * amp[493] - std::complex<double> (0, 1) *
      amp[494] - std::complex<double> (0, 1) * amp[496] - std::complex<double>
      (0, 1) * amp[497] - std::complex<double> (0, 1) * amp[498] +
      std::complex<double> (0, 1) * amp[499] + std::complex<double> (0, 1) *
      amp[500] + std::complex<double> (0, 1) * amp[502] + std::complex<double>
      (0, 1) * amp[507] + std::complex<double> (0, 1) * amp[508] +
      std::complex<double> (0, 1) * amp[510] + std::complex<double> (0, 1) *
      amp[511] + std::complex<double> (0, 1) * amp[513] + std::complex<double>
      (0, 1) * amp[514] + std::complex<double> (0, 1) * amp[515] +
      std::complex<double> (0, 1) * amp[519] + std::complex<double> (0, 1) *
      amp[520] + std::complex<double> (0, 1) * amp[522] - amp[525] +
      std::complex<double> (0, 1) * amp[526] - amp[527] - amp[529] - amp[530] -
      amp[531] + std::complex<double> (0, 1) * amp[532] + std::complex<double>
      (0, 1) * amp[533] + std::complex<double> (0, 1) * amp[534] - amp[535] +
      std::complex<double> (0, 1) * amp[536] - amp[551] - std::complex<double>
      (0, 1) * amp[553] + std::complex<double> (0, 1) * amp[626] -
      std::complex<double> (0, 1) * amp[628] + std::complex<double> (0, 1) *
      amp[629] + std::complex<double> (0, 1) * amp[634] - std::complex<double>
      (0, 1) * amp[636] + std::complex<double> (0, 1) * amp[637] -
      std::complex<double> (0, 1) * amp[639] + std::complex<double> (0, 1) *
      amp[640] + std::complex<double> (0, 1) * amp[641] + std::complex<double>
      (0, 1) * amp[642] + std::complex<double> (0, 1) * amp[646] -
      std::complex<double> (0, 1) * amp[648] + std::complex<double> (0, 1) *
      amp[649] + amp[665] - std::complex<double> (0, 1) * amp[666] + amp[667] +
      amp[669] + amp[670] + amp[671] - std::complex<double> (0, 1) * amp[672] -
      std::complex<double> (0, 1) * amp[673] - std::complex<double> (0, 1) *
      amp[674] + amp[675] - std::complex<double> (0, 1) * amp[676] -
      std::complex<double> (0, 1) * amp[677] - std::complex<double> (0, 1) *
      amp[755] - amp[759] - amp[760] - amp[761] - amp[763] - amp[765] +
      std::complex<double> (0, 1) * amp[783] - std::complex<double> (0, 1) *
      amp[785] + std::complex<double> (0, 1) * amp[789] - std::complex<double>
      (0, 1) * amp[791] + std::complex<double> (0, 1) * amp[792] -
      std::complex<double> (0, 1) * amp[794] + std::complex<double> (0, 1) *
      amp[795] - std::complex<double> (0, 1) * amp[797] - amp[813] + amp[811] -
      amp[816] + amp[814] - amp[819] + amp[817] - amp[822] + amp[820] -
      std::complex<double> (0, 1) * amp[826] + std::complex<double> (0, 1) *
      amp[824] - std::complex<double> (0, 1) * amp[828] - std::complex<double>
      (0, 1) * amp[829] - std::complex<double> (0, 1) * amp[834] -
      std::complex<double> (0, 1) * amp[835] - std::complex<double> (0, 1) *
      amp[837] - std::complex<double> (0, 1) * amp[838] - std::complex<double>
      (0, 1) * amp[840] - std::complex<double> (0, 1) * amp[841] + amp[844] +
      amp[843] + amp[847] + amp[846] + amp[850] + amp[849] + amp[853] +
      amp[852] - std::complex<double> (0, 1) * amp[870] - std::complex<double>
      (0, 1) * amp[869] - amp[887] + std::complex<double> (0, 1) * amp[888] -
      std::complex<double> (0, 1) * amp[934] + std::complex<double> (0, 1) *
      amp[974] - std::complex<double> (0, 1) * amp[964] + std::complex<double>
      (0, 1) * amp[970] + std::complex<double> (0, 1) * amp[975] +
      std::complex<double> (0, 1) * amp[965] - std::complex<double> (0, 1) *
      amp[977] - std::complex<double> (0, 1) * amp[967] - std::complex<double>
      (0, 1) * amp[978] - std::complex<double> (0, 1) * amp[963] +
      std::complex<double> (0, 1) * amp[969];
  jamp[32] = +2. * (-amp[459] - amp[472] - std::complex<double> (0, 1) *
      amp[478] - std::complex<double> (0, 1) * amp[479] - amp[524]);
  jamp[33] = -amp[238] - std::complex<double> (0, 1) * amp[239] - amp[240] -
      amp[242] - amp[243] - amp[244] - std::complex<double> (0, 1) * amp[245] -
      std::complex<double> (0, 1) * amp[246] - std::complex<double> (0, 1) *
      amp[247] - amp[248] - std::complex<double> (0, 1) * amp[249] - amp[257] -
      amp[259] - amp[261] - amp[262] - amp[263] - amp[267] +
      std::complex<double> (0, 1) * amp[270] + std::complex<double> (0, 1) *
      amp[271] + std::complex<double> (0, 1) * amp[272] + std::complex<double>
      (0, 1) * amp[273] + std::complex<double> (0, 1) * amp[274] +
      std::complex<double> (0, 1) * amp[275] - amp[279] - amp[280] - amp[281] -
      amp[283] - amp[285] - std::complex<double> (0, 1) * amp[310] + amp[314] +
      amp[313] + amp[317] + amp[316] + amp[320] + amp[319] + amp[324] +
      amp[323] + amp[481] + amp[485] + std::complex<double> (0, 1) * amp[486] +
      amp[490] + amp[491] + amp[492] + std::complex<double> (0, 1) * amp[493] +
      std::complex<double> (0, 1) * amp[494] + std::complex<double> (0, 1) *
      amp[496] + std::complex<double> (0, 1) * amp[497] + std::complex<double>
      (0, 1) * amp[498] - std::complex<double> (0, 1) * amp[499] -
      std::complex<double> (0, 1) * amp[500] - std::complex<double> (0, 1) *
      amp[502] - std::complex<double> (0, 1) * amp[507] - std::complex<double>
      (0, 1) * amp[508] - std::complex<double> (0, 1) * amp[510] -
      std::complex<double> (0, 1) * amp[511] - std::complex<double> (0, 1) *
      amp[513] - std::complex<double> (0, 1) * amp[514] - std::complex<double>
      (0, 1) * amp[515] - std::complex<double> (0, 1) * amp[519] -
      std::complex<double> (0, 1) * amp[520] - std::complex<double> (0, 1) *
      amp[522] + amp[525] - std::complex<double> (0, 1) * amp[526] + amp[527] +
      amp[529] + amp[530] + amp[531] - std::complex<double> (0, 1) * amp[532] -
      std::complex<double> (0, 1) * amp[533] - std::complex<double> (0, 1) *
      amp[534] + amp[535] - std::complex<double> (0, 1) * amp[536] + amp[551] +
      std::complex<double> (0, 1) * amp[553] + std::complex<double> (0, 1) *
      amp[627] + std::complex<double> (0, 1) * amp[628] + std::complex<double>
      (0, 1) * amp[630] + std::complex<double> (0, 1) * amp[635] +
      std::complex<double> (0, 1) * amp[636] + std::complex<double> (0, 1) *
      amp[638] + std::complex<double> (0, 1) * amp[639] + std::complex<double>
      (0, 1) * amp[643] + std::complex<double> (0, 1) * amp[644] +
      std::complex<double> (0, 1) * amp[645] + std::complex<double> (0, 1) *
      amp[647] + std::complex<double> (0, 1) * amp[648] + std::complex<double>
      (0, 1) * amp[650] + std::complex<double> (0, 1) * amp[678] +
      std::complex<double> (0, 1) * amp[679] + std::complex<double> (0, 1) *
      amp[681] + std::complex<double> (0, 1) * amp[682] + std::complex<double>
      (0, 1) * amp[686] + std::complex<double> (0, 1) * amp[687] +
      std::complex<double> (0, 1) * amp[689] + std::complex<double> (0, 1) *
      amp[690] + std::complex<double> (0, 1) * amp[692] + std::complex<double>
      (0, 1) * amp[693] + std::complex<double> (0, 1) * amp[694] +
      std::complex<double> (0, 1) * amp[695] + std::complex<double> (0, 1) *
      amp[696] + std::complex<double> (0, 1) * amp[697] + std::complex<double>
      (0, 1) * amp[698] + std::complex<double> (0, 1) * amp[699] +
      std::complex<double> (0, 1) * amp[701] + std::complex<double> (0, 1) *
      amp[702] - amp[704] - amp[706] - amp[708] - amp[709] - amp[710] -
      amp[714] - std::complex<double> (0, 1) * amp[729] + std::complex<double>
      (0, 1) * amp[828] + std::complex<double> (0, 1) * amp[827] +
      std::complex<double> (0, 1) * amp[834] + std::complex<double> (0, 1) *
      amp[833] + std::complex<double> (0, 1) * amp[837] + std::complex<double>
      (0, 1) * amp[836] + std::complex<double> (0, 1) * amp[840] +
      std::complex<double> (0, 1) * amp[839] - amp[843] - amp[842] - amp[846] -
      amp[845] - amp[849] - amp[848] - amp[852] - amp[851] +
      std::complex<double> (0, 1) * amp[869] + std::complex<double> (0, 1) *
      amp[868] + std::complex<double> (0, 1) * amp[871] + std::complex<double>
      (0, 1) * amp[872] + std::complex<double> (0, 1) * amp[875] +
      std::complex<double> (0, 1) * amp[874] + amp[887] - std::complex<double>
      (0, 1) * amp[888] + std::complex<double> (0, 1) * amp[939] +
      std::complex<double> (0, 1) * amp[940] + std::complex<double> (0, 1) *
      amp[945] + std::complex<double> (0, 1) * amp[946] + std::complex<double>
      (0, 1) * amp[948] + std::complex<double> (0, 1) * amp[949] +
      std::complex<double> (0, 1) * amp[951] + std::complex<double> (0, 1) *
      amp[952] - std::complex<double> (0, 1) * amp[974] - std::complex<double>
      (0, 1) * amp[975] - std::complex<double> (0, 1) * amp[981] +
      std::complex<double> (0, 1) * amp[977] - std::complex<double> (0, 1) *
      amp[982] + std::complex<double> (0, 1) * amp[962] + std::complex<double>
      (0, 1) * amp[957] + std::complex<double> (0, 1) * amp[963] +
      std::complex<double> (0, 1) * amp[958] - std::complex<double> (0, 1) *
      amp[979];
  jamp[34] = +2. * (+std::complex<double> (0, 1) * amp[38] +
      std::complex<double> (0, 1) * amp[76] - std::complex<double> (0, 1) *
      amp[528] - std::complex<double> (0, 1) * amp[541] + std::complex<double>
      (0, 1) * amp[707] + std::complex<double> (0, 1) * amp[720] +
      std::complex<double> (0, 1) * amp[854] + std::complex<double> (0, 1) *
      amp[867]);
  jamp[35] = +2. * (-std::complex<double> (0, 1) * amp[38] -
      std::complex<double> (0, 1) * amp[76] + std::complex<double> (0, 1) *
      amp[528] + std::complex<double> (0, 1) * amp[541] - std::complex<double>
      (0, 1) * amp[707] - std::complex<double> (0, 1) * amp[720] -
      std::complex<double> (0, 1) * amp[854] - std::complex<double> (0, 1) *
      amp[867]);
  jamp[36] = -std::complex<double> (0, 1) * amp[0] + std::complex<double> (0,
      1) * amp[2] + std::complex<double> (0, 1) * amp[4] - std::complex<double>
      (0, 1) * amp[8] + std::complex<double> (0, 1) * amp[10] -
      std::complex<double> (0, 1) * amp[11] + std::complex<double> (0, 1) *
      amp[13] + std::complex<double> (0, 1) * amp[17] + std::complex<double>
      (0, 1) * amp[18] + std::complex<double> (0, 1) * amp[19] -
      std::complex<double> (0, 1) * amp[20] + std::complex<double> (0, 1) *
      amp[22] + std::complex<double> (0, 1) * amp[24] - std::complex<double>
      (0, 1) * amp[25] - std::complex<double> (0, 1) * amp[27] -
      std::complex<double> (0, 1) * amp[28] - std::complex<double> (0, 1) *
      amp[29] - std::complex<double> (0, 1) * amp[30] + std::complex<double>
      (0, 1) * amp[32] - amp[33] - amp[34] + std::complex<double> (0, 1) *
      amp[35] + std::complex<double> (0, 1) * amp[36] + std::complex<double>
      (0, 1) * amp[37] - amp[39] - amp[40] - amp[41] + std::complex<double> (0,
      1) * amp[42] - amp[43] + amp[98] + amp[102] + std::complex<double> (0, 1)
      * amp[103] + amp[107] + amp[108] + amp[109] + std::complex<double> (0, 1)
      * amp[110] + amp[112] + amp[461] + amp[462] + amp[467] + amp[468] +
      amp[469] + amp[471] - std::complex<double> (0, 1) * amp[473] -
      std::complex<double> (0, 1) * amp[474] - std::complex<double> (0, 1) *
      amp[475] - std::complex<double> (0, 1) * amp[476] - std::complex<double>
      (0, 1) * amp[477] + std::complex<double> (0, 1) * amp[499] -
      std::complex<double> (0, 1) * amp[501] + std::complex<double> (0, 1) *
      amp[502] - std::complex<double> (0, 1) * amp[503] + std::complex<double>
      (0, 1) * amp[507] - std::complex<double> (0, 1) * amp[509] +
      std::complex<double> (0, 1) * amp[510] - std::complex<double> (0, 1) *
      amp[512] + std::complex<double> (0, 1) * amp[513] + std::complex<double>
      (0, 1) * amp[514] + std::complex<double> (0, 1) * amp[515] -
      std::complex<double> (0, 1) * amp[516] - std::complex<double> (0, 1) *
      amp[517] - std::complex<double> (0, 1) * amp[518] + std::complex<double>
      (0, 1) * amp[519] - std::complex<double> (0, 1) * amp[521] +
      std::complex<double> (0, 1) * amp[522] - std::complex<double> (0, 1) *
      amp[523] - amp[525] - amp[527] - amp[529] - amp[530] - amp[531] -
      amp[535] + std::complex<double> (0, 1) * amp[550] - std::complex<double>
      (0, 1) * amp[553] - std::complex<double> (0, 1) * amp[590] - amp[591] -
      amp[592] - amp[593] - amp[597] - amp[599] + std::complex<double> (0, 1) *
      amp[613] - amp[614] - amp[615] + std::complex<double> (0, 1) * amp[616] +
      std::complex<double> (0, 1) * amp[617] + std::complex<double> (0, 1) *
      amp[618] - amp[620] - amp[621] - amp[622] + std::complex<double> (0, 1) *
      amp[623] - amp[624] + std::complex<double> (0, 1) * amp[625] -
      std::complex<double> (0, 1) * amp[733] + std::complex<double> (0, 1) *
      amp[735] + std::complex<double> (0, 1) * amp[737] - std::complex<double>
      (0, 1) * amp[738] + std::complex<double> (0, 1) * amp[740] -
      std::complex<double> (0, 1) * amp[741] + std::complex<double> (0, 1) *
      amp[743] + std::complex<double> (0, 1) * amp[747] + std::complex<double>
      (0, 1) * amp[748] + std::complex<double> (0, 1) * amp[749] -
      std::complex<double> (0, 1) * amp[750] + std::complex<double> (0, 1) *
      amp[752] + std::complex<double> (0, 1) * amp[754] + std::complex<double>
      (0, 1) * amp[780] + amp[782] - std::complex<double> (0, 1) * amp[828] -
      std::complex<double> (0, 1) * amp[829] - std::complex<double> (0, 1) *
      amp[834] - std::complex<double> (0, 1) * amp[835] - std::complex<double>
      (0, 1) * amp[837] - std::complex<double> (0, 1) * amp[838] -
      std::complex<double> (0, 1) * amp[840] - std::complex<double> (0, 1) *
      amp[841] + amp[844] + amp[843] + amp[847] + amp[846] + amp[850] +
      amp[849] + amp[853] + amp[852] - std::complex<double> (0, 1) * amp[870] -
      std::complex<double> (0, 1) * amp[869] - amp[903] - amp[902] - amp[906] -
      amp[905] - amp[909] - amp[908] - amp[912] - amp[911] +
      std::complex<double> (0, 1) * amp[916] + std::complex<double> (0, 1) *
      amp[915] + std::complex<double> (0, 1) * amp[922] + std::complex<double>
      (0, 1) * amp[921] + std::complex<double> (0, 1) * amp[925] +
      std::complex<double> (0, 1) * amp[924] + std::complex<double> (0, 1) *
      amp[928] + std::complex<double> (0, 1) * amp[927] - std::complex<double>
      (0, 1) * amp[931] - std::complex<double> (0, 1) * amp[930] -
      std::complex<double> (0, 1) * amp[933] - std::complex<double> (0, 1) *
      amp[959] - std::complex<double> (0, 1) * amp[964] + std::complex<double>
      (0, 1) * amp[970] + std::complex<double> (0, 1) * amp[975] -
      std::complex<double> (0, 1) * amp[960] + std::complex<double> (0, 1) *
      amp[976] + std::complex<double> (0, 1) * amp[983] - std::complex<double>
      (0, 1) * amp[963] - std::complex<double> (0, 1) * amp[968] +
      std::complex<double> (0, 1) * amp[984];
  jamp[37] = +std::complex<double> (0, 1) * amp[0] - std::complex<double> (0,
      1) * amp[2] - std::complex<double> (0, 1) * amp[4] + std::complex<double>
      (0, 1) * amp[8] - std::complex<double> (0, 1) * amp[10] +
      std::complex<double> (0, 1) * amp[11] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[17] - std::complex<double>
      (0, 1) * amp[18] - std::complex<double> (0, 1) * amp[19] +
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[22] - std::complex<double> (0, 1) * amp[24] + std::complex<double>
      (0, 1) * amp[25] + std::complex<double> (0, 1) * amp[27] +
      std::complex<double> (0, 1) * amp[28] + std::complex<double> (0, 1) *
      amp[29] + std::complex<double> (0, 1) * amp[30] - std::complex<double>
      (0, 1) * amp[32] + amp[33] + amp[34] - std::complex<double> (0, 1) *
      amp[35] - std::complex<double> (0, 1) * amp[36] - std::complex<double>
      (0, 1) * amp[37] + amp[39] + amp[40] + amp[41] - std::complex<double> (0,
      1) * amp[42] + amp[43] - amp[98] - amp[102] - std::complex<double> (0, 1)
      * amp[103] - amp[107] - amp[108] - amp[109] - std::complex<double> (0, 1)
      * amp[110] - amp[112] + amp[257] + amp[259] + amp[261] + amp[262] +
      amp[263] + amp[267] - std::complex<double> (0, 1) * amp[270] -
      std::complex<double> (0, 1) * amp[271] - std::complex<double> (0, 1) *
      amp[272] - std::complex<double> (0, 1) * amp[273] - std::complex<double>
      (0, 1) * amp[274] - std::complex<double> (0, 1) * amp[275] - amp[276] -
      amp[277] - amp[278] - amp[282] - amp[284] + std::complex<double> (0, 1) *
      amp[298] - amp[299] - amp[300] + std::complex<double> (0, 1) * amp[301] +
      std::complex<double> (0, 1) * amp[302] + std::complex<double> (0, 1) *
      amp[303] - amp[305] - amp[306] - amp[307] + std::complex<double> (0, 1) *
      amp[308] - amp[309] - amp[313] - amp[312] - amp[316] - amp[315] -
      amp[319] - amp[318] - amp[323] - amp[322] + std::complex<double> (0, 1) *
      amp[325] - std::complex<double> (0, 1) * amp[678] - std::complex<double>
      (0, 1) * amp[679] - std::complex<double> (0, 1) * amp[681] -
      std::complex<double> (0, 1) * amp[682] - std::complex<double> (0, 1) *
      amp[686] - std::complex<double> (0, 1) * amp[687] - std::complex<double>
      (0, 1) * amp[689] - std::complex<double> (0, 1) * amp[690] -
      std::complex<double> (0, 1) * amp[692] - std::complex<double> (0, 1) *
      amp[693] - std::complex<double> (0, 1) * amp[694] - std::complex<double>
      (0, 1) * amp[695] - std::complex<double> (0, 1) * amp[696] -
      std::complex<double> (0, 1) * amp[697] - std::complex<double> (0, 1) *
      amp[698] - std::complex<double> (0, 1) * amp[699] - std::complex<double>
      (0, 1) * amp[701] - std::complex<double> (0, 1) * amp[702] + amp[704] +
      amp[706] + amp[708] + amp[709] + amp[710] + amp[714] +
      std::complex<double> (0, 1) * amp[729] + std::complex<double> (0, 1) *
      amp[733] + std::complex<double> (0, 1) * amp[734] + std::complex<double>
      (0, 1) * amp[736] + std::complex<double> (0, 1) * amp[738] +
      std::complex<double> (0, 1) * amp[739] + std::complex<double> (0, 1) *
      amp[741] + std::complex<double> (0, 1) * amp[742] + std::complex<double>
      (0, 1) * amp[744] + std::complex<double> (0, 1) * amp[745] +
      std::complex<double> (0, 1) * amp[746] + std::complex<double> (0, 1) *
      amp[750] + std::complex<double> (0, 1) * amp[751] + std::complex<double>
      (0, 1) * amp[753] - std::complex<double> (0, 1) * amp[780] - amp[782] -
      std::complex<double> (0, 1) * amp[827] + std::complex<double> (0, 1) *
      amp[829] - std::complex<double> (0, 1) * amp[833] + std::complex<double>
      (0, 1) * amp[835] - std::complex<double> (0, 1) * amp[836] +
      std::complex<double> (0, 1) * amp[838] - std::complex<double> (0, 1) *
      amp[839] + std::complex<double> (0, 1) * amp[841] - amp[844] + amp[842] -
      amp[847] + amp[845] - amp[850] + amp[848] - amp[853] + amp[851] +
      std::complex<double> (0, 1) * amp[870] - std::complex<double> (0, 1) *
      amp[868] - std::complex<double> (0, 1) * amp[872] + std::complex<double>
      (0, 1) * amp[873] - std::complex<double> (0, 1) * amp[875] -
      std::complex<double> (0, 1) * amp[876] - std::complex<double> (0, 1) *
      amp[939] - std::complex<double> (0, 1) * amp[938] - std::complex<double>
      (0, 1) * amp[945] - std::complex<double> (0, 1) * amp[944] -
      std::complex<double> (0, 1) * amp[948] - std::complex<double> (0, 1) *
      amp[947] - std::complex<double> (0, 1) * amp[951] - std::complex<double>
      (0, 1) * amp[950] + std::complex<double> (0, 1) * amp[964] +
      std::complex<double> (0, 1) * amp[985] - std::complex<double> (0, 1) *
      amp[970] + std::complex<double> (0, 1) * amp[981] - std::complex<double>
      (0, 1) * amp[956] + std::complex<double> (0, 1) * amp[982] -
      std::complex<double> (0, 1) * amp[962] - std::complex<double> (0, 1) *
      amp[957] - std::complex<double> (0, 1) * amp[983] + std::complex<double>
      (0, 1) * amp[968];
  jamp[38] = +2. * (-amp[31] - amp[44] - amp[82] + std::complex<double> (0, 1)
      * amp[95] + std::complex<double> (0, 1) * amp[96]);
  jamp[39] = -std::complex<double> (0, 1) * amp[1] - std::complex<double> (0,
      1) * amp[2] - std::complex<double> (0, 1) * amp[3] - std::complex<double>
      (0, 1) * amp[4] - std::complex<double> (0, 1) * amp[9] -
      std::complex<double> (0, 1) * amp[10] - std::complex<double> (0, 1) *
      amp[12] - std::complex<double> (0, 1) * amp[13] - std::complex<double>
      (0, 1) * amp[14] - std::complex<double> (0, 1) * amp[15] -
      std::complex<double> (0, 1) * amp[16] - std::complex<double> (0, 1) *
      amp[17] - std::complex<double> (0, 1) * amp[18] - std::complex<double>
      (0, 1) * amp[19] - std::complex<double> (0, 1) * amp[21] -
      std::complex<double> (0, 1) * amp[22] - std::complex<double> (0, 1) *
      amp[23] - std::complex<double> (0, 1) * amp[24] + amp[52] + amp[53] +
      amp[58] + amp[59] + amp[60] + amp[62] + std::complex<double> (0, 1) *
      amp[64] + std::complex<double> (0, 1) * amp[65] + std::complex<double>
      (0, 1) * amp[66] + std::complex<double> (0, 1) * amp[67] +
      std::complex<double> (0, 1) * amp[68] - amp[71] - amp[72] - amp[77] -
      amp[78] - amp[79] - amp[81] - std::complex<double> (0, 1) * amp[110] -
      std::complex<double> (0, 1) * amp[111] - amp[480] - amp[484] +
      std::complex<double> (0, 1) * amp[486] - amp[487] - amp[488] - amp[489] +
      std::complex<double> (0, 1) * amp[493] + std::complex<double> (0, 1) *
      amp[494] + std::complex<double> (0, 1) * amp[496] + std::complex<double>
      (0, 1) * amp[497] + std::complex<double> (0, 1) * amp[498] -
      std::complex<double> (0, 1) * amp[499] - std::complex<double> (0, 1) *
      amp[500] - std::complex<double> (0, 1) * amp[502] - std::complex<double>
      (0, 1) * amp[507] - std::complex<double> (0, 1) * amp[508] -
      std::complex<double> (0, 1) * amp[510] - std::complex<double> (0, 1) *
      amp[511] - std::complex<double> (0, 1) * amp[513] - std::complex<double>
      (0, 1) * amp[514] - std::complex<double> (0, 1) * amp[515] -
      std::complex<double> (0, 1) * amp[519] - std::complex<double> (0, 1) *
      amp[520] - std::complex<double> (0, 1) * amp[522] - amp[538] -
      std::complex<double> (0, 1) * amp[539] - amp[540] - amp[542] - amp[543] -
      amp[544] - std::complex<double> (0, 1) * amp[545] - std::complex<double>
      (0, 1) * amp[546] - std::complex<double> (0, 1) * amp[547] - amp[548] -
      std::complex<double> (0, 1) * amp[549] - amp[552] + std::complex<double>
      (0, 1) * amp[553] - std::complex<double> (0, 1) * amp[626] +
      std::complex<double> (0, 1) * amp[628] - std::complex<double> (0, 1) *
      amp[629] - std::complex<double> (0, 1) * amp[634] + std::complex<double>
      (0, 1) * amp[636] - std::complex<double> (0, 1) * amp[637] +
      std::complex<double> (0, 1) * amp[639] - std::complex<double> (0, 1) *
      amp[640] - std::complex<double> (0, 1) * amp[641] - std::complex<double>
      (0, 1) * amp[642] - std::complex<double> (0, 1) * amp[646] +
      std::complex<double> (0, 1) * amp[648] - std::complex<double> (0, 1) *
      amp[649] + amp[652] + std::complex<double> (0, 1) * amp[653] + amp[654] +
      amp[656] + amp[657] + amp[658] + std::complex<double> (0, 1) * amp[659] +
      std::complex<double> (0, 1) * amp[660] + std::complex<double> (0, 1) *
      amp[661] + amp[662] + std::complex<double> (0, 1) * amp[663] +
      std::complex<double> (0, 1) * amp[677] + std::complex<double> (0, 1) *
      amp[766] - amp[767] - amp[768] - amp[769] - amp[773] - amp[775] -
      std::complex<double> (0, 1) * amp[783] + std::complex<double> (0, 1) *
      amp[785] - std::complex<double> (0, 1) * amp[789] + std::complex<double>
      (0, 1) * amp[791] - std::complex<double> (0, 1) * amp[792] +
      std::complex<double> (0, 1) * amp[794] - std::complex<double> (0, 1) *
      amp[795] + std::complex<double> (0, 1) * amp[797] - amp[800] + amp[798] -
      amp[803] + amp[801] - amp[806] + amp[804] - amp[809] + amp[807] +
      std::complex<double> (0, 1) * amp[826] - std::complex<double> (0, 1) *
      amp[824] + std::complex<double> (0, 1) * amp[828] + std::complex<double>
      (0, 1) * amp[829] + std::complex<double> (0, 1) * amp[834] +
      std::complex<double> (0, 1) * amp[835] + std::complex<double> (0, 1) *
      amp[837] + std::complex<double> (0, 1) * amp[838] + std::complex<double>
      (0, 1) * amp[840] + std::complex<double> (0, 1) * amp[841] + amp[857] +
      amp[856] + amp[860] + amp[859] + amp[863] + amp[862] + amp[866] +
      amp[865] + std::complex<double> (0, 1) * amp[870] + std::complex<double>
      (0, 1) * amp[869] - amp[881] - std::complex<double> (0, 1) * amp[882] +
      std::complex<double> (0, 1) * amp[936] - std::complex<double> (0, 1) *
      amp[974] + std::complex<double> (0, 1) * amp[964] - std::complex<double>
      (0, 1) * amp[970] - std::complex<double> (0, 1) * amp[975] -
      std::complex<double> (0, 1) * amp[965] + std::complex<double> (0, 1) *
      amp[977] + std::complex<double> (0, 1) * amp[967] + std::complex<double>
      (0, 1) * amp[978] + std::complex<double> (0, 1) * amp[963] -
      std::complex<double> (0, 1) * amp[969];
  jamp[40] = +2. * (-amp[237] - amp[250] - std::complex<double> (0, 1) *
      amp[482] - std::complex<double> (0, 1) * amp[483] - amp[651]);
  jamp[41] = -amp[347] - amp[348] - amp[353] - amp[354] - amp[355] - amp[357] +
      std::complex<double> (0, 1) * amp[359] + std::complex<double> (0, 1) *
      amp[360] + std::complex<double> (0, 1) * amp[361] + std::complex<double>
      (0, 1) * amp[362] + std::complex<double> (0, 1) * amp[363] -
      std::complex<double> (0, 1) * amp[385] + std::complex<double> (0, 1) *
      amp[387] - std::complex<double> (0, 1) * amp[388] + std::complex<double>
      (0, 1) * amp[389] - std::complex<double> (0, 1) * amp[393] +
      std::complex<double> (0, 1) * amp[395] - std::complex<double> (0, 1) *
      amp[396] + std::complex<double> (0, 1) * amp[398] - std::complex<double>
      (0, 1) * amp[399] - std::complex<double> (0, 1) * amp[400] -
      std::complex<double> (0, 1) * amp[401] + std::complex<double> (0, 1) *
      amp[402] + std::complex<double> (0, 1) * amp[403] + std::complex<double>
      (0, 1) * amp[404] - std::complex<double> (0, 1) * amp[405] +
      std::complex<double> (0, 1) * amp[407] - std::complex<double> (0, 1) *
      amp[408] + std::complex<double> (0, 1) * amp[409] + amp[411] + amp[413] +
      amp[415] + amp[416] + amp[417] + amp[421] - std::complex<double> (0, 1) *
      amp[436] + std::complex<double> (0, 1) * amp[439] - std::complex<double>
      (0, 1) * amp[460] - amp[461] - amp[462] - std::complex<double> (0, 1) *
      amp[463] - std::complex<double> (0, 1) * amp[464] - std::complex<double>
      (0, 1) * amp[465] - amp[467] - amp[468] - amp[469] - std::complex<double>
      (0, 1) * amp[470] - amp[471] + amp[480] + amp[484] - std::complex<double>
      (0, 1) * amp[486] + amp[487] + amp[488] + amp[489] - std::complex<double>
      (0, 1) * amp[493] - std::complex<double> (0, 1) * amp[494] -
      std::complex<double> (0, 1) * amp[496] - std::complex<double> (0, 1) *
      amp[497] - std::complex<double> (0, 1) * amp[498] + std::complex<double>
      (0, 1) * amp[500] + std::complex<double> (0, 1) * amp[501] +
      std::complex<double> (0, 1) * amp[503] + std::complex<double> (0, 1) *
      amp[508] + std::complex<double> (0, 1) * amp[509] + std::complex<double>
      (0, 1) * amp[511] + std::complex<double> (0, 1) * amp[512] +
      std::complex<double> (0, 1) * amp[516] + std::complex<double> (0, 1) *
      amp[517] + std::complex<double> (0, 1) * amp[518] + std::complex<double>
      (0, 1) * amp[520] + std::complex<double> (0, 1) * amp[521] +
      std::complex<double> (0, 1) * amp[523] - std::complex<double> (0, 1) *
      amp[550] + amp[552] + std::complex<double> (0, 1) * amp[601] - amp[605] -
      amp[606] - amp[607] - amp[609] - amp[611] + std::complex<double> (0, 1) *
      amp[626] - std::complex<double> (0, 1) * amp[628] + std::complex<double>
      (0, 1) * amp[629] + std::complex<double> (0, 1) * amp[634] -
      std::complex<double> (0, 1) * amp[636] + std::complex<double> (0, 1) *
      amp[637] - std::complex<double> (0, 1) * amp[639] + std::complex<double>
      (0, 1) * amp[640] + std::complex<double> (0, 1) * amp[641] +
      std::complex<double> (0, 1) * amp[642] + std::complex<double> (0, 1) *
      amp[646] - std::complex<double> (0, 1) * amp[648] + std::complex<double>
      (0, 1) * amp[649] - amp[652] - std::complex<double> (0, 1) * amp[653] -
      amp[654] - amp[656] - amp[657] - amp[658] - std::complex<double> (0, 1) *
      amp[659] - std::complex<double> (0, 1) * amp[660] - std::complex<double>
      (0, 1) * amp[661] - amp[662] - std::complex<double> (0, 1) * amp[663] -
      std::complex<double> (0, 1) * amp[677] + std::complex<double> (0, 1) *
      amp[784] + std::complex<double> (0, 1) * amp[783] + std::complex<double>
      (0, 1) * amp[790] + std::complex<double> (0, 1) * amp[789] +
      std::complex<double> (0, 1) * amp[793] + std::complex<double> (0, 1) *
      amp[792] + std::complex<double> (0, 1) * amp[796] + std::complex<double>
      (0, 1) * amp[795] - amp[799] - amp[798] - amp[802] - amp[801] - amp[805]
      - amp[804] - amp[808] - amp[807] + std::complex<double> (0, 1) * amp[825]
      + std::complex<double> (0, 1) * amp[824] - std::complex<double> (0, 1) *
      amp[880] + amp[881] + amp[904] + amp[903] + amp[907] + amp[906] +
      amp[910] + amp[909] + amp[913] + amp[912] - std::complex<double> (0, 1) *
      amp[916] - std::complex<double> (0, 1) * amp[917] - std::complex<double>
      (0, 1) * amp[922] - std::complex<double> (0, 1) * amp[923] -
      std::complex<double> (0, 1) * amp[925] - std::complex<double> (0, 1) *
      amp[926] - std::complex<double> (0, 1) * amp[928] - std::complex<double>
      (0, 1) * amp[929] + std::complex<double> (0, 1) * amp[932] +
      std::complex<double> (0, 1) * amp[931] - std::complex<double> (0, 1) *
      amp[936] + std::complex<double> (0, 1) * amp[974] + std::complex<double>
      (0, 1) * amp[960] + std::complex<double> (0, 1) * amp[965] -
      std::complex<double> (0, 1) * amp[976] + std::complex<double> (0, 1) *
      amp[961] + std::complex<double> (0, 1) * amp[966] - std::complex<double>
      (0, 1) * amp[977] - std::complex<double> (0, 1) * amp[972] -
      std::complex<double> (0, 1) * amp[978] - std::complex<double> (0, 1) *
      amp[973];
  jamp[42] = +std::complex<double> (0, 1) * amp[1] + std::complex<double> (0,
      1) * amp[2] + std::complex<double> (0, 1) * amp[3] + std::complex<double>
      (0, 1) * amp[4] + std::complex<double> (0, 1) * amp[9] +
      std::complex<double> (0, 1) * amp[10] + std::complex<double> (0, 1) *
      amp[12] + std::complex<double> (0, 1) * amp[13] + std::complex<double>
      (0, 1) * amp[14] + std::complex<double> (0, 1) * amp[15] +
      std::complex<double> (0, 1) * amp[16] + std::complex<double> (0, 1) *
      amp[17] + std::complex<double> (0, 1) * amp[18] + std::complex<double>
      (0, 1) * amp[19] + std::complex<double> (0, 1) * amp[21] +
      std::complex<double> (0, 1) * amp[22] + std::complex<double> (0, 1) *
      amp[23] + std::complex<double> (0, 1) * amp[24] - amp[52] - amp[53] -
      amp[58] - amp[59] - amp[60] - amp[62] - std::complex<double> (0, 1) *
      amp[64] - std::complex<double> (0, 1) * amp[65] - std::complex<double>
      (0, 1) * amp[66] - std::complex<double> (0, 1) * amp[67] -
      std::complex<double> (0, 1) * amp[68] + amp[71] + amp[72] + amp[77] +
      amp[78] + amp[79] + amp[81] + std::complex<double> (0, 1) * amp[110] +
      std::complex<double> (0, 1) * amp[111] - amp[367] - amp[371] -
      std::complex<double> (0, 1) * amp[372] - amp[376] - amp[377] - amp[378] -
      std::complex<double> (0, 1) * amp[379] - std::complex<double> (0, 1) *
      amp[380] - std::complex<double> (0, 1) * amp[382] - std::complex<double>
      (0, 1) * amp[383] - std::complex<double> (0, 1) * amp[384] +
      std::complex<double> (0, 1) * amp[385] + std::complex<double> (0, 1) *
      amp[386] + std::complex<double> (0, 1) * amp[388] + std::complex<double>
      (0, 1) * amp[393] + std::complex<double> (0, 1) * amp[394] +
      std::complex<double> (0, 1) * amp[396] + std::complex<double> (0, 1) *
      amp[397] + std::complex<double> (0, 1) * amp[399] + std::complex<double>
      (0, 1) * amp[400] + std::complex<double> (0, 1) * amp[401] +
      std::complex<double> (0, 1) * amp[405] + std::complex<double> (0, 1) *
      amp[406] + std::complex<double> (0, 1) * amp[408] - amp[411] +
      std::complex<double> (0, 1) * amp[412] - amp[413] - amp[415] - amp[416] -
      amp[417] + std::complex<double> (0, 1) * amp[418] + std::complex<double>
      (0, 1) * amp[419] + std::complex<double> (0, 1) * amp[420] - amp[421] +
      std::complex<double> (0, 1) * amp[422] - amp[437] - std::complex<double>
      (0, 1) * amp[439] + std::complex<double> (0, 1) * amp[678] -
      std::complex<double> (0, 1) * amp[680] + std::complex<double> (0, 1) *
      amp[681] + std::complex<double> (0, 1) * amp[686] - std::complex<double>
      (0, 1) * amp[688] + std::complex<double> (0, 1) * amp[689] -
      std::complex<double> (0, 1) * amp[691] + std::complex<double> (0, 1) *
      amp[692] + std::complex<double> (0, 1) * amp[693] + std::complex<double>
      (0, 1) * amp[694] + std::complex<double> (0, 1) * amp[698] -
      std::complex<double> (0, 1) * amp[700] + std::complex<double> (0, 1) *
      amp[701] + amp[717] - std::complex<double> (0, 1) * amp[718] + amp[719] +
      amp[721] + amp[722] + amp[723] - std::complex<double> (0, 1) * amp[724] -
      std::complex<double> (0, 1) * amp[725] - std::complex<double> (0, 1) *
      amp[726] + amp[727] - std::complex<double> (0, 1) * amp[728] -
      std::complex<double> (0, 1) * amp[729] - std::complex<double> (0, 1) *
      amp[766] - amp[770] - amp[771] - amp[772] - amp[774] - amp[776] -
      std::complex<double> (0, 1) * amp[784] - std::complex<double> (0, 1) *
      amp[785] - std::complex<double> (0, 1) * amp[790] - std::complex<double>
      (0, 1) * amp[791] - std::complex<double> (0, 1) * amp[793] -
      std::complex<double> (0, 1) * amp[794] - std::complex<double> (0, 1) *
      amp[796] - std::complex<double> (0, 1) * amp[797] + amp[800] + amp[799] +
      amp[803] + amp[802] + amp[806] + amp[805] + amp[809] + amp[808] -
      std::complex<double> (0, 1) * amp[826] - std::complex<double> (0, 1) *
      amp[825] + std::complex<double> (0, 1) * amp[827] - std::complex<double>
      (0, 1) * amp[829] + std::complex<double> (0, 1) * amp[833] -
      std::complex<double> (0, 1) * amp[835] + std::complex<double> (0, 1) *
      amp[836] - std::complex<double> (0, 1) * amp[838] + std::complex<double>
      (0, 1) * amp[839] - std::complex<double> (0, 1) * amp[841] - amp[857] +
      amp[855] - amp[860] + amp[858] - amp[863] + amp[861] - amp[866] +
      amp[864] - std::complex<double> (0, 1) * amp[870] + std::complex<double>
      (0, 1) * amp[868] - amp[884] + std::complex<double> (0, 1) * amp[885] -
      std::complex<double> (0, 1) * amp[935] - std::complex<double> (0, 1) *
      amp[964] + std::complex<double> (0, 1) * amp[970] - std::complex<double>
      (0, 1) * amp[980] + std::complex<double> (0, 1) * amp[971] -
      std::complex<double> (0, 1) * amp[981] - std::complex<double> (0, 1) *
      amp[966] + std::complex<double> (0, 1) * amp[962] - std::complex<double>
      (0, 1) * amp[967] + std::complex<double> (0, 1) * amp[972] +
      std::complex<double> (0, 1) * amp[969];
  jamp[43] = +2. * (-amp[345] - amp[358] - std::complex<double> (0, 1) *
      amp[364] - std::complex<double> (0, 1) * amp[365] - amp[410]);
  jamp[44] = -amp[238] - amp[240] - amp[242] - amp[243] - amp[244] - amp[248] +
      std::complex<double> (0, 1) * amp[251] + std::complex<double> (0, 1) *
      amp[252] + std::complex<double> (0, 1) * amp[253] + std::complex<double>
      (0, 1) * amp[254] + std::complex<double> (0, 1) * amp[255] - amp[257] -
      std::complex<double> (0, 1) * amp[258] - amp[259] - amp[261] - amp[262] -
      amp[263] - std::complex<double> (0, 1) * amp[264] - std::complex<double>
      (0, 1) * amp[265] - std::complex<double> (0, 1) * amp[266] - amp[267] -
      std::complex<double> (0, 1) * amp[268] + std::complex<double> (0, 1) *
      amp[286] - amp[290] - amp[291] - amp[292] - amp[294] - amp[296] -
      std::complex<double> (0, 1) * amp[311] + amp[314] + amp[313] + amp[317] +
      amp[316] + amp[320] + amp[319] + amp[324] + amp[323] + amp[367] +
      amp[371] + std::complex<double> (0, 1) * amp[372] + amp[376] + amp[377] +
      amp[378] + std::complex<double> (0, 1) * amp[379] + std::complex<double>
      (0, 1) * amp[380] + std::complex<double> (0, 1) * amp[382] +
      std::complex<double> (0, 1) * amp[383] + std::complex<double> (0, 1) *
      amp[384] - std::complex<double> (0, 1) * amp[385] - std::complex<double>
      (0, 1) * amp[386] - std::complex<double> (0, 1) * amp[388] -
      std::complex<double> (0, 1) * amp[393] - std::complex<double> (0, 1) *
      amp[394] - std::complex<double> (0, 1) * amp[396] - std::complex<double>
      (0, 1) * amp[397] - std::complex<double> (0, 1) * amp[399] -
      std::complex<double> (0, 1) * amp[400] - std::complex<double> (0, 1) *
      amp[401] - std::complex<double> (0, 1) * amp[405] - std::complex<double>
      (0, 1) * amp[406] - std::complex<double> (0, 1) * amp[408] + amp[411] -
      std::complex<double> (0, 1) * amp[412] + amp[413] + amp[415] + amp[416] +
      amp[417] - std::complex<double> (0, 1) * amp[418] - std::complex<double>
      (0, 1) * amp[419] - std::complex<double> (0, 1) * amp[420] + amp[421] -
      std::complex<double> (0, 1) * amp[422] + amp[437] + std::complex<double>
      (0, 1) * amp[439] + std::complex<double> (0, 1) * amp[626] +
      std::complex<double> (0, 1) * amp[627] + std::complex<double> (0, 1) *
      amp[629] + std::complex<double> (0, 1) * amp[630] + std::complex<double>
      (0, 1) * amp[634] + std::complex<double> (0, 1) * amp[635] +
      std::complex<double> (0, 1) * amp[637] + std::complex<double> (0, 1) *
      amp[638] + std::complex<double> (0, 1) * amp[640] + std::complex<double>
      (0, 1) * amp[641] + std::complex<double> (0, 1) * amp[642] +
      std::complex<double> (0, 1) * amp[643] + std::complex<double> (0, 1) *
      amp[644] + std::complex<double> (0, 1) * amp[645] + std::complex<double>
      (0, 1) * amp[646] + std::complex<double> (0, 1) * amp[647] +
      std::complex<double> (0, 1) * amp[649] + std::complex<double> (0, 1) *
      amp[650] - amp[652] - amp[654] - amp[656] - amp[657] - amp[658] -
      amp[662] - std::complex<double> (0, 1) * amp[677] + std::complex<double>
      (0, 1) * amp[679] + std::complex<double> (0, 1) * amp[680] +
      std::complex<double> (0, 1) * amp[682] + std::complex<double> (0, 1) *
      amp[687] + std::complex<double> (0, 1) * amp[688] + std::complex<double>
      (0, 1) * amp[690] + std::complex<double> (0, 1) * amp[691] +
      std::complex<double> (0, 1) * amp[695] + std::complex<double> (0, 1) *
      amp[696] + std::complex<double> (0, 1) * amp[697] + std::complex<double>
      (0, 1) * amp[699] + std::complex<double> (0, 1) * amp[700] +
      std::complex<double> (0, 1) * amp[702] + std::complex<double> (0, 1) *
      amp[784] + std::complex<double> (0, 1) * amp[783] + std::complex<double>
      (0, 1) * amp[790] + std::complex<double> (0, 1) * amp[789] +
      std::complex<double> (0, 1) * amp[793] + std::complex<double> (0, 1) *
      amp[792] + std::complex<double> (0, 1) * amp[796] + std::complex<double>
      (0, 1) * amp[795] - amp[799] - amp[798] - amp[802] - amp[801] - amp[805]
      - amp[804] - amp[808] - amp[807] + std::complex<double> (0, 1) * amp[825]
      + std::complex<double> (0, 1) * amp[824] + std::complex<double> (0, 1) *
      amp[871] + std::complex<double> (0, 1) * amp[872] + std::complex<double>
      (0, 1) * amp[875] + std::complex<double> (0, 1) * amp[874] + amp[884] -
      std::complex<double> (0, 1) * amp[885] + std::complex<double> (0, 1) *
      amp[939] + std::complex<double> (0, 1) * amp[940] + std::complex<double>
      (0, 1) * amp[945] + std::complex<double> (0, 1) * amp[946] +
      std::complex<double> (0, 1) * amp[948] + std::complex<double> (0, 1) *
      amp[949] + std::complex<double> (0, 1) * amp[951] + std::complex<double>
      (0, 1) * amp[952] + std::complex<double> (0, 1) * amp[980] +
      std::complex<double> (0, 1) * amp[965] - std::complex<double> (0, 1) *
      amp[971] + std::complex<double> (0, 1) * amp[966] - std::complex<double>
      (0, 1) * amp[982] - std::complex<double> (0, 1) * amp[972] +
      std::complex<double> (0, 1) * amp[957] - std::complex<double> (0, 1) *
      amp[978] + std::complex<double> (0, 1) * amp[958] - std::complex<double>
      (0, 1) * amp[979];
  jamp[45] = +2. * (+std::complex<double> (0, 1) * amp[57] +
      std::complex<double> (0, 1) * amp[89] - std::complex<double> (0, 1) *
      amp[414] - std::complex<double> (0, 1) * amp[427] + std::complex<double>
      (0, 1) * amp[655] + std::complex<double> (0, 1) * amp[668] +
      std::complex<double> (0, 1) * amp[810] + std::complex<double> (0, 1) *
      amp[823]);
  jamp[46] = +2. * (-std::complex<double> (0, 1) * amp[57] -
      std::complex<double> (0, 1) * amp[89] + std::complex<double> (0, 1) *
      amp[414] + std::complex<double> (0, 1) * amp[427] - std::complex<double>
      (0, 1) * amp[655] - std::complex<double> (0, 1) * amp[668] -
      std::complex<double> (0, 1) * amp[810] - std::complex<double> (0, 1) *
      amp[823]);
  jamp[47] = +std::complex<double> (0, 1) * amp[0] + std::complex<double> (0,
      1) * amp[1] + std::complex<double> (0, 1) * amp[3] + std::complex<double>
      (0, 1) * amp[8] + std::complex<double> (0, 1) * amp[9] +
      std::complex<double> (0, 1) * amp[11] + std::complex<double> (0, 1) *
      amp[12] + std::complex<double> (0, 1) * amp[14] + std::complex<double>
      (0, 1) * amp[15] + std::complex<double> (0, 1) * amp[16] +
      std::complex<double> (0, 1) * amp[20] + std::complex<double> (0, 1) *
      amp[21] + std::complex<double> (0, 1) * amp[23] + std::complex<double>
      (0, 1) * amp[25] + std::complex<double> (0, 1) * amp[27] +
      std::complex<double> (0, 1) * amp[28] + std::complex<double> (0, 1) *
      amp[29] + std::complex<double> (0, 1) * amp[30] + std::complex<double>
      (0, 1) * amp[51] - amp[52] - amp[53] + std::complex<double> (0, 1) *
      amp[54] + std::complex<double> (0, 1) * amp[55] + std::complex<double>
      (0, 1) * amp[56] - amp[58] - amp[59] - amp[60] + std::complex<double> (0,
      1) * amp[61] - amp[62] - amp[98] - amp[102] - std::complex<double> (0, 1)
      * amp[103] - amp[107] - amp[108] - amp[109] + std::complex<double> (0, 1)
      * amp[111] - amp[112] + amp[347] + amp[348] + amp[353] + amp[354] +
      amp[355] + amp[357] - std::complex<double> (0, 1) * amp[359] -
      std::complex<double> (0, 1) * amp[360] - std::complex<double> (0, 1) *
      amp[361] - std::complex<double> (0, 1) * amp[362] - std::complex<double>
      (0, 1) * amp[363] + std::complex<double> (0, 1) * amp[385] -
      std::complex<double> (0, 1) * amp[387] + std::complex<double> (0, 1) *
      amp[388] - std::complex<double> (0, 1) * amp[389] + std::complex<double>
      (0, 1) * amp[393] - std::complex<double> (0, 1) * amp[395] +
      std::complex<double> (0, 1) * amp[396] - std::complex<double> (0, 1) *
      amp[398] + std::complex<double> (0, 1) * amp[399] + std::complex<double>
      (0, 1) * amp[400] + std::complex<double> (0, 1) * amp[401] -
      std::complex<double> (0, 1) * amp[402] - std::complex<double> (0, 1) *
      amp[403] - std::complex<double> (0, 1) * amp[404] + std::complex<double>
      (0, 1) * amp[405] - std::complex<double> (0, 1) * amp[407] +
      std::complex<double> (0, 1) * amp[408] - std::complex<double> (0, 1) *
      amp[409] - amp[411] - amp[413] - amp[415] - amp[416] - amp[417] -
      amp[421] + std::complex<double> (0, 1) * amp[436] - std::complex<double>
      (0, 1) * amp[439] - std::complex<double> (0, 1) * amp[601] - amp[602] -
      amp[603] - amp[604] - amp[608] - amp[610] - std::complex<double> (0, 1) *
      amp[613] + amp[614] + amp[615] - std::complex<double> (0, 1) * amp[616] -
      std::complex<double> (0, 1) * amp[617] - std::complex<double> (0, 1) *
      amp[618] + amp[620] + amp[621] + amp[622] - std::complex<double> (0, 1) *
      amp[623] + amp[624] - std::complex<double> (0, 1) * amp[625] +
      std::complex<double> (0, 1) * amp[733] - std::complex<double> (0, 1) *
      amp[735] - std::complex<double> (0, 1) * amp[737] + std::complex<double>
      (0, 1) * amp[738] - std::complex<double> (0, 1) * amp[740] +
      std::complex<double> (0, 1) * amp[741] - std::complex<double> (0, 1) *
      amp[743] - std::complex<double> (0, 1) * amp[747] - std::complex<double>
      (0, 1) * amp[748] - std::complex<double> (0, 1) * amp[749] +
      std::complex<double> (0, 1) * amp[750] - std::complex<double> (0, 1) *
      amp[752] - std::complex<double> (0, 1) * amp[754] + std::complex<double>
      (0, 1) * amp[781] - amp[782] - std::complex<double> (0, 1) * amp[784] -
      std::complex<double> (0, 1) * amp[785] - std::complex<double> (0, 1) *
      amp[790] - std::complex<double> (0, 1) * amp[791] - std::complex<double>
      (0, 1) * amp[793] - std::complex<double> (0, 1) * amp[794] -
      std::complex<double> (0, 1) * amp[796] - std::complex<double> (0, 1) *
      amp[797] + amp[800] + amp[799] + amp[803] + amp[802] + amp[806] +
      amp[805] + amp[809] + amp[808] - std::complex<double> (0, 1) * amp[826] -
      std::complex<double> (0, 1) * amp[825] - amp[904] + amp[902] - amp[907] +
      amp[905] - amp[910] + amp[908] - amp[913] + amp[911] -
      std::complex<double> (0, 1) * amp[915] + std::complex<double> (0, 1) *
      amp[917] - std::complex<double> (0, 1) * amp[921] + std::complex<double>
      (0, 1) * amp[923] - std::complex<double> (0, 1) * amp[924] +
      std::complex<double> (0, 1) * amp[926] - std::complex<double> (0, 1) *
      amp[927] + std::complex<double> (0, 1) * amp[929] - std::complex<double>
      (0, 1) * amp[932] + std::complex<double> (0, 1) * amp[930] +
      std::complex<double> (0, 1) * amp[933] + std::complex<double> (0, 1) *
      amp[959] - std::complex<double> (0, 1) * amp[961] - std::complex<double>
      (0, 1) * amp[966] - std::complex<double> (0, 1) * amp[967] +
      std::complex<double> (0, 1) * amp[972] - std::complex<double> (0, 1) *
      amp[983] + std::complex<double> (0, 1) * amp[968] + std::complex<double>
      (0, 1) * amp[973] - std::complex<double> (0, 1) * amp[984] +
      std::complex<double> (0, 1) * amp[969];
  jamp[48] = -std::complex<double> (0, 1) * amp[0] - std::complex<double> (0,
      1) * amp[1] - std::complex<double> (0, 1) * amp[3] - std::complex<double>
      (0, 1) * amp[8] - std::complex<double> (0, 1) * amp[9] -
      std::complex<double> (0, 1) * amp[11] - std::complex<double> (0, 1) *
      amp[12] - std::complex<double> (0, 1) * amp[14] - std::complex<double>
      (0, 1) * amp[15] - std::complex<double> (0, 1) * amp[16] -
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[21] - std::complex<double> (0, 1) * amp[23] - std::complex<double>
      (0, 1) * amp[25] - std::complex<double> (0, 1) * amp[27] -
      std::complex<double> (0, 1) * amp[28] - std::complex<double> (0, 1) *
      amp[29] - std::complex<double> (0, 1) * amp[30] - std::complex<double>
      (0, 1) * amp[51] + amp[52] + amp[53] - std::complex<double> (0, 1) *
      amp[54] - std::complex<double> (0, 1) * amp[55] - std::complex<double>
      (0, 1) * amp[56] + amp[58] + amp[59] + amp[60] - std::complex<double> (0,
      1) * amp[61] + amp[62] + amp[98] + amp[102] + std::complex<double> (0, 1)
      * amp[103] + amp[107] + amp[108] + amp[109] - std::complex<double> (0, 1)
      * amp[111] + amp[112] + amp[238] + amp[240] + amp[242] + amp[243] +
      amp[244] + amp[248] - std::complex<double> (0, 1) * amp[251] -
      std::complex<double> (0, 1) * amp[252] - std::complex<double> (0, 1) *
      amp[253] - std::complex<double> (0, 1) * amp[254] - std::complex<double>
      (0, 1) * amp[255] - std::complex<double> (0, 1) * amp[286] - amp[287] -
      amp[288] - amp[289] - amp[293] - amp[295] - std::complex<double> (0, 1) *
      amp[298] + amp[299] + amp[300] - std::complex<double> (0, 1) * amp[301] -
      std::complex<double> (0, 1) * amp[302] - std::complex<double> (0, 1) *
      amp[303] + amp[305] + amp[306] + amp[307] - std::complex<double> (0, 1) *
      amp[308] + amp[309] - amp[314] + amp[312] - amp[317] + amp[315] -
      amp[320] + amp[318] - amp[324] + amp[322] - std::complex<double> (0, 1) *
      amp[325] - std::complex<double> (0, 1) * amp[626] - std::complex<double>
      (0, 1) * amp[627] - std::complex<double> (0, 1) * amp[629] -
      std::complex<double> (0, 1) * amp[630] - std::complex<double> (0, 1) *
      amp[634] - std::complex<double> (0, 1) * amp[635] - std::complex<double>
      (0, 1) * amp[637] - std::complex<double> (0, 1) * amp[638] -
      std::complex<double> (0, 1) * amp[640] - std::complex<double> (0, 1) *
      amp[641] - std::complex<double> (0, 1) * amp[642] - std::complex<double>
      (0, 1) * amp[643] - std::complex<double> (0, 1) * amp[644] -
      std::complex<double> (0, 1) * amp[645] - std::complex<double> (0, 1) *
      amp[646] - std::complex<double> (0, 1) * amp[647] - std::complex<double>
      (0, 1) * amp[649] - std::complex<double> (0, 1) * amp[650] + amp[652] +
      amp[654] + amp[656] + amp[657] + amp[658] + amp[662] +
      std::complex<double> (0, 1) * amp[677] - std::complex<double> (0, 1) *
      amp[733] - std::complex<double> (0, 1) * amp[734] - std::complex<double>
      (0, 1) * amp[736] - std::complex<double> (0, 1) * amp[738] -
      std::complex<double> (0, 1) * amp[739] - std::complex<double> (0, 1) *
      amp[741] - std::complex<double> (0, 1) * amp[742] - std::complex<double>
      (0, 1) * amp[744] - std::complex<double> (0, 1) * amp[745] -
      std::complex<double> (0, 1) * amp[746] - std::complex<double> (0, 1) *
      amp[750] - std::complex<double> (0, 1) * amp[751] - std::complex<double>
      (0, 1) * amp[753] - std::complex<double> (0, 1) * amp[781] + amp[782] -
      std::complex<double> (0, 1) * amp[783] + std::complex<double> (0, 1) *
      amp[785] - std::complex<double> (0, 1) * amp[789] + std::complex<double>
      (0, 1) * amp[791] - std::complex<double> (0, 1) * amp[792] +
      std::complex<double> (0, 1) * amp[794] - std::complex<double> (0, 1) *
      amp[795] + std::complex<double> (0, 1) * amp[797] - amp[800] + amp[798] -
      amp[803] + amp[801] - amp[806] + amp[804] - amp[809] + amp[807] +
      std::complex<double> (0, 1) * amp[826] - std::complex<double> (0, 1) *
      amp[824] - std::complex<double> (0, 1) * amp[871] - std::complex<double>
      (0, 1) * amp[873] - std::complex<double> (0, 1) * amp[874] +
      std::complex<double> (0, 1) * amp[876] + std::complex<double> (0, 1) *
      amp[938] - std::complex<double> (0, 1) * amp[940] + std::complex<double>
      (0, 1) * amp[944] - std::complex<double> (0, 1) * amp[946] +
      std::complex<double> (0, 1) * amp[947] - std::complex<double> (0, 1) *
      amp[949] + std::complex<double> (0, 1) * amp[950] - std::complex<double>
      (0, 1) * amp[952] - std::complex<double> (0, 1) * amp[985] -
      std::complex<double> (0, 1) * amp[965] + std::complex<double> (0, 1) *
      amp[956] + std::complex<double> (0, 1) * amp[967] + std::complex<double>
      (0, 1) * amp[978] + std::complex<double> (0, 1) * amp[983] -
      std::complex<double> (0, 1) * amp[968] - std::complex<double> (0, 1) *
      amp[958] + std::complex<double> (0, 1) * amp[979] - std::complex<double>
      (0, 1) * amp[969];
  jamp[49] = +2. * (-amp[50] - amp[63] - amp[69] - std::complex<double> (0, 1)
      * amp[95] - std::complex<double> (0, 1) * amp[96]);

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



