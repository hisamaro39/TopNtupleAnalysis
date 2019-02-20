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
// Process: g g > t t~ c c~ NP<=2 @3

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
  jamp2[0] = new double[14]; 
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
  for(int i = 0; i < 14; i++ )
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
        t[0] = matrix_3_gg_ttxccx(); 

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
      t[0] = matrix_3_gg_ttxccx(); 

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
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]); 
  VVV2P0_1(w[0], w[1], pars->GC_28, pars->ZERO, pars->ZERO, w[6]); 
  FFV5P0_3(w[3], w[2], pars->GC_2, pars->ZERO, pars->ZERO, w[7]); 
  FFV1_1(w[4], w[6], pars->GC_7, pars->ZERO, pars->ZERO, w[8]); 
  FFV1_2(w[5], w[6], pars->GC_7, pars->ZERO, pars->ZERO, w[9]); 
  FFV5P0_3(w[3], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[10]); 
  FFV2_7_3(w[3], w[2], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[11]);
  VVV1P0_1(w[0], w[1], pars->GC_6, pars->ZERO, pars->ZERO, w[12]); 
  FFV1_1(w[4], w[12], pars->GC_7, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_2(w[5], w[12], pars->GC_7, pars->ZERO, pars->ZERO, w[14]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_170, pars->GC_91, pars->ZERO, pars->ZERO,
      w[15]);
  FFV3_8P0_3(w[3], w[2], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[16]);
  FFV3_8_3(w[3], w[2], pars->GC_166, pars->GC_81, pars->mdl_MZ, pars->mdl_WZ,
      w[17]);
  FFV2_3(w[3], w[2], pars->GC_138, pars->mdl_MZ, pars->mdl_WZ, w[18]); 
  VVS2_3(w[0], w[1], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[19]); 
  FFV2_7_3(w[5], w[4], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[20]);
  FFV1P0_3(w[5], w[4], pars->GC_7, pars->ZERO, pars->ZERO, w[21]); 
  FFS2_3(w[3], w[2], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[22]); 
  FFV1P0_3(w[5], w[4], pars->GC_2, pars->ZERO, pars->ZERO, w[23]); 
  FFS2_1(w[2], w[19], pars->GC_95, pars->mdl_MT, pars->mdl_WT, w[24]); 
  FFS2_2(w[3], w[19], pars->GC_95, pars->mdl_MT, pars->mdl_WT, w[25]); 
  FFV5_1(w[2], w[6], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[26]); 
  FFV5_2(w[3], w[6], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[27]); 
  FFV5_1(w[2], w[12], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[28]); 
  FFV3_8_1(w[2], w[12], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[29]);
  FFV5_2(w[3], w[12], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[30]); 
  FFV3_8_2(w[3], w[12], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[31]);
  FFFF4_5_1(w[2], w[3], w[4], pars->GC_17, pars->GC_16, pars->ZERO, pars->ZERO,
      w[32]);
  FFFF3_1(w[2], w[3], w[4], pars->GC_13, pars->ZERO, pars->ZERO, w[33]); 
  FFFF3_7_1(w[2], w[3], w[4], pars->GC_11, pars->GC_27, pars->ZERO, pars->ZERO,
      w[34]);
  FFFF4_5_4(w[5], w[2], w[3], pars->GC_17, pars->GC_16, pars->ZERO, pars->ZERO,
      w[35]);
  FFFF3_4(w[5], w[2], w[3], pars->GC_13, pars->ZERO, pars->ZERO, w[36]); 
  FFFF3_7_4(w[5], w[2], w[3], pars->GC_11, pars->GC_27, pars->ZERO, pars->ZERO,
      w[37]);
  FFFF4_5_3(w[5], w[2], w[4], pars->GC_17, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[38]);
  FFFF3_3(w[5], w[2], w[4], pars->GC_13, pars->mdl_MT, pars->mdl_WT, w[39]); 
  FFFF3_7_3(w[5], w[2], w[4], pars->GC_11, pars->GC_27, pars->mdl_MT,
      pars->mdl_WT, w[40]);
  FFFF4_5_2(w[5], w[3], w[4], pars->GC_17, pars->GC_16, pars->mdl_MT,
      pars->mdl_WT, w[41]);
  FFFF3_2(w[5], w[3], w[4], pars->GC_13, pars->mdl_MT, pars->mdl_WT, w[42]); 
  FFFF3_7_2(w[5], w[3], w[4], pars->GC_11, pars->GC_27, pars->mdl_MT,
      pars->mdl_WT, w[43]);
  FFV5_1(w[2], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[44]); 
  FFV5_2(w[3], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[45]); 
  FFV3_8_2(w[3], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[46]);
  FFV3_8_1(w[2], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[47]);
  FFV1_1(w[4], w[1], pars->GC_7, pars->ZERO, pars->ZERO, w[48]); 
  FFV5P0_3(w[3], w[44], pars->GC_2, pars->ZERO, pars->ZERO, w[49]); 
  FFV3_8P0_3(w[3], w[44], pars->GC_170, pars->GC_91, pars->ZERO, pars->ZERO,
      w[50]);
  FFV5P0_3(w[3], w[44], pars->GC_7, pars->ZERO, pars->ZERO, w[51]); 
  FFV3_8P0_3(w[3], w[44], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[52]);
  FFV3_8_3(w[3], w[44], pars->GC_166, pars->GC_81, pars->mdl_MZ, pars->mdl_WZ,
      w[53]);
  FFV2_7_3(w[3], w[44], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[54]);
  FFV2_3(w[3], w[44], pars->GC_138, pars->mdl_MZ, pars->mdl_WZ, w[55]); 
  FFV5P0_3(w[3], w[47], pars->GC_2, pars->ZERO, pars->ZERO, w[56]); 
  FFV5P0_3(w[3], w[47], pars->GC_7, pars->ZERO, pars->ZERO, w[57]); 
  FFV2_7_3(w[3], w[47], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[58]);
  FFV1_2(w[5], w[1], pars->GC_7, pars->ZERO, pars->ZERO, w[59]); 
  FFV5_1(w[44], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[60]); 
  FFV3_8_1(w[44], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[61]);
  FFS2_3(w[3], w[44], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[62]); 
  FFV5_1(w[47], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[63]); 
  FFV5_2(w[3], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[64]); 
  FFV5_1(w[2], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[65]); 
  FFV3_8_1(w[2], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[66]);
  FFV3_8_2(w[3], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[67]);
  FFV5P0_3(w[64], w[2], pars->GC_2, pars->ZERO, pars->ZERO, w[68]); 
  FFV3_8P0_3(w[64], w[2], pars->GC_170, pars->GC_91, pars->ZERO, pars->ZERO,
      w[69]);
  FFV5P0_3(w[64], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[70]); 
  FFV3_8P0_3(w[64], w[2], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[71]);
  FFV3_8_3(w[64], w[2], pars->GC_166, pars->GC_81, pars->mdl_MZ, pars->mdl_WZ,
      w[72]);
  FFV2_7_3(w[64], w[2], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[73]);
  FFV2_3(w[64], w[2], pars->GC_138, pars->mdl_MZ, pars->mdl_WZ, w[74]); 
  FFV5P0_3(w[67], w[2], pars->GC_2, pars->ZERO, pars->ZERO, w[75]); 
  FFV5P0_3(w[67], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[76]); 
  FFV2_7_3(w[67], w[2], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[77]);
  FFV5_2(w[64], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[78]); 
  FFV3_8_2(w[64], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[79]);
  FFS2_3(w[64], w[2], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[80]); 
  FFV5_2(w[67], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[81]); 
  FFV1_1(w[4], w[0], pars->GC_7, pars->ZERO, pars->ZERO, w[82]); 
  FFV1P0_3(w[5], w[82], pars->GC_2, pars->ZERO, pars->ZERO, w[83]); 
  FFV1P0_3(w[5], w[82], pars->GC_7, pars->ZERO, pars->ZERO, w[84]); 
  FFV2_7_3(w[5], w[82], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[85]);
  FFV1_1(w[82], w[1], pars->GC_7, pars->ZERO, pars->ZERO, w[86]); 
  FFVV1_2P0_3(w[3], w[2], w[1], pars->GC_146, pars->GC_85, pars->ZERO,
      pars->ZERO, w[87]);
  FFV1_2(w[5], w[0], pars->GC_7, pars->ZERO, pars->ZERO, w[88]); 
  FFV1P0_3(w[88], w[4], pars->GC_2, pars->ZERO, pars->ZERO, w[89]); 
  FFV1P0_3(w[88], w[4], pars->GC_7, pars->ZERO, pars->ZERO, w[90]); 
  FFV2_7_3(w[88], w[4], pars->GC_59, pars->GC_69, pars->mdl_MZ, pars->mdl_WZ,
      w[91]);
  FFV1_2(w[88], w[1], pars->GC_7, pars->ZERO, pars->ZERO, w[92]); 
  FFV5_1(w[65], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[93]); 
  FFV3_8_1(w[65], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[94]);
  VVS2_3(w[0], w[21], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[95]); 
  VVV2P0_1(w[0], w[21], pars->GC_28, pars->ZERO, pars->ZERO, w[96]); 
  VVV1P0_1(w[0], w[21], pars->GC_6, pars->ZERO, pars->ZERO, w[97]); 
  FFV5_1(w[66], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[98]); 
  FFV5_2(w[45], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[99]); 
  FFV3_8_2(w[45], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[100]);
  FFV5_2(w[46], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[101]); 
  VVS2P0_1(w[0], w[22], pars->GC_78, pars->ZERO, pars->ZERO, w[102]); 
  FFV1_1(w[48], w[0], pars->GC_7, pars->ZERO, pars->ZERO, w[103]); 
  VVV2P0_1(w[0], w[10], pars->GC_28, pars->ZERO, pars->ZERO, w[104]); 
  VVV1P0_1(w[0], w[10], pars->GC_6, pars->ZERO, pars->ZERO, w[105]); 
  VVV1P0_1(w[0], w[16], pars->GC_6, pars->ZERO, pars->ZERO, w[106]); 
  FFV1_2(w[59], w[0], pars->GC_7, pars->ZERO, pars->ZERO, w[107]); 
  FFVV1_2_1(w[2], w[0], w[1], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[108]);
  FFVV1_2_2(w[3], w[0], w[1], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[109]);
  FFVV1_2P0_3(w[3], w[2], w[0], pars->GC_146, pars->GC_85, pars->ZERO,
      pars->ZERO, w[110]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[5], w[8], w[7], pars->GC_2, amp[0]); 
  FFV1_0(w[9], w[4], w[7], pars->GC_2, amp[1]); 
  FFV1_0(w[5], w[8], w[10], pars->GC_7, amp[2]); 
  FFV1_0(w[9], w[4], w[10], pars->GC_7, amp[3]); 
  FFV2_7_0(w[5], w[8], w[11], pars->GC_59, pars->GC_69, amp[4]); 
  FFV2_7_0(w[9], w[4], w[11], pars->GC_59, pars->GC_69, amp[5]); 
  FFV1_0(w[5], w[13], w[7], pars->GC_2, amp[6]); 
  FFV1_0(w[14], w[4], w[7], pars->GC_2, amp[7]); 
  FFV1_0(w[5], w[13], w[15], pars->GC_2, amp[8]); 
  FFV1_0(w[14], w[4], w[15], pars->GC_2, amp[9]); 
  FFV1_0(w[5], w[13], w[10], pars->GC_7, amp[10]); 
  FFV1_0(w[14], w[4], w[10], pars->GC_7, amp[11]); 
  FFV1_0(w[5], w[13], w[16], pars->GC_7, amp[12]); 
  FFV1_0(w[14], w[4], w[16], pars->GC_7, amp[13]); 
  FFV2_7_0(w[5], w[13], w[17], pars->GC_59, pars->GC_69, amp[14]); 
  FFV2_7_0(w[14], w[4], w[17], pars->GC_59, pars->GC_69, amp[15]); 
  FFV2_7_0(w[5], w[13], w[11], pars->GC_59, pars->GC_69, amp[16]); 
  FFV2_7_0(w[14], w[4], w[11], pars->GC_59, pars->GC_69, amp[17]); 
  FFV2_7_0(w[5], w[13], w[18], pars->GC_59, pars->GC_69, amp[18]); 
  FFV2_7_0(w[14], w[4], w[18], pars->GC_59, pars->GC_69, amp[19]); 
  VVS1_0(w[11], w[20], w[19], pars->GC_93, amp[20]); 
  VVV1_0(w[6], w[10], w[21], pars->GC_6, amp[21]); 
  VVS2_0(w[12], w[21], w[22], pars->GC_78, amp[22]); 
  VVV2_0(w[12], w[10], w[21], pars->GC_28, amp[23]); 
  VVV1_0(w[12], w[10], w[21], pars->GC_6, amp[24]); 
  VVV1_0(w[12], w[16], w[21], pars->GC_6, amp[25]); 
  FFV5_0(w[3], w[24], w[23], pars->GC_2, amp[26]); 
  FFV5_0(w[25], w[2], w[23], pars->GC_2, amp[27]); 
  FFV5_0(w[3], w[24], w[21], pars->GC_7, amp[28]); 
  FFV5_0(w[25], w[2], w[21], pars->GC_7, amp[29]); 
  FFV2_7_0(w[3], w[24], w[20], pars->GC_59, pars->GC_69, amp[30]); 
  FFV2_7_0(w[25], w[2], w[20], pars->GC_59, pars->GC_69, amp[31]); 
  FFV5_0(w[3], w[26], w[23], pars->GC_2, amp[32]); 
  FFV5_0(w[27], w[2], w[23], pars->GC_2, amp[33]); 
  FFV5_0(w[3], w[26], w[21], pars->GC_7, amp[34]); 
  FFV5_0(w[27], w[2], w[21], pars->GC_7, amp[35]); 
  FFV2_7_0(w[3], w[26], w[20], pars->GC_59, pars->GC_69, amp[36]); 
  FFV2_7_0(w[27], w[2], w[20], pars->GC_59, pars->GC_69, amp[37]); 
  FFV5_0(w[3], w[28], w[23], pars->GC_2, amp[38]); 
  FFV3_8_0(w[3], w[28], w[23], pars->GC_170, pars->GC_91, amp[39]); 
  FFV5_0(w[3], w[29], w[23], pars->GC_2, amp[40]); 
  FFV5_0(w[30], w[2], w[23], pars->GC_2, amp[41]); 
  FFV3_8_0(w[30], w[2], w[23], pars->GC_170, pars->GC_91, amp[42]); 
  FFV5_0(w[31], w[2], w[23], pars->GC_2, amp[43]); 
  FFVV1_2_0(w[3], w[2], w[12], w[21], pars->GC_146, pars->GC_85, amp[44]); 
  FFV5_0(w[3], w[28], w[21], pars->GC_7, amp[45]); 
  FFV3_8_0(w[3], w[28], w[21], pars->GC_145, pars->GC_79, amp[46]); 
  FFV5_0(w[3], w[29], w[21], pars->GC_7, amp[47]); 
  FFV5_0(w[30], w[2], w[21], pars->GC_7, amp[48]); 
  FFV3_8_0(w[30], w[2], w[21], pars->GC_145, pars->GC_79, amp[49]); 
  FFV5_0(w[31], w[2], w[21], pars->GC_7, amp[50]); 
  FFV3_8_0(w[3], w[28], w[20], pars->GC_166, pars->GC_81, amp[51]); 
  FFV2_7_0(w[3], w[28], w[20], pars->GC_59, pars->GC_69, amp[52]); 
  FFV2_0(w[3], w[28], w[20], pars->GC_138, amp[53]); 
  FFV2_7_0(w[3], w[29], w[20], pars->GC_59, pars->GC_69, amp[54]); 
  FFV3_8_0(w[30], w[2], w[20], pars->GC_166, pars->GC_81, amp[55]); 
  FFV2_7_0(w[30], w[2], w[20], pars->GC_59, pars->GC_69, amp[56]); 
  FFV2_0(w[30], w[2], w[20], pars->GC_138, amp[57]); 
  FFV2_7_0(w[31], w[2], w[20], pars->GC_59, pars->GC_69, amp[58]); 
  FFV1_0(w[5], w[32], w[12], pars->GC_7, amp[59]); 
  FFV1_0(w[5], w[33], w[12], pars->GC_7, amp[60]); 
  FFV1_0(w[5], w[34], w[12], pars->GC_7, amp[61]); 
  FFV1_0(w[35], w[4], w[12], pars->GC_7, amp[62]); 
  FFV1_0(w[36], w[4], w[12], pars->GC_7, amp[63]); 
  FFV1_0(w[37], w[4], w[12], pars->GC_7, amp[64]); 
  FFV5_0(w[3], w[38], w[12], pars->GC_7, amp[65]); 
  FFV5_0(w[3], w[39], w[12], pars->GC_7, amp[66]); 
  FFV5_0(w[3], w[40], w[12], pars->GC_7, amp[67]); 
  FFV5_0(w[41], w[2], w[12], pars->GC_7, amp[68]); 
  FFV5_0(w[42], w[2], w[12], pars->GC_7, amp[69]); 
  FFV5_0(w[43], w[2], w[12], pars->GC_7, amp[70]); 
  FFFF4_5_0(w[5], w[44], w[45], w[4], pars->GC_17, pars->GC_16, amp[71]); 
  FFFF3_0(w[5], w[44], w[45], w[4], pars->GC_13, amp[72]); 
  FFFF3_7_0(w[5], w[44], w[45], w[4], pars->GC_11, pars->GC_27, amp[73]); 
  FFV5_0(w[45], w[44], w[23], pars->GC_2, amp[74]); 
  FFV3_8_0(w[45], w[44], w[23], pars->GC_170, pars->GC_91, amp[75]); 
  FFV5_0(w[45], w[44], w[21], pars->GC_7, amp[76]); 
  FFV3_8_0(w[45], w[44], w[21], pars->GC_145, pars->GC_79, amp[77]); 
  FFV3_8_0(w[45], w[44], w[20], pars->GC_166, pars->GC_81, amp[78]); 
  FFV2_7_0(w[45], w[44], w[20], pars->GC_59, pars->GC_69, amp[79]); 
  FFV2_0(w[45], w[44], w[20], pars->GC_138, amp[80]); 
  FFV5_0(w[46], w[44], w[23], pars->GC_2, amp[81]); 
  FFV5_0(w[46], w[44], w[21], pars->GC_7, amp[82]); 
  FFV2_7_0(w[46], w[44], w[20], pars->GC_59, pars->GC_69, amp[83]); 
  FFV5_0(w[45], w[47], w[23], pars->GC_2, amp[84]); 
  FFV5_0(w[45], w[47], w[21], pars->GC_7, amp[85]); 
  FFV2_7_0(w[45], w[47], w[20], pars->GC_59, pars->GC_69, amp[86]); 
  FFFF4_5_0(w[5], w[44], w[3], w[48], pars->GC_17, pars->GC_16, amp[87]); 
  FFFF3_0(w[5], w[44], w[3], w[48], pars->GC_13, amp[88]); 
  FFFF3_7_0(w[5], w[44], w[3], w[48], pars->GC_11, pars->GC_27, amp[89]); 
  FFV1_0(w[5], w[48], w[49], pars->GC_2, amp[90]); 
  FFV1_0(w[5], w[48], w[50], pars->GC_2, amp[91]); 
  FFV1_0(w[5], w[48], w[51], pars->GC_7, amp[92]); 
  FFV1_0(w[5], w[48], w[52], pars->GC_7, amp[93]); 
  FFV2_7_0(w[5], w[48], w[53], pars->GC_59, pars->GC_69, amp[94]); 
  FFV2_7_0(w[5], w[48], w[54], pars->GC_59, pars->GC_69, amp[95]); 
  FFV2_7_0(w[5], w[48], w[55], pars->GC_59, pars->GC_69, amp[96]); 
  FFV1_0(w[5], w[48], w[56], pars->GC_2, amp[97]); 
  FFV1_0(w[5], w[48], w[57], pars->GC_7, amp[98]); 
  FFV2_7_0(w[5], w[48], w[58], pars->GC_59, pars->GC_69, amp[99]); 
  FFFF4_5_0(w[59], w[44], w[3], w[4], pars->GC_17, pars->GC_16, amp[100]); 
  FFFF3_0(w[59], w[44], w[3], w[4], pars->GC_13, amp[101]); 
  FFFF3_7_0(w[59], w[44], w[3], w[4], pars->GC_11, pars->GC_27, amp[102]); 
  FFV1_0(w[59], w[4], w[49], pars->GC_2, amp[103]); 
  FFV1_0(w[59], w[4], w[50], pars->GC_2, amp[104]); 
  FFV1_0(w[59], w[4], w[51], pars->GC_7, amp[105]); 
  FFV1_0(w[59], w[4], w[52], pars->GC_7, amp[106]); 
  FFV2_7_0(w[59], w[4], w[53], pars->GC_59, pars->GC_69, amp[107]); 
  FFV2_7_0(w[59], w[4], w[54], pars->GC_59, pars->GC_69, amp[108]); 
  FFV2_7_0(w[59], w[4], w[55], pars->GC_59, pars->GC_69, amp[109]); 
  FFV1_0(w[59], w[4], w[56], pars->GC_2, amp[110]); 
  FFV1_0(w[59], w[4], w[57], pars->GC_7, amp[111]); 
  FFV2_7_0(w[59], w[4], w[58], pars->GC_59, pars->GC_69, amp[112]); 
  FFV5_0(w[3], w[60], w[23], pars->GC_2, amp[113]); 
  FFV3_8_0(w[3], w[60], w[23], pars->GC_170, pars->GC_91, amp[114]); 
  FFV5_0(w[3], w[61], w[23], pars->GC_2, amp[115]); 
  FFVV1_2_0(w[3], w[44], w[1], w[21], pars->GC_146, pars->GC_85, amp[116]); 
  FFV5_0(w[3], w[60], w[21], pars->GC_7, amp[117]); 
  FFV3_8_0(w[3], w[60], w[21], pars->GC_145, pars->GC_79, amp[118]); 
  FFV5_0(w[3], w[61], w[21], pars->GC_7, amp[119]); 
  VVS2_0(w[1], w[21], w[62], pars->GC_78, amp[120]); 
  VVV2_0(w[1], w[21], w[51], pars->GC_28, amp[121]); 
  VVV1_0(w[1], w[21], w[51], pars->GC_6, amp[122]); 
  VVV1_0(w[1], w[21], w[52], pars->GC_6, amp[123]); 
  FFV3_8_0(w[3], w[60], w[20], pars->GC_166, pars->GC_81, amp[124]); 
  FFV2_7_0(w[3], w[60], w[20], pars->GC_59, pars->GC_69, amp[125]); 
  FFV2_0(w[3], w[60], w[20], pars->GC_138, amp[126]); 
  FFV2_7_0(w[3], w[61], w[20], pars->GC_59, pars->GC_69, amp[127]); 
  FFV5_0(w[3], w[63], w[23], pars->GC_2, amp[128]); 
  FFV5_0(w[3], w[63], w[21], pars->GC_7, amp[129]); 
  VVV1_0(w[1], w[21], w[57], pars->GC_6, amp[130]); 
  FFV2_7_0(w[3], w[63], w[20], pars->GC_59, pars->GC_69, amp[131]); 
  FFV5_0(w[41], w[44], w[1], pars->GC_7, amp[132]); 
  FFV5_0(w[42], w[44], w[1], pars->GC_7, amp[133]); 
  FFV5_0(w[43], w[44], w[1], pars->GC_7, amp[134]); 
  FFFF4_5_0(w[5], w[65], w[64], w[4], pars->GC_17, pars->GC_16, amp[135]); 
  FFFF3_0(w[5], w[65], w[64], w[4], pars->GC_13, amp[136]); 
  FFFF3_7_0(w[5], w[65], w[64], w[4], pars->GC_11, pars->GC_27, amp[137]); 
  FFV5_0(w[64], w[65], w[23], pars->GC_2, amp[138]); 
  FFV3_8_0(w[64], w[65], w[23], pars->GC_170, pars->GC_91, amp[139]); 
  FFV5_0(w[64], w[65], w[21], pars->GC_7, amp[140]); 
  FFV3_8_0(w[64], w[65], w[21], pars->GC_145, pars->GC_79, amp[141]); 
  FFV3_8_0(w[64], w[65], w[20], pars->GC_166, pars->GC_81, amp[142]); 
  FFV2_7_0(w[64], w[65], w[20], pars->GC_59, pars->GC_69, amp[143]); 
  FFV2_0(w[64], w[65], w[20], pars->GC_138, amp[144]); 
  FFV5_0(w[64], w[66], w[23], pars->GC_2, amp[145]); 
  FFV5_0(w[64], w[66], w[21], pars->GC_7, amp[146]); 
  FFV2_7_0(w[64], w[66], w[20], pars->GC_59, pars->GC_69, amp[147]); 
  FFV5_0(w[67], w[65], w[23], pars->GC_2, amp[148]); 
  FFV5_0(w[67], w[65], w[21], pars->GC_7, amp[149]); 
  FFV2_7_0(w[67], w[65], w[20], pars->GC_59, pars->GC_69, amp[150]); 
  FFFF4_5_0(w[5], w[2], w[64], w[48], pars->GC_17, pars->GC_16, amp[151]); 
  FFFF3_0(w[5], w[2], w[64], w[48], pars->GC_13, amp[152]); 
  FFFF3_7_0(w[5], w[2], w[64], w[48], pars->GC_11, pars->GC_27, amp[153]); 
  FFV1_0(w[5], w[48], w[68], pars->GC_2, amp[154]); 
  FFV1_0(w[5], w[48], w[69], pars->GC_2, amp[155]); 
  FFV1_0(w[5], w[48], w[70], pars->GC_7, amp[156]); 
  FFV1_0(w[5], w[48], w[71], pars->GC_7, amp[157]); 
  FFV2_7_0(w[5], w[48], w[72], pars->GC_59, pars->GC_69, amp[158]); 
  FFV2_7_0(w[5], w[48], w[73], pars->GC_59, pars->GC_69, amp[159]); 
  FFV2_7_0(w[5], w[48], w[74], pars->GC_59, pars->GC_69, amp[160]); 
  FFV1_0(w[5], w[48], w[75], pars->GC_2, amp[161]); 
  FFV1_0(w[5], w[48], w[76], pars->GC_7, amp[162]); 
  FFV2_7_0(w[5], w[48], w[77], pars->GC_59, pars->GC_69, amp[163]); 
  FFFF4_5_0(w[59], w[2], w[64], w[4], pars->GC_17, pars->GC_16, amp[164]); 
  FFFF3_0(w[59], w[2], w[64], w[4], pars->GC_13, amp[165]); 
  FFFF3_7_0(w[59], w[2], w[64], w[4], pars->GC_11, pars->GC_27, amp[166]); 
  FFV1_0(w[59], w[4], w[68], pars->GC_2, amp[167]); 
  FFV1_0(w[59], w[4], w[69], pars->GC_2, amp[168]); 
  FFV1_0(w[59], w[4], w[70], pars->GC_7, amp[169]); 
  FFV1_0(w[59], w[4], w[71], pars->GC_7, amp[170]); 
  FFV2_7_0(w[59], w[4], w[72], pars->GC_59, pars->GC_69, amp[171]); 
  FFV2_7_0(w[59], w[4], w[73], pars->GC_59, pars->GC_69, amp[172]); 
  FFV2_7_0(w[59], w[4], w[74], pars->GC_59, pars->GC_69, amp[173]); 
  FFV1_0(w[59], w[4], w[75], pars->GC_2, amp[174]); 
  FFV1_0(w[59], w[4], w[76], pars->GC_7, amp[175]); 
  FFV2_7_0(w[59], w[4], w[77], pars->GC_59, pars->GC_69, amp[176]); 
  FFV5_0(w[78], w[2], w[23], pars->GC_2, amp[177]); 
  FFV3_8_0(w[78], w[2], w[23], pars->GC_170, pars->GC_91, amp[178]); 
  FFV5_0(w[79], w[2], w[23], pars->GC_2, amp[179]); 
  FFVV1_2_0(w[64], w[2], w[1], w[21], pars->GC_146, pars->GC_85, amp[180]); 
  FFV5_0(w[78], w[2], w[21], pars->GC_7, amp[181]); 
  FFV3_8_0(w[78], w[2], w[21], pars->GC_145, pars->GC_79, amp[182]); 
  FFV5_0(w[79], w[2], w[21], pars->GC_7, amp[183]); 
  VVS2_0(w[1], w[21], w[80], pars->GC_78, amp[184]); 
  VVV2_0(w[1], w[21], w[70], pars->GC_28, amp[185]); 
  VVV1_0(w[1], w[21], w[70], pars->GC_6, amp[186]); 
  VVV1_0(w[1], w[21], w[71], pars->GC_6, amp[187]); 
  FFV3_8_0(w[78], w[2], w[20], pars->GC_166, pars->GC_81, amp[188]); 
  FFV2_7_0(w[78], w[2], w[20], pars->GC_59, pars->GC_69, amp[189]); 
  FFV2_0(w[78], w[2], w[20], pars->GC_138, amp[190]); 
  FFV2_7_0(w[79], w[2], w[20], pars->GC_59, pars->GC_69, amp[191]); 
  FFV5_0(w[81], w[2], w[23], pars->GC_2, amp[192]); 
  FFV5_0(w[81], w[2], w[21], pars->GC_7, amp[193]); 
  VVV1_0(w[1], w[21], w[76], pars->GC_6, amp[194]); 
  FFV2_7_0(w[81], w[2], w[20], pars->GC_59, pars->GC_69, amp[195]); 
  FFV5_0(w[64], w[38], w[1], pars->GC_7, amp[196]); 
  FFV5_0(w[64], w[39], w[1], pars->GC_7, amp[197]); 
  FFV5_0(w[64], w[40], w[1], pars->GC_7, amp[198]); 
  FFFF4_5_0(w[5], w[65], w[3], w[82], pars->GC_17, pars->GC_16, amp[199]); 
  FFFF3_0(w[5], w[65], w[3], w[82], pars->GC_13, amp[200]); 
  FFFF3_7_0(w[5], w[65], w[3], w[82], pars->GC_11, pars->GC_27, amp[201]); 
  FFV5_0(w[3], w[65], w[83], pars->GC_2, amp[202]); 
  FFV3_8_0(w[3], w[65], w[83], pars->GC_170, pars->GC_91, amp[203]); 
  FFV5_0(w[3], w[65], w[84], pars->GC_7, amp[204]); 
  FFV3_8_0(w[3], w[65], w[84], pars->GC_145, pars->GC_79, amp[205]); 
  FFV3_8_0(w[3], w[65], w[85], pars->GC_166, pars->GC_81, amp[206]); 
  FFV2_7_0(w[3], w[65], w[85], pars->GC_59, pars->GC_69, amp[207]); 
  FFV2_0(w[3], w[65], w[85], pars->GC_138, amp[208]); 
  FFV5_0(w[3], w[66], w[83], pars->GC_2, amp[209]); 
  FFV5_0(w[3], w[66], w[84], pars->GC_7, amp[210]); 
  FFV2_7_0(w[3], w[66], w[85], pars->GC_59, pars->GC_69, amp[211]); 
  FFFF4_5_0(w[5], w[2], w[45], w[82], pars->GC_17, pars->GC_16, amp[212]); 
  FFFF3_0(w[5], w[2], w[45], w[82], pars->GC_13, amp[213]); 
  FFFF3_7_0(w[5], w[2], w[45], w[82], pars->GC_11, pars->GC_27, amp[214]); 
  FFV5_0(w[45], w[2], w[83], pars->GC_2, amp[215]); 
  FFV3_8_0(w[45], w[2], w[83], pars->GC_170, pars->GC_91, amp[216]); 
  FFV5_0(w[45], w[2], w[84], pars->GC_7, amp[217]); 
  FFV3_8_0(w[45], w[2], w[84], pars->GC_145, pars->GC_79, amp[218]); 
  FFV3_8_0(w[45], w[2], w[85], pars->GC_166, pars->GC_81, amp[219]); 
  FFV2_7_0(w[45], w[2], w[85], pars->GC_59, pars->GC_69, amp[220]); 
  FFV2_0(w[45], w[2], w[85], pars->GC_138, amp[221]); 
  FFV5_0(w[46], w[2], w[83], pars->GC_2, amp[222]); 
  FFV5_0(w[46], w[2], w[84], pars->GC_7, amp[223]); 
  FFV2_7_0(w[46], w[2], w[85], pars->GC_59, pars->GC_69, amp[224]); 
  FFFF4_5_0(w[59], w[2], w[3], w[82], pars->GC_17, pars->GC_16, amp[225]); 
  FFFF3_0(w[59], w[2], w[3], w[82], pars->GC_13, amp[226]); 
  FFFF3_7_0(w[59], w[2], w[3], w[82], pars->GC_11, pars->GC_27, amp[227]); 
  FFV1_0(w[59], w[82], w[7], pars->GC_2, amp[228]); 
  FFV1_0(w[59], w[82], w[15], pars->GC_2, amp[229]); 
  FFV1_0(w[59], w[82], w[10], pars->GC_7, amp[230]); 
  FFV1_0(w[59], w[82], w[16], pars->GC_7, amp[231]); 
  FFV2_7_0(w[59], w[82], w[17], pars->GC_59, pars->GC_69, amp[232]); 
  FFV2_7_0(w[59], w[82], w[11], pars->GC_59, pars->GC_69, amp[233]); 
  FFV2_7_0(w[59], w[82], w[18], pars->GC_59, pars->GC_69, amp[234]); 
  VVS2_0(w[1], w[84], w[22], pars->GC_78, amp[235]); 
  FFV1_0(w[5], w[86], w[7], pars->GC_2, amp[236]); 
  FFV1_0(w[5], w[86], w[15], pars->GC_2, amp[237]); 
  FFV1_0(w[5], w[86], w[10], pars->GC_7, amp[238]); 
  VVV2_0(w[1], w[10], w[84], pars->GC_28, amp[239]); 
  VVV1_0(w[1], w[10], w[84], pars->GC_6, amp[240]); 
  FFV1_0(w[5], w[86], w[16], pars->GC_7, amp[241]); 
  VVV1_0(w[1], w[16], w[84], pars->GC_6, amp[242]); 
  FFV2_7_0(w[5], w[86], w[17], pars->GC_59, pars->GC_69, amp[243]); 
  FFV2_7_0(w[5], w[86], w[11], pars->GC_59, pars->GC_69, amp[244]); 
  FFV2_7_0(w[5], w[86], w[18], pars->GC_59, pars->GC_69, amp[245]); 
  FFV1_0(w[5], w[82], w[87], pars->GC_7, amp[246]); 
  FFV1_0(w[35], w[82], w[1], pars->GC_7, amp[247]); 
  FFV1_0(w[36], w[82], w[1], pars->GC_7, amp[248]); 
  FFV1_0(w[37], w[82], w[1], pars->GC_7, amp[249]); 
  FFFF4_5_0(w[88], w[65], w[3], w[4], pars->GC_17, pars->GC_16, amp[250]); 
  FFFF3_0(w[88], w[65], w[3], w[4], pars->GC_13, amp[251]); 
  FFFF3_7_0(w[88], w[65], w[3], w[4], pars->GC_11, pars->GC_27, amp[252]); 
  FFV5_0(w[3], w[65], w[89], pars->GC_2, amp[253]); 
  FFV3_8_0(w[3], w[65], w[89], pars->GC_170, pars->GC_91, amp[254]); 
  FFV5_0(w[3], w[65], w[90], pars->GC_7, amp[255]); 
  FFV3_8_0(w[3], w[65], w[90], pars->GC_145, pars->GC_79, amp[256]); 
  FFV3_8_0(w[3], w[65], w[91], pars->GC_166, pars->GC_81, amp[257]); 
  FFV2_7_0(w[3], w[65], w[91], pars->GC_59, pars->GC_69, amp[258]); 
  FFV2_0(w[3], w[65], w[91], pars->GC_138, amp[259]); 
  FFV5_0(w[3], w[66], w[89], pars->GC_2, amp[260]); 
  FFV5_0(w[3], w[66], w[90], pars->GC_7, amp[261]); 
  FFV2_7_0(w[3], w[66], w[91], pars->GC_59, pars->GC_69, amp[262]); 
  FFFF4_5_0(w[88], w[2], w[45], w[4], pars->GC_17, pars->GC_16, amp[263]); 
  FFFF3_0(w[88], w[2], w[45], w[4], pars->GC_13, amp[264]); 
  FFFF3_7_0(w[88], w[2], w[45], w[4], pars->GC_11, pars->GC_27, amp[265]); 
  FFV5_0(w[45], w[2], w[89], pars->GC_2, amp[266]); 
  FFV3_8_0(w[45], w[2], w[89], pars->GC_170, pars->GC_91, amp[267]); 
  FFV5_0(w[45], w[2], w[90], pars->GC_7, amp[268]); 
  FFV3_8_0(w[45], w[2], w[90], pars->GC_145, pars->GC_79, amp[269]); 
  FFV3_8_0(w[45], w[2], w[91], pars->GC_166, pars->GC_81, amp[270]); 
  FFV2_7_0(w[45], w[2], w[91], pars->GC_59, pars->GC_69, amp[271]); 
  FFV2_0(w[45], w[2], w[91], pars->GC_138, amp[272]); 
  FFV5_0(w[46], w[2], w[89], pars->GC_2, amp[273]); 
  FFV5_0(w[46], w[2], w[90], pars->GC_7, amp[274]); 
  FFV2_7_0(w[46], w[2], w[91], pars->GC_59, pars->GC_69, amp[275]); 
  FFFF4_5_0(w[88], w[2], w[3], w[48], pars->GC_17, pars->GC_16, amp[276]); 
  FFFF3_0(w[88], w[2], w[3], w[48], pars->GC_13, amp[277]); 
  FFFF3_7_0(w[88], w[2], w[3], w[48], pars->GC_11, pars->GC_27, amp[278]); 
  FFV1_0(w[88], w[48], w[7], pars->GC_2, amp[279]); 
  FFV1_0(w[88], w[48], w[15], pars->GC_2, amp[280]); 
  FFV1_0(w[88], w[48], w[10], pars->GC_7, amp[281]); 
  FFV1_0(w[88], w[48], w[16], pars->GC_7, amp[282]); 
  FFV2_7_0(w[88], w[48], w[17], pars->GC_59, pars->GC_69, amp[283]); 
  FFV2_7_0(w[88], w[48], w[11], pars->GC_59, pars->GC_69, amp[284]); 
  FFV2_7_0(w[88], w[48], w[18], pars->GC_59, pars->GC_69, amp[285]); 
  VVS2_0(w[1], w[90], w[22], pars->GC_78, amp[286]); 
  FFV1_0(w[92], w[4], w[7], pars->GC_2, amp[287]); 
  FFV1_0(w[92], w[4], w[15], pars->GC_2, amp[288]); 
  FFV1_0(w[92], w[4], w[10], pars->GC_7, amp[289]); 
  VVV2_0(w[1], w[10], w[90], pars->GC_28, amp[290]); 
  VVV1_0(w[1], w[10], w[90], pars->GC_6, amp[291]); 
  FFV1_0(w[92], w[4], w[16], pars->GC_7, amp[292]); 
  VVV1_0(w[1], w[16], w[90], pars->GC_6, amp[293]); 
  FFV2_7_0(w[92], w[4], w[17], pars->GC_59, pars->GC_69, amp[294]); 
  FFV2_7_0(w[92], w[4], w[11], pars->GC_59, pars->GC_69, amp[295]); 
  FFV2_7_0(w[92], w[4], w[18], pars->GC_59, pars->GC_69, amp[296]); 
  FFV1_0(w[88], w[4], w[87], pars->GC_7, amp[297]); 
  FFV1_0(w[88], w[32], w[1], pars->GC_7, amp[298]); 
  FFV1_0(w[88], w[33], w[1], pars->GC_7, amp[299]); 
  FFV1_0(w[88], w[34], w[1], pars->GC_7, amp[300]); 
  FFV5_0(w[3], w[93], w[23], pars->GC_2, amp[301]); 
  FFV3_8_0(w[3], w[93], w[23], pars->GC_170, pars->GC_91, amp[302]); 
  FFV5_0(w[3], w[94], w[23], pars->GC_2, amp[303]); 
  FFVV1_2_0(w[3], w[65], w[0], w[21], pars->GC_146, pars->GC_85, amp[304]); 
  FFV5_0(w[3], w[93], w[21], pars->GC_7, amp[305]); 
  FFV3_8_0(w[3], w[93], w[21], pars->GC_145, pars->GC_79, amp[306]); 
  FFV5_0(w[3], w[94], w[21], pars->GC_7, amp[307]); 
  FFS2_0(w[3], w[65], w[95], pars->GC_95, amp[308]); 
  FFV5_0(w[3], w[65], w[96], pars->GC_7, amp[309]); 
  FFV5_0(w[3], w[65], w[97], pars->GC_7, amp[310]); 
  FFV3_8_0(w[3], w[65], w[97], pars->GC_145, pars->GC_79, amp[311]); 
  FFV3_8_0(w[3], w[93], w[20], pars->GC_166, pars->GC_81, amp[312]); 
  FFV2_7_0(w[3], w[93], w[20], pars->GC_59, pars->GC_69, amp[313]); 
  FFV2_0(w[3], w[93], w[20], pars->GC_138, amp[314]); 
  FFV2_7_0(w[3], w[94], w[20], pars->GC_59, pars->GC_69, amp[315]); 
  FFV5_0(w[3], w[98], w[23], pars->GC_2, amp[316]); 
  FFV5_0(w[3], w[98], w[21], pars->GC_7, amp[317]); 
  FFV5_0(w[3], w[66], w[97], pars->GC_7, amp[318]); 
  FFV2_7_0(w[3], w[98], w[20], pars->GC_59, pars->GC_69, amp[319]); 
  FFV5_0(w[41], w[65], w[0], pars->GC_7, amp[320]); 
  FFV5_0(w[42], w[65], w[0], pars->GC_7, amp[321]); 
  FFV5_0(w[43], w[65], w[0], pars->GC_7, amp[322]); 
  FFV5_0(w[99], w[2], w[23], pars->GC_2, amp[323]); 
  FFV3_8_0(w[99], w[2], w[23], pars->GC_170, pars->GC_91, amp[324]); 
  FFV5_0(w[100], w[2], w[23], pars->GC_2, amp[325]); 
  FFVV1_2_0(w[45], w[2], w[0], w[21], pars->GC_146, pars->GC_85, amp[326]); 
  FFV5_0(w[99], w[2], w[21], pars->GC_7, amp[327]); 
  FFV3_8_0(w[99], w[2], w[21], pars->GC_145, pars->GC_79, amp[328]); 
  FFV5_0(w[100], w[2], w[21], pars->GC_7, amp[329]); 
  FFS2_0(w[45], w[2], w[95], pars->GC_95, amp[330]); 
  FFV5_0(w[45], w[2], w[96], pars->GC_7, amp[331]); 
  FFV5_0(w[45], w[2], w[97], pars->GC_7, amp[332]); 
  FFV3_8_0(w[45], w[2], w[97], pars->GC_145, pars->GC_79, amp[333]); 
  FFV3_8_0(w[99], w[2], w[20], pars->GC_166, pars->GC_81, amp[334]); 
  FFV2_7_0(w[99], w[2], w[20], pars->GC_59, pars->GC_69, amp[335]); 
  FFV2_0(w[99], w[2], w[20], pars->GC_138, amp[336]); 
  FFV2_7_0(w[100], w[2], w[20], pars->GC_59, pars->GC_69, amp[337]); 
  FFV5_0(w[101], w[2], w[23], pars->GC_2, amp[338]); 
  FFV5_0(w[101], w[2], w[21], pars->GC_7, amp[339]); 
  FFV5_0(w[46], w[2], w[97], pars->GC_7, amp[340]); 
  FFV2_7_0(w[101], w[2], w[20], pars->GC_59, pars->GC_69, amp[341]); 
  FFV5_0(w[45], w[38], w[0], pars->GC_7, amp[342]); 
  FFV5_0(w[45], w[39], w[0], pars->GC_7, amp[343]); 
  FFV5_0(w[45], w[40], w[0], pars->GC_7, amp[344]); 
  FFV1_0(w[5], w[48], w[102], pars->GC_7, amp[345]); 
  FFV1_0(w[5], w[103], w[7], pars->GC_2, amp[346]); 
  FFV1_0(w[5], w[103], w[15], pars->GC_2, amp[347]); 
  FFV1_0(w[5], w[103], w[10], pars->GC_7, amp[348]); 
  FFV1_0(w[5], w[48], w[104], pars->GC_7, amp[349]); 
  FFV1_0(w[5], w[48], w[105], pars->GC_7, amp[350]); 
  FFV1_0(w[5], w[103], w[16], pars->GC_7, amp[351]); 
  FFV1_0(w[5], w[48], w[106], pars->GC_7, amp[352]); 
  FFV2_7_0(w[5], w[103], w[17], pars->GC_59, pars->GC_69, amp[353]); 
  FFV2_7_0(w[5], w[103], w[11], pars->GC_59, pars->GC_69, amp[354]); 
  FFV2_7_0(w[5], w[103], w[18], pars->GC_59, pars->GC_69, amp[355]); 
  FFV1_0(w[35], w[48], w[0], pars->GC_7, amp[356]); 
  FFV1_0(w[36], w[48], w[0], pars->GC_7, amp[357]); 
  FFV1_0(w[37], w[48], w[0], pars->GC_7, amp[358]); 
  FFV1_0(w[59], w[4], w[102], pars->GC_7, amp[359]); 
  FFV1_0(w[107], w[4], w[7], pars->GC_2, amp[360]); 
  FFV1_0(w[107], w[4], w[15], pars->GC_2, amp[361]); 
  FFV1_0(w[107], w[4], w[10], pars->GC_7, amp[362]); 
  FFV1_0(w[59], w[4], w[104], pars->GC_7, amp[363]); 
  FFV1_0(w[59], w[4], w[105], pars->GC_7, amp[364]); 
  FFV1_0(w[107], w[4], w[16], pars->GC_7, amp[365]); 
  FFV1_0(w[59], w[4], w[106], pars->GC_7, amp[366]); 
  FFV2_7_0(w[107], w[4], w[17], pars->GC_59, pars->GC_69, amp[367]); 
  FFV2_7_0(w[107], w[4], w[11], pars->GC_59, pars->GC_69, amp[368]); 
  FFV2_7_0(w[107], w[4], w[18], pars->GC_59, pars->GC_69, amp[369]); 
  FFV1_0(w[59], w[32], w[0], pars->GC_7, amp[370]); 
  FFV1_0(w[59], w[33], w[0], pars->GC_7, amp[371]); 
  FFV1_0(w[59], w[34], w[0], pars->GC_7, amp[372]); 
  VVVS1_0(w[0], w[1], w[21], w[22], pars->GC_84, amp[373]); 
  VVV1_0(w[1], w[21], w[102], pars->GC_6, amp[374]); 
  VVS2_0(w[1], w[97], w[22], pars->GC_78, amp[375]); 
  VVVV2_0(w[0], w[1], w[10], w[21], pars->GC_44, amp[376]); 
  VVVV8_0(w[0], w[1], w[10], w[21], pars->GC_44, amp[377]); 
  VVVV7_0(w[0], w[1], w[10], w[21], pars->GC_44, amp[378]); 
  VVVV1_0(w[0], w[1], w[10], w[21], pars->GC_8, amp[379]); 
  VVVV4_0(w[0], w[1], w[10], w[21], pars->GC_8, amp[380]); 
  VVVV5_0(w[0], w[1], w[10], w[21], pars->GC_8, amp[381]); 
  VVV1_0(w[1], w[21], w[104], pars->GC_6, amp[382]); 
  VVV2_0(w[1], w[21], w[105], pars->GC_28, amp[383]); 
  VVV1_0(w[1], w[21], w[105], pars->GC_6, amp[384]); 
  VVV1_0(w[1], w[10], w[96], pars->GC_6, amp[385]); 
  VVV2_0(w[1], w[10], w[97], pars->GC_28, amp[386]); 
  VVV1_0(w[1], w[10], w[97], pars->GC_6, amp[387]); 
  VVVV1_0(w[0], w[1], w[16], w[21], pars->GC_8, amp[388]); 
  VVVV4_0(w[0], w[1], w[16], w[21], pars->GC_8, amp[389]); 
  VVVV5_0(w[0], w[1], w[16], w[21], pars->GC_8, amp[390]); 
  VVV1_0(w[1], w[21], w[106], pars->GC_6, amp[391]); 
  VVV1_0(w[1], w[16], w[97], pars->GC_6, amp[392]); 
  FFV5_0(w[3], w[108], w[23], pars->GC_2, amp[393]); 
  FFV5_0(w[3], w[108], w[21], pars->GC_7, amp[394]); 
  FFV2_7_0(w[3], w[108], w[20], pars->GC_59, pars->GC_69, amp[395]); 
  FFV5_0(w[109], w[2], w[23], pars->GC_2, amp[396]); 
  FFV5_0(w[109], w[2], w[21], pars->GC_7, amp[397]); 
  FFV2_7_0(w[109], w[2], w[20], pars->GC_59, pars->GC_69, amp[398]); 
  FFV1_0(w[5], w[48], w[110], pars->GC_7, amp[399]); 
  FFV1_0(w[59], w[4], w[110], pars->GC_7, amp[400]); 
  VVV1_0(w[110], w[1], w[21], pars->GC_6, amp[401]); 
  VVV1_0(w[0], w[87], w[21], pars->GC_6, amp[402]); 

}
double CPPProcess::matrix_3_gg_ttxccx() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 403; 
  const int ncolor = 14; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1,
      1};
  static const double cf[ncolor][ncolor] = {{48, 16, 16, 6, 0, 16, -2, 0, -6,
      -2, -2, 6, 18, 6}, {16, 48, 6, 16, 16, 0, 0, -2, -2, -6, 6, -2, 6, 18},
      {16, 6, 48, 16, -2, 0, 0, 16, -2, 6, -6, -2, 6, 18}, {6, 16, 16, 48, 0,
      -2, 16, 0, 6, -2, -2, -6, 18, 6}, {0, 16, -2, 0, 48, 16, 16, 6, 0, -2,
      16, 0, 0, 6}, {16, 0, 0, -2, 16, 48, 6, 16, -2, 0, 0, 16, 6, 0}, {-2, 0,
      0, 16, 16, 6, 48, 16, 16, 0, 0, -2, 6, 0}, {0, -2, 16, 0, 6, 16, 16, 48,
      0, 16, -2, 0, 0, 6}, {-6, -2, -2, 6, 0, -2, 16, 0, 48, 16, 16, 6, 18, 6},
      {-2, -6, 6, -2, -2, 0, 0, 16, 16, 48, 6, 16, 6, 18}, {-2, 6, -6, -2, 16,
      0, 0, -2, 16, 6, 48, 16, 6, 18}, {6, -2, -2, -6, 0, 16, -2, 0, 6, 16, 16,
      48, 18, 6}, {6, 2, 2, 6, 0, 2, 2, 0, 6, 2, 2, 6, 18, 6}, {2, 6, 6, 2, 2,
      0, 0, 2, 2, 6, 6, 2, 6, 18}};

  // Calculate color flows
  jamp[0] = -std::complex<double> (0, 1) * amp[32] - std::complex<double> (0,
      1) * amp[33] + 1./6. * std::complex<double> (0, 1) * amp[34] + 1./6. *
      std::complex<double> (0, 1) * amp[35] - std::complex<double> (0, 1) *
      amp[36] - std::complex<double> (0, 1) * amp[37] - std::complex<double>
      (0, 1) * amp[38] - std::complex<double> (0, 1) * amp[39] -
      std::complex<double> (0, 1) * amp[40] - std::complex<double> (0, 1) *
      amp[41] - std::complex<double> (0, 1) * amp[42] - std::complex<double>
      (0, 1) * amp[43] + 1./6. * std::complex<double> (0, 1) * amp[45] + 1./6.
      * std::complex<double> (0, 1) * amp[46] + 1./6. * std::complex<double>
      (0, 1) * amp[47] + 1./6. * std::complex<double> (0, 1) * amp[48] + 1./6.
      * std::complex<double> (0, 1) * amp[49] + 1./6. * std::complex<double>
      (0, 1) * amp[50] - std::complex<double> (0, 1) * amp[51] -
      std::complex<double> (0, 1) * amp[52] - std::complex<double> (0, 1) *
      amp[53] - std::complex<double> (0, 1) * amp[54] - std::complex<double>
      (0, 1) * amp[55] - std::complex<double> (0, 1) * amp[56] -
      std::complex<double> (0, 1) * amp[57] - std::complex<double> (0, 1) *
      amp[58] + std::complex<double> (0, 1) * amp[66] - 1./6. *
      std::complex<double> (0, 1) * amp[67] + std::complex<double> (0, 1) *
      amp[69] - 1./6. * std::complex<double> (0, 1) * amp[70] - amp[72] + 1./6.
      * amp[73] + amp[74] + amp[75] - 1./6. * amp[76] - 1./6. * amp[77] +
      amp[78] + amp[79] + amp[80] + amp[81] - 1./6. * amp[82] + amp[83] +
      amp[84] - 1./6. * amp[85] + amp[86] + amp[113] + amp[114] + amp[115] -
      1./6. * amp[117] - 1./6. * amp[118] - 1./6. * amp[119] + amp[124] +
      amp[125] + amp[126] + amp[127] + amp[128] - 1./6. * amp[129] + amp[131] -
      amp[133] + 1./6. * amp[134] + amp[323] + amp[324] + amp[325] - 1./6. *
      amp[327] - 1./6. * amp[328] - 1./6. * amp[329] + amp[334] + amp[335] +
      amp[336] + amp[337] + amp[338] - 1./6. * amp[339] + amp[341] + 1./6. *
      amp[344] - amp[343] - std::complex<double> (0, 1) * amp[393] + 1./6. *
      std::complex<double> (0, 1) * amp[394] - std::complex<double> (0, 1) *
      amp[395] - std::complex<double> (0, 1) * amp[396] + 1./6. *
      std::complex<double> (0, 1) * amp[397] - std::complex<double> (0, 1) *
      amp[398];
  jamp[1] = -1./2. * std::complex<double> (0, 1) * amp[3] - 1./2. *
      std::complex<double> (0, 1) * amp[11] - 1./2. * std::complex<double> (0,
      1) * amp[13] + 1./2. * amp[21] + 1./2. * amp[23] + 1./2. * amp[24] +
      1./2. * amp[25] - 1./2. * std::complex<double> (0, 1) * amp[34] - 1./2. *
      amp[44] - 1./2. * std::complex<double> (0, 1) * amp[45] - 1./2. *
      std::complex<double> (0, 1) * amp[46] - 1./2. * std::complex<double> (0,
      1) * amp[47] + std::complex<double> (0, 1) * amp[59] + 1./2. *
      std::complex<double> (0, 1) * amp[61] + std::complex<double> (0, 1) *
      amp[68] + 1./2. * std::complex<double> (0, 1) * amp[70] - amp[100] -
      1./2. * amp[102] + 1./2. * amp[105] + 1./2. * amp[106] + 1./2. * amp[111]
      - 1./2. * std::complex<double> (0, 1) * amp[116] + 1./2. * amp[117] +
      1./2. * amp[118] + 1./2. * amp[119] - 1./2. * std::complex<double> (0, 1)
      * amp[121] - 1./2. * std::complex<double> (0, 1) * amp[122] - 1./2. *
      std::complex<double> (0, 1) * amp[123] + 1./2. * amp[129] - 1./2. *
      std::complex<double> (0, 1) * amp[130] - amp[132] - 1./2. * amp[134] +
      1./2. * amp[362] + 1./2. * std::complex<double> (0, 1) * amp[363] + 1./2.
      * std::complex<double> (0, 1) * amp[364] + 1./2. * amp[365] + 1./2. *
      std::complex<double> (0, 1) * amp[366] - 1./2. * amp[372] - amp[370] +
      1./2. * amp[376] + 1./2. * amp[377] + 1./2. * amp[379] + 1./2. * amp[380]
      + 1./2. * amp[382] + 1./2. * amp[383] + 1./2. * amp[384] + 1./2. *
      amp[388] + 1./2. * amp[389] + 1./2. * amp[391] - 1./2. *
      std::complex<double> (0, 1) * amp[394] + 1./2. * std::complex<double> (0,
      1) * amp[400] + 1./2. * amp[401];
  jamp[2] = -1./2. * std::complex<double> (0, 1) * amp[2] - 1./2. *
      std::complex<double> (0, 1) * amp[10] - 1./2. * std::complex<double> (0,
      1) * amp[12] - 1./2. * amp[21] - 1./2. * amp[23] - 1./2. * amp[24] -
      1./2. * amp[25] - 1./2. * std::complex<double> (0, 1) * amp[35] + 1./2. *
      amp[44] - 1./2. * std::complex<double> (0, 1) * amp[48] - 1./2. *
      std::complex<double> (0, 1) * amp[49] - 1./2. * std::complex<double> (0,
      1) * amp[50] + std::complex<double> (0, 1) * amp[62] + 1./2. *
      std::complex<double> (0, 1) * amp[64] + std::complex<double> (0, 1) *
      amp[65] + 1./2. * std::complex<double> (0, 1) * amp[67] - amp[212] -
      1./2. * amp[214] + 1./2. * amp[217] + 1./2. * amp[218] + 1./2. * amp[223]
      + 1./2. * amp[238] - 1./2. * std::complex<double> (0, 1) * amp[239] -
      1./2. * std::complex<double> (0, 1) * amp[240] + 1./2. * amp[241] - 1./2.
      * std::complex<double> (0, 1) * amp[242] - 1./2. * std::complex<double>
      (0, 1) * amp[246] - amp[247] - 1./2. * amp[249] + 1./2. *
      std::complex<double> (0, 1) * amp[326] + 1./2. * amp[327] + 1./2. *
      amp[328] + 1./2. * amp[329] + 1./2. * std::complex<double> (0, 1) *
      amp[331] + 1./2. * std::complex<double> (0, 1) * amp[332] + 1./2. *
      std::complex<double> (0, 1) * amp[333] + 1./2. * amp[339] + 1./2. *
      std::complex<double> (0, 1) * amp[340] - 1./2. * amp[344] - amp[342] -
      1./2. * amp[376] + 1./2. * amp[378] - 1./2. * amp[379] + 1./2. * amp[381]
      + 1./2. * amp[385] + 1./2. * amp[386] + 1./2. * amp[387] - 1./2. *
      amp[388] + 1./2. * amp[390] + 1./2. * amp[392] - 1./2. *
      std::complex<double> (0, 1) * amp[397] - 1./2. * amp[402];
  jamp[3] = -std::complex<double> (0, 1) * amp[0] - std::complex<double> (0, 1)
      * amp[1] + 1./6. * std::complex<double> (0, 1) * amp[2] + 1./6. *
      std::complex<double> (0, 1) * amp[3] - std::complex<double> (0, 1) *
      amp[4] - std::complex<double> (0, 1) * amp[5] - std::complex<double> (0,
      1) * amp[6] - std::complex<double> (0, 1) * amp[7] - std::complex<double>
      (0, 1) * amp[8] - std::complex<double> (0, 1) * amp[9] + 1./6. *
      std::complex<double> (0, 1) * amp[10] + 1./6. * std::complex<double> (0,
      1) * amp[11] + 1./6. * std::complex<double> (0, 1) * amp[12] + 1./6. *
      std::complex<double> (0, 1) * amp[13] - std::complex<double> (0, 1) *
      amp[14] - std::complex<double> (0, 1) * amp[15] - std::complex<double>
      (0, 1) * amp[16] - std::complex<double> (0, 1) * amp[17] -
      std::complex<double> (0, 1) * amp[18] - std::complex<double> (0, 1) *
      amp[19] - std::complex<double> (0, 1) * amp[22] + std::complex<double>
      (0, 1) * amp[60] - 1./6. * std::complex<double> (0, 1) * amp[61] +
      std::complex<double> (0, 1) * amp[63] - 1./6. * std::complex<double> (0,
      1) * amp[64] - amp[226] + 1./6. * amp[227] + amp[228] + amp[229] - 1./6.
      * amp[230] - 1./6. * amp[231] + amp[232] + amp[233] + amp[234] + amp[235]
      + amp[236] + amp[237] - 1./6. * amp[238] - 1./6. * amp[241] + amp[243] +
      amp[244] + amp[245] - amp[248] + 1./6. * amp[249] + amp[359] + amp[360] +
      amp[361] - 1./6. * amp[362] - 1./6. * amp[365] + amp[367] + amp[368] +
      amp[369] + 1./6. * amp[372] - amp[371] - std::complex<double> (0, 1) *
      amp[373] - std::complex<double> (0, 1) * amp[374] + std::complex<double>
      (0, 1) * amp[375];
  jamp[4] = -amp[88] + 1./6. * amp[89] + amp[90] + amp[91] - 1./6. * amp[92] -
      1./6. * amp[93] + amp[94] + amp[95] + amp[96] + amp[97] - 1./6. * amp[98]
      + amp[99] - amp[101] + 1./6. * amp[102] + amp[103] + amp[104] - 1./6. *
      amp[105] - 1./6. * amp[106] + amp[107] + amp[108] + amp[109] + amp[110] -
      1./6. * amp[111] + amp[112] + amp[120] - amp[152] + 1./6. * amp[153] +
      amp[154] + amp[155] - 1./6. * amp[156] - 1./6. * amp[157] + amp[158] +
      amp[159] + amp[160] + amp[161] - 1./6. * amp[162] + amp[163] - amp[165] +
      1./6. * amp[166] + amp[167] + amp[168] - 1./6. * amp[169] - 1./6. *
      amp[170] + amp[171] + amp[172] + amp[173] + amp[174] - 1./6. * amp[175] +
      amp[176] + amp[184];
  jamp[5] = -amp[71] - 1./2. * amp[73] + 1./2. * amp[76] + 1./2. * amp[77] +
      1./2. * amp[82] + 1./2. * amp[85] - amp[87] - 1./2. * amp[89] + 1./2. *
      amp[92] + 1./2. * amp[93] + 1./2. * amp[98] + 1./2. *
      std::complex<double> (0, 1) * amp[116] + 1./2. * std::complex<double> (0,
      1) * amp[121] + 1./2. * std::complex<double> (0, 1) * amp[122] + 1./2. *
      std::complex<double> (0, 1) * amp[123] + 1./2. * std::complex<double> (0,
      1) * amp[130] - amp[263] - 1./2. * amp[265] + 1./2. * amp[268] + 1./2. *
      amp[269] + 1./2. * amp[274] - amp[276] - 1./2. * amp[278] + 1./2. *
      amp[281] + 1./2. * amp[282] - 1./2. * std::complex<double> (0, 1) *
      amp[290] - 1./2. * std::complex<double> (0, 1) * amp[291] - 1./2. *
      std::complex<double> (0, 1) * amp[293] - 1./2. * std::complex<double> (0,
      1) * amp[297] - 1./2. * std::complex<double> (0, 1) * amp[326] - 1./2. *
      std::complex<double> (0, 1) * amp[331] - 1./2. * std::complex<double> (0,
      1) * amp[332] - 1./2. * std::complex<double> (0, 1) * amp[333] - 1./2. *
      std::complex<double> (0, 1) * amp[340] + 1./2. * std::complex<double> (0,
      1) * amp[349] + 1./2. * std::complex<double> (0, 1) * amp[350] + 1./2. *
      std::complex<double> (0, 1) * amp[352] - 1./2. * amp[377] - 1./2. *
      amp[378] - 1./2. * amp[380] - 1./2. * amp[381] - 1./2. * amp[382] - 1./2.
      * amp[383] - 1./2. * amp[384] - 1./2. * amp[385] - 1./2. * amp[386] -
      1./2. * amp[387] - 1./2. * amp[389] - 1./2. * amp[390] - 1./2. * amp[391]
      - 1./2. * amp[392] + 1./2. * std::complex<double> (0, 1) * amp[399] -
      1./2. * amp[401] + 1./2. * amp[402];
  jamp[6] = -amp[135] - 1./2. * amp[137] + 1./2. * amp[140] + 1./2. * amp[141]
      + 1./2. * amp[146] + 1./2. * amp[149] - amp[164] - 1./2. * amp[166] +
      1./2. * amp[169] + 1./2. * amp[170] + 1./2. * amp[175] - 1./2. *
      std::complex<double> (0, 1) * amp[180] - 1./2. * std::complex<double> (0,
      1) * amp[185] - 1./2. * std::complex<double> (0, 1) * amp[186] - 1./2. *
      std::complex<double> (0, 1) * amp[187] - 1./2. * std::complex<double> (0,
      1) * amp[194] - amp[199] - 1./2. * amp[201] + 1./2. * amp[204] + 1./2. *
      amp[205] + 1./2. * amp[210] - amp[225] - 1./2. * amp[227] + 1./2. *
      amp[230] + 1./2. * amp[231] + 1./2. * std::complex<double> (0, 1) *
      amp[239] + 1./2. * std::complex<double> (0, 1) * amp[240] + 1./2. *
      std::complex<double> (0, 1) * amp[242] + 1./2. * std::complex<double> (0,
      1) * amp[246] + 1./2. * std::complex<double> (0, 1) * amp[304] + 1./2. *
      std::complex<double> (0, 1) * amp[309] + 1./2. * std::complex<double> (0,
      1) * amp[310] + 1./2. * std::complex<double> (0, 1) * amp[311] + 1./2. *
      std::complex<double> (0, 1) * amp[318] - 1./2. * std::complex<double> (0,
      1) * amp[363] - 1./2. * std::complex<double> (0, 1) * amp[364] - 1./2. *
      std::complex<double> (0, 1) * amp[366] - 1./2. * amp[377] - 1./2. *
      amp[378] - 1./2. * amp[380] - 1./2. * amp[381] - 1./2. * amp[382] - 1./2.
      * amp[383] - 1./2. * amp[384] - 1./2. * amp[385] - 1./2. * amp[386] -
      1./2. * amp[387] - 1./2. * amp[389] - 1./2. * amp[390] - 1./2. * amp[391]
      - 1./2. * amp[392] - 1./2. * std::complex<double> (0, 1) * amp[400] -
      1./2. * amp[401] + 1./2. * amp[402];
  jamp[7] = -amp[200] + 1./6. * amp[201] + amp[202] + amp[203] - 1./6. *
      amp[204] - 1./6. * amp[205] + amp[206] + amp[207] + amp[208] + amp[209] -
      1./6. * amp[210] + amp[211] - amp[213] + 1./6. * amp[214] + amp[215] +
      amp[216] - 1./6. * amp[217] - 1./6. * amp[218] + amp[219] + amp[220] +
      amp[221] + amp[222] - 1./6. * amp[223] + amp[224] - amp[251] + 1./6. *
      amp[252] + amp[253] + amp[254] - 1./6. * amp[255] - 1./6. * amp[256] +
      amp[257] + amp[258] + amp[259] + amp[260] - 1./6. * amp[261] + amp[262] -
      amp[264] + 1./6. * amp[265] + amp[266] + amp[267] - 1./6. * amp[268] -
      1./6. * amp[269] + amp[270] + amp[271] + amp[272] + amp[273] - 1./6. *
      amp[274] + amp[275] + amp[308] + amp[330];
  jamp[8] = +std::complex<double> (0, 1) * amp[32] + std::complex<double> (0,
      1) * amp[33] - 1./6. * std::complex<double> (0, 1) * amp[34] - 1./6. *
      std::complex<double> (0, 1) * amp[35] + std::complex<double> (0, 1) *
      amp[36] + std::complex<double> (0, 1) * amp[37] + std::complex<double>
      (0, 1) * amp[38] + std::complex<double> (0, 1) * amp[39] +
      std::complex<double> (0, 1) * amp[40] + std::complex<double> (0, 1) *
      amp[41] + std::complex<double> (0, 1) * amp[42] + std::complex<double>
      (0, 1) * amp[43] - 1./6. * std::complex<double> (0, 1) * amp[45] - 1./6.
      * std::complex<double> (0, 1) * amp[46] - 1./6. * std::complex<double>
      (0, 1) * amp[47] - 1./6. * std::complex<double> (0, 1) * amp[48] - 1./6.
      * std::complex<double> (0, 1) * amp[49] - 1./6. * std::complex<double>
      (0, 1) * amp[50] + std::complex<double> (0, 1) * amp[51] +
      std::complex<double> (0, 1) * amp[52] + std::complex<double> (0, 1) *
      amp[53] + std::complex<double> (0, 1) * amp[54] + std::complex<double>
      (0, 1) * amp[55] + std::complex<double> (0, 1) * amp[56] +
      std::complex<double> (0, 1) * amp[57] + std::complex<double> (0, 1) *
      amp[58] - std::complex<double> (0, 1) * amp[66] + 1./6. *
      std::complex<double> (0, 1) * amp[67] - std::complex<double> (0, 1) *
      amp[69] + 1./6. * std::complex<double> (0, 1) * amp[70] - amp[136] +
      1./6. * amp[137] + amp[138] + amp[139] - 1./6. * amp[140] - 1./6. *
      amp[141] + amp[142] + amp[143] + amp[144] + amp[145] - 1./6. * amp[146] +
      amp[147] + amp[148] - 1./6. * amp[149] + amp[150] + amp[177] + amp[178] +
      amp[179] - 1./6. * amp[181] - 1./6. * amp[182] - 1./6. * amp[183] +
      amp[188] + amp[189] + amp[190] + amp[191] + amp[192] - 1./6. * amp[193] +
      amp[195] + 1./6. * amp[198] - amp[197] + amp[301] + amp[302] + amp[303] -
      1./6. * amp[305] - 1./6. * amp[306] - 1./6. * amp[307] + amp[312] +
      amp[313] + amp[314] + amp[315] + amp[316] - 1./6. * amp[317] + amp[319] -
      amp[321] + 1./6. * amp[322] + std::complex<double> (0, 1) * amp[393] -
      1./6. * std::complex<double> (0, 1) * amp[394] + std::complex<double> (0,
      1) * amp[395] + std::complex<double> (0, 1) * amp[396] - 1./6. *
      std::complex<double> (0, 1) * amp[397] + std::complex<double> (0, 1) *
      amp[398];
  jamp[9] = +1./2. * std::complex<double> (0, 1) * amp[3] + 1./2. *
      std::complex<double> (0, 1) * amp[11] + 1./2. * std::complex<double> (0,
      1) * amp[13] - 1./2. * amp[21] - 1./2. * amp[23] - 1./2. * amp[24] -
      1./2. * amp[25] + 1./2. * std::complex<double> (0, 1) * amp[34] + 1./2. *
      amp[44] + 1./2. * std::complex<double> (0, 1) * amp[45] + 1./2. *
      std::complex<double> (0, 1) * amp[46] + 1./2. * std::complex<double> (0,
      1) * amp[47] - std::complex<double> (0, 1) * amp[59] - 1./2. *
      std::complex<double> (0, 1) * amp[61] - std::complex<double> (0, 1) *
      amp[68] - 1./2. * std::complex<double> (0, 1) * amp[70] - amp[250] -
      1./2. * amp[252] + 1./2. * amp[255] + 1./2. * amp[256] + 1./2. * amp[261]
      + 1./2. * amp[289] + 1./2. * std::complex<double> (0, 1) * amp[290] +
      1./2. * std::complex<double> (0, 1) * amp[291] + 1./2. * amp[292] + 1./2.
      * std::complex<double> (0, 1) * amp[293] + 1./2. * std::complex<double>
      (0, 1) * amp[297] - 1./2. * amp[300] - amp[298] - 1./2. *
      std::complex<double> (0, 1) * amp[304] + 1./2. * amp[305] + 1./2. *
      amp[306] + 1./2. * amp[307] - 1./2. * std::complex<double> (0, 1) *
      amp[309] - 1./2. * std::complex<double> (0, 1) * amp[310] - 1./2. *
      std::complex<double> (0, 1) * amp[311] + 1./2. * amp[317] - 1./2. *
      std::complex<double> (0, 1) * amp[318] - amp[320] - 1./2. * amp[322] -
      1./2. * amp[376] + 1./2. * amp[378] - 1./2. * amp[379] + 1./2. * amp[381]
      + 1./2. * amp[385] + 1./2. * amp[386] + 1./2. * amp[387] - 1./2. *
      amp[388] + 1./2. * amp[390] + 1./2. * amp[392] + 1./2. *
      std::complex<double> (0, 1) * amp[394] - 1./2. * amp[402];
  jamp[10] = +1./2. * std::complex<double> (0, 1) * amp[2] + 1./2. *
      std::complex<double> (0, 1) * amp[10] + 1./2. * std::complex<double> (0,
      1) * amp[12] + 1./2. * amp[21] + 1./2. * amp[23] + 1./2. * amp[24] +
      1./2. * amp[25] + 1./2. * std::complex<double> (0, 1) * amp[35] - 1./2. *
      amp[44] + 1./2. * std::complex<double> (0, 1) * amp[48] + 1./2. *
      std::complex<double> (0, 1) * amp[49] + 1./2. * std::complex<double> (0,
      1) * amp[50] - std::complex<double> (0, 1) * amp[62] - 1./2. *
      std::complex<double> (0, 1) * amp[64] - std::complex<double> (0, 1) *
      amp[65] - 1./2. * std::complex<double> (0, 1) * amp[67] - amp[151] -
      1./2. * amp[153] + 1./2. * amp[156] + 1./2. * amp[157] + 1./2. * amp[162]
      + 1./2. * std::complex<double> (0, 1) * amp[180] + 1./2. * amp[181] +
      1./2. * amp[182] + 1./2. * amp[183] + 1./2. * std::complex<double> (0, 1)
      * amp[185] + 1./2. * std::complex<double> (0, 1) * amp[186] + 1./2. *
      std::complex<double> (0, 1) * amp[187] + 1./2. * amp[193] + 1./2. *
      std::complex<double> (0, 1) * amp[194] - 1./2. * amp[198] - amp[196] +
      1./2. * amp[348] - 1./2. * std::complex<double> (0, 1) * amp[349] - 1./2.
      * std::complex<double> (0, 1) * amp[350] + 1./2. * amp[351] - 1./2. *
      std::complex<double> (0, 1) * amp[352] - amp[356] - 1./2. * amp[358] +
      1./2. * amp[376] + 1./2. * amp[377] + 1./2. * amp[379] + 1./2. * amp[380]
      + 1./2. * amp[382] + 1./2. * amp[383] + 1./2. * amp[384] + 1./2. *
      amp[388] + 1./2. * amp[389] + 1./2. * amp[391] + 1./2. *
      std::complex<double> (0, 1) * amp[397] - 1./2. * std::complex<double> (0,
      1) * amp[399] + 1./2. * amp[401];
  jamp[11] = +std::complex<double> (0, 1) * amp[0] + std::complex<double> (0,
      1) * amp[1] - 1./6. * std::complex<double> (0, 1) * amp[2] - 1./6. *
      std::complex<double> (0, 1) * amp[3] + std::complex<double> (0, 1) *
      amp[4] + std::complex<double> (0, 1) * amp[5] + std::complex<double> (0,
      1) * amp[6] + std::complex<double> (0, 1) * amp[7] + std::complex<double>
      (0, 1) * amp[8] + std::complex<double> (0, 1) * amp[9] - 1./6. *
      std::complex<double> (0, 1) * amp[10] - 1./6. * std::complex<double> (0,
      1) * amp[11] - 1./6. * std::complex<double> (0, 1) * amp[12] - 1./6. *
      std::complex<double> (0, 1) * amp[13] + std::complex<double> (0, 1) *
      amp[14] + std::complex<double> (0, 1) * amp[15] + std::complex<double>
      (0, 1) * amp[16] + std::complex<double> (0, 1) * amp[17] +
      std::complex<double> (0, 1) * amp[18] + std::complex<double> (0, 1) *
      amp[19] + std::complex<double> (0, 1) * amp[22] - std::complex<double>
      (0, 1) * amp[60] + 1./6. * std::complex<double> (0, 1) * amp[61] -
      std::complex<double> (0, 1) * amp[63] + 1./6. * std::complex<double> (0,
      1) * amp[64] - amp[277] + 1./6. * amp[278] + amp[279] + amp[280] - 1./6.
      * amp[281] - 1./6. * amp[282] + amp[283] + amp[284] + amp[285] + amp[286]
      + amp[287] + amp[288] - 1./6. * amp[289] - 1./6. * amp[292] + amp[294] +
      amp[295] + amp[296] + 1./6. * amp[300] - amp[299] + amp[345] + amp[346] +
      amp[347] - 1./6. * amp[348] - 1./6. * amp[351] + amp[353] + amp[354] +
      amp[355] - amp[357] + 1./6. * amp[358] + std::complex<double> (0, 1) *
      amp[373] + std::complex<double> (0, 1) * amp[374] - std::complex<double>
      (0, 1) * amp[375];
  jamp[12] = +2. * amp[20] + 2. * amp[26] + 2. * amp[27] - 1./3. * amp[28] -
      1./3. * amp[29] + 2. * amp[30] + 2. * amp[31];
  jamp[13] = +amp[28] + amp[29]; 

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



