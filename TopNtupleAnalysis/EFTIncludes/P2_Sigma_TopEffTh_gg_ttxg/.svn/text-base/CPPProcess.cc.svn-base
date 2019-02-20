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
// Process: g g > t t~ g NP<=2 @2

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
  jamp2[0] = new double[11]; 
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
  for(int i = 0; i < 11; i++ )
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
        t[0] = matrix_2_gg_ttxg(); 

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
      t[0] = matrix_2_gg_ttxg(); 

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
  VVV2P0_1(w[0], w[1], pars->GC_28, pars->ZERO, pars->ZERO, w[5]); 
  FFV5P0_3(w[3], w[2], pars->GC_7, pars->ZERO, pars->ZERO, w[6]); 
  VVV1P0_1(w[0], w[1], pars->GC_6, pars->ZERO, pars->ZERO, w[7]); 
  FFS2_3(w[3], w[2], pars->GC_95, pars->mdl_MH, pars->mdl_WH, w[8]); 
  FFV3_8P0_3(w[3], w[2], pars->GC_145, pars->GC_79, pars->ZERO, pars->ZERO,
      w[9]);
  VVS2_3(w[0], w[1], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[10]); 
  FFV5_1(w[2], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[11]); 
  FFV3_8_1(w[2], w[4], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[12]);
  FFV5_2(w[3], w[4], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[13]); 
  FFV3_8_2(w[3], w[4], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[14]);
  FFV5_1(w[2], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[15]); 
  FFV5_2(w[3], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[16]); 
  FFV3_8_2(w[3], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[17]);
  FFV3_8_1(w[2], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[18]);
  VVS2_3(w[1], w[4], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[19]); 
  VVV2P0_1(w[1], w[4], pars->GC_28, pars->ZERO, pars->ZERO, w[20]); 
  VVV1P0_1(w[1], w[4], pars->GC_6, pars->ZERO, pars->ZERO, w[21]); 
  FFV5_2(w[3], w[0], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[22]); 
  FFV5_1(w[2], w[1], pars->GC_7, pars->mdl_MT, pars->mdl_WT, w[23]); 
  FFV3_8_1(w[2], w[1], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[24]);
  FFV3_8_2(w[3], w[0], pars->GC_145, pars->GC_79, pars->mdl_MT, pars->mdl_WT,
      w[25]);
  VVS2_3(w[0], w[4], pars->GC_78, pars->mdl_MH, pars->mdl_WH, w[26]); 
  VVV2P0_1(w[0], w[4], pars->GC_28, pars->ZERO, pars->ZERO, w[27]); 
  VVV1P0_1(w[0], w[4], pars->GC_6, pars->ZERO, pars->ZERO, w[28]); 
  FFVV1_2_1(w[2], w[0], w[1], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[29]);
  FFVV1_2_2(w[3], w[0], w[1], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[30]);
  VVVV2P0_1(w[0], w[1], w[4], pars->GC_44, pars->ZERO, pars->ZERO, w[31]); 
  VVVV8P0_1(w[0], w[1], w[4], pars->GC_44, pars->ZERO, pars->ZERO, w[32]); 
  VVVV7P0_1(w[0], w[1], w[4], pars->GC_44, pars->ZERO, pars->ZERO, w[33]); 
  VVVV1P0_1(w[0], w[1], w[4], pars->GC_8, pars->ZERO, pars->ZERO, w[34]); 
  VVVV4P0_1(w[0], w[1], w[4], pars->GC_8, pars->ZERO, pars->ZERO, w[35]); 
  VVVV5P0_1(w[0], w[1], w[4], pars->GC_8, pars->ZERO, pars->ZERO, w[36]); 
  VVVS1_4(w[0], w[1], w[4], pars->GC_84, pars->mdl_MH, pars->mdl_WH, w[37]); 
  FFVV1_2P0_3(w[3], w[2], w[0], pars->GC_146, pars->GC_85, pars->ZERO,
      pars->ZERO, w[38]);
  FFVV1_2_1(w[2], w[0], w[4], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[39]);
  FFVV1_2_2(w[3], w[0], w[4], pars->GC_146, pars->GC_85, pars->mdl_MT,
      pars->mdl_WT, w[40]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVV1_0(w[5], w[6], w[4], pars->GC_6, amp[0]); 
  VVS2_0(w[7], w[4], w[8], pars->GC_78, amp[1]); 
  VVV2_0(w[7], w[6], w[4], pars->GC_28, amp[2]); 
  VVV1_0(w[7], w[6], w[4], pars->GC_6, amp[3]); 
  VVV1_0(w[7], w[9], w[4], pars->GC_6, amp[4]); 
  FFS2_0(w[3], w[11], w[10], pars->GC_95, amp[5]); 
  FFV5_0(w[3], w[11], w[5], pars->GC_7, amp[6]); 
  FFV5_0(w[3], w[11], w[7], pars->GC_7, amp[7]); 
  FFV3_8_0(w[3], w[11], w[7], pars->GC_145, pars->GC_79, amp[8]); 
  FFV5_0(w[3], w[12], w[7], pars->GC_7, amp[9]); 
  FFS2_0(w[13], w[2], w[10], pars->GC_95, amp[10]); 
  FFV5_0(w[13], w[2], w[5], pars->GC_7, amp[11]); 
  FFV5_0(w[13], w[2], w[7], pars->GC_7, amp[12]); 
  FFV3_8_0(w[13], w[2], w[7], pars->GC_145, pars->GC_79, amp[13]); 
  FFV5_0(w[14], w[2], w[7], pars->GC_7, amp[14]); 
  FFVV1_2_0(w[3], w[2], w[4], w[7], pars->GC_146, pars->GC_85, amp[15]); 
  FFV5_0(w[16], w[15], w[4], pars->GC_7, amp[16]); 
  FFV3_8_0(w[16], w[15], w[4], pars->GC_145, pars->GC_79, amp[17]); 
  FFV5_0(w[17], w[15], w[4], pars->GC_7, amp[18]); 
  FFV5_0(w[16], w[18], w[4], pars->GC_7, amp[19]); 
  FFS2_0(w[3], w[15], w[19], pars->GC_95, amp[20]); 
  FFV5_0(w[3], w[15], w[20], pars->GC_7, amp[21]); 
  FFV5_0(w[3], w[15], w[21], pars->GC_7, amp[22]); 
  FFV3_8_0(w[3], w[15], w[21], pars->GC_145, pars->GC_79, amp[23]); 
  FFV5_0(w[3], w[18], w[21], pars->GC_7, amp[24]); 
  FFV5_0(w[13], w[15], w[1], pars->GC_7, amp[25]); 
  FFV3_8_0(w[13], w[15], w[1], pars->GC_145, pars->GC_79, amp[26]); 
  FFV5_0(w[14], w[15], w[1], pars->GC_7, amp[27]); 
  FFV5_0(w[13], w[18], w[1], pars->GC_7, amp[28]); 
  FFVV1_2_0(w[3], w[15], w[1], w[4], pars->GC_146, pars->GC_85, amp[29]); 
  FFV5_0(w[22], w[23], w[4], pars->GC_7, amp[30]); 
  FFV3_8_0(w[22], w[23], w[4], pars->GC_145, pars->GC_79, amp[31]); 
  FFV5_0(w[22], w[24], w[4], pars->GC_7, amp[32]); 
  FFV5_0(w[25], w[23], w[4], pars->GC_7, amp[33]); 
  FFS2_0(w[22], w[2], w[19], pars->GC_95, amp[34]); 
  FFV5_0(w[22], w[2], w[20], pars->GC_7, amp[35]); 
  FFV5_0(w[22], w[2], w[21], pars->GC_7, amp[36]); 
  FFV3_8_0(w[22], w[2], w[21], pars->GC_145, pars->GC_79, amp[37]); 
  FFV5_0(w[25], w[2], w[21], pars->GC_7, amp[38]); 
  FFV5_0(w[22], w[11], w[1], pars->GC_7, amp[39]); 
  FFV3_8_0(w[22], w[11], w[1], pars->GC_145, pars->GC_79, amp[40]); 
  FFV5_0(w[22], w[12], w[1], pars->GC_7, amp[41]); 
  FFV5_0(w[25], w[11], w[1], pars->GC_7, amp[42]); 
  FFVV1_2_0(w[22], w[2], w[1], w[4], pars->GC_146, pars->GC_85, amp[43]); 
  FFS2_0(w[3], w[23], w[26], pars->GC_95, amp[44]); 
  FFV5_0(w[3], w[23], w[27], pars->GC_7, amp[45]); 
  FFV5_0(w[3], w[23], w[28], pars->GC_7, amp[46]); 
  FFV3_8_0(w[3], w[23], w[28], pars->GC_145, pars->GC_79, amp[47]); 
  FFV5_0(w[3], w[24], w[28], pars->GC_7, amp[48]); 
  FFS2_0(w[16], w[2], w[26], pars->GC_95, amp[49]); 
  FFV5_0(w[16], w[2], w[27], pars->GC_7, amp[50]); 
  FFV5_0(w[16], w[2], w[28], pars->GC_7, amp[51]); 
  FFV3_8_0(w[16], w[2], w[28], pars->GC_145, pars->GC_79, amp[52]); 
  FFV5_0(w[17], w[2], w[28], pars->GC_7, amp[53]); 
  VVV1_0(w[27], w[1], w[6], pars->GC_6, amp[54]); 
  VVS2_0(w[28], w[1], w[8], pars->GC_78, amp[55]); 
  VVV2_0(w[28], w[1], w[6], pars->GC_28, amp[56]); 
  VVV1_0(w[28], w[1], w[6], pars->GC_6, amp[57]); 
  VVV1_0(w[28], w[1], w[9], pars->GC_6, amp[58]); 
  FFVV1_2_0(w[3], w[2], w[1], w[28], pars->GC_146, pars->GC_85, amp[59]); 
  FFV5_0(w[13], w[23], w[0], pars->GC_7, amp[60]); 
  FFV3_8_0(w[13], w[23], w[0], pars->GC_145, pars->GC_79, amp[61]); 
  FFV5_0(w[14], w[23], w[0], pars->GC_7, amp[62]); 
  FFV5_0(w[13], w[24], w[0], pars->GC_7, amp[63]); 
  FFV5_0(w[16], w[11], w[0], pars->GC_7, amp[64]); 
  FFV3_8_0(w[16], w[11], w[0], pars->GC_145, pars->GC_79, amp[65]); 
  FFV5_0(w[16], w[12], w[0], pars->GC_7, amp[66]); 
  FFV5_0(w[17], w[11], w[0], pars->GC_7, amp[67]); 
  VVV1_0(w[0], w[20], w[6], pars->GC_6, amp[68]); 
  VVS2_0(w[0], w[21], w[8], pars->GC_78, amp[69]); 
  VVV2_0(w[0], w[21], w[6], pars->GC_28, amp[70]); 
  VVV1_0(w[0], w[21], w[6], pars->GC_6, amp[71]); 
  VVV1_0(w[0], w[21], w[9], pars->GC_6, amp[72]); 
  FFV5_0(w[3], w[29], w[4], pars->GC_7, amp[73]); 
  FFV5_0(w[30], w[2], w[4], pars->GC_7, amp[74]); 
  FFV5_0(w[3], w[2], w[31], pars->GC_7, amp[75]); 
  FFV5_0(w[3], w[2], w[32], pars->GC_7, amp[76]); 
  FFV5_0(w[3], w[2], w[33], pars->GC_7, amp[77]); 
  FFV5_0(w[3], w[2], w[34], pars->GC_7, amp[78]); 
  FFV5_0(w[3], w[2], w[35], pars->GC_7, amp[79]); 
  FFV5_0(w[3], w[2], w[36], pars->GC_7, amp[80]); 
  FFV3_8_0(w[3], w[2], w[34], pars->GC_145, pars->GC_79, amp[81]); 
  FFV3_8_0(w[3], w[2], w[35], pars->GC_145, pars->GC_79, amp[82]); 
  FFV3_8_0(w[3], w[2], w[36], pars->GC_145, pars->GC_79, amp[83]); 
  FFS2_0(w[3], w[2], w[37], pars->GC_95, amp[84]); 
  VVV1_0(w[1], w[4], w[38], pars->GC_6, amp[85]); 
  FFV5_0(w[3], w[39], w[1], pars->GC_7, amp[86]); 
  FFV5_0(w[40], w[2], w[1], pars->GC_7, amp[87]); 

}
double CPPProcess::matrix_2_gg_ttxg() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 88; 
  const int ncolor = 11; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {9, 3, 9, 9, 3, 9, 3, 3, 9, 9, 3}; 
  static const double cf[ncolor][ncolor] = {{64, 24, -8, -8, -3, 1, 21, -6, 1,
      10, 24}, {8, 24, 8, -1, 3, 8, 0, 0, -1, 8, 3}, {-8, 24, 64, 1, 24, 10,
      -6, 21, -8, 1, -3}, {-8, -3, 1, 64, 24, -8, -6, 21, 10, 1, 24}, {-1, 3,
      8, 8, 24, 8, 0, 0, 8, -1, 3}, {1, 24, 10, -8, 24, 64, 21, -6, 1, -8, -3},
      {7, 0, -2, -2, 0, 7, 21, -6, 7, -2, 0}, {-2, 0, 7, 7, 0, -2, -6, 21, -2,
      7, 0}, {1, -3, -8, 10, 24, 1, 21, -6, 64, -8, 24}, {10, 24, 1, 1, -3, -8,
      -6, 21, -8, 64, 24}, {8, 3, -1, 8, 3, -1, 0, 0, 8, 8, 24}};

  // Calculate color flows
  jamp[0] = -amp[0] - amp[2] - amp[3] - amp[4] + std::complex<double> (0, 1) *
      amp[11] + std::complex<double> (0, 1) * amp[12] + std::complex<double>
      (0, 1) * amp[13] + std::complex<double> (0, 1) * amp[14] - amp[15] +
      std::complex<double> (0, 1) * amp[21] + std::complex<double> (0, 1) *
      amp[22] + std::complex<double> (0, 1) * amp[23] + std::complex<double>
      (0, 1) * amp[24] - amp[25] - amp[26] - amp[27] - amp[28] +
      std::complex<double> (0, 1) * amp[29] + amp[68] + amp[70] + amp[71] +
      amp[72] + std::complex<double> (0, 1) * amp[73] - amp[77] + amp[75] -
      amp[80] + amp[78] - amp[83] + amp[81] - amp[85];
  jamp[1] = +2. * (-amp[20] - amp[34]); 
  jamp[2] = -amp[16] - amp[17] - amp[18] - amp[19] - std::complex<double> (0,
      1) * amp[21] - std::complex<double> (0, 1) * amp[22] -
      std::complex<double> (0, 1) * amp[23] - std::complex<double> (0, 1) *
      amp[24] - std::complex<double> (0, 1) * amp[29] + std::complex<double>
      (0, 1) * amp[50] + std::complex<double> (0, 1) * amp[51] +
      std::complex<double> (0, 1) * amp[52] + std::complex<double> (0, 1) *
      amp[53] + amp[54] + amp[56] + amp[57] + amp[58] - amp[59] - amp[68] -
      amp[70] - amp[71] - amp[72] - amp[76] - amp[75] - amp[79] - amp[78] -
      amp[82] - amp[81] + amp[85] + std::complex<double> (0, 1) * amp[86];
  jamp[3] = +amp[0] + amp[2] + amp[3] + amp[4] - std::complex<double> (0, 1) *
      amp[11] - std::complex<double> (0, 1) * amp[12] - std::complex<double>
      (0, 1) * amp[13] - std::complex<double> (0, 1) * amp[14] + amp[15] +
      std::complex<double> (0, 1) * amp[45] + std::complex<double> (0, 1) *
      amp[46] + std::complex<double> (0, 1) * amp[47] + std::complex<double>
      (0, 1) * amp[48] - amp[54] - amp[56] - amp[57] - amp[58] + amp[59] -
      amp[60] - amp[61] - amp[62] - amp[63] - std::complex<double> (0, 1) *
      amp[73] + amp[77] + amp[76] + amp[80] + amp[79] + amp[83] + amp[82] +
      std::complex<double> (0, 1) * amp[87];
  jamp[4] = +2. * (-amp[44] - amp[49]); 
  jamp[5] = -amp[30] - amp[31] - amp[32] - amp[33] + std::complex<double> (0,
      1) * amp[35] + std::complex<double> (0, 1) * amp[36] +
      std::complex<double> (0, 1) * amp[37] + std::complex<double> (0, 1) *
      amp[38] + std::complex<double> (0, 1) * amp[43] - std::complex<double>
      (0, 1) * amp[45] - std::complex<double> (0, 1) * amp[46] -
      std::complex<double> (0, 1) * amp[47] - std::complex<double> (0, 1) *
      amp[48] + amp[54] + amp[56] + amp[57] + amp[58] - amp[59] - amp[68] -
      amp[70] - amp[71] - amp[72] - amp[76] - amp[75] - amp[79] - amp[78] -
      amp[82] - amp[81] + amp[85] - std::complex<double> (0, 1) * amp[87];
  jamp[6] = +2. * (+std::complex<double> (0, 1) * amp[1] - std::complex<double>
      (0, 1) * amp[55] + std::complex<double> (0, 1) * amp[69] +
      std::complex<double> (0, 1) * amp[84]);
  jamp[7] = +2. * (-std::complex<double> (0, 1) * amp[1] + std::complex<double>
      (0, 1) * amp[55] - std::complex<double> (0, 1) * amp[69] -
      std::complex<double> (0, 1) * amp[84]);
  jamp[8] = +amp[0] + amp[2] + amp[3] + amp[4] + std::complex<double> (0, 1) *
      amp[6] + std::complex<double> (0, 1) * amp[7] + std::complex<double> (0,
      1) * amp[8] + std::complex<double> (0, 1) * amp[9] + amp[15] -
      std::complex<double> (0, 1) * amp[50] - std::complex<double> (0, 1) *
      amp[51] - std::complex<double> (0, 1) * amp[52] - std::complex<double>
      (0, 1) * amp[53] - amp[54] - amp[56] - amp[57] - amp[58] + amp[59] -
      amp[64] - amp[65] - amp[66] - amp[67] + std::complex<double> (0, 1) *
      amp[74] + amp[77] + amp[76] + amp[80] + amp[79] + amp[83] + amp[82] -
      std::complex<double> (0, 1) * amp[86];
  jamp[9] = -amp[0] - amp[2] - amp[3] - amp[4] - std::complex<double> (0, 1) *
      amp[6] - std::complex<double> (0, 1) * amp[7] - std::complex<double> (0,
      1) * amp[8] - std::complex<double> (0, 1) * amp[9] - amp[15] -
      std::complex<double> (0, 1) * amp[35] - std::complex<double> (0, 1) *
      amp[36] - std::complex<double> (0, 1) * amp[37] - std::complex<double>
      (0, 1) * amp[38] - amp[39] - amp[40] - amp[41] - amp[42] -
      std::complex<double> (0, 1) * amp[43] + amp[68] + amp[70] + amp[71] +
      amp[72] - std::complex<double> (0, 1) * amp[74] - amp[77] + amp[75] -
      amp[80] + amp[78] - amp[83] + amp[81] - amp[85];
  jamp[10] = +2. * (-amp[5] - amp[10]); 

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



