//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "TopNtupleAnalysis/HelAmps_TopEffTh.h"
#include <complex> 
#include <cmath> 
#include <iostream> 
#include <cstdlib> 
using namespace std; 

namespace MG5_TopEffTh 
{


void txxxxx(double p[4], double tmass, int nhel, int nst, complex<double>
    tc[18])
{
  complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  sqh = pow(0.5, 0.5); 
  sqs = pow(0.5/3, 0.5); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], pow(pt2 + p[3] * p[3], 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 

  ft[4][0] = complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if(nhel >= 0)
  {
    if(pp == 0)
    {
      ep[0] = complex<double> (0, 0); 
      ep[1] = complex<double> (-sqh, 0); 
      ep[2] = complex<double> (0, nst * sqh); 
      ep[3] = complex<double> (0, 0); 
    }
    else
    {
      ep[0] = complex<double> (0, 0); 
      ep[3] = complex<double> (pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = complex<double> (-sqh, 0); 
        ep[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }

  }

  // construct eps-
  if(nhel <= 0)
  {
    if(pp == 0)
    {
      em[0] = complex<double> (0, 0); 
      em[1] = complex<double> (sqh, 0); 
      em[2] = complex<double> (0, nst * sqh); 
      em[3] = complex<double> (0, 0); 
    }
    else
    {
      em[0] = complex<double> (0, 0); 
      em[3] = complex<double> (-pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = complex<double> (sqh, 0); 
        em[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }
  }

  // construct eps0
  if(fabs(nhel) <= 1)
  {
    if(pp == 0)
    {
      e0[0] = complex<double> (0, 0); 
      e0[1] = complex<double> (0, 0); 
      e0[2] = complex<double> (0, 0); 
      e0[3] = complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = complex<double> (pp/tmass, 0); 
      e0[3] = complex<double> (p[3] * emp, 0); 

      if(pt != 0)
      {
        e0[1] = complex<double> (p[1] * emp, 0); 
        e0[2] = complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = complex<double> (0, 0); 
        e0[2] = complex<double> (0, 0); 
      }
    }
  }

  if(nhel == 2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if(nhel == -2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if(tmass == 0)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = 0; 
    }
  }
  else if(tmass != 0)
  {
    if(nhel == 1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if(nhel == 0)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqs * (ep[i] * em[j] + em[i] * ep[j]
         + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if(nhel == -1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      std::cerr <<  "Invalid helicity in txxxxx.\n"; 
      std::exit(1); 
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for(j = 0; j < 4; j++ )
  {
    for(i = 0; i < 4; i++ )
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], pow((pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2)), 0.5)); 
    if (pp == 0.0)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.0), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * pow(2.0 * p[0], 0.5), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = pow(0.5, 0.5); 
  hel = double(nhel); 
  nsvahl = nsv * abs(hel); 
  pt2 = pow(p[1], 2) + pow(p[2], 2); 
  pp = min(p[0], pow(pt2 + pow(p[3], 2), 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = pow(pow(p[1], 2) + pow(p[2], 2), 0.5); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
    if (pp == 0.000)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[abs(ip)]; 
      fo[3] = ip * nsf * sqm[abs(ip)]; 
      fo[4] = im * nsf * sqm[abs(im)]; 
      fo[5] = ip * sqm[abs(im)]; 
    }
    else
    {
      pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.00), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * pow(2.0 * p[0], 0.5); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void FFS2_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP39; 
  complex<double> TMP43; 
  TMP43 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  TMP39 = (F1[2] * F2[2] + F1[3] * F2[3]); 
  vertex = COUP * - S3[2] * (+cI * (TMP39 + TMP43)); 
}


void FFV3_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P3[4]; 
  complex<double> denom; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - 2. * cI * M1 * (F2[3] * (P3[0] * (V3[3] + cI * (V3[4])) +
      (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5]))
      + P3[3] * (V3[3] + cI * (V3[4]))))) + F2[2] * (V3[5] * P3[0] - cI *
      (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]));
  F1[3] = denom * - 2. * cI * M1 * (F2[2] * (P3[0] * (V3[3] - cI * (V3[4])) +
      (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + F2[3] * (V3[2] * P3[3] - cI * (V3[3]
      * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]));
  F1[4] = denom * 2. * cI * (F2[2] * (P1[1] * (P3[0] * (+cI * (V3[4]) - V3[3])
      + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[2]) + cI * (V3[5])) +
      P3[3] * (V3[3] - cI * (V3[4]))))) + (P1[2] * (P3[0] * - 1. * (V3[4] + cI
      * (V3[3])) + (P3[1] * (-cI * (V3[5]) + cI * (V3[2])) + (P3[2] * (V3[2] -
      V3[5]) + P3[3] * (V3[4] + cI * (V3[3]))))) + (P1[0] * (V3[5] * P3[0] - cI
      * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]) + P1[3] *
      (V3[2] * P3[3] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] *
      P3[0])))) + F2[3] * (P1[0] * (P3[0] * (V3[3] + cI * (V3[4])) + (P3[1] * -
      1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) + P3[3] *
      (V3[3] + cI * (V3[4]))))) + (P1[3] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) -
      P3[3] * (V3[3] + cI * (V3[4]))))) + (P1[1] * (V3[5] * P3[0] - cI * (V3[4]
      * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]) + P1[2] * (V3[4] * P3[1]
      - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[3] * P3[2])))));
  F1[5] = denom * 2. * cI * (F2[2] * (P1[0] * (P3[0] * (V3[3] - cI * (V3[4])) +
      (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + (P1[3] * (P3[0] * (V3[3] - cI *
      (V3[4])) + (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (+cI * (V3[4]) - V3[3])))) + (P1[1] * (V3[2] * P3[3] -
      cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]) + P1[2] *
      (V3[4] * P3[1] - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[3] *
      P3[2])))) + F2[3] * (P1[1] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P3[1] * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) - P3[3] *
      (V3[3] + cI * (V3[4]))))) + (P1[2] * (P3[0] * (+cI * (V3[3]) - V3[4]) +
      (P3[1] * - 1. * (+cI * (V3[2] + V3[5])) + (P3[2] * (V3[2] + V3[5]) +
      P3[3] * (+cI * (V3[3]) - V3[4])))) + (P1[0] * (V3[2] * P3[3] - cI *
      (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]) + P1[3] * (V3[2]
      * P3[3] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] *
      P3[0])))));
}

void FFV3_8_1(complex<double> F2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P3[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV3_1(F2, V3, COUP1, M1, W1, F1); 
  FFV8_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFFF2_3(complex<double> F1[], complex<double> F2[], complex<double> F4[],
    complex<double> COUP, double M3, double W3, complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  F3[0] = +F1[0] + F2[0] + F4[0]; 
  F3[1] = +F1[1] + F2[1] + F4[1]; 
  P3[0] = -F3[0].real(); 
  P3[1] = -F3[1].real(); 
  P3[2] = -F3[1].imag(); 
  P3[3] = -F3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  F3[2] = denom * 2. * cI * F1[3] * M3 * (F2[4] * F4[5] - F2[5] * F4[4]); 
  F3[3] = denom * 2. * cI * F1[2] * M3 * (F2[5] * F4[4] - F2[4] * F4[5]); 
  F3[4] = denom * - 2. * cI * (F1[2] * (F2[4] * F4[5] * (P3[1] + cI * (P3[2]))
      - F2[5] * F4[4] * (P3[1] + cI * (P3[2]))) + F1[3] * (F2[4] * F4[5] *
      (P3[0] - P3[3]) + F2[5] * F4[4] * (P3[3] - P3[0])));
  F3[5] = denom * - 2. * cI * (F1[2] * (F2[4] * - F4[5] * (P3[0] + P3[3]) +
      F2[5] * F4[4] * (P3[0] + P3[3])) + F1[3] * (F2[4] * F4[5] * (+cI *
      (P3[2]) - P3[1]) + F2[5] * F4[4] * (P3[1] - cI * (P3[2]))));
}


void FFV1_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * (+cI *
      (V3[4]) - V3[3]))));
  F2[3] = denom * - cI * (F1[2] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] * (V3[2] +
      V3[5]))));
  F2[4] = denom * - cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] +
      cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * - 1. *
      (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] *
      (V3[2] - V3[5]))));
}


void FFFF4_3(complex<double> F1[], complex<double> F2[], complex<double> F4[],
    complex<double> COUP, double M3, double W3, complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  complex<double> TMP43; 
  F3[0] = +F1[0] + F2[0] + F4[0]; 
  F3[1] = +F1[1] + F2[1] + F4[1]; 
  P3[0] = -F3[0].real(); 
  P3[1] = -F3[1].real(); 
  P3[2] = -F3[1].imag(); 
  P3[3] = -F3[0].imag(); 
  TMP43 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  F3[2] = denom * cI * F4[2] * TMP43 * M3; 
  F3[3] = denom * cI * F4[3] * TMP43 * M3; 
  F3[4] = denom * cI * TMP43 * (F4[2] * (P3[3] - P3[0]) + F4[3] * (P3[1] + cI *
      (P3[2])));
  F3[5] = denom * - cI * TMP43 * (F4[2] * (+cI * (P3[2]) - P3[1]) + F4[3] *
      (P3[0] + P3[3]));
}

void FFFF4_5_3(complex<double> F1[], complex<double> F2[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M3, double W3,
    complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  int i; 
  complex<double> Ftmp[6]; 
  FFFF4_3(F1, F2, F4, COUP1, M3, W3, F3); 
  FFFF5_3(F1, F2, F4, COUP2, M3, W3, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F3[i] = F3[i] + Ftmp[i]; 
    i++; 
  }
}

void FFV8_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P3[4]; 
  complex<double> denom; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * 2. * cI * (F2[4] * (P1[1] * (P3[0] * (+cI * (V3[4]) - V3[3])
      + (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + (P1[2] * (P3[0] * - 1. * (V3[4] + cI
      * (V3[3])) + (P3[1] * (+cI * (V3[2] + V3[5])) + (P3[2] * (V3[2] + V3[5])
      - P3[3] * (V3[4] + cI * (V3[3]))))) + (P1[0] * (V3[2] * P3[3] - cI *
      (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]) + P1[3] * (V3[2]
      * P3[3] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]))))
      + F2[5] * (P1[0] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) + (P3[1] *
      (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) + P3[3] *
      (V3[3] + cI * (V3[4]))))) + (P1[3] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P1[1] * (V3[5] * P3[0] -
      cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]) + P1[2] *
      (V3[3] * P3[2] - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[4] *
      P3[1])))));
  F1[3] = denom * 2. * cI * (F2[4] * (P1[0] * (P3[0] * (+cI * (V3[4]) - V3[3])
      + (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + (P1[3] * (P3[0] * (V3[3] - cI *
      (V3[4])) + (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] +
      V3[5])) + P3[3] * (V3[3] - cI * (V3[4]))))) + (P1[1] * (V3[2] * P3[3] -
      cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]) + P1[2] *
      (V3[3] * P3[2] - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[4] *
      P3[1])))) + F2[5] * (P1[1] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (V3[3] + cI * (V3[4]))))) + (P1[2] * (P3[0] * (+cI * (V3[3]) -
      V3[4]) + (P3[1] * (-cI * (V3[2]) + cI * (V3[5])) + (P3[2] * (V3[2] -
      V3[5]) + P3[3] * (V3[4] - cI * (V3[3]))))) + (P1[0] * (V3[5] * P3[0] - cI
      * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]) + P1[3] *
      (V3[2] * P3[3] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] *
      P3[0])))));
  F1[4] = denom * - 2. * cI * M1 * (F2[5] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + F2[4] * (V3[2] * P3[3] -
      cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]));
  F1[5] = denom * - 2. * cI * M1 * (F2[4] * (P3[0] * (+cI * (V3[4]) - V3[3]) +
      (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + F2[5] * (V3[5] * P3[0] - cI * (V3[3]
      * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]));
}


void VVV1_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP31; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> TMP33; 
  complex<double> TMP27; 
  complex<double> TMP29; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP33 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  vertex = COUP * (TMP10 * (-cI * (TMP37) + cI * (TMP33)) + (TMP12 * (-cI *
      (TMP29) + cI * (TMP30)) + TMP23 * (-cI * (TMP27) + cI * (TMP31))));
}


void FFFF1_0(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> F4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP51; 
  TMP51 = (F1[2] * F3[3] * (F2[5] * F4[4] - F2[4] * F4[5]) + F1[3] * F3[2] *
      (F2[4] * F4[5] - F2[5] * F4[4]));
  vertex = COUP * - 2. * cI * TMP51; 
}

void FFFF1_2_6_0(complex<double> F1[], complex<double> F2[], complex<double>
    F3[], complex<double> F4[], complex<double> COUP1, complex<double> COUP2,
    complex<double> COUP3, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFFF1_0(F1, F2, F3, F4, COUP1, vertex); 
  FFFF2_0(F1, F2, F3, F4, COUP2, tmp); 
  vertex = vertex + tmp; 
  FFFF6_0(F1, F2, F3, F4, COUP3, tmp); 
  vertex = vertex + tmp; 
}
void FFFF1_6_0(complex<double> F1[], complex<double> F2[], complex<double>
    F3[], complex<double> F4[], complex<double> COUP1, complex<double> COUP2,
    complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  complex<double> COUP3; 
  FFFF1_0(F1, F2, F3, F4, COUP1, vertex); 
  FFFF6_0(F1, F2, F3, F4, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV4_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP14; 
  complex<double> TMP13; 
  TMP14 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP13 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * (-cI * (TMP13) + 2. * cI * (TMP14)); 
}


void FFVV1_1(complex<double> F2[], complex<double> V3[], complex<double> V4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0] + V4[0]; 
  F1[1] = +F2[1] + V3[1] + V4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * 2. * cI * M1 * (F2[3] * (V3[2] * (V4[3] + cI * (V4[4])) +
      (V3[3] * - 1. * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5]))
      + V3[5] * (V4[3] + cI * (V4[4]))))) + F2[2] * (V3[2] * V4[5] - cI *
      (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[5] * V4[2]));
  F1[3] = denom * 2. * cI * M1 * (F2[2] * (V3[2] * (V4[3] - cI * (V4[4])) +
      (V3[3] * (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) +
      V3[5] * (+cI * (V4[4]) - V4[3])))) + F2[3] * (V3[5] * V4[2] - cI * (V3[4]
      * V4[3]) + cI * (V3[3] * V4[4]) - V3[2] * V4[5]));
  F1[4] = denom * - 2. * cI * (F2[2] * (P1[1] * (V3[2] * (+cI * (V4[4]) -
      V4[3]) + (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[2]) + cI *
      (V4[5])) + V3[5] * (V4[3] - cI * (V4[4]))))) + (P1[2] * (V3[2] * - 1. *
      (V4[4] + cI * (V4[3])) + (V3[3] * (-cI * (V4[5]) + cI * (V4[2])) + (V3[4]
      * (V4[2] - V4[5]) + V3[5] * (V4[4] + cI * (V4[3]))))) + (P1[0] * (V3[2] *
      V4[5] - cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[5] * V4[2]) +
      P1[3] * (V3[5] * V4[2] - cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) -
      V3[2] * V4[5])))) + F2[3] * (P1[0] * (V3[2] * (V4[3] + cI * (V4[4])) +
      (V3[3] * - 1. * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5]))
      + V3[5] * (V4[3] + cI * (V4[4]))))) + (P1[3] * (V3[2] * - 1. * (V4[3] +
      cI * (V4[4])) + (V3[3] * (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] +
      V4[5])) - V3[5] * (V4[3] + cI * (V4[4]))))) + (P1[1] * (V3[2] * V4[5] -
      cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[5] * V4[2]) + P1[2] *
      (V3[3] * V4[4] - cI * (V3[5] * V4[2]) + cI * (V3[2] * V4[5]) - V3[4] *
      V4[3])))));
  F1[5] = denom * - 2. * cI * (F2[2] * (P1[0] * (V3[2] * (V4[3] - cI * (V4[4]))
      + (V3[3] * (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) +
      V3[5] * (+cI * (V4[4]) - V4[3])))) + (P1[3] * (V3[2] * (V4[3] - cI *
      (V4[4])) + (V3[3] * (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI *
      (V4[2])) + V3[5] * (+cI * (V4[4]) - V4[3])))) + (P1[1] * (V3[5] * V4[2] -
      cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[2] * V4[5]) + P1[2] *
      (V3[3] * V4[4] - cI * (V3[5] * V4[2]) + cI * (V3[2] * V4[5]) - V3[4] *
      V4[3])))) + F2[3] * (P1[1] * (V3[2] * - 1. * (V4[3] + cI * (V4[4])) +
      (V3[3] * (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] + V4[5])) - V3[5] *
      (V4[3] + cI * (V4[4]))))) + (P1[2] * (V3[2] * (+cI * (V4[3]) - V4[4]) +
      (V3[3] * - 1. * (+cI * (V4[2] + V4[5])) + (V3[4] * (V4[2] + V4[5]) +
      V3[5] * (+cI * (V4[3]) - V4[4])))) + (P1[0] * (V3[5] * V4[2] - cI *
      (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[2] * V4[5]) + P1[3] * (V3[5]
      * V4[2] - cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[2] *
      V4[5])))));
}

void FFVV1_2_1(complex<double> F2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP1, complex<double> COUP2, double M1, double W1,
    complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFVV1_1(F2, V3, V4, COUP1, M1, W1, F1); 
  FFVV2_1(F2, V3, V4, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFFF3_2(complex<double> F1[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + F3[0] + F4[0]; 
  F2[1] = +F1[1] + F3[1] + F4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * 2. * cI * (F4[4] * (F1[2] * F3[3] * (P2[1] - cI * (P2[2])) +
      F1[3] * F3[2] * (+cI * (P2[2]) - P2[1])) + F4[5] * (F1[2] * F3[3] *
      (P2[0] - P2[3]) + F1[3] * F3[2] * (P2[3] - P2[0])));
  F2[3] = denom * 2. * cI * (F4[4] * (F1[2] * - F3[3] * (P2[0] + P2[3]) + F1[3]
      * F3[2] * (P2[0] + P2[3])) + F4[5] * (F1[2] * - F3[3] * (P2[1] + cI *
      (P2[2])) + F1[3] * F3[2] * (P2[1] + cI * (P2[2]))));
  F2[4] = denom * 2. * cI * F4[5] * M2 * (F1[2] * F3[3] - F1[3] * F3[2]); 
  F2[5] = denom * 2. * cI * F4[4] * M2 * (F1[3] * F3[2] - F1[2] * F3[3]); 
}

void FFFF3_7_2(complex<double> F1[], complex<double> F3[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M2, double W2,
    complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  double P2[4]; 
  complex<double> denom; 
  int i; 
  FFFF3_2(F1, F3, F4, COUP1, M2, W2, F2); 
  FFFF7_2(F1, F3, F4, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void VVVVV14P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP17; 
  complex<double> TMP20; 
  complex<double> TMP16; 
  complex<double> TMP21; 
  complex<double> denom; 
  double P5[4]; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P5[0] = V5[0].real(); 
  P5[1] = V5[1].real(); 
  P5[2] = V5[1].imag(); 
  P5[3] = V5[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP20 = (V2[2] * P5[0] - V2[3] * P5[1] - V2[4] * P5[2] - V2[5] * P5[3]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP17 = (V3[2] * P5[0] - V3[3] * P5[1] - V3[4] * P5[2] - V3[5] * P5[3]); 
  TMP16 = (V4[2] * P5[0] - V4[3] * P5[1] - V4[4] * P5[2] - V4[5] * P5[3]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P5[0] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP16 * (-cI * (TMP10 * V5[2]) + cI * (V2[2] * TMP19)) + (TMP17 * (-cI *
      (V2[2] * TMP21) + cI * (V4[2] * TMP18)) + TMP20 * (-cI * (V4[2] * TMP19)
      + cI * (V5[2] * TMP22)))));
  V1[3] = denom * (P5[1] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP16 * (-cI * (TMP10 * V5[3]) + cI * (V2[3] * TMP19)) + (TMP17 * (-cI *
      (V2[3] * TMP21) + cI * (V4[3] * TMP18)) + TMP20 * (-cI * (V4[3] * TMP19)
      + cI * (V5[3] * TMP22)))));
  V1[4] = denom * (P5[2] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP16 * (-cI * (TMP10 * V5[4]) + cI * (V2[4] * TMP19)) + (TMP17 * (-cI *
      (V2[4] * TMP21) + cI * (V4[4] * TMP18)) + TMP20 * (-cI * (V4[4] * TMP19)
      + cI * (V5[4] * TMP22)))));
  V1[5] = denom * (P5[3] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP16 * (-cI * (TMP10 * V5[5]) + cI * (V2[5] * TMP19)) + (TMP17 * (-cI *
      (V2[5] * TMP21) + cI * (V4[5] * TMP18)) + TMP20 * (-cI * (V4[5] * TMP19)
      + cI * (V5[5] * TMP22)))));
}


void FFFF5_2(complex<double> F1[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP38; 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + F3[0] + F4[0]; 
  F2[1] = +F1[1] + F3[1] + F4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  TMP38 = (F4[4] * F3[4] + F4[5] * F3[5]); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * TMP38 * F1[2] * M2; 
  F2[3] = denom * cI * TMP38 * F1[3] * M2; 
  F2[4] = denom * - cI * TMP38 * (F1[2] * - 1. * (P2[0] + P2[3]) + F1[3] * (+cI
      * (P2[2]) - P2[1]));
  F2[5] = denom * cI * TMP38 * (F1[2] * (P2[1] + cI * (P2[2])) + F1[3] * (P2[0]
      - P2[3]));
}


void FFVV1_2(complex<double> F1[], complex<double> V3[], complex<double> V4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0] + V4[0]; 
  F2[1] = +F1[1] + V3[1] + V4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * 2. * cI * M2 * (F1[3] * (V3[2] * (V4[3] - cI * (V4[4])) +
      (V3[3] * (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) +
      V3[5] * (+cI * (V4[4]) - V4[3])))) + F1[2] * (V3[2] * V4[5] - cI * (V3[3]
      * V4[4]) + cI * (V3[4] * V4[3]) - V3[5] * V4[2]));
  F2[3] = denom * 2. * cI * M2 * (F1[2] * (V3[2] * (V4[3] + cI * (V4[4])) +
      (V3[3] * - 1. * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5]))
      + V3[5] * (V4[3] + cI * (V4[4]))))) + F1[3] * (V3[5] * V4[2] - cI *
      (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[2] * V4[5]));
  F2[4] = denom * 2. * cI * (F1[2] * (P2[1] * (V3[2] * (V4[3] + cI * (V4[4])) +
      (V3[3] * - 1. * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5]))
      + V3[5] * (V4[3] + cI * (V4[4]))))) + (P2[2] * (V3[2] * (V4[4] - cI *
      (V4[3])) + (V3[3] * (+cI * (V4[2] + V4[5])) + (V3[4] * - 1. * (V4[2] +
      V4[5]) + V3[5] * (V4[4] - cI * (V4[3]))))) + (P2[0] * (V3[2] * V4[5] - cI
      * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[5] * V4[2]) + P2[3] *
      (V3[2] * V4[5] - cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[5] *
      V4[2])))) + F1[3] * (P2[0] * (V3[2] * (V4[3] - cI * (V4[4])) + (V3[3] *
      (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) + V3[5] * (+cI
      * (V4[4]) - V4[3])))) + (P2[3] * (V3[2] * (V4[3] - cI * (V4[4])) + (V3[3]
      * (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) + V3[5] *
      (+cI * (V4[4]) - V4[3])))) + (P2[1] * (V3[5] * V4[2] - cI * (V3[4] *
      V4[3]) + cI * (V3[3] * V4[4]) - V3[2] * V4[5]) + P2[2] * (V3[3] * V4[4] -
      cI * (V3[5] * V4[2]) + cI * (V3[2] * V4[5]) - V3[4] * V4[3])))));
  F2[5] = denom * 2. * cI * (F1[2] * (P2[0] * (V3[2] * (V4[3] + cI * (V4[4])) +
      (V3[3] * - 1. * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5]))
      + V3[5] * (V4[3] + cI * (V4[4]))))) + (P2[3] * (V3[2] * - 1. * (V4[3] +
      cI * (V4[4])) + (V3[3] * (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] +
      V4[5])) - V3[5] * (V4[3] + cI * (V4[4]))))) + (P2[1] * (V3[2] * V4[5] -
      cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[5] * V4[2]) + P2[2] *
      (V3[3] * V4[4] - cI * (V3[5] * V4[2]) + cI * (V3[2] * V4[5]) - V3[4] *
      V4[3])))) + F1[3] * (P2[1] * (V3[2] * (V4[3] - cI * (V4[4])) + (V3[3] *
      (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) + V3[5] * (+cI
      * (V4[4]) - V4[3])))) + (P2[2] * (V3[2] * (V4[4] + cI * (V4[3])) + (V3[3]
      * (-cI * (V4[2]) + cI * (V4[5])) + (V3[4] * (V4[5] - V4[2]) - V3[5] *
      (V4[4] + cI * (V4[3]))))) + (P2[0] * (V3[5] * V4[2] - cI * (V3[4] *
      V4[3]) + cI * (V3[3] * V4[4]) - V3[2] * V4[5]) + P2[3] * (V3[2] * V4[5] -
      cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[5] * V4[2])))));
}

void FFVV1_2_2(complex<double> F1[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP1, complex<double> COUP2, double M2, double W2,
    complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  double P2[4]; 
  complex<double> denom; 
  int i; 
  FFVV1_2(F1, V3, V4, COUP1, M2, W2, F2); 
  FFVV2_2(F1, V3, V4, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void VVVV1_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  complex<double> TMP10; 
  complex<double> TMP9; 
  TMP9 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  vertex = COUP * (-cI * (TMP9 * TMP10) + cI * (TMP11 * TMP12)); 
}


void FFV2_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])));
  F1[3] = denom * - cI * M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] * (V3[5]
      - V3[2]));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * -
      1. * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3]
      - cI * (V3[4]))))) + F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))));
}

void FFV2_4_1(complex<double> F2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV2_1(F2, V3, COUP1, M1, W1, F1); 
  FFV4_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_7_1(complex<double> F2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV2_1(F2, V3, COUP1, M1, W1, F1); 
  FFV7_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFVV2_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP56; 
  complex<double> TMP55; 
  TMP55 = (F1[4] * (F2[4] * (V3[2] * (V4[2] - V4[5]) + (V3[3] * - 1. * (V4[3] +
      cI * (V4[4])) + (V3[4] * (+cI * (V4[3]) - V4[4]) + V3[5] * (V4[2] -
      V4[5])))) + F2[5] * (V3[2] * - 1. * (V4[3] + cI * (V4[4])) + (V3[3] *
      (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) + V3[5] *
      (V4[3] + cI * (V4[4])))))) + F1[5] * (F2[4] * (V3[2] * (+cI * (V4[4]) -
      V4[3]) + (V3[3] * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] +
      V4[5])) + V3[5] * (+cI * (V4[4]) - V4[3])))) + F2[5] * (V3[2] * (V4[2] +
      V4[5]) + (V3[3] * (+cI * (V4[4]) - V4[3]) + (V3[4] * - 1. * (V4[4] + cI *
      (V4[3])) - V3[5] * (V4[2] + V4[5]))))));
  TMP56 = (F1[4] * (F2[4] * (V3[2] * (V4[2] + V4[5]) + (V3[3] * (+cI * (V4[4])
      - V4[3]) + (V3[4] * - 1. * (V4[4] + cI * (V4[3])) - V3[5] * (V4[2] +
      V4[5])))) + F2[5] * (V3[2] * (V4[3] + cI * (V4[4])) + (V3[3] * (V4[5] -
      V4[2]) + (V3[4] * (-cI * (V4[2]) + cI * (V4[5])) - V3[5] * (V4[3] + cI *
      (V4[4])))))) + F1[5] * (F2[4] * (V3[2] * (V4[3] - cI * (V4[4])) + (V3[3]
      * - 1. * (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] + V4[5])) + V3[5] *
      (V4[3] - cI * (V4[4]))))) + F2[5] * (V3[2] * (V4[2] - V4[5]) + (V3[3] * -
      1. * (V4[3] + cI * (V4[4])) + (V3[4] * (+cI * (V4[3]) - V4[4]) + V3[5] *
      (V4[2] - V4[5]))))));
  vertex = COUP * (-cI * (TMP55) + cI * (TMP56)); 
}


void VVVVV2P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP31; 
  complex<double> denom; 
  complex<double> TMP40; 
  complex<double> TMP24; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP40 = (P2[0] * V5[2] - P2[1] * V5[3] - P2[2] * V5[4] - P2[3] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P2[0] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP10 * (-cI * (V4[2] * TMP40) + cI * (V5[2] * TMP24)) + (TMP31 * (-cI *
      (TMP11 * V5[2]) + cI * (V4[2] * TMP18)) + V2[2] * (-cI * (TMP19 * TMP24)
      + cI * (TMP22 * TMP40)))));
  V1[3] = denom * (P2[1] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP10 * (-cI * (V4[3] * TMP40) + cI * (V5[3] * TMP24)) + (TMP31 * (-cI *
      (TMP11 * V5[3]) + cI * (V4[3] * TMP18)) + V2[3] * (-cI * (TMP19 * TMP24)
      + cI * (TMP22 * TMP40)))));
  V1[4] = denom * (P2[2] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP10 * (-cI * (V4[4] * TMP40) + cI * (V5[4] * TMP24)) + (TMP31 * (-cI *
      (TMP11 * V5[4]) + cI * (V4[4] * TMP18)) + V2[4] * (-cI * (TMP19 * TMP24)
      + cI * (TMP22 * TMP40)))));
  V1[5] = denom * (P2[3] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP10 * (-cI * (V4[5] * TMP40) + cI * (V5[5] * TMP24)) + (TMP31 * (-cI *
      (TMP11 * V5[5]) + cI * (V4[5] * TMP18)) + V2[5] * (-cI * (TMP19 * TMP24)
      + cI * (TMP22 * TMP40)))));
}


void FFV1_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI
      * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * - 1. * (V3[2] +
      V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] + cI *
      (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])))));
  F1[3] = denom * - cI * (F2[2] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] *
      - 1. * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] *
      (V3[5] - V3[2]))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (V3[5] - V3[2]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] +
      V3[5]))));
}


void FFV8_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  double P3[4]; 
  complex<double> denom; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * - 2. * cI * (F1[4] * (P2[1] * (P3[0] * (V3[3] + cI * (V3[4]))
      + (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[2]) + cI * (V3[5])) -
      P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[2] * (P3[0] * (V3[4] - cI *
      (V3[3])) + (P3[1] * (-cI * (V3[5]) + cI * (V3[2])) + (P3[2] * (V3[5] -
      V3[2]) + P3[3] * (+cI * (V3[3]) - V3[4])))) + (P2[0] * (V3[2] * P3[3] -
      cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]) + P2[3] *
      (V3[5] * P3[0] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] *
      P3[3])))) + F1[5] * (P2[0] * (P3[0] * (+cI * (V3[4]) - V3[3]) + (P3[1] *
      (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) + P3[3] * (+cI
      * (V3[4]) - V3[3])))) + (P2[3] * (P3[0] * (V3[3] - cI * (V3[4])) + (P3[1]
      * - 1. * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) + P3[3] *
      (V3[3] - cI * (V3[4]))))) + (P2[1] * (V3[2] * P3[3] - cI * (V3[4] *
      P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]) + P2[2] * (V3[3] * P3[2] -
      cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[4] * P3[1])))));
  F2[3] = denom * - 2. * cI * (F1[4] * (P2[0] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[3] * (P3[0] * - 1. *
      (V3[3] + cI * (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI *
      (V3[5]) + cI * (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[1] *
      (V3[5] * P3[0] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] *
      P3[3]) + P2[2] * (V3[3] * P3[2] - cI * (V3[2] * P3[3]) + cI * (V3[5] *
      P3[0]) - V3[4] * P3[1])))) + F1[5] * (P2[1] * (P3[0] * (V3[3] - cI *
      (V3[4])) + (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] +
      V3[5])) + P3[3] * (V3[3] - cI * (V3[4]))))) + (P2[2] * (P3[0] * (V3[4] +
      cI * (V3[3])) + (P3[1] * - 1. * (+cI * (V3[2] + V3[5])) + (P3[2] * - 1. *
      (V3[2] + V3[5]) + P3[3] * (V3[4] + cI * (V3[3]))))) + (P2[0] * (V3[5] *
      P3[0] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]) +
      P2[3] * (V3[5] * P3[0] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) -
      V3[2] * P3[3])))));
  F2[4] = denom * - 2. * cI * M2 * (F1[5] * (P3[0] * (+cI * (V3[4]) - V3[3]) +
      (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + F1[4] * (V3[2] * P3[3] - cI * (V3[4]
      * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]));
  F2[5] = denom * - 2. * cI * M2 * (F1[4] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + F1[5] * (V3[5] * P3[0] -
      cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]));
}


void FFS2_2(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + (F1[5] * (+cI *
      (P2[2]) - P2[1]) + F1[2] * M2));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] * -
      1. * (P2[0] + P2[3]) - F1[3] * M2));
  F2[4] = denom * - cI * S3[2] * (F1[2] * - 1. * (P2[0] + P2[3]) + (F1[3] *
      (+cI * (P2[2]) - P2[1]) - F1[4] * M2));
  F2[5] = denom * cI * S3[2] * (F1[2] * (P2[1] + cI * (P2[2])) + (F1[3] *
      (P2[0] - P2[3]) + F1[5] * M2));
}


void FFV5P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5]
      * F2[3]);
  V3[3] = denom * - cI * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3]
      * F2[4]);
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] *
      F2[4] + F1[4] * F2[3]));
  V3[5] = denom * - cI * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5]
      * F2[3]);
}


void VVS2_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP37; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP30; 
  complex<double> TMP47; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP47 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  vertex = COUP * S3[2] * (-cI * (TMP30 * TMP37) + cI * (TMP23 * TMP47)); 
}


void VVVVV15P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  complex<double> TMP17; 
  complex<double> TMP20; 
  complex<double> TMP16; 
  complex<double> TMP21; 
  complex<double> denom; 
  double P5[4]; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P5[0] = V5[0].real(); 
  P5[1] = V5[1].real(); 
  P5[2] = V5[1].imag(); 
  P5[3] = V5[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP20 = (V2[2] * P5[0] - V2[3] * P5[1] - V2[4] * P5[2] - V2[5] * P5[3]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP17 = (V3[2] * P5[0] - V3[3] * P5[1] - V3[4] * P5[2] - V3[5] * P5[3]); 
  TMP16 = (V4[2] * P5[0] - V4[3] * P5[1] - V4[4] * P5[2] - V4[5] * P5[3]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P5[0] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP16 * (-cI * (TMP10 * V5[2]) + cI * (V3[2] * TMP18)) + (TMP17 * (-cI *
      (V4[2] * TMP18) + cI * (TMP11 * V5[2])) + TMP20 * (-cI * (V3[2] * TMP21)
      + cI * (V4[2] * TMP19)))));
  V1[3] = denom * (P5[1] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP16 * (-cI * (TMP10 * V5[3]) + cI * (V3[3] * TMP18)) + (TMP17 * (-cI *
      (V4[3] * TMP18) + cI * (TMP11 * V5[3])) + TMP20 * (-cI * (V3[3] * TMP21)
      + cI * (V4[3] * TMP19)))));
  V1[4] = denom * (P5[2] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP16 * (-cI * (TMP10 * V5[4]) + cI * (V3[4] * TMP18)) + (TMP17 * (-cI *
      (V4[4] * TMP18) + cI * (TMP11 * V5[4])) + TMP20 * (-cI * (V3[4] * TMP21)
      + cI * (V4[4] * TMP19)))));
  V1[5] = denom * (P5[3] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP16 * (-cI * (TMP10 * V5[5]) + cI * (V3[5] * TMP18)) + (TMP17 * (-cI *
      (V4[5] * TMP18) + cI * (TMP11 * V5[5])) + TMP20 * (-cI * (V3[5] * TMP21)
      + cI * (V4[5] * TMP19)))));
}


void VVVVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP30; 
  complex<double> denom; 
  complex<double> TMP28; 
  complex<double> TMP27; 
  complex<double> TMP18; 
  complex<double> TMP19; 
  complex<double> TMP41; 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP41 = (P1[0] * V5[2] - P1[1] * V5[3] - P1[2] * V5[4] - P1[3] * V5[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP11 * (-cI * (V3[2] * TMP41) + cI * (V5[2] * TMP27)) +
      (TMP18 * (-cI * (V4[2] * TMP27) + cI * (V3[2] * TMP28)) + (TMP19 * (-cI *
      (V2[2] * TMP28) + cI * (V4[2] * TMP30)) + TMP22 * (-cI * (V5[2] * TMP30)
      + cI * (V2[2] * TMP41)))));
  V1[3] = denom * (TMP11 * (-cI * (V3[3] * TMP41) + cI * (V5[3] * TMP27)) +
      (TMP18 * (-cI * (V4[3] * TMP27) + cI * (V3[3] * TMP28)) + (TMP19 * (-cI *
      (V2[3] * TMP28) + cI * (V4[3] * TMP30)) + TMP22 * (-cI * (V5[3] * TMP30)
      + cI * (V2[3] * TMP41)))));
  V1[4] = denom * (TMP11 * (-cI * (V3[4] * TMP41) + cI * (V5[4] * TMP27)) +
      (TMP18 * (-cI * (V4[4] * TMP27) + cI * (V3[4] * TMP28)) + (TMP19 * (-cI *
      (V2[4] * TMP28) + cI * (V4[4] * TMP30)) + TMP22 * (-cI * (V5[4] * TMP30)
      + cI * (V2[4] * TMP41)))));
  V1[5] = denom * (TMP11 * (-cI * (V3[5] * TMP41) + cI * (V5[5] * TMP27)) +
      (TMP18 * (-cI * (V4[5] * TMP27) + cI * (V3[5] * TMP28)) + (TMP19 * (-cI *
      (V2[5] * TMP28) + cI * (V4[5] * TMP30)) + TMP22 * (-cI * (V5[5] * TMP30)
      + cI * (V2[5] * TMP41)))));
}


void FFFF4_1(complex<double> F2[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  complex<double> TMP48; 
  F1[0] = +F2[0] + F3[0] + F4[0]; 
  F1[1] = +F2[1] + F3[1] + F4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  TMP48 = (F4[2] * F3[2] + F4[3] * F3[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * TMP48 * (F2[4] * (P1[0] + P1[3]) + F2[5] * (P1[1] + cI
      * (P1[2])));
  F1[3] = denom * cI * TMP48 * (F2[4] * (+cI * (P1[2]) - P1[1]) + F2[5] *
      (P1[3] - P1[0]));
  F1[4] = denom * cI * TMP48 * F2[4] * M1; 
  F1[5] = denom * cI * TMP48 * F2[5] * M1; 
}

void FFFF4_5_1(complex<double> F2[], complex<double> F3[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M1, double W1,
    complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  int i; 
  complex<double> denom; 
  complex<double> Ftmp[6]; 
  FFFF4_1(F2, F3, F4, COUP1, M1, W1, F1); 
  FFFF5_1(F2, F3, F4, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void VVV2_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP31; 
  double P3[4]; 
  complex<double> TMP46; 
  complex<double> TMP30; 
  complex<double> TMP47; 
  complex<double> TMP33; 
  complex<double> TMP27; 
  complex<double> TMP29; 
  complex<double> TMP35; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP46 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  TMP47 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP33 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP35 = (P3[0] * P1[0] - P3[1] * P1[1] - P3[2] * P1[2] - P3[3] * P1[3]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  vertex = COUP * (TMP10 * (-cI * (TMP35 * TMP37) + cI * (TMP33 * TMP47)) +
      (TMP12 * (-cI * (TMP29 * TMP47) + cI * (TMP30 * TMP46)) + (TMP23 * (-cI *
      (TMP27 * TMP46) + cI * (TMP31 * TMP35)) + (-cI * (TMP30 * TMP31 * TMP33)
      + cI * (TMP27 * TMP29 * TMP37)))));
}


void FFFF1_2(complex<double> F1[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + F3[0] + F4[0]; 
  F2[1] = +F1[1] + F3[1] + F4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * 2. * cI * (F4[4] * (F1[2] * F3[3] * (+cI * (P2[2]) - P2[1]) +
      F1[3] * F3[2] * (P2[1] - cI * (P2[2]))) + F4[5] * (F1[2] * F3[3] * (P2[3]
      - P2[0]) + F1[3] * F3[2] * (P2[0] - P2[3])));
  F2[3] = denom * 2. * cI * (F4[4] * (F1[2] * F3[3] * (P2[0] + P2[3]) - F1[3] *
      F3[2] * (P2[0] + P2[3])) + F4[5] * (F1[2] * F3[3] * (P2[1] + cI *
      (P2[2])) - F1[3] * F3[2] * (P2[1] + cI * (P2[2]))));
  F2[4] = denom * 2. * cI * F4[5] * M2 * (F1[3] * F3[2] - F1[2] * F3[3]); 
  F2[5] = denom * 2. * cI * F4[4] * M2 * (F1[2] * F3[3] - F1[3] * F3[2]); 
}

void FFFF1_6_2(complex<double> F1[], complex<double> F3[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M2, double W2,
    complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  double P2[4]; 
  complex<double> denom; 
  int i; 
  FFFF1_2(F1, F3, F4, COUP1, M2, W2, F2); 
  FFFF6_2(F1, F3, F4, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFFF1_2_6_2(complex<double> F1[], complex<double> F3[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  double P2[4]; 
  complex<double> denom; 
  int i; 
  FFFF1_2(F1, F3, F4, COUP1, M2, W2, F2); 
  FFFF2_2(F1, F3, F4, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
  FFFF6_2(F1, F3, F4, COUP3, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFFF6_2(complex<double> F1[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + F3[0] + F4[0]; 
  F2[1] = +F1[1] + F3[1] + F4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * 2. * cI * F4[3] * M2 * (F1[5] * F3[4] - F1[4] * F3[5]); 
  F2[3] = denom * 2. * cI * F4[2] * M2 * (F1[4] * F3[5] - F1[5] * F3[4]); 
  F2[4] = denom * 2. * cI * (F4[2] * (F1[4] * F3[5] * (P2[1] - cI * (P2[2])) +
      F1[5] * F3[4] * (+cI * (P2[2]) - P2[1])) + F4[3] * (F1[4] * - F3[5] *
      (P2[0] + P2[3]) + F1[5] * F3[4] * (P2[0] + P2[3])));
  F2[5] = denom * 2. * cI * (F4[2] * (F1[4] * F3[5] * (P2[0] - P2[3]) + F1[5] *
      F3[4] * (P2[3] - P2[0])) + F4[3] * (F1[4] * - F3[5] * (P2[1] + cI *
      (P2[2])) + F1[5] * F3[4] * (P2[1] + cI * (P2[2]))));
}


void VVVVV9P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP31; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP40; 
  complex<double> TMP24; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP40 = (P2[0] * V5[2] - P2[1] * V5[3] - P2[2] * V5[4] - P2[3] * V5[5]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P2[0] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP18 * (-cI * (V3[2] * TMP24) + cI * (V4[2] * TMP31)) + (TMP40 * (-cI *
      (V4[2] * TMP10) + cI * (V3[2] * TMP11)) + V2[2] * (-cI * (TMP21 * TMP31)
      + cI * (TMP19 * TMP24)))));
  V1[3] = denom * (P2[1] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP18 * (-cI * (V3[3] * TMP24) + cI * (V4[3] * TMP31)) + (TMP40 * (-cI *
      (V4[3] * TMP10) + cI * (V3[3] * TMP11)) + V2[3] * (-cI * (TMP21 * TMP31)
      + cI * (TMP19 * TMP24)))));
  V1[4] = denom * (P2[2] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP18 * (-cI * (V3[4] * TMP24) + cI * (V4[4] * TMP31)) + (TMP40 * (-cI *
      (V4[4] * TMP10) + cI * (V3[4] * TMP11)) + V2[4] * (-cI * (TMP21 * TMP31)
      + cI * (TMP19 * TMP24)))));
  V1[5] = denom * (P2[3] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP18 * (-cI * (V3[5] * TMP24) + cI * (V4[5] * TMP31)) + (TMP40 * (-cI *
      (V4[5] * TMP10) + cI * (V3[5] * TMP11)) + V2[5] * (-cI * (TMP21 * TMP31)
      + cI * (TMP19 * TMP24)))));
}


void VVVVV6P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP32; 
  double P4[4]; 
  complex<double> TMP42; 
  complex<double> TMP25; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP42 = (V5[2] * P4[0] - V5[3] * P4[1] - V5[4] * P4[2] - V5[5] * P4[3]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P4[0] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP21 * (-cI * (V2[2] * TMP25) + cI * (V3[2] * TMP32)) + (TMP42 * (-cI *
      (V3[2] * TMP11) + cI * (V2[2] * TMP22)) + V4[2] * (-cI * (TMP19 * TMP32)
      + cI * (TMP18 * TMP25)))));
  V1[3] = denom * (P4[1] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP21 * (-cI * (V2[3] * TMP25) + cI * (V3[3] * TMP32)) + (TMP42 * (-cI *
      (V3[3] * TMP11) + cI * (V2[3] * TMP22)) + V4[3] * (-cI * (TMP19 * TMP32)
      + cI * (TMP18 * TMP25)))));
  V1[4] = denom * (P4[2] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP21 * (-cI * (V2[4] * TMP25) + cI * (V3[4] * TMP32)) + (TMP42 * (-cI *
      (V3[4] * TMP11) + cI * (V2[4] * TMP22)) + V4[4] * (-cI * (TMP19 * TMP32)
      + cI * (TMP18 * TMP25)))));
  V1[5] = denom * (P4[3] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP21 * (-cI * (V2[5] * TMP25) + cI * (V3[5] * TMP32)) + (TMP42 * (-cI *
      (V3[5] * TMP11) + cI * (V2[5] * TMP22)) + V4[5] * (-cI * (TMP19 * TMP32)
      + cI * (TMP18 * TMP25)))));
}


void FFFF4_0(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> F4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP48; 
  complex<double> TMP43; 
  TMP43 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  TMP48 = (F4[2] * F3[2] + F4[3] * F3[3]); 
  vertex = COUP * - cI * TMP48 * TMP43; 
}

void FFFF4_5_0(complex<double> F1[], complex<double> F2[], complex<double>
    F3[], complex<double> F4[], complex<double> COUP1, complex<double> COUP2,
    complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFFF4_0(F1, F2, F3, F4, COUP1, vertex); 
  FFFF5_0(F1, F2, F3, F4, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void VVVV5_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  complex<double> TMP22; 
  complex<double> TMP23; 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  vertex = COUP * (-cI * (TMP11 * TMP12) + cI * (TMP22 * TMP23)); 
}


void FFFF3_3(complex<double> F1[], complex<double> F2[], complex<double> F4[],
    complex<double> COUP, double M3, double W3, complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  F3[0] = +F1[0] + F2[0] + F4[0]; 
  F3[1] = +F1[1] + F2[1] + F4[1]; 
  P3[0] = -F3[0].real(); 
  P3[1] = -F3[1].real(); 
  P3[2] = -F3[1].imag(); 
  P3[3] = -F3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  F3[2] = denom * 2. * cI * F1[3] * M3 * (F2[5] * F4[4] - F2[4] * F4[5]); 
  F3[3] = denom * 2. * cI * F1[2] * M3 * (F2[4] * F4[5] - F2[5] * F4[4]); 
  F3[4] = denom * - 2. * cI * (F1[2] * (F2[4] * - F4[5] * (P3[1] + cI *
      (P3[2])) + F2[5] * F4[4] * (P3[1] + cI * (P3[2]))) + F1[3] * (F2[4] *
      F4[5] * (P3[3] - P3[0]) + F2[5] * F4[4] * (P3[0] - P3[3])));
  F3[5] = denom * - 2. * cI * (F1[2] * (F2[4] * F4[5] * (P3[0] + P3[3]) - F2[5]
      * F4[4] * (P3[0] + P3[3])) + F1[3] * (F2[4] * F4[5] * (P3[1] - cI *
      (P3[2])) + F2[5] * F4[4] * (+cI * (P3[2]) - P3[1])));
}

void FFFF3_7_3(complex<double> F1[], complex<double> F2[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M3, double W3,
    complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFFF3_3(F1, F2, F4, COUP1, M3, W3, F3); 
  FFFF7_3(F1, F2, F4, COUP2, M3, W3, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F3[i] = F3[i] + Ftmp[i]; 
    i++; 
  }
}

void VVVV4P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V4[2] * TMP10) + cI * (V2[2] * TMP22)); 
  V1[3] = denom * (-cI * (V4[3] * TMP10) + cI * (V2[3] * TMP22)); 
  V1[4] = denom * (-cI * (V4[4] * TMP10) + cI * (V2[4] * TMP22)); 
  V1[5] = denom * (-cI * (V4[5] * TMP10) + cI * (V2[5] * TMP22)); 
}


void FFFF3_0(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> F4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP52; 
  TMP52 = (F1[2] * F3[3] * (F2[4] * F4[5] - F2[5] * F4[4]) + F1[3] * F3[2] *
      (F2[5] * F4[4] - F2[4] * F4[5]));
  vertex = COUP * - 2. * cI * TMP52; 
}

void FFFF3_7_0(complex<double> F1[], complex<double> F2[], complex<double>
    F3[], complex<double> F4[], complex<double> COUP1, complex<double> COUP2,
    complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFFF3_0(F1, F2, F3, F4, COUP1, vertex); 
  FFFF7_0(F1, F2, F3, F4, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFFF4_4(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> COUP, double M4, double W4, complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P4[4]; 
  complex<double> TMP43; 
  F4[0] = +F1[0] + F2[0] + F3[0]; 
  F4[1] = +F1[1] + F2[1] + F3[1]; 
  P4[0] = -F4[0].real(); 
  P4[1] = -F4[1].real(); 
  P4[2] = -F4[1].imag(); 
  P4[3] = -F4[0].imag(); 
  TMP43 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  F4[2] = denom * cI * F3[2] * TMP43 * M4; 
  F4[3] = denom * cI * F3[3] * TMP43 * M4; 
  F4[4] = denom * - cI * TMP43 * (F3[2] * - 1. * (P4[0] + P4[3]) + F3[3] * (+cI
      * (P4[2]) - P4[1]));
  F4[5] = denom * cI * TMP43 * (F3[2] * (P4[1] + cI * (P4[2])) + F3[3] * (P4[0]
      - P4[3]));
}

void FFFF4_5_4(complex<double> F1[], complex<double> F2[], complex<double>
    F3[], complex<double> COUP1, complex<double> COUP2, double M4, double W4,
    complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> Ftmp[6]; 
  double P4[4]; 
  int i; 
  FFFF4_4(F1, F2, F3, COUP1, M4, W4, F4); 
  FFFF5_4(F1, F2, F3, COUP2, M4, W4, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F4[i] = F4[i] + Ftmp[i]; 
    i++; 
  }
}

void FFVV2_1(complex<double> F2[], complex<double> V3[], complex<double> V4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0] + V4[0]; 
  F1[1] = +F2[1] + V3[1] + V4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - 2. * cI * (F2[4] * (P1[1] * (V3[2] * (+cI * (V4[4]) -
      V4[3]) + (V3[3] * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] +
      V4[5])) + V3[5] * (+cI * (V4[4]) - V4[3])))) + (P1[2] * (V3[2] * - 1. *
      (V4[4] + cI * (V4[3])) + (V3[3] * (+cI * (V4[2] + V4[5])) + (V3[4] *
      (V4[2] + V4[5]) - V3[5] * (V4[4] + cI * (V4[3]))))) + (P1[0] * (V3[5] *
      V4[2] - cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[2] * V4[5]) +
      P1[3] * (V3[5] * V4[2] - cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) -
      V3[2] * V4[5])))) + F2[5] * (P1[0] * (V3[2] * - 1. * (V4[3] + cI *
      (V4[4])) + (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[5]) + cI *
      (V4[2])) + V3[5] * (V4[3] + cI * (V4[4]))))) + (P1[3] * (V3[2] * - 1. *
      (V4[3] + cI * (V4[4])) + (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI *
      (V4[5]) + cI * (V4[2])) + V3[5] * (V4[3] + cI * (V4[4]))))) + (P1[1] *
      (V3[2] * V4[5] - cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[5] *
      V4[2]) + P1[2] * (V3[4] * V4[3] - cI * (V3[5] * V4[2]) + cI * (V3[2] *
      V4[5]) - V3[3] * V4[4])))));
  F1[3] = denom * - 2. * cI * (F2[4] * (P1[0] * (V3[2] * (+cI * (V4[4]) -
      V4[3]) + (V3[3] * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] +
      V4[5])) + V3[5] * (+cI * (V4[4]) - V4[3])))) + (P1[3] * (V3[2] * (V4[3] -
      cI * (V4[4])) + (V3[3] * - 1. * (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2]
      + V4[5])) + V3[5] * (V4[3] - cI * (V4[4]))))) + (P1[1] * (V3[5] * V4[2] -
      cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[2] * V4[5]) + P1[2] *
      (V3[4] * V4[3] - cI * (V3[5] * V4[2]) + cI * (V3[2] * V4[5]) - V3[3] *
      V4[4])))) + F2[5] * (P1[1] * (V3[2] * - 1. * (V4[3] + cI * (V4[4])) +
      (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) +
      V3[5] * (V4[3] + cI * (V4[4]))))) + (P1[2] * (V3[2] * (+cI * (V4[3]) -
      V4[4]) + (V3[3] * (-cI * (V4[2]) + cI * (V4[5])) + (V3[4] * (V4[2] -
      V4[5]) + V3[5] * (V4[4] - cI * (V4[3]))))) + (P1[0] * (V3[2] * V4[5] - cI
      * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[5] * V4[2]) + P1[3] *
      (V3[5] * V4[2] - cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[2] *
      V4[5])))));
  F1[4] = denom * 2. * cI * M1 * (F2[5] * (V3[2] * - 1. * (V4[3] + cI *
      (V4[4])) + (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[5]) + cI *
      (V4[2])) + V3[5] * (V4[3] + cI * (V4[4]))))) + F2[4] * (V3[5] * V4[2] -
      cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[2] * V4[5]));
  F1[5] = denom * 2. * cI * M1 * (F2[4] * (V3[2] * (+cI * (V4[4]) - V4[3]) +
      (V3[3] * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5])) +
      V3[5] * (+cI * (V4[4]) - V4[3])))) + F2[5] * (V3[2] * V4[5] - cI * (V3[4]
      * V4[3]) + cI * (V3[3] * V4[4]) - V3[5] * V4[2]));
}


void FFFF7_3(complex<double> F1[], complex<double> F2[], complex<double> F4[],
    complex<double> COUP, double M3, double W3, complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  F3[0] = +F1[0] + F2[0] + F4[0]; 
  F3[1] = +F1[1] + F2[1] + F4[1]; 
  P3[0] = -F3[0].real(); 
  P3[1] = -F3[1].real(); 
  P3[2] = -F3[1].imag(); 
  P3[3] = -F3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  F3[2] = denom * - 2. * cI * (F1[4] * (F2[2] * F4[3] * (P3[1] + cI * (P3[2]))
      - F2[3] * F4[2] * (P3[1] + cI * (P3[2]))) + F1[5] * (F2[2] * - F4[3] *
      (P3[0] + P3[3]) + F2[3] * F4[2] * (P3[0] + P3[3])));
  F3[3] = denom * - 2. * cI * (F1[4] * (F2[2] * F4[3] * (P3[0] - P3[3]) + F2[3]
      * F4[2] * (P3[3] - P3[0])) + F1[5] * (F2[2] * F4[3] * (+cI * (P3[2]) -
      P3[1]) + F2[3] * F4[2] * (P3[1] - cI * (P3[2]))));
  F3[4] = denom * 2. * cI * F1[5] * M3 * (F2[3] * F4[2] - F2[2] * F4[3]); 
  F3[5] = denom * 2. * cI * F1[4] * M3 * (F2[2] * F4[3] - F2[3] * F4[2]); 
}


void VVVVV10P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP45; 
  double P3[4]; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP45 = (P3[0] * V5[2] - P3[1] * V5[3] - P3[2] * V5[4] - P3[3] * V5[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P3[0] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP19 * (-cI * (V2[2] * TMP26) + cI * (V4[2] * TMP29)) + (TMP45 * (-cI *
      (V4[2] * TMP10) + cI * (V2[2] * TMP22)) + V3[2] * (-cI * (TMP21 * TMP29)
      + cI * (TMP18 * TMP26)))));
  V1[3] = denom * (P3[1] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP19 * (-cI * (V2[3] * TMP26) + cI * (V4[3] * TMP29)) + (TMP45 * (-cI *
      (V4[3] * TMP10) + cI * (V2[3] * TMP22)) + V3[3] * (-cI * (TMP21 * TMP29)
      + cI * (TMP18 * TMP26)))));
  V1[4] = denom * (P3[2] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP19 * (-cI * (V2[4] * TMP26) + cI * (V4[4] * TMP29)) + (TMP45 * (-cI *
      (V4[4] * TMP10) + cI * (V2[4] * TMP22)) + V3[4] * (-cI * (TMP21 * TMP29)
      + cI * (TMP18 * TMP26)))));
  V1[5] = denom * (P3[3] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP19 * (-cI * (V2[5] * TMP26) + cI * (V4[5] * TMP29)) + (TMP45 * (-cI *
      (V4[5] * TMP10) + cI * (V2[5] * TMP22)) + V3[5] * (-cI * (TMP21 * TMP29)
      + cI * (TMP18 * TMP26)))));
}


void FFV7_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * 4. * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3]
      - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] -
      V3[2])))) + (+1./4. * (M1 * (F2[5] * (V3[3] + cI * (V3[4])) + 4. * (F2[4]
      * 1./4. * (V3[2] + V3[5])))) + F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) +
      (P1[1] * - 1. * (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5]))
      + P1[3] * (V3[3] + cI * (V3[4])))))));
  F1[3] = denom * 4. * cI * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (+1./4. * (M1 * (F2[5] * (V3[2] - V3[5]) + 4. *
      (F2[4] * 1./4. * (V3[3] - cI * (V3[4]))))) + F2[3] * (P1[0] * - 1. *
      (V3[2] + V3[5]) + (P1[1] * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI
      * (V3[3])) + P1[3] * (V3[2] + V3[5]))))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 4. * (V3[5] - V3[2]) + 4. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * 4. * (+cI * (V3[4]) - V3[3]) + 4. * (F2[3] *
      (V3[2] + V3[5])))));
}


void FFV2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP4; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP4 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP4); 
  V3[3] = denom * - cI * (-F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 * TMP4); 
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5]) + cI * (F1[3] * F2[4]) - P3[2]
      * OM3 * TMP4);
  V3[5] = denom * - cI * (F1[3] * F2[5] - F1[2] * F2[4] - P3[3] * OM3 * TMP4); 
}

void FFV2_4_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  double OM3; 
  complex<double> Vtmp[6]; 
  int i; 
  FFV2_3(F1, F2, COUP1, M3, W3, V3); 
  FFV4_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV2_7_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  double OM3; 
  complex<double> Vtmp[6]; 
  int i; 
  FFV2_3(F1, F2, COUP1, M3, W3, V3); 
  FFV7_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void FFS2_1(complex<double> F2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * S3[2] * (F2[4] * (P1[0] + P1[3]) + (F2[5] * (P1[1] +
      cI * (P1[2])) - F2[2] * M1));
  F1[3] = denom * cI * S3[2] * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) + F2[3] * M1));
  F1[4] = denom * cI * S3[2] * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] + cI
      * (P1[2])) + F2[4] * M1));
  F1[5] = denom * - cI * S3[2] * (F2[2] * (+cI * (P1[2]) - P1[1]) + (F2[3] *
      (P1[0] + P1[3]) - F2[5] * M1));
}


void FFV5_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP14; 
  complex<double> TMP13; 
  TMP14 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP13 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - 1. * (+cI * (TMP13 + TMP14)); 
}


void FFVV2_2(complex<double> F1[], complex<double> V3[], complex<double> V4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0] + V4[0]; 
  F2[1] = +F1[1] + V3[1] + V4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * 2. * cI * (F1[4] * (P2[1] * (V3[2] * (V4[3] + cI * (V4[4])) +
      (V3[3] * (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[2]) + cI * (V4[5])) -
      V3[5] * (V4[3] + cI * (V4[4]))))) + (P2[2] * (V3[2] * (V4[4] - cI *
      (V4[3])) + (V3[3] * (-cI * (V4[5]) + cI * (V4[2])) + (V3[4] * (V4[5] -
      V4[2]) + V3[5] * (+cI * (V4[3]) - V4[4])))) + (P2[0] * (V3[5] * V4[2] -
      cI * (V3[3] * V4[4]) + cI * (V3[4] * V4[3]) - V3[2] * V4[5]) + P2[3] *
      (V3[2] * V4[5] - cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[5] *
      V4[2])))) + F1[5] * (P2[0] * (V3[2] * (+cI * (V4[4]) - V4[3]) + (V3[3] *
      (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5])) + V3[5] * (+cI
      * (V4[4]) - V4[3])))) + (P2[3] * (V3[2] * (V4[3] - cI * (V4[4])) + (V3[3]
      * - 1. * (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] + V4[5])) + V3[5] *
      (V4[3] - cI * (V4[4]))))) + (P2[1] * (V3[5] * V4[2] - cI * (V3[3] *
      V4[4]) + cI * (V3[4] * V4[3]) - V3[2] * V4[5]) + P2[2] * (V3[4] * V4[3] -
      cI * (V3[5] * V4[2]) + cI * (V3[2] * V4[5]) - V3[3] * V4[4])))));
  F2[3] = denom * 2. * cI * (F1[4] * (P2[0] * (V3[2] * - 1. * (V4[3] + cI *
      (V4[4])) + (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[5]) + cI *
      (V4[2])) + V3[5] * (V4[3] + cI * (V4[4]))))) + (P2[3] * (V3[2] * - 1. *
      (V4[3] + cI * (V4[4])) + (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI *
      (V4[5]) + cI * (V4[2])) + V3[5] * (V4[3] + cI * (V4[4]))))) + (P2[1] *
      (V3[2] * V4[5] - cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[5] *
      V4[2]) + P2[2] * (V3[4] * V4[3] - cI * (V3[5] * V4[2]) + cI * (V3[2] *
      V4[5]) - V3[3] * V4[4])))) + F1[5] * (P2[1] * (V3[2] * (V4[3] - cI *
      (V4[4])) + (V3[3] * - 1. * (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] +
      V4[5])) + V3[5] * (V4[3] - cI * (V4[4]))))) + (P2[2] * (V3[2] * (V4[4] +
      cI * (V4[3])) + (V3[3] * - 1. * (+cI * (V4[2] + V4[5])) + (V3[4] * - 1. *
      (V4[2] + V4[5]) + V3[5] * (V4[4] + cI * (V4[3]))))) + (P2[0] * (V3[2] *
      V4[5] - cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[5] * V4[2]) +
      P2[3] * (V3[2] * V4[5] - cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) -
      V3[5] * V4[2])))));
  F2[4] = denom * 2. * cI * M2 * (F1[5] * (V3[2] * (+cI * (V4[4]) - V4[3]) +
      (V3[3] * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5])) +
      V3[5] * (+cI * (V4[4]) - V4[3])))) + F1[4] * (V3[5] * V4[2] - cI * (V3[3]
      * V4[4]) + cI * (V3[4] * V4[3]) - V3[2] * V4[5]));
  F2[5] = denom * 2. * cI * M2 * (F1[4] * (V3[2] * - 1. * (V4[3] + cI *
      (V4[4])) + (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[5]) + cI *
      (V4[2])) + V3[5] * (V4[3] + cI * (V4[4]))))) + F1[5] * (V3[2] * V4[5] -
      cI * (V3[4] * V4[3]) + cI * (V3[3] * V4[4]) - V3[5] * V4[2]));
}


void FFV5_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * (+cI *
      (V3[4]) - V3[3]))));
  F2[3] = denom * - cI * (F1[2] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] * (V3[2] +
      V3[5]))));
  F2[4] = denom * - cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] +
      cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * - 1. *
      (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] *
      (V3[2] - V3[5]))));
}


void FFV3_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  double P3[4]; 
  complex<double> denom; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * - 2. * cI * M2 * (F1[3] * (P3[0] * (V3[3] - cI * (V3[4])) +
      (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + F1[2] * (V3[5] * P3[0] - cI * (V3[4]
      * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]));
  F2[3] = denom * - 2. * cI * M2 * (F1[2] * (P3[0] * (V3[3] + cI * (V3[4])) +
      (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5]))
      + P3[3] * (V3[3] + cI * (V3[4]))))) + F1[3] * (V3[2] * P3[3] - cI *
      (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]));
  F2[4] = denom * - 2. * cI * (F1[2] * (P2[1] * (P3[0] * (V3[3] + cI * (V3[4]))
      + (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] +
      V3[5])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[2] * (P3[0] * (V3[4] -
      cI * (V3[3])) + (P3[1] * (+cI * (V3[2] + V3[5])) + (P3[2] * - 1. * (V3[2]
      + V3[5]) + P3[3] * (V3[4] - cI * (V3[3]))))) + (P2[0] * (V3[5] * P3[0] -
      cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]) + P2[3] *
      (V3[5] * P3[0] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] *
      P3[3])))) + F1[3] * (P2[0] * (P3[0] * (V3[3] - cI * (V3[4])) + (P3[1] *
      (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) + P3[3] * (+cI
      * (V3[4]) - V3[3])))) + (P2[3] * (P3[0] * (V3[3] - cI * (V3[4])) + (P3[1]
      * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) + P3[3] *
      (+cI * (V3[4]) - V3[3])))) + (P2[1] * (V3[2] * P3[3] - cI * (V3[3] *
      P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]) + P2[2] * (V3[4] * P3[1] -
      cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[3] * P3[2])))));
  F2[5] = denom * - 2. * cI * (F1[2] * (P2[0] * (P3[0] * (V3[3] + cI * (V3[4]))
      + (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] +
      V3[5])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[3] * (P3[0] * - 1. *
      (V3[3] + cI * (V3[4])) + (P3[1] * (V3[2] + V3[5]) + (P3[2] * (+cI *
      (V3[2] + V3[5])) - P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[1] * (V3[5] *
      P3[0] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]) +
      P2[2] * (V3[4] * P3[1] - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) -
      V3[3] * P3[2])))) + F1[3] * (P2[1] * (P3[0] * (V3[3] - cI * (V3[4])) +
      (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + (P2[2] * (P3[0] * (V3[4] + cI *
      (V3[3])) + (P3[1] * (-cI * (V3[2]) + cI * (V3[5])) + (P3[2] * (V3[5] -
      V3[2]) - P3[3] * (V3[4] + cI * (V3[3]))))) + (P2[0] * (V3[2] * P3[3] - cI
      * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]) + P2[3] *
      (V3[5] * P3[0] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] *
      P3[3])))));
}

void FFV3_8_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> denom; 
  int i; 
  FFV3_2(F1, V3, COUP1, M2, W2, F2); 
  FFV8_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFFF2_0(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> F4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP51; 
  TMP51 = (F1[2] * F3[3] * (F2[5] * F4[4] - F2[4] * F4[5]) + F1[3] * F3[2] *
      (F2[4] * F4[5] - F2[5] * F4[4]));
  vertex = COUP * - 2. * cI * TMP51; 
}


void VVS2P0_1(complex<double> V2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP30; 
  complex<double> TMP47; 
  complex<double> denom; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP47 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * S3[2] * (-cI * (P2[0] * TMP30) + cI * (V2[2] * TMP47)); 
  V1[3] = denom * S3[2] * (-cI * (P2[1] * TMP30) + cI * (V2[3] * TMP47)); 
  V1[4] = denom * S3[2] * (-cI * (P2[2] * TMP30) + cI * (V2[4] * TMP47)); 
  V1[5] = denom * S3[2] * (-cI * (P2[3] * TMP30) + cI * (V2[5] * TMP47)); 
}


void FFV8_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  complex<double> TMP0; 
  double P3[4]; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP1 = (F1[4] * (F2[4] * (P3[0] * (V3[2] + V3[5]) + (P3[1] * (+cI * (V3[4]) -
      V3[3]) + (P3[2] * - 1. * (V3[4] + cI * (V3[3])) - P3[3] * (V3[2] +
      V3[5])))) + F2[5] * (P3[0] * (V3[3] + cI * (V3[4])) + (P3[1] * (V3[5] -
      V3[2]) + (P3[2] * (-cI * (V3[2]) + cI * (V3[5])) - P3[3] * (V3[3] + cI *
      (V3[4])))))) + F1[5] * (F2[4] * (P3[0] * (V3[3] - cI * (V3[4])) + (P3[1]
      * - 1. * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) + P3[3] *
      (V3[3] - cI * (V3[4]))))) + F2[5] * (P3[0] * (V3[2] - V3[5]) + (P3[1] * -
      1. * (V3[3] + cI * (V3[4])) + (P3[2] * (+cI * (V3[3]) - V3[4]) + P3[3] *
      (V3[2] - V3[5]))))));
  TMP0 = (F1[4] * (F2[4] * (P3[0] * (V3[2] - V3[5]) + (P3[1] * - 1. * (V3[3] +
      cI * (V3[4])) + (P3[2] * (+cI * (V3[3]) - V3[4]) + P3[3] * (V3[2] -
      V3[5])))) + F2[5] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) + (P3[1] *
      (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) + P3[3] *
      (V3[3] + cI * (V3[4])))))) + F1[5] * (F2[4] * (P3[0] * (+cI * (V3[4]) -
      V3[3]) + (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] +
      V3[5])) + P3[3] * (+cI * (V3[4]) - V3[3])))) + F2[5] * (P3[0] * (V3[2] +
      V3[5]) + (P3[1] * (+cI * (V3[4]) - V3[3]) + (P3[2] * - 1. * (V3[4] + cI *
      (V3[3])) - P3[3] * (V3[2] + V3[5]))))));
  vertex = COUP * (-cI * (TMP1) + cI * (TMP0)); 
}


void FFFF6_1(complex<double> F2[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + F3[0] + F4[0]; 
  F1[1] = +F2[1] + F3[1] + F4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - 2. * cI * (F3[4] * (F2[2] * F4[3] * (P1[1] + cI * (P1[2]))
      - F2[3] * F4[2] * (P1[1] + cI * (P1[2]))) + F3[5] * (F2[2] * - F4[3] *
      (P1[0] + P1[3]) + F2[3] * F4[2] * (P1[0] + P1[3])));
  F1[3] = denom * - 2. * cI * (F3[4] * (F2[2] * F4[3] * (P1[0] - P1[3]) + F2[3]
      * F4[2] * (P1[3] - P1[0])) + F3[5] * (F2[2] * F4[3] * (+cI * (P1[2]) -
      P1[1]) + F2[3] * F4[2] * (P1[1] - cI * (P1[2]))));
  F1[4] = denom * 2. * cI * F3[5] * M1 * (F2[3] * F4[2] - F2[2] * F4[3]); 
  F1[5] = denom * 2. * cI * F3[4] * M1 * (F2[2] * F4[3] - F2[3] * F4[2]); 
}


void VVV2P0_1(complex<double> V2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP31; 
  double P3[4]; 
  complex<double> TMP46; 
  complex<double> TMP30; 
  complex<double> TMP47; 
  complex<double> denom; 
  complex<double> TMP27; 
  complex<double> TMP29; 
  complex<double> TMP35; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP46 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  TMP47 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP35 = (P3[0] * P1[0] - P3[1] * P1[1] - P3[2] * P1[2] - P3[3] * P1[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P2[0] * (-cI * (TMP10 * TMP35) + cI * (TMP27 * TMP29)) +
      (P3[0] * (-cI * (TMP30 * TMP31) + cI * (TMP10 * TMP47)) + (TMP46 * (-cI *
      (V2[2] * TMP27) + cI * (V3[2] * TMP30)) + (-cI * (V3[2] * TMP29 * TMP47)
      + cI * (V2[2] * TMP31 * TMP35)))));
  V1[3] = denom * (P2[1] * (-cI * (TMP10 * TMP35) + cI * (TMP27 * TMP29)) +
      (P3[1] * (-cI * (TMP30 * TMP31) + cI * (TMP10 * TMP47)) + (TMP46 * (-cI *
      (V2[3] * TMP27) + cI * (V3[3] * TMP30)) + (-cI * (V3[3] * TMP29 * TMP47)
      + cI * (V2[3] * TMP31 * TMP35)))));
  V1[4] = denom * (P2[2] * (-cI * (TMP10 * TMP35) + cI * (TMP27 * TMP29)) +
      (P3[2] * (-cI * (TMP30 * TMP31) + cI * (TMP10 * TMP47)) + (TMP46 * (-cI *
      (V2[4] * TMP27) + cI * (V3[4] * TMP30)) + (-cI * (V3[4] * TMP29 * TMP47)
      + cI * (V2[4] * TMP31 * TMP35)))));
  V1[5] = denom * (P2[3] * (-cI * (TMP10 * TMP35) + cI * (TMP27 * TMP29)) +
      (P3[3] * (-cI * (TMP30 * TMP31) + cI * (TMP10 * TMP47)) + (TMP46 * (-cI *
      (V2[5] * TMP27) + cI * (V3[5] * TMP30)) + (-cI * (V3[5] * TMP29 * TMP47)
      + cI * (V2[5] * TMP31 * TMP35)))));
}


void FFVV2P0_3(complex<double> F1[], complex<double> F2[], complex<double>
    V4[], complex<double> COUP, double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0] + V4[0]; 
  V3[1] = +F1[1] + F2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * 2. * cI * (F1[4] * (F2[5] * (V4[3] + cI * (V4[4])) + F2[4] *
      V4[5]) + F1[5] * (F2[4] * (V4[3] - cI * (V4[4])) - F2[5] * V4[5]));
  V3[3] = denom * 2. * cI * (F1[4] * (F2[5] * (V4[2] - V4[5]) - cI * (F2[4] *
      V4[4])) + F1[5] * (F2[4] * (V4[2] + V4[5]) + cI * (F2[5] * V4[4])));
  V3[4] = denom * 2. * (F1[4] * (F2[5] * (V4[5] - V4[2]) - F2[4] * V4[3]) +
      F1[5] * (F2[4] * (V4[2] + V4[5]) + F2[5] * V4[3]));
  V3[5] = denom * 2. * cI * (F1[4] * (F2[5] * (V4[3] + cI * (V4[4])) + F2[4] *
      V4[2]) + F1[5] * (F2[4] * (+cI * (V4[4]) - V4[3]) - F2[5] * V4[2]));
}


void FFV4_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP15; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP4; 
  complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP15 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP4 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * 2. * cI * (OM3 * 1./2. * P3[0] * (TMP4 - 2. * (TMP15)) +
      (-1./2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] *
      F2[3]));
  V3[3] = denom * 2. * cI * (OM3 * 1./2. * P3[1] * (TMP4 - 2. * (TMP15)) +
      (+1./2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] *
      F2[2]));
  V3[4] = denom * - 2. * cI * (OM3 * 1./2. * P3[2] * (+2. * (TMP15) - TMP4) +
      (-1./2. * cI * (F1[2] * F2[5]) + 1./2. * cI * (F1[3] * F2[4]) - cI *
      (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
  V3[5] = denom * - 2. * cI * (OM3 * 1./2. * P3[3] * (+2. * (TMP15) - TMP4) +
      (-1./2. * (F1[2] * F2[4]) + 1./2. * (F1[3] * F2[5]) - F1[4] * F2[2] +
      F1[5] * F2[3]));
}


void VVVV2_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP10; 
  double P3[4]; 
  complex<double> TMP47; 
  complex<double> TMP32; 
  complex<double> TMP28; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  complex<double> TMP37; 
  double P2[4]; 
  complex<double> TMP33; 
  complex<double> TMP53; 
  complex<double> TMP34; 
  complex<double> TMP24; 
  complex<double> TMP12; 
  double P1[4]; 
  complex<double> TMP23; 
  complex<double> TMP30; 
  complex<double> TMP27; 
  double P4[4]; 
  complex<double> TMP11; 
  complex<double> TMP22; 
  complex<double> TMP31; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP47 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP53 = (P3[0] * P4[0] - P3[1] * P4[1] - P3[2] * P4[2] - P3[3] * P4[3]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP33 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP34 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  vertex = COUP * (TMP10 * (TMP9 * - 1. * (+cI * (TMP47 + TMP53)) + (+cI *
      (TMP28 * TMP37 + TMP26 * TMP34))) + (TMP11 * (TMP12 * (+cI * (TMP47 +
      TMP53)) + (-cI * (TMP27 * TMP37 + TMP25 * TMP33))) + (TMP12 * - 1. * (+cI
      * (TMP24 * TMP30 + TMP26 * TMP32)) + (TMP22 * (-cI * (TMP29 * TMP34) + cI
      * (TMP32 * TMP33)) + (TMP23 * (-cI * (TMP28 * TMP31) + cI * (TMP24 *
      TMP27)) + TMP9 * (+cI * (TMP30 * TMP31 + TMP25 * TMP29)))))));
}


void VVVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V4[2] * TMP10) + cI * (V3[2] * TMP11)); 
  V1[3] = denom * (-cI * (V4[3] * TMP10) + cI * (V3[3] * TMP11)); 
  V1[4] = denom * (-cI * (V4[4] * TMP10) + cI * (V3[4] * TMP11)); 
  V1[5] = denom * (-cI * (V4[5] * TMP10) + cI * (V3[5] * TMP11)); 
}


void FFFF5_4(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> COUP, double M4, double W4, complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP39; 
  double P4[4]; 
  F4[0] = +F1[0] + F2[0] + F3[0]; 
  F4[1] = +F1[1] + F2[1] + F3[1]; 
  P4[0] = -F4[0].real(); 
  P4[1] = -F4[1].real(); 
  P4[2] = -F4[1].imag(); 
  P4[3] = -F4[0].imag(); 
  TMP39 = (F1[2] * F2[2] + F1[3] * F2[3]); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  F4[2] = denom * cI * TMP39 * (F3[4] * (P4[0] - P4[3]) + F3[5] * (+cI *
      (P4[2]) - P4[1]));
  F4[3] = denom * - cI * TMP39 * (F3[4] * (P4[1] + cI * (P4[2])) - F3[5] *
      (P4[0] + P4[3]));
  F4[4] = denom * cI * F3[4] * TMP39 * M4; 
  F4[5] = denom * cI * F3[5] * TMP39 * M4; 
}


void VVVVV12P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP32; 
  double P4[4]; 
  complex<double> TMP42; 
  complex<double> TMP25; 
  complex<double> TMP19; 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP42 = (V5[2] * P4[0] - V5[3] * P4[1] - V5[4] * P4[2] - V5[5] * P4[3]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P4[0] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP22 * (-cI * (V5[2] * TMP32) + cI * (V2[2] * TMP42)) + (TMP25 * (-cI *
      (V2[2] * TMP21) + cI * (TMP11 * V5[2])) + V4[2] * (-cI * (TMP10 * TMP42)
      + cI * (TMP19 * TMP32)))));
  V1[3] = denom * (P4[1] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP22 * (-cI * (V5[3] * TMP32) + cI * (V2[3] * TMP42)) + (TMP25 * (-cI *
      (V2[3] * TMP21) + cI * (TMP11 * V5[3])) + V4[3] * (-cI * (TMP10 * TMP42)
      + cI * (TMP19 * TMP32)))));
  V1[4] = denom * (P4[2] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP22 * (-cI * (V5[4] * TMP32) + cI * (V2[4] * TMP42)) + (TMP25 * (-cI *
      (V2[4] * TMP21) + cI * (TMP11 * V5[4])) + V4[4] * (-cI * (TMP10 * TMP42)
      + cI * (TMP19 * TMP32)))));
  V1[5] = denom * (P4[3] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP22 * (-cI * (V5[5] * TMP32) + cI * (V2[5] * TMP42)) + (TMP25 * (-cI *
      (V2[5] * TMP21) + cI * (TMP11 * V5[5])) + V4[5] * (-cI * (TMP10 * TMP42)
      + cI * (TMP19 * TMP32)))));
}


void FFVV1P0_3(complex<double> F1[], complex<double> F2[], complex<double>
    V4[], complex<double> COUP, double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0] + V4[0]; 
  V3[1] = +F1[1] + F2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - 2. * cI * (F1[2] * (F2[3] * (V4[3] + cI * (V4[4])) + F2[2]
      * V4[5]) + F1[3] * (F2[2] * (V4[3] - cI * (V4[4])) - F2[3] * V4[5]));
  V3[3] = denom * - 2. * cI * (F1[2] * (F2[3] * (V4[2] + V4[5]) + cI * (F2[2] *
      V4[4])) + F1[3] * (F2[2] * (V4[2] - V4[5]) - cI * (F2[3] * V4[4])));
  V3[4] = denom * 2. * (F1[2] * (F2[3] * (V4[2] + V4[5]) - F2[2] * V4[3]) +
      F1[3] * (F2[2] * (V4[5] - V4[2]) + F2[3] * V4[3]));
  V3[5] = denom * 2. * cI * (F1[2] * (F2[3] * (V4[3] + cI * (V4[4])) - F2[2] *
      V4[2]) + F1[3] * (F2[2] * (+cI * (V4[4]) - V4[3]) + F2[3] * V4[2]));
}

void FFVV1_2P0_3(complex<double> F1[], complex<double> F2[], complex<double>
    V4[], complex<double> COUP1, complex<double> COUP2, double M3, double W3,
    complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Vtmp[6]; 
  FFVV1P0_3(F1, F2, V4, COUP1, M3, W3, V3); 
  FFVV2P0_3(F1, F2, V4, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void FFFF6_4(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> COUP, double M4, double W4, complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P4[4]; 
  F4[0] = +F1[0] + F2[0] + F3[0]; 
  F4[1] = +F1[1] + F2[1] + F3[1]; 
  P4[0] = -F4[0].real(); 
  P4[1] = -F4[1].real(); 
  P4[2] = -F4[1].imag(); 
  P4[3] = -F4[0].imag(); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  F4[2] = denom * 2. * cI * F2[3] * M4 * (F1[4] * F3[5] - F1[5] * F3[4]); 
  F4[3] = denom * - 2. * cI * F2[2] * M4 * (F1[4] * F3[5] - F1[5] * F3[4]); 
  F4[4] = denom * 2. * cI * (F2[2] * (F1[4] * F3[5] * (+cI * (P4[2]) - P4[1]) +
      F1[5] * F3[4] * (P4[1] - cI * (P4[2]))) + F2[3] * (F1[4] * F3[5] * (P4[0]
      + P4[3]) - F1[5] * F3[4] * (P4[0] + P4[3])));
  F4[5] = denom * 2. * cI * (F2[2] * (F1[4] * F3[5] * (P4[3] - P4[0]) + F1[5] *
      F3[4] * (P4[0] - P4[3])) + F2[3] * (F1[4] * F3[5] * (P4[1] + cI *
      (P4[2])) - F1[5] * F3[4] * (P4[1] + cI * (P4[2]))));
}


void FFFF7_4(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> COUP, double M4, double W4, complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P4[4]; 
  F4[0] = +F1[0] + F2[0] + F3[0]; 
  F4[1] = +F1[1] + F2[1] + F3[1]; 
  P4[0] = -F4[0].real(); 
  P4[1] = -F4[1].real(); 
  P4[2] = -F4[1].imag(); 
  P4[3] = -F4[0].imag(); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  F4[2] = denom * 2. * cI * F2[3] * M4 * (F1[5] * F3[4] - F1[4] * F3[5]); 
  F4[3] = denom * 2. * cI * F2[2] * M4 * (F1[4] * F3[5] - F1[5] * F3[4]); 
  F4[4] = denom * 2. * cI * (F2[2] * (F1[4] * F3[5] * (P4[1] - cI * (P4[2])) +
      F1[5] * F3[4] * (+cI * (P4[2]) - P4[1])) + F2[3] * (F1[4] * - F3[5] *
      (P4[0] + P4[3]) + F1[5] * F3[4] * (P4[0] + P4[3])));
  F4[5] = denom * 2. * cI * (F2[2] * (F1[4] * F3[5] * (P4[0] - P4[3]) + F1[5] *
      F3[4] * (P4[3] - P4[0])) + F2[3] * (F1[4] * - F3[5] * (P4[1] + cI *
      (P4[2])) + F1[5] * F3[4] * (P4[1] + cI * (P4[2]))));
}


void VVVS1_4(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, double M4, double W4, complex<double> S4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP31; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> TMP33; 
  complex<double> denom; 
  complex<double> TMP27; 
  double P4[4]; 
  complex<double> TMP29; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  S4[0] = +V1[0] + V2[0] + V3[0]; 
  S4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -S4[0].real(); 
  P4[1] = -S4[1].real(); 
  P4[2] = -S4[1].imag(); 
  P4[3] = -S4[0].imag(); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP33 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  S4[2] = denom * (TMP10 * (-cI * (TMP33) + cI * (TMP37)) + (TMP12 * (-cI *
      (TMP30) + cI * (TMP29)) + TMP23 * (-cI * (TMP31) + cI * (TMP27))));
}


void FFV8P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * 2. * cI * (F1[4] * (F2[5] * (P3[1] + cI * (P3[2])) + P3[3] *
      F2[4]) + F1[5] * (F2[4] * (P3[1] - cI * (P3[2])) - P3[3] * F2[5]));
  V3[3] = denom * 2. * cI * (F1[4] * (F2[5] * (P3[0] - P3[3]) - cI * (P3[2] *
      F2[4])) + F1[5] * (F2[4] * (P3[0] + P3[3]) + cI * (P3[2] * F2[5])));
  V3[4] = denom * 2. * (F1[4] * (F2[5] * (P3[3] - P3[0]) - P3[1] * F2[4]) +
      F1[5] * (F2[4] * (P3[0] + P3[3]) + P3[1] * F2[5]));
  V3[5] = denom * 2. * cI * (F1[4] * (F2[5] * (P3[1] + cI * (P3[2])) + P3[0] *
      F2[4]) + F1[5] * (F2[4] * (+cI * (P3[2]) - P3[1]) - P3[0] * F2[5]));
}


void FFV2_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP3; 
  TMP3 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - cI * TMP3; 
}

void FFV2_4_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFV2_0(F1, F2, V3, COUP1, vertex); 
  FFV4_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV2_7_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFV2_0(F1, F2, V3, COUP1, vertex); 
  FFV7_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV5_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI
      * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * - 1. * (V3[2] +
      V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] + cI *
      (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])))));
  F1[3] = denom * - cI * (F2[2] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] *
      - 1. * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] *
      (V3[5] - V3[2]))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (V3[5] - V3[2]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] +
      V3[5]))));
}


void FFS1_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  complex<double> TMP54; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP54 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP54; 
}


void FFV1P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5]
      * F2[3]);
  V3[3] = denom * - cI * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3]
      * F2[4]);
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] *
      F2[4] + F1[4] * F2[3]));
  V3[5] = denom * - cI * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5]
      * F2[3]);
}


void VVVVS2_5(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, double M5, double W5,
    complex<double> S5[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP23; 
  complex<double> denom; 
  double P5[4]; 
  complex<double> TMP9; 
  S5[0] = +V1[0] + V2[0] + V3[0] + V4[0]; 
  S5[1] = +V1[1] + V2[1] + V3[1] + V4[1]; 
  P5[0] = -S5[0].real(); 
  P5[1] = -S5[1].real(); 
  P5[2] = -S5[1].imag(); 
  P5[3] = -S5[0].imag(); 
  TMP9 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/(pow(P5[0], 2) - pow(P5[1], 2) - pow(P5[2], 2) - pow(P5[3], 2) -
      M5 * (M5 - cI * W5));
  S5[2] = denom * (-cI * (TMP22 * TMP23) + cI * (TMP9 * TMP10)); 
}


void FFFF7_1(complex<double> F2[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + F3[0] + F4[0]; 
  F1[1] = +F2[1] + F3[1] + F4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - 2. * cI * (F3[4] * (F2[2] * - F4[3] * (P1[1] + cI *
      (P1[2])) + F2[3] * F4[2] * (P1[1] + cI * (P1[2]))) + F3[5] * (F2[2] *
      F4[3] * (P1[0] + P1[3]) - F2[3] * F4[2] * (P1[0] + P1[3])));
  F1[3] = denom * - 2. * cI * (F3[4] * (F2[2] * F4[3] * (P1[3] - P1[0]) + F2[3]
      * F4[2] * (P1[0] - P1[3])) + F3[5] * (F2[2] * F4[3] * (P1[1] - cI *
      (P1[2])) + F2[3] * F4[2] * (+cI * (P1[2]) - P1[1])));
  F1[4] = denom * 2. * cI * F3[5] * M1 * (F2[2] * F4[3] - F2[3] * F4[2]); 
  F1[5] = denom * 2. * cI * F3[4] * M1 * (F2[3] * F4[2] - F2[2] * F4[3]); 
}


void FFFF1_1(complex<double> F2[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + F3[0] + F4[0]; 
  F1[1] = +F2[1] + F3[1] + F4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * 2. * cI * F3[3] * M1 * (F2[5] * F4[4] - F2[4] * F4[5]); 
  F1[3] = denom * 2. * cI * F3[2] * M1 * (F2[4] * F4[5] - F2[5] * F4[4]); 
  F1[4] = denom * - 2. * cI * (F3[2] * (F2[4] * - F4[5] * (P1[1] + cI *
      (P1[2])) + F2[5] * F4[4] * (P1[1] + cI * (P1[2]))) + F3[3] * (F2[4] *
      F4[5] * (P1[3] - P1[0]) + F2[5] * F4[4] * (P1[0] - P1[3])));
  F1[5] = denom * - 2. * cI * (F3[2] * (F2[4] * F4[5] * (P1[0] + P1[3]) - F2[5]
      * F4[4] * (P1[0] + P1[3])) + F3[3] * (F2[4] * F4[5] * (P1[1] - cI *
      (P1[2])) + F2[5] * F4[4] * (+cI * (P1[2]) - P1[1])));
}

void FFFF1_2_6_1(complex<double> F2[], complex<double> F3[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFFF1_1(F2, F3, F4, COUP1, M1, W1, F1); 
  FFFF2_1(F2, F3, F4, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
  FFFF6_1(F2, F3, F4, COUP3, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}
void FFFF1_6_1(complex<double> F2[], complex<double> F3[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M1, double W1,
    complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  complex<double> COUP3; 
  int i; 
  complex<double> Ftmp[6]; 
  FFFF1_1(F2, F3, F4, COUP1, M1, W1, F1); 
  FFFF6_1(F2, F3, F4, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFS2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP39; 
  double P3[4]; 
  complex<double> TMP43; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP43 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  TMP39 = (F1[2] * F2[2] + F1[3] * F2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * (+cI * (TMP39 + TMP43)); 
}


void VVVV5P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V3[2] * TMP11) + cI * (V2[2] * TMP22)); 
  V1[3] = denom * (-cI * (V3[3] * TMP11) + cI * (V2[3] * TMP22)); 
  V1[4] = denom * (-cI * (V3[4] * TMP11) + cI * (V2[4] * TMP22)); 
  V1[5] = denom * (-cI * (V3[5] * TMP11) + cI * (V2[5] * TMP22)); 
}


void FFV2_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))));
  F2[4] = denom * - cI * M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]));
  F2[5] = denom * cI * M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] * (V3[2] -
      V3[5]));
}

void FFV2_7_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  double P2[4]; 
  complex<double> denom; 
  int i; 
  FFV2_2(F1, V3, COUP1, M2, W2, F2); 
  FFV7_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_4_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  double P2[4]; 
  complex<double> denom; 
  int i; 
  FFV2_2(F1, V3, COUP1, M2, W2, F2); 
  FFV4_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFV3_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP5; 
  double P3[4]; 
  complex<double> TMP6; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP5 = (F1[2] * (F2[2] * (P3[0] * (V3[2] + V3[5]) + (P3[1] * - 1. * (V3[3] +
      cI * (V3[4])) + (P3[2] * (+cI * (V3[3]) - V3[4]) - P3[3] * (V3[2] +
      V3[5])))) + F2[3] * (P3[0] * (V3[3] + cI * (V3[4])) + (P3[1] * - 1. *
      (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) + P3[3] *
      (V3[3] + cI * (V3[4])))))) + F1[3] * (F2[2] * (P3[0] * (V3[3] - cI *
      (V3[4])) + (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (+cI * (V3[4]) - V3[3])))) + F2[3] * (P3[0] * (V3[2] -
      V3[5]) + (P3[1] * (+cI * (V3[4]) - V3[3]) + (P3[2] * - 1. * (V3[4] + cI *
      (V3[3])) + P3[3] * (V3[2] - V3[5]))))));
  TMP6 = (F1[2] * (F2[2] * (P3[0] * (V3[2] - V3[5]) + (P3[1] * (+cI * (V3[4]) -
      V3[3]) + (P3[2] * - 1. * (V3[4] + cI * (V3[3])) + P3[3] * (V3[2] -
      V3[5])))) + F2[3] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) + (P3[1] *
      (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) - P3[3] * (V3[3] + cI
      * (V3[4])))))) + F1[3] * (F2[2] * (P3[0] * (+cI * (V3[4]) - V3[3]) +
      (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[2]) + cI * (V3[5])) +
      P3[3] * (V3[3] - cI * (V3[4]))))) + F2[3] * (P3[0] * (V3[2] + V3[5]) +
      (P3[1] * - 1. * (V3[3] + cI * (V3[4])) + (P3[2] * (+cI * (V3[3]) - V3[4])
      - P3[3] * (V3[2] + V3[5]))))));
  vertex = COUP * (-cI * (TMP6) + cI * (TMP5)); 
}

void FFV3_8_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> tmp; 
  FFV3_0(F1, F2, V3, COUP1, vertex); 
  FFV8_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void VVVVV13P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP32; 
  double P4[4]; 
  complex<double> TMP42; 
  complex<double> TMP25; 
  complex<double> TMP18; 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP42 = (V5[2] * P4[0] - V5[3] * P4[1] - V5[4] * P4[2] - V5[5] * P4[3]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P4[0] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP11 * (-cI * (V5[2] * TMP25) + cI * (V3[2] * TMP42)) + (TMP32 * (-cI *
      (V3[2] * TMP21) + cI * (V5[2] * TMP22)) + V4[2] * (-cI * (TMP10 * TMP42)
      + cI * (TMP18 * TMP25)))));
  V1[3] = denom * (P4[1] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP11 * (-cI * (V5[3] * TMP25) + cI * (V3[3] * TMP42)) + (TMP32 * (-cI *
      (V3[3] * TMP21) + cI * (V5[3] * TMP22)) + V4[3] * (-cI * (TMP10 * TMP42)
      + cI * (TMP18 * TMP25)))));
  V1[4] = denom * (P4[2] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP11 * (-cI * (V5[4] * TMP25) + cI * (V3[4] * TMP42)) + (TMP32 * (-cI *
      (V3[4] * TMP21) + cI * (V5[4] * TMP22)) + V4[4] * (-cI * (TMP10 * TMP42)
      + cI * (TMP18 * TMP25)))));
  V1[5] = denom * (P4[3] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP11 * (-cI * (V5[5] * TMP25) + cI * (V3[5] * TMP42)) + (TMP32 * (-cI *
      (V3[5] * TMP21) + cI * (V5[5] * TMP22)) + V4[5] * (-cI * (TMP10 * TMP42)
      + cI * (TMP18 * TMP25)))));
}


void FFFF2_2(complex<double> F1[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + F3[0] + F4[0]; 
  F2[1] = +F1[1] + F3[1] + F4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * 2. * cI * (F4[4] * (F1[2] * F3[3] * (+cI * (P2[2]) - P2[1]) +
      F1[3] * F3[2] * (P2[1] - cI * (P2[2]))) + F4[5] * (F1[2] * F3[3] * (P2[3]
      - P2[0]) + F1[3] * F3[2] * (P2[0] - P2[3])));
  F2[3] = denom * 2. * cI * (F4[4] * (F1[2] * F3[3] * (P2[0] + P2[3]) - F1[3] *
      F3[2] * (P2[0] + P2[3])) + F4[5] * (F1[2] * F3[3] * (P2[1] + cI *
      (P2[2])) - F1[3] * F3[2] * (P2[1] + cI * (P2[2]))));
  F2[4] = denom * 2. * cI * F4[5] * M2 * (F1[3] * F3[2] - F1[2] * F3[3]); 
  F2[5] = denom * 2. * cI * F4[4] * M2 * (F1[2] * F3[3] - F1[3] * F3[2]); 
}


void VVVV2P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP10; 
  double P3[4]; 
  complex<double> TMP47; 
  complex<double> TMP32; 
  complex<double> TMP28; 
  complex<double> TMP25; 
  complex<double> TMP11; 
  double P2[4]; 
  complex<double> TMP53; 
  complex<double> TMP24; 
  double P1[4]; 
  complex<double> TMP30; 
  complex<double> denom; 
  complex<double> TMP27; 
  double P4[4]; 
  complex<double> TMP22; 
  complex<double> TMP31; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP47 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP53 = (P3[0] * P4[0] - P3[1] * P4[1] - P3[2] * P4[2] - P3[3] * P4[3]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP10 * (V4[2] * - 1. * (+cI * (TMP47 + TMP53)) + (+cI *
      (P2[0] * TMP28 + P4[0] * TMP26))) + (TMP11 * (V3[2] * (+cI * (TMP47 +
      TMP53)) + (-cI * (P2[0] * TMP27 + P3[0] * TMP25))) + (TMP22 * (-cI *
      (P4[0] * TMP29) + cI * (P3[0] * TMP32)) + (TMP24 * (-cI * (V3[2] * TMP30)
      + cI * (V2[2] * TMP27)) + (TMP31 * (-cI * (V2[2] * TMP28) + cI * (V4[2] *
      TMP30)) + (-cI * (V3[2] * TMP26 * TMP32) + cI * (V4[2] * TMP25 *
      TMP29)))))));
  V1[3] = denom * (TMP10 * (V4[3] * - 1. * (+cI * (TMP47 + TMP53)) + (+cI *
      (P2[1] * TMP28 + P4[1] * TMP26))) + (TMP11 * (V3[3] * (+cI * (TMP47 +
      TMP53)) + (-cI * (P2[1] * TMP27 + P3[1] * TMP25))) + (TMP22 * (-cI *
      (P4[1] * TMP29) + cI * (P3[1] * TMP32)) + (TMP24 * (-cI * (V3[3] * TMP30)
      + cI * (V2[3] * TMP27)) + (TMP31 * (-cI * (V2[3] * TMP28) + cI * (V4[3] *
      TMP30)) + (-cI * (V3[3] * TMP26 * TMP32) + cI * (V4[3] * TMP25 *
      TMP29)))))));
  V1[4] = denom * (TMP10 * (V4[4] * - 1. * (+cI * (TMP47 + TMP53)) + (+cI *
      (P2[2] * TMP28 + P4[2] * TMP26))) + (TMP11 * (V3[4] * (+cI * (TMP47 +
      TMP53)) + (-cI * (P2[2] * TMP27 + P3[2] * TMP25))) + (TMP22 * (-cI *
      (P4[2] * TMP29) + cI * (P3[2] * TMP32)) + (TMP24 * (-cI * (V3[4] * TMP30)
      + cI * (V2[4] * TMP27)) + (TMP31 * (-cI * (V2[4] * TMP28) + cI * (V4[4] *
      TMP30)) + (-cI * (V3[4] * TMP26 * TMP32) + cI * (V4[4] * TMP25 *
      TMP29)))))));
  V1[5] = denom * (TMP10 * (V4[5] * - 1. * (+cI * (TMP47 + TMP53)) + (+cI *
      (P2[3] * TMP28 + P4[3] * TMP26))) + (TMP11 * (V3[5] * (+cI * (TMP47 +
      TMP53)) + (-cI * (P2[3] * TMP27 + P3[3] * TMP25))) + (TMP22 * (-cI *
      (P4[3] * TMP29) + cI * (P3[3] * TMP32)) + (TMP24 * (-cI * (V3[5] * TMP30)
      + cI * (V2[5] * TMP27)) + (TMP31 * (-cI * (V2[5] * TMP28) + cI * (V4[5] *
      TMP30)) + (-cI * (V3[5] * TMP26 * TMP32) + cI * (V4[5] * TMP25 *
      TMP29)))))));
}


void FFFF4_2(complex<double> F1[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  complex<double> TMP48; 
  F2[0] = +F1[0] + F3[0] + F4[0]; 
  F2[1] = +F1[1] + F3[1] + F4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  TMP48 = (F4[2] * F3[2] + F4[3] * F3[3]); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * TMP48 * (F1[4] * (P2[0] - P2[3]) + F1[5] * (+cI *
      (P2[2]) - P2[1]));
  F2[3] = denom * - cI * TMP48 * (F1[4] * (P2[1] + cI * (P2[2])) - F1[5] *
      (P2[0] + P2[3]));
  F2[4] = denom * cI * TMP48 * F1[4] * M2; 
  F2[5] = denom * cI * TMP48 * F1[5] * M2; 
}

void FFFF4_5_2(complex<double> F1[], complex<double> F3[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M2, double W2,
    complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  double P2[4]; 
  complex<double> denom; 
  int i; 
  FFFF4_2(F1, F3, F4, COUP1, M2, W2, F2); 
  FFFF5_2(F1, F3, F4, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFV3P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - 2. * cI * (F1[2] * (F2[3] * (P3[1] + cI * (P3[2])) + P3[3]
      * F2[2]) + F1[3] * (F2[2] * (P3[1] - cI * (P3[2])) - P3[3] * F2[3]));
  V3[3] = denom * - 2. * cI * (F1[2] * (F2[3] * (P3[0] + P3[3]) + cI * (P3[2] *
      F2[2])) + F1[3] * (F2[2] * (P3[0] - P3[3]) - cI * (P3[2] * F2[3])));
  V3[4] = denom * 2. * (F1[2] * (F2[3] * (P3[0] + P3[3]) - P3[1] * F2[2]) +
      F1[3] * (F2[2] * (P3[3] - P3[0]) + P3[1] * F2[3]));
  V3[5] = denom * 2. * cI * (F1[2] * (F2[3] * (P3[1] + cI * (P3[2])) - P3[0] *
      F2[2]) + F1[3] * (F2[2] * (+cI * (P3[2]) - P3[1]) + P3[0] * F2[3]));
}

void FFV3_8P0_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Vtmp[6]; 
  FFV3P0_3(F1, F2, COUP1, M3, W3, V3); 
  FFV8P0_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void FFFF6_3(complex<double> F1[], complex<double> F2[], complex<double> F4[],
    complex<double> COUP, double M3, double W3, complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  F3[0] = +F1[0] + F2[0] + F4[0]; 
  F3[1] = +F1[1] + F2[1] + F4[1]; 
  P3[0] = -F3[0].real(); 
  P3[1] = -F3[1].real(); 
  P3[2] = -F3[1].imag(); 
  P3[3] = -F3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  F3[2] = denom * - 2. * cI * (F1[4] * (F2[2] * - F4[3] * (P3[1] + cI *
      (P3[2])) + F2[3] * F4[2] * (P3[1] + cI * (P3[2]))) + F1[5] * (F2[2] *
      F4[3] * (P3[0] + P3[3]) - F2[3] * F4[2] * (P3[0] + P3[3])));
  F3[3] = denom * - 2. * cI * (F1[4] * (F2[2] * F4[3] * (P3[3] - P3[0]) + F2[3]
      * F4[2] * (P3[0] - P3[3])) + F1[5] * (F2[2] * F4[3] * (P3[1] - cI *
      (P3[2])) + F2[3] * F4[2] * (+cI * (P3[2]) - P3[1])));
  F3[4] = denom * 2. * cI * F1[5] * M3 * (F2[2] * F4[3] - F2[3] * F4[2]); 
  F3[5] = denom * 2. * cI * F1[4] * M3 * (F2[3] * F4[2] - F2[2] * F4[3]); 
}


void FFV4_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - 2. * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + (+1./2. * (M1 * (+2. * (F2[4] * - 1./2. * (V3[2] + V3[5]))
      - F2[5] * (V3[3] + cI * (V3[4])))) + F2[3] * (P1[0] * (V3[3] + cI *
      (V3[4])) + (P1[1] * - 1. * (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI *
      (V3[2] + V3[5])) + P1[3] * (V3[3] + cI * (V3[4])))))));
  F1[3] = denom * - 2. * cI * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1]
      * (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] *
      (+cI * (V3[4]) - V3[3])))) + (+1./2. * (M1 * (F2[5] * (V3[5] - V3[2]) +
      2. * (F2[4] * 1./2. * (+cI * (V3[4]) - V3[3])))) + F2[3] * (P1[0] * - 1.
      * (V3[2] + V3[5]) + (P1[1] * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] -
      cI * (V3[3])) + P1[3] * (V3[2] + V3[5]))))));
  F1[4] = denom * cI * (F2[4] * (P1[0] * - 1. * (V3[2] + V3[5]) + (P1[1] *
      (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[2]
      + V3[5])))) + (F2[5] * (P1[0] * - 1. * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 2. * (V3[5] - V3[2]) + 2. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * -
      1. * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3]
      - cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))) + M1 * (F2[2] * 2. * (+cI * (V3[4]) - V3[3]) + 2. *
      (F2[3] * (V3[2] + V3[5])))));
}


void FFFF1_4(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> COUP, double M4, double W4, complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P4[4]; 
  F4[0] = +F1[0] + F2[0] + F3[0]; 
  F4[1] = +F1[1] + F2[1] + F3[1]; 
  P4[0] = -F4[0].real(); 
  P4[1] = -F4[1].real(); 
  P4[2] = -F4[1].imag(); 
  P4[3] = -F4[0].imag(); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  F4[2] = denom * 2. * cI * (F2[4] * (F1[2] * F3[3] * (P4[1] - cI * (P4[2])) +
      F1[3] * F3[2] * (+cI * (P4[2]) - P4[1])) + F2[5] * (F1[2] * F3[3] *
      (P4[0] - P4[3]) + F1[3] * F3[2] * (P4[3] - P4[0])));
  F4[3] = denom * 2. * cI * (F2[4] * (F1[2] * - F3[3] * (P4[0] + P4[3]) + F1[3]
      * F3[2] * (P4[0] + P4[3])) + F2[5] * (F1[2] * - F3[3] * (P4[1] + cI *
      (P4[2])) + F1[3] * F3[2] * (P4[1] + cI * (P4[2]))));
  F4[4] = denom * 2. * cI * F2[5] * M4 * (F1[2] * F3[3] - F1[3] * F3[2]); 
  F4[5] = denom * 2. * cI * F2[4] * M4 * (F1[3] * F3[2] - F1[2] * F3[3]); 
}

void FFFF1_6_4(complex<double> F1[], complex<double> F2[], complex<double>
    F3[], complex<double> COUP1, complex<double> COUP2, double M4, double W4,
    complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  complex<double> denom; 
  double P4[4]; 
  int i; 
  FFFF1_4(F1, F2, F3, COUP1, M4, W4, F4); 
  FFFF6_4(F1, F2, F3, COUP2, M4, W4, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F4[i] = F4[i] + Ftmp[i]; 
    i++; 
  }
}
void FFFF1_2_6_4(complex<double> F1[], complex<double> F2[], complex<double>
    F3[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    double M4, double W4, complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  complex<double> denom; 
  double P4[4]; 
  int i; 
  FFFF1_4(F1, F2, F3, COUP1, M4, W4, F4); 
  FFFF2_4(F1, F2, F3, COUP2, M4, W4, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F4[i] = F4[i] + Ftmp[i]; 
    i++; 
  }
  FFFF6_4(F1, F2, F3, COUP3, M4, W4, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F4[i] = F4[i] + Ftmp[i]; 
    i++; 
  }
}

void VVVVV8P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP31; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP40; 
  complex<double> TMP24; 
  complex<double> TMP18; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP40 = (P2[0] * V5[2] - P2[1] * V5[3] - P2[2] * V5[4] - P2[3] * V5[5]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P2[0] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP11 * (-cI * (V3[2] * TMP40) + cI * (V5[2] * TMP31)) + (TMP24 * (-cI *
      (TMP10 * V5[2]) + cI * (V3[2] * TMP18)) + V2[2] * (-cI * (TMP21 * TMP31)
      + cI * (TMP22 * TMP40)))));
  V1[3] = denom * (P2[1] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP11 * (-cI * (V3[3] * TMP40) + cI * (V5[3] * TMP31)) + (TMP24 * (-cI *
      (TMP10 * V5[3]) + cI * (V3[3] * TMP18)) + V2[3] * (-cI * (TMP21 * TMP31)
      + cI * (TMP22 * TMP40)))));
  V1[4] = denom * (P2[2] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP11 * (-cI * (V3[4] * TMP40) + cI * (V5[4] * TMP31)) + (TMP24 * (-cI *
      (TMP10 * V5[4]) + cI * (V3[4] * TMP18)) + V2[4] * (-cI * (TMP21 * TMP31)
      + cI * (TMP22 * TMP40)))));
  V1[5] = denom * (P2[3] * (-cI * (TMP18 * TMP22) + cI * (TMP10 * TMP21)) +
      (TMP11 * (-cI * (V3[5] * TMP40) + cI * (V5[5] * TMP31)) + (TMP24 * (-cI *
      (TMP10 * V5[5]) + cI * (V3[5] * TMP18)) + V2[5] * (-cI * (TMP21 * TMP31)
      + cI * (TMP22 * TMP40)))));
}


void VVVVS3_5(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, double M5, double W5,
    complex<double> S5[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  complex<double> TMP22; 
  complex<double> TMP23; 
  complex<double> denom; 
  double P5[4]; 
  S5[0] = +V1[0] + V2[0] + V3[0] + V4[0]; 
  S5[1] = +V1[1] + V2[1] + V3[1] + V4[1]; 
  P5[0] = -S5[0].real(); 
  P5[1] = -S5[1].real(); 
  P5[2] = -S5[1].imag(); 
  P5[3] = -S5[0].imag(); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/(pow(P5[0], 2) - pow(P5[1], 2) - pow(P5[2], 2) - pow(P5[3], 2) -
      M5 * (M5 - cI * W5));
  S5[2] = denom * (-cI * (TMP22 * TMP23) + cI * (TMP11 * TMP12)); 
}


void VVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP31; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> denom; 
  complex<double> TMP27; 
  complex<double> TMP29; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP10 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI *
      (TMP27) + cI * (TMP31)) + V3[2] * (-cI * (TMP29) + cI * (TMP30))));
  V1[3] = denom * (TMP10 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI *
      (TMP27) + cI * (TMP31)) + V3[3] * (-cI * (TMP29) + cI * (TMP30))));
  V1[4] = denom * (TMP10 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI *
      (TMP27) + cI * (TMP31)) + V3[4] * (-cI * (TMP29) + cI * (TMP30))));
  V1[5] = denom * (TMP10 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI *
      (TMP27) + cI * (TMP31)) + V3[5] * (-cI * (TMP29) + cI * (TMP30))));
}


void FFFF7_2(complex<double> F1[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + F3[0] + F4[0]; 
  F2[1] = +F1[1] + F3[1] + F4[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * 2. * cI * F4[3] * M2 * (F1[4] * F3[5] - F1[5] * F3[4]); 
  F2[3] = denom * 2. * cI * F4[2] * M2 * (F1[5] * F3[4] - F1[4] * F3[5]); 
  F2[4] = denom * 2. * cI * (F4[2] * (F1[4] * F3[5] * (+cI * (P2[2]) - P2[1]) +
      F1[5] * F3[4] * (P2[1] - cI * (P2[2]))) + F4[3] * (F1[4] * F3[5] * (P2[0]
      + P2[3]) - F1[5] * F3[4] * (P2[0] + P2[3])));
  F2[5] = denom * 2. * cI * (F4[2] * (F1[4] * F3[5] * (P2[3] - P2[0]) + F1[5] *
      F3[4] * (P2[0] - P2[3])) + F4[3] * (F1[4] * F3[5] * (P2[1] + cI *
      (P2[2])) - F1[5] * F3[4] * (P2[1] + cI * (P2[2]))));
}


void FFV3_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP7; 
  double P3[4]; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP7 = (F1[2] * - F2[2] * (P3[1] * P3[1] + P3[2] * P3[2] + P3[3] * P3[3] -
      P3[0] * P3[0]) - F1[3] * F2[3] * (P3[1] * P3[1] + P3[2] * P3[2] + P3[3] *
      P3[3] - P3[0] * P3[0]));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - 2. * cI * (F1[2] * (F2[3] * (P3[1] + cI * (P3[2])) + P3[3]
      * F2[2]) + F1[3] * (F2[2] * (P3[1] - cI * (P3[2])) - P3[3] * F2[3]));
  V3[3] = denom * - 2. * cI * (F1[2] * (F2[3] * (P3[0] + P3[3]) + cI * (P3[2] *
      F2[2])) + F1[3] * (F2[2] * (P3[0] - P3[3]) - cI * (P3[2] * F2[3])));
  V3[4] = denom * 2. * (F1[2] * (F2[3] * (P3[0] + P3[3]) - P3[1] * F2[2]) +
      F1[3] * (F2[2] * (P3[3] - P3[0]) + P3[1] * F2[3]));
  V3[5] = denom * 2. * cI * (F1[2] * (F2[3] * (P3[1] + cI * (P3[2])) - P3[0] *
      F2[2]) + F1[3] * (F2[2] * (+cI * (P3[2]) - P3[1]) + P3[0] * F2[3]));
}

void FFV3_8_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  int i; 
  complex<double> Vtmp[6]; 
  FFV3_3(F1, F2, COUP1, M3, W3, V3); 
  FFV8_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void FFV7_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * 4. * (V3[2] - V3[5]) + 4. * (F1[5]
      * (+cI * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))) + M2 * (F1[4] * - 4. * (V3[3] + cI * (V3[4])) + 4. *
      (F1[5] * (V3[2] + V3[5])))));
  F2[4] = denom * - 4. * cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + (+1./4. * (M2 * (F1[3] * (+cI * (V3[4]) - V3[3]) + 4. *
      (F1[2] * - 1./4. * (V3[2] + V3[5])))) + F1[5] * (P2[0] * (V3[3] - cI *
      (V3[4])) + (P2[1] * - 1. * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] +
      V3[5])) + P2[3] * (V3[3] - cI * (V3[4])))))));
  F2[5] = denom * - 4. * cI * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./4. * (M2 * (F1[3] * (V3[5] - V3[2]) + 4.
      * (F1[2] * - 1./4. * (V3[3] + cI * (V3[4]))))) + F1[5] * (P2[0] * - 1. *
      (V3[2] + V3[5]) + (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI
      * (V3[3])) + P2[3] * (V3[2] + V3[5]))))));
}


void VVVVV3P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP45; 
  double P3[4]; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP45 = (P3[0] * V5[2] - P3[1] * V5[3] - P3[2] * V5[4] - P3[3] * V5[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P3[0] * (-cI * (TMP11 * TMP19) + cI * (TMP18 * TMP22)) +
      (TMP10 * (-cI * (V4[2] * TMP45) + cI * (V5[2] * TMP26)) + (TMP29 * (-cI *
      (V5[2] * TMP22) + cI * (V4[2] * TMP19)) + V3[2] * (-cI * (TMP18 * TMP26)
      + cI * (TMP11 * TMP45)))));
  V1[3] = denom * (P3[1] * (-cI * (TMP11 * TMP19) + cI * (TMP18 * TMP22)) +
      (TMP10 * (-cI * (V4[3] * TMP45) + cI * (V5[3] * TMP26)) + (TMP29 * (-cI *
      (V5[3] * TMP22) + cI * (V4[3] * TMP19)) + V3[3] * (-cI * (TMP18 * TMP26)
      + cI * (TMP11 * TMP45)))));
  V1[4] = denom * (P3[2] * (-cI * (TMP11 * TMP19) + cI * (TMP18 * TMP22)) +
      (TMP10 * (-cI * (V4[4] * TMP45) + cI * (V5[4] * TMP26)) + (TMP29 * (-cI *
      (V5[4] * TMP22) + cI * (V4[4] * TMP19)) + V3[4] * (-cI * (TMP18 * TMP26)
      + cI * (TMP11 * TMP45)))));
  V1[5] = denom * (P3[3] * (-cI * (TMP11 * TMP19) + cI * (TMP18 * TMP22)) +
      (TMP10 * (-cI * (V4[5] * TMP45) + cI * (V5[5] * TMP26)) + (TMP29 * (-cI *
      (V5[5] * TMP22) + cI * (V4[5] * TMP19)) + V3[5] * (-cI * (TMP18 * TMP26)
      + cI * (TMP11 * TMP45)))));
}


void VVVV8P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP36; 
  double P3[4]; 
  complex<double> TMP32; 
  complex<double> TMP28; 
  complex<double> TMP25; 
  complex<double> TMP11; 
  double P2[4]; 
  complex<double> TMP10; 
  complex<double> TMP24; 
  double P1[4]; 
  complex<double> TMP30; 
  complex<double> denom; 
  complex<double> TMP27; 
  double P4[4]; 
  complex<double> TMP35; 
  complex<double> TMP22; 
  complex<double> TMP31; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP36 = (P2[0] * P4[0] - P2[1] * P4[1] - P2[2] * P4[2] - P2[3] * P4[3]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP35 = (P3[0] * P1[0] - P3[1] * P1[1] - P3[2] * P1[2] - P3[3] * P1[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP10 * (V4[2] * - 1. * (+cI * (TMP35 + TMP36)) + (+cI *
      (P3[0] * TMP28 + P4[0] * TMP24))) + (TMP22 * (V2[2] * (+cI * (TMP35 +
      TMP36)) + (-cI * (P3[0] * TMP30 + P2[0] * TMP32))) + (TMP11 * (-cI *
      (P4[0] * TMP31) + cI * (P2[0] * TMP25)) + (TMP26 * (-cI * (V2[2] * TMP27)
      + cI * (V3[2] * TMP30)) + (TMP29 * (-cI * (V3[2] * TMP28) + cI * (V4[2] *
      TMP27)) + (-cI * (V2[2] * TMP24 * TMP25) + cI * (V4[2] * TMP31 *
      TMP32)))))));
  V1[3] = denom * (TMP10 * (V4[3] * - 1. * (+cI * (TMP35 + TMP36)) + (+cI *
      (P3[1] * TMP28 + P4[1] * TMP24))) + (TMP22 * (V2[3] * (+cI * (TMP35 +
      TMP36)) + (-cI * (P3[1] * TMP30 + P2[1] * TMP32))) + (TMP11 * (-cI *
      (P4[1] * TMP31) + cI * (P2[1] * TMP25)) + (TMP26 * (-cI * (V2[3] * TMP27)
      + cI * (V3[3] * TMP30)) + (TMP29 * (-cI * (V3[3] * TMP28) + cI * (V4[3] *
      TMP27)) + (-cI * (V2[3] * TMP24 * TMP25) + cI * (V4[3] * TMP31 *
      TMP32)))))));
  V1[4] = denom * (TMP10 * (V4[4] * - 1. * (+cI * (TMP35 + TMP36)) + (+cI *
      (P3[2] * TMP28 + P4[2] * TMP24))) + (TMP22 * (V2[4] * (+cI * (TMP35 +
      TMP36)) + (-cI * (P3[2] * TMP30 + P2[2] * TMP32))) + (TMP11 * (-cI *
      (P4[2] * TMP31) + cI * (P2[2] * TMP25)) + (TMP26 * (-cI * (V2[4] * TMP27)
      + cI * (V3[4] * TMP30)) + (TMP29 * (-cI * (V3[4] * TMP28) + cI * (V4[4] *
      TMP27)) + (-cI * (V2[4] * TMP24 * TMP25) + cI * (V4[4] * TMP31 *
      TMP32)))))));
  V1[5] = denom * (TMP10 * (V4[5] * - 1. * (+cI * (TMP35 + TMP36)) + (+cI *
      (P3[3] * TMP28 + P4[3] * TMP24))) + (TMP22 * (V2[5] * (+cI * (TMP35 +
      TMP36)) + (-cI * (P3[3] * TMP30 + P2[3] * TMP32))) + (TMP11 * (-cI *
      (P4[3] * TMP31) + cI * (P2[3] * TMP25)) + (TMP26 * (-cI * (V2[5] * TMP27)
      + cI * (V3[5] * TMP30)) + (TMP29 * (-cI * (V3[5] * TMP28) + cI * (V4[5] *
      TMP27)) + (-cI * (V2[5] * TMP24 * TMP25) + cI * (V4[5] * TMP31 *
      TMP32)))))));
}


void FFFF3_4(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> COUP, double M4, double W4, complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P4[4]; 
  F4[0] = +F1[0] + F2[0] + F3[0]; 
  F4[1] = +F1[1] + F2[1] + F3[1]; 
  P4[0] = -F4[0].real(); 
  P4[1] = -F4[1].real(); 
  P4[2] = -F4[1].imag(); 
  P4[3] = -F4[0].imag(); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  F4[2] = denom * 2. * cI * (F2[4] * (F1[2] * F3[3] * (+cI * (P4[2]) - P4[1]) +
      F1[3] * F3[2] * (P4[1] - cI * (P4[2]))) + F2[5] * (F1[2] * F3[3] * (P4[3]
      - P4[0]) + F1[3] * F3[2] * (P4[0] - P4[3])));
  F4[3] = denom * 2. * cI * (F2[4] * (F1[2] * F3[3] * (P4[0] + P4[3]) - F1[3] *
      F3[2] * (P4[0] + P4[3])) + F2[5] * (F1[2] * F3[3] * (P4[1] + cI *
      (P4[2])) - F1[3] * F3[2] * (P4[1] + cI * (P4[2]))));
  F4[4] = denom * 2. * cI * F2[5] * M4 * (F1[3] * F3[2] - F1[2] * F3[3]); 
  F4[5] = denom * 2. * cI * F2[4] * M4 * (F1[2] * F3[3] - F1[3] * F3[2]); 
}

void FFFF3_7_4(complex<double> F1[], complex<double> F2[], complex<double>
    F3[], complex<double> COUP1, complex<double> COUP2, double M4, double W4,
    complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
  complex<double> denom; 
  double P4[4]; 
  int i; 
  FFFF3_4(F1, F2, F3, COUP1, M4, W4, F4); 
  FFFF7_4(F1, F2, F3, COUP2, M4, W4, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F4[i] = F4[i] + Ftmp[i]; 
    i++; 
  }
}

void FFFF5_0(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> F4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP39; 
  complex<double> TMP38; 
  TMP39 = (F1[2] * F2[2] + F1[3] * F2[3]); 
  TMP38 = (F4[4] * F3[4] + F4[5] * F3[5]); 
  vertex = COUP * - cI * TMP38 * TMP39; 
}


void FFFF5_3(complex<double> F1[], complex<double> F2[], complex<double> F4[],
    complex<double> COUP, double M3, double W3, complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP39; 
  double P3[4]; 
  F3[0] = +F1[0] + F2[0] + F4[0]; 
  F3[1] = +F1[1] + F2[1] + F4[1]; 
  P3[0] = -F3[0].real(); 
  P3[1] = -F3[1].real(); 
  P3[2] = -F3[1].imag(); 
  P3[3] = -F3[0].imag(); 
  TMP39 = (F1[2] * F2[2] + F1[3] * F2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  F3[2] = denom * - cI * TMP39 * (F4[4] * (P3[0] + P3[3]) + F4[5] * (P3[1] + cI
      * (P3[2])));
  F3[3] = denom * cI * TMP39 * (F4[4] * (+cI * (P3[2]) - P3[1]) + F4[5] *
      (P3[3] - P3[0]));
  F3[4] = denom * cI * F4[4] * TMP39 * M3; 
  F3[5] = denom * cI * F4[5] * TMP39 * M3; 
}


void VVVVS1_5(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, double M5, double W5,
    complex<double> S5[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  complex<double> TMP10; 
  complex<double> denom; 
  double P5[4]; 
  complex<double> TMP9; 
  S5[0] = +V1[0] + V2[0] + V3[0] + V4[0]; 
  S5[1] = +V1[1] + V2[1] + V3[1] + V4[1]; 
  P5[0] = -S5[0].real(); 
  P5[1] = -S5[1].real(); 
  P5[2] = -S5[1].imag(); 
  P5[3] = -S5[0].imag(); 
  TMP9 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/(pow(P5[0], 2) - pow(P5[1], 2) - pow(P5[2], 2) - pow(P5[3], 2) -
      M5 * (M5 - cI * W5));
  S5[2] = denom * (-cI * (TMP11 * TMP12) + cI * (TMP9 * TMP10)); 
}


void FFFF2_1(complex<double> F2[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + F3[0] + F4[0]; 
  F1[1] = +F2[1] + F3[1] + F4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * 2. * cI * F3[3] * M1 * (F2[5] * F4[4] - F2[4] * F4[5]); 
  F1[3] = denom * 2. * cI * F3[2] * M1 * (F2[4] * F4[5] - F2[5] * F4[4]); 
  F1[4] = denom * - 2. * cI * (F3[2] * (F2[4] * - F4[5] * (P1[1] + cI *
      (P1[2])) + F2[5] * F4[4] * (P1[1] + cI * (P1[2]))) + F3[3] * (F2[4] *
      F4[5] * (P1[3] - P1[0]) + F2[5] * F4[4] * (P1[0] - P1[3])));
  F1[5] = denom * - 2. * cI * (F3[2] * (F2[4] * F4[5] * (P1[0] + P1[3]) - F2[5]
      * F4[4] * (P1[0] + P1[3])) + F3[3] * (F2[4] * F4[5] * (P1[1] - cI *
      (P1[2])) + F2[5] * F4[4] * (+cI * (P1[2]) - P1[1])));
}


void FFV1_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP8; 
  TMP8 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  vertex = COUP * - cI * TMP8; 
}


void FFFF2_4(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> COUP, double M4, double W4, complex<double> F4[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P4[4]; 
  F4[0] = +F1[0] + F2[0] + F3[0]; 
  F4[1] = +F1[1] + F2[1] + F3[1]; 
  P4[0] = -F4[0].real(); 
  P4[1] = -F4[1].real(); 
  P4[2] = -F4[1].imag(); 
  P4[3] = -F4[0].imag(); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  F4[2] = denom * 2. * cI * (F2[4] * (F1[2] * F3[3] * (P4[1] - cI * (P4[2])) +
      F1[3] * F3[2] * (+cI * (P4[2]) - P4[1])) + F2[5] * (F1[2] * F3[3] *
      (P4[0] - P4[3]) + F1[3] * F3[2] * (P4[3] - P4[0])));
  F4[3] = denom * 2. * cI * (F2[4] * (F1[2] * - F3[3] * (P4[0] + P4[3]) + F1[3]
      * F3[2] * (P4[0] + P4[3])) + F2[5] * (F1[2] * - F3[3] * (P4[1] + cI *
      (P4[2])) + F1[3] * F3[2] * (P4[1] + cI * (P4[2]))));
  F4[4] = denom * 2. * cI * F2[5] * M4 * (F1[2] * F3[3] - F1[3] * F3[2]); 
  F4[5] = denom * 2. * cI * F2[4] * M4 * (F1[3] * F3[2] - F1[2] * F3[3]); 
}


void FFV8_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP2 = (F1[4] * - F2[4] * (P3[1] * P3[1] + P3[2] * P3[2] + P3[3] * P3[3] -
      P3[0] * P3[0]) - F1[5] * F2[5] * (P3[1] * P3[1] + P3[2] * P3[2] + P3[3] *
      P3[3] - P3[0] * P3[0]));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * 2. * cI * (F1[4] * (F2[5] * (P3[1] + cI * (P3[2])) + P3[3] *
      F2[4]) + F1[5] * (F2[4] * (P3[1] - cI * (P3[2])) - P3[3] * F2[5]));
  V3[3] = denom * 2. * cI * (F1[4] * (F2[5] * (P3[0] - P3[3]) - cI * (P3[2] *
      F2[4])) + F1[5] * (F2[4] * (P3[0] + P3[3]) + cI * (P3[2] * F2[5])));
  V3[4] = denom * 2. * (F1[4] * (F2[5] * (P3[3] - P3[0]) - P3[1] * F2[4]) +
      F1[5] * (F2[4] * (P3[0] + P3[3]) + P3[1] * F2[5]));
  V3[5] = denom * 2. * cI * (F1[4] * (F2[5] * (P3[1] + cI * (P3[2])) + P3[0] *
      F2[4]) + F1[5] * (F2[4] * (+cI * (P3[2]) - P3[1]) - P3[0] * F2[5]));
}


void VVVV7_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP10; 
  double P3[4]; 
  complex<double> TMP32; 
  complex<double> TMP28; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  complex<double> TMP37; 
  double P2[4]; 
  complex<double> TMP46; 
  complex<double> TMP33; 
  complex<double> TMP34; 
  complex<double> TMP24; 
  complex<double> TMP49; 
  complex<double> TMP12; 
  double P1[4]; 
  complex<double> TMP23; 
  complex<double> TMP30; 
  complex<double> TMP27; 
  double P4[4]; 
  complex<double> TMP11; 
  complex<double> TMP22; 
  complex<double> TMP31; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP46 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP49 = (P1[0] * P4[0] - P1[1] * P4[1] - P1[2] * P4[2] - P1[3] * P4[3]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP33 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP34 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  vertex = COUP * (TMP11 * (TMP12 * - 1. * (+cI * (TMP46 + TMP49)) + (+cI *
      (TMP27 * TMP34 + TMP31 * TMP33))) + (TMP22 * (TMP23 * (+cI * (TMP46 +
      TMP49)) + (-cI * (TMP30 * TMP34 + TMP29 * TMP37))) + (TMP10 * (-cI *
      (TMP24 * TMP33) + cI * (TMP26 * TMP37)) + (TMP12 * (+cI * (TMP28 * TMP32
      + TMP24 * TMP29)) + (TMP23 * - 1. * (+cI * (TMP25 * TMP28 + TMP26 *
      TMP31)) + TMP9 * (-cI * (TMP27 * TMP32) + cI * (TMP25 * TMP30)))))));
}


void FFV4_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * 2. * (V3[5] - V3[2]) + 2. * (F1[5]
      * (V3[3] - cI * (V3[4]))))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))) + M2 * (F1[4] * 2. * (V3[3] + cI * (V3[4])) - 2. *
      (F1[5] * (V3[2] + V3[5])))));
  F2[4] = denom * 2. * cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3]
      + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (+1./2. * (M2 * (F1[3] * (V3[3] - cI * (V3[4])) + 2. * (F1[2]
      * 1./2. * (V3[2] + V3[5])))) + F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) +
      (P2[1] * - 1. * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) +
      P2[3] * (V3[3] - cI * (V3[4])))))));
  F2[5] = denom * 2. * cI * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./2. * (M2 * (F1[3] * (V3[2] - V3[5]) + 2.
      * (F1[2] * 1./2. * (V3[3] + cI * (V3[4]))))) + F1[5] * (P2[0] * - 1. *
      (V3[2] + V3[5]) + (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI
      * (V3[3])) + P2[3] * (V3[2] + V3[5]))))));
}


void VVS2_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP37; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> TMP47; 
  complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP47 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * (-cI * (TMP23 * TMP47) + cI * (TMP30 * TMP37)); 
}


void VVVV7P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP10; 
  double P3[4]; 
  complex<double> TMP32; 
  complex<double> TMP28; 
  complex<double> TMP25; 
  complex<double> TMP11; 
  double P2[4]; 
  complex<double> TMP46; 
  complex<double> TMP24; 
  complex<double> TMP49; 
  double P1[4]; 
  complex<double> TMP30; 
  complex<double> denom; 
  complex<double> TMP27; 
  double P4[4]; 
  complex<double> TMP22; 
  complex<double> TMP31; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP46 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP49 = (P1[0] * P4[0] - P1[1] * P4[1] - P1[2] * P4[2] - P1[3] * P4[3]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP11 * (V3[2] * - 1. * (+cI * (TMP46 + TMP49)) + (+cI *
      (P4[0] * TMP27 + P3[0] * TMP31))) + (TMP22 * (V2[2] * (+cI * (TMP46 +
      TMP49)) + (-cI * (P4[0] * TMP30 + P2[0] * TMP29))) + (TMP10 * (-cI *
      (P3[0] * TMP24) + cI * (P2[0] * TMP26)) + (TMP25 * (-cI * (V2[2] * TMP28)
      + cI * (V4[2] * TMP30)) + (TMP32 * (-cI * (V4[2] * TMP27) + cI * (V3[2] *
      TMP28)) + (-cI * (V2[2] * TMP26 * TMP31) + cI * (V3[2] * TMP24 *
      TMP29)))))));
  V1[3] = denom * (TMP11 * (V3[3] * - 1. * (+cI * (TMP46 + TMP49)) + (+cI *
      (P4[1] * TMP27 + P3[1] * TMP31))) + (TMP22 * (V2[3] * (+cI * (TMP46 +
      TMP49)) + (-cI * (P4[1] * TMP30 + P2[1] * TMP29))) + (TMP10 * (-cI *
      (P3[1] * TMP24) + cI * (P2[1] * TMP26)) + (TMP25 * (-cI * (V2[3] * TMP28)
      + cI * (V4[3] * TMP30)) + (TMP32 * (-cI * (V4[3] * TMP27) + cI * (V3[3] *
      TMP28)) + (-cI * (V2[3] * TMP26 * TMP31) + cI * (V3[3] * TMP24 *
      TMP29)))))));
  V1[4] = denom * (TMP11 * (V3[4] * - 1. * (+cI * (TMP46 + TMP49)) + (+cI *
      (P4[2] * TMP27 + P3[2] * TMP31))) + (TMP22 * (V2[4] * (+cI * (TMP46 +
      TMP49)) + (-cI * (P4[2] * TMP30 + P2[2] * TMP29))) + (TMP10 * (-cI *
      (P3[2] * TMP24) + cI * (P2[2] * TMP26)) + (TMP25 * (-cI * (V2[4] * TMP28)
      + cI * (V4[4] * TMP30)) + (TMP32 * (-cI * (V4[4] * TMP27) + cI * (V3[4] *
      TMP28)) + (-cI * (V2[4] * TMP26 * TMP31) + cI * (V3[4] * TMP24 *
      TMP29)))))));
  V1[5] = denom * (TMP11 * (V3[5] * - 1. * (+cI * (TMP46 + TMP49)) + (+cI *
      (P4[3] * TMP27 + P3[3] * TMP31))) + (TMP22 * (V2[5] * (+cI * (TMP46 +
      TMP49)) + (-cI * (P4[3] * TMP30 + P2[3] * TMP29))) + (TMP10 * (-cI *
      (P3[3] * TMP24) + cI * (P2[3] * TMP26)) + (TMP25 * (-cI * (V2[5] * TMP28)
      + cI * (V4[5] * TMP30)) + (TMP32 * (-cI * (V4[5] * TMP27) + cI * (V3[5] *
      TMP28)) + (-cI * (V2[5] * TMP26 * TMP31) + cI * (V3[5] * TMP24 *
      TMP29)))))));
}


void VVVVV7P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP17; 
  complex<double> TMP20; 
  complex<double> TMP16; 
  complex<double> TMP21; 
  complex<double> denom; 
  double P5[4]; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P5[0] = V5[0].real(); 
  P5[1] = V5[1].real(); 
  P5[2] = V5[1].imag(); 
  P5[3] = V5[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP20 = (V2[2] * P5[0] - V2[3] * P5[1] - V2[4] * P5[2] - V2[5] * P5[3]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP17 = (V3[2] * P5[0] - V3[3] * P5[1] - V3[4] * P5[2] - V3[5] * P5[3]); 
  TMP16 = (V4[2] * P5[0] - V4[3] * P5[1] - V4[4] * P5[2] - V4[5] * P5[3]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P5[0] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP16 * (-cI * (V2[2] * TMP19) + cI * (V3[2] * TMP18)) + (TMP17 * (-cI *
      (TMP11 * V5[2]) + cI * (V2[2] * TMP21)) + TMP20 * (-cI * (V3[2] * TMP21)
      + cI * (V5[2] * TMP22)))));
  V1[3] = denom * (P5[1] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP16 * (-cI * (V2[3] * TMP19) + cI * (V3[3] * TMP18)) + (TMP17 * (-cI *
      (TMP11 * V5[3]) + cI * (V2[3] * TMP21)) + TMP20 * (-cI * (V3[3] * TMP21)
      + cI * (V5[3] * TMP22)))));
  V1[4] = denom * (P5[2] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP16 * (-cI * (V2[4] * TMP19) + cI * (V3[4] * TMP18)) + (TMP17 * (-cI *
      (TMP11 * V5[4]) + cI * (V2[4] * TMP21)) + TMP20 * (-cI * (V3[4] * TMP21)
      + cI * (V5[4] * TMP22)))));
  V1[5] = denom * (P5[3] * (-cI * (TMP18 * TMP22) + cI * (TMP11 * TMP19)) +
      (TMP16 * (-cI * (V2[5] * TMP19) + cI * (V3[5] * TMP18)) + (TMP17 * (-cI *
      (TMP11 * V5[5]) + cI * (V2[5] * TMP21)) + TMP20 * (-cI * (V3[5] * TMP21)
      + cI * (V5[5] * TMP22)))));
}


void VVVV4_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP23; 
  complex<double> TMP9; 
  TMP9 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  vertex = COUP * (-cI * (TMP9 * TMP10) + cI * (TMP22 * TMP23)); 
}


void FFFF6_0(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> F4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP50; 
  TMP50 = (F1[4] * F3[5] * (F2[3] * F4[2] - F2[2] * F4[3]) + F1[5] * F3[4] *
      (F2[2] * F4[3] - F2[3] * F4[2]));
  vertex = COUP * - 2. * cI * TMP50; 
}


void FFVV1_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP58; 
  complex<double> TMP57; 
  TMP58 = (F1[2] * (F2[2] * (V3[2] * (V4[2] - V4[5]) + (V3[3] * (+cI * (V4[4])
      - V4[3]) + (V3[4] * - 1. * (V4[4] + cI * (V4[3])) + V3[5] * (V4[2] -
      V4[5])))) + F2[3] * (V3[2] * - 1. * (V4[3] + cI * (V4[4])) + (V3[3] *
      (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] + V4[5])) - V3[5] * (V4[3] + cI
      * (V4[4])))))) + F1[3] * (F2[2] * (V3[2] * (+cI * (V4[4]) - V4[3]) +
      (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[2]) + cI * (V4[5])) +
      V3[5] * (V4[3] - cI * (V4[4]))))) + F2[3] * (V3[2] * (V4[2] + V4[5]) +
      (V3[3] * - 1. * (V4[3] + cI * (V4[4])) + (V3[4] * (+cI * (V4[3]) - V4[4])
      - V3[5] * (V4[2] + V4[5]))))));
  TMP57 = (F1[2] * (F2[2] * (V3[2] * (V4[2] + V4[5]) + (V3[3] * - 1. * (V4[3] +
      cI * (V4[4])) + (V3[4] * (+cI * (V4[3]) - V4[4]) - V3[5] * (V4[2] +
      V4[5])))) + F2[3] * (V3[2] * (V4[3] + cI * (V4[4])) + (V3[3] * - 1. *
      (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5])) + V3[5] *
      (V4[3] + cI * (V4[4])))))) + F1[3] * (F2[2] * (V3[2] * (V4[3] - cI *
      (V4[4])) + (V3[3] * (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI *
      (V4[2])) + V3[5] * (+cI * (V4[4]) - V4[3])))) + F2[3] * (V3[2] * (V4[2] -
      V4[5]) + (V3[3] * (+cI * (V4[4]) - V4[3]) + (V3[4] * - 1. * (V4[4] + cI *
      (V4[3])) + V3[5] * (V4[2] - V4[5]))))));
  vertex = COUP * (-cI * (TMP57) + cI * (TMP58)); 
}

void FFVV1_2_0(complex<double> F1[], complex<double> F2[], complex<double>
    V3[], complex<double> V4[], complex<double> COUP1, complex<double> COUP2,
    complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFVV1_0(F1, F2, V3, V4, COUP1, vertex); 
  FFVV2_0(F1, F2, V3, V4, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void VVVVV4P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP30; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP28; 
  complex<double> TMP27; 
  complex<double> TMP18; 
  complex<double> TMP41; 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP41 = (P1[0] * V5[2] - P1[1] * V5[3] - P1[2] * V5[4] - P1[3] * V5[5]); 
  TMP18 = (V2[2] * V5[2] - V2[3] * V5[3] - V2[4] * V5[4] - V2[5] * V5[5]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP10 * (-cI * (V4[2] * TMP41) + cI * (V5[2] * TMP28)) +
      (TMP18 * (-cI * (V3[2] * TMP28) + cI * (V4[2] * TMP27)) + (TMP21 * (-cI *
      (V2[2] * TMP27) + cI * (V3[2] * TMP30)) + TMP22 * (-cI * (V5[2] * TMP30)
      + cI * (V2[2] * TMP41)))));
  V1[3] = denom * (TMP10 * (-cI * (V4[3] * TMP41) + cI * (V5[3] * TMP28)) +
      (TMP18 * (-cI * (V3[3] * TMP28) + cI * (V4[3] * TMP27)) + (TMP21 * (-cI *
      (V2[3] * TMP27) + cI * (V3[3] * TMP30)) + TMP22 * (-cI * (V5[3] * TMP30)
      + cI * (V2[3] * TMP41)))));
  V1[4] = denom * (TMP10 * (-cI * (V4[4] * TMP41) + cI * (V5[4] * TMP28)) +
      (TMP18 * (-cI * (V3[4] * TMP28) + cI * (V4[4] * TMP27)) + (TMP21 * (-cI *
      (V2[4] * TMP27) + cI * (V3[4] * TMP30)) + TMP22 * (-cI * (V5[4] * TMP30)
      + cI * (V2[4] * TMP41)))));
  V1[5] = denom * (TMP10 * (-cI * (V4[5] * TMP41) + cI * (V5[5] * TMP28)) +
      (TMP18 * (-cI * (V3[5] * TMP28) + cI * (V4[5] * TMP27)) + (TMP21 * (-cI *
      (V2[5] * TMP27) + cI * (V3[5] * TMP30)) + TMP22 * (-cI * (V5[5] * TMP30)
      + cI * (V2[5] * TMP41)))));
}


void FFV7_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP15; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP4; 
  complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP15 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP4 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - 4. * cI * (OM3 * - 1./4. * P3[0] * (TMP4 + 4. * (TMP15)) +
      (+1./4. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] *
      F2[3]));
  V3[3] = denom * - 4. * cI * (OM3 * - 1./4. * P3[1] * (TMP4 + 4. * (TMP15)) +
      (-1./4. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] *
      F2[2]));
  V3[4] = denom * 4. * cI * (OM3 * 1./4. * P3[2] * (TMP4 + 4. * (TMP15)) +
      (+1./4. * cI * (F1[2] * F2[5]) - 1./4. * cI * (F1[3] * F2[4]) - cI *
      (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
  V3[5] = denom * 4. * cI * (OM3 * 1./4. * P3[3] * (TMP4 + 4. * (TMP15)) +
      (+1./4. * (F1[2] * F2[4]) - 1./4. * (F1[3] * F2[5]) - F1[4] * F2[2] +
      F1[5] * F2[3]));
}


void VVVVV11P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  complex<double> TMP45; 
  double P3[4]; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  complex<double> TMP19; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP45 = (P3[0] * V5[2] - P3[1] * V5[3] - P3[2] * V5[4] - P3[3] * V5[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P3[0] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP22 * (-cI * (V2[2] * TMP45) + cI * (V5[2] * TMP29)) + (TMP26 * (-cI *
      (TMP10 * V5[2]) + cI * (V2[2] * TMP19)) + V3[2] * (-cI * (TMP21 * TMP29)
      + cI * (TMP11 * TMP45)))));
  V1[3] = denom * (P3[1] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP22 * (-cI * (V2[3] * TMP45) + cI * (V5[3] * TMP29)) + (TMP26 * (-cI *
      (TMP10 * V5[3]) + cI * (V2[3] * TMP19)) + V3[3] * (-cI * (TMP21 * TMP29)
      + cI * (TMP11 * TMP45)))));
  V1[4] = denom * (P3[2] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP22 * (-cI * (V2[4] * TMP45) + cI * (V5[4] * TMP29)) + (TMP26 * (-cI *
      (TMP10 * V5[4]) + cI * (V2[4] * TMP19)) + V3[4] * (-cI * (TMP21 * TMP29)
      + cI * (TMP11 * TMP45)))));
  V1[5] = denom * (P3[3] * (-cI * (TMP11 * TMP19) + cI * (TMP10 * TMP21)) +
      (TMP22 * (-cI * (V2[5] * TMP45) + cI * (V5[5] * TMP29)) + (TMP26 * (-cI *
      (TMP10 * V5[5]) + cI * (V2[5] * TMP19)) + V3[5] * (-cI * (TMP21 * TMP29)
      + cI * (TMP11 * TMP45)))));
}


void FFFF5_1(complex<double> F2[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP38; 
  complex<double> denom; 
  double P1[4]; 
  F1[0] = +F2[0] + F3[0] + F4[0]; 
  F1[1] = +F2[1] + F3[1] + F4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  TMP38 = (F4[4] * F3[4] + F4[5] * F3[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * TMP38 * F2[2] * M1; 
  F1[3] = denom * cI * TMP38 * F2[3] * M1; 
  F1[4] = denom * cI * TMP38 * (F2[2] * (P1[3] - P1[0]) + F2[3] * (P1[1] + cI *
      (P1[2])));
  F1[5] = denom * - cI * TMP38 * (F2[2] * (+cI * (P1[2]) - P1[1]) + F2[3] *
      (P1[0] + P1[3]));
}


void FFFF3_1(complex<double> F2[], complex<double> F3[], complex<double> F4[],
    complex<double> COUP, double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + F3[0] + F4[0]; 
  F1[1] = +F2[1] + F3[1] + F4[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * 2. * cI * F3[3] * M1 * (F2[4] * F4[5] - F2[5] * F4[4]); 
  F1[3] = denom * 2. * cI * F3[2] * M1 * (F2[5] * F4[4] - F2[4] * F4[5]); 
  F1[4] = denom * - 2. * cI * (F3[2] * (F2[4] * F4[5] * (P1[1] + cI * (P1[2]))
      - F2[5] * F4[4] * (P1[1] + cI * (P1[2]))) + F3[3] * (F2[4] * F4[5] *
      (P1[0] - P1[3]) + F2[5] * F4[4] * (P1[3] - P1[0])));
  F1[5] = denom * - 2. * cI * (F3[2] * (F2[4] * - F4[5] * (P1[0] + P1[3]) +
      F2[5] * F4[4] * (P1[0] + P1[3])) + F3[3] * (F2[4] * F4[5] * (+cI *
      (P1[2]) - P1[1]) + F2[5] * F4[4] * (P1[1] - cI * (P1[2]))));
}

void FFFF3_7_1(complex<double> F2[], complex<double> F3[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M1, double W1,
    complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFFF3_1(F2, F3, F4, COUP1, M1, W1, F1); 
  FFFF7_1(F2, F3, F4, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFFF1_3(complex<double> F1[], complex<double> F2[], complex<double> F4[],
    complex<double> COUP, double M3, double W3, complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  F3[0] = +F1[0] + F2[0] + F4[0]; 
  F3[1] = +F1[1] + F2[1] + F4[1]; 
  P3[0] = -F3[0].real(); 
  P3[1] = -F3[1].real(); 
  P3[2] = -F3[1].imag(); 
  P3[3] = -F3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  F3[2] = denom * 2. * cI * F1[3] * M3 * (F2[4] * F4[5] - F2[5] * F4[4]); 
  F3[3] = denom * 2. * cI * F1[2] * M3 * (F2[5] * F4[4] - F2[4] * F4[5]); 
  F3[4] = denom * - 2. * cI * (F1[2] * (F2[4] * F4[5] * (P3[1] + cI * (P3[2]))
      - F2[5] * F4[4] * (P3[1] + cI * (P3[2]))) + F1[3] * (F2[4] * F4[5] *
      (P3[0] - P3[3]) + F2[5] * F4[4] * (P3[3] - P3[0])));
  F3[5] = denom * - 2. * cI * (F1[2] * (F2[4] * - F4[5] * (P3[0] + P3[3]) +
      F2[5] * F4[4] * (P3[0] + P3[3])) + F1[3] * (F2[4] * F4[5] * (+cI *
      (P3[2]) - P3[1]) + F2[5] * F4[4] * (P3[1] - cI * (P3[2]))));
}

void FFFF1_2_6_3(complex<double> F1[], complex<double> F2[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    double M3, double W3, complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFFF1_3(F1, F2, F4, COUP1, M3, W3, F3); 
  FFFF2_3(F1, F2, F4, COUP2, M3, W3, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F3[i] = F3[i] + Ftmp[i]; 
    i++; 
  }
  FFFF6_3(F1, F2, F4, COUP3, M3, W3, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F3[i] = F3[i] + Ftmp[i]; 
    i++; 
  }
}
void FFFF1_6_3(complex<double> F1[], complex<double> F2[], complex<double>
    F4[], complex<double> COUP1, complex<double> COUP2, double M3, double W3,
    complex<double> F3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  complex<double> COUP3; 
  int i; 
  complex<double> Ftmp[6]; 
  FFFF1_3(F1, F2, F4, COUP1, M3, W3, F3); 
  FFFF6_3(F1, F2, F4, COUP2, M3, W3, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F3[i] = F3[i] + Ftmp[i]; 
    i++; 
  }
}

void VVVV8_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> V4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP10; 
  double P3[4]; 
  complex<double> TMP32; 
  complex<double> TMP28; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  complex<double> TMP11; 
  double P2[4]; 
  complex<double> TMP36; 
  complex<double> TMP33; 
  complex<double> TMP34; 
  complex<double> TMP24; 
  complex<double> TMP12; 
  double P1[4]; 
  complex<double> TMP23; 
  complex<double> TMP30; 
  complex<double> TMP27; 
  double P4[4]; 
  complex<double> TMP37; 
  complex<double> TMP35; 
  complex<double> TMP22; 
  complex<double> TMP31; 
  complex<double> TMP26; 
  complex<double> TMP29; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP25 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP26 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP22 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP9 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP33 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP32 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP36 = (P2[0] * P4[0] - P2[1] * P4[1] - P2[2] * P4[2] - P2[3] * P4[3]); 
  TMP35 = (P3[0] * P1[0] - P3[1] * P1[1] - P3[2] * P1[2] - P3[3] * P1[3]); 
  TMP34 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  vertex = COUP * (TMP10 * (TMP9 * - 1. * (+cI * (TMP35 + TMP36)) + (+cI *
      (TMP28 * TMP33 + TMP24 * TMP34))) + (TMP22 * (TMP23 * (+cI * (TMP35 +
      TMP36)) + (-cI * (TMP30 * TMP33 + TMP32 * TMP37))) + (TMP11 * (-cI *
      (TMP31 * TMP34) + cI * (TMP25 * TMP37)) + (TMP12 * (-cI * (TMP28 * TMP29)
      + cI * (TMP26 * TMP30)) + (TMP23 * - 1. * (+cI * (TMP24 * TMP25 + TMP26 *
      TMP27)) + TMP9 * (+cI * (TMP27 * TMP29 + TMP31 * TMP32)))))));
}


void VVVS1_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> S4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP31; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> TMP33; 
  complex<double> TMP27; 
  complex<double> TMP29; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP29 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP33 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP31 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP37 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP12 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  vertex = COUP * S4[2] * (TMP10 * (-cI * (TMP37) + cI * (TMP33)) + (TMP12 *
      (-cI * (TMP29) + cI * (TMP30)) + TMP23 * (-cI * (TMP27) + cI * (TMP31))));
}


void FFFF7_0(complex<double> F1[], complex<double> F2[], complex<double> F3[],
    complex<double> F4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP44; 
  TMP44 = (F1[4] * F3[5] * (F2[2] * F4[3] - F2[3] * F4[2]) + F1[5] * F3[4] *
      (F2[3] * F4[2] - F2[2] * F4[3]));
  vertex = COUP * - 2. * cI * TMP44; 
}


void VVS1_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP23; 
  TMP23 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  vertex = COUP * - cI * TMP23 * S3[2]; 
}


void FFV7_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP14; 
  complex<double> TMP13; 
  TMP14 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP13 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - 1. * (+cI * (TMP13) + 4. * cI * (TMP14)); 
}


void VVVVV5P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> V5[], complex<double> COUP, double M1, double W1,
    complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  complex<double> TMP30; 
  complex<double> TMP21; 
  complex<double> denom; 
  complex<double> TMP28; 
  complex<double> TMP41; 
  complex<double> TMP19; 
  complex<double> TMP27; 
  V1[0] = +V2[0] + V3[0] + V4[0] + V5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + V5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP41 = (P1[0] * V5[2] - P1[1] * V5[3] - P1[2] * V5[4] - P1[3] * V5[5]); 
  TMP21 = (V4[2] * V5[2] - V4[3] * V5[3] - V4[4] * V5[4] - V4[5] * V5[5]); 
  TMP28 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP27 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP19 = (V3[2] * V5[2] - V3[3] * V5[3] - V3[4] * V5[4] - V3[5] * V5[5]); 
  TMP30 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP10 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP10 * (-cI * (V4[2] * TMP41) + cI * (V5[2] * TMP28)) +
      (TMP11 * (-cI * (V5[2] * TMP27) + cI * (V3[2] * TMP41)) + (TMP19 * (-cI *
      (V2[2] * TMP28) + cI * (V4[2] * TMP30)) + TMP21 * (-cI * (V3[2] * TMP30)
      + cI * (V2[2] * TMP27)))));
  V1[3] = denom * (TMP10 * (-cI * (V4[3] * TMP41) + cI * (V5[3] * TMP28)) +
      (TMP11 * (-cI * (V5[3] * TMP27) + cI * (V3[3] * TMP41)) + (TMP19 * (-cI *
      (V2[3] * TMP28) + cI * (V4[3] * TMP30)) + TMP21 * (-cI * (V3[3] * TMP30)
      + cI * (V2[3] * TMP27)))));
  V1[4] = denom * (TMP10 * (-cI * (V4[4] * TMP41) + cI * (V5[4] * TMP28)) +
      (TMP11 * (-cI * (V5[4] * TMP27) + cI * (V3[4] * TMP41)) + (TMP19 * (-cI *
      (V2[4] * TMP28) + cI * (V4[4] * TMP30)) + TMP21 * (-cI * (V3[4] * TMP30)
      + cI * (V2[4] * TMP27)))));
  V1[5] = denom * (TMP10 * (-cI * (V4[5] * TMP41) + cI * (V5[5] * TMP28)) +
      (TMP11 * (-cI * (V5[5] * TMP27) + cI * (V3[5] * TMP41)) + (TMP19 * (-cI *
      (V2[5] * TMP28) + cI * (V4[5] * TMP30)) + TMP21 * (-cI * (V3[5] * TMP30)
      + cI * (V2[5] * TMP27)))));
}


}  // end namespace $(namespace)s_TopEffTh

