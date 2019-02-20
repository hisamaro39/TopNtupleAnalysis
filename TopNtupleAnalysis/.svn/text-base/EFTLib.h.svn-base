#ifndef EFTLIB_H
#define EFTLIB_H

#include <iostream>
#include <string>
#include "TLorentzVector.h"
#include "TopNtupleAnalysis/Parameters_TopEffTh.h"
#include "TopNtupleAnalysis/Parameters_sm.h"

void initPDFForReweighting(const std::string &s, int id);
double pdfAlphaS(double Q2);

void initEFTModels(float eftLambda, float eftCvv);
double getEFTWeight(int i1_pid, int i2_pid, std::vector<int> f_pid, TLorentzVector i1, TLorentzVector i2, TLorentzVector t, TLorentzVector tbar, std::vector<TLorentzVector> f, double aS, double cvv = 1.0);
double getSMWeight(int i1_pid, int i2_pid, std::vector<int> f_pid, TLorentzVector i1, TLorentzVector i2, TLorentzVector t, TLorentzVector tbar, std::vector<TLorentzVector> f, double aS);

#endif

