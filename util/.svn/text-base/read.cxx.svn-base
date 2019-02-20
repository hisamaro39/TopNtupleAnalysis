/**
 * @brief Main executable to run over input mini flat ntuples,
 * read events in an object-oriented way and
 * run the event analysis class that derives from Analysis.
 * The main objective is to be minimal and to allow users to run this part
 * of the code in their laptops.
 *
 * @author Danilo Enoque Ferreira de Lima <dferreir@mail.cern.ch>
 */

#include "TopNtupleAnalysis/Event.h"
#include "TChain.h"
#include "TopNtupleAnalysis/MiniTree.h"

#include <iostream>

//#include "TopNtupleAnalysis/EventCount.h"

#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresSL.h"
#include "TopNtupleAnalysis/AnaTtresSLMtt.h"
#include "TopNtupleAnalysis/AnaTtresQCD.h"
#include "TopNtupleAnalysis/AnaTtresMM.h"

#include "TopDataPreparation/SampleXsection.h"

#include "TopNtupleAnalysis/ParseUtils.h"
#include "TopNtupleAnalysis/MMUtils.h"

#include <fstream>

#include "TROOT.h"
#include "TInterpreter.h"

#include <sstream>

#include "TopNtupleAnalysis/WeakCorrScaleFactorParam.h"

double wjetsSF(const Event &sel, const std::string &syst) {
  bool isWjets = false;
  if (sel.channelNumber() >= 363330 && sel.channelNumber() <= 363354) {
    isWjets = true;
  }
  if (sel.channelNumber() >= 363436 && sel.channelNumber() <= 363483) {
    isWjets = true;
  }
  if (!isWjets) return 1.0;
  const unsigned int runNumber = sel.runNumber_or_RandomRunNumber();

  int year = 2015;
  if (runNumber > 297730) {
    year = 2016;
  }

  int nj = 2;
  if (sel.jet().size() >= 4) {
    nj = 4;
  }

  std::vector<std::string> flavours;
  flavours.push_back("bb");
  flavours.push_back("cc");
  flavours.push_back("c");
  flavours.push_back("l");

  double hfweight = 1.0;
  double f_ca = 1.0;
  std::map<int, std::map<int, std::map<std::string, std::map<std::string, float> > > > frac;
  frac[2015] = std::map<int, std::map<std::string, std::map<std::string, float> > >();
  frac[2016] = std::map<int, std::map<std::string, std::map<std::string, float> > >();
  frac[2015][2] = std::map<std::string, std::map<std::string, float> >();
  frac[2016][2] = std::map<std::string, std::map<std::string, float> >();
  frac[2015][4] = std::map<std::string, std::map<std::string, float> >();
  frac[2016][4] = std::map<std::string, std::map<std::string, float> >();
  frac[2015][2]["e"]  = std::map<std::string, float>();
  frac[2015][2]["mu"] = std::map<std::string, float>();
  frac[2016][2]["e"]  = std::map<std::string, float>();
  frac[2016][2]["mu"] = std::map<std::string, float>();
  frac[2015][4]["e"]  = std::map<std::string, float>();
  frac[2015][4]["mu"] = std::map<std::string, float>();
  frac[2016][4]["e"]  = std::map<std::string, float>();
  frac[2016][4]["mu"] = std::map<std::string, float>();

  frac[2015][2]["e"]["bb"] = 0.091; frac[2015][2]["e"]["cc"] = 0.070; frac[2015][2]["e"]["c"] = 0.226; frac[2015][2]["e"]["l"] = 0.613;
  frac[2015][2]["mu"]["bb"] = 0.088; frac[2015][2]["mu"]["cc"] = 0.065; frac[2015][2]["mu"]["c"] = 0.216; frac[2015][2]["mu"]["l"] = 0.631;

  frac[2016][2]["e"]["bb"] = 0.090; frac[2016][2]["e"]["cc"] = 0.067; frac[2016][2]["e"]["c"] = 0.232; frac[2016][2]["e"]["l"] = 0.611;
  frac[2016][2]["mu"]["bb"] = 0.089; frac[2016][2]["mu"]["cc"] = 0.065; frac[2016][2]["mu"]["c"] = 0.221; frac[2016][2]["mu"]["l"] = 0.625;

  frac[2015][4]["e"]["bb"] = 0.140; frac[2015][4]["e"]["cc"] = 0.157; frac[2015][4]["e"]["c"] = 0.220; frac[2015][4]["e"]["l"] = 0.483;
  frac[2015][4]["mu"]["bb"] = 0.153; frac[2015][4]["mu"]["cc"] = 0.155; frac[2015][4]["mu"]["c"] = 0.210; frac[2015][4]["mu"]["l"] = 0.482;

  frac[2016][4]["e"]["bb"] = 0.148; frac[2016][4]["e"]["cc"] = 0.155; frac[2016][4]["e"]["c"] = 0.221; frac[2016][4]["e"]["l"] = 0.487;
  frac[2016][4]["mu"]["bb"] = 0.146; frac[2016][4]["mu"]["cc"] = 0.149; frac[2016][4]["mu"]["c"] = 0.213; frac[2016][4]["mu"]["l"] = 0.492;

  std::map<int, std::map<int, std::map<std::string, float> > > f_ca_map;
  f_ca_map[2015] = std::map<int, std::map<std::string, float> >();
  f_ca_map[2016] = std::map<int, std::map<std::string, float> >();
  f_ca_map[2015][2] = std::map<std::string, float>();
  f_ca_map[2016][2] = std::map<std::string, float>();
  f_ca_map[2015][4] = std::map<std::string, float>();
  f_ca_map[2016][4] = std::map<std::string, float>();

  f_ca_map[2015][2]["e"] = 1.01; f_ca_map[2015][2]["mu"] = 1.02;
  f_ca_map[2016][2]["e"] = 0.91; f_ca_map[2016][2]["mu"] = 1.11;

  f_ca_map[2015][4]["e"] = 0.76; f_ca_map[2015][4]["mu"] = 0.94;
  f_ca_map[2016][4]["e"] = 0.90; f_ca_map[2016][4]["mu"] = 0.82;

  std::map<int, std::map<int, std::map<std::string, std::map<std::string, float> > > > flav_map;
  flav_map[2015] = std::map<int, std::map<std::string, std::map<std::string, float> > >();
  flav_map[2016] = std::map<int, std::map<std::string, std::map<std::string, float> > >();
  flav_map[2015][2] = std::map<std::string, std::map<std::string, float> >();
  flav_map[2016][2] = std::map<std::string, std::map<std::string, float> >();
  flav_map[2015][4] = std::map<std::string, std::map<std::string, float> >();
  flav_map[2016][4] = std::map<std::string, std::map<std::string, float> >();
  flav_map[2015][2]["e"] = std::map<std::string, float>();
  flav_map[2015][2]["mu"] = std::map<std::string, float>();
  flav_map[2016][2]["e"] = std::map<std::string, float>();
  flav_map[2016][2]["mu"] = std::map<std::string, float>();
  flav_map[2015][4]["e"] = std::map<std::string, float>();
  flav_map[2015][4]["mu"] = std::map<std::string, float>();
  flav_map[2016][4]["e"] = std::map<std::string, float>();
  flav_map[2016][4]["mu"] = std::map<std::string, float>();

  flav_map[2015][2]["e"]["bb"] = 1.57; flav_map[2015][2]["e"]["cc"] = 1.57; flav_map[2015][2]["e"]["c"] = 1.00; flav_map[2015][2]["e"]["l"] = 0.85;
  flav_map[2016][2]["e"]["bb"] = 2.06; flav_map[2016][2]["e"]["cc"] = 2.06; flav_map[2016][2]["e"]["c"] = 1.00; flav_map[2016][2]["e"]["l"] = 0.73;

  flav_map[2015][2]["mu"]["bb"] = 1.49; flav_map[2015][2]["mu"]["cc"] = 1.49; flav_map[2015][2]["mu"]["c"] = 1.00; flav_map[2015][2]["mu"]["l"] = 0.88;
  flav_map[2016][2]["mu"]["bb"] = 1.70; flav_map[2016][2]["mu"]["cc"] = 1.70; flav_map[2016][2]["mu"]["c"] = 1.00; flav_map[2016][2]["mu"]["l"] = 0.83;

  flav_map[2015][4]["e"]["bb"] = 1.43; flav_map[2015][4]["e"]["cc"] = 1.43; flav_map[2015][4]["e"]["c"] = 0.91; flav_map[2015][4]["e"]["l"] = 0.78;
  flav_map[2016][4]["e"]["bb"] = 1.73; flav_map[2016][4]["e"]["cc"] = 1.73; flav_map[2016][4]["e"]["c"] = 0.84; flav_map[2016][4]["e"]["l"] = 0.61;

  flav_map[2015][4]["mu"]["bb"] = 1.37; flav_map[2015][4]["mu"]["cc"] = 1.37; flav_map[2015][4]["mu"]["c"] = 0.91; flav_map[2015][4]["mu"]["l"] = 0.80;
  flav_map[2016][4]["mu"]["bb"] = 1.52; flav_map[2016][4]["mu"]["cc"] = 1.52; flav_map[2016][4]["mu"]["c"] = 0.89; flav_map[2016][4]["mu"]["l"] = 0.74;

  //for (auto &c : frac) {
  //    for (auto &f : frac[c.first]) {
  for (std::map<std::string, std::map<std::string, float> >::iterator it = frac[year][nj].begin(); it != frac[year][nj].end(); ++it) {
  	for (std::map<std::string, float>::iterator jt = frac[year][nj][it->first].begin(); jt != frac[year][nj][it->first].end(); ++jt) {
	  //frac[c.first][f.first] *= flav_map[c.first][f.first];
	  frac[year][nj][it->first][jt->first] *= flav_map[year][nj][it->first][jt->first];
	}
  }

  std::map<int, std::map<std::string, std::map<std::string, float> > > flav_map_unc;
  flav_map_unc[2015] = std::map<std::string, std::map<std::string, float> >();
  flav_map_unc[2016] = std::map<std::string, std::map<std::string, float> >();
  flav_map_unc[2015]["e"]["bb"] = 0.253; flav_map_unc[2015]["e"]["cc"] = 0.253; flav_map_unc[2015]["e"]["c"] = 0.500; flav_map_unc[2015]["e"]["l"] = 0.123;
  flav_map_unc[2015]["mu"]["bb"] = 0.100; flav_map_unc[2015]["mu"]["cc"] = 0.100; flav_map_unc[2015]["mu"]["c"] = 0.500; flav_map_unc[2015]["mu"]["l"] = 0.041;
  flav_map_unc[2016]["e"]["bb"] = 0.043; flav_map_unc[2016]["e"]["cc"] = 0.043; flav_map_unc[2016]["e"]["c"] = 0.500; flav_map_unc[2016]["e"]["l"] = 0.031;
  flav_map_unc[2016]["mu"]["bb"] = 0.046; flav_map_unc[2016]["mu"]["cc"] = 0.046; flav_map_unc[2016]["mu"]["c"] = 0.500; flav_map_unc[2016]["mu"]["l"] = 0.023;

  std::string chan = "";
  if (sel.electron().size() == 1) {
    chan = "e";
  } else if (sel.muon().size() == 1) {
    chan = "mu";
  }

  std::string flav = "";
  int flag = sel.Wfilter_Sherpa_nT();
  if (flag == 3 || flag == 4) {
    flav = "bb";
  } else if (flag == 1) {
    flav = "cc";
  } else if (flag == 2) {
    flav = "c";
  } else if (flag == 5) {
    flav = "l";
  }

  f_ca = f_ca_map[year][nj][chan];
  double norm = 1.0;
  hfweight = flav_map[year][nj][chan][flav];

  std::string flavunc = "";
  if (syst.find("wbb__") != std::string::npos) {
    flavunc = "bb";
  } else if (syst.find("wcc__") != std::string::npos) {
    flavunc = "cc";
  } else if (syst.find("wc__") != std::string::npos) {
    flavunc = "c";
  } else if (syst.find("wl__") != std::string::npos) {
    flavunc = "l";
  }
  double updown = 0.0;
  if (syst.find("1up") != std::string::npos) {
    updown = 1.0;
  } else if (syst.find("1down") != std::string::npos) {
    updown = -1.0;
  }

  norm  = 0;
  //for (auto &f : flavours) {
  for (std::vector<std::string>::iterator it = flavours.begin(); it != flavours.end(); ++it) {
    std::string f = *it;
    if (flavunc == f) {
	  norm += (1.0 + updown*flav_map_unc[year][chan][f])*frac[year][nj][chan][f];
	  hfweight *= (1.0 + updown*flav_map_unc[year][chan][f]);
	} else {
	  norm += frac[year][nj][chan][f];
	}
  }
  if (syst == "wnorm__1up") {
    f_ca *= 1.10;
  } else if (syst == "wnorm__1down") {
    f_ca *= 0.90;
  }
  hfweight /= norm;
  return f_ca*hfweight;
}

double ttWeight(Event &sel) {
  double w1 = 1;
  double w2 = 1;

  TLorentzVector t    = sel.MC_t(); 
  TLorentzVector tbar = sel.MC_tbar(); 

  float pT_ttbar = (t + tbar).Perp();
  float pT_top = t.Perp();
  static Float_t ttbarPt_NOMINAL[4] = { 0.0409119,  -0.0121258,  -0.188448,  -0.316183 };
  static Float_t ttbarPtUpBins[4]= { 40.e3, 170.e3, 340.e3, 1000.e3 };
  int index = 3;
  for (unsigned int i = 0; i < 4; i++) {
    if (pT_ttbar < ttbarPtUpBins[i]) {
      index = i;
      break;
    }
  }
  w1 = 1 + ttbarPt_NOMINAL[index];
  static Float_t sequential_topPt_NOMINAL[7] = { 0.0139729, 0.0128711, 0.00951743, 0.00422323, -0.0352631, -0.0873852, -0.120025, };
  static Float_t sequential_topPtUpBins[7]= { 50.e3, 100.e3, 150.e3, 200.e3, 250.e3, 350.e3, 800.e3 };
  int index2 = 6;
  for (unsigned int i = 0; i < 7; i++) {
    if (pT_top < sequential_topPtUpBins[i]) {
      index2=i;
      break;
    }
  }
  w2 = 1 + sequential_topPt_NOMINAL[index2];
  return w1*w2;
}

int main(int argc, char **argv) {

  //ROOT::Cintex::Cintex::Enable();

  gROOT->ProcessLine("#include <vector>");
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");

  //initEventCount();

  // input files
  int help = 0;
  int isData = 0;
  int isAtlFastII = 0;
  std::string outFile = "hist_re.root,hist_rmu.root,hist_be.root,hist_bmu.root";
  std::string files = "output.root";
  std::string input_fullFileList = "";
  std::string analysis = "AnaTtresSL";
  std::string systs = "nominal";
  int doWeightSystematics = 0;
  int loose = 0;
  int _nentries = -1;
  int _btags = 1;
  int removeOverlapHighMtt = 0;
  int runMM = 0;
  int runMM_StatErr = 0;
  int applyPtRew = 0;
  std::string pdf = "";
  int applyWSF = 0;
  int applyEWK = 1;
  std::string only = "";
  int opt_dsid(-1);
  float opt_dsidSoW(0);
  std::string sumWeightsFile = "";

  static struct extendedOption extOpt[] = {
        {"help",                  no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"data",                  required_argument,     0, 'd', "Is this data?", &isData, extendedOption::eOTInt},
        {"atlFastII",             required_argument,     0, 'a', "Is this AtlFastII? (0/1)", &isAtlFastII, extendedOption::eOTInt},
        {"files",                 required_argument,     0, 'f', "Input list of comma-separated files to apply the selection on. If argument has input_ and .txt in the name, it is assumed that a file list is inside the text file provided as argument.", &files, extendedOption::eOTString},
        {"fullFiles",             required_argument,     0, 'F', "Full list of input files in this sample to be used to calculate the sum of weights for the normalisation. If left untouched, the contents of --files will be assumed to represent the full sample.", &input_fullFileList, extendedOption::eOTString},
        {"sumWeights",            required_argument,     0, 'w', "File with sum of weights.", &sumWeightsFile, extendedOption::eOTString},
        {"analysis",              required_argument,     0, 'A', "Analysis to run. Choices: AnaTtresSL", &analysis, extendedOption::eOTString},
        {"output",                required_argument,     0, 'o', "Comma-separated list of output files.", &outFile, extendedOption::eOTString},
        {"systs",                 required_argument,     0, 's', "Comma-separated list of systematics.", &systs, extendedOption::eOTString},
        {"loose",                 required_argument,     0, 'l', "Should I run over the loose TTree too?", &loose, extendedOption::eOTInt},
        {"nentries",              required_argument,     0, 'N', "Run over only the first entries if > 0.", &_nentries, extendedOption::eOTInt},
        {"btags",                 required_argument,     0, 'B', "Add cut on b-tagged jets >= abs(X). If negative use track-jet b-tagging.", &_btags, extendedOption::eOTInt},
        {"removeOverlapHighMtt",  required_argument,     0, 'R', "Veto events with true mtt > 1.1 TeV in the 410000 sample only (to be activated if one wnats to use the mtt sliced samples).", &removeOverlapHighMtt, extendedOption::eOTInt},
        {"doWeightSystematics",   required_argument,     0, 'S', "Include the variation of the systematics in the SFs.", &doWeightSystematics, extendedOption::eOTInt},
        {"runMM",                 required_argument,     0, 'M', "Implement the QCD weiths to data: runMM 1 for nominal QCD and runMM 2 for W+Jets", &runMM, extendedOption::eOTInt},
	{"runMM_StatErr",         required_argument,     0, 'm', "Shift the real/fake efficiencies to estimate the stattistical error", &runMM_StatErr, extendedOption::eOTInt},
        {"applyPtRew",            required_argument,     0, 'P', "Apply pt reweighting.", &applyPtRew, extendedOption::eOTInt},
        {"pdf",                   required_argument,     0, 'p', "Only run PDF variations.", &pdf, extendedOption::eOTString},
        {"applyWSF",              required_argument,     0, 'W', "Apply W SF.", &applyWSF, extendedOption::eOTInt},
        {"applyEWK",              required_argument,     0, 'E', "Apply electroweak correction.", &applyEWK, extendedOption::eOTInt},
        {"onlyChannel",           required_argument,     0, 'O', "Run code only for comma-separated channels.", &only, extendedOption::eOTString},
        {"dsid",                  required_argument,     0, 'D', "dsid", &opt_dsid, extendedOption::eOTInt},
        {"dsidsumOfWeight",       required_argument,     0, 'Z', "dsid 's Sum of Weight (when one needs to split jobs for 1 DSID)", &opt_dsidSoW, extendedOption::eOTFloat},

        {0, 0, 0, 0, 0, 0, extendedOption::eOTInt}
      };

  if (!parseArguments(argc, argv, extOpt) || help) {
    dumpHelp("read", extOpt, "read\nRead selected events after preselection and generate histograms.\n");
    return 0;
  } else {
    std::cout << "Dumping options:" << std::endl;
    dumpOptions(extOpt);
  }

  std::vector<int> onlyChannelList;
  std::string ocStr = only;
  for (size_t i = 0,n; i <= ocStr.length(); i=n+1) {
    n = ocStr.find_first_of(',',i);
    if (n == std::string::npos)
      n = ocStr.length();
    std::string tmp = ocStr.substr(i,n-i);
    if (tmp.length() > 0)
      onlyChannelList.push_back(std::atoi(tmp.c_str()));
  }

  // parse file list
  std::vector<std::string> fileList;
  std::vector<std::string> fullFileList;
  if ( ((files.find("input") != std::string::npos) && (files.find(".txt") != std::string::npos)) ) {
    std::cout << "Using file given as text list." << std::endl;
    ifstream f(files.c_str());
    if (!f) {
      std::cout << "Cannot open " << files << std::endl;
      std::exit(-2);
    }
    std::string thePathStr;
    while (std::getline(f, thePathStr)) {
      if (thePathStr != "") {
        size_t idx = std::string::npos;
        idx = thePathStr.find("\n");
        std::string aFile = thePathStr.substr(0, idx);
        fileList.push_back(aFile);
      }
    }
    for (std::vector<std::string>::const_iterator it = fileList.begin(); it != fileList.end(); ++it) {
      std::cout << "Input file \""<<*it<<"\""<< std::endl;
    }
  } else {
    // split by ','
    std::string argStr = files;
    for (size_t i = 0,n; i <= argStr.length(); i=n+1) {
      n = argStr.find_first_of(',',i);
      if (n == std::string::npos)
        n = argStr.length();
      std::string tmp = argStr.substr(i,n-i);
      fileList.push_back(tmp);
    }
  }
  if (fileList.size() == 0) {
    std::cout << "ERROR: You must input at least one file." << std::endl;
    std::exit(-1);
  }
  if (input_fullFileList == "") {
    fullFileList = fileList;
  } else {
    if ( ((input_fullFileList.find("input") != std::string::npos) && (input_fullFileList.find(".txt") != std::string::npos)) ) {
      std::cout << "Using file given as text list." << std::endl;
      ifstream f(input_fullFileList.c_str());
      if (!f) {
        std::cout << "Cannot open " << input_fullFileList << std::endl;
        std::exit(-2);
      }
      std::string thePathStr;
      while (std::getline(f, thePathStr)) {
        if (thePathStr != "") {
          size_t idx = std::string::npos;
          idx = thePathStr.find("\n");
          std::string aFile = thePathStr.substr(0, idx);
          fullFileList.push_back(aFile);
        }
      }
      for (std::vector<std::string>::const_iterator it = fullFileList.begin(); it != fullFileList.end(); ++it) {
        std::cout << "(full input for norm.) Input file \""<<*it<<"\""<< std::endl;
      }
    } else {
      // split by ','
      std::string argStr = input_fullFileList;
      for (size_t i = 0,n; i <= argStr.length(); i=n+1) {
        n = argStr.find_first_of(',',i);
        if (n == std::string::npos)
          n = argStr.length();
        std::string tmp = argStr.substr(i,n-i);
        fullFileList.push_back(tmp);
      }
    }
  }

  // get output files
  std::string outStr = outFile;
  std::vector<std::string> outList;
  for (size_t i = 0,n; i <= outStr.length(); i=n+1) {
    n = outStr.find_first_of(',',i);
    if (n == std::string::npos)
      n = outStr.length();
    std::string tmp = outStr.substr(i,n-i);
    outList.push_back(tmp);
  }

  // split systs by comma
  std::vector<std::string> systsList;
  std::vector<std::string> systsListWithBlankNominal;
  if (systs.find(".txt") != std::string::npos) {
    // systs is a file with the list of systematics
    ifstream f(systs.c_str());
    if (!f) {
      std::cout << "Cannot open " << systs << std::endl;
      std::exit(-2);
    }
    std::string finalSysts = "";
    std::string theLineStr;
    while (std::getline(f, theLineStr)) {
      if (theLineStr != "") {
        size_t idx = std::string::npos;
        idx = theLineStr.find("\n");
        std::string aLine = theLineStr.substr(0, idx);
        if (finalSysts.size() != 0 && finalSysts[finalSysts.size()-1] != ',') finalSysts += ",";
        finalSysts += aLine;
      }
    }
    systs = finalSysts;
  }
  for (size_t i = 0,n; i <= systs.length(); i=n+1) {
    n = systs.find_first_of(',',i);
    if (n == std::string::npos)
      n = systs.length();
    std::string tmp = systs.substr(i,n-i);
    std::string tmpWithBlankNominal = tmp;
    if (tmpWithBlankNominal == "nominal") tmpWithBlankNominal = "";
    systsList.push_back(tmp);
    systsListWithBlankNominal.push_back(tmpWithBlankNominal);
    if (loose) {
      systsListWithBlankNominal.push_back(tmpWithBlankNominal+std::string("_Loose"));
    }
  }

  MMUtils * MM_0b_re      = NULL;
  MMUtils * MM_0b_rmu     = NULL;
  MMUtils * MM_1b_re      = NULL;
  MMUtils * MM_1b_rmu     = NULL;
    
  if (runMM){ // >= 2 jets
  	// MMUtils(isBtagged, real2015, fake2015, real2016, fake2016)
	MM_0b_re  = new MMUtils(0, "scripts/QCDestimation/RATES_2015/resolved_e_eff_ttbar.root", // real2015
				   "scripts/QCDestimation/RATES_2015/resolved_e_btag0_fake.root",	// fake2015
				
		         	   "scripts/QCDestimation/RATES_2016/resolved_e_eff_ttbar.root",	// real2016
				   "scripts/QCDestimation/RATES_2016/resolved_e_btag0_fake.root");	// fake2016
				
	MM_1b_re  = new MMUtils(1, "scripts/QCDestimation/RATES_2015/resolved_e_eff_ttbar.root", // real2015  
				   "scripts/QCDestimation/RATES_2015/resolved_e_btag1_fake.root",	// fake2015
				
				   "scripts/QCDestimation/RATES_2016/resolved_e_eff_ttbar.root",   // real2016 
				   "scripts/QCDestimation/RATES_2016/resolved_e_btag1_fake.root"); // fake2016

	MM_0b_rmu = new MMUtils(0, "scripts/QCDestimation/RATES_2015/resolved_mu_eff_ttbar.root",// real2015  
				   "scripts/QCDestimation/RATES_2015/resolved_mu_btag0_fake.root",	// fake2015
				
				  "scripts/QCDestimation/RATES_2016/resolved_mu_eff_ttbar.root",	// real2016  
				  "scripts/QCDestimation/RATES_2016/resolved_mu_btag0_fake.root");// fake2016
				
	MM_1b_rmu = new MMUtils(1, "scripts/QCDestimation/RATES_2015/resolved_mu_eff_ttbar.root",// real2015  
				   "scripts/QCDestimation/RATES_2015/resolved_mu_btag1_fake.root",// fake2015
				
				   "scripts/QCDestimation/RATES_2016/resolved_mu_eff_ttbar.root",	// real2016
				   "scripts/QCDestimation/RATES_2016/resolved_mu_btag1_fake.root");// fake2016
	  
  } 

  std::vector<std::string> pdfList;
  for (size_t i = 0,n; i <= pdf.length(); i=n+1) {
    n = pdf.find_first_of(',',i);
    if (n == std::string::npos)
      n = pdf.length();
    std::string tmp = pdf.substr(i,n-i);
    pdfList.push_back(tmp);
  }


  // retrieve, list of sum of weights
  std::map<int, float> sumOfWeights;
  std::map<int, std::map<std::string, std::vector<float> > > PDFsumOfWeights;
  TChain t_sumWeights("sumWeights");
  TChain t_PDFsumWeights("PDFsumWeights");
  if(!isData && sumWeightsFile == "") {
    for (int k = 0; k < fullFileList.size(); ++k) {
      t_sumWeights.Add(fullFileList[k].c_str());
      if (pdf != "") t_PDFsumWeights.Add(fullFileList[k].c_str());
    }
    int dsid;
    float value;
    int dsidPdf;
    std::map<std::string, std::vector<float> *> valuePdf;
    t_sumWeights.SetBranchAddress("dsid", &dsid);
    t_sumWeights.SetBranchAddress("totalEventsWeighted", &value);
    t_PDFsumWeights.SetBranchAddress("dsid", &dsidPdf);
    for (int k = 0; k < pdfList.size(); ++k) {
      valuePdf[pdfList[k].c_str()] = 0;
    }
    for (int k = 0; k < pdfList.size(); ++k) {
      t_PDFsumWeights.SetBranchAddress(pdfList[k].c_str(), &valuePdf[pdfList[k]]);
    }
    for (int k = 0; k < t_sumWeights.GetEntries(); ++k) {
      t_sumWeights.GetEntry(k);
      if (sumOfWeights.find(dsid) == sumOfWeights.end()) sumOfWeights[dsid] = 0;
        sumOfWeights[dsid] += value;
    }
    if(opt_dsid>0 && opt_dsid==dsid){
       cout << "For DSID= " << opt_dsid << " : Overwriting sumOfWeights=" << sumOfWeights[dsid] << " by " << opt_dsidSoW << endl;
       sumOfWeights[dsid]=opt_dsidSoW;
    }
    cout << "sumOfWeights[ "<< dsid <<" ] = " << sumOfWeights[dsid] << endl;
    
    for (int k = 0; k < t_PDFsumWeights.GetEntries(); ++k) {
      t_PDFsumWeights.GetEntry(k);
      if (PDFsumOfWeights.find(dsid) == PDFsumOfWeights.end()) PDFsumOfWeights[dsid] = std::map<std::string, std::vector<float> >();
      for (int l = 0; l < pdfList.size(); ++l) {
        if (PDFsumOfWeights[dsid].find(pdfList[l]) == PDFsumOfWeights[dsid].end()) PDFsumOfWeights[dsid][pdfList[l]] = std::vector<float>();
        if (PDFsumOfWeights[dsid][pdfList[l]].size() == 0)
          PDFsumOfWeights[dsid][pdfList[l]].resize(valuePdf[pdfList[l]]->size());

        for (int m = 0; m < PDFsumOfWeights[dsid][pdfList[l]].size(); ++m)
          PDFsumOfWeights[dsid][pdfList[l]][m] += valuePdf[pdfList[l]]->at(m);
      }
    }
  }
  if (!isData && sumWeightsFile != "") {
    ifstream f(sumWeightsFile.c_str());
    if (!f) {
      std::cout << "Cannot open " << sumWeightsFile << std::endl;
      std::exit(-2);
    }
    std::string thePathStr;
    while (std::getline(f, thePathStr)) {
      if (thePathStr != "") {
        size_t idx = std::string::npos;
        idx = thePathStr.find("\n");
        std::string aFile = thePathStr.substr(0, idx);
		std::stringstream ss(aFile);
		int dsid;
		float sw;
		ss >> dsid >> sw;
		sumOfWeights[dsid] = sw;
		//std::cout << "sumOfWeights " << dsid << " = " << sumOfWeights[dsid]  << std::endl;
      }
    }
  }

  
  int n_eigenvars_b = 0;
  int n_eigenvars_c = 0;
  int n_eigenvars_l = 0;
  //std::vector<std::string> trackjetSFs;
  //std::string trackjet_pre = "trackjet_btagSF_70_eigenvars";
  //size_t trackjet_presize = std::string("trackjet_btagSF_70_eigenvars").size();

  if (doWeightSystematics && pdf == "") { // include SF systs. too
    std::cout << "adding more systematics" << std::endl;
    systsListWithBlankNominal.push_back("eTrigSF__1up");
    systsListWithBlankNominal.push_back("eTrigSF__1down");
    systsListWithBlankNominal.push_back("eRecoSF__1up");
    systsListWithBlankNominal.push_back("eRecoSF__1down");
    systsListWithBlankNominal.push_back("eIDSF__1up");
    systsListWithBlankNominal.push_back("eIDSF__1down");
    systsListWithBlankNominal.push_back("eIsolSF__1up");
    systsListWithBlankNominal.push_back("eIsolSF__1down");

    systsListWithBlankNominal.push_back("muTrigStatSF__1up");
    systsListWithBlankNominal.push_back("muTrigStatSF__1down");
    systsListWithBlankNominal.push_back("muIDStatSF__1up");
    systsListWithBlankNominal.push_back("muIDStatSF__1down");
    systsListWithBlankNominal.push_back("muIsolStatSF__1up");
    systsListWithBlankNominal.push_back("muIsolStatSF__1down");

    systsListWithBlankNominal.push_back("muTrigSystSF__1up");
    systsListWithBlankNominal.push_back("muTrigSystSF__1down");
    systsListWithBlankNominal.push_back("muIDSystSF__1up");
    systsListWithBlankNominal.push_back("muIDSystSF__1down");
    systsListWithBlankNominal.push_back("muIsolSystSF__1up");
    systsListWithBlankNominal.push_back("muIsolSystSF__1down");

    systsListWithBlankNominal.push_back("wnorm__1up");
    systsListWithBlankNominal.push_back("wnorm__1down");

    systsListWithBlankNominal.push_back("wbb__1up");
    systsListWithBlankNominal.push_back("wbb__1down");

    systsListWithBlankNominal.push_back("wcc__1up");
    systsListWithBlankNominal.push_back("wcc__1down");

    systsListWithBlankNominal.push_back("wc__1up");
    systsListWithBlankNominal.push_back("wc__1down");

    systsListWithBlankNominal.push_back("wl__1up");
    systsListWithBlankNominal.push_back("wl__1down");

    systsListWithBlankNominal.push_back("ttEWK__1up");
    systsListWithBlankNominal.push_back("ttEWK__1down");

    // count b-tagging variations
    MiniTree mt(false, fileList[0].c_str(), "nominal");
    Event sel;
    mt.read(0, sel);

    if (_btags > 0) {
      n_eigenvars_b = mt.vf("weight_bTagSF_70_eigenvars_B_up")->size();
      n_eigenvars_c = mt.vf("weight_bTagSF_70_eigenvars_C_up")->size();
      n_eigenvars_l = mt.vf("weight_bTagSF_70_eigenvars_Light_up")->size();

      for (int i = 0; i < n_eigenvars_b; ++i) {
        std::stringstream ss;
	ss << i;
        systsListWithBlankNominal.push_back("btagbSF_"+ss.str()+"__1up");
        systsListWithBlankNominal.push_back("btagbSF_"+ss.str()+"__1down");
      }
      for (int i = 0; i < n_eigenvars_c; ++i) {
        std::stringstream ss;
	ss << i;
        systsListWithBlankNominal.push_back("btagcSF_"+ss.str()+"__1up");
        systsListWithBlankNominal.push_back("btagcSF_"+ss.str()+"__1down");
      }
      for (int i = 0; i < n_eigenvars_l; ++i) {
        std::stringstream ss;
	ss << i;
        systsListWithBlankNominal.push_back("btaglSF_"+ss.str()+"__1up");
        systsListWithBlankNominal.push_back("btaglSF_"+ss.str()+"__1down");
      }
      systsListWithBlankNominal.push_back("btageSF_0__1up");
      systsListWithBlankNominal.push_back("btageSF_0__1down");
      systsListWithBlankNominal.push_back("btageSF_1__1up");
      systsListWithBlankNominal.push_back("btageSF_1__1down");
      if (loose) {
        for (int i = 0; i < n_eigenvars_b; ++i) {
          std::stringstream ss;
  	  ss << i;
          systsListWithBlankNominal.push_back("btagbSF_"+ss.str()+"__1up_Loose");
          systsListWithBlankNominal.push_back("btagbSF_"+ss.str()+"__1down_Loose");
        }
        for (int i = 0; i < n_eigenvars_c; ++i) {
          std::stringstream ss;
  	  ss << i;
          systsListWithBlankNominal.push_back("btagcSF_"+ss.str()+"__1up_Loose");
          systsListWithBlankNominal.push_back("btagcSF_"+ss.str()+"__1down_Loose");
        }
        for (int i = 0; i < n_eigenvars_l; ++i) {
          std::stringstream ss;
  	  ss << i;
          systsListWithBlankNominal.push_back("btaglSF_"+ss.str()+"__1up_Loose");
          systsListWithBlankNominal.push_back("btaglSF_"+ss.str()+"__1down_Loose");
        }
        systsListWithBlankNominal.push_back("btageSF_0__1up_Loose");
        systsListWithBlankNominal.push_back("btageSF_0__1down_Loose");
        systsListWithBlankNominal.push_back("btageSF_1__1up_Loose");
        systsListWithBlankNominal.push_back("btageSF_1__1down_Loose");
      }
    }
    if (_btags <= 0) {
      /*
      for (std::map<std::string, MiniTree::MTType>::const_iterator it = mt.m_brs.begin(); it != mt.m_brs.end(); ++it) {
        if (it->first.find(trackjet_pre) != std::string::npos) {
          systsListWithBlankNominal.push_back(std::string("bsf_")+it->first.substr(trackjet_presize));
          if (loose) systsListWithBlankNominal.push_back(std::string("bsf_")+it->first.substr(trackjet_presize)+std::string("_Loose"));
          trackjetSFs.push_back(std::string("bsf_")+it->first.substr(trackjet_presize));
          if (loose) trackjetSFs.push_back(std::string("bsf_")+it->first.substr(trackjet_presize)+std::string("_Loose"));
        }
      }
      */

      n_eigenvars_b = mt.vf("weight_trackjet_bTagSF_70_eigenvars_B_up")->size();
      n_eigenvars_c = mt.vf("weight_trackjet_bTagSF_70_eigenvars_C_up")->size();
      n_eigenvars_l = mt.vf("weight_trackjet_bTagSF_70_eigenvars_Light_up")->size();

      for (int i = 0; i < n_eigenvars_b; ++i) {
        std::stringstream ss;
  	ss << i;
        systsListWithBlankNominal.push_back("btagbSF_"+ss.str()+"__1up");
        systsListWithBlankNominal.push_back("btagbSF_"+ss.str()+"__1down");
      }
      for (int i = 0; i < n_eigenvars_c; ++i) {
        std::stringstream ss;
  	ss << i;
        systsListWithBlankNominal.push_back("btagcSF_"+ss.str()+"__1up");
        systsListWithBlankNominal.push_back("btagcSF_"+ss.str()+"__1down");
      }
      for (int i = 0; i < n_eigenvars_l; ++i) {
        std::stringstream ss;
  	ss << i;
        systsListWithBlankNominal.push_back("btaglSF_"+ss.str()+"__1up");
        systsListWithBlankNominal.push_back("btaglSF_"+ss.str()+"__1down");
      }

      systsListWithBlankNominal.push_back("btagbSF_0_pt1__1up");
      systsListWithBlankNominal.push_back("btagbSF_0_pt1__1down");
      systsListWithBlankNominal.push_back("btagbSF_0_pt2__1up");
      systsListWithBlankNominal.push_back("btagbSF_0_pt2__1down");
      systsListWithBlankNominal.push_back("btagbSF_0_pt3__1up");
      systsListWithBlankNominal.push_back("btagbSF_0_pt3__1down");

      systsListWithBlankNominal.push_back("btagcSF_0_pt1__1up");
      systsListWithBlankNominal.push_back("btagcSF_0_pt1__1down");
      systsListWithBlankNominal.push_back("btagcSF_0_pt2__1up");
      systsListWithBlankNominal.push_back("btagcSF_0_pt2__1down");
      systsListWithBlankNominal.push_back("btagcSF_0_pt3__1up");
      systsListWithBlankNominal.push_back("btagcSF_0_pt3__1down");

      systsListWithBlankNominal.push_back("btaglSF_0_pt1__1up");
      systsListWithBlankNominal.push_back("btaglSF_0_pt1__1down");
      systsListWithBlankNominal.push_back("btaglSF_0_pt2__1up");
      systsListWithBlankNominal.push_back("btaglSF_0_pt2__1down");
      systsListWithBlankNominal.push_back("btaglSF_0_pt3__1up");
      systsListWithBlankNominal.push_back("btaglSF_0_pt3__1down");



      systsListWithBlankNominal.push_back("btageSF_0__1up");
      systsListWithBlankNominal.push_back("btageSF_0__1down");
      systsListWithBlankNominal.push_back("btageSF_1__1up");
      systsListWithBlankNominal.push_back("btageSF_1__1down");
      if (loose) {
        for (int i = 0; i < n_eigenvars_b; ++i) {
        std::stringstream ss;
  	ss << i;
          systsListWithBlankNominal.push_back("btagbSF_"+ss.str()+"__1up_Loose");
          systsListWithBlankNominal.push_back("btagbSF_"+ss.str()+"__1down_Loose");
        }
        for (int i = 0; i < n_eigenvars_c; ++i) {
        std::stringstream ss;
  	ss << i;
          systsListWithBlankNominal.push_back("btagcSF_"+ss.str()+"__1up_Loose");
          systsListWithBlankNominal.push_back("btagcSF_"+ss.str()+"__1down_Loose");
        }
        for (int i = 0; i < n_eigenvars_l; ++i) {
        std::stringstream ss;
  	ss << i;
          systsListWithBlankNominal.push_back("btaglSF_"+ss.str()+"__1up_Loose");
          systsListWithBlankNominal.push_back("btaglSF_"+ss.str()+"__1down_Loose");
        }
        systsListWithBlankNominal.push_back("btageSF_0__1up_Loose");
        systsListWithBlankNominal.push_back("btageSF_0__1down_Loose");
        systsListWithBlankNominal.push_back("btageSF_1__1up_Loose");
        systsListWithBlankNominal.push_back("btageSF_1__1down_Loose");
      }
    }
    if (loose) {
      systsListWithBlankNominal.push_back("eTrigSF__1up_Loose");
      systsListWithBlankNominal.push_back("eTrigSF__1down_Loose");
      systsListWithBlankNominal.push_back("eRecoSF__1up_Loose");
      systsListWithBlankNominal.push_back("eRecoSF__1down_Loose");
      systsListWithBlankNominal.push_back("eIDSF__1up_Loose");
      systsListWithBlankNominal.push_back("eIDSF__1down_Loose");
      systsListWithBlankNominal.push_back("eIsolSF__1up_Loose");
      systsListWithBlankNominal.push_back("eIsolSF__1down_Loose");

      systsListWithBlankNominal.push_back("muTrigStatSF__1up_Loose");
      systsListWithBlankNominal.push_back("muTrigStatSF__1down_Loose");
      systsListWithBlankNominal.push_back("muIDStatSF__1up_Loose");
      systsListWithBlankNominal.push_back("muIDStatSF__1down_Loose");
      systsListWithBlankNominal.push_back("muIsolStatSF__1up_Loose");
      systsListWithBlankNominal.push_back("muIsolStatSF__1down_Loose");

      systsListWithBlankNominal.push_back("muTrigSystSF__1up_Loose");
      systsListWithBlankNominal.push_back("muTrigSystSF__1down_Loose");
      systsListWithBlankNominal.push_back("muIDSystSF__1up_Loose");
      systsListWithBlankNominal.push_back("muIDSystSF__1down_Loose");
      systsListWithBlankNominal.push_back("muIsolSystSF__1up_Loose");
      systsListWithBlankNominal.push_back("muIsolSystSF__1down_Loose");

      systsListWithBlankNominal.push_back("boostedWSF__1up_Loose");
      systsListWithBlankNominal.push_back("boostedWSF__1down_Loose");

      systsListWithBlankNominal.push_back("ttEWK__1up_Loose");
      systsListWithBlankNominal.push_back("ttEWK__1down_Loose");
    }
  } else if (pdf != "") {
    for (int m = 0; m < pdfList.size(); ++m) {
      int nvar = (*PDFsumOfWeights.begin()).second[pdfList[m]].size();
      for (int l = 0; l < nvar; ++l) {
        std::string s = "pdf_";
        s += pdfList[m];
        s += "_";
	std::stringstream ss;
	ss << l;
        s += ss.str();
        systsListWithBlankNominal.push_back(s);
      }
    }
  }

  std::vector<Analysis *> vec_analysis;
  
  if (analysis == "AnaTtresSL") {
    vec_analysis.push_back(new AnaTtresSL(outList[0], true,  false, systsListWithBlankNominal)); // resolved electron
    vec_analysis.push_back(new AnaTtresSL(outList[1], false, false, systsListWithBlankNominal)); // resolved muon    
    vec_analysis.push_back(new AnaTtresSL(outList[2], true,  true,  systsListWithBlankNominal)); // boosted  electron
    vec_analysis.push_back(new AnaTtresSL(outList[3], false, true,  systsListWithBlankNominal)); // boosted  muon
  } else if(analysis == "AnaTtresSL_QCDVR1_2j_2015" || analysis == "AnaTtresSL_QCDVR2_2j_2015" ||analysis == "AnaTtresSL_QCDCR2j_2015" || analysis == "AnaTtresSL_QCDCR4j_2015" || analysis == "AnaTtresSL_WjetsCR2j_2015" || analysis == "AnaTtresSL_SR4j_2015" || analysis == "AnaTtresSL_QCDVR1_4j_2015" || analysis == "AnaTtresSL_QCDVR2_4j_2015" 
         || analysis == "AnaTtresSL_QCDVR1_2j_2016" || analysis == "AnaTtresSL_QCDVR2_2j_2016" ||analysis == "AnaTtresSL_QCDCR2j_2016" || analysis == "AnaTtresSL_QCDSR2j_2016" || analysis == "AnaTtresSL_QCDCR4j_2016" || analysis == "AnaTtresSL_QCDVR1_4j_2016" || analysis == "AnaTtresSL_QCDVR2_4j_2016"){
    vec_analysis.push_back(new AnaTtresMM(outList[0], true,  false, systsListWithBlankNominal, opt_dsid)); // resolved electron
    vec_analysis.push_back(new AnaTtresMM(outList[1], false, false, systsListWithBlankNominal, opt_dsid)); // resolved muon    
    vec_analysis.push_back(new AnaTtresMM(outList[2], true,  true,  systsListWithBlankNominal, opt_dsid)); // boosted  electron
    vec_analysis.push_back(new AnaTtresMM(outList[3], false, true,  systsListWithBlankNominal, opt_dsid)); // boosted  muon 
  } 
    else if(analysis == "AnaTtresWQCDreal_2015"|| analysis == "AnaTtresQCDreal_2015"|| analysis == "AnaTtresQCDfake_2015"||analysis == "AnaTtresWQCDfake_2015"
         || analysis == "AnaTtresWQCDreal_2016"|| analysis == "AnaTtresQCDreal_2016"|| analysis == "AnaTtresQCDfake_2016"||analysis == "AnaTtresWQCDfake_2016") {
    vec_analysis.push_back(new AnaTtresQCD(outList[0], true,  false, systsListWithBlankNominal, opt_dsid) ); //resolved electron
    vec_analysis.push_back(new AnaTtresQCD(outList[1], false, false, systsListWithBlankNominal, opt_dsid) ); // resolved muon
    vec_analysis.push_back(new AnaTtresQCD(outList[2], true,  true,  systsListWithBlankNominal, opt_dsid) ); // boosted  electron
    vec_analysis.push_back(new AnaTtresQCD(outList[3], false, true,  systsListWithBlankNominal, opt_dsid) ); // boosted  muon
  }

    else if (analysis == "AnaTtresSLMtt") {
    vec_analysis.push_back(new AnaTtresSLMtt(outList[0], true,  false, systsListWithBlankNominal)); // resolved electron
    vec_analysis.push_back(new AnaTtresSLMtt(outList[1], false, false, systsListWithBlankNominal)); // resolved muon    
    vec_analysis.push_back(new AnaTtresSLMtt(outList[2], true,  true,  systsListWithBlankNominal)); // boosted  electron
    vec_analysis.push_back(new AnaTtresSLMtt(outList[3], false, true,  systsListWithBlankNominal)); // boosted  muon
  }

  SampleXsection sampleXsection;
  sampleXsection.readFromFile("../TopDataPreparation/data/XSection-MC15-13TeV.data");

  WeakCorr::WeakCorrScaleFactorParam ewkTool("share/EWcorr_param.root");

  unsigned long procEvents = 0;

  // systsList contains the list of TTrees representing systematics
  // given by the user
  for (size_t systIdx = 0; systIdx < systsList.size(); ++systIdx) {
    std::cout<< "* Running on syst: " << systIdx << std::endl;
    // it does not have the xxx_Loose TTrees, so we add it ourselves and run over it too
    std::vector<std::string> lepton_modes;
    
    if (!runMM && loose){
        lepton_modes.push_back(systsList[systIdx]+std::string("_Loose"));
    } else if (!runMM && !loose){
        lepton_modes.push_back(systsList[systIdx]);
    } else {
        //lepton_modes.push_back(systsList[systIdx]);
        lepton_modes.push_back(systsList[systIdx]+std::string("_Loose"));
    }
    
    for (size_t lmodeIdx = 0; lmodeIdx < lepton_modes.size(); ++lmodeIdx) {
      std::cout<< "** Running on lepton mode: " << lmodeIdx << std::endl;

      // now we are looping over all possible TTrees with all systematic uncertainties
      // Call MiniTree to open the files and read that TTree
      std::string tname = lepton_modes[lmodeIdx];
      MiniTree mt(false, fileList[0], tname.c_str());
      for (int k = 1; k < fileList.size(); ++k) {
        mt.addFileToRead(fileList[k], tname.c_str());
      }

      // for the nominal "systematic unc.", use an empty string as suffix to be backward compatible
      // also, it would be horrible to plot histograms with the nominal suffix all the time ...
      std::string systSuffixForHistograms = tname;
      if (systSuffixForHistograms == "nominal")
        systSuffixForHistograms = "";
      else if (systSuffixForHistograms == "nominal_Loose")
        systSuffixForHistograms = "_Loose";

      Event sel; // selected objects

      // now loop over all entries for this systematic uncertainty
      int nentries = mt.GetEntries();
      if (_nentries < nentries && _nentries > 0)
        nentries = _nentries;
      for (int k = 0; k < nentries; ++k) {
        if (k % 1000 == 0)
          std::cout << "("<< tname << ") Entry " << k << "/" << nentries << std::endl;

        bool isTight = false;
    
        // read the event into an Event object to be sent to the Analysis code later
        mt.read(k, sel);
        int channel = sel.channelNumber();

        if (sel.electron().size() == 1 && sel.muon().size() == 0)
            isTight = sel.electron()[0].isTightPP();
        else if (sel.electron().size() == 0 && sel.muon().size() == 1)
            isTight = sel.muon()[0].isTight();

        if (!runMM && loose) {
            // for the fake rate or efficiency estimate, just separate loose and tight
            if (isTight)
                systSuffixForHistograms = "";
            else
                systSuffixForHistograms = "_Loose";
        } else if (runMM && loose) { // run the QCD matrix method in the SR adding weights
            systSuffixForHistograms = "";
        }

        if (channel == 410000 && removeOverlapHighMtt) // for SM ttbar, we have mtt sliced samples above 1.1 TeV
          if (sel.MC_ttbar_beforeFSR().M() > 1.1e6)
            continue;

    
        // this is list just contains systSuffixForHistograms for all systematics
        // except for the nominal and nominal_Loose, on which it contains
        // the electron SF, muon SF and b-tagging SF systematics
        // these systematics do not show up as separate TTrees, so they need special treatment
        std::vector<std::string> weightSystematics;
        if ( !(doWeightSystematics || pdf != "") || (systSuffixForHistograms != "" && systSuffixForHistograms != "_Loose") ) {
          weightSystematics.push_back(systSuffixForHistograms);
        } else if (pdf != "") {
          weightSystematics.push_back(systSuffixForHistograms);
          if (pdfList.size() != 0) {
            for (int m = 0; m < pdfList.size(); ++m) {
              int nvar = (*PDFsumOfWeights.begin()).second[pdfList[m]].size();
              for (int l = 0; l < nvar; ++l) {
                std::string s = "pdf_";
                s += pdfList[m];
                s += "_";
		std::stringstream ss;
		ss << l;
                s += ss.str();
                weightSystematics.push_back(s);
              }
            }
          }
        } else { // apply variations on the nominal
          weightSystematics.push_back(systSuffixForHistograms);
          weightSystematics.push_back(std::string("eTrigSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("eTrigSF__1down")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("eRecoSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("eRecoSF__1down")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("eIDSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("eIDSF__1down")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("eIsolSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("eIsolSF__1down")+systSuffixForHistograms);

          weightSystematics.push_back(std::string("muTrigStatSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muTrigStatSF__1down")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muIDStatSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muIDStatSF__1down")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muIsolStatSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muIsolStatSF__1down")+systSuffixForHistograms);

          weightSystematics.push_back(std::string("muTrigSystSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muTrigSystSF__1down")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muIDSystSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muIDSystSF__1down")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muIsolSystSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("muIsolSystSF__1down")+systSuffixForHistograms);

          weightSystematics.push_back(std::string("boostedWSF__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("boostedWSF__1down")+systSuffixForHistograms);

          weightSystematics.push_back(std::string("ttEWK__1up")+systSuffixForHistograms);
          weightSystematics.push_back(std::string("ttEWK__1down")+systSuffixForHistograms);

          //if (_btags > 0) {
            for (int i = 0; i < n_eigenvars_b; ++i) {
              std::stringstream ss;
     	      ss << i;
              weightSystematics.push_back(std::string("btagbSF_"+ss.str()+"__1up")+systSuffixForHistograms);
              weightSystematics.push_back(std::string("btagbSF_"+ss.str()+"__1down")+systSuffixForHistograms);
            }
            weightSystematics.push_back(std::string("btagbSF_0_pt1__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagbSF_0_pt1__1down")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagbSF_0_pt2__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagbSF_0_pt2__1down")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagbSF_0_pt3__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagbSF_0_pt3__1down")+systSuffixForHistograms);
            for (int i = 0; i < n_eigenvars_c; ++i) {
              std::stringstream ss;
     	      ss << i;
              weightSystematics.push_back(std::string("btagcSF_"+ss.str()+"__1up")+systSuffixForHistograms);
              weightSystematics.push_back(std::string("btagcSF_"+ss.str()+"__1down")+systSuffixForHistograms);
            }
            weightSystematics.push_back(std::string("btagcSF_0_pt1__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagcSF_0_pt1__1down")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagcSF_0_pt2__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagcSF_0_pt2__1down")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagcSF_0_pt3__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btagcSF_0_pt3__1down")+systSuffixForHistograms);
            for (int i = 0; i < n_eigenvars_l; ++i) {
              std::stringstream ss;
     	      ss << i;
              weightSystematics.push_back(std::string("btaglSF_"+ss.str()+"__1up")+systSuffixForHistograms);
              weightSystematics.push_back(std::string("btaglSF_"+ss.str()+"__1down")+systSuffixForHistograms);
            }
            weightSystematics.push_back(std::string("btaglSF_0_pt1__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btaglSF_0_pt1__1down")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btaglSF_0_pt2__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btaglSF_0_pt2__1down")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btaglSF_0_pt3__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btaglSF_0_pt3__1down")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btageSF_0__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btageSF_0__1down")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btageSF_1__1up")+systSuffixForHistograms);
            weightSystematics.push_back(std::string("btageSF_1__1down")+systSuffixForHistograms);
          //} else if (_btags < 0) {
          //  for (size_t i = 0; i < trackjetSFs.size(); ++i) {
          //    weightSystematics.push_back(trackjetSFs[i]+systSuffixForHistograms);
          //  }
          //}
        }
	
	const unsigned int runNumber = sel.runNumber_or_RandomRunNumber();
	//std::cout << runNumber << std::endl;
	
        // loop over weight systematics
        for (size_t wIdx = 0; wIdx < weightSystematics.size(); ++wIdx) {
          std::string &suffix = weightSystematics[wIdx];

          TLorentzVector lep;
          if (sel.electron().size() == 1 && sel.muon().size() == 0) {
            lep = sel.electron()[0].mom();
          } else if (sel.electron().size() == 0 && sel.muon().size() == 1) {
            lep = sel.muon()[0].mom();
          }


          double weight = 1;
          if (!isData) {
            weight *= sel.weight_mc()*sel.weight_pileup();	    	    
            weight *= sampleXsection.getXsection(channel);
	    
            if(sel.weight_Sherpa22_corr())	weight *= sel.weight_Sherpa22_corr();

			//weight *= wjetsSF(sel, suffix); 
	    
            double pdfw = 1.0;
            bool isPdf = false;
            std::string pdfname = "";
            int pdfvar = -1;
            if (suffix.find("pdf_") != std::string::npos) {
              isPdf = true;
              size_t last = suffix.rfind("_");
              size_t first = std::string("pdf_").size();
              pdfname = suffix.substr(first, last - first);
              pdfvar = atoi(suffix.substr(last+1).c_str());
              pdfw = mt.vf(pdfname)->at(pdfvar);
            }
            weight *= pdfw;

            if ( ( (channel == 410000) || (channel == 301528) || (channel == 301529) ||
                   (channel == 301530) || (channel == 301531) || (channel == 301532) ) && applyPtRew) {
              weight *= ttWeight(sel);
            }

            if ( ( (channel == 410000) || (channel == 301528) || (channel == 301529) ||
                   (channel == 301530) || (channel == 301531) || (channel == 301532) ) && applyEWK) {
              float MC_t_pt = sel.MC_t().Perp();
              float MC_t_eta = sel.MC_t().Eta();
              float MC_t_phi = sel.MC_t().Phi();
              float MC_t_m = sel.MC_t().M();
              float MC_tbar_pt = sel.MC_tbar().Perp();
              float MC_tbar_eta = sel.MC_tbar().Eta();
              float MC_tbar_phi = sel.MC_tbar().Phi();
              float MC_tbar_m = sel.MC_tbar().M();
              int initial_type = mt.i("initial_type");
              WeakCorr::ScaleFactor sf; 
              sf = ewkTool.getScaleFactor(MC_t_pt,    MC_t_eta,    MC_t_phi,    MC_t_m, \
                                          MC_tbar_pt, MC_tbar_eta, MC_tbar_phi, MC_tbar_m, \
                                          initial_type);
              double ewkW = sf.nominal;
              if (suffix == "ttEWK__1down" || suffix == "ttEWK__1down_Loose") {
                ewkW = sf.down;
              } else if (suffix == "ttEWK__1up" || suffix == "ttEWK__1up_Loose") {
                ewkW = sf.up;
              }
              weight *= ewkW;
            }
		
            double btagsf = 1.0;
            std::string pref = "weight_bTagSF_70";
            std::string pref_pt = "tjet_bTagSF_70";
            //std::string prefe = "weight_bTagSF_70_env";
            std::string prefe = "weight_bTagSF_70";
            if (_btags <= 0 && _btags > -10) {
              pref = "weight_trackjet_bTagSF_70";
              //prefe = "weight_trackjet_bTagSF_70_env";
              prefe = "weight_trackjet_bTagSF_70";
            } else if (_btags < 0 && _btags <= -10) { // tag leading track jet only
              pref = "tjet_bTagSF_70";
              prefe = "tjet_bTagSF_70";
            } else if (_btags > 0 && _btags >= 10) { // tag leading calo jet only
              pref = "jet_bTagSF_70";
              prefe = "jet_bTagSF_70";
            }
            if (std::fabs(_btags) < 10) {
              size_t last = suffix.find("__1");
              btagsf = 1;
              if (suffix.find("btag") != std::string::npos && suffix.find("pt") != std::string::npos) {
                // split it in pt
                if (suffix.find("btagbSF_") != std::string::npos) {
                  size_t first = std::string("btagbSF_").size();
                  size_t first_pt = suffix.find("pt");
                  int eig = atoi(suffix.substr(first, last - first - 4).c_str());
                  int ptr = atoi(suffix.substr(first_pt+2, last - first_pt - 2).c_str());
                  for (size_t bidx = 0; bidx < mt.vf("tjet_mv2c10")->size(); ++bidx) {
                    TLorentzVector pb;
                    pb.SetPtEtaPhiE(mt.vf("tjet_pt")->at(bidx), mt.vf("tjet_eta")->at(bidx), mt.vf("tjet_phi")->at(bidx), mt.vf("tjet_e")->at(bidx));
                    if (mt.vf("tjet_pt")->at(bidx) > 10e3 &&
                      std::fabs(mt.vf("tjet_eta")->at(bidx)) < 2.5 && mt.vi("tjet_numConstituents")->at(bidx) >= 2) {
                      int ptbin = 1;
                      if (mt.vf("tjet_pt")->at(bidx) < 50e3) {
                        ptbin = 1;
                      } else if (mt.vf("tjet_pt")->at(bidx) < 100e3) {
                        ptbin = 2;
                      } else {
                        ptbin = 3;
                      }
                      if (ptr == ptbin) {
                        if (suffix.find("1up") != std::string::npos) {
                          btagsf *= mt.vvf(pref_pt+"_eigen_B_up")->at(bidx).at(eig);
                        } else if (suffix.find("1down") != std::string::npos) {
                          btagsf *= mt.vvf(pref_pt+"_eigen_B_down")->at(bidx).at(eig);
                        } else {
                          btagsf *= mt.vf(pref_pt)->at(bidx);
                        }
                      } else {
			btagsf *= mt.vf(pref_pt)->at(bidx);
		      }
                    }
                  }
                } else if (suffix.find("btagcSF_") != std::string::npos) {
                  size_t first = std::string("btagcSF_").size();
                  size_t first_pt = suffix.find("pt");
                  int eig = atoi(suffix.substr(first, last - first - 4).c_str());
                  int ptr = atoi(suffix.substr(first_pt+2, last - first_pt - 2).c_str());
                  for (size_t bidx = 0; bidx < mt.vf("tjet_mv2c10")->size(); ++bidx) {
                    TLorentzVector pb;
                    pb.SetPtEtaPhiE(mt.vf("tjet_pt")->at(bidx), mt.vf("tjet_eta")->at(bidx), mt.vf("tjet_phi")->at(bidx), mt.vf("tjet_e")->at(bidx));
                    if (mt.vf("tjet_pt")->at(bidx) > 10e3 &&
                      std::fabs(mt.vf("tjet_eta")->at(bidx)) < 2.5 && mt.vi("tjet_numConstituents")->at(bidx) >= 2) {
                      int ptbin = 1;
                      if (mt.vf("tjet_pt")->at(bidx) < 50e3) {
                        ptbin = 1;
                      } else if (mt.vf("tjet_pt")->at(bidx) < 100e3) {
                        ptbin = 2;
                      } else {
                        ptbin = 3;
                      }
                      if (ptr == ptbin) {
                        if (suffix.find("1up") != std::string::npos) {
                          btagsf *= mt.vvf(pref_pt+"_eigen_C_up")->at(bidx).at(eig);
                        } else if (suffix.find("1down") != std::string::npos) {
                          btagsf *= mt.vvf(pref_pt+"_eigen_C_down")->at(bidx).at(eig);
                        } else {
                          btagsf *= mt.vf(pref_pt)->at(bidx);
                        }
                      } else {
			btagsf *= mt.vf(pref_pt)->at(bidx);
		      }
                    }
                  }
                } else if (suffix.find("btaglSF_") != std::string::npos) {
                  size_t first = std::string("btaglSF_").size();
                  size_t first_pt = suffix.find("pt");
                  int eig = atoi(suffix.substr(first, last - first - 4).c_str());
                  int ptr = atoi(suffix.substr(first_pt+2, last - first_pt - 2).c_str());
                  for (size_t bidx = 0; bidx < mt.vf("tjet_mv2c10")->size(); ++bidx) {
                    TLorentzVector pb;
                    pb.SetPtEtaPhiE(mt.vf("tjet_pt")->at(bidx), mt.vf("tjet_eta")->at(bidx), mt.vf("tjet_phi")->at(bidx), mt.vf("tjet_e")->at(bidx));
                    if (mt.vf("tjet_pt")->at(bidx) > 10e3 &&  
                      std::fabs(mt.vf("tjet_eta")->at(bidx)) < 2.5 && mt.vi("tjet_numConstituents")->at(bidx) >= 2) {
                      int ptbin = 1;
                      if (mt.vf("tjet_pt")->at(bidx) < 50e3) {
                        ptbin = 1;
                      } else if (mt.vf("tjet_pt")->at(bidx) < 100e3) {
                        ptbin = 2;
                      } else {
                        ptbin = 3;
                      }
                      if (ptr == ptbin) {
                        if (suffix.find("1up") != std::string::npos) {
                          btagsf *= mt.vvf(pref_pt+"_eigen_Light_up")->at(bidx).at(eig);
                        } else if (suffix.find("1down") != std::string::npos) {
                          btagsf *= mt.vvf(pref_pt+"_eigen_Light_down")->at(bidx).at(eig);
                        } else {
                          btagsf *= mt.vf(pref_pt)->at(bidx);
                        }
                      } else {
			btagsf *= mt.vf(pref_pt)->at(bidx);
		      }
                    }
                  }
                } else {
                  btagsf = 0;
                }
              } else if (suffix.find("btagbSF_") != std::string::npos) {
                size_t first = std::string("btagbSF_").size();
                int eig = atoi(suffix.substr(first, last - first).c_str());
                for (size_t bidx = 0; bidx < mt.vf("tjet_mv2c10")->size(); ++bidx) {
                  TLorentzVector pb;
                  pb.SetPtEtaPhiE(mt.vf("tjet_pt")->at(bidx), mt.vf("tjet_eta")->at(bidx), mt.vf("tjet_phi")->at(bidx), mt.vf("tjet_e")->at(bidx));
                  if (mt.vf("tjet_pt")->at(bidx) > 10e3 &&
                    std::fabs(mt.vf("tjet_eta")->at(bidx)) < 2.5 && mt.vi("tjet_numConstituents")->at(bidx) >= 2) {
                    if (suffix.find("1up") != std::string::npos) {
                      btagsf *= mt.vvf(pref_pt+"_eigen_B_up")->at(bidx).at(eig);
                    } else if (suffix.find("1down") != std::string::npos) {
                      btagsf *= mt.vvf(pref_pt+"_eigen_B_down")->at(bidx).at(eig);
                    } else {
                      btagsf *= mt.vf(pref_pt)->at(bidx);
                    }
                  }
                }
              } else if (suffix.find("btagcSF_") != std::string::npos) {
                size_t first = std::string("btagcSF_").size();
                int eig = atoi(suffix.substr(first, last - first).c_str());
                for (size_t bidx = 0; bidx < mt.vf("tjet_mv2c10")->size(); ++bidx) {
                  TLorentzVector pb;
                  pb.SetPtEtaPhiE(mt.vf("tjet_pt")->at(bidx), mt.vf("tjet_eta")->at(bidx), mt.vf("tjet_phi")->at(bidx), mt.vf("tjet_e")->at(bidx));
                  if (mt.vf("tjet_pt")->at(bidx) > 10e3 && 
                    std::fabs(mt.vf("tjet_eta")->at(bidx)) < 2.5 && mt.vi("tjet_numConstituents")->at(bidx) >= 2) {
                    if (suffix.find("1up") != std::string::npos) {
                      btagsf *= mt.vvf(pref_pt+"_eigen_C_up")->at(bidx).at(eig);
                    } else if (suffix.find("1down") != std::string::npos) {
                      btagsf *= mt.vvf(pref_pt+"_eigen_C_down")->at(bidx).at(eig);
                    } else {
                      btagsf *= mt.vf(pref_pt)->at(bidx);
                    }
                  }
                }
              } else if (suffix.find("btaglSF_") != std::string::npos) {
                size_t first = std::string("btaglSF_").size();
                int eig = atoi(suffix.substr(first, last - first).c_str());
                for (size_t bidx = 0; bidx < mt.vf("tjet_mv2c10")->size(); ++bidx) {
                  TLorentzVector pb;
                  pb.SetPtEtaPhiE(mt.vf("tjet_pt")->at(bidx), mt.vf("tjet_eta")->at(bidx), mt.vf("tjet_phi")->at(bidx), mt.vf("tjet_e")->at(bidx));
                  if (mt.vf("tjet_pt")->at(bidx) > 10e3 &&
                    std::fabs(mt.vf("tjet_eta")->at(bidx)) < 2.5 && mt.vi("tjet_numConstituents")->at(bidx) >= 2) {
                    if (suffix.find("1up") != std::string::npos) {
                      btagsf *= mt.vvf(pref_pt+"_eigen_Light_up")->at(bidx).at(eig);
                    } else if (suffix.find("1down") != std::string::npos) {
                      btagsf *= mt.vvf(pref_pt+"_eigen_Light_down")->at(bidx).at(eig);
                    } else {
                      btagsf *= mt.vf(pref_pt)->at(bidx);
                    }
                  }
                }
              } else if (suffix.find("btageSF_0") != std::string::npos) {
                if (suffix.find("1up") != std::string::npos) {
                  btagsf = mt.f(prefe+"_extrapolation_up");
                } else {
                  btagsf = mt.f(prefe+"_extrapolation_down");
                }
              } else if (suffix.find("btageSF_1") != std::string::npos) {
                if (suffix.find("1up") != std::string::npos) {
                  btagsf = mt.f(prefe+"_extrapolation_from_charm_up");
                } else {
                  btagsf = mt.f(prefe+"_extrapolation_from_charm_down");
                }
              } else {
                //btagsf = mt.f(pref);
                for (size_t bidx = 0; bidx < mt.vf("tjet_mv2c10")->size(); ++bidx) {
                  TLorentzVector pb;
                  pb.SetPtEtaPhiE(mt.vf("tjet_pt")->at(bidx), mt.vf("tjet_eta")->at(bidx), mt.vf("tjet_phi")->at(bidx), mt.vf("tjet_e")->at(bidx));
                  if (mt.vf("tjet_pt")->at(bidx) > 10e3 &&
                    std::fabs(mt.vf("tjet_eta")->at(bidx)) < 2.5 && mt.vi("tjet_numConstituents")->at(bidx) >= 2) {
                    btagsf *= mt.vf(pref_pt)->at(bidx);
                  }
                }
              }
            } else if (std::fabs(_btags) >= 10) {
              if (mt.vf("tjet_pt")->size() > 0) {
                size_t last = suffix.find("__1");
                if (suffix.find("btagbSF_") != std::string::npos) {
                  size_t first = std::string("btagbSF_").size();
                  int eig = atoi(suffix.substr(first, last - first).c_str());
                  if (suffix.find("1up") != std::string::npos) {
                    btagsf = mt.vvf(pref+"_eigen_B_up")->at(0)[eig];
                  } else {
                    btagsf = mt.vvf(pref+"_eigen_B_down")->at(0)[eig];
                  }
                } else if (suffix.find("btagcSF_") != std::string::npos) {
                  size_t first = std::string("btagcSF_").size();
                  int eig = atoi(suffix.substr(first, last - first).c_str());
                  if (suffix.find("1up") != std::string::npos) {
                    btagsf = mt.vvf(pref+"_eigen_C_up")->at(0)[eig];
                  } else {
                    btagsf = mt.vvf(pref+"_eigen_C_down")->at(0)[eig];
                  }
                } else if (suffix.find("btaglSF_") != std::string::npos) {
                  size_t first = std::string("btaglSF_").size();
                  int eig = atoi(suffix.substr(first, last - first).c_str());
                  if (suffix.find("1up") != std::string::npos) {
                    btagsf = mt.vvf(pref+"_eigen_Light_up")->at(0)[eig];
                  } else {
                    btagsf = mt.vvf(pref+"_eigen_Light_down")->at(0)[eig];
                  }
                } else if (suffix.find("btageSF_0") != std::string::npos) {
                  if (suffix.find("1up") != std::string::npos) {
                    btagsf = mt.vf(prefe+"_syst_extrapolation_up")->at(0);
                  } else {
                    btagsf = mt.vf(prefe+"_syst_extrapolation_down")->at(0);
                  }
                } else if (suffix.find("btageSF_1") != std::string::npos) {
                  if (suffix.find("1up") != std::string::npos) {
                    btagsf = mt.vf(prefe+"_syst_extrapolation_from_charm_up")->at(0);
                  } else {
                    btagsf = mt.vf(prefe+"_syst_extrapolation_from_charm_down")->at(0);
                  }
                } else {
                  btagsf = mt.vf(pref)->at(0);
                }
              }
            }

            weight *= btagsf;
            if (suffix == "eTrigSF__1up" || suffix == "eTrigSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_EL_SF_Trigger_UP");
            } else if (suffix == "eTrigSF__1down" || suffix == "eTrigSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_EL_SF_Trigger_DOWN");
            } else if (suffix == "eRecoSF__1up" || suffix == "eRecoSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_EL_SF_Reco_UP");
            } else if (suffix == "eRecoSF__1down" || suffix == "eRecoSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_EL_SF_Reco_DOWN");
            } else if (suffix == "eIDSF__1up" || suffix == "eIDSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_EL_SF_ID_UP");
            } else if (suffix == "eIDSF__1down" || suffix == "eIDSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_EL_SF_ID_DOWN");
            } else if (suffix == "eIsolSF__1up" || suffix == "eIsolSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_EL_SF_Isol_UP");
            } else if (suffix == "eIsolSF__1down" || suffix == "eIsolSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_EL_SF_Isol_DOWN");
            } else if (suffix == "muTrigStatSF__1up" || suffix == "muTrigStatSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_Trigger_STAT_UP");
            } else if (suffix == "muTrigStatSF__1down" || suffix == "muTrigStatSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_Trigger_STAT_DOWN");
            } else if (suffix == "muIDStatSF__1up" || suffix == "muIDStatSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_ID_STAT_UP");
            } else if (suffix == "muIDStatSF__1down" || suffix == "muIDStatSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_ID_STAT_DOWN");
            } else if (suffix == "muIsolStatSF__1up" || suffix == "muIsolStatSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_Isol_STAT_UP");
            } else if (suffix == "muIsolStatSF__1down" || suffix == "muIsolStatSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_Isol_STAT_DOWN");
            } else if (suffix == "muTrigSystSF__1up" || suffix == "muTrigSystSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_Trigger_SYST_UP");
            } else if (suffix == "muTrigSystSF__1down" || suffix == "muTrigSystSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_Trigger_SYST_DOWN");
            } else if (suffix == "muIDSystSF__1up" || suffix == "muIDSystSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_ID_SYST_UP");
            } else if (suffix == "muIDSystSF__1down" || suffix == "muIDSystSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_ID_SYST_DOWN");
            } else if (suffix == "muIsolSystSF__1up" || suffix == "muIsolSystSF__1up_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_Isol_SYST_UP");
            } else if (suffix == "muIsolSystSF__1down" || suffix == "muIsolSystSF__1down_Loose") {
              weight *= mt.f("weight_leptonSF_MU_SF_Isol_SYST_DOWN");
            } else {
              weight *= sel.weight_leptonSF();
            }

            if (sumOfWeights[channel] != 0 && !isPdf) {
              weight /= sumOfWeights[channel];
            } else if (isPdf) {
              weight /= PDFsumOfWeights[channel][pdfname][pdfvar];
            }
		
	  }//!isData
      
          // this applies b-tagging early
          int nBtagged = 0; 
          int nLeadTagged = 0;
          if (_btags > 0) { //calo b-tagging
            for (size_t bidx = 0; bidx < sel.jet().size(); ++bidx) {
              if (sel.jet()[bidx].btag_mv2c20_70()) {
                nBtagged += 1;	
              }
            }
            if (sel.jet()[0].btag_mv2c20_70()) nLeadTagged++;
          } else if (_btags <= 0) { //track b-tagging
            for (size_t bidx = 0; bidx < mt.vf("tjet_mv2c10")->size(); ++bidx) {
              TLorentzVector pb;
              pb.SetPtEtaPhiE(mt.vf("tjet_pt")->at(bidx), mt.vf("tjet_eta")->at(bidx), mt.vf("tjet_phi")->at(bidx), mt.vf("tjet_e")->at(bidx));
              if (mt.vf("tjet_mv2c10")->at(bidx) > 0.6455 && mt.vf("tjet_pt")->at(bidx) > 10e3 && 
                  std::fabs(mt.vf("tjet_eta")->at(bidx)) < 2.5 && mt.vi("tjet_numConstituents")->at(bidx) >= 2) {
                nBtagged += 1;
              }
            }
            if (mt.vf("tjet_pt")->size() > 0) {
              if (mt.vf("tjet_mv2c10")->at(0) > 0.6455 && mt.vf("tjet_pt")->at(0) > 10e3 &&
                  std::fabs(mt.vf("tjet_eta")->at(0)) < 2.5 && mt.vi("tjet_numConstituents")->at(0) >= 2) {
                nLeadTagged += 1;
              }
            }
          }
          if (_btags != 0 && std::fabs(_btags) < 10) {
            if (nBtagged < abs(_btags))
              continue;
          } else if (_btags != 0 && std::fabs(_btags) >= 10) {
            if (nLeadTagged == 0)
              continue;
          }
	  
          procEvents++;
	  if (analysis=="AnaTtresSL") { // signal region
	    for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresSL*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresSL*>(vec_analysis[iAna]))->isBoosted());
	        if (runMM) {

		   if(nBtagged==0){		     
	  	     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);// use w+jets rates (>=2 jets) everywhere for the time being
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);// use w+jets rates (>=2 jets) everywhere for the time being

		   }else{		      
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);// use w+jets rates (>=2 jets) everywhere for the time being
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);// use w+jets rates (>=2 jets) everywhere for the time being
		     
		   }//if
		}//runMM

                vec_analysis[iAna]->run(sel, weight, suffix);
	    }
          } 
	    else if (analysis=="AnaTtresSL_QCDCR2j_2015") { // QCD control region, >=2 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDCR2j_2015(sel, weight, suffix);
            }
          	  
	  }  
            else if (analysis=="AnaTtresSL_QCDCR2j_2016") { // QCD control region, >=2 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDCR2j_2016(sel, weight, suffix);
            }
          	  
	  }           
           else if (analysis=="AnaTtresSL_QCDVR1_2j_2015") { // QCD vvalidation region for muon channel, >=2 jets; 3< |d_0Sig| < 5
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDVR1_2j_2015(sel, weight, suffix);
            }
	 
	  }

          else if (analysis=="AnaTtresSL_QCDVR2_2j_2015") { // QCD vvalidation region for muon channel, >=2 jets; 3< |d_0Sig| < 5
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
                bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
                bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());

                if (runMM) {

                   if(nBtagged==0){
                     if(iselect)        weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }else{
                     if(iselect)        weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }//if  
               }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDVR2_2j_2015(sel, weight, suffix);
            }

          }
          
           else if (analysis=="AnaTtresSL_QCDVR1_2j_2016") { // QCD vvalidation region, >=2 jets; 3< |d_0Sig| < 5
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDVR1_2j_2016(sel, weight, suffix);
            }
	 
	  }
         
          else if (analysis=="AnaTtresSL_QCDVR2_2j_2016") { // QCD vvalidation region, >=2 jets; 3< |d_0Sig| < 5
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
                bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
                bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());

                if (runMM) {

                   if(nBtagged==0){
                     if(iselect)        weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }else{
                     if(iselect)        weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }//if  
               }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDVR2_2j_2016(sel, weight, suffix);
            }

          }


            else if (analysis=="AnaTtresSL_WjetsCR2j_2015") { // W+jets enriched control region, >=2 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_WjetsCR2j_2015(sel, weight, suffix);
            }

	  } 

           else if (analysis=="AnaTtresSL_QCDSR2j_2016") { // W+jets enriched control region, >=2 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDSR2j_2016(sel, weight, suffix);
            }

	  }    
           else if (analysis=="AnaTtresSL_QCDCR4j_2016") { // resolved control region, >=4 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDCR4j_2016(sel, weight, suffix);
            }

 	  } 

        else if (analysis=="AnaTtresSL_QCDCR4j_2015") { // resolved control region, >=4 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
                bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
                bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());

                if (runMM) {

                   if(nBtagged==0){
                     if(iselect)        weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }else{
                     if(iselect)        weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }//if  
               }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDCR4j_2015(sel, weight, suffix);
            }

          }
              

            else if (analysis=="AnaTtresSL_SR4j_2015") { // resolved channel (4jets, no chi2) and boosted channel (large jet Pt > 300 GeV)
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_SR4j_2015(sel, weight, suffix);
            }	 
	  
	  } 
            else if (analysis=="AnaTtresSL_QCDVR1_4j_2015") { // QCD vvalidation region, >=4 jets; 3< |d_0Sig| < 5
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDVR1_4j_2015(sel, weight, suffix);
            }	
           }

           else if (analysis=="AnaTtresSL_QCDVR2_4j_2015") { // QCD vvalidation region, >=4 jets; 3< |d_0Sig| < 5
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
                bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
                bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());

                if (runMM) {

                   if(nBtagged==0){
                     if(iselect)        weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }else{
                     if(iselect)        weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }//if  
               }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDVR2_4j_2015(sel, weight, suffix);
            }
 	  
	  }

            else if (analysis=="AnaTtresSL_QCDVR1_4j_2016") { // QCD vvalidation region, >=4 jets; 3< |d_0Sig| < 5
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
	        bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
	        bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());
	        
	        if (runMM) {
                   
		   if(nBtagged==0){
		     if(iselect)	weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
	  	     
		  }else{  
		     if(iselect)	weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     else		weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
		     
		  }//if  
	       }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDVR1_4j_2016(sel, weight, suffix);
            }	 
	  
	  }

          else if (analysis=="AnaTtresSL_QCDVR2_4j_2016") { // QCD vvalidation region, >=4 jets; 3< |d_0Sig| < 5
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
                bool iselect((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isElectron());
                bool isboost((dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->isBoosted());

                if (runMM) {

                   if(nBtagged==0){
                     if(iselect)        weight = MM_0b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_0b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }else{
                     if(iselect)        weight = MM_1b_re->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);
                     else               weight = MM_1b_rmu->getMMweights(sel, runMM_StatErr, iselect, 0, runNumber);

                  }//if  
               }//runMM

               (dynamic_cast<AnaTtresMM*>(vec_analysis[iAna]))->runMatrixMethod_QCDVR2_4j_2016(sel, weight, suffix);
            }

          }
          

	    else if (analysis=="AnaTtresWQCDreal_2015") { // to generate the real rates, >=2 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) 
		(dynamic_cast<AnaTtresQCD*>(vec_analysis[iAna]))->runRealRateWQCDCR_2015(sel, weight, suffix);
	  }
            else if (analysis=="AnaTtresQCDreal_2015") { // to generate the real rates, >=4 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna)
                (dynamic_cast<AnaTtresQCD*>(vec_analysis[iAna]))->runRealRateQCDCR_2015(sel, weight, suffix);
          }
            else if (analysis=="AnaTtresWQCDreal_2016") { // to generate the real rates, >=2 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) 
		(dynamic_cast<AnaTtresQCD*>(vec_analysis[iAna]))->runRealRateWQCDCR_2016(sel, weight, suffix);
	  }
            else if (analysis=="AnaTtresQCDreal_2016") { // to generate the real rates, >=4 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna)
                (dynamic_cast<AnaTtresQCD*>(vec_analysis[iAna]))->runRealRateQCDCR_2016(sel, weight, suffix);
          }
            else if (analysis=="AnaTtresQCDfake_2015") { // to generate the fake rates, >=4 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) 
		(dynamic_cast<AnaTtresQCD*>(vec_analysis[iAna]))->runFakeRateQCDCR_2015(sel, weight, suffix);
          }
            else if (analysis=="AnaTtresQCDfake_2016") { // to generate the fake rates, >=4 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) 
		(dynamic_cast<AnaTtresQCD*>(vec_analysis[iAna]))->runFakeRateQCDCR_2016(sel, weight, suffix);
          }
            else if (analysis=="AnaTtresWQCDfake_2015") { // to generate the fake rates, >=2 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) 
		(dynamic_cast<AnaTtresQCD*>(vec_analysis[iAna]))->runFakeRateWQCDCR_2015(sel, weight, suffix);
          } 
            else if (analysis=="AnaTtresWQCDfake_2016") { // to generate the fake rates, >=2 jets
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) 
		(dynamic_cast<AnaTtresQCD*>(vec_analysis[iAna]))->runFakeRateWQCDCR_2016(sel, weight, suffix);
          } 
            else if (analysis=="AnaTtresSLMtt") {
            for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) 
                vec_analysis[iAna]->run(sel, weight, suffix);
          } //if  
        } // end of loop over weight systematics
      } // end of loop over entries
    } // end of loop over lepton modes (tight and loose)
  } // end of loop over systematics

  for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
    std::cout << "Removed " << vec_analysis[iAna]->getNduplicate() << " duplicate entries" << std::endl;
    vec_analysis[iAna]->terminate();
    delete vec_analysis[iAna];
  }

  vec_analysis.clear();

  std::cout << "Called analysis code " << procEvents << " times." << std::endl;

  return 0;
}


