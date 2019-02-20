/**
 * @brief Main executable to run over input mini flat ntuples,
 * read events in an object-oriented way and
 * run the event analysis class that derives from Analysis.
 * The main objective is to be minimal and to allow users to run this part
 * of the code in their laptops or using the E4 MapReduce hadoop cluster.
 *
 * @author Danilo Enoque Ferreira de Lima <dferreir@mail.cern.ch>
 * @author Hendrik Esch <esch@.cern.ch>
 */

#include "TopNtupleAnalysis/Event.h"
#include "TChain.h"
#include "TopNtupleAnalysis/MiniTree.h"

#include <iostream>

#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresSL.h"
#include "TopNtupleAnalysis/AnaTuDoSL.h"
#include "TopNtupleAnalysis/AnaTuDoTtresResolved.h"
#include "TopNtupleAnalysis/AnaTuDoTtresBoosted.h"

// #include "TopDataPreparation/SampleXsection.h"

#include "TopNtupleAnalysis/ParseUtils.h"

#include <fstream>

#include "TROOT.h"
#include "TString.h"
#include "TInterpreter.h"

#include <sstream>

int main(int argc, char **argv) {

//   ROOT::Cintex::Cintex::Enable();

  gROOT->ProcessLine("#include <vector>");
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");

//   initEventCount();

  // input files
  int help = 0;
  int systFlag = 0;
  int isAtlFastII = 0;
  std::string outFile = "";
  std::string files = "output.root";
  std::string analysis = "AnaTuDoSL";
  std::string channel = "mu";

  static struct extendedOption extOpt[] = {
        {"help",          no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"systFlag",         required_argument,     0, 's', "Systematics tag. 0 means nominal, many other integer values are possible.", &systFlag, extendedOption::eOTInt},
        {"atlFastII",         required_argument,     0, 'a', "Is this AtlFastII? (0/1)", &isAtlFastII, extendedOption::eOTInt},
        {"files",         required_argument,     0, 'f', "Input list of comma-separated D3PD files to apply the selection on", &files, extendedOption::eOTString},
        {"analysis",   required_argument,     0, 'A', "Analysis to run. Choices: AnaTuDoSL", &analysis, extendedOption::eOTString},
        {"channel",   required_argument,     0, 'c', "Decay channel. Choices: el, mu", &channel, extendedOption::eOTString},
        {"output",   required_argument,     0, 'o', "Output file (optional)", &outFile, extendedOption::eOTString},

        {0, 0, 0, 0, 0, 0, extendedOption::eOTInt}
      };

    bool printhelp = parseArguments(argc, argv, extOpt);
    
    bool isData = false;
    if (systFlag == -1) {
        isData = true;
    }
    
    if (outFile == "") {
        TString outFileName = files;
        int rootIndex = outFileName.Index(".root");
        outFileName.Insert(rootIndex, Form("_%s_output_%i", channel.c_str(), systFlag));
        outFile = outFileName;
    }

    if (!printhelp || help) {
        dumpHelp("read", extOpt, "read\nRead selected events after preselection and generate histograms.\n");
        return 0;
    } else {
        std::cout << "Dumping options:" << std::endl;
        dumpOptions(extOpt);
    }

  // parse file list
  std::vector<std::string> fileList;
  if ( ((files.find("input") != std::string::npos) && (files.find(".txt") != std::string::npos)) ) {
    std::cout << "Using file given as text list." << std::endl;
    std::vector<std::string> inputList;
    // split by ','
    std::string argStr = files;
    for (size_t i = 0,n; i <= argStr.length(); i=n+1) {
      n = argStr.find_first_of(',',i);
      if (n == std::string::npos)
        n = argStr.length();
      std::string tmp = argStr.substr(i,n-i);
      if (tmp != "")
        inputList.push_back(tmp);
    }

    for (int k = 0; k < inputList.size(); ++k) {

      ifstream f(inputList[k].c_str());
      std::string thePathStr;
      while (f.good()) {
        std::stringstream ss;
        f.get(*ss.rdbuf(), '\n');
        if (f.get() == EOF)
          break;

       thePathStr = ss.str();
         if (thePathStr != "") {
          // bug in Ganga/Grid nodes: \\n is substituted instead of a newline
          size_t idx = std::string::npos;
          do {
            idx = thePathStr.find("\\n");
            std::string aFile = thePathStr.substr(0, idx);
            if (aFile.find(".root") != std::string::npos)
              fileList.push_back(aFile);
            if (idx != std::string::npos) {
              thePathStr = thePathStr.substr(idx+2);
            }
          } while (idx != std::string::npos);
        }
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
  }

 
  std::string treeName = "nominal";
  
  //TODO
  // include possibility to switch from "nominal" tree to tree with systematic variation
  
  MiniTree *mt = new MiniTree(false, fileList[0].c_str(), treeName);
  for (int k = 1; k < fileList.size(); ++k) {
    mt->addFileToRead(fileList[k]);
  }

//   // get output files
//   std::string outStr = outFile;
//   std::vector<std::string> outList;
//   for (size_t i = 0,n; i <= outStr.length(); i=n+1) {
//     n = outStr.find_first_of(',',i);
//     if (n == std::string::npos)
//       n = outStr.length();
//     std::string tmp = outStr.substr(i,n-i);
//     outList.push_back(tmp);
//   }

    std::vector<Analysis *> vec_analysis;
    if (analysis == "AnaTuDoSL") {
        if (channel == "el") {
            vec_analysis.push_back(new AnaTuDoSL(outFile, true, mt )); // electron channel
        } else if (channel == "mu") {
            vec_analysis.push_back(new AnaTuDoSL(outFile, false, mt )); // muon channel
        } else {
            std::cerr << "Error! Decay channel not known. Choose either el or mu. Aborting!" << std::endl;
            exit(0);
        }
    } else if (analysis == "AnaTuDoTtresResolved") {
        if (channel == "el") {
            vec_analysis.push_back(new AnaTuDoTtresResolved(outFile, true, mt )); // electron channel
        } else if (channel == "mu") {
            vec_analysis.push_back(new AnaTuDoTtresResolved(outFile, false, mt )); // muon channel
        } else {
            std::cerr << "Error! Decay channel not known. Choose either el or mu. Aborting!" << std::endl;
            exit(0);
        }
    } else if (analysis == "AnaTuDoTtresBoosted") {
        if (channel == "el") {
            vec_analysis.push_back(new AnaTuDoTtresBoosted(outFile, true, mt )); // electron channel
        } else if (channel == "mu") {
            vec_analysis.push_back(new AnaTuDoTtresBoosted(outFile, false, mt )); // muon channel
        } else {
            std::cerr << "Error! Decay channel not known. Choose either el or mu. Aborting!" << std::endl;
            exit(0);
        }
    }

  Event sel; // selected objects

//   SampleXsection sampleXsection;
//   sampleXsection.readFromFile("../TopDataPreparation/data/XSection-MC15-13TeV-fromSusyGrp.data");

  // retrieve, list of sum of weights
  std::map<int, float> sumOfWeights;
  TChain t_sumWeights("sumWeights");
  for (int k = 0; k < fileList.size(); ++k) {
    t_sumWeights.Add(fileList[k].c_str());
  }
  int dsid;
  float value;
  t_sumWeights.SetBranchAddress("dsid", &dsid);
  t_sumWeights.SetBranchAddress("totalEventsWeighted", &value);
  for (int k = 0; k < t_sumWeights.GetEntries(); ++k) {
    t_sumWeights.GetEntry(k);
    if (sumOfWeights.find(dsid) == sumOfWeights.end())
      sumOfWeights[dsid] = 0;
    sumOfWeights[dsid] += value;
  }
 
  int channelNumber = -99;
 
  for (int k = 0; k < mt->GetEntries(); ++k) {
    if (k % 1000 == 0)
      std::cout << "Entry " << k << "/" << mt->GetEntries() << std::endl;

    mt->read(k, sel);

    channelNumber = sel.channelNumber();
    
    double weight = 1;
    if (!isData) {
//       weight *= sel.mcWeight() *sel.pileupWeight();
//       weight *= sampleXsection.getXsection(channel);
//       //weight /= getEventCountBeforeSkimming(channel);
//       if (sumOfWeights[channel] != 0)
//         weight /= sumOfWeights[channel]; // this will be the correct way of doing this
      // but keeping this commented as it ihas only been added in the trunk of AnalysisTop now
      // if you use a recent version of AnalysisTop, uncomment the last line
    }
    
    for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
      vec_analysis[iAna]->setIsData(isData);  
      vec_analysis[iAna]->run(sel, weight);
    }

//     std::cout << "channel " << channel << " and sumOfWeights " << sumOfWeights[channel] << " and cross-section " << sampleXsection.getXsection(channel) << std::endl;

  } // event loop ends here
  
  for (size_t iAna = 0; iAna < vec_analysis.size(); ++iAna) {
    vec_analysis[iAna]->terminate();
    // INFO: If using the HistoList tool, the following lines have to be commented out. Otherwise 
    // both HistogramSvc tools will try to close the same file --> Crashing code 
    delete vec_analysis[iAna];
  }
  // INFO: If using the HistoList tool, the following lines have to be commented out. Otherwise 
  // both HistogramSvc tools will try to close the same file --> Crashing code 
  vec_analysis.clear();

  return 0;
}

