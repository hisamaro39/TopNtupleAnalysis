#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
#include <iterator>

#include "TFile.h"
#include "TH1F.h"
#include "TH1.h"

#include <iomanip>
#include <sstream>

#include <signal.h>

#include "Hist.h"
#include "SystematicImplementation.h"
#include "utils.h"
#include "SampleSet.h"
#include "SystematicCalculation.h"

#include "ParseUtils.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TGraphErrors.h"

using namespace std;

int main(int argc, char **argv) {
  signal(SIGSEGV, handler);

  try {

    int help = 0;
    string h_input = "llargeJetPt";
    string h_files = "out_ttbarMatched_e.root,out_mcgen3m_e.root";
    string h_titles = "t#bar{t} matched,MC gen var. matched";
    string _outfile = "";
    int underflow = 0;
    string _extraText = "";
    string yTitle = "";
    float xMax = -9999;
    float xMin = -9999;
    float lumi = 5;
    std::string config = "";

    static struct extendedOption extOpt[] = {
        {"help",            no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"histogram",       required_argument,     0, 'h', "Histograms to compare", &h_input, extendedOption::eOTString},
        {"files",           required_argument,     0, 'f', "Comma-separated list of files to compare.", &h_files, extendedOption::eOTString},
        {"titles",          required_argument,     0, 't', "Comma-separated list of titles corresponding to files.", &h_titles, extendedOption::eOTString},
        {"outfile",         required_argument,     0, 'o', "Output file.", &_outfile, extendedOption::eOTString},
        {"underflow",       required_argument,     0, 'U', "Include underflow (0/1)", &underflow, extendedOption::eOTInt},
        {"extraText",       required_argument,     0, 'T', "Extra text to add in the plot.", &_extraText, extendedOption::eOTString},
        {"yTitle",          required_argument,     0, 'y', "Y title.", &yTitle, extendedOption::eOTString},
        {"xMax",            required_argument,     0, 'X', "Limit X maximum value.", &xMax, extendedOption::eOTFloat},
        {"xMin",            required_argument,     0, 'x', "Limit X minimum value.", &xMin, extendedOption::eOTFloat},
        {"lumi",            required_argument,     0, 'l', "Luminosity value to show.", &lumi, extendedOption::eOTFloat},
        {"config",          required_argument,     0, 'C', "Config file.", &config, extendedOption::eOTString},

        {0, 0, 0, 0, 0, 0, extendedOption::eOTInt}
      };

    
    if (!parseArguments(argc, argv, extOpt) || help) {
      dumpHelp("compare", extOpt, "plot\nPlot two histograms for comparison.\n");
      return 0;
    } else {
      std::cout << "Dumping options:" << std::endl;
      dumpOptions(extOpt);
    }

    lumi_scale = lumi;
    if (config != "")
      loadConfig(config.c_str());
    else
      loadConfig(std::string(argv[0]).substr(0, std::string(argv[0]).rfind('/'))+"/config.txt");

    vector<string> items;
    vector<string> titles;
    split(h_files, ',', items);
    split(h_titles, ',', titles);
    string histogram = h_input;

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(1);

    // for Data/MC comparison
    SampleSetConfiguration stackConfig;

    Color_t colorList[] = {kBlack, kBlue, kRed, kGreen, kYellow, kCyan, kMagenta};
    for (int i = 0; i < items.size(); ++i) {
      stackConfig.addType(Form("%d items[i]", i));
      stackConfig.add(Form("%d items[i]", i), items[i], titles[i], titles[i], titles[i],
                                              "PL", 1, colorList[i], 1,      0, 20, 1, "e1");
    }
    SystematicCalculator systCalc(stackConfig);
    systCalc.calculate(histogram);

    if (underflow) stackConfig.showUnderflow();
    if (xMax > -998.0) stackConfig.limitMaxX(xMax);
    if (xMin > -998.0) stackConfig.limitMinX(xMin);

    cout << "Yields:" << endl << endl;
    systCalc.printYields(stackConfig);
  
    vector<string> extraText;
    string outfile = _outfile;
    string cname = "";
    if (items[0].find("_e.root") != string::npos) {
      cname = "e";
      extraText.push_back("e channel");
    } else if (items[0].find("_mu.root") != string::npos) {
      cname = "mu";
      extraText.push_back("#mu channel");
    } else if (items[0].find("_comb.root") != string::npos) {
      cname = "comb";
      extraText.push_back("e,#mu-channel");
    }
    if (outfile == "") {
      outfile = string("compare_")+histogram+string("_")+cname+string(".pdf");
    }
    extraText.push_back(_extraText);
    drawCompare(stackConfig, extraText, outfile, true, lumi);

  } catch (string s) {
    cout << "Crashed with exception: " << s << endl;
    dumpTrace();
  }

  return 0;
}

