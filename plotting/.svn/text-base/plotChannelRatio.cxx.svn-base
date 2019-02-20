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
  //signal(SIGSEGV, handler);

  try {

    int help = 0;
    string channel = "e";
    string prefix = "boosted";
    string h_input = "lepPt";
    int mcOnly = 0;
    string _outfile = "";
    int verbose = 0;
    int underflow = 0;
    string _extraText = "boosted";
    string yTitle = "";
    string xTitle = "";
    float xMax = -9999;
    float xMin = -9999;
    int mustBeBigger = 0;
    int posLegend = 0;
    int rebin = 1;
    float yMax = -1;
    float yMin = -1;
    int arrow = 0;
    int stamp = 0;
    float lumi = 0;
    std::string config = "";
    std::string saveTH1 = "";
    int normBinWidth = 0;
    int alternateStyle = 0;
    int _logY = 0;

    static struct extendedOption extOpt[] = {
        {"help",            no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"channel",         required_argument,     0, 'c', "Channel: one of e, mu, comb", &channel, extendedOption::eOTString},
        {"prefix",          required_argument,     0, 'p', "Prefix.", &prefix, extendedOption::eOTString},
        {"histogram",       required_argument,     0, 'h', "Histogram name.", &h_input, extendedOption::eOTString},
        {"mcOnly",          required_argument,     0, 'm', "Only look for the ttbarMatched histograms? (0/1)", &mcOnly, extendedOption::eOTInt},
        {"outfile",         required_argument,     0, 'o', "Output file.", &_outfile, extendedOption::eOTString},
        {"verbose",         required_argument,     0, 'v', "Verbose. (0/1)", &verbose, extendedOption::eOTInt},
        {"underflow",       required_argument,     0, 'U', "Include underflow (0/1)", &underflow, extendedOption::eOTInt},
        {"extraText",       required_argument,     0, 'T', "Extra text to add in the plot.", &_extraText, extendedOption::eOTString},
        {"yTitle",          required_argument,     0, 'E', "Y title.", &yTitle, extendedOption::eOTString},
        {"xTitle",          required_argument,     0, 't', "X title.", &xTitle, extendedOption::eOTString},
        {"xMax",            required_argument,     0, 'X', "Limit X maximum value.", &xMax, extendedOption::eOTFloat},
        {"xMin",            required_argument,     0, 'x', "Limit X minimum value.", &xMin, extendedOption::eOTFloat},
        {"mustBeBigger",    required_argument,     0, 'M', "Widen range of the Y axis in the ratio plots.", &mustBeBigger, extendedOption::eOTInt},
        {"posLegend",       required_argument,     0, 'L', "Move legend to the left.", &posLegend, extendedOption::eOTInt},
        {"rebin",           required_argument,     0, 'r', "Rebin by this factor.", &rebin, extendedOption::eOTInt},
        {"yMax",            required_argument,     0, 'Y', "Maximum of the Y axis.", &yMax, extendedOption::eOTFloat},
        {"yMin",            required_argument,     0, 'y', "Minimum of the Y axis.", &yMax, extendedOption::eOTFloat},
        {"arrow",           required_argument,     0, 'a', "Draw arrow.", &arrow, extendedOption::eOTInt},
        {"stamp",           required_argument,     0, 's', "0 = ATLAS Internal, 1 = ATLAS Preliminary.", &stamp, extendedOption::eOTInt},
        {"lumi",            required_argument,     0, 'l', "Luminosity value to show", &lumi, extendedOption::eOTFloat},
        {"config",          required_argument,     0, 'C', "Configuration file for items.", &config, extendedOption::eOTString},
        {"smoothen",        required_argument,     0, 'k', "Smoothen systematics.", &smooth, extendedOption::eOTInt},
        {"saveTH1",         required_argument,     0, 'R', "Save as ROOT files called hist_[name][sufix].root with all systs. Provide the sufix.", &saveTH1, extendedOption::eOTString},
        {"normBinWidth",    required_argument,     0, 'b', "Divide bin content by bin width?", &normBinWidth, extendedOption::eOTInt},
        {"alternateStyle",  required_argument,     0, 'A', "Alternative style", &alternateStyle, extendedOption::eOTInt},
        {"logY",            required_argument,     0, 'g', "Y axis in log?", &_logY, extendedOption::eOTInt},

        {0, 0, 0, 0, 0, 0, extendedOption::eOTInt}
      };

    
    if (!parseArguments(argc, argv, extOpt) || help) {
      dumpHelp(std::string(argv[0]), extOpt, "plotRatioChannels\nCalculate systematic uncertainties and plot ratio of the histogram in two channels.\n");
      return 0;
    } else {
      std::cout << "Dumping options:" << std::endl;
      dumpOptions(extOpt);
    }

    logY = _logY;
    lumi_scale = lumi;
    if (config != "")
      loadConfig(config.c_str());
    else
      loadConfig(std::string(argv[0]).substr(0, std::string(argv[0]).rfind('/'))+"/config.txt");

    _stamp = stamp;

    std::string histogram = "";
    histogram = h_input;

    vector<string> h_items;
    split(channel, '/', h_items);

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(1);

    shared_ptr<SystematicRatioCalculator> systRatCalcD;
    shared_ptr<SystematicRatioCalculator> systRatCalc;
    shared_ptr<SystematicRatioCalculator> systRatCalcMCAtNLO;
    shared_ptr<SystematicRatioRatioCalculator> systRatRatCalc;
    shared_ptr<SystematicRatioRatioCalculator> systRatRatCalcMCAtNLO;

    if (!mcOnly) {
      // for data eff. numerator
      SampleSetConfiguration stackConfigNumD = makeConfigurationDataOnly(prefix, h_items[0]);
      SystematicCalculator systCalcNumD(stackConfigNumD);
      systCalcNumD.calculate(histogram);
      if (xMax > -998.0) stackConfigNumD.limitMaxX(xMax);
      if (xMin > -998.0) stackConfigNumD.limitMinX(xMin);

      if (verbose) {
        cout << "Data numerator yields:" << endl << endl;
        systCalcNumD.printYields(stackConfigNumD);
      }

      // for data eff. denominator
      SampleSetConfiguration stackConfigDenD = makeConfigurationDataOnly(prefix, h_items[1]);
      SystematicCalculator systCalcDenD(stackConfigDenD);
      systCalcDenD.calculate(histogram);
      if (xMax > -998.0) stackConfigDenD.limitMaxX(xMax);
      if (xMin > -998.0) stackConfigDenD.limitMinX(xMin);

      if (verbose) {
        cout << "Data denominator yields:" << endl << endl;
        systCalcDenD.printYields(stackConfigDenD);
      }

      // calculate ratio
      systRatCalcD.reset(new SystematicRatioCalculator(systCalcNumD, systCalcDenD, false));
      systRatCalcD->calculate(histogram, histogram, false);

      cout << "Ratio nominal in data:" << systRatCalcD->_sr["Ratio"]._item[0].nominal << endl;
    }
  
    // for MC eff. numerator using PP
    SampleSetConfiguration stackConfigNum = makeConfigurationMCOnly(prefix, h_items[0]);
    SystematicCalculator systCalcNum(stackConfigNum);
    addAllSystematics(systCalcNum, prefix, h_items[0]);

    systCalcNum.calculate(histogram);

    if (xMax > -998.0) stackConfigNum.limitMaxX(xMax);
    if (xMin > -998.0) stackConfigNum.limitMinX(xMin);

    if (verbose) {
      cout << "MC numerator systs. (standard):" << endl << endl;
      systCalcNum.printSysts(stackConfigNum["MC"]);
      cout << "MC numerator yields (standard):" << endl << endl;
      systCalcNum.printYields(stackConfigNum);
    }
  
    // for MC eff. denominator
    SampleSetConfiguration stackConfigDen = makeConfigurationMCOnly(prefix, h_items[1]);
    SystematicCalculator systCalcDen(stackConfigDen);
    addAllSystematics(systCalcDen, prefix, h_items[1]);
  
    systCalcDen.calculate(histogram);

    if (xMax > -998.0) stackConfigDen.limitMaxX(xMax);
    if (xMin > -998.0) stackConfigDen.limitMinX(xMin);

    if (verbose) {
      cout << "MC denominator systs. (standard):" << endl << endl;
      systCalcDen.printSysts(stackConfigDen["MC"]);
      cout << "MC denominator yields (standard):" << endl << endl;
      systCalcDen.printYields(stackConfigDen);
    }

    // calculate ratio
    systRatCalc.reset(new SystematicRatioCalculator(systCalcNum, systCalcDen, false));
    systRatCalc->calculate(histogram, histogram, false);
    cout.precision(3);
    cout << "Ratio nominal in MC (Standard):" << systRatCalc->_sr["Ratio"]._item[0].nominal << endl;
    cout.precision(1);
    cout << "Systematic uncertainties for MC eff. (standard):" << endl << endl << endl;
    systRatCalc->printAverageSysts(systRatCalc->_sr["Ratio"]);

    systRatCalc->printBinSysts(systRatCalc->_sr["Ratio"]);

    systRatRatCalc.reset(new SystematicRatioRatioCalculator(*systRatCalcD.get(), *systRatCalc.get(), false));
    if (!mcOnly) {
      systRatRatCalc->calculate(histogram, histogram, false);

      std::cout << "Systematics on ratio of ratio: " << std::endl;
      systRatRatCalc->printBinSysts(systRatRatCalc->_sr["Ratio"]);
    }

    vector<string> extraText;
    vector<string> split_extraText;
    split(_extraText, ';', split_extraText);
    for (vector<string>::iterator i = split_extraText.begin(); i!=split_extraText.end();++i) extraText.push_back(*i);
    string outfile = _outfile;
    if (outfile == "") {
      outfile = string("chratio_")+prefix+"_"+histogram;
      if (stamp == 1) outfile += "_ATLASPrelim";
      if (stamp == 2) outfile += "_ATLAS";
      outfile += ".pdf";
    }
    if (!mcOnly) {
      drawChannelRatio(&(systRatCalc->_sr["Ratio"]), extraText, outfile, yTitle, &(systRatCalcD->_sr["Ratio"]), true, mustBeBigger, yMax, xTitle, &(systRatRatCalc->_sr["Ratio"]), lumi);
    } else {
      drawChannelRatio(&(systRatCalc->_sr["Ratio"]), extraText, outfile, yTitle, 0, true, mustBeBigger, yMax, xTitle, 0, lumi);
    }

  } catch (string s) {
    cout << "Crashed with exception: " << s << endl;
    dumpTrace();
  }

  return 0;
}

