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
    string prefix = "";
    string channel = "mu";
    string h_input = "lepPt";
    string h_other = "";
    string h_titles = "";
    string h_syst = "";
    string h_systTitles = "";
    string _outfile = "";
    int underflow = 0;
    string _extraText = "";
    string yTitle = "";
    string xTitle = "";
    float xMax = -9999;
    float xMin = -9999;
    float lumi = 0;
    std::string config = "";
    int mcOnly = 0;
    smooth = 40;
    int _logY = 0;

    static struct extendedOption extOpt[] = {
        {"help",            no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"channel",         required_argument,     0, 'c', "Channel: one of e, mu, comb.", &channel, extendedOption::eOTString},
        {"prefix",          required_argument,     0, 'p', "Prefix.", &prefix, extendedOption::eOTString},
        {"histogram",       required_argument,     0, 'h', "Histogram.", &h_input, extendedOption::eOTString},
        {"other",           required_argument,     0, 'O', "Other file.", &h_other, extendedOption::eOTString},
        {"titles",          required_argument,     0, 'z', "Other titles.", &h_titles, extendedOption::eOTString},
        {"syst",            required_argument,     0, 's', "Systematics list.", &h_syst, extendedOption::eOTString},
        {"systTitles",      required_argument,     0, 'S', "Systematic list titles.", &h_systTitles, extendedOption::eOTString},
        {"outfile",         required_argument,     0, 'o', "Output file.", &_outfile, extendedOption::eOTString},
        {"underflow",       required_argument,     0, 'U', "Include underflow (0/1)", &underflow, extendedOption::eOTInt},
        {"extraText",       required_argument,     0, 'T', "Extra text to add in the plot.", &_extraText, extendedOption::eOTString},
        {"yTitle",          required_argument,     0, 'y', "Y title.", &yTitle, extendedOption::eOTString},
        {"xTitle",          required_argument,     0, 't', "X title.", &xTitle, extendedOption::eOTString},
        {"xMax",            required_argument,     0, 'X', "Limit X maximum value.", &xMax, extendedOption::eOTFloat},
        {"xMin",            required_argument,     0, 'x', "Limit X minimum value.", &xMin, extendedOption::eOTFloat},
        {"lumi",            required_argument,     0, 'l', "Luminosity value to show.", &lumi, extendedOption::eOTFloat},
        {"config",          required_argument,     0, 'C', "Configuration file.", &config, extendedOption::eOTString},
        {"mcOnly",          required_argument,     0, 'm', "Do not include data? (0/1)", &mcOnly, extendedOption::eOTInt},
        {"smoothen",        required_argument,     0, 'k', "Smoothen systematics.", &smooth, extendedOption::eOTInt},
        {"logY",            required_argument,     0, 'g', "Y axis in log?", &_logY, extendedOption::eOTInt},

        {0, 0, 0, 0, 0, 0, extendedOption::eOTInt}
      };

    
    if (!parseArguments(argc, argv, extOpt) || help) {
      dumpHelp("plotCompareNominal", extOpt, "plot\nCalculate systematic uncertainties and make histograms with them.\n");
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

    std::string histogram = "";
    histogram = h_input;

    vector<string> other_items;
    split(h_other, ',', other_items);

    vector<string> other_titles;
    split(h_titles, ',', other_titles);

    vector<string> syst_items;
    split(h_syst, ',', syst_items);

    vector<string> syst_titles;
    split(h_systTitles, ',', syst_titles);

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(1);

    // for Data/MC comparison
    //SampleSetConfiguration stackConfig = makeConfigurationPlotsCompare(prefix, channel, other_items, other_titles, mcOnly);
    std::vector<std::string> n;
    SampleSetConfiguration stackConfig = makeConfigurationPlotsCompare(prefix, channel, n, n, mcOnly);
    SystematicCalculator systCalc(stackConfig);
    std::vector<std::string> nsyst_items,nsyst_titles;
    for (int z = 0; z < syst_items.size(); ++z) {
      if (syst_items[z].find("pdf_") == 0) {
        if (smooth) {
          nsyst_items.push_back(syst_items[z]+std::string("_smooth")); nsyst_titles.push_back(syst_titles[z]+std::string(" smooth"));
        }
        nsyst_items.push_back(syst_items[z]); nsyst_titles.push_back(syst_titles[z]);

        std::vector<std::string> patterns;
        patterns.push_back("pdf_PDF4LHC15_nlo_30_0");
        patterns.push_back(syst_items[z]);

        std::vector<std::string> filenam;
        std::vector<std::string> sample;
        filenam.push_back(Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), "MC15_13TeV_25ns_FS_EXOT4_ttbaraMcAtNlo_PDF"));
        sample.push_back("ttbar");
        if (smooth) {
          systCalc.add(syst_items[z]+std::string("_smooth"), new NotData(new HistDiffMany(filenam, patterns, sample, true)), syst_items[z]+std::string(" smooth"));
        }
        systCalc.add(syst_items[z], new NotData(new HistDiffMany(filenam, patterns, sample, false)), syst_items[z]);
      } else {
        if (smooth) {
          nsyst_items.push_back(syst_items[z]+std::string("_smooth")); nsyst_titles.push_back(syst_titles[z]+std::string(" smooth"));
        }
        nsyst_items.push_back(syst_items[z]); nsyst_titles.push_back(syst_titles[z]);
        if (smooth) {
          systCalc.add(syst_items[z]+std::string("_smooth"), new NotData(new HistDiff(syst_items[z].c_str(), "", true)), syst_items[z]+std::string(" smooth"));
        }
        systCalc.add(syst_items[z], new NotData(new HistDiff(syst_items[z].c_str(), "", false)), syst_items[z]);
      }
    }
    if (other_items.size() == 2) {
      std::vector<std::string> pat;
      pat.push_back("ttbar");
      if (smooth) {
        std::string s = other_items[0];
        s += std::string("_smooth");
        std::string st = other_titles[0];
        st += std::string(" smooth");
        nsyst_items.push_back(s); nsyst_titles.push_back(st);
        systCalc.add(s, new RelativeISRFSR(Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), other_items[0].c_str()), \
                                           Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), other_items[1].c_str()), \
                                           pat, true), st);
      }
      nsyst_items.push_back(other_items[0]); nsyst_titles.push_back(other_titles[0]);
      systCalc.add(other_items[0], new RelativeISRFSR(Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), other_items[0].c_str()), \
                                                      Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), other_items[1].c_str()), \
                                                      pat, false), other_titles[0]);
      if (smooth) {
        std::string s = other_items[1];
        s += std::string("_smooth");
        std::string st = other_titles[1];
        st += std::string(" smooth");
        nsyst_items.push_back(s); nsyst_titles.push_back(st);
        systCalc.add(s, new RelativeISRFSR(Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), other_items[1].c_str()), \
                                           Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), other_items[0].c_str()), \
                                           pat, true), st);
      }
      nsyst_items.push_back(other_items[1]); nsyst_titles.push_back(other_titles[1]);
      systCalc.add(other_items[1], new RelativeISRFSR(Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), other_items[1].c_str()), \
                                                      Form("%s_%s_%s.root", prefix.c_str(), channel.c_str(), other_items[0].c_str()), \
                                                      pat, false), other_titles[1]);
    }
    //addAllSystematics(systCalc, prefix, channel, false);
    systCalc.calculate(histogram);

    if (underflow) stackConfig.showUnderflow();
    if (xMax > -998.0) stackConfig.limitMaxX(xMax);
    if (xMin > -998.0) stackConfig.limitMinX(xMin);

    cout << "Yields:" << endl << endl;
    systCalc.printYields(stackConfig);
  
    vector<string> extraText;
    string outfile = _outfile;
    if (outfile == "") {
	  if (prefix != "") {
	    outfile = prefix+"_";
	  }
      outfile += histogram;
      outfile += string("_");
      outfile += channel;
      outfile += ".pdf";
    }
    if (channel == "e") {
      extraText.push_back("e channel");
    } else if (channel == "mu") {
      extraText.push_back("#mu channel");
    } else {
      extraText.push_back(channel);
	}
    vector<string> split_extraText;
    split(_extraText, ';', split_extraText);
    for (vector<string>::iterator i = split_extraText.begin(); i!=split_extraText.end();++i) extraText.push_back(*i);
    drawDataMCCompare(stackConfig, extraText, outfile, true, xTitle, nsyst_items, nsyst_titles, lumi);
  } catch (string s) {
    cout << "Crashed with exception: " << s << endl;
    dumpTrace();
  }

  return 0;
}

