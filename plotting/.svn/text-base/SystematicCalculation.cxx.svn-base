#include "SystematicCalculation.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

SystematicCalculatorBase::SystematicCalculatorBase() { }

void SystematicCalculatorBase::add(const string &name, Syst *s, const string &title) {
  _syst.insert(make_pair<string, unique_ptr<Syst> >(string(name), unique_ptr<Syst>(s)));
  string ntitle = title;
  if (ntitle == "") ntitle = name;
  _syst_title.insert(make_pair<string, string>(string(name), string(ntitle)));
}

double SystematicCalculatorBase::totalSystInSample(Sample &item) {
  double total = 0;
  for (map<string, Hist>::const_iterator k = item.syst.begin(); k != item.syst.end(); ++k) {
    string systName = k->first;
    total += pow(item.syst[systName].yield(), 2);
  }
  return sqrt(total);
}

/*
 * Returns the total systematic unc. of name systName
 * over all samples in the stack st, treated in full correlation
 */
double SystematicCalculatorBase::totalSyst(SampleSet &st, string systName) {
  double total = 0;
  for (int j = 0; j < st._item.size(); ++j) {
    Sample &one_item = st._item[j];
    total += one_item.syst[systName].yield();
  }
  return total;
}

/*
 * Returns the average systematic unc. of name systName
 * over all samples in the stack st, treated in full correlation, bin-by-bin
 */
double SystematicCalculatorBase::averageSyst(SampleSet &st, string systName) {
  Hist total;
  Hist total_syst;
  for (int j = 0; j < st._item.size(); ++j) {
    Sample &one_item = st._item[j];
    total_syst += one_item.syst[systName];
    total += one_item.nominal;
  }
  Hist average = total_syst/total;
  return average.meanY();
}

double SystematicCalculatorBase::totalSyst(SampleSet &st) {
  double total = 0;
  // for each systematic unc.
  for (map<string, Hist>::const_iterator k = st._item[0].syst.begin(); k != st._item[0].syst.end(); ++k) {
    // get the syst. unc in all samples and add in quad.
    total += pow(totalSyst(st, k->first), 2);
  }
  return sqrt(total);
}

double SystematicCalculatorBase::totalYield(SampleSet &st) {
  double total = 0;
  for (int j = 0; j < st._item.size(); ++j) {
    Sample &one_item = st._item[j];
    total += one_item.nominal.yield();
  }
  return total;
}

void SystematicCalculatorBase::printYields(SampleSetConfiguration &sc) {
  cout << "\\begin{tabular}{c|c}" << endl;
  cout << "\\hline" << endl;
  cout << "\\hline" << endl;
  cout << setw(30) << "Type"<< " & " << setw(20) << "Yield" << "\\\\" << endl;
  cout << "\\hline" << endl;
  
  float nData(0.);
  float nDataErr(0.);
  double total = 0;
  double total_stat = 0;
  double total_syst = 0;
  
  for (map<string, SampleSet>::iterator i = sc._stack.begin(); i != sc._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    std::string lcaseName = i->first;
    std::transform(lcaseName.begin(), lcaseName.end(), lcaseName.begin(), ::tolower);
    // i->second is a stack
    // i->second._item is the list of items of the stack
    bool isLast = (next(i) == sc._stack.end());
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      bool lastItem = (isLast && (j == i->second._item.size()-1));
      if (lcaseName.find("data") != std::string::npos) {
        nData=one_item.nominal.yield();
        nDataErr=one_item.nominal.yieldstat();
        break;
      }
      cout << setw(30) << one_item.name_latex << " & " << setw(20)  << one_item.nominal.yield() << " $\\pm$ " << setw(20)  << std::sqrt(std::pow(one_item.nominal.yieldstat(), 2) + std::pow(totalSystInSample(one_item), 2)) << " \\\\ \\hline" << endl;
      total += one_item.nominal.yield();
      total_stat += pow(one_item.nominal.yieldstat(), 2);
    }
    for (map<string, Hist>::iterator k = i->second._item[0].syst.begin(); k != i->second._item[0].syst.end(); ++k) {
      string systName = k->first;
      double s = fabs(totalSyst(i->second, systName));
      total_syst += s*s;
    }
  }
  cout << setw(30) << "Total" << " & " << setw(20)  << total<< " $\\pm$ " << setw(20)  << sqrt(total_stat + total_syst) << " \\\\ \\hline" << endl;
  cout << setw(30) << "Data" << " & " << setw(20)  << nData  << " $\\pm$ " << setw(20) << nDataErr <<  " \\\\ \\hline" << endl;   
  cout << "\\end{tabular}" << endl;
}

void SystematicCalculatorBase::printSysts(SampleSet &st) {
  cout << "\\begin{tabular}{|c|c|}" << endl;
  cout << "\\hline" << endl;
  cout << setw(50) << "Systematics" <<" & " << setw(20) << "Percentual variation" << " \\\\" << endl;
  cout << "\\hline" << endl;

  double total_nom = totalYield(st);
  double total_syst = 0;
  for (map<string, Hist>::iterator k = st._item[0].syst.begin(); k != st._item[0].syst.end(); ++k) {
    string systName = k->first;
    double s = fabs(totalSyst(st, systName));
    string systTitle = _syst_title[systName];
    cout << setw(50) << systTitle << " & " << setw(20) << s*100.0/total_nom;
    cout << " \\\\" << endl;
    total_syst += s*s;
  }
  cout << "\\hline" << endl;
  cout << setw(50) << "Total" << " & " << setw(20) << sqrt(total_syst)*100/total_nom << "\\\\" << endl;
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
}

void SystematicCalculatorBase::printAverageSysts(SampleSet &st) {
  cout << "\\begin{tabular}{|c|c|}" << endl;
  cout << "\\hline" << endl;
  cout << setw(50) << "Systematics" <<" & " << setw(20) << "Percentual variation" << " \\\\" << endl;
  cout << "\\hline" << endl;

  double total_syst = 0;
  for (map<string, Hist>::iterator k = st._item[0].syst.begin(); k != st._item[0].syst.end(); ++k) {
    string systName = k->first;
    double s = fabs(averageSyst(st, systName));
    cout << setw(50) << systName << " & " << setw(20) << s*100.0;
    cout << " \\\\" << endl;
    total_syst += s*s;
  }
  cout << "\\hline" << endl;
  cout << setw(50) << "Total" << " & " << setw(20) << sqrt(total_syst)*100 << "\\\\" << endl;
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
}

void SystematicCalculatorBase::printBinSysts(SampleSet &st) {
  cout << "\\begin{tabular}{|c|";
  for (int z = 1; z < st._item[0].nominal._size-1; ++z) {
    cout << "c|";
  }
  cout << "}" << endl;
  cout << "\\hline" << endl;
  cout << setw(50) << "Systematics" <<" & ";
  for (int z = 1; z < st._item[0].nominal._size-1; ++z) {
    cout << setw(3) << st._item[0].nominal.x(z) << "-" << setw(3) << st._item[0].nominal.x(z+1) << " GeV";
    if (z != st._item[0].nominal._size-1-1) cout << " & ";
  }
  cout << " \\\\" << endl;
  cout << "\\hline" << endl;

  Hist total_syst(st._item[0].nominal);
  for (int z = 0; z < total_syst._size; ++z) total_syst[z] = 0;

  for (map<string, Hist>::iterator k = st._item[0].syst.begin(); k != st._item[0].syst.end(); ++k) {
    string systName = k->first;
    cout << setw(50) << systName << " & ";
    for (int z = 1; z < st._item[0].nominal._size-1; ++z) {
      double s = 0;
      if (st._item[0].nominal[z] != 0)
        s = fabs(st._item[0].syst[systName][z]/st._item[0].nominal[z]);
      cout << setw(15) << s*100.0;
      if (z != st._item[0].nominal._size-1-1) cout << " & ";
      total_syst[z] += pow(st._item[0].syst[systName][z], 2);
    }
    cout << " \\\\" << endl;
  }
  total_syst.squareRoot();
  cout << "\\hline" << endl;
  cout << setw(50) << "Total" << " & ";
  for (int z = 1; z < st._item[0].nominal._size-1; ++z) {
    double s = 0;
    if (st._item[0].nominal[z] != 0)
      s = total_syst[z]/st._item[0].nominal[z];
    cout << setw(15) << s*100.0;
    if (z != st._item[0].nominal._size-1-1) cout << " & ";
  }
  cout << "\\\\" << endl;
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
}

void SystematicCalculatorBase::printBigTable(SampleSetConfiguration &sc) {
  cout << "\\begin{tabular}{|c|";
  for (map<string, SampleSet>::const_iterator i = sc._stack.begin(); i != sc._stack.end(); ++i) cout << "c|";
  cout << "}" << endl;
  cout << "\\hline" << endl;
  cout << setw(50) << "Systematics"<<" & ";
  for (map<string, SampleSet>::const_iterator i = sc._stack.begin(); i != sc._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    bool isLast = (next(i) == sc._stack.end());
    for (int j = 0; j < i->second._item.size(); ++j) {
      const Sample &one_item = i->second._item[j];
      bool lastItem = (isLast && (j == i->second._item.size()-1));
      cout << setw(20) << one_item.name_latex;
      if (!lastItem) cout << " & ";
      else cout << " \\\\" << endl;
    }
  }
  cout << "\\hline" << endl;

  map<string, Hist> &firstListSyst = sc._stack.begin()->second._item[0].syst;
  for (map<string, Hist >::const_iterator k = firstListSyst.begin(); k != firstListSyst.end(); ++k) {
    string systName = k->first;
    cout << setw(50) << systName << " & ";
    for (map<string, SampleSet>::iterator i = sc._stack.begin(); i != sc._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
      // i->second is a stack
      // i->second._item is the list of items of the stack
      bool isLast = (next(i) == sc._stack.end());
      for (int j = 0; j < i->second._item.size(); ++j) {
        Sample &one_item = i->second._item[j];
        bool lastItem = (isLast && (j == i->second._item.size()-1));
        if (one_item.nominal.yield() != 0)
          cout << setw(20)  << 100.0*fabs(one_item.syst[systName].yield())/one_item.nominal.yield();
        else
          cout << setw(20)  << 0;
        if (!lastItem) cout << " & ";
        else cout << " \\\\" << endl;
      }
    }
  }

  cout << "\\hline" << endl;
  cout << setw(50) << "Total" << " & ";
  for (map<string, SampleSet>::iterator i = sc._stack.begin(); i != sc._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    bool isLast = (next(i) == sc._stack.end());
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      bool lastItem = (isLast && (j == i->second._item.size()-1));
      if (one_item.nominal.yield() != 0)
        cout << setw(20)  << 100.0*totalSystInSample(one_item)/one_item.nominal.yield();
      else
        cout << setw(20)  << 0;
      if (!lastItem) cout << " & ";
      else cout << " \\\\" << endl;
    }
  }
  cout << "\\hline" << endl;

  cout << "\\end{tabular}" << endl;
}

SystematicCalculator::SystematicCalculator(SampleSetConfiguration &sc)
  : _sc(sc) {
}

void SystematicCalculator::calculate(const std::string &histogram) {
  for (map<string, SampleSet>::iterator i = _sc._stack.begin(); i != _sc._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      //one_item.nominal and one_item.syst(map of string, Hist)
      one_item.nominal = Hist(histogram, "", one_item.fname);
      for (map<string, unique_ptr<Syst> >::const_iterator k = _syst.begin(); k != _syst.end(); ++k) {
        const string &systName = k->first;
        Syst *thisSyst = k->second.get();
        one_item.syst[systName] = thisSyst->get(histogram, one_item.fname);
      }
    }
  }
}



SystematicRatioCalculator::SystematicRatioCalculator(SystematicCalculator &sca, SystematicCalculator &scb, bool binomial)
  :_sca(sca), _scb(scb), _binomial(binomial) {
  _sr.addType("Ratio");
  _sr.add("Ratio", "", "", "", "");
}

void SystematicRatioCalculator::dataMinusBkg(SampleSetConfiguration &s, Sample &result, bool subtractBkg) {
  Hist nom;
  for (map<string, SampleSet>::iterator i = s._stack.begin(); i != s._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    double s = 1.0;
    if (i->first != "Data" && subtractBkg) s = -1;

    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      //one_item.nominal and one_item.syst(map of string, Hist)
      nom += one_item.nominal*s;
    }
  }

  for (map<string, SampleSet>::iterator i = s._stack.begin(); i != s._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    double s = 1.0;
    if (i->first != "Data" && subtractBkg) s = -1;
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      for (map<string, Hist>::const_iterator k = one_item.syst.begin(); k != one_item.syst.end(); ++k) {
        const string &systName = k->first;
        if (result.syst.find(systName) == result.syst.end())
          result.syst[systName] = nom - nom; // create a zero histogram with the correct binning
        result.syst[systName] += one_item.syst[systName]*s;
      }
    }
  }

  result.nominal = nom;
}

void SystematicRatioCalculator::calculate(const std::string &hnum, const string &hden, bool subtractBkg) {
  // calculate systematic uncertainties in numerator
  // and denominator
  // It should be done beforehand
  //_sca.calculate(hnum);
  //_scb.calculate(hden);

  Sample &ratio = _sr._stack["Ratio"]._item[0];
  Sample numerator;
  Sample denominator;
  dataMinusBkg(_sca._sc, numerator, subtractBkg);
  dataMinusBkg(_scb._sc, denominator, subtractBkg);

  // now one needs to propagate the nominal and syst into the ratio
  vector<string> systNames;
  for (map<string, Hist>::const_iterator k = numerator.syst.begin(); k != numerator.syst.end(); ++k) {
    systNames.push_back(k->first);
  }
  //ratio.nominal and ratio.syst(map of string, Hist)
  if (_binomial) {
    ratio.nominal = numerator.nominal/denominator.nominal;
  } else {
    ratio.nominal = numerator.nominal;
    ratio.nominal /= denominator.nominal;
  }
  //ratio.nominal.clipErrorForEfficiency();
  for (int k = 0; k < systNames.size(); ++k) {
    string name = systNames[k];
    Hist ap = numerator.nominal + numerator.syst[name];
    Hist bp = denominator.nominal + denominator.syst[name];
    ratio.syst[name] = (ap/bp) - ratio.nominal;
  }
}

SystematicRatioRatioCalculator::SystematicRatioRatioCalculator(SystematicRatioCalculator &sca, SystematicRatioCalculator &scb, bool binomial)
  :_sca(sca), _scb(scb), _binomial(binomial) {
  _sr.addType("Ratio");
  _sr.add("Ratio", "", "", "", "");
}

void SystematicRatioRatioCalculator::dataMinusBkg(SampleSetConfiguration &s, Sample &result, bool subtractBkg) {
  Hist nom;
  for (map<string, SampleSet>::iterator i = s._stack.begin(); i != s._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    double s = 1.0;

    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      //one_item.nominal and one_item.syst(map of string, Hist)
      nom += one_item.nominal*s;
    }
  }

  for (map<string, SampleSet>::iterator i = s._stack.begin(); i != s._stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    double s = 1.0;
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      for (map<string, Hist>::const_iterator k = one_item.syst.begin(); k != one_item.syst.end(); ++k) {
        const string &systName = k->first;
        if (result.syst.find(systName) == result.syst.end())
          result.syst[systName] = nom - nom; // create a zero histogram with the correct binning
        result.syst[systName] += one_item.syst[systName]*s;
      }
    }
  }

  result.nominal = nom;
}

void SystematicRatioRatioCalculator::calculate(const std::string &hnum, const string &hden, bool subtractBkg) {
  // calculate systematic uncertainties in numerator
  // and denominator
  // It should be done beforehand

  Sample &ratio = _sr._stack["Ratio"]._item[0];
  Sample numerator;
  Sample denominator;
  dataMinusBkg(_sca._sr, numerator, subtractBkg);
  dataMinusBkg(_scb._sr, denominator, subtractBkg);

  // now one needs to propagate the nominal and syst into the ratio
  vector<string> systNames;
  for (map<string, Hist>::const_iterator k = numerator.syst.begin(); k != numerator.syst.end(); ++k) {
    systNames.push_back(k->first);
  }
  for (map<string, Hist>::const_iterator k = denominator.syst.begin(); k != denominator.syst.end(); ++k) {
    if (std::find(systNames.begin(), systNames.end(), k->first) == systNames.end()) {
      systNames.push_back(k->first);
    }
  }
  //ratio.nominal and ratio.syst(map of string, Hist)
  if (_binomial) {
    ratio.nominal = numerator.nominal/denominator.nominal;
  } else {
    ratio.nominal = numerator.nominal;
    ratio.nominal /= denominator.nominal;
  }
  //ratio.nominal.clipErrorForEfficiency();
  for (int k = 0; k < systNames.size(); ++k) {
    string name = systNames[k];
    Hist ap = numerator.nominal + numerator.syst[name];
    Hist bp = denominator.nominal + denominator.syst[name];
    ratio.syst[name] = (ap/bp) - ratio.nominal;
  }
}


