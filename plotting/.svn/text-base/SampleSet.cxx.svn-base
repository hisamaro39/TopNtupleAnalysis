#include "SampleSet.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

Sample::Sample(const string &_fname, const string &_legname, const string &_latex, const string &_plot,
               const string &_legendstyle, int _linestyle, int _linecolor, int _fillstyle, int _fillcolor, int _markerstyle,
               float _markersize, const string &_option)
  : fname(_fname), legname(_legname), name_latex(_latex), name_plot(_plot),
    legendstyle(_legendstyle), linestyle(_linestyle), linecolor(_linecolor), fillstyle(_fillstyle), fillcolor(_fillcolor), markerstyle(_markerstyle),
    markersize(_markersize), option(_option) {
  
}

shared_ptr<TH1D> Sample::makeTH1(const string &name, const string &systName) {
  shared_ptr<TH1D> me;
  if (systName == "")
    me = nominal.makeTH1(name);
  else if (syst.find(systName) != syst.end()) {
    Hist thisVar = nominal;
    thisVar += syst[systName];
    me = thisVar.makeTH1(name);
  } else {
    for (map<string, Hist>::const_iterator i = syst.begin(); i != syst.end(); ++i) {
      cout << "NOTE: Systematic uncertainty " << i->first << " is available." << endl;
    }
    throw string("Sample::makeTH1 - Did not find systematic uncertainty ")+systName+(" for sample named ")+legname;
  }

  me->SetLineWidth(2);
  me->SetLineStyle(linestyle);
  me->SetLineColor(linecolor);
  me->SetFillStyle(fillstyle);
  me->SetFillColor(fillcolor);
  me->SetMarkerStyle(markerstyle);
  me->SetMarkerSize(markersize);
  me->SetMarkerColor(linecolor);
  return me;
}

SampleSet::SampleSet() {
}

void SampleSet::add(const string &_fname, const string &_legname, const string &_latex, const string &_plot,
                    const string &_legendstyle, int _linestyle, int _linecolor, int _fillstyle, int _fillcolor, int _markerstyle,
                    float _markersize, const string &_option) {
  _item.push_back(Sample(_fname, _legname, _latex, _plot, _legendstyle, _linestyle, _linecolor, _fillstyle, _fillcolor, _markerstyle, _markersize, _option));
}

shared_ptr<TGraphAsymmErrors> SampleSet::makeBand(bool isRatio) {
  shared_ptr<TGraphAsymmErrors> g(new TGraphAsymmErrors());
  g->Set((int) (_item[0].nominal._size - 2));
  Hist total;
  Hist syst;
  vector<string> sNames;
  for (map<string, Hist>::iterator k = _item[0].syst.begin(); k != _item[0].syst.end(); ++k) {
    string systName = k->first;
    sNames.push_back(systName);
  }

  for (int j = 0; j < _item.size(); ++j) {
    Sample &one_item = _item[j];
    total += one_item.nominal;
  }

  for (int k = 0; k < sNames.size(); ++k) {
    Hist this_syst;
    for (int j = 0; j < _item.size(); ++j) {
      Sample &one_item = _item[j];
      this_syst += one_item.syst[sNames[k]];
    }
    syst += this_syst*this_syst;
  }
  syst.squareRoot();
  for (int k = 1; k < syst._size-1; ++k) {
    double delta = (total.x(k+1)-total.x(k))*0.5;
    g->SetPoint(k-1, total.x(k)+delta, total[k]);
    g->SetPointEXhigh(k-1, delta);
    g->SetPointEXlow(k-1, delta);
    double eyh = syst[k];
    double eyl = syst[k];
    if (total[k] - eyl <= 0) eyl = total[k]*0.9999;
    if (isRatio && total[k] + eyh > 1) eyh = 1 - total[k];
    g->SetPointEYhigh(k-1, eyh);
    g->SetPointEYlow(k-1, eyl);
  }
  return g;
}

shared_ptr<TH1D> SampleSet::makeBandLine(const string &name, bool isRatio, bool up) {
  shared_ptr<TH1D> g(_item[0].makeTH1(Form("%s%d", name.c_str(), (int) up)));
  Hist total;
  Hist syst;
  vector<string> sNames;
  for (map<string, Hist>::iterator k = _item[0].syst.begin(); k != _item[0].syst.end(); ++k) {
    string systName = k->first;
    sNames.push_back(systName);
  }

  for (int j = 0; j < _item.size(); ++j) {
    Sample &one_item = _item[j];
    total += one_item.nominal;
  }

  for (int k = 0; k < sNames.size(); ++k) {
    Hist this_syst;
    for (int j = 0; j < _item.size(); ++j) {
      Sample &one_item = _item[j];
      this_syst += one_item.syst[sNames[k]];
    }
    syst += this_syst*this_syst;
  }
  syst.squareRoot();
  for (int k = 1; k < syst._size-1; ++k) {
    double factor = 1;
    if (!up) factor = -1;
    double ll = total[k] + factor*syst[k];
    if (!up && ll <= 0) ll = 0.0001;
    g->SetBinContent(k, total[k] + factor*syst[k]);
  }
  return g;
}

shared_ptr<THStack> SampleSet::makeStack(const string &name, shared_ptr<TLegend> legend, vector<shared_ptr<TH1D> > &vechist) {
  shared_ptr<THStack> stack(new THStack(name.c_str(), name.c_str()));
  vector<string> names;
  vector<string> style;
  for (int j = _item.size()-1; j >= 0; --j) { // this inverts the stack
    Sample &one_item = _item[j];
    vechist.push_back(one_item.makeTH1(Form("%s%d", name.c_str(), j)));
    stack->Add(vechist[vechist.size()-1].get(), one_item.option.c_str());
    names.push_back(one_item.name_plot);
    style.push_back(one_item.legendstyle);
  }
  for (int j = vechist.size()-1; j >= 0; --j) {
    if (style[j] != "")
      legend->AddEntry(vechist[j].get(), names[j].c_str(), style[j].c_str());
  }
  return stack;
}

shared_ptr<TH1D> SampleSet::makeTH1(const string &name, const string &systName) {
  shared_ptr<TH1D> mySum(_item[0].makeTH1(name, systName));
  for (int j = 1; j < _item.size(); ++j) {
    shared_ptr<TH1D> t(_item[j].makeTH1(Form("%s_%d", name.c_str(), j), systName));
    mySum->Add(t.get());
  }
  return shared_ptr<TH1D>(mySum);
}

void SampleSet::saveTH1(const std::string &s) {
  for (int j = 0; j < _item.size(); ++j) {
    TFile f((string("hist_")+_item[j].legname+string(".root")).c_str(), "update");
    TH1D *t = new TH1D(*_item[j].makeTH1(Form("x%s%d", s.c_str(), j), "").get());
    f.cd();
    t->Write(Form("x%s", s.c_str()), TObject::kOverwrite);
    delete t;
    for (map<string, Hist>::const_iterator it = _item[j].syst.begin(); it != _item[j].syst.end(); ++it) {
      TH1D *ts = new TH1D(*_item[j].makeTH1(Form("x%s%s%d", s.c_str(), it->first.c_str(), j), it->first).get());
      f.cd();
      ts->Write(Form("x%s%s", s.c_str(), it->first.c_str()), TObject::kOverwrite);
      delete ts;
    }
    f.Close();
  }
}


void SampleSetConfiguration::addType(const std::string &type) {
  _stack.insert(std::make_pair(type, SampleSet()));
}

void SampleSetConfiguration::add(const std::string &type, const std::string &fname, const std::string &legname, const string &_latex, const string &_plot,
                                 const string &_legendstyle, int _linestyle, int _linecolor, int _fillstyle, int _fillcolor, int _markerstyle,
                                 float _markersize, const string &_option) {

  _stack[type].add(fname, legname, _latex, _plot, _legendstyle, _linestyle, _linecolor, _fillstyle, _fillcolor, _markerstyle, _markersize, _option);
}

SampleSet &SampleSetConfiguration::operator [](const string &name) {
  return _stack[name];
}


void SampleSetConfiguration::showUnderflow() {
  for (map<string, SampleSet>::iterator i = _stack.begin(); i != _stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      //one_item.nominal and one_item.syst(map of string, Hist)
      one_item.nominal.showUnderflow();
      for (map<string, Hist>::const_iterator k = one_item.syst.begin(); k != one_item.syst.end(); ++k) {
        const string &systName = k->first;
        one_item.syst[systName].showUnderflow();
      }
    }
  }
}

void SampleSetConfiguration::limitMaxX(double xMax, bool addOverflow) {
  for (map<string, SampleSet>::iterator i = _stack.begin(); i != _stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      //one_item.nominal and one_item.syst(map of string, Hist)
      one_item.nominal.limitMaxX(xMax, addOverflow);
      for (map<string, Hist>::const_iterator k = one_item.syst.begin(); k != one_item.syst.end(); ++k) {
        const string &systName = k->first;
        one_item.syst[systName].limitMaxX(xMax, addOverflow);
      }
    }
  }
}

void SampleSetConfiguration::rebin(int r) {
  for (map<string, SampleSet>::iterator i = _stack.begin(); i != _stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      //one_item.nominal and one_item.syst(map of string, Hist)
      one_item.nominal.rebin(r);
      for (map<string, Hist>::const_iterator k = one_item.syst.begin(); k != one_item.syst.end(); ++k) {
        const string &systName = k->first;
        one_item.syst[systName].rebin(r);
      }
    }
  }
}

void SampleSetConfiguration::normBinWidth(float s) {
  for (map<string, SampleSet>::iterator i = _stack.begin(); i != _stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      //one_item.nominal and one_item.syst(map of string, Hist)
      one_item.nominal.normBinWidth(s);
      for (map<string, Hist>::const_iterator k = one_item.syst.begin(); k != one_item.syst.end(); ++k) {
        const string &systName = k->first;
        one_item.syst[systName].normBinWidth(s);
      }
    }
  }
}

void SampleSetConfiguration::limitMinX(double xMin) {
  for (map<string, SampleSet>::iterator i = _stack.begin(); i != _stack.end(); ++i) { // for every element in the stack (one for data, one for MC)
    // i->second is a stack
    // i->second._item is the list of items of the stack
    for (int j = 0; j < i->second._item.size(); ++j) {
      Sample &one_item = i->second._item[j];
      //one_item.nominal and one_item.syst(map of string, Hist)
      one_item.nominal.limitMinX(xMin);
      for (map<string, Hist>::const_iterator k = one_item.syst.begin(); k != one_item.syst.end(); ++k) {
        const string &systName = k->first;
        one_item.syst[systName].limitMinX(xMin);
      }
    }
  }
}

void SampleSetConfiguration::normToData() {
  Hist &myData = _stack["Data"]._item[0].nominal;
  double data_yield = myData.yield();

  Sample &firstMC = _stack["MC"]._item[0];

  Hist fullMC = _stack["MC"]._item[0].nominal;
  for (int j = 1; j < _stack["MC"]._item.size(); ++j) {
    fullMC += _stack["MC"]._item[j].nominal;
  }

  for (map<string, Hist>::const_iterator k = firstMC.syst.begin(); k != firstMC.syst.end(); ++k) {
    const string &systName = k->first;

    Hist fullMCSyst = _stack["MC"]._item[0].syst[systName];
    for (int j = 1; j < _stack["MC"]._item.size(); ++j) {
      fullMCSyst += _stack["MC"]._item[j].syst[systName];
    }

    for (int j = 0; j < _stack["MC"]._item.size(); ++j) {
      Sample &myMC = _stack["MC"]._item[j];
      Hist a = myMC.nominal + myMC.syst[systName];
      a *= 1.0/(fullMC + fullMCSyst).yield();
      Hist b = myMC.nominal;
      b *= 1.0/fullMC.yield();
      a -= b;
      myMC.syst[systName] = a*data_yield;
    }
  }
  for (int j = 0; j < _stack["MC"]._item.size(); ++j) {
    Sample &myMC = _stack["MC"]._item[j];
    myMC.nominal *= data_yield/fullMC.yield();
  }
}


int SampleSetConfiguration::n() {
  return _stack.size();
}


