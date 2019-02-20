
#include <vector>
#include <string>
#include <map>

#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>

#include <algorithm>

#include "TGraphAsymmErrors.h"
#include <cmath>

#include "TLatex.h"
#include <sstream>

bool includeUF = false;
std::string extraText = "";

TStyle *AtlasStyle()  {
  TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);
  atlasStyle->SetStatW(0.15); 
  atlasStyle->SetStatH(0.11);
  atlasStyle->SetLegendBorderSize(icol);
  //atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);

  // set margin sizes
  atlasStyle->SetPadTopMargin(0.06);
  atlasStyle->SetPadRightMargin(0.08);
  atlasStyle->SetPadBottomMargin(0.15);
  atlasStyle->SetPadLeftMargin(0.16);
  
  // Set Canvas sizes
  atlasStyle->SetCanvasDefH(500);
  atlasStyle->SetCanvasDefW(800);
  atlasStyle->SetCanvasDefH(600);
  atlasStyle->SetCanvasDefW(600);

  // set title offsets (for axis label)
  atlasStyle->SetTitleXOffset(1.2);
  atlasStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  atlasStyle->SetTextFont(font);

  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");
  
  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth(2);
  atlasStyle->SetHatchesLineWidth(1);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //atlasStyle->SetErrorX(0.001);
  // get rid of error bar caps
  atlasStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  atlasStyle->SetOptStat(1111);
  //atlasStyle->SetOptStat(0);
  atlasStyle->SetOptFit(1111);
  //atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);
  
  // PALETTE
  atlasStyle->SetPalette(1,0);

  return atlasStyle;
}

void stampATLAS(const std::string &text, float x, float y) {
  TLatex l;
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(1);
  l.SetTextSize(0.05);
  l.DrawLatex(x, y, "ATLAS");
  TLatex p;
  p.SetNDC();
  p.SetTextFont(42);
  p.SetTextColor(1);
  p.SetTextSize(0.05);
  p.DrawLatex(x+0.12, y, text.c_str());
}
void stampLumi(float lumi, float x, float y) {
  std::stringstream ss;
  ss << "#int L dt = " << lumi << " fb^{-1}";
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextColor(1);
  l.SetTextSize(0.04);
  l.DrawLatex(x, y, ss.str().c_str());
}
void stampText(std::string &text, float x, float y) {
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextColor(1);
  l.SetTextSize(0.05);
  l.DrawLatex(x, y, text.c_str());
}

struct PlotConfig {
  std::string xTitle;
  std::string yTitle;
  float rebin;
};

struct StackElement {
  std::string filen;
  std::string legend;
  std::string legend_style;
  Style_t linestyle;
  Color_t linecolor;
  Style_t fillstyle;
  Color_t fillcolor;
  Style_t markerstyle;
  Int_t markersize;
  std::string option;
  TFile *file;
  int start;
  int stop;

  std::map<std::string, TH1F*> hist;

  StackElement() {
    filen = "";
    legend = "";
    legend_style = "";
    fillcolor = 0;
    fillstyle = 0;
    linecolor = 0;
    linestyle = 0;
    markerstyle = 0;
    markersize = 0;
    option = "";
    file = 0;
    start = 0;
    stop = 0;
  }
  StackElement(const std::string &_file, const std::string &_legend, const std::string &_legend_style,
               Style_t _linestyle, Color_t _linecolor, Style_t _fillstyle, Color_t _fillcolor,
               Style_t _markerstyle, int _markersize, const std::string &_option) {
    filen = _file;
    legend = _legend;
    legend_style = _legend_style;
    fillcolor = _fillcolor;
    fillstyle = _fillstyle;
    linecolor = _linecolor;
    linestyle = _linestyle;
    markerstyle = _markerstyle;
    markersize = _markersize;
    option = _option;
    file = new TFile(_file.c_str());
  }
  StackElement(int _start, int _stop, const std::string &_file, const std::string &_legend, const std::string &_legend_style,
               Style_t _linestyle, Color_t _linecolor, Style_t _fillstyle, Color_t _fillcolor,
               Style_t _markerstyle, int _markersize, const std::string &_option) {
    start = _start;
    stop = _stop;
    filen = _file;
    legend = _legend;
    legend_style = _legend_style;
    fillcolor = _fillcolor;
    fillstyle = _fillstyle;
    linecolor = _linecolor;
    linestyle = _linestyle;
    markerstyle = _markerstyle;
    markersize = _markersize;
    option = _option;
    file = new TFile(_file.c_str());
  }
  TH1F *get(const std::string &hname, const PlotConfig &pc, const std::string &s = "") {
    if (hist.find(s) != hist.end()) {
      return hist[s];
    }
    std::string hs = hname;
    hs += s;
    TH1F *hf = (TH1F *) file->Get(hs.c_str());
    if (pc.rebin != 1 && pc.rebin != 0) hf->Rebin(pc.rebin);
    std::string nname = hf->GetName();
    nname += "_stack_";
    nname += hs;
    TH1F *h = (TH1F *) hf->Clone(nname.c_str());
    if (includeUF && h->GetBinContent(0) != 0) {
      h->SetBinContent(1, h->GetBinContent(0) + h->GetBinContent(1));
      h->SetBinError(1, std::sqrt(std::pow(h->GetBinError(0),2) + std::pow(h->GetBinError(1),2)));
      h->SetBinContent(0, 0);
      h->SetBinError(0, 0);
    }
    //if (h->GetBinContent(h->GetNbinsX()+1) != 0) {
    //  h->SetBinContent(h->GetNbinsX(), h->GetBinContent(h->GetNbinsX()) + h->GetBinContent(h->GetNbinsX()+1));
    //  h->SetBinError(h->GetNbinsX(), std::sqrt(std::pow(h->GetBinError(h->GetNbinsX()),2) + std::pow(h->GetBinError(h->GetNbinsX()+1),2)));
    //  h->SetBinContent(h->GetNbinsX()+1, 0);
    //  h->SetBinError(h->GetNbinsX()+1, 0);
    //}
    h->SetDirectory(0);
    h->SetLineColor(linecolor);
    h->SetFillColor(fillcolor);
    h->SetFillStyle(fillstyle);
    h->SetLineStyle(linestyle);
    h->SetMarkerStyle(markerstyle);
    h->SetMarkerSize(markersize);
    h->SetLineWidth(2);
    //h->SetTitle("");
    //h->GetXaxis()->SetTitle(pc.xTitle.c_str());
    //h->GetYaxis()->SetTitle(pc.yTitle.c_str());
    hist.insert(make_pair(s,h));
    return h;
  }
};

struct StackConfig {
  bool log;

  std::vector<StackElement> main;
  // overlay
  std::vector<StackElement> overlay;
  // standalone plots
  std::vector<StackElement> standalone;

  // for ratios
  int size_ratio;
  std::vector<int> numerator;
  std::vector<int> numerator_type; // 0 = normal stack, 1 = overlay, 2 = standalone
  int denominator;
  int den_type;
  std::string ratio_title;

  // legend
  float leg_x1;
  float leg_y1;
  float leg_x2;
  float leg_y2;
  int   leg_nc;

  StackConfig();

  void setRatio(int _denom, const std::string &_ratio_title);
  void setRatio(int _denom, int _type, const std::string &_ratio_title);
  void addRatio(int _num, int _type);
  
};

void StackConfig::setRatio(int _denom, const std::string &_ratio_title) {
  denominator = _denom;
  ratio_title = _ratio_title;
  den_type = 0;
}
void StackConfig::setRatio(int _denom, int _type, const std::string &_ratio_title) {
  denominator = _denom;
  ratio_title = _ratio_title;
  den_type = _type;
}

void StackConfig::addRatio(int _num, int _type) {
  size_ratio++;
  numerator.push_back(_num);
  numerator_type.push_back(_type);
}

StackConfig::StackConfig() {
  log = false;
  size_ratio = 0;
  denominator = -1;

  leg_x1 = 0.55;
  leg_y1 = 0.76;
  leg_x2 = 0.90;
  leg_y2 = 0.88;
  leg_nc = 2;
}

std::map<std::string, PlotConfig> plotConfig;

void addPlotConfig(const char *dir, const std::string &name, const std::string &xTitle, const std::string &yTitle, float rebin = 1) {
  std::string fname;
  if (dir != "") {
    fname = "/";
    fname += dir;
    fname += "/";
  }
  fname += name;
  PlotConfig p;
  p.xTitle = xTitle;
  p.yTitle = yTitle;
  p.rebin = rebin;
  plotConfig.insert(make_pair(fname, p));
}

PlotConfig &getPlotConfig(const std::string &x) {
  size_t ending = std::string::npos;
  if (ending == std::string::npos) {
    ending = x.find("El");
  }
  if (ending == std::string::npos) {
    ending = x.find("Mu");
  }
  if (ending == std::string::npos) {
    ending = x.find("MJ");
  }
  std::string pref = x.substr(0,ending);
  return plotConfig[pref];
}


std::vector<THStack *> makeStack(const PlotConfig &pc, StackConfig &s, const std::string &hname) {
  std::vector<THStack *> vs;
  THStack *stack = new THStack("stack_master", "");
  for (int i = s.main.size()-1; i >= 0; --i) {
    stack->Add(s.main[i].get(hname, pc), s.main[i].option.c_str());
  }
  vs.push_back(stack);
  for (size_t k = 0; k < s.overlay.size(); ++k) {
    std::string name = "stack_ov";
    name += k;
    THStack *stacko = new THStack(name.c_str(), "");
    for (int i = s.overlay[k].start-1; i >= s.overlay[k].stop; --i) {
      stacko->Add(s.main[i].get(hname, pc), s.main[i].option.c_str());
    }
    stacko->Add(s.overlay[k].get(hname, pc), s.overlay[k].option.c_str());
    vs.push_back(stacko);
  }
  for (size_t k = 0; k < s.standalone.size(); ++k) {
    std::string name = "stack_st";
    name += k;
    THStack *stacks = new THStack(name.c_str(), "");
    stacks->Add(s.standalone[k].get(hname, pc), s.standalone[k].option.c_str());
    vs.push_back(stacks);
  }
  return vs;
}

void zero(TGraphAsymmErrors *syst) {
  for (int i = 0; i < syst->GetN(); ++i) {
    double x, y;
    syst->GetPoint(i, x, y);
    syst->SetPoint(i, x, 0);
  }
}

void diff(TGraphAsymmErrors *syst, TGraphAsymmErrors *nominal) {
  for (int i = 0; i < syst->GetN(); ++i) {
    double x, y;
    syst->GetPoint(i, x, y);
    double xn, yn;
    nominal->GetPoint(i, xn, yn);
    syst->SetPoint(i, x, y - yn);
  }
}

void sum(TGraphAsymmErrors *syst, TGraphAsymmErrors *nominal) {
  for (int i = 0; i < syst->GetN(); ++i) {
    double x, y;
    syst->GetPoint(i, x, y);
    double xn, yn;
    nominal->GetPoint(i, xn, yn);
    syst->SetPoint(i, x, y + yn);
  }
}

std::vector<float> convert(TH1F *h) {
  std::vector<float> x;
  for (int k = 0; k < h->GetNbinsX()+1; ++k) {
    x.push_back(h->GetBinContent(k));
  }
  return x;
}


std::vector<float> zero(TH1F *h) {
  std::vector<float> x;
  for (int k = 0; k < h->GetNbinsX()+1; ++k) {
    x.push_back(0);
  }
  return x;
}

struct SystException {
  SystException() { }
  virtual void load(std::vector<float> &, const std::string &, StackElement &, PlotConfig &, TH1F *) { }
  virtual ~SystException() { }
};

struct SystDummy : public SystException {
  SystDummy() { }
  virtual void load(std::vector<float> &ret, const std::string &/*hname*/, StackElement &, PlotConfig &, TH1F *nominal) {
    ret.clear();
    ret = convert(nominal);
  }
  virtual ~SystDummy() { }
};

struct SystNorm : public SystException {
  double value;
  std::vector<std::string> exclude;
  std::vector<std::string> include;
  SystNorm(double _value = 0.039) : value(_value) { }
  virtual void load(std::vector<float> &ret, const std::string &, StackElement &s, PlotConfig &, TH1F *nominal) {
    ret.clear();
    ret = convert(nominal);
    if ( (std::find(exclude.begin(), exclude.end(), s.legend) == exclude.end()) )
      if ( (include.size() == 0) || (std::find(include.begin(), include.end(), s.legend) != include.end()) )
        for (size_t k = 0; k < ret.size(); ++k) {
          ret[k] *= 1 + value;
          if (ret[k] < 0) ret[k] = 0;
        }
  }
  virtual ~SystNorm() { }
};

struct SystNormAdd : public SystException {
  double value;
  std::vector<std::string> include;
  SystNormAdd(double _value = 7.58*20.2769) : value(_value) { }
  virtual void load(std::vector<float> &ret, const std::string &, StackElement &s, PlotConfig &, TH1F *nominal) {
    ret.clear();
    ret = convert(nominal);
    if (std::find(include.begin(), include.end(), s.legend) != include.end())
      for (size_t k = 0; k < ret.size(); ++k) {
        ret[k] += value;
        if (ret[k] < 0) ret[k] = 0;
      }
  }
  virtual ~SystNormAdd() { }
};

std::map<std::string, SystException *> systExceptions;

void loadExceptions() {
  systExceptions.insert(make_pair("_pdfrw", new SystDummy));
  systExceptions.insert(make_pair("_mcgen", new SystDummy));
  systExceptions.insert(make_pair("_topmass", new SystDummy));
  systExceptions.insert(make_pair("_pshower", new SystDummy));
  systExceptions.insert(make_pair("_fakeNormSystUp", new SystDummy));
  systExceptions.insert(make_pair("_fakeNormSystDown", new SystDummy));
  systExceptions.insert(make_pair("_isrfsr", new SystDummy));
  systExceptions.insert(make_pair("_lumiUp", new SystNorm(0.039)));
  systExceptions.insert(make_pair("_lumiDown", new SystNorm(-0.039)));
  //((SystNorm *) systExceptions["_lumiUp"])->exclude.push_back("W + jets");
  //((SystNorm *) systExceptions["_lumiUp"])->exclude.push_back("QCD multi-jets");
  //((SystNorm *) systExceptions["_lumiDown"])->exclude.push_back("W + jets");
  //((SystNorm *) systExceptions["_lumiDown"])->exclude.push_back("QCD multi-jets");

  systExceptions.insert(make_pair("_ttmxsecUp", new SystNorm(7.58/252.89)));
  systExceptions.insert(make_pair("_ttmxsecDown", new SystNorm(-7.33/252.89)));
  ((SystNorm *) systExceptions["_ttmxsecUp"])->include.push_back("t#bar{t}");
  ((SystNorm *) systExceptions["_ttmxsecDown"])->include.push_back("t#bar{t}");

  systExceptions.insert(make_pair("_ttxsecUp", new SystNorm(13.30/252.89)));
  systExceptions.insert(make_pair("_ttxsecDown", new SystNorm(-14.52/252.89)));
  ((SystNorm *) systExceptions["_ttxsecUp"])->include.push_back("t#bar{t}");
  ((SystNorm *) systExceptions["_ttxsecDown"])->include.push_back("t#bar{t}");

  //systExceptions.insert(make_pair("_zjxsecUp", new SystNorm(0.48)));
  //systExceptions.insert(make_pair("_zjxsecDown", new SystNorm(-0.48)));
  //((SystNorm *) systExceptions["_zjxsecUp"])->include.push_back("Z + jets");
  //((SystNorm *) systExceptions["_zjxsecDown"])->include.push_back("Z + jets");

  //systExceptions.insert(make_pair("_stxsecUp", new SystNorm(0.10)));
  //systExceptions.insert(make_pair("_stxsecDown", new SystNorm(-0.10)));
  //((SystNorm *) systExceptions["_stxsecUp"])->include.push_back("single top");
  //((SystNorm *) systExceptions["_stxsecDown"])->include.push_back("single top");
}

struct SystItem {
  std::string name;
  TH1F *nominal;
  std::map<std::string, std::vector<float> > systs;
  std::vector<float> maxsyst;
  float all_nom;
  float all_syst;

  SystItem() 
    : name(""), nominal(0) {
    all_nom = 0;
    all_syst = 0;
  }
  SystItem(const std::string &_name, TH1F *_nominal, const std::vector<std::string> &_systsName)
    : name(_name), nominal(_nominal) {
    all_nom = 0;
    all_syst = 0;
    for(size_t k = 0; k < _systsName.size(); ++k) {
      systs.insert(make_pair(_systsName[k], std::vector<float>()));
    }
  }
  SystItem(const SystItem &si)
    : name(si.name), nominal(si.nominal) {
    all_nom = 0;
    all_syst = 0;
    for(std::map<std::string, std::vector<float> >::const_iterator it = si.systs.begin(); it != si.systs.end(); ++it) {
      systs.insert(make_pair(it->first, std::vector<float>()));
    }
  }

  void load(const std::string &hname, StackElement &s, PlotConfig &pc) {
    // DEAL WITH EXCEPTIONS!!
    for (std::map<std::string, std::vector<float> >::iterator it = systs.begin(); it != systs.end(); ++it) {
      bool isExcpt = false;
      for (std::map<std::string, SystException *>::iterator jt = systExceptions.begin(); jt != systExceptions.end(); ++jt) {
        if (it->first == jt->first) {
          isExcpt = true;
          jt->second->load(it->second, hname, s, pc, nominal);
          break;
        }
      }
      if (!isExcpt) {
        it->second = convert(s.get(hname, pc, it->first));
      }
    }
  }

  void loadzero(const std::string &hname, StackElement &s, PlotConfig &pc) {
    for (std::map<std::string, std::vector<float> >::iterator it = systs.begin(); it != systs.end(); ++it) {
      it->second = zero(s.get(hname, pc, it->first));
    }
  }

  void calculateMax() {
    all_nom = 0;
    for (int z = 0; z < nominal->GetNbinsX()+1; ++z) all_nom += nominal->GetBinContent(z);
    all_syst = 0;
    std::vector<float> all_syst_m;
    maxsyst.clear();
    for (int i = 0; i < nominal->GetNbinsX()+1; ++i) {
      maxsyst.push_back(0);
    }
    for (std::map<std::string, std::vector<float> >::const_iterator it = systs.begin(); it != systs.end(); ++it) {
      float ss = 0;
      for (size_t z = 0; z < it->second.size(); ++z) ss += it->second[z];
      all_syst_m.push_back(ss);
      for (int i = 0; i < nominal->GetNbinsX()+1; ++i) {
        float x = std::fabs(it->second[i] - nominal->GetBinContent(i));
        if (x > std::fabs(maxsyst[i])) {
          maxsyst[i] = x;
        }
      }
    }
    all_syst = all_nom;
    for (size_t z = 0; z < all_syst_m.size(); ++z) {
      if (std::fabs(all_nom - all_syst_m[z]) > std::fabs(all_nom - all_syst)) all_syst = all_syst_m[z];
    }
  }
};

struct SystConfig {
  TH1F *nominal;
  std::vector<SystItem> systs;
  std::vector<float> total;
  SystConfig(TH1F *_nominal)
    : nominal(0) {
    nominal = _nominal;
  }
  SystConfig(const SystConfig &s)
    : nominal(s.nominal) {
    for (size_t k = 0; k < s.systs.size(); ++k) {
      systs.push_back(s.systs[k]);
      systs[k].nominal = nominal;
    }
  }
  void addSystItem(const std::string &name, const std::vector<std::string> &systsv) {
    systs.push_back(SystItem(name, nominal, systsv));
  }
  void addSystItem(const std::string &name, const char *s[]) {
    std::vector<std::string> systsv;
    for (size_t k = 0; k < sizeof(s)/sizeof(const char *); ++k) {
      systsv.push_back(s[k]);
    }
    addSystItem(name, systsv);
  }
};

std::string bla = "";

void calculateSyst(PlotConfig &pc, StackConfig &s, const std::string &hname, std::vector<SystConfig> &systConfig_per_sample, TGraphAsymmErrors *&error) {
  bla += "a";
  std::string name_nom = "nominalForSyst_";
  name_nom += hname;
  name_nom += bla;
  TH1F *nominal = (TH1F *) s.main[0].get(hname, pc)->Clone(name_nom.c_str());
  for (int i = 0; i < nominal->GetNbinsX()+1; ++i) {
    nominal->SetBinContent(i, 0);
    nominal->SetBinError(i, 0);
  }
  for (size_t i = 0; i < s.main.size(); ++i) {
    nominal->Add(s.main[i].get(hname, pc));
  }
  for (size_t i = 0; i < systConfig_per_sample.size(); ++i) {
    systConfig_per_sample[i].nominal = s.main[i].get(hname, pc);
    for (size_t j = 0; j < systConfig_per_sample[i].systs.size(); ++j) {
      systConfig_per_sample[i].systs[j].nominal = s.main[i].get(hname, pc);
      systConfig_per_sample[i].systs[j].load(hname, s.main[i], pc);
    }
  }

  // adds all nominal inside total_syst
  SystConfig total_syst = systConfig_per_sample[0];
  std::string name_nom1 = total_syst.nominal->GetName();
  name_nom1 += bla;
  total_syst.nominal = (TH1F *) systConfig_per_sample[0].nominal->Clone(name_nom1.c_str());
  for (size_t i = 1; i < systConfig_per_sample.size(); ++i) {
    total_syst.nominal->Add(systConfig_per_sample[i].nominal);
  }
  for (size_t j = 0; j < total_syst.systs.size(); ++j) {
    total_syst.systs[j].nominal = total_syst.nominal;
  }
  // adds all syst. variations individually in the total_syst
  for (size_t i = 0; i < systConfig_per_sample.size(); ++i) {
    for (size_t j = 0; j < systConfig_per_sample[i].systs.size(); ++j) {
      for (std::map<std::string, std::vector<float> >::iterator it = systConfig_per_sample[i].systs[j].systs.begin(); it != systConfig_per_sample[i].systs[j].systs.end(); ++it) {
         std::vector<float> &thisSyst = total_syst.systs[j].systs[it->first];
         if (thisSyst.size() == 0) thisSyst = zero(nominal); 
         std::vector<float> &toAdd = systConfig_per_sample[i].systs[j].systs[it->first];
         for (size_t m = 0; m < thisSyst.size(); ++m) thisSyst[m] += toAdd[m];
      }
    }
  }
  // calculate maximum in groups
  error = new TGraphAsymmErrors(nominal);
  std::vector<float> total(nominal->GetNbinsX()+1, 0);
  total_syst.total.clear();
  total_syst.total.resize(nominal->GetNbinsX()+1, 0);
  for (size_t l = 0; l < total_syst.systs.size(); ++l) {
    total_syst.systs[l].calculateMax();
    for (size_t k = 0; k < total_syst.systs[l].maxsyst.size(); ++k) {
      total[k] += std::pow(total_syst.systs[l].maxsyst[k], 2);
      total_syst.total[k] += std::pow(total_syst.systs[l].maxsyst[k], 2);
    }
  }
  for (size_t k = 0; k < total_syst.total.size(); ++k) {
    total_syst.total[k] = std::sqrt(total_syst.total[k]);
  }
  for (size_t k = 0; k < total.size(); ++k) {
    error->SetPoint(k, nominal->GetBinCenter(k), nominal->GetBinContent(k));
    error->SetPointEXlow(k, nominal->GetBinWidth(k)/2.0);
    error->SetPointEXhigh(k, nominal->GetBinWidth(k)/2.0);
    error->SetPointEYlow(k, std::sqrt(total[k]));
    error->SetPointEYhigh(k, std::sqrt(total[k]));
  }

  for (size_t i = 0; i < systConfig_per_sample.size(); ++i) { // sample
    systConfig_per_sample[i].total.clear();
    systConfig_per_sample[i].total.resize(nominal->GetNbinsX()+1, 0);
    for (size_t l = 0; l < systConfig_per_sample[i].systs.size(); ++l) {
      systConfig_per_sample[i].systs[l].calculateMax();
      for (size_t k = 0; k < systConfig_per_sample[i].systs[l].maxsyst.size(); ++k) {
        systConfig_per_sample[i].total[k] += std::pow(systConfig_per_sample[i].systs[l].maxsyst[k], 2);
      }
    }
    for (size_t k = 0; k < systConfig_per_sample[i].total.size(); ++k) {
      systConfig_per_sample[i].total[k] = std::sqrt(systConfig_per_sample[i].total[k]);
    }
  }
}

void dumpSysts(StackConfig &s, std::vector<SystConfig> &sc, PlotConfig &pc, const std::string &hname) {
  int float_prec = 3;
  float units = 100;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\centering" << std::endl;
  std::cout << "\\caption{Systematic uncertainties in percentages in the muon channel.\
 The maximum between up and down systematic variation in each bin is summed in absolute value to calculate the impact of each uncertainty separately in each sample shown.\
 The total systematic uncertainty in the last column shows the effect of a single systematic variation in the sum of all simulation samples, and it is calculated similarly,\
by summing the absolute value of the variations in each bin.\
 The total systematic uncertainty is the sum, in quadrature, of each element.}" << std::endl;
  std::cout << "\\begin{tabular}{|c|";
  for (size_t i = 0; i < s.main.size(); ++i) {
    std::cout << "c|";
  }
  std::cout << "c|}" << std::endl;

  std::cout << "\\hline" << std::endl;
  std::cout << std::setw(20) << "Syst.";
  float nominal = 0;
  for (size_t i = 0; i < s.main.size(); ++i) {
    std::cout << " & " << std::setw(20) << s.main[i].legend << " ";
    nominal += s.main[i].get(hname, pc)->Integral(0, s.main[i].get(hname, pc)->GetNbinsX()+1);
  }
  std::cout << std::setw(20) << " & Total \\\\" << std::endl;

  std::cout << "\\hline" << std::endl;

  // adds all nominal inside total_syst
  SystConfig total_syst = sc[0];
  total_syst.nominal = (TH1F *) sc[0].nominal->Clone("blablablabla");
  for (size_t i = 1; i < sc.size(); ++i) {
    total_syst.nominal->Add(sc[i].nominal);
  }
  for (size_t j = 0; j < total_syst.systs.size(); ++j) {
    total_syst.systs[j].nominal = total_syst.nominal;
  }
  // adds all syst. variations individually in the total_syst
  for (size_t i = 0; i < sc.size(); ++i) {
    for (size_t j = 0; j < sc[i].systs.size(); ++j) {
      for (std::map<std::string, std::vector<float> >::iterator it = sc[i].systs[j].systs.begin(); it != sc[i].systs[j].systs.end(); ++it) {
         std::vector<float> &thisSyst = total_syst.systs[j].systs[it->first];
         if (thisSyst.size() == 0) thisSyst = zero(total_syst.nominal); 
         std::vector<float> &toAdd = sc[i].systs[j].systs[it->first];
         for (size_t m = 0; m < thisSyst.size(); ++m) thisSyst[m] += toAdd[m];
      }
    }
  }
  total_syst.total.clear();
  total_syst.total.resize(total_syst.nominal->GetNbinsX()+1, 0);
  for (size_t l = 0; l < total_syst.systs.size(); ++l) {
    total_syst.systs[l].calculateMax();
    for (size_t k = 0; k < total_syst.systs[l].maxsyst.size(); ++k) {
      total_syst.total[k] += std::pow(total_syst.systs[l].maxsyst[k], 2);
    }
  }

  std::vector<std::vector<float> > sumit;
  std::vector<float> totsum;
  for (size_t i = 0; i < sc[0].systs.size(); ++i) {
    sumit.push_back(std::vector<float>());
    std::cout << std::setw(20) << sc[0].systs[i].name << " ";
    for (size_t l = 0; l < s.main.size(); ++l) {
      sumit[i].push_back(0);
      for (size_t k = 0; k < sc[l].systs[i].maxsyst.size(); ++k) {
        sumit[i][l] += sc[l].systs[i].maxsyst[k];
      }
      if (nominal != 0)
        std::cout << " & " << std::setw(20) << std::setprecision(float_prec) << units*sumit[i][l]/nominal << " ";
      else
        std::cout << " & " << std::setw(20) << std::setprecision(float_prec) << 0 << " ";
    }
    totsum.push_back(0);
    for (size_t k = 0; k < total_syst.systs[i].maxsyst.size(); ++k) {
      totsum[i] += total_syst.systs[i].maxsyst[k];
    }
    totsum[i] /= nominal;
    std::cout << " & " << std::setw(20) << std::setprecision(float_prec) << units*totsum[i] << "\\\\" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << std::setw(20) << "Total ";
  std::vector<float> sumtot_s;
  for (size_t l = 0; l < s.main.size(); ++l) {
    sumtot_s.push_back(0);
    for (size_t i = 0; i < sc[0].systs.size(); ++i) {
      sumtot_s[l] += std::pow(sumit[i][l]/nominal, 2);
    }
    sumtot_s[l] = std::sqrt(sumtot_s[l]);
    std::cout << " & " << std::setw(20) << std::setprecision(float_prec) << units*sumtot_s[l] << " ";
  }
  float tottot = 0;
  for (size_t i = 0; i < sc[0].systs.size(); ++i) {
    tottot += std::pow(totsum[i], 2);
  }
  tottot = std::sqrt(tottot);
  std::cout << " & " << std::setw(20) << std::setprecision(float_prec) << units*tottot << "\\\\" << std::endl;

  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{table}" << std::endl;
  std::cout << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\centering" << std::endl;
  std::cout << "\\caption{Systematic uncertainties effect in the nominal yields, in percentages in the muon channel.\
 The total systematic uncertainty is the sum, in quadrature, of each element.}" << std::endl;
  std::cout << "\\begin{tabular}{|c|c|}" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << std::setw(20) << "Syst. name " << std::setw(20) << "& Uncertainty \\\\" << std::endl;
  std::cout << "\\hline" << std::endl;
  float s_total_syst = 0;
  for (size_t l = 0; l < total_syst.systs.size(); ++l) {
    std::cout << std::setw(20) << sc[0].systs[l].name << " ";
    if (nominal != 0) {
      std::cout << " & " << std::setw(20) << std::setprecision(float_prec) << units*std::fabs(1 - total_syst.systs[l].all_syst/nominal) << " ";
      s_total_syst += std::pow(1 - total_syst.systs[l].all_syst/nominal, 2);
    } else {
      std::cout << " & " << std::setw(20) << std::setprecision(float_prec) << 0 << " ";
    }
    std::cout << " \\\\" << std::endl;
  }
  s_total_syst = std::sqrt(s_total_syst);
  std::cout << "\\hline" << std::endl;
  std::cout << std::setw(20) << "Total ";
  std::cout << " & " << std::setw(20) << std::setprecision(float_prec) << units*s_total_syst << " \\\\" <<std::endl;;
  std::cout << "\\hline" << std::endl;

  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{table}" << std::endl;
}

void dumpYields(StackConfig &s, std::vector<SystConfig> &sc, PlotConfig &pc, const std::string &hname) {
  int float_prec = 4;

  // calculate systs first
  float nominal = 0;
  for (size_t i = 0; i < s.main.size(); ++i) {
    nominal += s.main[i].get(hname, pc)->Integral(0,  s.main[i].get(hname, pc)->GetNbinsX()+1);
    //nominal += s.main[i].get(hname, pc)->Integral(1,  s.main[i].get(hname, pc)->GetNbinsX());
  }
  std::vector<std::vector<float> > sumit;
  std::vector<float> totsum;
  for (size_t i = 0; i < sc[0].systs.size(); ++i) {
    sumit.push_back(std::vector<float>());
    for (size_t l = 0; l < s.main.size(); ++l) {
      sumit[i].push_back(0);
      sumit[i][l] = std::fabs(sc[l].systs[i].all_syst - sc[l].systs[i].all_nom); // DANILO
      //for (size_t k = 0; k < sc[l].systs[i].maxsyst.size(); ++k) {
      //  sumit[i][l] += sc[l].systs[i].maxsyst[k];
      //}
    }
    totsum.push_back(0);
    for (size_t l = 0; l < s.main.size(); ++l) {
      totsum[i] += sumit[i][l]/nominal;
    }
  }
  std::vector<float> sumtot_s;
  for (size_t l = 0; l < s.main.size(); ++l) {
    sumtot_s.push_back(0);
    for (size_t i = 0; i < sc[0].systs.size(); ++i) {
      sumtot_s[l] += std::pow(sumit[i][l]/nominal, 2);
    }
    sumtot_s[l] = std::sqrt(sumtot_s[l]);
  }
  float tottot = 0;
  for (size_t i = 0; i < sc[0].systs.size(); ++i) {
    tottot += std::pow(totsum[i], 2);
  }
  tottot = std::sqrt(tottot);

  // now dump it
  std::cout << "Dump of yields" << std::endl;

  std::cout << "\\begin{table}" << std::endl;
  std::cout << "\\centering" << std::endl;
  std::cout << "\\caption{Expected event yields in simulation and data.}" << std::endl;
  std::cout << "\\begin{tabular}{|c|ccc|}" << std::endl;

  std::cout << std::setw(20) << "\\hline" << std::endl;

  float tot = 0;
  float tot_errs = 0;
  for (size_t i = 0; i < s.main.size(); ++i) {
    std::cout << std::setw(20) << s.main[i].legend << " & ";
    //float err = sumtot_s[i]*nominal; // should be equivalent
    float err = 0;
    for (size_t k = 0; k < sc[i].total.size(); ++k) {
      err += sc[i].total[k];
    }
    float val = 0;
    double errs = 0;
    val = s.main[i].get(hname, pc)->IntegralAndError(0,s.main[i].get(hname, pc)->GetNbinsX()+1, errs);
    //val = s.main[i].get(hname, pc)->IntegralAndError(1,s.main[i].get(hname, pc)->GetNbinsX(), errs);
    tot += val;
    tot_errs += std::pow(errs, 2);
    std::cout << std::setw(20) << std::setprecision(float_prec) << val  << " & $\\pm$ " << std::setw(20) << std::setprecision(float_prec) << err << " (syst.) & $\\pm$ " << std::setw(20) << std::setprecision(float_prec) << errs << " (stat.) \\\\" << std::endl;
  }
  tot_errs = std::sqrt(tot_errs);
  std::cout << "\\hline" << std::endl;
  std::cout << std::setw(20) << "Expectation & " << std::setw(20) << std::setprecision(float_prec) << tot << " & $\\pm$ " << std::setw(20) << std::setprecision(float_prec) << tottot*nominal << " (syst.) & $\\pm$ " << std::setw(20) << std::setprecision(float_prec) << tot_errs << " (stat.) \\\\" << std::endl;
  std::cout << "\\hline" << std::endl;
  for (size_t i = 0; i < s.standalone.size(); ++i) {
    std::cout << std::setw(20) << s.standalone[i].legend << " & ";
    float val = 0;
    double errs = 0;
    val = s.standalone[i].get(hname, pc)->IntegralAndError(0,s.standalone[i].get(hname, pc)->GetNbinsX()+1, errs);
    //val = s.standalone[i].get(hname, pc)->IntegralAndError(1,s.standalone[i].get(hname, pc)->GetNbinsX(), errs);
    std::cout << std::setw(20) << std::setprecision(float_prec) << val  << " & & $\\pm$ " << std::setw(20) << std::setprecision(float_prec) << errs  << " (stat.) \\\\" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{table}" << std::endl;
  std::cout << std::endl;
}

std::string replaceSlash(std::string s) {
  size_t pos = 0;
  while((pos = s.find("/", pos)) != std::string::npos)
  {
     if (pos == 0) {
       s.replace(pos, 1, "");
     } else {
       s.replace(pos, 1, "_");
       pos += 1;
     }
  }
  return s;
}

void draw(const std::string &hname, StackConfig &s, std::vector<THStack *> &vstack, TCanvas *c, PlotConfig &pc, TGraphAsymmErrors *error = 0, TGraphAsymmErrors *errorO = 0, const std::string &suf = "_mu") {
  // apply style
  TStyle *atlasStyle = AtlasStyle();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
  std::cout << "Using AtlasStyle (" << atlasStyle << ")" << std::endl;

  // do ratio?
  TPad *pad_ratio = 0;
  if (s.size_ratio != 0) {
    c->Divide(1, 2, 0.015, 0.01);
    pad_ratio = (TPad *) c->cd(2);
    pad_ratio->SetPad(0,0,1,0.30);
    pad_ratio->SetTopMargin(0.01);
    pad_ratio->SetBottomMargin(0.45);
    c->cd(1)->SetPad(0,0.2,1,1);
    c->cd(1)->SetBottomMargin(0.13);
  } else {
    c->SetBottomMargin(0.14);
  }

  if (s.log)
    c->cd(1)->SetLogy();

  c->cd(1);
  int nitems = 0;
  nitems += vstack[0]->GetStack()->GetEntries();
  nitems += vstack.size()-1;
  TLegend *leg = new TLegend(s.leg_x1, s.leg_y1 + 0.05 - 0.03*nitems, s.leg_x2, s.leg_y2);
  leg->SetNColumns(s.leg_nc);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  for (int i = (int) s.main.size()-1; i >= 0; --i) {
    leg->AddEntry(s.main[i].get(hname, pc), s.main[i].legend.c_str(), s.main[i].legend_style.c_str());
  }
  for (size_t i = 0; i < s.overlay.size(); ++i) {
    leg->AddEntry(s.overlay[i].get(hname, pc), s.overlay[i].legend.c_str(), s.overlay[i].legend_style.c_str());
  }
  for (size_t i = 0; i < s.standalone.size(); ++i) {
    leg->AddEntry(s.standalone[i].get(hname, pc), s.standalone[i].legend.c_str(), s.standalone[i].legend_style.c_str());
  }

  float maxi = 0;
  if (s.log)
    maxi = 2e-1;
  for (size_t k = 0; k < vstack.size(); ++k) {
    float thisM = 0;
    for (size_t l = 1; l < (size_t) ((TH1 *) vstack[k]->GetStack()->Last())->GetNbinsX(); ++l)
      if (( (TH1 *) vstack[k]->GetStack()->Last())->GetBinContent(l) > thisM)
        thisM = ((TH1 *) vstack[k]->GetStack()->Last())->GetBinContent(l);
    if (thisM > maxi) {
      maxi = thisM;
    }
  }

  // now draw
  for (size_t k = 0; k < vstack.size(); ++k) {
    if (s.log)
      vstack[k]->SetMaximum(10*maxi);
    else
      vstack[k]->SetMaximum(1.7*maxi);
    if (k == 0) vstack[k]->Draw();
    else        vstack[k]->Draw("same");
    vstack[k]->GetXaxis()->SetTitleSize(0.07);
    vstack[k]->GetXaxis()->SetLabelSize(0.07);
    vstack[k]->GetYaxis()->SetTitleSize(0.08);
    vstack[k]->GetYaxis()->SetLabelSize(0.06);
    vstack[k]->GetXaxis()->SetTitleOffset(0.9);
    vstack[k]->GetYaxis()->SetTitleOffset(1.0);
    vstack[k]->GetYaxis()->SetTitleFont(42);
    vstack[k]->GetYaxis()->SetLabelFont(42);
    vstack[k]->GetXaxis()->SetTitleFont(42);
    vstack[k]->GetXaxis()->SetLabelFont(42);
    //vstack[k]->GetYaxis()->SetTitle(pc.yTitle.c_str());
    vstack[k]->GetYaxis()->SetTitle(s.main[0].get(hname, pc)->GetYaxis()->GetTitle());

    if (pad_ratio)
      vstack[k]->GetXaxis()->SetTitle("");
    else
      vstack[k]->GetXaxis()->SetTitle(s.main[0].get(hname, pc)->GetXaxis()->GetTitle());
  }
  c->Update();
  leg->Draw();
  if (errorO) {
    errorO->SetFillColor(kGreen-4);
    errorO->SetFillStyle(3244);
    errorO->SetLineWidth(0);
    errorO->Draw("2 ][ SAME");
  }
  if (error) {
    error->SetFillColor(kBlack);
    error->SetFillStyle(3354);
    error->SetLineWidth(0);
    error->Draw("2 ][ SAME");
  }
  for (size_t k = 0; k < vstack.size(); ++k) {
    vstack[k]->Draw("same");
  }

  if (pad_ratio != 0) {
    c->cd(2);
    TH1F *den = 0;
    if (s.den_type == 0)
      den = (TH1F*) vstack[0]->GetStack()->At(vstack[0]->GetStack()->LastIndex() - s.denominator);
    else
      den = (TH1F*) vstack[1+s.overlay.size()+s.denominator]->GetStack()->At(vstack[1+s.overlay.size()+s.denominator]->GetStack()->LastIndex());

    for (int k = 0; k < s.size_ratio; ++k) {
      TGraphAsymmErrors *er = 0;
      TGraphAsymmErrors *erO = 0;
      TH1F *num = 0;
      if (s.numerator_type[k] == 0) {
        num = (TH1F *) vstack[0]->GetStack()->At(vstack[0]->GetStack()->LastIndex() - s.numerator[k]);
        if (error) er = new TGraphAsymmErrors(*error);
        if (errorO) erO = new TGraphAsymmErrors(*errorO);
      } else if (s.numerator_type[k] == 1) {
        num = (TH1F *) vstack[1+s.numerator[k]]->GetStack()->Last();
        if (error) er = new TGraphAsymmErrors(*error);
        if (errorO) erO = new TGraphAsymmErrors(*errorO);
      } else if (s.numerator_type[k] == 2) {
        num = (TH1F *) vstack[1+s.overlay.size()+s.numerator[k]]->GetStack()->Last();
        if (error) er = new TGraphAsymmErrors(*error);
        if (errorO) erO = new TGraphAsymmErrors(*errorO);
      } else continue;
      std::string namer = num->GetName();
      namer += "ratio";
      namer += k;
      TH1F *rat = (TH1F *) num->Clone(namer.c_str());
      rat->SetDirectory(0);
      rat->Divide(den);
      if (er) {
        for (int m = 0; m < rat->GetNbinsX()+1; ++m) {
          er->SetPoint(m, rat->GetBinCenter(m), 1.0);
          er->SetPointEYlow(m, er->GetErrorYhigh(m)/num->GetBinContent(m));
          er->SetPointEYhigh(m, er->GetErrorYhigh(m)/num->GetBinContent(m));
        }
      }
      if (erO) {
        for (int m = 0; m < rat->GetNbinsX()+1; ++m) {
          erO->SetPoint(m, rat->GetBinCenter(m), 1.0);
          erO->SetPointEYlow(m, erO->GetErrorYhigh(m)/num->GetBinContent(m));
          erO->SetPointEYhigh(m, erO->GetErrorYhigh(m)/num->GetBinContent(m));
        }
      }
      rat->SetStats(0);
      rat->GetYaxis()->SetRangeUser(0.3, 1.7);
      rat->SetTitle("");
      //rat->GetXaxis()->SetTitle(pc.xTitle.c_str());
      rat->GetXaxis()->SetTitle(s.main[0].get(hname, pc)->GetXaxis()->GetTitle());
      rat->GetYaxis()->SetTitle(s.ratio_title.c_str());
      rat->GetYaxis()->SetNdivisions(3, 0, 2);
      rat->GetXaxis()->SetLabelFont(42);
      rat->GetXaxis()->SetTitleFont(42);
      rat->GetYaxis()->SetLabelFont(42);
      rat->GetYaxis()->SetTitleFont(42);
      rat->GetYaxis()->SetLabelSize(0.18);
      rat->GetYaxis()->SetTitleOffset(0.45);
      rat->GetYaxis()->SetLabelOffset(0.02);
      rat->GetYaxis()->SetTitleSize(0.18);
      rat->GetXaxis()->SetTitleSize(0.2);
      rat->GetXaxis()->SetLabelSize(0.18);
      rat->GetXaxis()->SetTitleSize(0.18);
      rat->GetXaxis()->SetTitleOffset(1.1);
      rat->SetTitle("");
      std::string option;
      if (s.numerator_type[k] == 0) {
        option = "e1";//s.main[s.numerator[k]].option;
      } else if (s.numerator_type[k] == 1) {
        option = "e1";//s.overlay[s.numerator[k]].option;
      } else if (s.numerator_type[k] == 2) {
        option = s.standalone[s.numerator[k]].option;
      } else continue;
      if (k != 0)
        option += " same";
      rat->Draw(option.c_str());
      if (erO) {
        erO->SetFillColor(kGreen-4);
        erO->SetFillStyle(3244);
        erO->SetLineWidth(0);
        erO->Draw("2 ][ same");
      }
      if (er) {
        er->SetFillColor(kBlack);
        er->SetFillStyle(3354);
        er->SetLineWidth(0);
        er->Draw("2 ][ same");
      }
      option += " same";
      rat->Draw(option.c_str());
    }
    c->cd(1);
  }
  // major hack to remove the zero label in the main plot, which is cut in half by the ratio pad
  if (!s.log && pad_ratio) {
    TPad *pe = new TPad("pe", "pe", 0, 0, 0.99*c->cd(1)->GetLeftMargin(), 0.18);
    pe->SetFillColor(c->cd(1)->GetFillColor());
    pe->SetBorderMode(0);
    pe->Draw();
  }
  stampATLAS("Internal", 0.15, 0.85);
  stampLumi(1, 0.15, 0.77);
  stampText(extraText, 0.15, 0.70);
  std::string sqrts = "#sqrt{s} = 13 TeV";
  stampText(sqrts, 0.15, 0.65);
  std::string hname1 = replaceSlash(hname);
  hname1 += suf;
  std::string outfileEps = hname1;
  outfileEps += ".eps";
  std::string outfileC = hname1;
  outfileC += ".C";
  c->SaveAs(outfileEps.c_str());
  c->SaveAs(outfileC.c_str());
}

std::vector<SystConfig> buildSystConfig(StackConfig *s, SystConfig &temp) {
  std::vector<SystConfig> x;
  for (size_t k = 0; k < s->main.size(); ++k) {
    x.push_back(SystConfig(temp));
  }
  return x;
}


void initPlotConfig() {
  const char *dirs[] = {""};
  for (size_t k = 0; k < sizeof(dirs)/sizeof(const char *); ++k) {
    addPlotConfig(dirs[k], "lepPt", "Lepton p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lepEta", "Lepton #eta", "Events");
    addPlotConfig(dirs[k], "lepPhi", "Lepton #phi", "Events");
    addPlotConfig(dirs[k], "jetPt", "Small-R jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "jetEta", "Small-R jet #eta", "Events");
    addPlotConfig(dirs[k], "jetPhi", "Small-R jet #phi", "Events");

    addPlotConfig(dirs[k], "largeJetPt", "Large-R jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "largeJetEta", "Large-R jet #eta", "Events");
    addPlotConfig(dirs[k], "largeJetPhi", "Large-R jet #phi", "Events");
    addPlotConfig(dirs[k], "largeJetM", "Large-R jet mass [GeV]", "Events");
    addPlotConfig(dirs[k], "subjetPt", "C/A R=0.2 sub-jets p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "subjetN", "C/A R=0.2 sub-jet multiplicity", "Events");

    addPlotConfig(dirs[k], "subjet0Pt", "C/A R=0.2 leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "subjet1Pt", "C/A R=0.2 sub-leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "subjet2Pt", "C/A R=0.2 third leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "subjet3Pt", "C/A R=0.2 fourth leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "subjet4Pt", "C/A R=0.2 fifth leading sub-jet p_{T} [GeV]", "Events");

    addPlotConfig(dirs[k], "chi", "SD log(#chi)", "Events");
    addPlotConfig(dirs[k], "chi_pt550", "SD log(#chi) for p_{T} > 550 GeV", "Events");
    addPlotConfig(dirs[k], "chi_400pt550", "SD log(#chi) for 400 GeV < p_{T} < 550 GeV", "Events");
    addPlotConfig(dirs[k], "chi_300pt400", "SD log(#chi) for 300 GeV < p_{T} < 400 GeV", "Events");

    addPlotConfig(dirs[k], "sd_sig", "SD log(signal weight)", "Events");
    addPlotConfig(dirs[k], "sd_bkg", "SD log(background weight)", "Events");

    addPlotConfig(dirs[k], "llargeJetPtPos", "Large-R jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "llargeJetPtPosChi40", "Large-R jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "llargeJetPtPosChib40", "Large-R jet p_{T} [GeV]", "Events");

    addPlotConfig(dirs[k], "llargeJetMPos", "Large-R jet mass [GeV]", "Events");
    addPlotConfig(dirs[k], "llargeJetMPosChi40", "Large-R jet mass [GeV]", "Events");
    addPlotConfig(dirs[k], "llargeJetMPosChib40", "Large-R jet mass [GeV]", "Events");

    addPlotConfig(dirs[k], "llargeJetSDMChi40", "SD-weighted large-R jet mass [GeV]", "Events");
    addPlotConfig(dirs[k], "llargeJetSDMChib40", "SD-weighted large-R jet mass [GeV]", "Events");

    addPlotConfig(dirs[k], "llargeJetPt", "Large-R jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "llargeJetEta", "Large-R jet #eta", "Events");
    addPlotConfig(dirs[k], "llargeJetPhi", "Large-R jet #phi", "Events");
    addPlotConfig(dirs[k], "llargeJetM", "Large-R jet mass [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjetPt", "C/A R=0.2 sub-jets p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjetN", "C/A R=0.2 sub-jet multiplicity", "Events");

    addPlotConfig(dirs[k], "lsubjet0Pt", "C/A R=0.2 leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet1Pt", "C/A R=0.2 sub-leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet2Pt", "C/A R=0.2 third leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet3Pt", "C/A R=0.2 fourth leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet4Pt", "C/A R=0.2 fifth leading sub-jet p_{T} [GeV]", "Events");

    addPlotConfig(dirs[k], "lsubjet0Eta", "C/A R=0.2 leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet1Eta", "C/A R=0.2 sub-leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet2Eta", "C/A R=0.2 third leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet3Eta", "C/A R=0.2 fourth leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet4Eta", "C/A R=0.2 fifth leading sub-jet #eta", "Events");

    addPlotConfig(dirs[k], "lsubjet0EtaChi40", "C/A R=0.2 leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet1EtaChi40", "C/A R=0.2 sub-leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet2EtaChi40", "C/A R=0.2 third leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet3EtaChi40", "C/A R=0.2 fourth leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet4EtaChi40", "C/A R=0.2 fifth leading sub-jet #eta", "Events");

    addPlotConfig(dirs[k], "lsubjet0EtaChib40", "C/A R=0.2 leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet1EtaChib40", "C/A R=0.2 sub-leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet2EtaChib40", "C/A R=0.2 third leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet3EtaChib40", "C/A R=0.2 fourth leading sub-jet #eta", "Events");
    addPlotConfig(dirs[k], "lsubjet4EtaChib40", "C/A R=0.2 fifth leading sub-jet #eta", "Events");

    addPlotConfig(dirs[k], "lsubjet0PtChi40", "C/A R=0.2 leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet1PtChi40", "C/A R=0.2 sub-leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet2PtChi40", "C/A R=0.2 third leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet3PtChi40", "C/A R=0.2 fourth leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet4PtChi40", "C/A R=0.2 fifth leading sub-jet p_{T} [GeV]", "Events");

    addPlotConfig(dirs[k], "lsubjet0PtChib40", "C/A R=0.2 leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet1PtChib40", "C/A R=0.2 sub-leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet2PtChib40", "C/A R=0.2 third leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet3PtChib40", "C/A R=0.2 fourth leading sub-jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "lsubjet4PtChib40", "C/A R=0.2 fifth leading sub-jet p_{T} [GeV]", "Events");

    addPlotConfig(dirs[k], "lchi", "SD log(#chi)", "Events");
    addPlotConfig(dirs[k], "lchib", "SD log(#chi) with sub-jet b-tagging", "Events");
    addPlotConfig(dirs[k], "lchi_pt550", "SD log(#chi) for p_{T} > 550 GeV", "Events");
    addPlotConfig(dirs[k], "lchi_400pt550", "SD log(#chi) for 400 GeV < p_{T} < 550 GeV", "Events");
    addPlotConfig(dirs[k], "lchi_300pt400", "SD log(#chi) for 300 GeV < p_{T} < 400 GeV", "Events");

    addPlotConfig(dirs[k], "lsd_sig", "SD log(signal weight)", "Events");
    addPlotConfig(dirs[k], "lsd_bkg", "SD log(background weight)", "Events");

    addPlotConfig(dirs[k], "lsd_sigb", "SD log(signal weight)", "Events");
    addPlotConfig(dirs[k], "lsd_bkgb", "SD log(background weight)", "Events");

    addPlotConfig(dirs[k], "roc_largeJetPt", "Large-R jet p_{T} [GeV]", "Events");
    addPlotConfig(dirs[k], "roc_largeJetEta", "Large-R jet #eta", "Events");
    addPlotConfig(dirs[k], "roc_largeJetPhi", "Large-R jet #phi", "Events");
    addPlotConfig(dirs[k], "roc_largeJetM", "Large-R jet mass [GeV]", "Events");

    addPlotConfig(dirs[k], "roc_sd_sig", "SD log(signal weight)", "Events");
    addPlotConfig(dirs[k], "roc_sd_bkg", "SD log(background weight)", "Events");

    addPlotConfig(dirs[k], "mu", "Average #mu", "Events");
    addPlotConfig(dirs[k], "npv", "Number of Primary Vertices", "Events");

    addPlotConfig(dirs[k], "ldpfjl", "#Delta #phi(large-R jet, lepton)", "Events");
    addPlotConfig(dirs[k], "ldrfjb", "#Delta R(large-R jet, b-jet)", "Events");

    addPlotConfig(dirs[k], "lm12", "Mass of the two leading sub-jets [GeV]", "Events");
    addPlotConfig(dirs[k], "lm23", "Mass of the second and third leading sub-jets [GeV]", "Events");
    addPlotConfig(dirs[k], "lm123", "Mass of the three leading sub-jets [GeV]", "Events");
  }
}

//  StackElement(const std::string &_file, const std::string &_legend, const std::string &_legend_style,
//               Style_t _linestyle, Color_t _linecolor, Style_t _fillstyle, Color_t _fillcolor,
//               Style_t _markerstyle, int _markersize, const std::string &_option);

const std::string concat(const std::string &dir, const std::string &file, const std::string &end = "") {
  std::string ret = dir;
  ret += file;
  ret += end;
  return ret;
}

StackConfig *initStackConfig(const std::string &dir, const std::string &channel = "_mu.root") {
  StackConfig *s = new StackConfig;
  int elements = 0;
  s->main.push_back(StackElement(concat(dir, "13TeV_FS_ZprimePythia3000_e4", channel), "Z' 3 TeV (#times 10^{3})", "F",
                                 1, kBlack, 1001, 5,
                                 1, 0, "hist"));
  elements++;
  s->main.push_back(StackElement(concat(dir, "13TeV_FS_ttbarPowhegPythia_e4", channel), "t#bar{t}", "F",
                                 1, kBlack, 1001, 0,
                                 1, 0, "hist"));
  elements++;
  s->main.push_back(StackElement(concat(dir, "13TeV_FS_SingleTopPowhegPythia_e4", channel), "single top", "F",
                                 1, kBlack, 1001, 62,
                                 1, 0, "hist"));
  elements++;
  s->main.push_back(StackElement(concat(dir, "13TeV_FS_WmassivebcSherpa_e4", channel), "W + jets", "F",
                                 1, kBlack, 1001, 92,
                                 1, 0, "hist"));
  elements++;
  s->main.push_back(StackElement(concat(dir, "13TeV_FS_ZmassivebcSherpa_e4", channel), "Z + jets", "F",
                                 1, kBlack, 1001, 95,
                                 1, 0, "hist"));
  elements++;
  //s->main.push_back(StackElement(concat(dir, "out_VV", channel), "Diboson", "F",
  //                               1, kBlack, 1001, 5,
  //                               1, 0, "hist"));
  //elements++;

  //if (channel == "_mu.root")
  //  s->standalone.push_back(StackElement(0, 0, concat(dir, "out_DataMuon", channel), "Data", "PL",
  //                                       1, kBlack, 1, 0,
  //                                       20, 1, "e1"));
  //else if (channel == "_e.root")
  //  s->standalone.push_back(StackElement(0, 0, concat(dir, "out_DataElectron", channel), "Data", "PL",
  //                                       1, kBlack, 1, 0,
  //                                       20, 1, "e1"));
  //s->setRatio(0, "Data/Sim."); // 0 = index in main stack for denominator, followed by Y axis text
  //s->addRatio(0, 2); // 0 = index in vector for numerator, 2 = standalone
  //s->addRatio(0, 1); // 0 = index in vector for numerator, 1 = overlay
  //s->setRatio(0, 2, "Sim./Data"); // 0 = index in stack for denominator, 2 = use standalone (data) for denominator, followed by Y axis text
  //s->addRatio(0, 0); // 0 = index in vector for numerator, 0 = main stack
  //s->addRatio(0, 1); // 0 = index in vector for numerator, 1 = overlay
  s->log = false;
  return s;
}

SystConfig initSystConfig() {
  loadExceptions();
  SystConfig x(0);
  //const char *lumi[] = {"_lumiUp", "_lumiDown"}; x.addSystItem("Luminosity", lumi);
  //const char *jee[] = {"jee"}; x.addSystItem("JEE", jee);
  //const char *jer[] = {"jer"}; x.addSystItem("JER", jer);
  //const char *jes[] = {"jesUp", "jesDown"}; x.addSystItem("JES", jes);
  //const char *bjes[] = {"boostedJESUp", "_boostedJESDown"}; x.addSystItem("Boosted JES", bjes);
  //const char *bjms[] = {"boostedJMSUp", "_boostedJMSDown"}; x.addSystItem("Boosted JMS", bjms);
  //const char *msmear1[] = {"muSmearMSLOW", "muSmearMSUP"}; x.addSystItem("Muon res. (MS)", msmear1);
  //const char *msmear2[] = {"muSmearIDLOW", "muSmearIDUP"}; x.addSystItem("Muon res. (ID)", msmear2);
  //const char *mscale[] = {"muSmearSCALEUP"}; x.addSystItem("Muon scale", mscale);
  //const char *esmear[] = {"eSmearUp", "eSmearDown"}; x.addSystItem("Electron res.", esmear);
  //const char *ees[] = {"eRescaleUp", "eRescaleDown"}; x.addSystItem("Electron scale", ees);
  //const char *mets[] = {"metResoSoftTermsUp", "metResoSoftTermsDown"}; x.addSystItem("MET Soft Jet scale", mets);
  //const char *metr[] = {"metScaleSoftTermsUp", "metScaleSoftTermsUp"}; x.addSystItem("MET Soft Jet res.", metr);
  //const char *bt1[] = {"_bt1Up", "_bt1Down"}; x.addSystItem("b-tagging efficiency", bt1);
  //const char *bt2[] = {"_bt2Up", "_bt2Down"}; x.addSystItem("c-jet mistag rate", bt2);
  //const char *bt3[] = {"_bt3Up", "_bt3Down"}; x.addSystItem("l-jet mistag rate", bt3);
  //const char *esf[] = {"_esfUp", "_esfDown"}; x.addSystItem("Electron SF", esf);
  //const char *msf[] = {"_musfUp", "_musfDown"}; x.addSystItem("Muon SF", msf);
  //const char *ttx[] = {"_ttxsecUp", "_ttxsecDown"}; x.addSystItem("ttbar cross section", ttx);
  //const char *ttm[] = {"_ttmxsecUp", "_ttmxsecDown"}; x.addSystItem("ttbar mass", ttm);
  //const char *st[] = {"_singletop"}; x.addSystItem("single top cross section", st);
  return x;
}

SystConfig initSystConfigO() {
  SystConfig x(0);
  return x;
}

void plot(const std::string &h = "lepPt", const std::string &dir="./", const std::string &channel = "_mu.root", bool _includeUF = false) {
  if (channel == "_mu.root") {
    extraText = "#mu channel";
  } else if (channel == "_e.root") {
    extraText = "e channel";
  } else if (channel == "_be.root") {
    extraText = "boosted e channel";
  } else if (channel == "_bmu.root") {
    extraText = "boosted #mu channel";
  } else if (channel == "_re.root") {
    extraText = "resolved e channel";
  } else if (channel == "_rmu.root") {
    extraText = "resolved #mu channel";
  }
  includeUF = _includeUF;
  initPlotConfig();
  PlotConfig pc;// = getPlotConfig(h);
  pc.xTitle = "";
  pc.yTitle = "";
  pc.rebin = 1;
  StackConfig *s = initStackConfig(dir, channel);
  SystConfig temp = initSystConfig();
  SystConfig tempO = initSystConfigO();

  std::vector<THStack *> vstack = makeStack(pc, *s, h);
  std::vector<SystConfig> sconfig = buildSystConfig(s, temp);
  std::vector<SystConfig> sconfig_o = buildSystConfig(s, tempO);

  TGraphAsymmErrors *error = 0;
  calculateSyst(pc, *s, h, sconfig, error);

  TGraphAsymmErrors *errorO = 0;
  calculateSyst(pc, *s, h, sconfig_o, errorO);

  dumpSysts(*s, sconfig, pc, h);
  dumpYields(*s, sconfig, pc, h);

  TCanvas *c = new TCanvas("c", "", 800, 600);
  std::string suf = "";
  if (channel == "_mu.root") {
    suf = "_mu";
  } else if (channel == "_e.root") {
    suf = "_e";
  } else if (channel == "_be.root") {
    suf = "_be";
  } else if (channel == "_bmu.root") {
    suf = "_bmu";
  } else if (channel == "_re.root") {
    suf = "_re";
  } else if (channel == "_rmu.root") {
    suf = "_rmu";
  }
  draw(h, *s, vstack, c, pc, error, errorO, suf);
}

