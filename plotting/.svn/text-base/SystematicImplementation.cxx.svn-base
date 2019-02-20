#include "SystematicImplementation.h"

#include "Hist.h"
#include <string>
#include <iostream>

using namespace std;

Syst::Syst() {
}

HistDiff::HistDiff(const string &a, const string &b, int smoothLevel, const std::vector<std::string> toExclude) {
  _a = a;
  _b = b;
  _smoothLevel = smoothLevel;
  _toExclude = toExclude;
}

Hist HistDiff::get(const string &name, const string &fname) {
  bool excluded = false;
  for (int k = 0; k < _toExclude.size(); ++k) {
    if (fname.find(_toExclude[k]) != std::string::npos) {
      excluded = true;
      break;
    }
  }
  if (excluded) {
    return Hist();
  }
  Hist ha(name, _a, fname);
  Hist hb(name, _b, fname);
  if (_smoothLevel <= 0)
    return (ha - hb);
  Hist ha_smooth = ha.smoothRun1(hb, _smoothLevel);
  return (ha_smooth - hb);
}

HistDiffMany::HistDiffMany(std::vector<std::string> &file, std::vector<string> &list, std::vector<std::string> &sample, int smoothLevel) {
  _file = file;
  _list = list;
  _sample = sample;
  _smoothLevel = smoothLevel;
}

Hist HistDiffMany::get(const string &name, const string &fname) {
  int idx = -1;
  for (int i = 0; i < _sample.size(); ++i) {
    if (fname.find(_sample[i]) != std::string::npos) {
      idx = i;
      break;
    }
  }
  if (idx < 0)
    return Hist();

  Hist hb(name, _list[0], _file[idx]);
  Hist ret = hb;
  for (size_t i = 0; i < ret._size; ++i) ret[i] = 0;
  for (size_t z = 1; z < _list.size(); ++z) {
    Hist ha(name, _list[z], _file[idx]);
    if (_smoothLevel <= 0) {
      ret += (ha - hb).power(2.0);
    } else {
      Hist ha_smooth = ha.smoothRun1(hb, _smoothLevel);
      ret += (ha_smooth - hb).power(2.0);
    }
  }
  ret.squareRoot();
  return ret;
}

HistNorm::HistNorm(double n) {
  _n = n;
}

HistNorm::HistNorm(double n, const vector<string> &only) {
  _n = n;
  _only = only;
}

Hist HistNorm::get(const string &name, const string &fname) {
  bool good = false;
  for (int i = 0; i < _only.size(); ++i) {
    if (fname.find(_only[i]) != string::npos) good = true;
  }
  if (!good && _only.size() != 0) return Hist();
  Hist ha(name, "", fname);
  ha *= _n;
  return ha;
}

Symm::Symm(Syst *a, Syst *b) {
  _a.reset(a);
  if (b) _b.reset(b);
}

Hist Symm::get(const string &name, const string &fname) {
  Hist ha = _a->get(name, fname);
  if (!_b) {
    return ha;
  }
  Hist hb = _b->get(name, fname);
  Hist hc;
  hc.max(ha, hb);
  return hc;
}

NotData::NotData(Syst *a) {
  _a.reset(a);
}

Hist NotData::get(const string &name, const string &fname) {
  if (fname.find("Data") != string::npos || fname.find("data") != string::npos) return Hist();
  if (fname.find("QCD") != string::npos || fname.find("qcd") != string::npos) return Hist();
  Hist ha = _a->get(name, fname);
  return ha;
}

SystCorr::SystCorr(Syst *a, Syst *b) {
  _a.reset(a);
  _b.reset(b);
}

Hist SystCorr::get(const string &name, const string &fname) {
  Hist ha = _a->get(name, fname);
  Hist hb = _b->get(name, fname);
  Hist hc = ha + hb;
  return hc;
}

Relative::Relative(const string &a, const string &b, const vector<string> &only, double scale, int smoothLevel) {
  _a = a;
  _b = b;
  _only = only;
  _scale = scale;
  _smoothLevel = smoothLevel;
}

Hist Relative::get(const string &name, const string &fname) {
  bool good = false;
  for (int i = 0; i < _only.size(); ++i) {
    if (fname.find(_only[i]) != string::npos) good = true;
  }
  if (!good && _only.size() != 0) return Hist();

  Hist ha(name, "", _a);
  Hist hb(name, "", _b);
  Hist hnom(name, "", fname);
  if (_smoothLevel <= 0)
    return hnom*hb*_scale/ha - hnom;
  Hist ha_smooth = (hnom*hb*_scale/ha).smoothRun1(hnom, _smoothLevel);
  return ha_smooth - hnom;
  //return (ha-hb)*0.5;
}

RelativeISRFSR::RelativeISRFSR(const string &a, const string &b, const vector<string> &only, int smoothLevel) {
  _a = a;
  _b = b;
  _only = only;
  _smoothLevel = smoothLevel;
}

Hist RelativeISRFSR::get(const string &name, const string &fname) {
  bool good = false;
  for (int i = 0; i < _only.size(); ++i) {
    if (fname.find(_only[i]) != string::npos) good = true;
  }
  if (!good) return Hist();

  Hist ha(name, "", _a);
  Hist hb(name, "", _b);
  Hist hnom(name, "", fname);
  Hist hc = ha + hb;
  Hist hd = ha - hb;
  hc *= 0.5;
  if (_smoothLevel <= 0)
    return (hnom*hd/hc); // == (hnom*hd/hc + hnom) - hnom, where the first term can be seen as the variation histogram
  //return (hnom*hd/hc).abs();
  Hist ha_smooth = (hnom*hd/hc + hnom).smoothRun1(hnom, _smoothLevel);
  return ha_smooth - hnom;
}

