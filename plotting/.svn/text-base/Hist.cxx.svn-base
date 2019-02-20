#include "Hist.h"
#include <cmath>
#include <iostream>
#include "utils.h"

using namespace std;

Hist::Hist(int nbins, double value) {
  _size = nbins;
  for (int i = 0; i < _size; ++i) {
    _x.push_back(value);
    _y.push_back(value);
    _ye.push_back(value);
  }
  _e = 0;
  yTitle = "Entries";
  xTitle = "x";
}

Hist::Hist(const string &name, const string &syst, const string &file) {
  TFile *f = new TFile(file.c_str());
  if (!f || !f->IsOpen()) {
    throw string("Failed to open file ")+file+string(", trying to get ")+name+string(" for syst ")+syst;
  }
  string n = name;
  n += syst;
  TH1 *h = (TH1 *) f->Get(n.c_str());
  if (!h) {
    throw string("Failed to get histogram in file ")+file+string(", trying to get ")+name+string("for syst ")+syst;
  }
  if (file.find("data") == std::string::npos && file.find("Data") == std::string::npos &&
      file.find("QCD") == std::string::npos && file.find("qcd") == std::string::npos && lumi_scale > 0)
    h->Scale(lumi_scale*1000);
  _size = 0;
  for (int i = 0; i <= h->GetNbinsX()+1; ++i) {
    _size++;
    _x.push_back(h->GetXaxis()->GetBinLowEdge(i));
    _y.push_back(h->GetBinContent(i));
    _ye.push_back(h->GetBinError(i));
  }
  h->IntegralAndError(0,h->GetNbinsX()+1, _e);
  xTitle = h->GetXaxis()->GetTitle();
  yTitle = h->GetYaxis()->GetTitle();
  delete f;
}

Hist::Hist(const Hist &other)
  : _size(other._size), _x(other._x), _y(other._y), _ye(other._ye), _e(other._e),
    xTitle(other.xTitle),yTitle(other.yTitle) {
}

double &Hist::operator[](int i) { return _y.at(i); }
double Hist::operator[](int i) const { return _y.at(i); }
double &Hist::x(int i) { return _x.at(i); }
double Hist::x(int i) const { return _x.at(i); }
double &Hist::e(int i) { return _ye[i]; }
double Hist::e(int i) const { return _ye[i]; }

Hist &Hist::operator =(const Hist &other) {
  _x = other._x;
  _y = other._y;
  _size = other._size;
  _ye = other._ye;
  _e = other._e;
  yTitle = other.yTitle;
  xTitle = other.xTitle;
  return *this;
}

Hist Hist::operator +(Hist a) const {
  if (_size == 0) return a;
  if (a._size == 0) return *this;
  if (a._size != _size) throw string("Trying to add histograms with different sizes.");
  Hist me(*this);
  for (int i = 0; i < _size; ++i) {
    me[i] += a[i];
    me.e(i) = sqrt(me.e(i)*me.e(i) + a.e(i)*a.e(i));
  }
  return me;
}
Hist Hist::operator -(Hist a) const {
  if (_size == 0) return a;
  if (a._size == 0) return *this;
  if (a._size != _size) throw string("Trying to add histograms with different sizes.");
  Hist me(*this);
  for (int i = 0; i < _size; ++i) {
    me[i] -= a[i];
    me.e(i) = sqrt(me.e(i)*me.e(i) + a.e(i)*a.e(i));
  }
  return me;
}

Hist Hist::smooth(int nbins) {
  //if (_size-2 < 3) return *this;
  // need at least 3 bins to smoothen it

  Hist me(*this);
  double pre_yield = me.yield();

  Double_t *xx = new Double_t[_size];
  for (unsigned int i = 0; i < _size; ++i)
    xx[i] = 100 + me[i];

  TH1::SmoothArray(_size, xx, nbins);

  for (int i = 0; i < _size; ++i) {
    me[i] = xx[i] - 100;
  }
  delete [] xx;

  // keep the normalisation fixed before and after smoothing
  //double nb = _size;
  //double perbin = (pre_yield - me.yield())/nb;
  //for (unsigned int i = 0; i < _size; ++i)
  //  me[i] += perbin;
  if (me.yield() != 0)
    me *= pre_yield/me.yield();

  return me;
}

int Hist::GetMaxStatErrBin(Hist &h0, Hist &h1,  double thr, bool skipZero) {
  int iMin = -999;
  double vStatMax = 0;
  for (unsigned int i = 1; i < h0._size-1; ++i) {   // ignore under/overflow bins for now
    double val = h0[i];
    if (skipZero && val == 0)
      continue;

    double err0 = h0.e(i);
    double err1 = h1.e(i);
    double errR = std::fabs(std::sqrt(err0*err0+err1*err1)/val);
    if (errR < thr && errR > vStatMax) {
      vStatMax = errR;
      iMin = i;
    }
  }
  return iMin;
}

/*
 * Start from the bin with lowest stat. 
 * If the syst variation is smaller than stat uncertainty, merge it with the adjacent bin with lower stat.
 * Else check the next lowest bin.
 * Restart after merging is done, until no more merging is needed (either single-bin or all bins with syst larger than stat)
 */
Hist Hist::smoothRun1(Hist &nom, int smoothLevel) {
  if (smoothLevel <= 0) return *this;

  Hist var = *this;

  const float systOverStatThreshold = 2.;
  bool verbose = true;
  verbose = false;
  Hist nom_merged = nom;
  Hist var_merged = var;
  while (nom_merged._size-2 > 1) {    // iterate until single-bin, or nothing to merge
    int nbins = nom_merged._size-1;

    int binToMerge = -999;
    int iMerge = -999;
    double maxErrR = 1e9;
    while (iMerge == -999) {  // Loop over all bins, starting from the one with largest stat error. Stop if merge happens.

      int iMin = Hist::GetMaxStatErrBin(nom_merged, var_merged, maxErrR);
      if (verbose)
        std::cout << "Check " << iMin << "-th bin out of "<< var_merged._size-2 << " bins" << std::endl;
      if (iMin<0)
        break;
      double var = var_merged[iMin];
      double nom = nom_merged[iMin];
      double stat = std::sqrt(std::pow(var_merged.e(iMin), 2) + std::pow(nom_merged.e(iMin), 2));
      maxErrR = stat/std::fabs(nom);
      if (verbose)
        std::cout << Form("yields=%f(%f)\tRelSyst=%f\tRelStat=%f", nom, nom_merged.e(iMin), (var - nom)/nom, maxErrR) << std::endl;

      if (std::fabs(var-nom)/stat < systOverStatThreshold) {
        
        double nomL  = nom_merged[iMin-1];
        double systL = var_merged[iMin-1] - nomL;
        double statL = std::sqrt(std::pow(var_merged.e(iMin-1), 2) + std::pow(nom_merged.e(iMin-1), 2));
        double nomR  = nom_merged[iMin+1];
        double systR = var_merged[iMin+1] - nomR;
        double statR = std::sqrt(std::pow(var_merged.e(iMin+1), 2) + std::pow(nom_merged.e(iMin+1), 2));

        //if (systL*systR<0) continue; // Skip the bins where systematics turning opposite signs (these bins are often merged by following criteria, which remove the shape information)

        bool  mergeL = (iMin != 1 && nomL != 0);
        bool  mergeR = (iMin != nbins && nomR != 0);
        // The neighbor bin is considered not good for merge, if its systematics is statistically significant enough, and can get this bin's systematics over-estimated after merging.
        if (std::fabs(systL)/statL > systOverStatThreshold && std::fabs(systL + var - nom)/(nomL+nom)*nom/stat > systOverStatThreshold) mergeL = false;
        if (std::fabs(systR)/statR > systOverStatThreshold && std::fabs(systR + var - nom)/(nomR+nom)*nom/stat > systOverStatThreshold) mergeR = false;

        if (!mergeL && !mergeR)
          continue;  // Skip bin if both side not appropriate for merge

        iMerge=iMin; // else this bin will be marked as to be merged

        if (mergeL && mergeR) {  // Choose the side with worse statistics to merge, if both sides are appropriate
          double vL = statL/nomL;
          double vR = statR/nomR;
          if (vR > vL) mergeL=false;
          else mergeR=false;
        }

        if (verbose) {
          std::cout << Form("L:%d\tN=%f\tsyst=%f\tstat=%f\tR:%d\tN=%f\tsyst=%f\tstat=%f", mergeL, nomL, systL/nomL, statL/nomL, mergeR, nomR, systR/nomR, statR/nomR) << std::endl;
        }

        
        if (mergeR) binToMerge = iMerge;
        if (mergeL) binToMerge = iMerge-1;
      }
    }
    if (iMerge == -999) break;

    // Rebin and replace 
    Hist nom_new(nom_merged._size-1);
    Hist var_new(nom_merged._size-1);
    int ol = 1;
    nom_new[0] = nom_merged[0];
    nom_new.x(0) = nom_merged.x(0);
    var_new[0] = var_merged[0];
    var_new.x(0) = var_merged.x(0);
    for (int l = 1; l < nom_merged._size; ++l) {
      if (l == binToMerge) {
        nom_new[ol] = nom_merged[l] + nom_merged[l+1];
        nom_new.e(ol) = std::sqrt(std::pow(nom_merged.e(l), 2) + std::pow(nom_merged.e(l+1), 2));
        nom_new.x(ol) = nom_merged.x(l);
        var_new[ol] = var_merged[l] + var_merged[l+1];
        var_new.e(ol) = std::sqrt(std::pow(var_merged.e(l), 2) + std::pow(var_merged.e(l+1), 2));
        var_new.x(ol) = var_merged.x(l);
        l++;
      } else {
        nom_new[ol] = nom_merged[l];
        nom_new.e(ol) = nom_merged.e(l);
        nom_new.x(ol) = nom_merged.x(l);
        var_new[ol] = var_merged[l];
        var_new.e(ol) = var_merged.e(l);
        var_new.x(ol) = var_merged.x(l);
      }
      ol++;
    }
    nom_merged = nom_new;
    var_merged = var_new;

  }
  // Take relative variation from merged hists
  var_merged -= nom_merged;
  var_merged /= nom_merged;

  if (verbose)
    std::cout << Form("Merged %d bins into %d bins.", var._size-2, nom_merged._size-2)<<endl;
  
  for (int i = 1; i < var._size-1; ++i) {
    double vN = nom[i];
    double vX = var.x(i)+0.5*(var.x(i+1)-var.x(i));
    int    iX = 0;
    for (iX = 1; iX < var_merged._size; ++iX) {
      if (var_merged.x(iX) > vX) {
        iX--;
        break;
      }
    }
    if (iX >= var_merged._size) iX = var_merged._size-1;
    double vY = var_merged[iX];
    if (verbose)
      std::cout<< Form("bin %d (%d-th in merged binning): vN=%.3f vOld=%.3f vNew=%.3f relU=%.5f", i, iX, vN, var[i], (1+vY)*vN, vY) << std::endl;
    var[i] = (1+vY)*vN;
  }
  return var;
}


Hist Hist::operator *(Hist a) const {
  if (_size == 0) return a;
  if (a._size == 0) return *this;
  if (a._size != _size) throw string("Trying to add histograms with different sizes.");
  Hist me(*this);
  for (int i = 0; i < _size; ++i) {
    me.e(i) = sqrt(pow(me[i]*a.e(i), 2) + pow(a[i]*me.e(i), 2));
    me[i] *= a[i];
  }
  return me;
}
Hist Hist::operator /(Hist a) const {
  if (_size == 0) return a;
  if (a._size == 0) return *this;
  if (a._size != _size) throw string("Trying to add histograms with different sizes.");
  Hist me(*this);
  for (int i = 0; i < _size; ++i) {
    if (a[i] != 0) {
      me[i] /= a[i];
      // binomial if using a/b
      me.e(i) = sqrt(fabs(((1.-2.*me[i])*me.e(i)*me.e(i) + me[i]*me[i]*a.e(i)*a.e(i) )/(a[i]*a[i]) ));
      //me.e(i) = sqrt(pow((me[i]/pow(a[i], 2))*a.e(i), 2) + pow(me.e(i)/a[i], 2));
    } else if (a[i] == 0 && me[i] == 0) {
      me.e(i) = 1;
      me[i] = 1;
    } else {
      me.e(i) = 0;
      me[i] = 0;
    }
  }
  return me;
}

Hist &Hist::operator +=(Hist a) {
  if (_size == 0) {
    (*this) = a;
    return *this;
  }
  if (a._size == 0) {
    return *this;
  }
  if (a._size != _size) throw string("Trying to add histograms with different sizes.");
  for (int i = 0; i < _size; ++i) {
    (*this).e(i) = sqrt((*this).e(i)*(*this).e(i) + a.e(i)*a.e(i));
    (*this)[i] += a[i];
  }
  return *this;
}


Hist &Hist::operator -=(Hist a) {
  if (_size == 0) {
    (*this) = a;
    return *this;
  }
  if (a._size == 0) {
    return *this;
  }
  if (a._size != _size) throw string("Trying to add histograms with different sizes.");
  for (int i = 0; i < _size; ++i) {
    (*this).e(i) = sqrt((*this).e(i)*(*this).e(i) + a.e(i)*a.e(i));
    (*this)[i] -= a[i];
  }
  return *this;
}

Hist &Hist::operator *=(Hist a) {
  if (_size == 0) {
    (*this) = a;
    return *this;
  }
  if (a._size == 0) {
    return *this;
  }
  
  if (a._size != _size) throw string("Trying to add histograms with different sizes.");
  for (int i = 0; i < _size; ++i) {
    (*this).e(i) = sqrt((*this).e(i)*(*this).e(i) + a.e(i)*a.e(i));
    (*this)[i] *= a[i];
    (*this).e(i) = sqrt((*this).e(i)*(*this).e(i) + a.e(i)*a.e(i));
    (*this)[i] *= a[i];
  }
  return *this;
}

Hist &Hist::operator /=(Hist a) {
  if (_size == 0) {
    (*this) = a;
    return *this;
  }
  if (a._size == 0) {
    return *this;
  }
  
  if (a._size != _size) throw string("Trying to add histograms with different sizes.");
  for (int i = 0; i < _size; ++i) {
    // not binomial if using /=!
    double b1sq = std::pow((*this)[i], 2);
    double b2sq = std::pow(a[i], 2);
    double e1sq = std::pow((*this).e(i), 2);
    double e2sq = std::pow(a.e(i), 2);
    (*this).e(i) = b2sq!=0?std::sqrt((e1sq/b2sq + e2sq*b1sq/(b2sq*b2sq))):0;
    if (a[i] != 0)
      (*this)[i] /= a[i];
    else
      (*this)[i] = 0;
  }
  return *this;
}

Hist Hist::operator *(double a) const {
  Hist me(*this);
  for (int i = 0; i < _size; ++i) {
    me[i] *= a;
    me.e(i) *= a;
  }
  return me;
}
Hist &Hist::operator *=(double f) {
  for (int i = 0; i < _size; ++i) (*this)[i] *= f;
  for (int i = 0; i < _size; ++i) (*this).e(i) *= f;
}


void Hist::max(Hist a, Hist b) {
  *this = a;
  for (int i = 0; i < _size; ++i) {
    if (std::fabs(b[i]) > std::fabs(a[i])) (*this)[i] = std::fabs(b[i]);
    else (*this)[i] = std::fabs(a[i]);
  }
}

double Hist::yield() const {
  double y = 0;
  for (int i = 0; i < _size; ++i) {
    y += (*this)[i];
  }
  return y;
}

double Hist::meanY() const {
  double y = 0;
  for (int i = 1; i < _size-1; ++i) {
    y += (*this)[i]*((*this).x(i+1) - (*this).x(i));
  }
  y /= ((*this).x(_size-1) - (*this).x(1));
  return y;
}

double Hist::yieldstat() const {
  return _e;
}

void Hist::squareRoot() {
  for (int i = 0; i < _size; ++i) {
    (*this).e(i) = sqrt((*this).e(i));
    (*this)[i] = sqrt((*this)[i]);
  }
}

Hist Hist::power(double ex) {
  Hist ret = *this;
  for (int i = 0; i < _size; ++i) {
    ret.e(i) = pow(ret.e(i), ex);
    ret[i] = pow(ret[i], ex);
  }
  return ret;
}

void Hist::showUnderflow() {
  if (_size == 0) return;
  _x.insert(_x.begin(), _x[0]-1);

  _y.insert(_y.begin(), 0);
  _size++;
  _ye.insert(_ye.begin(), 0);
}

void Hist::limitMaxX(double xMax, bool addOverflow) {
  if (_size <= 2) return;
  size_t iMax = 0;
  while (_x[iMax] < xMax) {
    if (iMax >= _x.size()-1) {
      break;
    } else {
      ++iMax;
    }
  }
  if (iMax >= _x.size()-1) return;

  double cumInt = 0;
  double cumIntError = 0;
  xMax = _x[iMax];
  for (size_t i = iMax; i < _x.size(); ++i) {
    cumInt += _y[i];
    cumIntError += pow(_y[i], 2);
    if (!addOverflow) break;
  }
  cumIntError = sqrt(cumIntError);
  _x.erase(_x.begin()+iMax, _x.end());
  _y.erase(_y.begin()+iMax, _y.end());
  _ye.erase(_ye.begin()+iMax, _ye.end());

  _x.push_back(xMax);
  _y.push_back(cumInt);
  _ye.push_back(cumIntError);

  _size = _x.size();
}

void Hist::limitMinX(double xMin) {
  if (_size <= 2) return;
  size_t iMin = _x.size()-1;
  while (_x[iMin] > xMin) {
    if (iMin == 0) {
      break;
    } else {
      --iMin;
    }
  }
  if (iMin <= 0) return;

  double cumInt = 0;
  double cumIntError = 0;
  xMin = _x[iMin];
  for (size_t i = 0; i < iMin; ++i) {
    cumInt += _y[i];
    cumIntError += pow(_y[i], 2);
  }
  cumIntError = sqrt(cumIntError);
  _x.erase(_x.begin(), _x.begin()+iMin);
  _y.erase(_y.begin(), _y.begin()+iMin);
  _ye.erase(_ye.begin(), _ye.begin()+iMin);

  _x.insert(_x.begin(), xMin-1);
  _y.insert(_y.begin(), cumInt);
  _ye.insert(_ye.begin(), cumIntError);

  _size = _x.size();
}

void Hist::normBinWidth(float s) {
  if (_size <= 2) return;

  for (size_t i = 1; i < _size-1; ++i) {
    double binwidth = _x[i+1] - _x[i];
    _y[i] *= s/binwidth;
    _ye[i] *= s/binwidth;
  }
}

void Hist::rebin(int r) {
  if (_size <= 2) return;

  vector<double> nx;
  vector<double> ny;
  vector<double> nye;
  nx.push_back(_x[0]);
  ny.push_back(_y[0]);
  nye.push_back(_ye[0]);
  int count = 0;
  double addX = _x[1];
  double oldX = _x[1];
  double addY = 0;
  double addYE = 0;
  for (size_t i = 1; i < _size-1; ++i) {
    addX = _x[i];
    if (count < r) {
      addY += _y[i];
      addYE += pow(_ye[i], 2);
    } else {
      nx.push_back(oldX);
      ny.push_back(addY);
      nye.push_back(sqrt(addYE));

      oldX  = _x[i];
      addY  = _y[i];
      addYE = pow(_ye[i], 2);
      count = 0;
    }
    count++;
  }
  if (count != 0) {
    nx.push_back(oldX);
    ny.push_back(addY);
    nye.push_back(sqrt(addYE));
  }
  nx.push_back(_x[_size-1]);
  ny.push_back(_y[_size-1]);
  nye.push_back(_ye[_size-1]);
  _x = nx;
  _y = ny;
  _ye = nye;

  _size = _x.size();
}

void Hist::showOverflow() {
  _x.insert(_x.end(), _x[_size]+1);

  _y.insert(_y.end(), 0);
  _size++;
  _ye.insert(_ye.end(), 0);
}

shared_ptr<TH1D> Hist::makeTH1(const string name) {
  shared_ptr<TH1D> h;
  int mySize = _size-2; // skip underflow
  double *myX = new double[mySize+1];
  for (int p = 0; p <= mySize; ++p) {
    myX[p] = (*this).x(p+1);
  }
  h.reset(new TH1D(name.c_str(), name.c_str(), mySize, myX));
  delete [] myX;
  h->Sumw2();
  for (int p = 0; p < _size; ++p) {
    h->SetBinContent(p, (*this)[p]);
    h->SetBinError(p, (*this).e(p));
  }
  h->GetXaxis()->SetTitle(xTitle.c_str());
  h->GetYaxis()->SetTitle(yTitle.c_str());
  return h;
}

Hist Hist::abs() {
  Hist h(*this);
  for (size_t i = 0; i < _size; ++i) h[i] = fabs(h[i]);
  return h;
}

ostream &operator<<(ostream &o, Hist &h) {
  o << "Histogram with x title '" << h.xTitle << "', size = " << h._size << ", yield = " << h.yield() << ": ";
  for (int k = 0; k < h._size; ++k) o << "(" << h.x(k) << ", " << h[k] << ", " << h.e(k) << ") ";
  return o;
}
ostream &operator<<(ostream &o, TH1D &h) {
  o << "TH1D with x title " << h.GetXaxis()->GetTitle() << ", size = " << h.GetNbinsX() << ": ";
  for (int k = 0; k <= h.GetNbinsX()+1; ++k) o << "(" << h.GetXaxis()->GetBinLowEdge(k) << ", " << h.GetBinContent(k) <<") ";
  return o;
}

void Hist::clipErrorForEfficiency() {
  for (int i = 0; i < _size; ++i) {
    if ((*this).e(i) + (*this)[i] > 1) {
      (*this).e(i) = 1 - (*this)[i];
    }
  }
}
