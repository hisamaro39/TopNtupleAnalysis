#ifndef HIST_H
#define HIST_H

#include <string>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include <cmath>
#include <memory>

using namespace std;

/*
 * Class that implements a histogram
 * in a very simple way
 * with auxiliary operators to make the systematics handling easier.
 */
class Hist {
  public:
    vector<double> _x;  // x axis values
    vector<double> _y;  // y axis values
    vector<double> _ye; // y axis error bars
    int _size;          // number of bins
    double _e;          // integral error

    // save when reading TH1
    string xTitle;
    string yTitle;

    /*
     * Constructor.
     * If size == 0, the histogram is treated
     * as a dummy and all operations with it
     * (such as addition or multiplication)
     * react as the identity transformation.
     */
    Hist(int nbins = 0, double value = 0);
    /*
     * Constructor to read the histogram
     * values from a TH1 in a file.
     */
    Hist(const string &name, const string &syst, const string &file);
    /*
     * Copy constructor.
     */
    Hist(const Hist &other);

    /*
     * Acessors.
     */ 
    double &operator[](int i);
    double operator[](int i) const;
    double &x(int i);
    double x(int i) const;
    double &e(int i);
    double e(int i) const;

    /*
     * Copy operation.
     */
    Hist &operator =(const Hist &other);

    /*
     * Operations.
     * Includes addition, multiplication, subtraction and division.
     * Also includes taking the maximum bin-by-bin.
     */
    Hist operator +(Hist a) const;
    Hist operator -(Hist a) const;
    Hist operator *(Hist a) const;
    Hist operator /(Hist a) const;
    Hist &operator +=(Hist a);
    Hist &operator -=(Hist a);
    Hist &operator *=(Hist a);
    Hist &operator /=(Hist a);


    /* 
     * Smoothen it.
     *
     */
    Hist smooth(int smoothLevel = 1);
    Hist smoothRun1(Hist &nominal, int smoothLevel = 1);

    static int GetMaxStatErrBin(Hist &h0, Hist &h1,  double thr, bool skipZero = true);

    /*
     * Multiply by scalar.
     */
    Hist operator *(double a) const;
    Hist &operator *=(double f);

    /*
     * Take maximum bin-by-bin, keeping sign.
     */ 
    void max(Hist a, Hist b);

    /*
     * Return total histogram yield.
     */
    double yield() const;

    /*
     * Return mean in the Y axis.
     */
    double meanY() const;

    /*
     * Return histogram error when reading from TH1.
     */
    double yieldstat() const;

    /*
     * Calculate square root bin-by-bin.
     */
    void squareRoot();

    /*
     * Return histogram to the power of e bin-by-bin.
     */
    Hist power(double e);

    /*
     * Functions that add one extra empty bin in the beginning or end of Hist.
     * This is understood as the underflow or overflow bin, so this allows
     * the existing underflow and overflow bins to be plotted.
     */
    void showUnderflow();
    void showOverflow();

    /*
     * Put entries for x > xMax (x < xMin) in the overflow (underflow) bin and
     * limit the histogram at xMax (xMin) in the X axis.
     */
    void limitMaxX(double xMax, bool addOverflow = false);
    void limitMinX(double xMin);

    /*
     * Rebin plots with factor r.
     */
    void rebin(int r);

    /*
     * Normalise by bin width.
     */
    void normBinWidth(float s = 1.0);

    /*
     * Clip error bars at 1.
     */
    void clipErrorForEfficiency();

    /*
     * Make TH1D from histogram.
     */
    shared_ptr<TH1D> makeTH1(const string name = "myName");

    /*
     *
     * Take absolute value of all bins in histogram.
     */
    Hist abs();

};

ostream &operator<<(ostream &, Hist &h);
ostream &operator<<(ostream &, TH1D &h);

#endif

