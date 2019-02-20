#ifndef SAMPLESET_H
#define SAMPLESET_H

#include "Hist.h"
#include <string>
#include <vector>
#include <map>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"

using namespace std;

class Sample {
  public:
    string fname;
    string legname;
    string name_latex;
    string name_plot;

    string legendstyle;
    int linestyle;
    int linecolor;
    int fillstyle;
    int fillcolor;
    int markerstyle;
    float markersize;
    string option;

    Hist nominal;
    map<string, Hist> syst;


    Sample(const string &_fname = "", const string &_legname = "", const string &_latex = "", const string &_plot = "",
           const string &_legendstyle = "F", int _linestyle = 1, int _linecolor = kBlack, int _fillstyle = 1001, int _fillcolor = 0,
           int _markerstyle = 1,
           float _markersize = 0, const string &_option = "hist");

    /*
     * Return TH1D with the proper style.
     */
    shared_ptr<TH1D> makeTH1(const std::string &name = "", const std::string &syst = "");

};

class SampleSet {
  public:
    vector<Sample> _item;

    SampleSet();

    void add(const string &_fname = "", const string &_legname = "", const string &_latex = "", const string &_plot = "",
             const string &_legendstyle = "F", int _linestyle = 1, int _linecolor = kBlack, int _fillstyle = 1001, int _fillcolor = 0,
             int _markerstyle = 1,
             float _markersize = 0, const string &_option = "hist");

    /*
     * Return TGraphErrors with systematic error band of SampleSet st.
     */
    shared_ptr<TGraphAsymmErrors> makeBand(bool isRatio = false);

    /*
     * Return TH1D with systematic error band as histogram for up or down variation of SampleSet st.
     */
    shared_ptr<TH1D> makeBandLine(const string &name = "banderror", bool isRatio = false, bool up = true);

    /*
     * Return THStack with all histograms in the SampleSet st.
     */
    shared_ptr<THStack> makeStack(const string &name, shared_ptr<TLegend> legend, vector<shared_ptr<TH1D> > &vechist);

    /*
     * Return TH1D with the sum of samples.
     * The style will be the style of the first histogram in the sample.
     */
    shared_ptr<TH1D> makeTH1(const string &name = "Data", const string &syst = "");

    /*
     *
     * Save TH1D and variations collected in set of files named hist_[legname].root for each sub sample.
     */
    void saveTH1(const std::string &s);
};

class SampleSetConfiguration {
  public:
    std::map<std::string, SampleSet> _stack;

    void addType(const std::string &type);
    void add(const std::string &type, const std::string &fname, const std::string &legname, const string &_latex, const string &_plot,
             const string &_legendstyle = "LP", int _linestyle = 1, int _linecolor = kBlack, int _fillstyle = 0, int _fillcolor = 0, int _markerstyle = 1,
             float _markersize = 1.0, const string &_option = "e1");
    SampleSet &operator [](const string &name);
    int n();

    /*
     * Move underflow bin to the first visible bin in all histograms.
     */
    void showUnderflow();

    /*
     * Limit maximum/minimum X axis value, moving everything else to the overflow/underflow bin.
     */
    void limitMaxX(double xMax, bool addOverflow = false);
    void limitMinX(double xMin);

    /*
     * Rebin plot with factor r.
     */
    void rebin(int r);

    /*
     * Normalise by bin width.
     */
    void normBinWidth(float s = 1.0);

    /*
     * Assume two SampleSets exist: MC and Data.
     * Normalise all MC histograms (nominal and systs.) to the data yield.
     */
    void normToData();
};

#endif

