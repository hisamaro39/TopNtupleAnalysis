#ifndef SYSTEMATICIMPLEMENTATION_H
#define SYSTEMATICIMPLEMENTATION_H

#include "Hist.h"
#include <string>
#include <memory>

using namespace std;

/*
 * Base class for all systematic calculation classes.
 * All classes willing to implement a syst. variation,
 * must implement a function get() that returns the Hist with
 * the variations and receives the name of the histogram and the file
 * name of the sample the histogram is in.
 */
class Syst {
  public:
  Syst();
  virtual Hist get(const string &name, const string &fname) = 0;
};

/*
 * Calculates the syst. var. caused by taking the difference between two histograms.
 * The constructor receives the histogram suffixes to be concatenated to histogram name
 * when trying to find the syst. variations in the file.
 */
class HistDiff : public Syst{
  public:
  string _a;
  string _b;
  int _smoothLevel;
  std::vector<std::string> _toExclude;
  HistDiff(const string &a, const string &b, int smoothLevel = 0, const std::vector<std::string> toExclude = std::vector<std::string>());
  Hist get(const string &name, const string &fname);
};

/*
 * Calculates the syst. var. caused by taking the difference between N-1 histograms and a reference histogram.
 * The constructor receives the N histograms suffixes to be concatenated to histogram name
 * when trying to find the syst. variations in the file. The first histogram is the reference.
 */
class HistDiffMany : public Syst{
  public:
  vector<string> _file;
  vector<string> _list;
  vector<string> _sample;
  int _smoothLevel;
  HistDiffMany(vector<string> &_file, vector<string> &list, vector<std::string> &sample, int smoothLevel = 0);
  Hist get(const string &name, const string &fname);
};

/*
 * Calculates the syst. var. caused by a flat normalisation
 * applied to the nominal histogram.
 * The constructor receives the fractional normalisation uncertainty.
 * It can also receive a list of file names for exclusive samples where it should be applied.
 */
class HistNorm : public Syst{
  public:
  double _n;
  std::vector<std::string> _only;
  HistNorm(double n = 0);
  HistNorm(double n, const std::vector<std::string> &only);
  Hist get(const string &name, const string &fname);
};

/*
 * Takes two syst. var. calculation objects
 * and symmetrises the resulting systematic variation,
 * by taking the maximum of the absolute value of the
 * two input variations.
 * The constructor receives the two systematic variation
 * calculation object pointers (to be dynamically allocated).
 * The class takes ownership of the pointers and deletes them.
 */
class Symm : public Syst {
  public:
  unique_ptr<Syst> _a;
  unique_ptr<Syst> _b;
  Symm(Syst *a, Syst *b = 0);
  Hist get(const string &name, const string &fname);
};

/*
 * Systematic variation calculation object that
 * returns a zero histogram in data and it is
 * the identity for everything else.
 * The constructor receives all other systematic uncertainties
 * to be called if the sample is MC.
 */
class NotData : public Syst {
  public:
  unique_ptr<Syst> _a;
  NotData(Syst *a);
  Hist get(const string &name, const string &fname);
};

/*
 * Systematic variation class that adds up two
 * sources of systematic uncertainty in full correlation.
 * The constructor receives the pointers to dynamically
 * allocated systematic uncertainties and the class
 * takes ownership of them, deleting them as necessary.
 */
class SystCorr : public Syst {
  public:
  unique_ptr<Syst> _a;
  unique_ptr<Syst> _b;
  SystCorr(Syst *a, Syst *b);
  Hist get(const string &name, const string &fname);
};

/* Systematic variation class that takes the relative
 * uncertainty betweem the nominal histogram in two different files
 * (whose names are specified in the constructor) and
 * returns them bin-by-bin.
 */
class Relative : public Syst{
  public:
  string _a;
  string _b;
  vector<string> _only;
  double _scale;
  int _smoothLevel;
  Relative(const string &a, const string &b, const vector<string> &only, double scale = 1.0, int smoothLevel = 0);
  Hist get(const string &name, const string &fname);
};

/*
 * Systematic variation class
 * that takes the nominal histogram in two files
 * a and b, given in the constructor, and returns, bin by bin,
 * (a-b)/(a+b).
 * If the vector only is non-empty, it specifies that this is only to be applied
 * to the sample whose file name is given in the contents of "only".
 */
class RelativeISRFSR : public Syst{
  public:
  string _a;
  string _b;
  vector<string> _only;
  int _smoothLevel;
  RelativeISRFSR(const string &a, const string &b, const vector<string> &only, int smoothLevel = 0);
  Hist get(const string &name, const string &fname);
};

#endif

