#ifndef SYSTEMATICCALCULATION_H
#define SYSTEMATICCALCULATION_H

#include <map>
#include <string>
#include "SystematicImplementation.h"
#include "SampleSet.h"
#include <vector>

#include <iostream>

#include "TGraphErrors.h"

using namespace std;

/*
 * Base class for all others that calculate
 * systematic uncertainties.
 * It has common utilities to print out the systematic uncertainties.
 */
class SystematicCalculatorBase {
  public:
    // maps systematic unc. names to a Syst class that implements
    // it, given a histogram name and a file that identifies the sample to vary
    // Each implementation of this will have its own systematic map
    map<string, unique_ptr<Syst> > _syst;

    map<string, string> _syst_title;

    /* 
     * Constructor. Does nothing
     */
    SystematicCalculatorBase();

    // adds syst. uncertainty
    void add(const string &name, Syst *s, const string &title = "");

    /*
     * Returns sum in quadrature of systematics  in
     * the sample item
     */
    double totalSystInSample(Sample &item);

    /*
     * Returns the total systematic unc. of name systName
     * over all samples in the stack st, treated in full correlation
     */
    double totalSyst(SampleSet &st, string systName);

    /*
     * Returns the total systematic unc. in the samples of st
     * Variations of the same syst. unc. in different samples of st
     * are treated in full correlation.
     * Different systematic variations are treated as independent sources.
     */
    double totalSyst(SampleSet &st);

    /*
     * Returns total yield of samples in st
     */
    double totalYield(SampleSet &st);

    /*
     * Prints yields for each sample of each set in sc
     */
    void printYields(SampleSetConfiguration &sc);

    /*
     * Returns the average systematic for the sample, averaging over different bins.
     */
    double averageSyst(SampleSet &st, string systName);

    /*
     * Prints total systematic uncertainty in st
     * discriminated by sources of syst. unc.
     * but with the contribution of all samples added up.
     */
    void printSysts(SampleSet &st);

    /*
     * Prints bin-by-bin averaged systematic uncertainty in st
     * discriminated by sources of syst. unc.
     * but with the contribution of all samples added up.
     */
    void printAverageSysts(SampleSet &st);

    /*
     * Prints a huge table of all systematic impact in each sample
     */
    void printBigTable(SampleSetConfiguration &sc);


    /*
     * Prints the systematic table per bin.
     */
    void printBinSysts(SampleSet &st);
};

/*
 * Class that implements calculation of the systematic map for a data/MC
 * comparison plot.
 * The resulting systematic uncertainties are stored
 * as histograms (Hist) inside the syst maps of the Sample
 */
class SystematicCalculator : public SystematicCalculatorBase {
  public:
    SampleSetConfiguration &_sc;
    
    SystematicCalculator(SampleSetConfiguration &sc);
    void calculate(const std::string &histogram);
};

/*
 * Implements the uncertainty calculation
 * of a ratio of SystematicCalculators.
 * Aimed at calculating data1 - bkg1 / data2 - bkg2
 * to extract data efficiency.
 */

class SystematicRatioCalculator : public SystematicCalculatorBase {
  public:
    SystematicCalculator &_sca; // numerator configuration
    SystematicCalculator &_scb; // denominator configuration

    SampleSetConfiguration _sr; // ratio setup

    bool _binomial; // use binomial errors?

    SystematicRatioCalculator(SystematicCalculator &sca, SystematicCalculator &scb, bool binomial = true);

    /*
     * Starts from the full data+MC setup in the SampleSetConfiguration s
     * and subtracts the set of samples in the set called "MC" from everything
     * else.
     * Inverts the systematic uncertainties for the set "MC" and adds it to
     * all others.
     * If subtractBkg is false, it only sums all samples in the SampleSetConfiguration and returns it in result
     */
    void dataMinusBkg(SampleSetConfiguration &s, Sample &result, bool subtractBkg = true);

    void calculate(const std::string &hnum, const string &hden, bool subtractBkg = true);
};

class SystematicRatioRatioCalculator : public SystematicCalculatorBase {
  public:
    SystematicRatioCalculator &_sca; // numerator configuration
    SystematicRatioCalculator &_scb; // denominator configuration

    SampleSetConfiguration _sr; // ratio setup

    bool _binomial; // use binomial errors

    SystematicRatioRatioCalculator(SystematicRatioCalculator &sca, SystematicRatioCalculator &scb, bool binomial = true);

    /*
     * Starts from the full data+MC setup in the SampleSetConfiguration s
     * and subtracts the set of samples in the set called "MC" from everything
     * else.
     * Inverts the systematic uncertainties for the set "MC" and adds it to
     * all others.
     * If subtractBkg is false, it only sums all samples in the SampleSetConfiguration and returns it in result
     */
    void dataMinusBkg(SampleSetConfiguration &s, Sample &result, bool subtractBkg = true);

    void calculate(const std::string &hnum, const string &hden, bool subtractBkg = true);
};
#endif

