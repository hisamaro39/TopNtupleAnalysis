//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_uux_ttxuux_H
#define MG5_Sigma_sm_uux_ttxuux_H

//#include <complex> 
//#include <vector> 

//#include "Parameters_sm.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > t t~ u u~ WEIGHTED<=4 @3
// Process: c c~ > t t~ c c~ WEIGHTED<=4 @3
// Process: d d~ > t t~ d d~ WEIGHTED<=4 @3
// Process: s s~ > t t~ s s~ WEIGHTED<=4 @3
//--------------------------------------------------------------------------

class CPPProcess
{
  public:

    // Constructor.
    CPPProcess() {}

    // Initialize process.
    virtual void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual string name() const {return "u u~ > t t~ u u~ (sm)";}

    virtual int code() const {return 3;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 6; 
    static const int nprocesses = 2; 

    // Pointer to the model parameters
    Parameters_sm * pars; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 22; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 14; 
    std::complex<double> amp[namplitudes]; 
    double matrix_3_uux_ttxuux(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_sm_uux_ttxuux_H
