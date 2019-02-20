/**
 * @brief Analysis base class
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#ifndef ANALYSISTRIG_H
#define ANALYSISTRIG_H

#include <iostream>
#include <string>
#include "TH1F.h"
#include "TFile.h"
#include "TopNtupleAnalysis/Event.h"
#include "TopNtupleAnalysis/HistogramService.h"
#include <set>

class AnalysisTrig {
  public:
    AnalysisTrig(const std::string &filename, std::vector<std::string> &systematicList);
    AnalysisTrig(const std::string &filename);
    virtual ~AnalysisTrig();

    virtual void run(const Event &e, double weight, const std::string &systUnc = "", int is2016run=0, int trthr=0) = 0;
    virtual void terminate() = 0;
    virtual void setIsData(bool isData) = 0;
    virtual void clearDuplicateList();
    virtual unsigned int getNduplicate(){return m_Nduplicate;};
    
  protected:
    std::string m_filename;
    HistogramService m_hSvc;
    //std::set< std::pair<unsigned int, unsigned int> > m_runEventPair;
    std::set< std::pair<unsigned int,  std::pair<unsigned int,unsigned int> > > m_runEventPair;
    bool isDuplicateEvent(unsigned int runNumber, unsigned int eventNumber, double leptPt);
    
    unsigned int m_Nduplicate;
};

#endif

