/**
 * @brief Analysis base class
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include "TopNtupleAnalysis/AnalysisTrig.h"
#include "TopNtupleAnalysis/Event.h"

AnalysisTrig::AnalysisTrig(const std::string &filename, std::vector<std::string> &systematicList)
  : m_filename(filename), m_hSvc(filename), m_Nduplicate(0) {
  for (size_t i = 0; i < systematicList.size(); ++i) {
    m_hSvc.addSystematics(systematicList[i]);
  }
  m_hSvc.addTrigger("");
}

AnalysisTrig::AnalysisTrig(const std::string &filename)
  : m_filename(filename), m_hSvc(filename), m_Nduplicate(0) {
  m_hSvc.addSystematics("");
  m_hSvc.addTrigger("");
}

AnalysisTrig::~AnalysisTrig() {
}


void AnalysisTrig::clearDuplicateList(){
   m_runEventPair.clear();
}

bool AnalysisTrig::isDuplicateEvent(unsigned int runNumber, unsigned int eventNumber, double leptPt){
  std::pair<unsigned int, unsigned int>  temp( eventNumber, leptPt);
  std::pair<unsigned int, std::pair<unsigned int, unsigned int> > runEvent(runNumber, temp);
  if( m_runEventPair.end() == m_runEventPair.find(runEvent) ) { 
    m_runEventPair.insert(runEvent);
    return false;
  }
  
  return true;
}
