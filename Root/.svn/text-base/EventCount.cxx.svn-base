/**
 * @brief Reads sum of initial events from EventCount.root.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/EventCount.h"

std::map<std::string, TFile *> f = std::map<std::string, TFile *>();
std::map<std::string, TTree *> t = std::map<std::string, TTree *>();

float getEventCountBeforeSkimming(int mc_channel_number, const std::string &suff) {
  int type;
  float value;
  t[suff]->SetBranchAddress("type", &type);
  t[suff]->SetBranchAddress("value", &value);
  for (int i = 0; i < t[suff]->GetEntries(); ++i) {
    t[suff]->GetEntry(i);
    if (type == mc_channel_number) {
      return value;
    }
  }
  return 1e10;
}

void initEventCount() {
  // use different keys when using PDF systematics .. for now just this key
  f.insert(std::pair<std::string, TFile *>("", new TFile("EventCount.root")));
  for (std::map<std::string, TFile *>::const_iterator it = f.begin(); it != f.end(); ++it) {
    t.insert(std::pair<std::string, TTree *>(it->first, (TTree *) it->second->Get("count")));
  }
}

