/**
 * @brief Reads sum of initial events from EventCount.root.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef EVENTCOUNT_H
#define EVENTCOUNT_H

#include "TFile.h"
#include "TTree.h"
#include <map>
#include <string>

float getEventCountBeforeSkimming(int mc_channel_number, const std::string &suff = "");
void initEventCount();

#endif
