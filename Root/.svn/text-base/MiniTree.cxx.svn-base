/**
 * @brief Class that reads information from the input file and puts it in a object-oriented format in Event.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "TopNtupleAnalysis/MiniTree.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <vector>
#include <string>
#include <iostream>
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include <algorithm>

#include "TTree.h"

MiniTree::MiniTree(bool toWrite, const std::string &file, const std::string &name)
  : m_file(0), m_chain(0), m_name(name) {

  //m_file = TFile::Open(file.c_str());
  m_chain = new TChain(name.c_str());
  ((TChain *) m_chain)->Add(file.c_str());

  m_sumWeights = 0;

  prepareBranches();
}

MiniTree::~MiniTree() {
  if (m_chain) delete m_chain;
  if (m_file) delete m_file;
}

void MiniTree::read(int event, Event &e) {
  e.clear();
  m_chain->GetEntry(event);
  e.channelNumber() = ui("mcChannelNumber");
  e.isData()= e.channelNumber()==0;
  e.eventNumber() = ul64("eventNumber"); 
  e.runNumber() = ui("runNumber"); 
  e.randomRunNumber() = ui("randomRunNumber"); 

  for (std::map<std::string, char>::iterator it = m_c.begin(); it != m_c.end(); ++it) {
    if (it->first.find("HLT_") == 0) {
      e.setTrigger(it->first, it->second != 0);
    }
  }

  e.npv() = ui("npv"); 
  e.vtxz() = f("vtxz");
  e.mu() = f("mu");
  e.mu_original() = f("mu_original");
  e.weight_mc() = f("weight_mc");
  e.weight_pileup() = f("weight_pileup");
  e.weight_bTagSF() = f("weight_bTagSF_70");
  e.weight_leptonSF() = f("weight_leptonSF");
  e.weight_Sherpa22_corr() = f("weight_Sherpa22_corr");
  e.Wfilter_Sherpa_nT() = i("Wfilter_Sherpa_nT");
  
  // adding the truth information into the event        
  if (f("MC_w1h_pt") > 0)	e.MC_w1h().SetPtEtaPhiM(f("MC_w1h_pt"), f("MC_w1h_eta"), f("MC_w1h_phi"), f("MC_w1h_m"));
  else			e.MC_w1h().SetPtEtaPhiM(2000, -9., -9., 0);  
  e.MC_w1h_pdgId() 	= i("MC_w1h_pdgId");
  
  if (f("MC_w2h_pt") > 0)	e.MC_w2h().SetPtEtaPhiM(f("MC_w2h_pt"), f("MC_w2h_eta"), f("MC_w2h_phi"), f("MC_w2h_m"));
  else			e.MC_w2h().SetPtEtaPhiM(2000, -9., -9., 0);  
  e.MC_w2h_pdgId() 	= i("MC_w2h_pdgId");
  
  if (f("MC_bh_pt") > 0)	e.MC_bh().SetPtEtaPhiM(f("MC_bh_pt"), f("MC_bh_eta"), f("MC_bh_phi"), f("MC_bh_m"));
  else			e.MC_bh().SetPtEtaPhiM(2000, -9., -9., 0);  
  
  if (f("MC_w1l_pt") > 0)	e.MC_w1l().SetPtEtaPhiM(f("MC_w1l_pt"), f("MC_w1l_eta"), f("MC_w1l_phi"), f("MC_w1l_m"));
  else			e.MC_w1l().SetPtEtaPhiM(2000, -9., 0., 0.);  
  e.MC_w1l_pdgId() 	= i("MC_w1l_pdgId");
  
  if (f("MC_w2l_pt") > 0)	e.MC_w2l().SetPtEtaPhiM(f("MC_w2l_pt"), f("MC_w2l_eta"), f("MC_w2l_phi"), f("MC_w2l_m"));
  else			e.MC_w2l().SetPtEtaPhiM(2000, -9., 0., 0.);
  e.MC_w2l_pdgId() 	= i("MC_w2l_pdgId");
  
  if (f("MC_bl_pt") > 0)	e.MC_bl().SetPtEtaPhiM(f("MC_bl_pt"), f("MC_bl_eta"), f("MC_bl_phi"), f("MC_bl_m"));
  else			e.MC_bl().SetPtEtaPhiM(2000, -9., 0., 0.);
  
  if (f("MC_ttbar_beforeFSR_pt") > 0)	e.MC_ttbar_beforeFSR().SetPtEtaPhiM(f("MC_ttbar_beforeFSR_pt"), f("MC_ttbar_beforeFSR_eta"), f("MC_ttbar_beforeFSR_phi"), f("MC_ttbar_beforeFSR_m"));
  else				e.MC_ttbar_beforeFSR().SetPtEtaPhiM(2000, -9., 0., 0.);

  if (f("MC_t_pt") > 0)	e.MC_t().SetPtEtaPhiM(f("MC_t_pt"), f("MC_t_eta"), f("MC_t_phi"), f("MC_t_m"));
  else				e.MC_t().SetPtEtaPhiM(2000, -9., 0., 0.);

  if (f("MC_tbar_pt") > 0)	e.MC_tbar().SetPtEtaPhiM(f("MC_tbar_pt"), f("MC_tbar_eta"), f("MC_tbar_phi"), f("MC_tbar_m"));
  else				e.MC_tbar().SetPtEtaPhiM(2000, -9., 0., 0.);
  
  if (i("MC_ttbar_type") > 0)	e.MC_ttbar_type() = i("MC_ttbar_type");
  
  
    // adding the truth information into the event        
  if (f("MC_Wdecay1_from_t_pt") > 0)	e.MC_Wdecay1_from_t().SetPtEtaPhiM(f("MC_Wdecay1_from_t_pt"), f("MC_Wdecay1_from_t_eta"), f("MC_Wdecay1_from_t_phi"), f("MC_Wdecay1_from_t_m"));
  else					e.MC_Wdecay1_from_t().SetPtEtaPhiM(2000, -9., -9., 0);  
  e.MC_Wdecay1_from_t_pdgId() 	= i("MC_Wdecay1_from_t_pdgId");
  
  if (f("MC_Wdecay2_from_t_pt") > 0)	e.MC_Wdecay2_from_t().SetPtEtaPhiM(f("MC_Wdecay2_from_t_pt"), f("MC_Wdecay2_from_t_eta"), f("MC_Wdecay2_from_t_phi"), f("MC_Wdecay2_from_t_m"));
  else					e.MC_Wdecay2_from_t().SetPtEtaPhiM(2000, -9., -9., 0);  
  e.MC_Wdecay2_from_t_pdgId() 	= i("MC_Wdecay2_from_t_pdgId");
  
  if (f("MC_b_from_t_pt") > 0)	e.MC_b_from_t().SetPtEtaPhiM(f("MC_b_from_t_pt"), f("MC_b_from_t_eta"), f("MC_b_from_t_phi"), f("MC_b_from_t_m"));
  else				e.MC_b_from_t().SetPtEtaPhiM(2000, -9., -9., 0);  
    
  if (f("MC_Wdecay1_from_tbar_pt") > 0)	e.MC_Wdecay1_from_tbar().SetPtEtaPhiM(f("MC_Wdecay1_from_tbar_pt"), f("MC_Wdecay1_from_tbar_eta"), f("MC_Wdecay1_from_tbar_phi"), f("MC_Wdecay1_from_tbar_m"));
  else					e.MC_Wdecay1_from_tbar().SetPtEtaPhiM(2000, -9., -9., 0);  
  e.MC_Wdecay1_from_tbar_pdgId() 	= i("MC_Wdecay1_from_tbar_pdgId");
  
  if (f("MC_Wdecay2_from_tbar_pt") > 0)	e.MC_Wdecay2_from_tbar().SetPtEtaPhiM(f("MC_Wdecay2_from_tbar_pt"), f("MC_Wdecay2_from_tbar_eta"), f("MC_Wdecay2_from_tbar_phi"), f("MC_Wdecay2_from_tbar_m"));
  else					e.MC_Wdecay2_from_tbar().SetPtEtaPhiM(2000, -9., -9., 0);  
  e.MC_Wdecay2_from_tbar_pdgId() 	= i("MC_Wdecay2_from_tbar_pdgId");
  
  if (f("MC_b_from_tbar_pt") > 0)	e.MC_b_from_tbar().SetPtEtaPhiM(f("MC_b_from_tbar_pt"), f("MC_b_from_tbar_eta"), f("MC_b_from_tbar_phi"), f("MC_b_from_tbar_m"));
  else					e.MC_b_from_tbar().SetPtEtaPhiM(2000, -9., -9., 0);  
      
  // adding the matched objets into the event 
  if (f("MA_w1h_pt") > 0)	e.MA_w1h().SetPtEtaPhiM(f("MA_w1h_pt"), f("MA_w1h_eta"), f("MA_w1h_phi"), f("MA_w1h_m"));
  else   		e.MA_w1h().SetPtEtaPhiM(2000, -9., 0., 0.);
  e.MA_w1h_pdgId() 	= i("MA_w1h_pdgId");
  
  if (f("MA_w2h_pt") > 0)	e.MA_w2h().SetPtEtaPhiM(f("MA_w2h_pt"), f("MA_w2h_eta"), f("MA_w2h_phi"), f("MA_w2h_m"));
  else   		e.MA_w2h().SetPtEtaPhiM(2000, -9., 0., 0.);
  e.MA_w2h_pdgId() 	= i("MA_w2h_pdgId");
  
  if (f("MA_bh_pt") > 0)	e.MA_bh().SetPtEtaPhiM(f("MA_bh_pt"), f("MA_bh_eta"), f("MA_bh_phi"), f("MA_bh_m"));
  else   		e.MA_bh().SetPtEtaPhiM(2000, -9., 0., 0.);
  
  if (f("MA_w1l_pt") > 0)	e.MA_w1l().SetPtEtaPhiM(f("MA_w1l_pt"), f("MA_w1l_eta"), f("MA_w1l_phi"), f("MA_w1l_m"));
  else   		e.MA_w1l().SetPtEtaPhiM(2000, -9., 0., 0.);
  e.MA_w1l_pdgId() 	= i("MA_w1l_pdgId");
  
  if (f("MA_w2l_pt") > 0)	e.MA_w2l().SetPtEtaPhiM(f("MA_w2l_pt"), f("MA_w2l_eta"), f("MA_w2l_phi"), f("MA_w2l_m"));
  else			e.MA_w2l().SetPtEtaPhiM(2000, -9., 0., 0.);
  e.MA_w2l_pdgId() 	= i("MA_w2l_pdgId");
  
  if (f("MA_bl_pt") > 0)	e.MA_bl().SetPtEtaPhiM(f("MA_bl_pt"), f("MA_bl_eta"), f("MA_bl_phi"), f("MA_bl_m"));  
  else   		e.MA_bl().SetPtEtaPhiM(2000, -9., 0., 0.);

  for (int k = 0; k < vf("el_pt")->size(); ++k) {
    e.electron().push_back(Electron());
    e.electron()[k].mom().SetPtEtaPhiE(vf("el_pt")->at(k), vf("el_eta")->at(k), vf("el_phi")->at(k), vf("el_e")->at(k));
    //e.electron()[k].setMI(vf("el_miniiso")->at(k));
    e.electron()[k].setTightPP(false);
    if (vc("el_isTight"))	e.electron()[k].setTightPP(vc("el_isTight")->at(k));
    e.electron()[k].caloMom() = e.electron()[k].mom();
    e.electron()[k].trkMom() = e.electron()[k].mom();
    e.electron()[k].Dz0() = vf("el_delta_z0_sintheta")->at(k);
    e.electron()[k].d0() = vf("el_d0")->at(k);
    e.electron()[k].sd0() = vf("el_d0sig")->at(k); 
    e.electron()[k].ptvarcone20() = vf("el_ptvarcone20")->at(k);
    e.electron()[k].topoetcone20()= vf("el_topoetcone20")->at(k);  
    e.electron()[k].HLT_e24_lhmedium_iloose_L1EM20VH()= 0;
    e.electron()[k].HLT_e24_lhmedium_L1EM18VH() = 0;
    e.electron()[k].HLT_e24_lhmedium_L1EM20VH() = 0;
    e.electron()[k].HLT_e60_lhmedium() 		= 0;
    e.electron()[k].HLT_e120_lhloose() 		= 0;
    if (vc("el_trigMatch_HLT_e24_lhmedium_iloose_L1EM20VH"))	e.electron()[k].HLT_e24_lhmedium_iloose_L1EM20VH()= vc("el_trigMatch_HLT_e24_lhmedium_iloose_L1EM20VH")->at(k);
    if (vc("el_trigMatch_HLT_e24_lhmedium_L1EM18VH"))		e.electron()[k].HLT_e24_lhmedium_L1EM18VH() = vc("el_trigMatch_HLT_e24_lhmedium_L1EM18VH")->at(k); 
    if (vc("el_trigMatch_HLT_e24_lhmedium_L1EM20VH"))		e.electron()[k].HLT_e24_lhmedium_L1EM20VH() = vc("el_trigMatch_HLT_e24_lhmedium_L1EM20VH")->at(k); 
    e.electron()[k].HLT_e60_lhmedium() 		= vc("el_trigMatch_HLT_e60_lhmedium")->at(k); 
    e.electron()[k].HLT_e120_lhloose() 		= vc("el_trigMatch_HLT_e120_lhloose")->at(k); 
    e.electron()[k].author() = 1;
  }
  for (int k = 0; k < vf("mu_pt")->size(); ++k) {
    e.muon().push_back(Muon());
    e.muon()[k].mom().SetPtEtaPhiE(vf("mu_pt")->at(k), vf("mu_eta")->at(k), vf("mu_phi")->at(k), vf("mu_e")->at(k));
    e.muon()[k].setMI(0);
    e.muon()[k].setTight(false);
    if (vc("mu_isTight"))	e.muon()[k].setTight(vc("mu_isTight")->at(k));
    e.muon()[k].Dz0() = vf("mu_delta_z0_sintheta")->at(k);
    e.muon()[k].d0() = vf("mu_d0")->at(k);
    e.muon()[k].sd0() = vf("mu_d0sig")->at(k);
    e.muon()[k].ptvarcone30() = vf("mu_ptvarcone30")->at(k);
    e.muon()[k].topoetcone20()= vf("mu_topoetcone20")->at(k); 
    if (vc("mu_trigMatch_HLT_mu20_L1MU15"))		e.muon()[k].HLT_mu20_L1MU15() 	 = vc("mu_trigMatch_HLT_mu20_L1MU15")->at(k);
    if (vc("mu_trigMatch_HLT_mu20_iloose_L1MU15"))	e.muon()[k].HLT_mu20_iloose_L1MU15() = vc("mu_trigMatch_HLT_mu20_iloose_L1MU15")->at(k); 
    e.muon()[k].HLT_mu50() 		 = vc("mu_trigMatch_HLT_mu50")->at(k);   
    e.muon()[k].author() = 0;
    e.muon()[k].passTrkCuts() = true;
  }

  for (int k = 0; k < vf("jet_pt")->size(); ++k) {
    e.jet().push_back(Jet());
    e.jet()[k].mom().SetPtEtaPhiE(vf("jet_pt")->at(k), vf("jet_eta")->at(k), vf("jet_phi")->at(k), vf("jet_e")->at(k));
    e.jet()[k].trueFlavour() = vi("jet_trueflav")->at(k); // TODO jet_trueflav==0?-99:jet_trueflav->at(k);
    e.jet()[k].mv1() = vf("jet_mv1")==0?-99:vf("jet_mv1")->at(k);
    e.jet()[k].ip3dsv1() = vf("jet_ip3dsv1")==0?-99:vf("jet_ip3dsv1")->at(k);
    e.jet()[k].mv2c20() = vf("jet_mv2c20")==0?-99:vf("jet_mv2c20")->at(k);
    e.jet()[k].jvt() = vf("jet_jvt")==0?-99:vf("jet_jvt")->at(k);
    e.jet()[k].closeToLepton() = vi("jet_closeToLepton")==0?-99:vi("jet_closeToLepton")->at(k);
  }

  for (int k = 0; k < vf("tjet_pt")->size(); ++k) {
    e.tjet().push_back(Jet());
    e.tjet()[k].mom().SetPtEtaPhiE(vf("tjet_pt")->at(k), vf("tjet_eta")->at(k), vf("tjet_phi")->at(k), vf("tjet_e")->at(k));
    e.tjet()[k].trueFlavour() = vi("tjet_label")->at(k);
    e.tjet()[k].mv1() = -99;
    e.tjet()[k].ip3dsv1() = -99;
    e.tjet()[k].mv2c10() = vf("tjet_mv2c10")==0?-99:vf("tjet_mv2c10")->at(k);
    e.tjet()[k].jvt() = -99;
    e.tjet()[k].closeToLepton() = -99;
    e.tjet()[k].numConstituents() = vi("tjet_numConstituents")->at(k);
  }
   
  for (int k = 0; k < vf("ljet_pt")->size(); ++k) {
    e.largeJet().push_back(LargeJet());
    e.largeJet()[k].mom().SetPtEtaPhiE(vf("ljet_pt")->at(k), vf("ljet_eta")->at(k), vf("ljet_phi")->at(k), vf("ljet_e")->at(k));
    e.largeJet()[k].split12() = vf("ljet_sd12")->at(k);
    e.largeJet()[k].setGood((vi("ljet_good")->at(k) == 1)?true:false);
    e.largeJet()[k].trueFlavour() = 0; //TODO ljet_trueflav==0?-99:ljet_trueflav->at(k);
    e.largeJet()[k].subs("tau32") = vf("ljet_tau32")->at(k);
    e.largeJet()[k].subs("tau32_wta") = vf("ljet_tau32_wta")->at(k);
    e.largeJet()[k].subs("tau21") = vf("ljet_tau21")->at(k);
    e.largeJet()[k].subs("tau21_wta") = vf("ljet_tau21_wta")->at(k);
  }
  
  e.met(f("met_met")*std::cos(f("met_phi")), f("met_met")*std::sin(f("met_phi")));

  e.passes().clear();

  if (i("bejets")) e.passes().push_back("bejets");
  if (i("bmujets")) e.passes().push_back("bmujets");
  if (i("bejets_2015")) e.passes().push_back("bejets_2015");
  if (i("bmujets_2015")) e.passes().push_back("bmujets_2015");
  if (i("bejets_2016")) e.passes().push_back("bejets_2016");
  if (i("bmujets_2016")) e.passes().push_back("bmujets_2016");
  if (i("bmu2jets")) e.passes().push_back("bmu2jets");

  if (i("rejets")) e.passes().push_back("rejets");
  if (i("rmujets")) e.passes().push_back("rmujets");
  if (i("rejets_2015")) e.passes().push_back("rejets_2015");
  if (i("rmujets_2015")) e.passes().push_back("rmujets_2015");
  if (i("rejets_2016")) e.passes().push_back("rejets_2016");
  if (i("rmujets_2016")) e.passes().push_back("rmujets_2016");
  if (i("rmu2jets")) e.passes().push_back("rmu2jets");
  
  if (i("rejetsQCDCR_2015")) 	e.passes().push_back("rejetsQCDCR_2015");
  if (i("bejetsQCDCR_2015")) 	e.passes().push_back("bejetsQCDCR_2015");
  if (i("rmujetsQCDCR_2015")) 	e.passes().push_back("rmujetsQCDCR_2015");
  if (i("bmujetsQCDCR_2015")) 	e.passes().push_back("bmujetsQCDCR_2015");  
  if (i("rejetsIncluR_2015")) 	e.passes().push_back("rejetsIncluR_2015");
  if (i("bejetsIncluR_2015")) 	e.passes().push_back("bejetsIncluR_2015");
  if (i("rejetsIncluR_2016"))   e.passes().push_back("rejetsIncluR_2016");
  if (i("bejetsIncluR_2016"))   e.passes().push_back("bejetsIncluR_2016");

  if (i("rejetsWCR_2015")) 	e.passes().push_back("rejetsWCR_2015");
  if (i("bejetsWCR_2015")) 	e.passes().push_back("bejetsWCR_2015");
  if (i("rmujetsWCR_2015")) 	e.passes().push_back("rmujetsWCR_2015");
  if (i("bmujetsWCR_2015")) 	e.passes().push_back("bmujetsWCR_2015");
  if (i("rejetsWCR_2016")) 	e.passes().push_back("rejetsWCR_2016");
  if (i("bejetsWCR_2016")) 	e.passes().push_back("bejetsWCR_2016");
  if (i("rmujetsWCR_2016")) 	e.passes().push_back("rmujetsWCR_2016");
  if (i("bmujetsWCR_2016")) 	e.passes().push_back("bmujetsWCR_2016");

  if (i("rmujetsQCDCR_2016")) 	e.passes().push_back("rmujetsQCDCR_2016");
  if (i("bmujetsQCDCR_2016")) 	e.passes().push_back("bmujetsQCDCR_2016");
  if (i("rejetsQCDCR_2016")) 	e.passes().push_back("rejetsQCDCR_2016");
  if (i("bejetsQCDCR_2016")) 	e.passes().push_back("bejetsQCDCR_2016");

  if (i("ejets")) e.passes().push_back("ejets");
  if (i("mujets")) e.passes().push_back("mujets");
  if (i("ee")) e.passes().push_back("ee");
  if (i("emu")) e.passes().push_back("emu");
  if (i("mumu")) e.passes().push_back("mumu");
  
  
  
  
}

ULong64_t          &MiniTree::ul64(const std::string &n) {
  return m_ul64[n];
}
unsigned int         &MiniTree::ui(const std::string &n) {
  return m_ui[n];
}
int                  &MiniTree::i(const std::string &n) {
  return m_i[n];
}
float                &MiniTree::f(const std::string &n) {
  return m_f[n];
}
char                 &MiniTree::c(const std::string &n) {
  return m_c[n];
}

std::vector<std::vector<float> > *MiniTree::vvf(const std::string &n) {
  if (m_vvf.find(n) == m_vvf.end()) return 0;
  return m_vvf[n];
}

std::vector<std::vector<int> > *MiniTree::vvi(const std::string &n) {
  if (m_vvi.find(n) == m_vvi.end()) return 0;
  return m_vvi[n];
}
std::vector<int>     *MiniTree::vi(const std::string &n) {
  if (m_vi.find(n) == m_vi.end()) return 0;
  return m_vi[n];
}
std::vector<char>    *MiniTree::vc(const std::string &n) {
  if (m_vc.find(n) == m_vc.end()) return 0;
  return m_vc[n];
}
std::vector<float>   *MiniTree::vf(const std::string &n) {
  if (m_vf.find(n) == m_vf.end()) return 0;

  return m_vf[n];
}

double &MiniTree::sumWeights() {
  return m_sumWeights;
}

void MiniTree::addFileToRead(const std::string &fname) {
  ((TChain *) m_chain)->Add(fname.c_str());
}

void MiniTree::addFileToRead(const std::string &fname, const std::string &treeName) {
  ((TChain *) m_chain)->AddFile(fname.c_str(), -1, treeName.c_str());
}

double MiniTree::getSumWeights() {
  return m_sumWeights;
}

int MiniTree::GetEntries() {
  return m_chain->GetEntries();
}

void MiniTree::prepareBranches() {
  TObjArray *l = m_chain->GetListOfBranches();
  for (size_t z = 0; z < l->GetEntries(); ++z) {
    std::string name = l->At(z)->GetName();
    std::string type = ((TBranch *) l->At(z))->GetLeaf(name.c_str())->GetTypeName();
    //std::cout<< name << "\t" << type << std::endl;
    if (type == "Float_t" || type == "float") {
      m_f[name] = -50000;
    } else if (type == "Int_t" || type == "int") {
      m_i[name] = -50000;
    } else if (type == "ULong64_t") {
      m_ul64[name] = 0;
    } else if (type == "UInt_t" || type == "unsigned int") {
      m_ui[name] = 0;
    } else if (type == "char" || type == "Char_t") {
      m_c[name] = 0;
    } else if (type == "vector<char>" || type == "std::vector<char>" || type == "vector<Char_t>" || type == "std::vector<Char_t>") {
      m_vc[name] = 0;
    } else if (type == "vector<float>" || type == "std::vector<float>" || type == "vector<Float_t>" || type == "std::vector<Float_t>") {
      m_vf[name] = 0;
    } else if (type == "vector<int>" || type == "std::vector<int>" || type == "vector<Int_t>" || type == "std::vector<Int_t>") {
      m_vi[name] = 0;
    } else if (type == "vector<vector<int> >" || type == "std::vector<std::vector<int> >" || type == "vector<vector<Int_t> >" || type == "std::vector<std::vector<Int_t> >") {
      m_vvi[name] = 0;
    } else if (type == "vector<vector<float> >" || type == "std::vector<std::vector<float> >" || type == "vector<vector<Float_t> >" || type == "std::vector<std::vector<Float_t> >") {
      m_vvf[name] = 0;
    } else {
      std::cout << "ERROR: I could not figure out the type of this branch! Name = " << name << ", type = " << type << std::endl;
    }
  }

  //for (auto& it : m_f) {
  for (std::map<std::string, float>::iterator it = m_f.begin(); it != m_f.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtFloat;
  }
  //for (auto& it : m_ui) {
  for (std::map<std::string, unsigned int>::iterator it = m_ui.begin(); it != m_ui.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtUint;
  }
  //for (auto& it : m_ul64) {
  for (std::map<std::string, ULong64_t>::iterator it = m_ul64.begin(); it != m_ul64.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtULong64;
  }
  //for (auto& it : m_i) {
  for (std::map<std::string, int>::iterator it = m_i.begin(); it != m_i.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtInt;
  }
  //for (auto& it : m_c) {
  for (std::map<std::string, char>::iterator it = m_c.begin(); it != m_c.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtChar;
  }
  //for (auto& it : m_vc) {
  for (std::map<std::string, std::vector<char> *>::iterator it = m_vc.begin(); it != m_vc.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtVChar;
  }
  //for (auto& it : m_vf) {
  for (std::map<std::string, std::vector<float> *>::iterator it = m_vf.begin(); it != m_vf.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtVFloat;
  }
  //for (auto& it : m_vi) {
  for (std::map<std::string, std::vector<int> *>::iterator it = m_vi.begin(); it != m_vi.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtVInt;
  }
  //for (auto& it : m_vvi) {
  for (std::map<std::string, std::vector<std::vector<int> > *>::iterator it = m_vvi.begin(); it != m_vvi.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtVVInt;
  }
  //for (auto& it : m_vvf) {
  for (std::map<std::string, std::vector<std::vector<float> > *>::iterator it = m_vvf.begin(); it != m_vvf.end(); ++it) {
    m_chain->SetBranchAddress(it->first.c_str(), &(it->second));
    m_brs[it->first] = mtVVFloat;
  }
}

