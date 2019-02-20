/**
 * @brief Analysis class for tt resonances.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */

#include "TopNtupleAnalysis/Analysis.h"
#include "TopNtupleAnalysis/AnaTtresSLMtt.h"
#include "TopNtupleAnalysis/Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "TopNtupleAnalysis/HistogramService.h"

AnaTtresSLMtt::AnaTtresSLMtt(const std::string &filename, bool electron, bool boosted, std::vector<std::string> &systList)
  : Analysis(filename, systList), m_electron(electron), m_boosted(boosted),
    m_neutrinoBuilder("MeV"), m_chi2("MeV") {

  m_chi2.Init(TtresChi2::DATA2015_MC15C);

  m_hSvc.m_tree->Branch("truemtt",    &_tree_truemtt);
  m_hSvc.m_tree->Branch("mtt",    &_tree_mtt);
  m_hSvc.m_tree->Branch("weight", &_tree_weight);
  m_hSvc.m_tree->Branch("cat",    &_tree_cat);
  m_hSvc.m_tree->Branch("syst",   &_tree_syst);
}

AnaTtresSLMtt::~AnaTtresSLMtt() {
}

void AnaTtresSLMtt::run(const Event &evt, double weight, const std::string &s) {
  // check channel
  if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
    return;

  if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
    return;

  if (m_boosted)
    if (!(evt.passes("bejets") || evt.passes("bmujets")))
      return;

  if (!m_boosted)
    if (!(evt.passes("rejets") || evt.passes("rmujets")))
      return;

  if (evt.channelNumber() == 410000) // for SM ttbar, we have mtt sliced samples above 1.1 TeV
    if (evt.MC_ttbar_beforeFSR().M() > 1.1e6)
      return;
    
  HistogramService *h = &m_hSvc;
  TLorentzVector l;
  if (m_electron) {
    l = evt.electron()[0].mom();
  } else {
    l = evt.muon()[0].mom();
  }

  _tree_truemtt = evt.MC_ttbar_beforeFSR().M()*1e-3;
  float mtt = -1;
  if (m_boosted && (evt.passes("bejets") || evt.passes("bmujets"))) {
    
    size_t close_idx = 0;
    for (; close_idx < evt.jet().size(); ++close_idx)
      if (evt.jet()[close_idx].closeToLepton())
        break;
    const TLorentzVector &sj = evt.jet()[close_idx].mom();
    
    size_t goodljet_idx = 0;
    for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx)
      if (evt.largeJet()[goodljet_idx].good())
        break;
    
    const TLorentzVector &lj = evt.largeJet()[goodljet_idx].mom();
    
    // recalc. mtt
    // lepton = l
    // large-R jet = hadronic top = lj
    // selected jet = leptonic top's b-jet = sj
    // neutrino px, py = met
    std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, evt.met().Perp(), evt.met().Phi(), true);   
    TLorentzVector nu(0,0,0,0);
    if (vec_nu.size() > 0) {
      nu = *(vec_nu[0]);
      for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
      vec_nu.clear();
    }
    
    if (evt.largeJet().size()!=0)    mtt = (lj+sj+nu+l).M();
    _tree_mtt = mtt*1e-3;
    _tree_weight = weight;
    _tree_cat = -1;
    if (evt.passes("bejets") && m_boosted && m_electron) _tree_cat = 0;
    if (evt.passes("bmujets") && m_boosted && !m_electron) _tree_cat = 1;
    _tree_syst = s;
    h->m_tree->Fill();
    
  } else if (!m_boosted && (evt.passes("rejets") || evt.passes("rmujets"))) {
    
    // inputs 
    // LEPTON --> TLorentzVector for your lepton
    // vjets -->  std::vector<TLorentzVector*> for the jets
    // vjets_btagged --> std::vector<bool> to say if the jets are btagged or not
    // met --> TLorentzVector for your MET

    // outputs, they will be filled by the TTBarLeptonJetsBuilder_chi2
    int  igj3, igj4; // index for the Whad
    int igb3, igb4; // index for the b's
    int  ign1;  // index for the neutrino (because chi2 can test both pz solution)
    double chi2ming1, chi2ming1H, chi2ming1L;

    std::vector<TLorentzVector *> vjets;
    std::vector<bool> vjets_btagged;
    for (size_t z = 0; z < evt.jet().size(); ++z) {
      vjets.push_back(new TLorentzVector(0,0,0,0));
      vjets[z]->SetPtEtaPhiE(evt.jet()[z].mom().Perp(), evt.jet()[z].mom().Eta(), evt.jet()[z].mom().Phi(), evt.jet()[z].mom().E());
      // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BTagingxAODEDM
      // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks
      vjets_btagged.push_back(evt.jet()[z].btag_mv2c20_60());
    }
    TLorentzVector met = evt.met();
    bool status = m_chi2.findMinChiSquare(&l, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L); 
    
    float chi2Value = 1000000; // log10(1000000) = 6
    float mwh = -1;
    float mtl = -1;
    float mth = -1;

    if (status){
        chi2Value = chi2ming1;
	mwh = m_chi2.getResult_Mwh();
	mtl = m_chi2.getResult_Mtl();
	mth = m_chi2.getResult_Mth();
	mtt = m_chi2.getResult_Mtt();
    }
    
    for (size_t z = 0; z < vjets.size(); ++z) {
      delete vjets[z];
    }
    vjets.clear();
    vjets_btagged.clear();

    //Fill histograms
    _tree_mtt = mtt*1e-3;
    _tree_weight = weight;
    _tree_syst = s;
    _tree_cat = -1;
    if (evt.passes("rejets") && !m_boosted && m_electron) _tree_cat = 2;
    if (evt.passes("rmujets") && !m_boosted && !m_electron) _tree_cat = 3;
    h->m_tree->Fill();
  }

}

