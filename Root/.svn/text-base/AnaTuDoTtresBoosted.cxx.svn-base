/**
 * @brief Analysis class for TU Dortmund to run ttbar (boosted) analysis
 * @author Hendrik Esch <hendrik.esch@tu-dortmund.de>
 */

#include "TopNtupleAnalysis/AnaTuDoTtresBoosted.h"
#include <vector>
#include <string>


AnaTuDoTtresBoosted::AnaTuDoTtresBoosted(const std::string &filename, bool electron, MiniTree *mini)
  : Analysis(filename), m_electron(electron), m_isData(false),
    m_neutrinoBuilder("MeV") {

    TH1::SetDefaultSumw2();
        
//     m_histoSvc = HistoListSvc::svc(filename, "foo.xml", mini->m_chain); 
    m_mini = mini;

//     m_leppt = m_histoSvc->th1F("leppt", "p_T of lepton", 20, 0, 200);
//     m_lepeta = m_histoSvc->th1F("lepeta", "#eta of lepton", 20, -5, 5);
//     m_lepphi = m_histoSvc->th1F("lepphi", "#phi of lepton", 14, -3.5, 3.5);
//     m_lepcharge = m_histoSvc->th1F("lepcharge", "charge of lepton", 5, -2, 2);
//     
//     m_leadJetPt = m_histoSvc->th1F("leadJetPt", "; Leading Jet p_{T} ; Events", 50, 0, 500);
//     m_leadbJetPt = m_histoSvc->th1F("leadbJetPt", "; Leading b-jet p_{T} ; Events", 50, 0, 500);
//     
//     m_met = m_histoSvc->th1F("met", "; Missing E_{T} ; Events", 50, 0, 500);
//     m_met_phi = m_histoSvc->th1F("met_phi", "; Missing E_{T} #phi; Events", 64, -3.2, 3.2);
// 
//     m_closeJetPt = m_histoSvc->th1F("closeJetPt", "; Selected Jet p_{T} ; Events", 25, 0, 500);
//     m_largeJetPt = m_histoSvc->th1F("largeJetPt", "; Large jet p_{T} ; Events", 20, 300, 700);
//     m_largeJetM = m_histoSvc->th1F("largeJetM", "; Large jet M ; Events", 15, 0, 300);
//     m_largeJetSd12 = m_histoSvc->th1F("largeJetSd12", "; Large jet #sqrt{d_{12}} ; Events", 30, 0, 300);
    
//     m_mtlep_boo = m_histoSvc->th1F("mtlep_boo", "; m_{t,lep} ; Events", 40, 0, 400);
    
//     m_mtt = m_histoSvc->th1F("mtt", "; m_{t#bar{t}} ; Events", 60, 0, 6000);
     
    h_TuDoCutFlow = new TH1F("TuDoCutFlow", "TuDoCutFlow", 14, 0, 14);
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(1,"INITIAL");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(2,"HFOR veto");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(3,"Preselection");

}

AnaTuDoTtresBoosted::~AnaTuDoTtresBoosted() {
}

void AnaTuDoTtresBoosted::setIsData(bool isData){
    m_isData = isData;
}

void AnaTuDoTtresBoosted::run(const Event &evt, double weight, const std::string &s) {
    
    h_TuDoCutFlow -> Fill(0);
    h_TuDoCutFlow -> Fill(1);

    // check channel
    if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
        return;

    if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
        return;
    
    if (!(evt.passes("bejets") || evt.passes("bmujets")))
        return;

    h_TuDoCutFlow -> Fill(2);
    
    TLorentzVector lepton;
    float lepcharge = -99;
    if (m_electron) {
        lepton = evt.electron()[0].mom();
//         lepcharge = evt.electron()[0].charge();
    } else {
        lepton = evt.muon()[0].mom();
//         lepcharge = evt.electron()[0].charge();
    }
    
    h_TuDoCutFlow -> Fill(3);

//     m_mini->binning01 = -1;
//     m_mini->binning02 = abs(lepton.Eta());
    
    const TLorentzVector &leadingjet = evt.jet()[0].mom();
    int index_leadingjet = 0;
    
//     m_leppt->Fill(lepton.Perp()*1e-3, weight, index_leadingjet);
//     m_lepeta->Fill(lepton.Eta(), weight, index_leadingjet);
//     m_lepphi->Fill(lepton.Phi(), weight, index_leadingjet);
//     m_lepcharge->Fill(lepcharge, weight, index_leadingjet);
    
//     m_met->Fill(evt.met().Perp()*1e-3, weight, index_leadingjet);
//     m_met_phi->Fill(evt.met().Phi(), weight, index_leadingjet);

//     m_leadJetPt->Fill(leadingjet.Perp()*1e-3, weight, index_leadingjet);
    size_t bidx = 0;
    for (; bidx < evt.jet().size(); ++bidx)
        if (evt.jet()[bidx].btag_mv2c20_60())
        break;
//     m_leadbJetPt->Fill(evt.jet()[bidx].mom().Perp()*1e-3, weight, index_leadingjet);

    float mtt = 0;
    

    size_t close_idx = 0;
    for (; close_idx < evt.jet().size(); ++close_idx)
      if (evt.jet()[close_idx].closeToLepton())
        break;
      
    // selected jet = leptonic top's b-jet = sj
    const TLorentzVector &sj = evt.jet()[close_idx].mom();
    
//     m_closeJetPt->Fill(sj.Perp()*1e-3, weight, index_leadingjet);
    
    size_t goodljet_idx = 0;
    for (; goodljet_idx < evt.largeJet().size(); ++goodljet_idx)
      if (evt.largeJet()[goodljet_idx].good())
        break;
    // large-R jet = hadronic top = lj
    const TLorentzVector &lj = evt.largeJet()[goodljet_idx].mom();
    
//     m_largeJetPt->Fill(lj.Perp()*1e-3, weight, index_leadingjet);
//     m_largeJetM->Fill(lj.M()*1e-3, weight, index_leadingjet);
//     m_largeJetSd12->Fill(evt.largeJet()[goodljet_idx].split12()*1e-3, weight, index_leadingjet);

    // recalc. mtt
    // lepton = l
    // neutrino px, py = met
    std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&lepton, evt.met().Perp(), evt.met().Phi(), true);   
    TLorentzVector nu(0,0,0,0);
    if (vec_nu.size() > 0) {
      nu = *(vec_nu[0]);
      for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
      vec_nu.clear();
    }
    mtt = (lj+sj+nu+lepton).M();
    
//     m_mtlep_boo->Fill((sj+nu+lepton).M()*1e-3, weight, index_leadingjet);
//     m_mtt->Fill(mtt*1e-3, weight, index_leadingjet);
}

void AnaTuDoTtresBoosted::terminate() {
    
    h_TuDoCutFlow->SetDirectory(0);
    h_TuDoCutFlow->Write();
    
    TH1D* cutflow;
    if(!m_isData) {
        if(m_electron) {
            cutflow = (TH1D*) m_mini->m_file->Get("bejets/cutflow_mc_pu_zvtx")->Clone("cutflow_mc_pu_zvtx");
        } else {
            cutflow = (TH1D*) m_mini->m_file->Get("bmujets/cutflow_mc_pu_zvtx")->Clone("cutflow_mc_pu_zvtx");
        }
    } else {
        if(m_electron) {
            cutflow = (TH1D*) m_mini->m_file->Get("bejets/cutflow")->Clone("cutflow");
        } else {
            cutflow = (TH1D*) m_mini->m_file->Get("bmujets/cutflow")->Clone("cutflow");
        }

    }
    
    cutflow->SetDirectory(0);
    
    float number_events = 0;
    number_events = cutflow->GetBinContent(1);
    
    TH1D *nEvt = new TH1D("event_counter", "event_counter", 1, 0, 1);
    nEvt->SetBinContent(1, number_events);
    nEvt->Write();

    cutflow->Write();
    
//     std::cout << "number_events = " << number_events << " ; event_number_TRC = " << event_number_TRC << std::endl;
    
//     HistoListSvc::cleanup();
}
