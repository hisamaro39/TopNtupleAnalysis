/**
 * @brief Analysis class for TU Dortmund to run ttbar (resolved) analysis
 * @author Hendrik Esch <hendrik.esch@tu-dortmund.de>
 */

#include "TopNtupleAnalysis/AnaTuDoTtresResolved.h"
#include <vector>
#include <string>


AnaTuDoTtresResolved::AnaTuDoTtresResolved(const std::string &filename, bool electron, MiniTree *mini)
  : Analysis(filename), m_electron(electron), m_isData(false),
    m_neutrinoBuilder("MeV"), m_chi2("MeV") {

    TH1::SetDefaultSumw2();
        
//     m_histoSvc = HistoListSvc::svc(filename, "foo.xml", mini->m_chain); 
    m_mini = mini;
    
    m_chi2.Init(TtresChi2::DATA2015_MC15C);

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
//     m_mtlep_res = m_histoSvc->th1F("mtlep_res", "; m_{t,lep} ; Events", 40, 0, 400);
//     m_mthad_res = m_histoSvc->th1F("mthad_res", "; m_{t,had} ; Events", 40, 0, 400);
//     m_mwhad_res = m_histoSvc->th1F("mwhad_res", "; m_{W,had} ; Events", 40, 0, 400);
//     m_hl_chi2 = m_histoSvc->th1F("chi2", "; #chi^{2} ; Events", 40, 0, 400);
//     
//     m_mtt = m_histoSvc->th1F("mtt", "; m_{t#bar{t}} ; Events", 60, 0, 6000);

    h_TuDoCutFlow = new TH1F("TuDoCutFlow", "TuDoCutFlow", 14, 0, 14);
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(1,"INITIAL");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(2,"HFOR veto");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(3,"Preselection");
}

AnaTuDoTtresResolved::~AnaTuDoTtresResolved() {
}

void AnaTuDoTtresResolved::setIsData(bool isData){
    m_isData = isData;
}

void AnaTuDoTtresResolved::run(const Event &evt, double weight, const std::string &s) {
    
    h_TuDoCutFlow -> Fill(0);
    h_TuDoCutFlow -> Fill(1);

    // check channel
    if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
        return;

    if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
        return;
    
    if (!(evt.passes("rejets") || evt.passes("rmujets")))
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

    m_mini->f("binning01") = -1;
    m_mini->f("binning02") = abs(lepton.Eta());
    
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

    std::vector<TLorentzVector *> jets;
    std::vector<bool> jets_btagged;
    for (size_t it_jet = 0; it_jet < evt.jet().size(); ++it_jet) {
        jets.push_back(new TLorentzVector(0,0,0,0));
        jets[it_jet]->SetPtEtaPhiE(evt.jet()[it_jet].mom().Perp(), evt.jet()[it_jet].mom().Eta(), evt.jet()[it_jet].mom().Phi(), evt.jet()[it_jet].mom().E());
        // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BTagingxAODEDM
        // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks
        jets_btagged.push_back(evt.jet()[it_jet].btag_mv2c20_60());
    }

    TLorentzVector met = evt.met();
    
    bool status = m_chi2.findMinChiSquare(&lepton, &jets, &jets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L); 
    // status has to be true

    if (status) mtt = m_chi2.getResult_Mtt();
    for (size_t z = 0; z < jets.size(); ++z) {
      delete jets[z];
    }
    jets.clear();
    jets_btagged.clear();
    
//     m_mtlep_res->Fill(m_chi2.getResult_Mtl()*1e-3, weight, index_leadingjet);
//     m_mthad_res->Fill(m_chi2.getResult_Mth()*1e-3, weight, index_leadingjet);
//     m_mwhad_res->Fill(m_chi2.getResult_Mwh()*1e-3, weight, index_leadingjet);
//     m_hl_chi2->Fill(chi2ming1, weight, index_leadingjet);
//     m_mtt->Fill(mtt*1e-3, weight, index_leadingjet);
}

void AnaTuDoTtresResolved::terminate() {
    
    h_TuDoCutFlow->SetDirectory(0);
    h_TuDoCutFlow->Write();
    
    TH1D* cutflow;
    if(!m_isData) {
        if(m_electron) {
            cutflow = (TH1D*) m_mini->m_file->Get("rejets/cutflow_mc_pu_zvtx")->Clone("cutflow_mc_pu_zvtx");
        } else {
            cutflow = (TH1D*) m_mini->m_file->Get("rmujets/cutflow_mc_pu_zvtx")->Clone("cutflow_mc_pu_zvtx");
        }
    } else {
        if(m_electron) {
            cutflow = (TH1D*) m_mini->m_file->Get("rejets/cutflow")->Clone("cutflow");
        } else {
            cutflow = (TH1D*) m_mini->m_file->Get("rmujets/cutflow")->Clone("cutflow");
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
