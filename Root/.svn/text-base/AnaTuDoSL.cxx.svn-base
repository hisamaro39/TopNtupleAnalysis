/**
 * @brief Analysis class for TU Dortmund to run analysis with 2 jets and 1 lepton.
 * @author Hendrik Esch <hendrik.esch@tu-dortmund.de>
 */

#include "TopNtupleAnalysis/AnaTuDoSL.h"
#include <vector>
#include <string>


AnaTuDoSL::AnaTuDoSL(const std::string &filename, bool electron, MiniTree *mini)
  : Analysis(filename), m_electron(electron), m_isData(false),
    m_neutrinoBuilder("MeV") {

    TH1::SetDefaultSumw2();
        
//     m_histoSvc = HistoListSvc::svc(filename, "foo.xml", mini->m_chain); 
    m_mini = mini;
    

//     m_lep_DeltaPhi = m_histoSvc->th2F("hl_lep_DeltaPhi", "lep_DeltaPhi (hl)",140,-3.5,3.5,100,0,200);
//     m_lep_DeltaPhiCut = m_histoSvc->th2F("hl_lep_DeltaPhiCut", "lep_DeltaPhi (hl,cut)",140,-3.5,3.5,100,0,200);
    
//     m_leppt = m_histoSvc->th1F("leppt", "p_T of lepton", 20, 0, 200);
//     m_lepeta = m_histoSvc->th1F("lepeta", "#eta of lepton", 20, -5, 5);
//     m_lepphi = m_histoSvc->th1F("lepphi", "#phi of lepton", 14, -3.5, 3.5);
//     m_lepcharge = m_histoSvc->th1F("lepcharge", "charge of lepton", 5, -2, 2);
//     
//     m_nupt = m_histoSvc->th1F("nupt", "p_T of neutrino", 20, 0, 200);
//     m_nueta = m_histoSvc->th1F("nueta", "#eta of neutrino", 20, -5, 5);
//     m_nuphi = m_histoSvc->th1F("nuphi", "#phi of neutrino", 14, -3.5, 3.5);
//     
//     m_jpt_1 = m_histoSvc->th1F("jpt_1", "p_T of first jet", 20, 0, 200);
//     m_jeta_1 = m_histoSvc->th1F("jeta_1", "#eta of first jet", 20, -5, 5);
//     m_jphi_1 = m_histoSvc->th1F("jphi_1", "#phi of first jet", 14, -3.5, 3.5);
//     m_jpt_2 = m_histoSvc->th1F("jpt_2", "p_T of second jet", 20, 0, 200);
//     m_jeta_2 = m_histoSvc->th1F("jeta_2", "#eta of second jet", 20, -5, 5);
//     m_jphi_2 = m_histoSvc->th1F("jphi_2", "#phi of second jet", 14, -3.5, 3.5);
//     
//     m_bpt = m_histoSvc->th1F("bpt", "p_T of b-jet", 20, 0, 200);
//     m_beta = m_histoSvc->th1F("beta", "#eta of b-jet", 20, -5, 5);
//     m_bphi = m_histoSvc->th1F("bphi", "#phi of b-jet", 14, -3.5, 3.5);
//     m_lightpt = m_histoSvc->th1F("lightpt", "p_T of light-jet", 20, 0, 200);
//     m_lighteta = m_histoSvc->th1F("lighteta", "#eta of light-jet", 20, -5, 5);
//     m_lightphi = m_histoSvc->th1F("lightphi", "#phi of light-jet", 14, -3.5, 3.5);
// 
//     m_met = m_histoSvc->th1F("met", "Missing Energy", 40, 0, 200);
//     m_metx = m_histoSvc->th1F("metx", "Missing Energy (x-component)", 40, -100, 100);
//     m_mety = m_histoSvc->th1F("mety", "Missing Energy (y-component)", 40, -100, 100);
//     m_metphi = m_histoSvc->th1F("metphi", "Missing Energy (#phi)", 14, -3.5, 3.5);
//     m_mtw = m_histoSvc->th1F("mtw", "Transverse W Mass", 30, 0, 150);
// 
//     m_deltaR_l_lnub = m_histoSvc->th1F("deltaR_l_lnub", "deltaR between lepton and top",40, 0, 8);
//     m_pt_lnu = m_histoSvc->th1F("pt_lnu", "p_T of W", 60,0,300);
//     m_pt_lnub = m_histoSvc->th1F("pt_lnub", "p_T of top", 80,0,400);
//     m_mnu = m_histoSvc->th1F("mnu", "m(nu)", 100, -500, 500);
//     m_costheta_lj_top = m_histoSvc->th1F("costheta_lj_top", "polarization", 40, -1.0, 1.0);
//     m_WHelicity = m_histoSvc->th1F("WHelicity", "helicity of W", 40, -1.0, 1.0);
//     m_aplanarity = m_histoSvc->th1F("aplanarity", "aplanarity", 45, 0.0, 0.45);
//     m_sphericity = m_histoSvc->th1F("sphericity", "sphericity", 50, 0.0, 1.0);
//     m_mlnubj = m_histoSvc->th1F("mlnubj", "m(lnubj)", 150, 0, 1500);
//     m_deltaEta_l_b = m_histoSvc->th1F("deltaEta_l_b", "deltaEta between lepton and b-jet",80, 0, 4);
//     m_deltaPhi_l_b = m_histoSvc->th1F("deltaPhi_l_b", "deltaPhi between lepton and b-jet",80, -4, 4);
//     m_deltaR_l_b = m_histoSvc->th1F("deltaR_l_b", "deltaR between lepton and b-jet",120, -6, 6);
//     m_mlnub = m_histoSvc->th1F("mlnub", "m(lnub)", 250, 0, 500);
//     m_mb = m_histoSvc->th1F("mb", "m(b)", 14, 0, 21);
//     m_mjb = m_histoSvc->th1F("mjb", "m(jb)", 25, 0, 500);
//     m_etaj = m_histoSvc->th1F("etaj", "|#eta| of not tagged jet", 25, 0, 5);
//     m_etalnu = m_histoSvc->th1F("etalnu", "#eta of l+nu", 20, -5, 5);
//     m_etalnub = m_histoSvc->th1F("etalnub", "#eta of l+nu+b", 20, -5, 5);
//     m_mlb = m_histoSvc->th1F("mlb", "m(lb)", 120, 0, 300);
//     m_deltaR_j1_lnu = m_histoSvc->th1F("deltaR_j1_lnu", "deltaR(j1,lnu)",40, 0, 8);
//     m_deltaR_j_b = m_histoSvc->th1F("deltaR_j_b", "deltaR(j,b)",40, 0, 8);
//     m_ht = m_histoSvc->th1F("ht","ht of the event", 20,0,600);
// //     m_hfor_pretag = m_histoSvc->th1F("hfor_pretag", "hfor", 6, -1, 5);
// //     m_hfor_tag = m_histoSvc->th1F("hfor_tag", "hfor", 6, -1, 5);
//     m_lepmetphi = m_histoSvc->th1F("lepmetphi", "#Delta #phi of lepton and met", 33, 0, 3.25);

    h_TuDoCutFlow = new TH1F("TuDoCutFlow", "TuDoCutFlow", 14, 0, 14);
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(1,"INITIAL");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(2,"HFOR veto");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(3,"Preselection");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(4,"==1 Lepton 25 GeV");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(5,"==2 Jets 30 GeV");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(6,"==1 BTAG mv2c20@60");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(7,"MWT cut");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(8,"ptLep(DeltaPhi) +");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(9,"ptLep(DeltaPhi) -");
    h_TuDoCutFlow->GetXaxis()->SetBinLabel(10,"MET cut");
}

AnaTuDoSL::~AnaTuDoSL() {
}

void AnaTuDoSL::setIsData(bool isData){
    m_isData = isData;
}

void AnaTuDoSL::run(const Event &evt, double weight, const std::string &s) {
    
    h_TuDoCutFlow -> Fill(0);
    h_TuDoCutFlow -> Fill(1);

    // check channel
    if (m_electron && (evt.electron().size() != 1 || evt.muon().size() != 0))
        return;

    if (!m_electron && (evt.electron().size() != 0 || evt.muon().size() != 1))
        return;
    
    if (!(evt.passes("ejets") || evt.passes("mujets")))
        return;

    h_TuDoCutFlow -> Fill(2);
    
//     HistogramService *h = &m_hSvc;
    TLorentzVector lepton, leptonGeV;
    float lepcharge = -99;
    if (m_electron) {
        lepton = evt.electron()[0].mom();
//         lepcharge = evt.electron()[0].charge();
    } else {
        lepton = evt.muon()[0].mom();
//         lepcharge = evt.electron()[0].charge();
    }
    
    if(lepton.Pt() / 1000.0 < 25.0) // cut on exactly one lepton with 25 GeV
        return;
    
    h_TuDoCutFlow -> Fill(3);
    
    leptonGeV.SetPtEtaPhiE(lepton.Perp() / 1000.0, lepton.Eta(), lepton.Phi(), lepton.E()/1000.0);
//     m_mini->binning01 = -1;
//     m_mini->binning02 = abs(lepton.Eta());

    std::vector<TLorentzVector *> jets;
    std::vector<bool> jets_btagged;
    for (size_t it_jet = 0; it_jet < evt.jet().size(); ++it_jet) {
        jets.push_back(new TLorentzVector(0,0,0,0));
        jets[it_jet]->SetPtEtaPhiE(evt.jet()[it_jet].mom().Perp(), evt.jet()[it_jet].mom().Eta(), evt.jet()[it_jet].mom().Phi(), evt.jet()[it_jet].mom().E());
        // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BTagingxAODEDM
        // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks
        jets_btagged.push_back(evt.jet()[it_jet].btag_mv2c20_60());
    }
  
    // calculate neutrino px, py and pz
    std::vector<TLorentzVector*> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&lepton, evt.met().Perp(), evt.met().Phi(), true);   
    TLorentzVector neutrino(0,0,0,0), neutrinoGeV(0,0,0,0);
    if (vec_nu.size() > 0) {
        neutrino = *(vec_nu[0]);
        for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
        vec_nu.clear();
    }
  
    neutrinoGeV.SetPtEtaPhiE(neutrino.Perp() / 1000.0, neutrino.Eta(), neutrino.Phi(), neutrino.E()/1000.0);
    leptonGeV.SetPtEtaPhiE(lepton.Perp() / 1000.0, lepton.Eta(), lepton.Phi(), lepton.E()/1000.0);
  
    if(jets.size() != 2) // cut on exactly two jets
        return;
    
    h_TuDoCutFlow -> Fill(4);
    
    int ntags = 0;
    int index_taggedjet = -99;
    int index_untaggedjet = -99;
    int index_leadingjet = -99;
    int index_secondleadingjet = -99;
    for (size_t it_jet = 0; it_jet < jets.size(); ++it_jet) {
        if (jets_btagged.at(it_jet) == true) {
            ntags++;
            index_taggedjet = it_jet;
        } else {
            index_untaggedjet = it_jet;
        }
    }
    
    if(index_taggedjet < index_untaggedjet) {
        index_leadingjet = index_taggedjet;
        index_secondleadingjet = index_untaggedjet;
        
    } else {
        index_leadingjet = index_untaggedjet;
        index_secondleadingjet = index_taggedjet;
    }
  
    if(ntags != 1) // cut on exactly one tagged jet
        return;
    
    h_TuDoCutFlow -> Fill(5);

    if (index_taggedjet == index_untaggedjet) {
        std::cout << "Error, index_taggedjet == index_untaggedjet = " << index_taggedjet << " should never happen. Aborting!" << std::endl;
        exit(0);
    }
    
    TLorentzVector bjet, nonbjet;
    bjet.SetPtEtaPhiE(jets.at(index_taggedjet)->Pt() / 1000.0,jets.at(index_taggedjet)->Eta(),jets.at(index_taggedjet)->Phi(),jets.at(index_taggedjet)->E() / 1000.0);
    nonbjet.SetPtEtaPhiE(jets.at(index_untaggedjet)->Pt() / 1000.0,jets.at(index_untaggedjet)->Eta(),jets.at(index_untaggedjet)->Phi(),jets.at(index_untaggedjet)->E() / 1000.0);
//     bjet = *jets.at(index_taggedjet);
//     nonbjet = *jets.at(index_untaggedjet);
    
    TLorentzVector leadingjet, secondleadingjet;
    leadingjet.SetPtEtaPhiE(jets.at(index_leadingjet)->Pt() / 1000.0,jets.at(index_leadingjet)->Eta(),jets.at(index_leadingjet)->Phi(),jets.at(index_leadingjet)->E() / 1000.0);
    secondleadingjet.SetPtEtaPhiE(jets.at(index_secondleadingjet)->Pt() / 1000.0,jets.at(index_secondleadingjet)->Eta(),jets.at(index_secondleadingjet)->Phi(),jets.at(index_secondleadingjet)->E() / 1000.0);
//     leadingjet = *jets.at(index_leadingjet);
//     secondleadingjet = *jets.at(index_secondleadingjet);
        
    float mtw = KinematicUtils::transMass(lepton.Pt() / 1000.0, lepton.Phi(), evt.met().Perp() / 1000.0, evt.met().Phi());
    
    if(mtw < 50.0) // cut on transverse w-mass (in GeV!!!)
        return;
    
    h_TuDoCutFlow -> Fill(6);
  
    float deltaPhi_j1_lep = 0.0;
    deltaPhi_j1_lep = TuDoAtlas::delta_phi(leadingjet,lepton);
//     std::cout << lepton.Pt() / 1000.0 << " " << deltaPhi_j1_lep << " " << leadingjet.Phi() << " " << lepton.Phi() << std::endl;
    
    bool triangularcut = false;
    // leppt > (40 + (40/2.14 * (DeltaPhi - 3.14))) || leppt > (40 - (40/2.14 * (DeltaPhi + 3.14)))
    if (   (  lepton.Pt() / 1000.0 < (40+ ((40/(2.14))*(deltaPhi_j1_lep-3.14))) )
        || (  lepton.Pt() / 1000.0 < (40- ((40/(2.14))*(deltaPhi_j1_lep+3.14))) )) {
        triangularcut = true;
    }

    h_TuDoCutFlow -> Fill(7);
    h_TuDoCutFlow -> Fill(8);
    
    if(evt.met().Perp() / 1000.0 > 30.0) {
        h_TuDoCutFlow -> Fill(9);
    }

//     m_lep_DeltaPhi->Fill(deltaPhi_j1_lep,lepton.Pt() / 1000.0,weight,index_leadingjet);
    
    if(triangularcut) // triangular cut rejecting QCD
        return;
    
//     m_lep_DeltaPhiCut->Fill(deltaPhi_j1_lep,lepton.Pt() / 1000.0,weight,index_leadingjet);
    
    TLorentzVector lb = leptonGeV + bjet;
    TLorentzVector lnu = leptonGeV + neutrinoGeV;
    TLorentzVector lnub = lnu + bjet;
    TLorentzVector jb = bjet + nonbjet;
    TLorentzVector lnubj = lnub + nonbjet;
    
    double deltaR_jb = KinematicUtils::deltaR(bjet.Eta(), nonbjet.Eta(),bjet.Phi(), nonbjet.Phi());
    double deltaR_j1lnu = KinematicUtils::deltaR(leadingjet.Eta(), lnu.Eta(),leadingjet.Phi(),lnu.Phi());
    float ht = TuDoAtlas::ht2(bjet,nonbjet,leptonGeV,neutrinoGeV);
    float deltaR_l_lnub = KinematicUtils::deltaR(leptonGeV.Eta(), lnub.Eta(),leptonGeV.Phi(),lnub.Phi());
    float pt_lnu = lnu.Pt();
    float pt_lnub = lnub.Pt();
    float mnu = TuDoAtlas::m_nu(leptonGeV,neutrinoGeV);
    float costheta_lj_top = TuDoAtlas::Polarization(nonbjet,leptonGeV,lnub);
    float WHelicity = TuDoAtlas::WHelicity(lnub,lnu,leptonGeV);
    float aplanarity = TuDoAtlas::aplanarity3(bjet,nonbjet,leptonGeV,neutrinoGeV);
    float sphericity = TuDoAtlas::spherisity3(bjet,nonbjet,leptonGeV,neutrinoGeV);
    
//     InputArray[0] = ht;
//     InputArray[1] = mtw;
//     InputArray[2] = lb.M();
//     InputArray[3] = jb.M();
//     InputArray[4] = lnub.M();
//     InputArray[5] = std::fabs(nonbjet.Eta());
//     InputArray[6] = lnu.Eta();
//     InputArray[7] = lnub.Eta();
//     InputArray[8] = bjet.Eta();
//     InputArray[9]  = costheta_lj_top;
//     InputArray[10] = evt.met().Perp() / 1000.0;
//     InputArray[11] = pt_lnu;
//     InputArray[12] = pt_lnub;
//     InputArray[13] = deltaR_l_lnub;
//     
//     float NNoutput = NBExpert->nb_expert(InputArray);
//     float NNoutput_transformed = (NNoutput + 1.0)/2.0;
//     m_mini->binning03 = NNoutput_transformed;

    float deltaEta_l_b = TuDoAtlas::delta_eta(leptonGeV, bjet);
    float deltaPhi_l_b = TuDoAtlas::delta_phi(leptonGeV, bjet);
    float deltaR_l_b = TuDoAtlas::delta_r(leptonGeV, bjet);
    
//     m_leppt->Fill(lepton.Pt() / 1000.0, weight, index_leadingjet);
//     m_lepeta->Fill(lepton.Eta(), weight, index_leadingjet);
//     m_lepphi->Fill(lepton.Phi(), weight, index_leadingjet);
//     m_lepcharge->Fill(lepcharge, weight, index_leadingjet);
//     m_lepmetphi->Fill(abs(evt.met().Phi() - lepton.Phi()), weight,index_leadingjet);
//     
//     m_nupt->Fill(neutrino.Pt() / 1000.0, weight, index_leadingjet);
//     m_nueta->Fill(neutrino.Eta(), weight, index_leadingjet);
//     m_nuphi->Fill(neutrino.Phi(), weight, index_leadingjet);
// 
//     m_jpt_1->Fill(jets.at(index_leadingjet)->Pt() / 1000.0, weight, index_leadingjet);
//     m_jeta_1->Fill(jets.at(index_leadingjet)->Eta(), weight, index_leadingjet);
//     m_jphi_1->Fill(jets.at(index_leadingjet)->Phi(), weight, index_leadingjet);
//     m_jpt_2->Fill(jets.at(index_secondleadingjet)->Pt() / 1000.0, weight, index_leadingjet);
//     m_jeta_2->Fill(jets.at(index_secondleadingjet)->Eta(), weight, index_leadingjet);
//     m_jphi_2->Fill(jets.at(index_secondleadingjet)->Phi(), weight, index_leadingjet);
//     
//     m_bpt->Fill(jets.at(index_taggedjet)->Pt() / 1000.0, weight, index_leadingjet);
//     m_beta->Fill(jets.at(index_taggedjet)->Eta(), weight, index_leadingjet);
//     m_bphi->Fill(jets.at(index_taggedjet)->Phi(), weight, index_leadingjet);
//     m_lightpt->Fill(jets.at(index_untaggedjet)->Pt() / 1000.0, weight, index_leadingjet);
//     m_lighteta->Fill(jets.at(index_untaggedjet)->Eta(), weight, index_leadingjet);
//     m_lightphi->Fill(jets.at(index_untaggedjet)->Phi(), weight, index_leadingjet);
// 
//     m_met->Fill(evt.met().Perp() / 1000.0, weight,index_leadingjet);
//     m_metx->Fill(evt.met().Perp()*std::cos(evt.met().Phi()) / 1000.0, weight,index_leadingjet);
//     m_mety->Fill(evt.met().Perp()*std::sin(evt.met().Phi()) / 1000.0, weight,index_leadingjet);
//     m_metphi->Fill(evt.met().Phi(), weight,index_leadingjet);
//     m_mtw->Fill(mtw, weight,index_leadingjet);
// 
//     m_deltaR_j1_lnu->Fill(deltaR_j1lnu, weight, index_leadingjet);
//     m_ht->Fill(ht,weight, index_leadingjet);
//     m_deltaR_j_b->Fill(deltaR_jb, weight, index_leadingjet);
//     m_mlnub->Fill(lnub.M(), weight, index_leadingjet);
//     m_mb->Fill(bjet.M(), weight, index_leadingjet);
//     m_mjb->Fill(jb.M(), weight, index_leadingjet);
//     m_etalnu->Fill(lnu.Eta(),weight, index_leadingjet);
//     m_etalnub->Fill(lnub.Eta(),weight, index_leadingjet);
//     m_mlb->Fill(lb.M(), weight, index_leadingjet);
//     m_deltaR_l_lnub->Fill(deltaR_l_lnub, weight, index_leadingjet);
//     m_pt_lnu->Fill(pt_lnu, weight, index_leadingjet);
//     m_pt_lnub->Fill(pt_lnub, weight, index_leadingjet);
//     m_mnu->Fill(mnu, weight, index_leadingjet);
//     m_costheta_lj_top->Fill(costheta_lj_top, weight, index_leadingjet);
//     m_WHelicity->Fill(WHelicity, weight, index_leadingjet);
//     m_aplanarity->Fill(aplanarity, weight, index_leadingjet);
//     m_sphericity->Fill(sphericity, weight, index_leadingjet);
//     m_mlnubj->Fill(lnubj.M(), weight, index_leadingjet);
//     m_etaj->Fill(std::fabs(nonbjet.Eta()), weight, index_leadingjet);
//     m_deltaEta_l_b->Fill(deltaEta_l_b, weight, index_leadingjet);
//     m_deltaPhi_l_b->Fill(deltaPhi_l_b, weight, index_leadingjet);
//     m_deltaR_l_b->Fill(deltaR_l_b, weight, index_leadingjet);
}

void AnaTuDoSL::terminate() {
    
    h_TuDoCutFlow->SetDirectory(0);
    h_TuDoCutFlow->Write();
    
    TH1D* cutflow;
    if(!m_isData) {
        if(m_electron) {
            cutflow = (TH1D*) m_mini->m_file->Get("ejets/cutflow_mc_pu_zvtx")->Clone("cutflow_mc_pu_zvtx");
        } else {
            cutflow = (TH1D*) m_mini->m_file->Get("mujets/cutflow_mc_pu_zvtx")->Clone("cutflow_mc_pu_zvtx");
        }
    } else {
        if(m_electron) {
            cutflow = (TH1D*) m_mini->m_file->Get("ejets/cutflow")->Clone("cutflow");
        } else {
            cutflow = (TH1D*) m_mini->m_file->Get("mujets/cutflow")->Clone("cutflow");
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
