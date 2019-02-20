#include "TopNtupleAnalysis/UserFunktions.h"
// #include <cmath>

namespace TuDoAtlas {
    
  double angle(TLorentzVector v1, TLorentzVector v2){
    return ROOT::Math::VectorUtil::Angle<TLorentzVector, TLorentzVector>(v1,v2);
  }

  double cos_theta(TLorentzVector v1, TLorentzVector v2){
    return ROOT::Math::VectorUtil::CosTheta<TLorentzVector, TLorentzVector>(v1,v2);
  }

  double delta_phi(TLorentzVector v1, TLorentzVector v2){
    return ROOT::Math::VectorUtil::DeltaPhi<TLorentzVector, TLorentzVector>(v1,v2);
  }
  double invariant_mass(TLorentzVector v1, TLorentzVector v2){
    return ROOT::Math::VectorUtil::InvariantMass<TLorentzVector, TLorentzVector>(v1,v2);
  }

  double delta_r(TLorentzVector v1, TLorentzVector v2){
    return ROOT::Math::VectorUtil::DeltaR<TLorentzVector, TLorentzVector>(v1,v2);
  }

  double delta_eta(TLorentzVector v1, TLorentzVector v2){
    return fabs(v1.Eta()-v2.Eta());
  }

  double WHelicity(TLorentzVector top4, TLorentzVector W4, TLorentzVector l4)
  {
    double costheta;

    // boost into rest frame
    l4.Boost( -W4.BoostVector() );
    W4.Boost( -top4.BoostVector() );

    TVector3 W3rest(W4.Px(), W4.Py(), W4.Pz());
    TVector3 l3rest(l4.Px(), l4.Py(), l4.Pz());

    costheta = W3rest * l3rest / (W3rest.Mag() * l3rest.Mag());
    return costheta;
  }
  double Polarization(TLorentzVector lqjet, TLorentzVector lepton, TLorentzVector top)
  {
    double pol=-999;

    // Polarisation of Top quark
    // cos Theta of Angle between Lepton fom W_top and Light Quark Jet in the
    // Top Quark rest Frame
    TLorentzVector lqj = lqjet;
    TLorentzVector l4 = lepton;

    TVector3 top_boost = top.BoostVector();
    l4.Boost(-top_boost);
    lqj.Boost(-top_boost);
    pol = l4.Vect()*lqj.Vect()/(l4.Vect().Mag()*lqj.Vect().Mag() );
    return pol;
  }

  double mT_W(TLorentzVector l4, TLorentzVector nu4)
  {
    double mt = sqrt(2*(l4.Pt()*nu4.Pt() - l4.Px()*nu4.Px() - l4.Py()*nu4.Py()));
    return mt;
  }

  double m_nu(TLorentzVector l4, TLorentzVector nu4)
  {

    double mt = sqrt(2*l4.Pt()*nu4.Pt()*(1-cos(delta_phi(l4,nu4))));
    double costheta_lep=l4.Pz()/sqrt(l4.Pt()*l4.Pt()+l4.Pz()*l4.Pt());
    double Wet=sqrt(mt*mt+80.376*80.376*costheta_lep*costheta_lep);
    double Wet_corr=Wet-l4.Pt();
    double TMpseudo_nu_METbase_lep_M=Wet_corr*Wet_corr-(nu4.Px()*nu4.Px()+nu4.Py()*nu4.Py());
    double c=1;
    if(TMpseudo_nu_METbase_lep_M<0.) {c=-1.;}
    return c*sqrt(fabs(TMpseudo_nu_METbase_lep_M));
  }

  //triangular cut
  bool PassSlideMeT(double dPhiAnchor, double MeTAnchor, double eldPhiMet, double met){
    bool pass = true;
    if ( (eldPhiMet) < ((-dPhiAnchor/MeTAnchor)*met + dPhiAnchor) ) pass = false;
    return pass;
  }
  
  //aplanarity + spherisity
  std::pair<double, double> calc_apl_shere(std::vector<TLorentzVector> jets, std::vector<TLorentzVector> leps, std::vector<TLorentzVector> mets){

    double denom =0;
    double numXX =0;
    double numYY =0;
    double numZZ =0;
    double numXY =0;
    double numXZ =0;
    double numYZ =0;
   
    for(unsigned int j=0; j<jets.size();++j){
      denom += jets[j].Vect().Mag2();
      numXX += jets[j].X() * jets[j].X();
      numYY += jets[j].Y() * jets[j].Y();
      numZZ += jets[j].Z() * jets[j].Z();
      numXY += jets[j].X() * jets[j].Y();
      numXZ += jets[j].X() * jets[j].Z();
      numYZ += jets[j].Y() * jets[j].Z();
    }

    for(unsigned int j=0; j<leps.size();++j){
      denom += leps[j].Vect().Mag2();
      numXX += leps[j].X() * leps[j].X();
      numYY += leps[j].Y() * leps[j].Y();
      numZZ += leps[j].Z() * leps[j].Z();
      numXY += leps[j].X() * leps[j].Y();
      numXZ += leps[j].X() * leps[j].Z();
      numYZ += leps[j].Y() * leps[j].Z();
    }

    for(unsigned int j=0; j<mets.size();++j){
      denom += mets[j].Vect().Mag2();
      numXX += mets[j].X() * mets[j].X();
      numYY += mets[j].Y() * mets[j].Y();
      numZZ += mets[j].Z() * mets[j].Z();
      numXY += mets[j].X() * mets[j].Y();
      numXZ += mets[j].X() * mets[j].Z();
      numYZ += mets[j].Y() * mets[j].Z();
    }

    TMatrix spherTens(3,3);
    spherTens(0,0) = numXX/denom;
    spherTens(1,1) = numYY/denom;
    spherTens(2,2) = numZZ/denom;
    spherTens(0,1) = numXY/denom;
    spherTens(1,0) = numXY/denom;
    spherTens(0,2) = numXZ/denom;
    spherTens(2,0) = numXZ/denom;
    spherTens(1,2) = numYZ/denom;
    spherTens(2,1) = numYZ/denom;

    TVector eigenval(3);
    spherTens.EigenVectors(eigenval);

    double sphericity=0;
    double aplanarity=0;

    if (eigenval(0)<=eigenval(1) && eigenval(1)<=eigenval(2)){
      sphericity = 3.0*(eigenval(0)+eigenval(1))/2.0;
      aplanarity = 3.0*eigenval(0)/2.0;
    }
    if (eigenval(0)<=eigenval(2) && eigenval(2)<=eigenval(1)){
      sphericity = 3.0*(eigenval(0)+eigenval(2))/2.0;
      aplanarity = 3.0*eigenval(0)/2.0;
    }
    if (eigenval(1)<=eigenval(0) && eigenval(0)<=eigenval(2)){
      sphericity = 3.0*(eigenval(1)+eigenval(0))/2.0;
      aplanarity = 3.0*eigenval(1)/2.0;
    }
    if (eigenval(1)<=eigenval(2) && eigenval(2)<=eigenval(0)){
      sphericity = 3.0*(eigenval(1)+eigenval(2))/2.0;
      aplanarity = 3.0*eigenval(1)/2.0;
    }
    if (eigenval(2)<=eigenval(0) && eigenval(0)<=eigenval(1)){
      sphericity = 3.0*(eigenval(2)+eigenval(0))/2.0;
      aplanarity = 3.0*eigenval(2)/2.0;
    }
    if (eigenval(2)<=eigenval(1) && eigenval(1)<=eigenval(0)){
      sphericity = 3.0*(eigenval(2)+eigenval(1))/2.0;
      aplanarity = 3.0*eigenval(2)/2.0;
    }

    return std::make_pair(aplanarity, sphericity);
  }

  double aplanarity(std::vector<TLorentzVector> jets, std::vector<TLorentzVector> leps, std::vector<TLorentzVector> mets){
    return calc_apl_shere(jets, leps, mets).first;
  }

  double aplanarity2(std::vector<TLorentzVector> jets, TLorentzVector lep, TLorentzVector met){
    std::vector<TLorentzVector> leps,mets;
    leps.push_back(lep);
    mets.push_back(met);
    return calc_apl_shere(jets, leps, mets).first;
  }
  
  double aplanarity3(TLorentzVector bjet, TLorentzVector nonbjet, TLorentzVector lep, TLorentzVector met){
    std::vector<TLorentzVector> jets,leps,mets;
    leps.push_back(lep);
    mets.push_back(met);
    jets.push_back(bjet);
    jets.push_back(nonbjet);
    return calc_apl_shere(jets, leps, mets).first;
  }

  double spherisity(std::vector<TLorentzVector> jets, std::vector<TLorentzVector> leps, std::vector<TLorentzVector> mets){
    return calc_apl_shere(jets, leps, mets).second;
  }

  double spherisity2(std::vector<TLorentzVector> jets, TLorentzVector lep, TLorentzVector met){
    std::vector<TLorentzVector> leps,mets;
    leps.push_back(lep);
    mets.push_back(met);

    return calc_apl_shere(jets, leps, mets).second;
  }
  
  double spherisity3(TLorentzVector bjet, TLorentzVector nonbjet, TLorentzVector lep, TLorentzVector met){
    std::vector<TLorentzVector> jets,leps,mets;
    leps.push_back(lep);
    mets.push_back(met);
    jets.push_back(bjet);
    jets.push_back(nonbjet);
    return calc_apl_shere(jets, leps, mets).second;
  }
  
  double jetprobRND()
  {
    double weight=-999;
    TRandom3 rnd(0);
    double random = rnd.Rndm();
    if(-log(random)>4.26)
     weight=-log(random);
    return weight;
  }
  
  double ht(std::vector<TLorentzVector>& leps, std::vector<TLorentzVector>& jets, std::vector<TLorentzVector>& met){

    double _ht = 0; 

    for(unsigned int i=0; i<leps.size(); ++i)
      _ht += leps[i].Pt();

    for(unsigned int i=0; i<jets.size(); ++i)
      _ht += jets[i].Pt();
      
    for(unsigned int i=0; i<met.size(); ++i)
      _ht += met[i].Pt();

    return _ht;
  }
  
  double ht2(TLorentzVector bjet, TLorentzVector nonbjet, TLorentzVector lep, TLorentzVector met){
    std::vector<TLorentzVector> jets,leps,mets;
    leps.push_back(lep);
    mets.push_back(met);
    jets.push_back(bjet);
    jets.push_back(nonbjet);
    return ht(jets, leps, mets);
  }
    

//   //find truth particle closest to 4 vector
//   int matchClosestTruthParticle(MCEvent* mcevent,TLorentzVector v4,int status, double dr)
//   {
//     int pdg_id=-999;
//     for(unsigned int i=0;i<mcevent->particles.size();i++){
//       if(fabs(mcevent->particles[i].status==3)){
//         if(mcevent->particles[i].v4.DeltaR(v4)<dr){
//           pdg_id=mcevent->particles[i].pdg_id;
//           dr=mcevent->particles[i].v4.DeltaR(v4);
//         }
//       }
//     }
//     return pdg_id;
//   }

//   int MuonOrElectronInJet(MCEvent* mcevent,TLorentzVector v4,double dr)
//   {
//     int pdg_id=-999;
//     for(unsigned int i=0;i<mcevent->particles.size();i++){
//       if(abs(mcevent->particles[i].pdg_id==11) || abs(mcevent->particles[i].pdg_id==13)){
//         if(mcevent->particles[i].v4.DeltaR(v4)<dr){
//           return mcevent->particles[i].pdg_id;
//         }
//       }
//     }
//     return pdg_id;
// 
//   }

//   //find truth particle closest to 4 vector
//   double matchClosestDr(MCEvent* mcevent,TLorentzVector v4,int status, double dr)
//   {
//     for(unsigned int i=0;i<mcevent->particles.size();i++){
//       if(fabs(mcevent->particles[i].status==3)){
//         if(mcevent->particles[i].v4.DeltaR(v4)<dr){
//           dr=mcevent->particles[i].v4.DeltaR(v4);
//         }
//       }
//     }
//     return dr;
//   }




//   //check for duplicate events
//   std::map<unsigned int, std::set<unsigned int> > _duplicate_check;
//   bool _total_events_resetted = false;
//   int _ev = 0;
//   bool check_duplicate(SampleInfo& sa_info, EventInfo& ev_info, EventLog& ev_log){
// 
//     if(!_total_events_resetted){
//       sa_info.total_events -= ev_log.weight_sum;
//       _total_events_resetted = true;
//     }
// 
//     unsigned int run_num = ev_info.run_number;
//     unsigned int event_num = ev_info.event_number;
//     bool dup = false;
//     std::cout << "\n" << ++_ev << "\t" << run_num << "\t" << event_num << "\t\t";
// 
//     std::map<unsigned int, std::set<unsigned int> >::iterator it = _duplicate_check.find(run_num);
//     if(it == _duplicate_check.end()){
//       std::cout << " inserting run :" << run_num << "\t";
//       _duplicate_check.insert(std::make_pair(run_num, std::set<unsigned int>()));
//       _duplicate_check[run_num].insert(event_num);
//     }
//     else{
//       std::set<unsigned int>::iterator it2 = it->second.find(event_num);
//       if(it2 == it->second.end()){
//         std::cout << " inserting ev :" << event_num << "\t\t";
//         _duplicate_check[run_num].insert(event_num);
//       }
//       else{
//         dup = true;
//         std::cout << "DUP!!!!!!!!!!!!!!!!!!!!!!";
//       }
//     }
//     //std::cout << "dup = " << dup << " \t Run number= "<< run_num << " \t Event number= " << event_num << "\n";
//     return dup;
//   }
//   bool clean_up(){
//     _duplicate_check.clear();
//     return true;
//   }
// 

  
//   int getNbtags(std::vector<Particle>& jets, std::string btagger, double wp){
// 
//     int nbtags = 0;
// 
//     for(unsigned int i=0; i<jets.size();++i){
//       double jeteta = fabs(jets[i].v4.Eta());
//       if (jets[i].GetDouble("EMScale_Eta") != -999)
//      jeteta =  fabs(jets[i].GetDouble("EMScale_Eta")+jets[i].GetDouble("EMJES_EtaCorr"));
//       
//       if ( jeteta < 2.5 ) {
//      double weight = jets[i].GetDouble(btagger);
//              if (weight > wp)
//        nbtags++;
// 
//      //std::cout<<jets[i].v4.Eta()<<" "<<weight<<" "<<nbtags<<std::endl;
//       }
//     }
// 
//     return nbtags;
//   }
  
//   double getCharge(std::vector<Particle>& ele, std::vector<Particle>& muo){
// 
//     double charge = 0;
// 
//     if (ele.size() == 1) 
//       charge = ele[0].GetDouble("charge");
//     else if (muo.size() == 1) 
//       charge = muo[0].GetDouble("charge");
//     else
//       std::cout<<"getCharge: too many leptons in the event!"<<std::endl;
// 
//     return charge;
//   }
  
//   bool hasTagableJet(std::vector<Particle>& jets){
// 
//     bool istagable = 0;
// 
//     for(unsigned int i=0; i<jets.size();++i){
//       double jeteta = fabs(jets[i].v4.Eta());
//       if (jets[i].GetDouble("EMScale_Eta") != -999)
//              jeteta =  fabs(jets[i].GetDouble("EMScale_Eta")+jets[i].GetDouble("EMJES_EtaCorr"));
//       
//       if ( jeteta < 2.5 ) {
//         istagable = 1;
//       }
//     }
// 
//     return istagable;
//   }

//   bool vetoBtag(std::vector<Particle>& jets,std::string taggername, double wp){
// 
//     bool hasBtag = 0;
// 
//     for(unsigned int i=0; i<jets.size();++i){
//       double jeteta = fabs(jets[i].v4.Eta());
//       if (jets[i].GetDouble("EMScale_Eta") != -999)
//              jeteta =  fabs(jets[i].GetDouble("EMScale_Eta")+jets[i].GetDouble("EMJES_EtaCorr"));
//       
//       if ( jeteta < 2.5 && jets[i].GetDouble(taggername) > wp) {
//         hasBtag = 1;
//       }
//     }
// 
//     return hasBtag; // true, if there is a btag
//   }

//   int getTruthFlavour(std::vector<Particle>& jets){
// 
//     int max_truth_flavour = 0;
// 
//     for(unsigned int i=0; i<jets.size();++i){
//       int truth_flavour = jets[i].GetInt("TruthFlavour");
//       
//       if ((truth_flavour > max_truth_flavour) && (truth_flavour < 6)) {
//         max_truth_flavour = truth_flavour;
//       }
//     }
// 
//     return max_truth_flavour; 
//   }

//   double getMomentumFraction_lightjet(MCEvent* mcevent)
//   {
//     double x = 0.0;
//     if((fabs(mcevent->GetInt("pdf_id1"))!= 5) && (fabs(mcevent->GetInt("pdf_id1"))!= 21)) {
//       x = mcevent->GetDouble("pdf_x1");
//     }
//     else {
//       x = mcevent->GetDouble("pdf_x2");
//     }
// 
//     return x;
//   }

//   bool TruthJetKinCuts(std::vector<Particle>& truth_jets, double pt_cut, double eta_cut, int jetbin){
// 
//     int nSelected = 0;
//     bool selected = false;
// 
//     for(unsigned int i=0; i<truth_jets.size();++i){
//       double jeteta = fabs(truth_jets[i].v4.Eta());
//       double jetpt = truth_jets[i].v4.Pt();
//        
//       if ( (jeteta < eta_cut) && (jetpt > pt_cut)) {
//         nSelected++;
//       }
//     }
//     if(jetbin < 0 ){
//       if(nSelected >= fabs(jetbin)) {
//         selected = true;
//       }
//     }
//     else {
//       if(nSelected == jetbin) {
//         selected = true;
//       }
//     }
// 
//     return selected;
//   }

//   TLorentzVector isLightTruthJet(std::vector<Particle>& truth_jets){
// 
//     TLorentzVector v4 = truth_jets[0].v4;
// 
//     for(unsigned int i=0; i<truth_jets.size();++i){
//       if(truth_jets[i].GetInt("isBHadron") == 0) {
//         v4 = truth_jets[i].v4;
//       }
//     }
//  
//     return v4;
//   }

//   bool matchSelectedLightQuark(std::vector<Particle>& truth_jets, TLorentzVector v4){
//     bool matchSelected = false;
//     double dr = 999.0;
//     unsigned int index = 99; 
//     
//     for(unsigned int i=0; i<truth_jets.size();++i){
//       if(truth_jets[i].v4.DeltaR(v4)<dr){
//         index=i;
//         dr=truth_jets[i].v4.DeltaR(v4);
//       }
//     }
//     if ( (fabs(truth_jets[index].v4.Eta()) < 4.5) && (truth_jets[index].v4.Pt() > 30.)) {
//       matchSelected = true;
//     }
//     return matchSelected;
//   }
  //analysis dependent functions________________________________________________________

//   #include "vbf/include/vbfFunktions.h"
//   #include "zprime/include/ZprimeFunctions.h"  
    
}