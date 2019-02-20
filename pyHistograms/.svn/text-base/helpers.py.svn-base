
import ROOT
from array import array
import math
import sys
sys.path.append('2HDM')
#import T2HDM


def initBinds():
    bind_neutrino = '''

#include "TopNtupleAnalysis/TtresNeutrinoBuilder.h"
#include "TopNtupleAnalysis/TtresChi2.h"
#include "TopNtupleAnalysis/WeakCorrScaleFactorParam.h"
#include "TopNtupleAnalysis/MMUtils.h"
#include "TopNtupleAnalysis/Event.h"
#include "TopNtupleAnalysis/Electron.h"
#include "TopNtupleAnalysis/Muon.h"
#include "TopNtupleAnalysis/Jet.h"
#include "TopNtupleAnalysis/LargeJet.h"
#include "TopNtupleAnalysis/MObject.h"

/*
#include "xAODBTaggingEfficiency/BTaggingEfficiencyTool.h"
#include "CalibrationDataInterface/CalibrationDataVariables.h"
#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"
#include "PATInterfaces/SystematicRegistry.h"
*/

#include "TLorentzVector.h"
#include <vector>
#include <cstdlib>
#include <sstream>

#include "TopNtupleAnalysis/EFTLib.h"

class WrapperExtras {
  public:
    TtresNeutrinoBuilder m_neutrinoBuilder;
    TtresChi2 m_chi2;
    WeakCorr::WeakCorrScaleFactorParam m_ewkTool;
    //PMGCorrsAndSysts m_pmg;
    MMUtils m_mm_b0_boosted_e;
    MMUtils m_mm_b1_boosted_e;
    MMUtils m_mm_b0_res_e;
    MMUtils m_mm_b1_res_e;
    MMUtils m_mm_b0_boosted_mu;
    MMUtils m_mm_b1_boosted_mu;
    MMUtils m_mm_b0_res_mu;
    MMUtils m_mm_b1_res_mu;

    //BTaggingEfficiencyTool m_btageff;

    WrapperExtras()
                    : m_neutrinoBuilder("MeV"), m_chi2("MeV"), m_ewkTool("../share/EWcorr_param.root"),
                      m_mm_b0_boosted_e(0, "../scripts/QCDestimation/RATES_2015/resolved_e_eff_ttbar.root", // real2015
		                           "../scripts/QCDestimation/RATES_2015/resolved_e_btag0_fake.root",       // fake2015
					   "../scripts/QCDestimation/RATES_2016/resolved_e_eff_ttbar.root",        // real2016
					   "../scripts/QCDestimation/RATES_2016/resolved_e_btag0_fake.root"),      // fake2016
                      m_mm_b1_boosted_e(1, "../scripts/QCDestimation/RATES_2015/resolved_e_eff_ttbar.root", // real2015
		                           "../scripts/QCDestimation/RATES_2015/resolved_e_btag1_fake.root",       // fake2015
					   "../scripts/QCDestimation/RATES_2016/resolved_e_eff_ttbar.root",   // real2016
					   "../scripts/QCDestimation/RATES_2016/resolved_e_btag1_fake.root"), // fake2016

                      m_mm_b0_res_e(0, "../scripts/QCDestimation/RATES_2015/resolved_e_eff_ttbar.root", // real2015
		                           "../scripts/QCDestimation/RATES_2015/resolved_e_btag0_fake.root",       // fake2015
					   "../scripts/QCDestimation/RATES_2016/resolved_e_eff_ttbar.root",        // real2016
					   "../scripts/QCDestimation/RATES_2016/resolved_e_btag0_fake.root"),      // fake2016
                      m_mm_b1_res_e(1, "../scripts/QCDestimation/RATES_2015/resolved_e_eff_ttbar.root", // real2015
		                           "../scripts/QCDestimation/RATES_2015/resolved_e_btag1_fake.root",       // fake2015
					   "../scripts/QCDestimation/RATES_2016/resolved_e_eff_ttbar.root",   // real2016
					   "../scripts/QCDestimation/RATES_2016/resolved_e_btag1_fake.root"), // fake2016
                      m_mm_b0_boosted_mu(0, "../scripts/QCDestimation/RATES_2015/resolved_mu_eff_ttbar.root",// real2015
		                            "../scripts/QCDestimation/RATES_2015/resolved_mu_btag0_fake.root",      // fake2015
					    "../scripts/QCDestimation/RATES_2016/resolved_mu_eff_ttbar.root",        // real2016
					    "../scripts/QCDestimation/RATES_2016/resolved_mu_btag0_fake.root"), // fake2016
                      m_mm_b1_boosted_mu(1, "../scripts/QCDestimation/RATES_2015/resolved_mu_eff_ttbar.root",// real2015
		                            "../scripts/QCDestimation/RATES_2015/resolved_mu_btag1_fake.root",// fake2015
					    "../scripts/QCDestimation/RATES_2016/resolved_mu_eff_ttbar.root",       // real2016
					    "../scripts/QCDestimation/RATES_2016/resolved_mu_btag1_fake.root"), // fake2016

                      m_mm_b0_res_mu(0, "../scripts/QCDestimation/RATES_2015/resolved_mu_eff_ttbar.root",// real2015
		                            "../scripts/QCDestimation/RATES_2015/resolved_mu_btag0_fake.root",      // fake2015
					    "../scripts/QCDestimation/RATES_2016/resolved_mu_eff_ttbar.root",        // real2016
					    "../scripts/QCDestimation/RATES_2016/resolved_mu_btag0_fake.root"), // fake2016
                      m_mm_b1_res_mu(1, "../scripts/QCDestimation/RATES_2015/resolved_mu_eff_ttbar.root",// real2015
		                            "../scripts/QCDestimation/RATES_2015/resolved_mu_btag1_fake.root",// fake2015
					    "../scripts/QCDestimation/RATES_2016/resolved_mu_eff_ttbar.root",       // real2016
					    "../scripts/QCDestimation/RATES_2016/resolved_mu_btag1_fake.root") // fake2016
                      //m_btageff("BTaggingEfficiencyTool")
    {
      m_chi2.Init(TtresChi2::DATA2015_MC15C);
      /*
      if (!m_btageff.setProperty("TaggerName", "MV2c10")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("OperatingPoint", "FixedCutBEff_70")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("JetAuthor", "AntiKt2PV0TrackJets")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("EfficiencyFileName", "../share/2016-20_7-13TeV-MC15-CDI-June27_v1.root")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("ScaleFactorFileName", "../share/2016-20_7-13TeV-MC15-CDI-June27_v1.root")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("ScaleFactorBCalibration", "default")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("ScaleFactorCCalibration", "default")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("ScaleFactorTCalibration", "default")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("ScaleFactorLightCalibration", "default")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("ExcludeFromEigenVectorTreatment", "")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("EfficiencyBCalibrations", "410000;410004;410006;410187")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("EfficiencyCCalibrations", "410000;410004;410006;410187")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("EfficiencyTCalibrations", "410000;410004;410006;410187")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.setProperty("EfficiencyLightCalibrations", "410000;410004;410006;410187")) std::cout << "Failed to set b-tagging property." << std::endl;
      if (!m_btageff.initialize()) std::cout << "Failed to initialize b-tagging eff." << std::endl;
      CP::SystematicSet s = m_btageff.affectingSystematics();
      //for (CP::SystematicSet::const_iterator it = s.begin(); it != s.end(); it++) {
      //  std::cout << "Syst. " << it->name() << std::endl;
      //}
      */
    }
    TLorentzVector getNu(TLorentzVector l, double met, double met_phi) {
      std::vector<TLorentzVector *> vec_nu = m_neutrinoBuilder.candidatesFromWMass_Rotation(&l, met, met_phi, true);
      TLorentzVector nu;
      if (vec_nu.size() > 0) {
        nu = *(vec_nu[0]);
        for (size_t z = 0; z < vec_nu.size(); ++z) delete vec_nu[z];
        vec_nu.clear();
      }
      return nu;
    }

    std::map<std::string, double> getMtt(TLorentzVector lep, std::vector<TLorentzVector> jets, std::vector<bool> btag, TLorentzVector met) {
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
      for (size_t z = 0; z < jets.size(); ++z) {
        vjets.push_back(new TLorentzVector(0,0,0,0));
        vjets[z]->SetPtEtaPhiE(jets[z].Perp(), jets[z].Eta(), jets[z].Phi(), jets[z].E());
        // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BTagingxAODEDM
        // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks
        vjets_btagged.push_back(btag[z]);
      }
      bool status = m_chi2.findMinChiSquare(&lep, &vjets, &vjets_btagged, &met, igj3, igj4, igb3, igb4, ign1, chi2ming1, chi2ming1H, chi2ming1L); 
    
      float chi2Value = 1000000; // log10(1000000) = 6
      float mwh = -1;
      float mtl = -1;
      float mth = -1;
      float mtt = -1;

      if (status) {
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
      std::map<std::string, double> ret;
      ret["mtt"] = mtt;
      ret["mtl"] = mtl;
      ret["mth"] = mth;
      ret["mwh"] = mwh;
      ret["chi2"] = chi2Value;
      return ret;
    }
    double getEWK(TLorentzVector top, TLorentzVector topbar, int initial_type, int var = 0) {
      WeakCorr::ScaleFactor sf; 
      float t_pt = top.Perp();
      float t_eta = top.Eta();
      float t_phi = top.Phi();
      float t_m = top.M();

      float tbar_pt = topbar.Perp();
      float tbar_eta = topbar.Eta();
      float tbar_phi = topbar.Phi();
      float tbar_m = topbar.M();
      sf = m_ewkTool.getScaleFactor(t_pt, t_eta, t_phi, t_m, tbar_pt, tbar_eta, tbar_phi, tbar_m, initial_type);
      if (var == 1) {
        return sf.up;
      } else if (var == 2) {
        return sf.down;
      }
      return sf.nominal;
    }
	/*
    float getBtaggingSF(std::vector<TLorentzVector> jet, std::vector<int> flavour, std::vector<float> mv2c10, std::vector<bool> vetoSyst, std::string syst = "") {
      m_btageff.setMapIndex("B", 0);
      m_btageff.setMapIndex("C", 0);
      m_btageff.setMapIndex("T", 0);
      m_btageff.setMapIndex("Light", 0);

      float btagsf = 1.0;
      for (int k = 0; k < jet.size(); ++k) {
	    float sf = 1;
            Analysis::CalibrationDataVariables v;
            v.jetPt = jet[k].Perp();
            v.jetEta = jet[k].Eta();
            v.jetTagWeight = mv2c10[k];
            v.jetAuthor = "AntiKt2PV0TrackJets";
            if (vetoSyst[k]) {
              if (!m_btageff.applySystematicVariation(CP::SystematicSet())) {
                std::cout << "Could not apply the sys variation " << "" << std::endl;
                exit(-1);
              }
            } else {
              if (!m_btageff.applySystematicVariation(CP::SystematicSet(syst))) {
                std::cout << "Could not apply the sys variation " << syst << std::endl;
                exit(-1);
              }
            }

            if (mv2c10[k] > 0.6455) {
              if (!m_btageff.getScaleFactor(flavour[k], v, sf)) {
                std::cout << "Could not get btag SF for " << syst << ", jet pt, eta, weight = " << v.jetPt << ", " << v.jetEta << ", " << v.jetTagWeight << std::endl;
                exit(-1);
              }
            } else {
              if (!m_btageff.getInefficiencyScaleFactor(flavour[k], v, sf)) {
                std::cout << "Could not get btag ineff. SF for " << syst << ", jet pt, eta, weight = " << v.jetPt << ", " << v.jetEta << ", " << v.jetTagWeight << std::endl;                exit(-1);
              }
            }
            btagsf *= sf;
      }
      return btagsf;
    }
    */
    double getQCDWeight(int btags, int boosted, TLorentzVector met, TLorentzVector lep, int isTight, std::vector<TLorentzVector> jet, float sd0, int isElectron, int muonTrigger, float topoetcone20, int runNumber) {
      Event e;
      e.met(met.Px(), met.Py());
      if (isElectron) {
        e.electron().push_back(Electron(lep));
        e.electron()[0].setTightPP(isTight);
        e.electron()[0].topoetcone20() = topoetcone20;
        e.electron()[0].sd0() = sd0;
      } else {
        e.muon().push_back(Muon(lep));
    	e.muon()[0].sd0() = sd0;
    	e.muon()[0].setTight(isTight);
    	e.muon()[0].HLT_mu20_iloose_L1MU15() = muonTrigger;
    	e.muon()[0].HLT_mu50() = muonTrigger;
    	e.muon()[0].topoetcone20() = topoetcone20;
      }
      for (int k = 0; k < jet.size(); ++k) {
        e.jet().push_back(Jet(jet[k]));
      }
      double w = 0;
      if (!btags && boosted && isElectron) {
        w = m_mm_b0_boosted_e.getMMweights(e, 0, isElectron, boosted, runNumber);
      } else if (btags && boosted && isElectron) {
        w = m_mm_b1_boosted_e.getMMweights(e, 0, isElectron, boosted, runNumber);
      } else if (!btags && !boosted && isElectron) {
        w = m_mm_b0_res_e.getMMweights(e, 0, isElectron, boosted, runNumber);
      } else if (btags && !boosted && isElectron) {
        w = m_mm_b1_res_e.getMMweights(e, 0, isElectron, boosted, runNumber);
      } else if (!btags && boosted && !isElectron) {
        w = m_mm_b0_boosted_mu.getMMweights(e, 0, isElectron, boosted, runNumber);
      } else if (btags && boosted && !isElectron) {
        w = m_mm_b1_boosted_mu.getMMweights(e, 0, isElectron, boosted, runNumber);
      } else if (!btags && !boosted && !isElectron) {
        w = m_mm_b0_res_mu.getMMweights(e, 0, isElectron, boosted, runNumber);
      } else if (btags && !boosted && !isElectron) {
        w = m_mm_b1_res_mu.getMMweights(e, 0, isElectron, boosted, runNumber);
      }
      return w;
    }
    void initPDF(const std::string &s = "NNPDF21_lo_as_0130_100") {
      initPDFForReweighting(s.c_str(), 0);
    }
    double alphaS(double Q2) {
      return pdfAlphaS(Q2);
    }
    void setEFT(float eftLambda, float eftCvv) {
      initEFTModels(eftLambda, eftCvv);
    }
    double getEFTSMWeight(int i1_pid, int i2_pid, std::vector<int> f_pid, TLorentzVector i1, TLorentzVector i2, TLorentzVector t, TLorentzVector tbar, std::vector<TLorentzVector> f, double Q2) {
      double eftw = getEFTWeight(i1_pid, i2_pid, f_pid, i1, i2, t, tbar, f, Q2, 1.0);
      double smw  = getSMWeight(i1_pid, i2_pid, f_pid, i1, i2, t, tbar, f, Q2);
      //double smw  = getEFTWeight(i1_pid, i2_pid, f_pid, i1, i2, t, tbar, f, Q2, 0.0);
      //double s = (t+tbar).M2();
      //double tm = (i1-t).M2();
      //double u = (i1-tbar).M2();
      //std::stringstream str;
      //for (int i = 0; i < f_pid.size(); ++i) { str << f_pid[i] << " "; }
      //double test = s/(4e3*4e3*4*3.14159*alphaS(Q2));
      //if (i1_pid != 21 && f_pid.size() == 0 && std::sqrt(s) > 2000)
      //  std::cout << "EFT/SM = " << eftw << "/" << smw << "-1 = " << eftw/smw-1<<", i1, i2 = " << i1_pid<< ", " << i2_pid << ", s,t,u = " << s << ", " << tm << ", " << u << ", Q2 = " << Q2 << ", aS = " << alphaS(Q2) << ", nj = " << f_pid.size() << str.str() << ", analyt = " << test << "  -> " << (eftw/smw-1)/test << std::endl;
      return eftw/smw - 1.0;
    }
};
'''
    ROOT.gInterpreter.ProcessLine("gSystem->AddIncludePath(\" -I../../RootCoreBin/include \") ")
    ROOT.gInterpreter.ProcessLine("gSystem->AddIncludePath(\" -I../ \") ")
    ROOT.gInterpreter.ProcessLine(".include ../../RootCoreBin/include ")
    ROOT.gInterpreter.ProcessLine(".include .. ")
    ROOT.gInterpreter.ProcessLine("gROOT->LoadMacro(\"/cvmfs/atlas.cern.ch/repo/sw/ASG/AnalysisTop/2.4.12/RootCore/scripts/load_packages.C\")")
    ROOT.gInterpreter.ProcessLine("load_packages()");
    ROOT.gInterpreter.ProcessLine(bind_neutrino)
    return ROOT.WrapperExtras()

wrapperC = initBinds()

## Initialise the T2HDM class and load the precompiled modules
nameX = ""
noX   = ""
mX    = -1
tanb  = -1
stanb = "-1"
itanb = -1
typeX = -1
sba   = -999
cutsX = ""
def init2HDM(MH,MA,SBA,TANB,TYPE):
  import T2HDM
  ### Needed to modify global copies of these vars:
  global T2HDM
  global nameX
  global noX
  global mX
  global tanb
  global stanb
  global itanb
  global typeX
  global sba
  global cutsX
  ### No need for global declaration to read value of these vars
  ############
  nameX = "H" if(MH>0) else "A"
  noX   = "h" if(nameX=="A") else "h1"
  mX    = MH if(MH>0) else MA
  tanb  = TANB
  stanb = '%.2f' % tanb
  typeX = T2HDM.model.typeX ## can be also modified but should match the generation step
  sba   = T2HDM.model.sba   ## can be also modified but should match the generation step
  cutsX = T2HDM.model.cuts  ## can be also modified but should match the generation step
  T2HDM.setParameters(nameX,mX,cutsX,typeX,sba)
  itanb = T2HDM.getItanb(tanb)
  print "itanb=",itanb
  T2HDM.setModules(nameX,"All",itanb)
  print "cuts:",cutsX
  print T2HDM.modules



def loadXsec(m, fName):
    f = open(fName)
    for l in f.readlines():
        if l[0] == '#':
	    continue
        lsplit = l.split()
	if len(lsplit) == 0:
	    continue
	if len(lsplit) < 3:
	    print "I cannot understand this line:", lsplit
	    continue
	m[int(lsplit[0])] = float(lsplit[1])*float(lsplit[2])

def addFilesInChain(c, txtFileOption, n = -1):
    counter = 0
    for f in txtFileOption.split(','):
        txtf = open(f)
        for l in txtf.readlines():
	    if l[-1] == '\n':
	        l = l[0:-1]
            if counter > n and n > 0:
                break
            c.Add(l)
            counter += 1
        if counter > n:
            break

class Event:
    def __init__(self):
        pass

MV2C10_CUT = 0.6455
# systematic variations are specified according to this syntax:
# btag[A]SF_[eig][_pt[X]]__1[up/down]
# where:
# [A] = b,c,l,e for b-jet, c-jet, light-jet or extrapolation variations
# [eig] is the number of the eigen vector
# [X] is the number of the pt bin if we want to vary only jets within a pt bin
# if there is no 'btag' in the name of the variation, then the nominal SF is used
def applyBtagSF(sel, s):
    if sel.mcChannelNumber == 0:
        return 1
    ptbinInS = -1
    if 'btag' in s and 'pt' in s:
       ptbinInS = int(s.split('pt')[1][0])

    jetList = ROOT.vector('TLorentzVector')()
    jetFlavour = ROOT.vector('int')()
    jetWeight = ROOT.vector('float')()
    vetoSyst = ROOT.vector('bool')()
    for k in range(0, len(sel.tjet_pt)):
        jetList.push_back(ROOT.TLorentzVector(sel.tjet_pt[k], sel.tjet_eta[k], sel.tjet_phi[k], sel.tjet_e[k]))
        jetFlavour.push_back(int(sel.tjet_label[k]))
        jetWeight.push_back(float(sel.tjet_mv2c10[k]))
        ptbin = -1
        if ptbinInS >= 0:
            if jetList[k].Perp() < 50e3:
                ptbin = 1
            elif jetList[k].Perp() < 100e3:
                ptbin = 2
            else:
                ptbin = 3
        if ptbinInS == ptbin:
            # if this jet is in the correct pt bin, apply the syst. variation
            vetoSyst.push_back(False)
        else:
            # otherwise, take the nominal SF for this jet
            # even if we are supposed to vary it for any other criteria
            vetoSyst.push_back(True)
    syst = ""
    if 'down' in s:
       direction = 'down'
    else:
       direction = 'up'
    eig = -1
    if 'btagbSF_' in s:
        eig = int(s.split('_')[1])
        syst = "FT_EFF_Eigen_B_%d__1%s" % (eig, direction)
    elif 'btagcSF_' in s:
        eig = int(s.split('_')[1])
        syst = "FT_EFF_Eigen_C_%d__1%s" % (eig, direction)
    elif 'btaglSF_' in s:
        eig = int(s.split('_')[1])
        syst = "FT_EFF_Eigen_Light_%d__1%s" % (eig, direction)
    elif 'btageSF_0' in s:
        syst = "FT_EFF_Eigen_extrapolation__1%s" % (direction)
    elif 'btageSF_1' in s:
        syst = "FT_EFF_Eigen_extrapolation from charm__1%s" % (direction)
    return wrapperC.getBtaggingSF(jetList, jetFlavour, jetWeight, vetoSyst, syst)

# same as above, but it reads the SFs from the file
def applyBtagSFFromFile(sel, s):
    pref = 'tjet_bTagSF_70'
    varName = ''
    nomName = pref
    ptbinInS = 0
    eig = -1
    # if this is a b-tagging variation
    # then check which branch to take to apply the variation
    if 'btag' in s:
        # if we are also varying only jets in a specific pt bin, the pt bin is identified as "ptX", where X is an integer
        # get "X" as an integer from the systematic name
        if 'pt' in s:
            ptbinInS = int(s.split('pt')[1][0])
            # later, if this is != -1, we will use it to apply the variation only at the jets within this pt bin
        # get direction
        direction = 'up'
        if 'down' in s:
            direction = 'down'
        # check if it is a b,c,light or extrapolation variation and set name
        if 'btagbSF_' in s:
            varName = pref+'_eigen_B_'+direction
            # get eigenvector index
            eig = int(s.split('_')[1])
        if 'btagcSF_' in s:
            varName = pref+'_eigen_C_'+direction
            # get eigenvector index
            eig = int(s.split('_')[1])
        if 'btaglSF_' in s:
            varName = pref+'_eigen_Light_'+direction
            # get eigenvector index
            eig = int(s.split('_')[1])
        if 'btageSF_0' in s:
            varName = pref+'_syst_extrapolation_'+direction
            eig = -1
        if 'btageSF_1' in s:
            varName = pref+'_syst_extrapolation_from_charm_'+direction
            eig = -1
    else: # not a b-tagging variation, so take the nominal
        varName = pref

    btagsf = 1.0
    # loop over track jets
    for bidx in range(0, len(sel.tjet_mv2c20)):
        pb = ROOT.TLorentzVector(sel.tjet_pt[bidx], sel.tjet_eta[bidx], sel.tjet_phi[bidx], sel.tjet_e[bidx])
        # and apply object definition cuts (should already have been applied)
        if pb.Perp() > 10e3 and math.fabs(pb.Eta()) < 2.5 and sel.tjet_numConstituents[bidx] >= 2:
            # if we requested to vary only jets in a specific pt bin, check if this jet is in the pt bin
            vetoed = False
            if ptbinInS != 0:
                ptbin = 0
                if pb.Perp() < 50e3:
                    ptbin = 1
                elif pb.Perp() < 100e3:
                    ptbin = 2
                else:
                    ptbin = 3
                if ptbinInS == ptbin:
                    # if this jet is in the correct pt bin, apply the syst. variation
                    theWeight = getattr(sel, varName)
                else:
                    # otherwise, take the nominal SF for this jet
                    # even if we are supposed to vary it for any other criteria
                    theWeight = getattr(sel, nomName)
                    vetoed = True
            else: # we have not requested the pt binned variation
                # so we always vary all jets, regardless of the pt bin
                # if no variation was requested, varName will be set to the nominal above
                theWeight = getattr(sel, varName)
            if eig > -1 and not vetoed: # branches in b,c,light variation have a second index related to the eigenvector number
                btagsf *= theWeight[bidx][eig]
            else: # extrapolation branches or the nominal, have no second index
                btagsf *= theWeight[bidx]
    return btagsf

# TODO
def readEvent(mt):
    return mt

listEWK = [410000, 301528, 301529, 301530, 301531, 301532]
def applyEWK(sel, s):
    if sel.mcChannelNumber == 0:
        return 1
    top = ROOT.TLorentzVector()
    top.SetPtEtaPhiM(sel.MC_t_pt, sel.MC_t_eta, sel.MC_t_phi, sel.MC_t_m)
    topbar = ROOT.TLorentzVector()
    topbar.SetPtEtaPhiM(sel.MC_tbar_pt, sel.MC_tbar_eta, sel.MC_tbar_phi, sel.MC_tbar_m)
    if s == 'ttEWK__1up':
        w = wrapperC.getEWK(top, topbar, sel.initial_type, 1)
    elif s == 'ttEWK__1down':
        w = wrapperC.getEWK(top, topbar, sel.initial_type, 2)
    else:
        w = wrapperC.getEWK(top, topbar, sel.initial_type, 0)
    return w


'''
SM topo's         |  H/A matched topo's  | H/A partially matched topo's
------------------+----------------------+-------------------------------
P0_gg_ttx         |  P0_gg_ttx           |  
P0_uux_ttx        |  P0_uux_ttx          |  P0_bbx_ttx
P1_uux_ttxg       |  P1_uux_ttxg         |  P1_bbx_ttxg
P1_gux_ttxux      |  P1_gux_ttxux        |  P1_gbx_ttxbx 
P1_gu_ttxu        |  P1_gu_ttxu          |  P1_gb_ttxb
P1_gg_ttxg        |  P1_gg_ttxg          |  
P2_gu_ttxgu       |  P2_gu_ttxgu         |  P2_gb_ttxgb
P2_uxux_ttxuxux   |  P2_uxux_ttxuxux     |  P2_bxbx_ttxbxbx
P2_uxcx_ttxuxcx   |  P2_uxcx_ttxuxcx     |  P2_uxbx_ttxuxbx
P2_uux_ttxuux     |  P2_uux_ttxuux       |  P2_bbx_ttxbbx
P2_uux_ttxgg      |  P2_uux_ttxgg        |  P2_bbx_ttxgg
P2_uux_ttxccx     |  P2_uux_ttxccx       |  P2_uux_ttxbbx  P2_bbx_ttxuux(?)
P2_uu_ttxuu       |  P2_uu_ttxuu         |  P2_bb_ttxbb
P2_ucx_ttxucx     |  P2_ucx_ttxucx       |  P2_ubx_ttxubx   
P2_uc_ttxuc       |  P2_uc_ttxuc         |  P2_ub_ttxub
P2_gux_ttxgux     |  P2_gux_ttxgux       |  P2_gbx_ttxgbx
P2_gg_ttxgg       |  P2_gg_ttxgg         |     
P2_gg_ttxuux      |  P2_gg_ttxuux        |  P2_gg_ttxbbx

H/A unmatched topo's:
P2_uxb_ttxuxb (?)
'''

def getTopology(prod,id1,id2=999):
    isProd = True if(prod!="") else False
    topology = ""
    GLU=21
    TOP=6
    CHM=4
    idsum = id1+id2
    absidsum = abs(id1)+abs(id2)
    absidsub = abs(abs(id1)-abs(id2))
    if(abs(id2)>GLU):
        if(id1==GLU):                            topology = "g"
        if(abs(id1)<TOP and id1>0):              topology = "u"
        if(abs(id1)<TOP and id1<0):              topology = "ux"
    elif(idsum==2*GLU):                          topology = "gg"
    elif(idsum==0 and id1<TOP):                  topology = "uux"
    elif(idsum==0 and id1>=CHM and prod=="uux"): topology = "ccx"
    elif(idsum>GLU and idsum<2*GLU):             topology = "gu"
    elif(idsum<GLU and idsum>=GLU-TOP):          topology = "gux"
    elif(absidsum<2*TOP):                    
        if(id1>0           and absidsub==0):     topology = "uu"
        if(id1<0           and absidsub==0):     topology = "uxux"
        if(id1>0 and id2>0 and absidsub>0):      topology = "uc"
        if(id1*id2<0       and absidsub>0):      topology = "ucx"
        if(id1<0 and id2<0 and absidsub>0):      topology = "uxcx"
    else:
        print "Cannot understand the topology - quitting"
        quit()
    return topology

def combineTopology(ids):
    topology = ""
    n = len(ids)
    if(n==4): topology += "P0_"
    if(n==5): topology += "P1_"
    if(n==6): topology += "P2_"
    prod = getTopology("",ids[0],ids[1])
    topology += prod
    topology += "_ttx" ## should be always there...
    decay = ""
    if(n==5): decay = getTopology(prod,ids[4])
    if(n==6): decay = getTopology(prod,ids[4],ids[5])
    topology += decay
    return topology

def stripTopology(topology):
    topology = topology.replace("P0_","")
    topology = topology.replace("P1_","")
    topology = topology.replace("P2_","")
    topology = topology.replace("_","")
    return topology

def correct2HDMTopology(sel,topology):
    n = len(sel.MC_id_me)
    topoparts = topology.split("_")
    if(len(topoparts)<3):
        print "Error: topoparts=",topoparts
        print "Too few topo parts - quitting"
        quit()
    BOT=5
    if(n==4 and abs(sel.MC_id_me[0])!=BOT and abs(sel.MC_id_me[1])!=BOT): return topology
    if(n==5 and abs(sel.MC_id_me[0])!=BOT and abs(sel.MC_id_me[1])!=BOT and abs(sel.MC_id_me[4])!=BOT): return topology
    if(n==6 and abs(sel.MC_id_me[0])!=BOT and abs(sel.MC_id_me[1])!=BOT and abs(sel.MC_id_me[4])!=BOT and abs(sel.MC_id_me[5])!=BOT): return topology
    prod  = topoparts[1]
    decay = topoparts[2].replace("ttx","")
    if(abs(sel.MC_id_me[0])==BOT or abs(sel.MC_id_me[1])==BOT):
        if(prod=="uux"):                                               prod = "bbx"
        elif(prod=="gu" or prod=="gux" or prod=="uu" or prod=="uxux"): prod = prod.replace("u","b")
        elif(prod=="uc" or prod=="ucx" or prod=="uxcx"):               prod = prod.replace("c","b")
        else:
            print "Cannot understand the production with b's - quitting:",prod
            print "Topology is:",topology
            quit()
    if(n==5 and abs(sel.MC_id_me[4])==BOT):
        if(decay=="u" or decay=="ux"): decay = decay.replace("u","b")
        else:
            print "Cannot understand the decay with b's (n=5)- quitting:",decay
            quit()
    if(n==6 and (abs(sel.MC_id_me[4])==BOT or abs(sel.MC_id_me[5])==BOT)):
        if(decay=="uux"):                                                    decay = "bbx"
        elif(decay=="gu" or decay=="gux" or decay=="uu"   or decay=="uxux"): decay = decay.replace("u","b")
        elif(decay=="uc" or decay=="ucx" or decay=="uxcx" or decay=="ccx"):  decay = decay.replace("c","b")
        else:
            print "Cannot understand the decay with b's (n=6) - quitting:",decay
            quit()
    topology = topoparts[0]+"_"+prod+"_ttx"+decay
    for i in xrange(3,len(topoparts)): topology += "_"+topoparts[i]
    return topology

def getMomenta(sel,topology):
    n = len(sel.MC_id_me)
    topoparts = topology.split("_")
    if(len(topoparts)<3):
        print "Error: topoparts=",topoparts
        print "Too few topo parts - quitting"
        quit()
    prod  = topoparts[1]
    decay = topoparts[2].replace("ttx","")
    GLU=21
    TOP=6
    CHM=4
    me = {}
    ### Fill production partons' momenta ###
    if(prod=="gg" or prod=="uu" or prod=="uxux"):
        me.update({"p1":{"id":sel.MC_id_me[0],"idx":0}})
        me.update({"p2":{"id":sel.MC_id_me[1],"idx":1}})
    elif(prod=="uux"):
        q  = 0 if(sel.MC_id_me[0]>0) else 1
        qx = 1 if(sel.MC_id_me[0]>0) else 0
        me.update({"p1":{"id":sel.MC_id_me[q], "idx":q}})
        me.update({"p2":{"id":sel.MC_id_me[qx],"idx":qx}})
    elif(prod=="gu" or prod=="gux"):
        g = 0 if(sel.MC_id_me[0]==GLU) else 1
        q = 1 if(sel.MC_id_me[0]==GLU) else 0
        me.update({"p1":{"id":sel.MC_id_me[g],"idx":g}})
        me.update({"p2":{"id":sel.MC_id_me[q],"idx":q}})
    elif(prod=="uc" or prod=="ucx" or prod=="uxcx"):
        qLight = 0 if(abs(sel.MC_id_me[0])<abs(sel.MC_id_me[1])) else 1
        qHeavy = 1 if(abs(sel.MC_id_me[0])<abs(sel.MC_id_me[1])) else 0
        me.update({"p1":{"id":sel.MC_id_me[qLight],"idx":qLight}})
        me.update({"p2":{"id":sel.MC_id_me[qHeavy],"idx":qHeavy}})
    else:
        print "Cannot understand the production - quitting"
        quit()
    ### Fill top quarks' momenta ###
    t  = 2 if(sel.MC_id_me[2]>0) else 3
    tx = 3 if(sel.MC_id_me[2]>0) else 2
    me.update({"t1":{"id":sel.MC_id_me[t], "idx":t}})
    me.update({"t2":{"id":sel.MC_id_me[tx],"idx":tx}}) 
    ### Fill extra partons' momenta ###
    if(n==5):
        if(decay=="g" or decay=="u" or decay=="ux"):
            me.update({"j1":{"id":sel.MC_id_me[4],"idx":4}})
    if(n==6):
        if(decay=="gg" or decay=="uu" or decay=="uxux"):
            me.update({"j1":{"id":sel.MC_id_me[4],"idx":4}})
            me.update({"j2":{"id":sel.MC_id_me[5],"idx":5}})
        elif(decay=="uux" or decay=="ccx"):
            q  = 4 if(sel.MC_id_me[4]>0) else 5
            qx = 5 if(sel.MC_id_me[4]>0) else 4
            me.update({"j1":{"id":sel.MC_id_me[q], "idx":q}})
            me.update({"j2":{"id":sel.MC_id_me[qx],"idx":qx}})
        elif(decay=="gu" or decay=="gux"):
            g = 4 if(sel.MC_id_me[4]==GLU) else 5
            q = 5 if(sel.MC_id_me[4]==GLU) else 4
            me.update({"j1":{"id":sel.MC_id_me[g],"idx":g}})
            me.update({"j2":{"id":sel.MC_id_me[q],"idx":q}})
        elif(decay=="uc" or decay=="ucx" or decay=="uxcx"):
            qLight = 4 if(abs(sel.MC_id_me[4])<abs(sel.MC_id_me[5])) else 5
            qHeavy = 5 if(abs(sel.MC_id_me[4])<abs(sel.MC_id_me[5])) else 4
            me.update({"j1":{"id":sel.MC_id_me[qLight],"idx":qLight}})
            me.update({"j2":{"id":sel.MC_id_me[qHeavy],"idx":qHeavy}})
        else:
            print "Cannot understand the decay - quitting"
            quit()
    ### Fill the momenta array
    p = [[ sel.MC_e_me[me["p1"]["idx"]], sel.MC_px_me[me["p1"]["idx"]], sel.MC_py_me[me["p1"]["idx"]], sel.MC_pz_me[me["p1"]["idx"]] ],
         [ sel.MC_e_me[me["p2"]["idx"]], sel.MC_px_me[me["p2"]["idx"]], sel.MC_py_me[me["p2"]["idx"]], sel.MC_pz_me[me["p2"]["idx"]] ],
         [ sel.MC_e_me[me["t1"]["idx"]], sel.MC_px_me[me["t1"]["idx"]], sel.MC_py_me[me["t1"]["idx"]], sel.MC_pz_me[me["t1"]["idx"]] ],
         [ sel.MC_e_me[me["t2"]["idx"]], sel.MC_px_me[me["t2"]["idx"]], sel.MC_py_me[me["t2"]["idx"]], sel.MC_pz_me[me["t2"]["idx"]] ]]
    if(n>4): p.append([ sel.MC_e_me[me["j1"]["idx"]], sel.MC_px_me[me["j1"]["idx"]], sel.MC_py_me[me["j1"]["idx"]], sel.MC_pz_me[me["j1"]["idx"]] ])
    if(n>5): p.append([ sel.MC_e_me[me["j2"]["idx"]], sel.MC_px_me[me["j2"]["idx"]], sel.MC_py_me[me["j2"]["idx"]], sel.MC_pz_me[me["j2"]["idx"]] ])
    return p

def getEFTSMWeight(sel, s = ""):
    #double getEFTSMWeight(int i1_pid, int i2_pid, std::vector<int> f_pid, TLorentzVector i1, TLorentzVector i2, TLorentzVector t, TLorentzVector tbar, std::vector<TLorentzVector> f, double Q2) {
    if sel.mcChannelNumber == 0:
        return 1
    #if not sel.mcChannelNumber in [407200, 407201, 407202, 407203, 407204]:
    #    return 1
    i1 = ROOT.TLorentzVector()
    i2 = ROOT.TLorentzVector()
    top = ROOT.TLorentzVector()
    topbar = ROOT.TLorentzVector()
    i1.SetPxPyPzE(sel.MC_px_me[0], sel.MC_py_me[0], sel.MC_pz_me[0], sel.MC_e_me[0])
    i1_pid = sel.MC_id_me[0]
    i2.SetPxPyPzE(sel.MC_px_me[1], sel.MC_py_me[1], sel.MC_pz_me[1], sel.MC_e_me[1])
    i2_pid = sel.MC_id_me[1]
    if sel.MC_id_me[2] > 0:
        top.SetPxPyPzE(sel.MC_px_me[2], sel.MC_py_me[2], sel.MC_pz_me[2], sel.MC_e_me[2])
        topbar.SetPxPyPzE(sel.MC_px_me[3], sel.MC_py_me[3], sel.MC_pz_me[3], sel.MC_e_me[3])
    else:
        top.SetPxPyPzE(sel.MC_px_me[3], sel.MC_py_me[3], sel.MC_pz_me[3], sel.MC_e_me[3])
        topbar.SetPxPyPzE(sel.MC_px_me[2], sel.MC_py_me[2], sel.MC_pz_me[2], sel.MC_e_me[2])
    f = ROOT.vector('TLorentzVector')()
    f_pid = ROOT.vector('int')()
    for idx in range(4,len(sel.MC_px_me)):
        f_pid.push_back(sel.MC_id_me[idx])
        f.push_back(ROOT.TLorentzVector(sel.MC_px_me[idx], sel.MC_py_me[idx], sel.MC_pz_me[idx], sel.MC_e_me[idx]))
    Q2 = (1e-3*sel.MC_Q_me)**2
    if 'eftScale' in s and 'up' in s:
    	Q2 *= 4
    elif 'eftScale' in s and 'down' in s:
    	Q2 *= 0.25
    w = wrapperC.getEFTSMWeight(i1_pid, i2_pid, f_pid, i1, i2, top, topbar, f, Q2)
    return w


def print2HDM(sel,topology,alphaS,ids,pCpp,me2XX,me2SM):
    pt1 = ROOT.TLorentzVector()
    pt2 = ROOT.TLorentzVector()
    pt1.SetPxPyPzE(sel.MC_px_me[2],sel.MC_py_me[2],sel.MC_pz_me[2],sel.MC_e_me[2])
    pt2.SetPxPyPzE(sel.MC_px_me[3],sel.MC_py_me[3],sel.MC_pz_me[3],sel.MC_e_me[3])
    if(abs((pt1+pt2).M()-mX)<100):
       weightX = me2XX/me2SM
       print "# 2HDM setup: mX=%g, sba=%g, tanb=%g, type=%g --> topo=%s, mtt=%g, me2XX/me2SM=%g/%g=%g" % (mX,sba,tanb,typeX,topology,(pt1+pt2).M(),me2XX,me2SM,weightX)
       print "# ME ids are: ",ids
       print "# Q =",sel.MC_Q_me/1000.
       print "alphaS =",alphaS
       print "p=["
       for x in xrange(len(ids)):
          print "  [%g,%g,%g,%g]," % (sel.MC_e_me[x],sel.MC_px_me[x],sel.MC_py_me[x],sel.MC_pz_me[x])
       print "]"
       print "used p =", pCpp



def get2HDMWeight(sel):
    if sel.mcChannelNumber == 0:
        return 1
    if not sel.mcChannelNumber in [407200, 407201, 407202, 407203, 407204]:
        return 1
    topologySM_name = combineTopology(sel.MC_id_me)
    topologyXX_name = topologySM_name+"_no_"+noX
    topologyXX_name = correct2HDMTopology(sel,topologyXX_name)
    topologySM_lib = 'matrix2SM'+stripTopology(topologySM_name)+'py'
    topologyXX_lib = 'matrix2'+nameX+str(itanb)+stripTopology(topologyXX_name)+'py'
    ids = []
    for i in xrange(sel.MC_id_me.size()): ids.append(sel.MC_id_me[i])
    if(topologySM_lib not in T2HDM.modules or topologyXX_lib not in T2HDM.modules):
        print "topology %s or %s (or both) is not in T2HDM.modules" % (topologySM_lib,topologyXX_lib) 
        print "ME id's are: ",ids
        return 1

    pCpp = getMomenta(sel,topologySM_name)
    pPython = T2HDM.invert_momenta(pCpp)
    Q2 = (sel.MC_Q_me/1000.)*(sel.MC_Q_me/1000.)
    alphaS = wrapperC.alphaS(Q2)
    nhel = 0 # sum over all helicities
    me2SM = T2HDM.modules[topologySM_lib].get_me(pPython,alphaS,nhel) ### calculate the SM ME^2
    me2XX = T2HDM.modules[topologyXX_lib].get_me(pPython,alphaS,nhel) ### calculate the X ME^2
    weightX = me2XX/me2SM                                             ### calculate the weight
    #print2HDM(sel,topologySM_name,alphaS,ids,pCpp,me2XX,me2SM)        ### print basic info
    
    return weightX


listWjets22 = []
listWjets22.extend(range(363330, 363354+1))
listWjets22.extend(range(363436, 363483+1))
#def applyWjets22Weight(sel):
#    if not sel.mcChannelNumber in listWjets22:
#        return 1.0
#    njet = 0
#    for i in range(0, len(sel.akt4truthjet_pt)):
#        if sel.akt4truthjet_pt[i] > 20e3 and fabs(sel.akt4truthjet_eta[i]) < 4.5:
#	    njet += 1
#    return wrapperC.getWjets22Weight(njet)

# list of systematics related to changes in the weights only
weightSF = {'' : ['pileup', 'leptonSF', 'jvt'],
            'eTrigSF__1up': ['pileup', 'leptonSF_EL_SF_Trigger_UP', 'jvt'],
            'eTrigSF__1down': ['pileup', 'leptonSF_EL_SF_Trigger_DOWN', 'jvt'],
            'eIDSF__1up': ['pileup', 'leptonSF_EL_SF_ID_UP', 'jvt'],
            'eIDSF__1down': ['pileup', 'leptonSF_EL_SF_ID_DOWN', 'jvt'],
            'eRecoSF__1up': ['pileup', 'leptonSF_EL_SF_Reco_UP', 'jvt'],
            'eRecoSF__1down': ['pileup', 'leptonSF_EL_SF_Reco_DOWN', 'jvt'],
            'eIsolSF__1up': ['pileup', 'leptonSF_EL_SF_Isol_UP', 'jvt'],
            'eIsolSF__1down': ['pileup', 'leptonSF_EL_SF_Isol_DOWN', 'jvt'],
            'muTrigStatSF__1up': ['pileup', 'leptonSF_MU_SF_Trigger_STAT_UP', 'jvt'],
            'muTrigStatSF__1down': ['pileup', 'leptonSF_MU_SF_Trigger_STAT_DOWN', 'jvt'],
            'muTrigSystSF__1up': ['pileup', 'leptonSF_MU_SF_Trigger_SYST_UP', 'jvt'],
            'muTrigSystSF__1down': ['pileup', 'leptonSF_MU_SF_Trigger_SYST_DOWN', 'jvt'],
            'muIDStatSF__1up': ['pileup', 'leptonSF_MU_SF_ID_STAT_UP', 'jvt'],
            'muIDStatSF__1down': ['pileup', 'leptonSF_MU_SF_ID_STAT_DOWN', 'jvt'],
            'muIDSystSF__1up': ['pileup', 'leptonSF_MU_SF_ID_SYST_UP', 'jvt'],
            'muIDSystSF__1down': ['pileup', 'leptonSF_MU_SF_ID_SYST_DOWN', 'jvt'],
            'muIsolStatSF__1up': ['pileup', 'leptonSF_MU_SF_Isol_STAT_UP', 'jvt'],
            'muIsolStatSF__1down': ['pileup', 'leptonSF_MU_SF_Isol_STAT_DOWN', 'jvt'],
            'muIsolSystSF__1up': ['pileup', 'leptonSF_MU_SF_Isol_SYST_UP', 'jvt'],
            'muIsolSystSF__1down': ['pileup', 'leptonSF_MU_SF_Isol_SYST_DOWN', 'jvt'],
            'pileupSF__1up': ['pileup_UP', 'leptonSF', 'jvt'],
            'pileupSF__1down': ['pileup_DOWN', 'leptonSF', 'jvt'],
            'jvtSF__1up': ['pileup', 'leptonSF', 'jvt_UP'],
            'jvtSF__1down': ['pileup', 'leptonSF', 'jvt_DOWN'],
	   }

weightChangeSystematics = weightSF.keys()

