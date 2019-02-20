import helpers
import ROOT
import math
from array import array


class Analysis:
	ch = ''
	fi = None
	histSuffixes = [] # systematic copies of histograms
	h = {} # map of histogram names to map of systematics to histograms
	noMttSlices = False
	applyMET = 0
	eftLambda = -1
	eftCvv = 0
        w2HDM = 1
	def __init__(self, channel, suf, outputFile):
		self.fi = ROOT.TFile(outputFile, "recreate")
		self.ch = channel
		self.histSuffixes = suf
		self.noMttSlices = False
		self.applyMET = 0
		self.eftLambda = -1
		self.eftCvv = -1
		self.w2HDM = 1
		self.h = {}
		self.keep = '' # can be 'bb', 'cc', 'c' or 'l' and only applies to W+jets

	def add(self, hName, nBins, xLow, xHigh):
		self.h[hName] = {}
		#self.fi.cd()
		for s in self.histSuffixes:
			#print "adding histogram with name ", hName+self.ch+s
			self.h[hName][s] = ROOT.TH1D(hName+self.ch+s, "", nBins, xLow, xHigh)
			self.h[hName][s].Sumw2()
			self.h[hName][s].SetDirectory(0)

	def addVar(self, hName, nBinsList):
		ar = array("d", nBinsList)
		#self.fi.cd()
		self.h[hName] = {}
		for s in self.histSuffixes:
			self.h[hName][s] = ROOT.TH1D(hName+self.ch+s, "", len(nBinsList) - 1, ar)
			self.h[hName][s].Sumw2()
			self.h[hName][s].SetDirectory(0)

	def add2D(self, hName, nBins, xLow, xHigh, nBinsY, yLow, yHigh):
		self.h[hName] = {}
		#self.fi.cd()
		for s in self.histSuffixes:
			#print "adding histogram with name ", hName+self.ch+s
			self.h[hName][s] = ROOT.TH2D(hName+self.ch+s, "", nBins, xLow, xHigh, nBinsY, yLow, yHigh)
			self.h[hName][s].Sumw2()
			self.h[hName][s].SetDirectory(0)

	def write(self):
		self.fi.cd()
		for hName in self.h:
			for s in self.histSuffixes:
				#print "writing histogram with name ", hName+s, " in file ",self.fi.GetName()
				self.h[hName][s].Write(hName+s)
		self.fi.Close()

	def getWeight(self, sel, s):
		# this applies all the weights that come out of the box
		if sel.mcChannelNumber == 0:
			return 1.0
		k = s
		if not s in helpers.weightSF:
			k = ''
		listWeights = helpers.weightSF[k]
		weight = 1.0
		for item in listWeights:
			wItem = getattr(sel, 'weight_'+item)
			weight *= wItem

		# this applies the EWK weight
		channel = sel.mcChannelNumber
		if channel in helpers.listEWK and not self.noMttSlices:
			weight *= helpers.applyEWK(sel, s)

		# this applies the W+jets Sherpa 2.2 nJets reweighting correction
		#weight *= helpers.applyWjets22Weight(sel)
		weight *= sel.weight_Sherpa22_corr

		return weight


class AnaTtresSL(Analysis):
	def __init__(self, channel, suf, outputFile):
		Analysis.__init__(self, channel, suf, outputFile)
		self.applyQCD = False
		self.noMttSlices = False
		self.applyMET = 0
		self.eftLambda = -1
		self.eftCvv = 0
		self.w2HDM = 1
		# make histograms
		self.add("yields", 1, 0.5, 1.5)
		self.add("yieldsPos", 1, 0.5, 1.5)
		self.add("yieldsNeg", 1, 0.5, 1.5)
		self.add("runNumber", 24647, 276261.5, 300908.5)
		self.add("lepPt", 100, 25, 525)
		self.add("lepEta", 20, -2.5, 2.5)
		self.add("lepPhi", 32, -3.2, 3.2)
		self.add("nJets", 10, -0.5, 9.5)
		self.add("nTrkBtagJets", 10, -0.5, 9.5)
		self.addVar("MET", [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 340, 380, 450, 500])
		self.add("MET_phi", 32, -3.2, 3.2)
		self.add("mwt", 20, 0, 200)
		self.add("closestJetDr", 20, 0, 2.0)
		self.add("closestJetPt", 20, 0, 200)
		self.add("mu", 100, 0, 100)
		self.add("vtxz", 40, -400, 400)
		self.add("npv", 50, 0, 50)
		self.addVar("closeJetPt", [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 340, 380, 450, 500])
		self.addVar("largeJetPt", [300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 540, 580, 620, 660, 700, 800, 1e3, 1.2e3, 1.5e3])
		self.add("largeJetM", 30, 0, 300)
		self.add("largeJetPtMtt", 50, 0, 1)
		self.add("largeJetEta", 20, -2., 2.)
		self.add("largeJetPhi", 32, -3.2, 3.2)
		self.addVar("mtlep_boo", [80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500])
		self.addVar("mtlep_res", [80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500])
		self.addVar("mthad_res", [80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500])
		self.add("mwhad_res", 40, 0, 400)
		self.add("chi2", 50, -3, 7)
                self.add("trueMttLow1", 20,300,1000)
                self.add("trueMttLow1r", 20,300,1000)
                self.add("trueMttLow2", 25,300,800)
                self.add("trueMttLow2r", 25,300,800)
                self.add("mttLow1", 20,300,1000)
                self.add("mttLow1r", 20,300,1000)
                self.add("mttLow2", 25,300,800)
                self.add("mttLow2r", 25,300,800)
		self.addVar("mtt", [0, 80, 160, 240, 320, 400, 480, 560,640,720,800,920,1040,1160,1280,1400,1550,1700,2000,2300,2600,2900,3200,3600,4100,4600,5100,6000])
		self.addVar("mttr", [0, 80, 160, 240, 320, 400, 480, 560,640,720,800,920,1040,1160,1280,1400,1550,1700,2000,2300,2600,2900,3200,3600,4100,4600,5100,6000])
		self.addVar("mttPos", [0, 80, 160, 240, 320, 400, 480, 560,640,720,800,920,1040,1160,1280,1400,1550,1700,2000,2300,2600,2900,3200,3600,4100,4600,5100,6000])
		self.addVar("mttNeg", [0, 80, 160, 240, 320, 400, 480, 560,640,720,800,920,1040,1160,1280,1400,1550,1700,2000,2300,2600,2900,3200,3600,4100,4600,5100,6000])
		self.addVar("trueMtt", [0, 80, 160, 240, 320, 400, 480, 560,640,720,800,920,1040,1160,1280,1400,1550,1700,2000,2300,2600,2900,3200,3600,4100,4600,5100,6000])
		self.addVar("trueMttr", [0, 80, 160, 240, 320, 400, 480, 560,640,720,800,920,1040,1160,1280,1400,1550,1700,2000,2300,2600,2900,3200,3600,4100,4600,5100,6000])
		self.add("btagSF", 50, 0, 2)

		self.add("largeJet_tau32_wta", 20, 0, 1)
		self.add("largeJet_tau21_wta", 20, 0, 1)

	# only apply the reco weights
	def getWeight(self, sel, s):
		if sel.mcChannelNumber == 0:
  			return 1.0
		weight = Analysis.getWeight(self, sel, s)

		# just add the btagging SFs on top of those, as this Analysis implementation applies b-tagging
		#### warning: disabled for now in mc15c
		btagsf = helpers.applyBtagSFFromFile(sel, s)
		weight *= btagsf

		# for EFT
		if self.eftLambda > 0:
			weight *= helpers.getEFTSMWeight(sel, s)

                # for 2HDM
                self.w2HDM = 1
                if(helpers.nameX != ""):
                        self.w2HDM = helpers.get2HDMWeight(sel)
                        weight *= self.w2HDM

		# W+jets C/A and HF syst. variations
		# assuming b-tagging
		hfweight = 1.0
		f_ca = 1.0

		nj = 4
		if len(sel.jet_pt) < 4:
			nj = 2
		year = "2015"
		if sel.runNumber > 297730:
			year = "2016"
		frac = {}
		frac["2015"] = {}
		frac["2016"] = {}
		frac["2015"][2] = {}
		frac["2016"][2] = {}
		frac["2015"][4] = {}
		frac["2016"][4] = {}

		frac["2015"][2]['e'] = {'bb': 0.091, 'cc': 0.070, 'c': 0.226, 'l': 0.613}
		frac["2015"][2]['mu'] = {'bb': 0.088, 'cc': 0.065, 'c': 0.216, 'l': 0.631}

		frac["2016"][2]['e'] = {'bb': 0.090, 'cc': 0.067, 'c': 0.232, 'l': 0.611}
		frac["2016"][2]['mu'] = {'bb': 0.089, 'cc': 0.065, 'c': 0.221, 'l': 0.625}

		frac["2015"][4]['e'] = {'bb': 0.140, 'cc': 0.157, 'c': 0.220, 'l': 0.483}
		frac["2015"][4]['mu'] = {'bb': 0.153, 'cc': 0.155, 'c': 0.210, 'l': 0.482}

		frac["2016"][4]['e'] = {'bb': 0.148, 'cc': 0.155, 'c': 0.211, 'l': 0.487}
		frac["2016"][4]['mu'] = {'bb': 0.146, 'cc': 0.149, 'c': 0.213, 'l': 0.492}

		f_ca_map = {}
		f_ca_map["2015"] = {}
		f_ca_map["2016"] = {}

		f_ca_map["2015"][2] = {'e': 1.01, 'mu': 1.02}
		f_ca_map["2016"][2] = {'e': 0.91, 'mu': 1.11}

		f_ca_map["2015"][4] = {'e': 0.76, 'mu': 0.94}
		f_ca_map["2016"][4] = {'e': 0.90, 'mu': 0.82}

		flav_map = {}
		flav_map["2015"] = {}
		flav_map["2016"] = {}
		flav_map["2015"][2] = {}
		flav_map["2016"][2] = {}
		flav_map["2015"][4] = {}
		flav_map["2016"][4] = {}

		flav_map["2015"][2]["e"] = {'bb': 1.57, 'cc': 1.57, 'c': 1.00, 'l': 0.85}
		flav_map["2016"][2]["e"] = {'bb': 2.06, 'cc': 2.06, 'c': 1.00, 'l': 0.73}

		flav_map["2015"][2]["mu"] = {'bb': 1.49, 'cc': 1.49, 'c': 1.00, 'l': 0.88}
		flav_map["2016"][2]["mu"] = {'bb': 1.70, 'cc': 1.70, 'c': 1.00, 'l': 0.83}

		flav_map["2015"][4]["e"] = {'bb': 1.43, 'cc': 1.43, 'c': 0.91, 'l': 0.78}
		flav_map["2016"][4]["e"] = {'bb': 1.73, 'cc': 1.73, 'c': 0.84, 'l': 0.61}

		flav_map["2015"][4]["mu"] = {'bb': 1.37, 'cc': 1.37, 'c': 0.91, 'l': 0.80}
		flav_map["2016"][4]["mu"] = {'bb': 1.52, 'cc': 1.52, 'c': 0.89, 'l': 0.74}
		for c in frac[year][nj]:
			for f in frac[year][nj][c]:
				frac[year][nj][c][f] *= flav_map[year][nj][c][f]
		flav_map_unc = {}
		flav_map_unc["2015"] = {}
		flav_map_unc["2016"] = {}
		flav_map_unc["2015"]["e"] = {}
		flav_map_unc["2016"]["e"] = {}
		flav_map_unc["2015"]["mu"] = {}
		flav_map_unc["2016"]["mu"] = {}
		flav_map_unc["2015"]['e'] = {'bb': 0.253, 'cc': 0.253, 'c': 0.50, 'l': 0.123}
		flav_map_unc["2015"]['mu'] = {'bb': 0.100, 'cc': 0.100, 'c': 0.50, 'l': 0.041}
		flav_map_unc["2016"]['e'] = {'bb': 0.043, 'cc': 0.043, 'c': 0.50, 'l': 0.031}
		flav_map_unc["2016"]['mu'] = {'bb': 0.046, 'cc': 0.046, 'c': 0.50, 'l': 0.023}
		#flav_map_unc["2015"]['e'] = {'bb': 0.500, 'cc': 0.500, 'c': 0.50, 'l': 0.500}
		#flav_map_unc["2015"]['mu'] = {'bb': 0.500, 'cc': 0.500, 'c': 0.50, 'l': 0.500}
		#flav_map_unc["2016"]['e'] = {'bb': 0.500, 'cc': 0.500, 'c': 0.50, 'l': 0.500}
		#flav_map_unc["2016"]['mu'] = {'bb': 0.500, 'cc': 0.500, 'c': 0.50, 'l': 0.500}

		chan = ''
		if len(sel.el_pt) == 1:
			chan = 'e'
		elif len(sel.mu_pt) == 1:
			chan = 'mu'
		flav = ''
		if sel.mcChannelNumber in helpers.listWjets22:
			flag = sel.Wfilter_Sherpa_nT
			if flag == 3 or flag == 4:
				flav = 'bb'
			elif flag == 1:
				flav = 'cc'
			elif flag == 2:
				flav = 'c'
			elif flag == 5:
				flav = 'l'
			f_ca = f_ca_map[year][nj][chan]
			norm = 1.0
			hfweight = flav_map[year][nj][chan][flav]
			flavunc = ''
			if 'wbb__' in s:
				flavunc = 'bb'
			elif 'wcc__' in s:
				flavunc = 'cc'
			elif 'wc__' in s:
				flavunc = 'c'
			elif 'wl__' in s:
				flavunc = 'l'
			updown = 0
			if '1up' in s:
				updown = 1
			elif '1down' in s:
				updown = -1
			norm = 0
			for f in ['bb', 'cc', 'c', 'l']:
				if flavunc == f:
					norm += (1+updown*flav_map_unc[year][chan][f])*frac[year][nj][chan][f]
					hfweight *= (1+updown*flav_map_unc[year][chan][f])
				else:
					norm += frac[year][nj][chan][f]
			if s == "wnorm__1up":
				f_ca *= 1.10
			elif s == "wnorm__1down":
				f_ca *= 0.90
			hfweight /= norm
			weight *= f_ca*hfweight
		return weight

	def qcdWeight(self, sel):
		if sel.mcChannelNumber != 0:
			return 1

		nBtag = 0
		for bidx in range(0, len(sel.tjet_mv2c10)):
			pb = ROOT.TLorentzVector(sel.tjet_pt[bidx], sel.tjet_eta[bidx], sel.tjet_phi[bidx], sel.tjet_e[bidx])
			if pb.Perp() > 10e3 and math.fabs(pb.Eta()) < 2.5 and sel.tjet_numConstituents[bidx] >= 2:
				if sel.tjet_mv2c10[bidx] > helpers.MV2C10_CUT:
					nBtag += 1
		isBoosted = 0
		#if 'be' in self.ch or 'bmu' in self.ch:
		#	isBoosted = 1

		met = ROOT.TLorentzVector(sel.met_met, 0, sel.met_phi, sel.met_met)
		l = ROOT.TLorentzVector()
		lisTight = 0
		lsd0 = 0
		isElectron = 0
		muonTrigger = 0
		topoetcone20 = 0
		runNumber = sel.runNumber
		if len(sel.el_pt) == 1:
			l.SetPtEtaPhiE(sel.el_pt[0], sel.el_eta[0], sel.el_phi[0], sel.el_e[0])
			lisTight = sel.el_isTight[0]
			lisTight = bool(int(lisTight.encode('hex'), 16))
			lsd0 = sel.el_d0sig[0]
			topoetcone20 = sel.el_topoetcone20[0]
			isElectron = 1
		elif len(sel.mu_pt) == 1:
			l.SetPtEtaPhiE(sel.mu_pt[0], sel.mu_eta[0], sel.mu_phi[0], sel.mu_e[0])
			lisTight = sel.mu_isTight[0]
			lisTight = bool(int(lisTight.encode('hex'), 16))
			lsd0 = sel.mu_d0sig[0]
			muonTrigger = bool(int(sel.mu_trigMatch_HLT_mu20_iloose_L1MU15[0].encode('hex'), 16)) or bool(int(sel.mu_trigMatch_HLT_mu50[0].encode('hex'), 16))
			topoetcone20 = sel.mu_topoetcone20[0]
		jets = ROOT.vector('TLorentzVector')()
		for k in range(0, len(sel.jet_pt)):
			jets.push_back(ROOT.TLorentzVector(sel.jet_pt[k], sel.jet_eta[k], sel.jet_phi[k], sel.jet_e[k]))
		w = helpers.wrapperC.getQCDWeight(nBtag, isBoosted, met, l, lisTight, jets, lsd0, isElectron, muonTrigger, topoetcone20, runNumber)
		return w

	def selectChannel(self, sel, syst, w):
		# OR all channels in the comma-separated list
		mapSel = {
					'be': 'bejets_2015,bejets_2016',
					'bmu': 'bmujets_2015,bmujets_2016',
					're': 'rejets_2015,rejets_2016',
					'rmu': 'rmujets_2015,rmujets_2016',
					'be0': 'bejets_2015,bejets_2016',
					'bmu0': 'bmujets_2015,bmujets_2016',
					're0': 'rejets_2015,rejets_2016',
					'rmu0': 'rmujets_2015,rmujets_2016',
					'be1': 'bejets_2015,bejets_2016',
					'bmu1': 'bmujets_2015,bmujets_2016',
					're1': 'rejets_2015,rejets_2016',
					'rmu1': 'rmujets_2015,rmujets_2016',
					'be2': 'bejets_2015,bejets_2016',
					'bmu2': 'bmujets_2015,bmujets_2016',
					're2': 'rejets_2015,rejets_2016',
					'rmu2': 'rmujets_2015,rmujets_2016',
					'be3': 'bejets_2015,bejets_2016',
					'bmu3': 'bmujets_2015,bmujets_2016',
					're3': 'rejets_2015,rejets_2016',
					'rmu3': 'rmujets_2015,rmujets_2016',
					'be2015': 'bejets_2015',
					'bmu2015': 'bmujets_2015',
					're2015': 'rejets_2015',
					'rmu2015': 'rmujets_2015',
					'be2016': 'bejets_2016',
					'bmu2016': 'bmujets_2016',
					're2016': 'rejets_2016',
					'rmu2016': 'rmujets_2016'}
		passSel = {}
		for i in mapSel:
			passSel[i] = True
			listSel = mapSel[i].split(',')

			passORChannels = False
			for item in listSel:
				if self.applyQCD == "e" and "mujets" in item:
					continue
				elif self.applyQCD == "mu" and "ejets" in item:
					continue
				hardPass = True
				if 'be' in i or 'bmu' in i:
					hardPass = False
					for k in range(0, len(sel.ljet_pt)):
						if sel.ljet_good[k] and sel.ljet_pt[k] > 300e3:
							hardPass = True
							break

				passChannel = getattr(sel, item)
				if passChannel and hardPass:
					passORChannels = True
					break
			passSel[i] = passORChannels

		if not passSel[self.ch]:
			return False

		# bugfix for no large-R jet requirement in QCD
		if 'be' in self.ch or 'bmu' in self.ch:
			passes = False
			for i in range(0, len(sel.ljet_pt)):
				if sel.ljet_good[i]:
					passes = True
					break
			if not passes:
				return False

		# veto resolved event if it passes the boosted channel
		if 're' in self.ch or 'rmu' in self.ch:
			#if sel.bejets or sel.bmujets or sel.bmujetsmet80:
			if passSel['be'] or passSel['bmu']:
				return False

		# apply b-tagging cut
		nBtag = 0
		for bidx in range(0, len(sel.tjet_mv2c10)):
			pb = ROOT.TLorentzVector(sel.tjet_pt[bidx], sel.tjet_eta[bidx], sel.tjet_phi[bidx], sel.tjet_e[bidx])
			if pb.Perp() > 10e3 and math.fabs(pb.Eta()) < 2.5 and sel.tjet_numConstituents[bidx] >= 2:
				if sel.tjet_mv2c10[bidx] > helpers.MV2C10_CUT:
					nBtag += 1
		#for bidx in range(0, len(sel.jet_mv2c10)):
		#	if sel.jet_mv2c10[bidx] > 0.8244273:
		#		nBtag += 1
		if nBtag < 1:
			return False

		if self.ch in ['be0', 'bmu0', 're0', 'rmu0']:
			if sel.Btagcat != 0:
				return False
		if self.ch in ['be1', 'bmu1', 're1', 'rmu1']:
			if sel.Btagcat != 1:
				return False
		if self.ch in ['be2', 'bmu2', 're2', 'rmu2']:
			if sel.Btagcat != 2:
				return False
		if self.ch in ['be3', 'bmu3', 're3', 'rmu3']:
			if sel.Btagcat != 3:
				return False

		# veto events in nominal ttbar overlapping with the mtt sliced samples
		# commented now as it is not available in mc15c
		if sel.mcChannelNumber == 410000 and hasattr(sel, "MC_ttbar_beforeFSR_m") and not self.noMttSlices:
			if sel.MC_ttbar_beforeFSR_m > 1.1e6:
				return False
		return True


	def run(self, sel, syst, wo, wTruth):
		if sel.mcChannelNumber in helpers.listWjets22:
			flag = sel.Wfilter_Sherpa_nT
			if self.keep == 'bb':
				if flag != 3 and flag != 4:
					return
			if self.keep == 'cc':
				if flag != 1:
					return
			if self.keep == 'c':
				if flag != 2:
					return
			if self.keep == 'l':
				if flag != 5:
					return

		w = wo
		if not self.selectChannel(sel, syst, w):
			return


		if self.applyQCD:
			w *= self.qcdWeight(sel)
		if self.applyMET > 0 and not ('be' in self.ch or 'bmu' in self.ch):
			if sel.met_met*1e-3 < self.applyMET:
				return

		if sel.mcChannelNumber != 0 and hasattr(sel, "MC_ttbar_beforeFSR_m"):
                        w0 = w/self.w2HDM
			self.h["trueMtt"][syst].Fill(sel.MC_ttbar_beforeFSR_m*1e-3, w)
			self.h["trueMttr"][syst].Fill(sel.MC_ttbar_beforeFSR_m*1e-3, w0*(self.w2HDM-1.))
                        self.h["trueMttLow1"][syst].Fill(sel.MC_ttbar_beforeFSR_m*1e-3, w)
                        self.h["trueMttLow1r"][syst].Fill(sel.MC_ttbar_beforeFSR_m*1e-3, w0*(self.w2HDM-1.))
                        self.h["trueMttLow2"][syst].Fill(sel.MC_ttbar_beforeFSR_m*1e-3, w)
                        self.h["trueMttLow2r"][syst].Fill(sel.MC_ttbar_beforeFSR_m*1e-3, w0*(self.w2HDM-1.))
		self.h["yields"][syst].Fill(1, w)
		self.h["runNumber"][syst].Fill(sel.runNumber, w)
		l = ROOT.TLorentzVector()
		lQ = 0
		if len(sel.el_pt) == 1:
			l.SetPtEtaPhiE(sel.el_pt[0], sel.el_eta[0], sel.el_phi[0], sel.el_e[0])
			lQ = sel.el_charge[0]
		elif len(sel.mu_pt) == 1:
			l.SetPtEtaPhiE(sel.mu_pt[0], sel.mu_eta[0], sel.mu_phi[0], sel.mu_e[0])
			lQ = sel.mu_charge[0]
		if lQ > 0:
			self.h["yieldsPos"][syst].Fill(1, w)
		elif lQ < 0:
			self.h["yieldsNeg"][syst].Fill(1, w)
		self.h["lepPt"][syst].Fill(l.Perp()*1e-3, w)
		self.h["lepEta"][syst].Fill(l.Eta(), w)
		self.h["lepPhi"][syst].Fill(l.Phi(), w)
		self.h["MET_phi"][syst].Fill(sel.met_phi, w)
		self.h["MET"][syst].Fill(sel.met_met*1e-3, w)
		self.h["nJets"][syst].Fill(len(sel.jet_pt), w)
		closestJetIdx = -1
		closestJetPt = 0
		closestJetDr = 99
		for i in range(0, len(sel.jet_pt)):
			cj = ROOT.TLorentzVector()
			cj.SetPtEtaPhiE(sel.jet_pt[i], sel.jet_eta[i], sel.jet_phi[i], sel.jet_e[i])
			dy = (cj.Rapidity() - l.Rapidity())
			dp = cj.DeltaPhi(l)
			dr = (dy**2 + dp**2)**0.5
			if dr < closestJetDr:
				closestJetIdx = i
				closestJetDr = dr
				closestJetPt = cj.Perp()

		self.h["closestJetDr"][syst].Fill(closestJetDr, w)
		self.h["closestJetPt"][syst].Fill(closestJetPt*1e-3, w)
		nBtags = 0
		tjets = []
		tb = []
		for bidx in range(0, len(sel.tjet_mv2c10)):
			pb = ROOT.TLorentzVector(sel.tjet_pt[bidx], sel.tjet_eta[bidx], sel.tjet_phi[bidx], sel.tjet_e[bidx])
			if pb.Perp() > 10e3 and math.fabs(pb.Eta()) < 2.5 and sel.tjet_numConstituents[bidx] >= 2:
				tjets.append(pb)
				b = False
				if sel.tjet_mv2c10[bidx] > helpers.MV2C10_CUT:
					nBtags += 1
					b = True
				tb.append(b)
		self.h["nTrkBtagJets"][syst].Fill(nBtags, w)
		self.h["mwt"][syst].Fill(math.sqrt(2*l.Perp()*sel.met_met*(1 - math.cos(l.Phi() - sel.met_phi)))*1e-3, w)
		self.h["mu"][syst].Fill(sel.mu, w)
		self.h["npv"][syst].Fill(sel.npv, w)
		self.h["vtxz"][syst].Fill(sel.vtxz, w)
		##TLorentzVector getNu(TLorentzVector l, double met, double met_phi) {
		##double getMtt(TLorentzVector lep, std::vector<TLorentzVector> jets, std::vector<bool> btag, TLorentzVector met) {
		nu = helpers.wrapperC.getNu(l, sel.met_met, sel.met_phi)
		if 'be' in self.ch or 'bmu' in self.ch:
			goodJetIdx = -1
			closeJetIdx = -1
			for i in range(0, len(sel.jet_pt)):
				if sel.jet_closeToLepton[i]:
					break
				#pt = 0
				#cj = ROOT.TLorentzVector()
				#cj.SetPtEtaPhiM(sel.jet_pt[i], sel.jet_eta[i], sel.jet_phi[i], sel.jet_m[i])
				#if cj.DeltaR(l) < 1.5:
				#	if cj.Perp() > pt:
				#		closeJetIdx = i
				#		pt = cj.Perp()
			for i in range(0, len(sel.ljet_pt)):
				if sel.ljet_good[i]:
					goodJetIdx = i
					break
			lj = ROOT.TLorentzVector()
			lj.SetPtEtaPhiM(sel.ljet_pt[goodJetIdx], sel.ljet_eta[goodJetIdx], sel.ljet_phi[goodJetIdx], sel.ljet_m[goodJetIdx])
			closeJet = ROOT.TLorentzVector()
			closeJet.SetPtEtaPhiE(sel.jet_pt[closeJetIdx], sel.jet_eta[closeJetIdx], sel.jet_phi[closeJetIdx], sel.jet_e[closeJetIdx])
			
			self.h["largeJetPt"][syst].Fill(lj.Perp()*1e-3, w)
			self.h["largeJetM"][syst].Fill(lj.M()*1e-3, w)
			self.h["largeJetEta"][syst].Fill(lj.Eta(), w)
			self.h["largeJetPhi"][syst].Fill(lj.Phi(), w)
	  		self.h["largeJet_tau32_wta"][syst].Fill(sel.ljet_tau32_wta[goodJetIdx], w)
			self.h["largeJet_tau21_wta"][syst].Fill(sel.ljet_tau21_wta[goodJetIdx], w)
			self.h["mtlep_boo"][syst].Fill((closeJet+nu+l).M()*1e-3, w)
			self.h["mtt"][syst].Fill((closeJet+nu+l+lj).M()*1e-3, w)
                        w0 = w/self.w2HDM
			self.h["mttr"][syst].Fill((closeJet+nu+l+lj).M()*1e-3, w0*(self.w2HDM-1.))
                        self.h["mttLow1"][syst].Fill((closeJet+nu+l+lj).M()*1e-3, w)
                        self.h["mttLow1r"][syst].Fill((closeJet+nu+l+lj).M()*1e-3, w0*(self.w2HDM-1.))
                        self.h["mttLow2"][syst].Fill((closeJet+nu+l+lj).M()*1e-3, w)
                        self.h["mttLow2r"][syst].Fill((closeJet+nu+l+lj).M()*1e-3, w0*(self.w2HDM-1.))
			if lQ > 0:
				self.h["mttPos"][syst].Fill((closeJet+nu+l+lj).M()*1e-3, w)
			elif lQ < 0:
				self.h["mttNeg"][syst].Fill((closeJet+nu+l+lj).M()*1e-3, w)
			self.h["largeJetPtMtt"][syst].Fill(lj.Perp()/(closeJet+nu+l+lj).M(), w)
			mtt = (closeJet+nu+l+lj).M()*1e-3
			if mtt > 1.1e3:
				self.h["btagSF"][syst].Fill(helpers.applyBtagSFFromFile(sel, syst), w)
		elif 're' in self.ch or 'rmu' in self.ch:
			jets = ROOT.vector('TLorentzVector')()
			#jets = []
			btag = ROOT.vector('bool')()
			#btag = []
			for k in range(0, len(sel.jet_pt)):
				jets.push_back(ROOT.TLorentzVector(sel.jet_pt[k], sel.jet_eta[k], sel.jet_phi[k], sel.jet_e[k]))
				tagged = False
				for t in range(0, len(tjets)):
					dr = jets[k].DeltaR(tjets[t])
					if dr < 0.4 and tb[t]:
						tagged = True
					btag.push_back(tagged)
			met = ROOT.TLorentzVector(sel.met_met, 0, sel.met_phi, sel.met_met)
			res_info = helpers.wrapperC.getMtt(l, jets, btag, met)
			mtt = res_info["mtt"]
			mtl = res_info["mtl"]
			mth = res_info["mth"]
			mwh = res_info["mwh"]
			chi2 = res_info["chi2"]
                        w0 = w/self.w2HDM
			self.h["mtt"][syst].Fill(mtt*1e-3, w)
			self.h["mttr"][syst].Fill(mtt*1e-3, w0*(self.w2HDM-1.))
			self.h["mttLow1"][syst].Fill(mtt*1e-3, w)
			self.h["mttLow1r"][syst].Fill(mtt*1e-3, w0*(self.w2HDM-1.))
			self.h["mttLow2"][syst].Fill(mtt*1e-3, w)
			self.h["mttLow2r"][syst].Fill(mtt*1e-3, w0*(self.w2HDM-1.))
			if lQ > 0:
				self.h["mttPos"][syst].Fill(mtt*1e-3, w)
			elif lQ < 0:
				self.h["mttNeg"][syst].Fill(mtt*1e-3, w)
			if mtt > 1.1e3:
				self.h["btagSF"][syst].Fill(helpers.applyBtagSFFromFile(sel, syst), w)
			self.h["mtlep_res"][syst].Fill(mtl*1e-3, w)
			self.h["mthad_res"][syst].Fill(mth*1e-3, w)
			self.h["mwhad_res"][syst].Fill(mwh*1e-3, w)
			self.h["chi2"][syst].Fill(chi2, w)
	def end(self):
		#print "Yield for channel ", self.ch, self.h["yields"][""].GetBinContent(1)
		self.write()

