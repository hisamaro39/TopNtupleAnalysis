#!/usr/bin/env python

from helpers import *
import math
from ROOT import *
from analysis import *
from optparse import OptionParser

def main():

	parser = OptionParser()
	parser.add_option("-d", "--data",
							 action="store_true", dest="data", default=False,
							 help="Is this data?", metavar="BOOL")
	parser.add_option("-f", "--files",
							 dest="files", default="input.txt",
				  help="Text file with list of input files.", metavar="FILE")
	parser.add_option("-A", "--analysis",
							 dest="analysis", default="AnaTtresSL",
				  help="Analysis code to run.", metavar="ANALYSIS")
	parser.add_option("-o", "--output",
							 dest="output", default="re:hist_re.root,rmu:hist_rmu.root,be:hist_be.root,bmu:hist_bmu.root",
				  help="Comma-separated list of output ROOT files.", metavar="FILES")
	parser.add_option("-s", "--systs",
							 dest="systs", default="nominal",
				  help="Comma-separated list of systematic uncertainties in TTrees in the input file. Use 'all' to run over all the default ones.", metavar="SYSTEMATICS")
	parser.add_option("-W", "--WjetsHF",
							 dest="WjetsHF", default="all",
				  help="Which W+jets HF to keep. Can be all, bb, cc, c or l.", metavar="FLAVOURS")
	parser.add_option("-P", "--pdf",
							 dest="pdf", default="",
				  help="Which PDFs to reweight to.", metavar="PDFS")
	parser.add_option("-Q", "--qcd",
							 dest="qcd", default="False",
							 help="Apply QCD weights?", metavar="CHANNEL")
	parser.add_option("-N", "--noMttSlices",
							 action='store_true', dest="noMttSlices", default=False,
				  help="If set, stop vetoing high mtt events in 410000 sample.", metavar="BOOL")
	parser.add_option("-M", "--applyMET",
							 dest="applyMET", default="0",
				  help="Extra MET cut to be applied.", metavar="CUT")
        parser.add_option("-S", "--SCALAR",
                                                         dest="SCALAR", default="",
                                  help="Parameters to use when reweighting LO ttbar+jets to a scalar 2HDM setup with a configuration of {mH,mA,sin(b-a), tan(b) and the model type}.", metavar="MH,MA,SBA,TANB,TYPE")
	parser.add_option("-E", "--EFT",
							 dest="EFT", default="",
				  help="Parameters to use when reweighting LO ttbar to an EFT setup with a lambda and cvv configuration. Set lambda to a negative value to disable this.", metavar="LAMBDA,CVV")
	parser.add_option("-p", "--pdfForWeight",
							 dest="pdfForWeight", default="NNPDF30_nlo_as_0118", # this is the PDF used for LO, so it should be used for the alphaS
				  help="PDF to use to get alpha_S when doing either the EFT or the scalar model reweighting.", metavar="PDF")
	 
	(options, args) = parser.parse_args()
	 
	pdfList = options.pdf.split(',')
	Xsec = {}
	 
	sumOfWeights = {} # map of DSID to sum of weights
	# TODO
	pdfSumOfWeights = {} # map of DSID to map of PDF variation names to sum of weights
	
	if not options.data:
		if options.pdf != "":
			pfs = open("sumOfWeightspdf_new.txt")
			for line in pfs.readlines():
				line_spl = line.split()
				if not int(line_spl[0]) in pdfSumOfWeights:
					pdfSumOfWeights[int(line_spl[0])] = {}
				if not line_spl[1] in pdfSumOfWeights[int(line_spl[0])]:
					pdfSumOfWeights[int(line_spl[0])][line_spl[1]] = {}
				pdfSumOfWeights[int(line_spl[0])][line_spl[1]][int(line_spl[2])] = float(line_spl[3])
			pfs.close()

		fs = open("sumOfWeights_new.txt")
		for line in fs.readlines():
			line_spl = line.split()
			sumOfWeights[int(line_spl[0])] = float(line_spl[1])
		fs.close()
		fs = open("sumOfWeightssyst_new.txt")
		for line in fs.readlines():
			line_spl = line.split()
			sumOfWeights[int(line_spl[0])] = float(line_spl[1])
		fs.close()

	#print sumOfWeights
	 
	loadXsec(Xsec, "../scripts/XSection-MC15-13TeV-ttres.data")
	loadXsec(Xsec, "../../TopDataPreparation/data/XSection-MC15-13TeV.data")
	#loadXsec(Xsec, "../share/MC15c-SherpaWZ.data")

	# check if there is any W+jets sample there
	isWjets = False
	doEWK = False

	mt_load = TChain("nominal")
	addFilesInChain(mt_load, options.files, 20)
	ent = mt_load.GetEntries()
	for k in range(0, ent):
		mt_load.GetEntry(k)
		sel = readEvent(mt_load)
		if sel.mcChannelNumber in helpers.listWjets22:
			isWjets = True
		if sel.mcChannelNumber in helpers.listEWK:
			doEWK = True
		if isWjets or doEWK:
			break

	# systematics list
	if options.systs == 'all':
		systList = []
		systList.append('nominal')
		if isWjets:
			systList.append('wnorm__1up')
			systList.append('wnorm__1down')
			systList.append('wbb__1up')
			systList.append('wbb__1down')
			systList.append('wcc__1up')
			systList.append('wcc__1down')
			systList.append('wc__1up')
			systList.append('wc__1down')
			systList.append('wl__1up')
			systList.append('wl__1down')
		for i in range(0, 4):
			systList.append('btagbSF_'+str(i)+'__1up')
			systList.append('btagbSF_'+str(i)+'__1down')
			if i == 0:
				for j in range(1, 4):
					systList.append('btagbSF_'+str(i)+'_pt'+str(j)+'__1up')
					systList.append('btagbSF_'+str(i)+'_pt'+str(j)+'__1down')
		for i in range(0, 4):
			systList.append('btagcSF_'+str(i)+'__1up')
			systList.append('btagcSF_'+str(i)+'__1down')
			if i == 0:
				for j in range(1, 4):
					systList.append('btagcSF_'+str(i)+'_pt'+str(j)+'__1up')
					systList.append('btagcSF_'+str(i)+'_pt'+str(j)+'__1down')
		for i in range(0, 11):
			systList.append('btaglSF_'+str(i)+'__1up')
			systList.append('btaglSF_'+str(i)+'__1down')
			if i == 0:
				for j in range(1, 4):
					systList.append('btaglSF_'+str(i)+'_pt'+str(j)+'__1up')
					systList.append('btaglSF_'+str(i)+'_pt'+str(j)+'__1down')
		systList.append('btageSF_0__1up')
		systList.append('btageSF_0__1down')
		systList.append('btageSF_1__1up')
		systList.append('btageSF_1__1down')
		systList.extend(weightChangeSystematics)
		systList.remove('')
		if doEWK:
			systList.append('ttEWK__1up')
			systList.append('ttEWK__1down')
		if options.EFT != '':
			systList.append("eftScale__1up")
			systList.append("eftScale__1down")
		systematics  = 'EG_RESOLUTION_ALL__1down,EG_RESOLUTION_ALL__1up,EG_SCALE_ALL__1down,EG_SCALE_ALL__1up'
		systematics += ',JET_JER_SINGLE_NP__1up'
		# 3NP for the akt4 jets
		#systematics += ',JET_NPScenario1_JET_EtaIntercalibration_NonClosure__1down,JET_NPScenario1_JET_EtaIntercalibration_NonClosure__1up'
		#systematics += ',JET_NPScenario1_JET_GroupedNP_2__1down,JET_NPScenario1_JET_GroupedNP_2__1up,JET_NPScenario1_JET_GroupedNP_3__1down,JET_NPScenario1_JET_GroupedNP_3__1up,JET_NPScenario1_JET_GroupedNP_1__1up,JET_NPScenario1_JET_GroupedNP_1__1down'
		# 19 NP for the akt4 jets
		systematics += ',JET_19NP_JET_EffectiveNP_3__1down,JET_19NP_JET_EffectiveNP_6restTerm__1up,JET_19NP_JET_EffectiveNP_6restTerm__1down,JET_19NP_JET_Pileup_RhoTopology__1down,JET_19NP_JET_Pileup_OffsetNPV__1down,JET_19NP_JET_BJES_Response__1up,JET_19NP_JET_Pileup_RhoTopology__1up,JET_19NP_JET_EtaIntercalibration_TotalStat__1up,JET_19NP_JET_EffectiveNP_1__1up,JET_19NP_JET_EtaIntercalibration_NonClosure__1up,JET_19NP_JET_EffectiveNP_1__1down,JET_19NP_JET_BJES_Response__1down,JET_19NP_JET_Flavor_Response__1up,JET_19NP_JET_Pileup_OffsetMu__1up,JET_19NP_JET_EffectiveNP_4__1down,JET_19NP_JET_EffectiveNP_5__1up,JET_19NP_JET_EffectiveNP_5__1down,JET_19NP_JET_EffectiveNP_2__1down,JET_19NP_JET_PunchThrough_MC15__1down,JET_19NP_JET_EffectiveNP_2__1up,JET_19NP_JET_SingleParticle_HighPt__1up,JET_19NP_JET_EffectiveNP_3__1up,JET_19NP_JET_SingleParticle_HighPt__1down,JET_19NP_JET_Flavor_Composition__1up,JET_19NP_JET_Pileup_PtTerm__1down,JET_19NP_JET_PunchThrough_MC15__1up,JET_19NP_JET_Flavor_Response__1down,JET_19NP_JET_EtaIntercalibration_Modelling__1down,JET_19NP_JET_Flavor_Composition__1down,JET_19NP_JET_Pileup_PtTerm__1up,JET_19NP_JET_EtaIntercalibration_Modelling__1up,JET_19NP_JET_EffectiveNP_4__1up,JET_19NP_JET_Pileup_OffsetMu__1down,JET_19NP_JET_EtaIntercalibration_TotalStat__1down,JET_19NP_JET_Pileup_OffsetNPV__1up,JET_19NP_JET_EtaIntercalibration_NonClosure__1down'
		systematics += ',MET_SoftTrk_ResoPara,MET_SoftTrk_ResoPerp,MET_SoftTrk_ScaleDown,MET_SoftTrk_ScaleUp'
		# muon scale and resolution
		systematics += ',MUONS_ID__1down,MUONS_ID__1up,MUONS_MS__1down,MUONS_MS__1up,MUONS_SCALE__1down,MUONS_SCALE__1up'
		# strong systematics for large-R jets
		systematics += ',LARGERJET_Strong_JET_Rtrk_Modelling_All__1up,LARGERJET_Strong_JET_Rtrk_Modelling_All__1down,LARGERJET_Strong_JET_Rtrk_Baseline_All__1down,LARGERJET_Strong_JET_Rtrk_Baseline_All__1up,LARGERJET_Strong_JET_Rtrk_Tracking_All__1down,LARGERJET_Strong_JET_Rtrk_Tracking_All__1up,LARGERJET_Strong_JET_Rtrk_TotalStat_All__1down,LARGERJET_Strong_JET_Rtrk_TotalStat_All__1up'
		# weak systematics for large-R jets
		#systematics += ',LARGERJET_Weak_JET_Rtrk_TotalStat_mass__1up,LARGERJET_Weak_JET_Rtrk_Modelling_mass__1down,LARGERJET_Weak_JET_Rtrk_Modelling_pT__1up,LARGERJET_Weak_JET_Rtrk_Modelling_Tau32__1up,LARGERJET_Weak_JET_Rtrk_Tracking_mass__1up,LARGERJET_Weak_JET_Rtrk_Baseline_Tau32__1up,LARGERJET_Weak_JET_Rtrk_TotalStat_Tau32__1up,LARGERJET_Weak_JET_Rtrk_Tracking_pT__1up,LARGERJET_Weak_JET_Rtrk_Modelling_Tau32__1down,LARGERJET_Weak_JET_Rtrk_Modelling_pT__1down,LARGERJET_Weak_JET_Rtrk_Tracking_pT__1down,LARGERJET_Weak_JET_Rtrk_Baseline_Tau32__1down,LARGERJET_Weak_JET_Rtrk_TotalStat_mass__1down,LARGERJET_Weak_JET_Rtrk_Baseline_mass__1up,LARGERJET_Weak_JET_Rtrk_Baseline_pT__1down,LARGERJET_Weak_JET_Rtrk_Tracking_Tau32__1up,LARGERJET_Weak_JET_Rtrk_Baseline_pT__1up,LARGERJET_Weak_JET_Rtrk_Modelling_mass__1up,LARGERJET_Weak_JET_Rtrk_Tracking_Tau32__1down,LARGERJET_Weak_JET_Rtrk_TotalStat_Tau32__1down,LARGERJET_Weak_JET_Rtrk_TotalStat_pT__1down,LARGERJET_Weak_JET_Rtrk_TotalStat_pT__1up,LARGERJET_Weak_JET_Rtrk_Baseline_mass__1down,LARGERJET_Weak_JET_Rtrk_Tracking_mass__1down'
		# medium systematics for large-R jets
		#systematics += ',LARGERJET_Medium_JET_Rtrk_Modelling_Tau32__1up,LARGERJET_Medium_JET_Rtrk_TotalStat_Kin__1down,LARGERJET_Medium_JET_Rtrk_Tracking_Tau32__1down,LARGERJET_Medium_JET_Rtrk_Modelling_Kin__1down,LARGERJET_Medium_JET_Rtrk_Tracking_Kin__1down,LARGERJET_Medium_JET_Rtrk_Baseline_Kin__1down,LARGERJET_Medium_JET_Rtrk_Baseline_Kin__1up,LARGERJET_Medium_JET_Rtrk_Tracking_Kin__1up,LARGERJET_Medium_JET_Rtrk_TotalStat_Tau32__1up,LARGERJET_Medium_JET_Rtrk_TotalStat_Tau32__1down,LARGERJET_Medium_JET_Rtrk_TotalStat_Kin__1up,LARGERJET_Medium_JET_Rtrk_Baseline_Tau32__1down,LARGERJET_Medium_JET_Rtrk_Modelling_Tau32__1down,LARGERJET_Medium_JET_Rtrk_Baseline_Tau32__1up,LARGERJET_Medium_JET_Rtrk_Tracking_Tau32__1up,LARGERJET_Medium_JET_Rtrk_Modelling_Kin__1up'
		systList.extend(systematics.split(','))
	elif options.systs == 'pdf':
		systList = []
		systList.append('nominal')
		for m in pdfList:
			nvar = len(pdfSumOfWeights[pdfSumOfWeights.keys()[0]][m])
			for k in range(0, nvar):
				systList.append('pdf_%s_%d' % (m, k))
	else:
		systList = options.systs.split(',')
	 
	# load analysis code
	histSuffixes = []
	for item in systList:
		if item == 'nominal':
			histSuffixes.append('')
		else:
			histSuffixes.append(item)
	channels = {}
	for k in options.output.split(','):
		channels[k.split(':')[0]] = k.split(':')[1]
	analysisCode = {}
	import analysis
	#print "Systematics: ", histSuffixes
	#print "To loop over: ", systList
	anaClass = getattr(analysis, options.analysis) 
	eftLambda = -1
	eftCvv = 0
        scalarMH   = -1
        scalarMA   = -1
        scalarSBA  = -999
        scalarTANB = -1
        scalarTYPE = -1
	if options.EFT != "":
		eftStr = options.EFT.split(",")
		eftLambda = float(eftStr[0])
		eftCvv = float(eftStr[1])
		helpers.wrapperC.initPDF(options.pdfForWeight)
		helpers.wrapperC.setEFT(eftLambda, eftCvv)
        if options.SCALAR != "":
                scalarStr  = options.SCALAR.split(",")
                scalarMH   = float(scalarStr[0])
                scalarMA   = float(scalarStr[1])
                scalarSBA  = float(scalarStr[2])
                scalarTANB = float(scalarStr[3])
                scalarTYPE = int(scalarStr[4])
                helpers.wrapperC.initPDF(options.pdfForWeight)
                helpers.init2HDM(scalarMH,scalarMA,scalarSBA,scalarTANB,scalarTYPE)
                print "2HDM setup: mH=%g, mA=%g, sba=%g, tanb=%g, type=%g" % (scalarMH, scalarMA, scalarSBA, scalarTANB, scalarTYPE)
	for k in channels:
		analysisCode[k] = anaClass(k, histSuffixes, channels[k])
		analysisCode[k].keep = options.WjetsHF
		analysisCode[k].applyQCD = False
		if options.qcd != "False":
			analysisCode[k].applyQCD = options.qcd
		if options.noMttSlices:
			analysisCode[k].noMttSlices = True
		if options.applyMET != "":
			analysisCode[k].applyMET = float(options.applyMET)
		if options.EFT != "":
			eftStr = options.EFT.split(",")
			analysisCode[k].eftLambda = eftLambda
			analysisCode[k].eftCvv = eftCvv
                if options.SCALAR != "":
                        scalarStr = options.SCALAR.split(",")
                        analysisCode[k].scalarMH   = scalarMH
                        analysisCode[k].scalarMA   = scalarMA
                        analysisCode[k].scalarSBA  = scalarSBA
                        analysisCode[k].scalarTANB = scalarTANB
                        analysisCode[k].scalarTYPE = scalarTYPE
		print k, analysisCode[k], channels[k]
	 
	for s in systList:
		# s is nominal, or the name of systematic
		treeName = s # systematic name is the same as the TTree name
		if treeName in weightChangeSystematics or 'btag' in treeName or 'wnorm' in treeName or 'wbb_' in treeName or 'wcc_' in treeName or 'wc_' in treeName or 'wl_' in treeName or 'ttEWK_' in treeName or 'pdf_' in treeName:
			treeName = 'nominal'
		if options.qcd != "False":
			treeName += '_Loose'
		mt = TChain(treeName)
		suffix = s
		if suffix == 'nominal':
			suffix = ''
		addFilesInChain(mt, options.files)

		ent = mt.GetEntries()
		for k in range(0, ent):
			mt.GetEntry(k)
			if k % 10000 == 0:
				print "(tree = ",treeName,", syst = ",suffix,") Entry ", k, "/", mt.GetEntries()
			sel = readEvent(mt)

			# common part of the weight
			weight = 1
			if not options.data:
				weight *= sel.weight_mc
				channel = sel.mcChannelNumber
				weight *= Xsec[channel]
				if not 'pdf_' in suffix:
					if not channel in sumOfWeights:
						print "Could not find DSID ",channel, " in sum of weights."
						weight = 0
					else:
						weight /= sumOfWeights[channel]
				else:
					pdfName = (suffix.split('_', 1)[1]).rsplit('_', 1)[0]
					pdfNumber = int(suffix.rsplit('_', 1)[1])
					if not channel in pdfSumOfWeights:
						print "Could not find DSID ",channel, " in sum of weights."
						weight = 0
					else:
						weight /= pdfSumOfWeights[channel][pdfName][pdfNumber]
					

			for ana in analysisCode:
				weight_reco = analysisCode[ana].getWeight(sel, suffix)
				if 'pdf_' in suffix:
					pdfName = (suffix.split('_', 1)[1]).rsplit('_', 1)[0]
					pdfNumber = int(suffix.rsplit('_', 1)[1])
					pdfAttr = getattr(sel, pdfName)
					weight_reco *= pdfAttr[pdfNumber]
				analysisCode[ana].run(sel, suffix, weight*weight_reco, weight)
	 
	for k in analysisCode:
		analysisCode[k].end()

if __name__ == "__main__":
	main()

