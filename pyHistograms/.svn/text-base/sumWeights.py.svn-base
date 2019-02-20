#!/usr/bin/env python

import math
from ROOT import *
from optparse import OptionParser

def addFilesInChain(c, txtFileOption):
    for f in txtFileOption.split(','):
        txtf = open(f)
        for l in txtf.readlines():
	    if l[-1] == '\n':
	        l = l[0:-1]
            c.Add(l)

pdfList = ['PDF4LHC15_nlo_30']

def main():

	parser = OptionParser()
	parser.add_option("-f", "--files",
							 #dest="files", default="/nfs/dust/atlas/user/danilo/hists_sr2416/inputSW_all.txt,/nfs/dust/atlas/user/danilo/hists_sr2416/inputSW_syst.txt,/nfs/dust/atlas/user/danilo/hists_sr2416/inputSW_pdf.txt",
							 dest="files", default="inputSW_all.txt,inputSW_syst.txt,inputSW_pdf.txt",
				  help="Text file with list of input files.", metavar="FILELIST")
	parser.add_option("-t", "--types",
							 dest="types", default=",syst,pdf",
				  help="Categories they should be added in.", metavar="LIST")
	parser.add_option("-s", "--suffix",
							 dest="suffix", default="_new",
				  help="Suffix to be added in the output file.", metavar="SUFFIX")

	(options, args) = parser.parse_args()

	suf = options.suffix
	files = options.files.split(',')
	types = options.types.split(',')
	for idx in range(0, len(files)):
		f = files[idx]
		t = types[idx]

		if t == "pdf":
			pdfSumOfWeights = {} # map of DSID to map of PDF variation names to sum of weights
			t_pdfSumWeights = TChain("PDFsumWeights")
			addFilesInChain(t_pdfSumWeights, f)
			for k in range(0, t_pdfSumWeights.GetEntries()):
				t_pdfSumWeights.GetEntry(k)
				if not t_pdfSumWeights.dsid in pdfSumOfWeights:
					pdfSumOfWeights[t_pdfSumWeights.dsid] = {}
				for l in pdfList:
					if l not in pdfSumOfWeights[t_pdfSumWeights.dsid]:
						pdfSumOfWeights[t_pdfSumWeights.dsid][l] = []
					pdfAttr = getattr(t_pdfSumWeights, l)
					if len(pdfSumOfWeights[t_pdfSumWeights.dsid][l]) == 0:
						pdfSumOfWeights[t_pdfSumWeights.dsid][l] = [0]*len(pdfAttr)
					for m in range(0, len(pdfAttr)):
						pdfSumOfWeights[t_pdfSumWeights.dsid][l][m] += pdfAttr[m]

			pfs = open("sumOfWeights%s%s.txt" % (t, suf), "w")
			for channel in sorted(pdfSumOfWeights.keys()):
				for pdfName in sorted(pdfSumOfWeights[channel].keys()):
					for pdfNumber in range(0, len(pdfSumOfWeights[channel][pdfName])):
						pfs.write("%20d%20s%20s%20s%20d%20s%20f\n" % (channel, "", pdfName, "", pdfNumber, "", pdfSumOfWeights[channel][pdfName][pdfNumber]))
			pfs.close()
		else:
			sumOfWeights = {} # map of DSID to sum of weights
			t_sumWeights = TChain("sumWeights")
			addFilesInChain(t_sumWeights, f)

			for k in range(0, t_sumWeights.GetEntries()):
				t_sumWeights.GetEntry(k)
				if not t_sumWeights.dsid in sumOfWeights:
					sumOfWeights[t_sumWeights.dsid] = 0
				sumOfWeights[t_sumWeights.dsid] += t_sumWeights.totalEventsWeighted

			fs = open("sumOfWeights%s%s.txt" % (t, suf), "w")
			for channel in sorted(sumOfWeights.keys()):
				fs.write("%20d%20s%20f\n" % (channel, "", sumOfWeights[channel]))
			fs.close()


if __name__ == "__main__":
	main()

