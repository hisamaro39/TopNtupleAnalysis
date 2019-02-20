#!/usr/bin/env python

import HQTTtResonancesTools.DC15MC13TeV_25ns_mc15c_EXOT4
import HQTTtResonancesTools.DC15Data13TeV_25ns_207_EXOT4

def main():
	# input directory
	ntuplesDir = '/nfs/dust/atlas/user/danilo/02062016WjetsCRv1'
	
	# output directory
	#outputDir = '/afs/desy.de/user/d/danilo/private/topana/Top246/TopNtupleAnalysis/pyHistograms/hists_WjetsMETCorrHF'
	outputDir = '/afs/desy.de/user/d/danilo/private/topana/Top246/TopNtupleAnalysis/pyHistograms/hists_WjetsHF3'

	# number of files per job
	nFilesPerJob = 80

	# use it to setup AnalysisTop
	rundir = '/afs/desy.de/user/d/danilo/private/topana/Top246'

	# email to use to tell us when the job is done
	email = 'dferreir@cern.ch'

	# queue to submit to
	#queue = '1nd'
	queue = 'default.q'

	# the default is AnaTtresSL, which produces many control pltos for tt res.
	# The Mtt version produces a TTree to do the limit setting
	# the QCD version aims at plots for QCD studies using the matrix method
	# look into read.cxx to see what is available
	# create yours, if you wish
	analysisType='AnaWjetsCR'
	
	# leave it for nominal to run only the nominal
	systematics = 'nominal'
	
	# 25 ns datasets
	names   = []
	names  += ['wbbjets']
	names  += ['wccjets']
	names  += ['wcjets']
	names  += ['wljets']
	names  += ['zjets']
	names  += ['tt']

	names  += ["data"]
	names  += ['singletop']
	names  += ['vv']

	mapToSamples = {
					'wbbjets': 'MC15c_13TeV_25ns_FS_EXOT4_Wjets22',
					'wccjets': 'MC15c_13TeV_25ns_FS_EXOT4_Wjets22',
					'wcjets': 'MC15c_13TeV_25ns_FS_EXOT4_Wjets22',
					'wljets': 'MC15c_13TeV_25ns_FS_EXOT4_Wjets22',
					'data': 'Data15_13TeV_25ns_207_EXOT4',
					'tt':'MC15c_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia',
					'singletop':'MC15c_13TeV_25ns_FS_EXOT4_singletop',
					'zjets':'MC15c_13TeV_25ns_FS_EXOT4_Zjets22',
					'vv': 'MC15c_13TeV_25ns_FS_EXOT4_VV',
		   }
	
	import TopExamples.grid
	
	
	import glob
	import os
	import sys
	# get list of processed datasets
	dirs = glob.glob(ntuplesDir+'/*')
	
	# each "sample" below means an item in the list names above
	# there may contain multiple datasets
	# for each sample we want to read
	for sn in names:
		sname = mapToSamples[sn]
		samples = TopExamples.grid.Samples([sname])
		sample = samples[0]

		forJobs = {}
		nJob = 0
		nFile = 0

		# write list of files to be read when processing this sample
		f = open(outputDir+"/input_"+sn+'.txt', 'w')
		# output file after running read
		outfile = sn
		
		# go over all directories in the ntuplesDir
		for d in dirs:
			# remove path and get only dir name in justfile
			justfile = d.split('/')[-1]
			dsid_dir = justfile.split('.')[2] # get the DSID of the directory
			# this will include all directories, so check if this director is in the sample
	
			# now go over the list of datasets in sample
			# and check if this DSID is there
			for s in sample.datasets:
				if len(s.split(':')) > 1:
					s = s.split(':')[1] # skip mc15_13TeV
				dsid_sample = s.split('.')[1] # get DSID
				if dsid_dir == dsid_sample: # this dataset belongs in the sample in the big for loop
					# get all files in the directory
					files = glob.glob(d+'/*.root*')
					# and write it in ht elist of input files to process
					for item in files:
						if not '.part' in item:
							f.write(item+'\n')
							if nFile == nFilesPerJob:
								nFile = 0
								nJob = nJob + 1
							if nFile == 0:
								forJobs[str(nJob)] = []
							forJobs[str(nJob)].append(item+'\n')
							nFile = nFile + 1
					# go to the next directory in the same sample
					break
		f.close()
		theSysts = systematics
		isData = ''
		extra = ''
		if "data" in sn:
			theSysts = "nominal"
			isData = ' -d '
		if "wbbjets" in sn:
			extra = ' --WjetsHF bb '
		elif 'wccjets' in sn:
			extra = ' --WjetsHF cc'
		elif 'wcjets' in sn:
			extra = ' --WjetsHF c'
		elif 'wljets' in sn:
			extra = ' --WjetsHF l'
	
		#os.system('./makeHistograms.py - '+isData+'  '+extra+' --files '+outputDir+"/input_"+sn+'.txt --analysis '+analysisType+' --output Wpre_resjets_el:'+outputDir+'/Wpre_resjets_el_'+outfile+'.root,Wpre_resjets_mu:'+outputDir+'/Wpre_resjets_mu_'+outfile+'.root,Wtag_resjets_el:'+outputDir+'/Wtag_resjets_el_'+outfile+'.root,Wtag_resjets_mu:'+outputDir+'/Wtag_resjets_mu_'+outfile+'.root --systs '+theSysts)

		for job in forJobs:
			jobName = sn+'_'+job
			infile = outputDir+"/input_"+jobName+'.txt'
			infullfile = outputDir+"/input_"+sn+'.txt'
			f = open(infile, 'w')
			for item in forJobs[job]:
				f.write(item)
			f.close()
			logfile = outputDir+"/stdout_"+jobName+'.txt'
			runfile = outputDir+"/run_"+jobName+'.sh'
			fr = open(runfile, 'w')
			fr.write('#!/bin/sh\n')
			fr.write('#$ -cwd\n')
			fr.write('#$ -j y\n')
			fr.write('#$ -l cvmfs\n')
			fr.write('#$ -l distro=sld6\n')
			fr.write('#$ -o '+logfile+'\n')
			fr.write('#$ -q '+queue+'\n')
			fr.write('#$ -m '+'eas'+'\n')
			fr.write('#$ -M '+email+'\n')
			fr.write('#$ -N tnapy_'+jobName+'\n')
			# bsub options
			#fr.write('#$ -e '+logfile+'\n')
			#fr.write('#$ -o '+logfile+'\n')
			#fr.write('#$ -q '+queue+'\n')
			#fr.write('#$ -N '+''+'\n')
			#fr.write('#$ -u '+email+'\n')
			#fr.write('#$ -J tnapy_'+jobName+'\n')
			fr.write('cd '+rundir+'\n')
			fr.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
			fr.write('export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK \n')
			fr.write('source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n')
			fr.write('lsetup rcsetup\n')
			fr.write('cd TopNtupleAnalysis/pyHistograms\n')
			#out = 'Wpre_resjets_el:'+outputDir+'/Wpre_resjets_el_'+jobName+'.root,Wpre_resjets_mu:'+outputDir+'/Wpre_resjets_mu_'+jobName+'.root,Wtag_resjets_el:'+outputDir+'/Wtag_resjets_el_'+jobName+'.root,Wtag_resjets_mu:'+outputDir+'/Wtag_resjets_mu_'+jobName+'.root'
			out = 'Wpre_resjets_el:'+outputDir+'/Wpre_resjets_el_'+jobName+'.root,Wpre_resjets_mu:'+outputDir+'/Wpre_resjets_mu_'+jobName+'.root,Wtag_resjets_el:'+outputDir+'/Wtag_resjets_el_'+jobName+'.root,Wtag_resjets_mu:'+outputDir+'/Wtag_resjets_mu_'+jobName+'.root,Wtag2_resjets_el:'+outputDir+'/Wtag2_resjets_el_'+jobName+'.root,Wtag2_resjets_mu:'+outputDir+'/Wtag2_resjets_mu_'+jobName+'.root'
			fr.write('./makeHistograms.py - '+isData+'   '+extra+'  --files '+infile+' --fullFiles '+infullfile+' --analysis '+analysisType+' --output '+out+'   --systs '+theSysts+'\n')
			fr.close()
			os.system('chmod a+x '+runfile)
			subcmd = 'qsub '+runfile
			os.system(subcmd)
			#sys.exit(0)
	
if __name__ == '__main__':
	main()

