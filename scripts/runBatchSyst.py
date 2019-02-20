
import HQTTtResonancesTools.DC15MC13TeV_25ns_EXOT4
import HQTTtResonancesTools.DC15Data13TeV_25ns_EXOT4
from subprocess import Popen,PIPE

# lxbatch queue to submit to
queue = '1nd'

# email to use to tell us when the job is done
email = 'dferreir@cern.ch'

# directory with the base of the RootCore stuff
rundir = '/afs/cern.ch/work/d/dferreir/private/topana/Top2337'

# number of files per job
nFilesPerJob = 2

# input directory
ntuplesDir = '/eos/atlas/user/d/dferreir/topana/01122015TT'

#eosrun = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine.user/bin/eos.select'
eosrun='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'
entry = 'root://eosatlas.cern.ch/'

# output directory
outputDir = rundir+'/TopNtupleAnalysis/ttres33'


# the default is AnaTtresSL, which produces many control pltos for tt res.
# The Mtt version produces a TTree to do the limit setting
# the QCD version aims at plots for QCD studies using the matrix method
# look into read.cxx to see what is available
# create yours, if you wish
#analysisType='AnaTtresSLMtt'
analysisType='AnaTtresSL'

# leave it for nominal to run only the nominal
systematics = 'nominal'

# set to 1 to run the loose selection for QCD
loose = 0

# number of btags (negative for track jet btagging)
btags = -1

names   = []
names += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarMCAtNLOHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythiaAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadLo', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadHi']
names += ['MC15_13TeV_25ns_FS_EXOT4_ttbarSherpaAF2']

import TopExamples.grid

import glob
#files = glob.glob(ntuplesDir+'/*.root')
samples = TopExamples.grid.Samples(names)

import glob
import os
# get list of processed datasets
#dirs = glob.glob(ntuplesDir+'/*')
dirs = Popen([eosrun, "ls", ntuplesDir], stdout=PIPE).communicate()[0].split('\n')

# each "sample" below means an item in the list names above
# there may contain multiple datasets
# for each sample we want to read
for sample in samples:
    # write list of files to be read when processing this sample
    f = open(outputDir+"/input_"+sample.name+'.txt', 'w')
    # output file after running read
    #outfile = outputDir+"/"+sample.name
    outfile = sample.name

    forJobs = {}
    nJob = 0
    nFile = 0
    
    # go over all directories in the ntuplesDir
    for d in dirs:
      if d == '':
         continue

      # remove path and get only dir name in justfile
      #justfile = d.split('/')[-1]
      #dsid_dir = justfile.split('.')[2] # get the DSID of the directory
      dsid_dir = d.split('.')[2] # get the DSID of the directory
      # this will include all directories, so check if this director is in the sample

      # now go over the list of datasets in sample
      # and check if this DSID is there
      for s in sample.datasets:
        if len(s.split(':')) > 1:
          s = s.split(':')[1] # skip mc15_13TeV
        dsid_sample = s.split('.')[1] # get DSID
        if dsid_dir == dsid_sample: # this dataset belongs in the sample in the big for loop
          # get all files in the directory
          #files = glob.glob(d+'/*.root*')
          files = Popen([eosrun, "ls", ntuplesDir+"/"+d], stdout=PIPE).communicate()[0].split('\n')
          # and write it in ht elist of input files to process
          for item in files:
              if not '.part' in item and item != '':
                  #fullname = item
                  fullname = entry+ntuplesDir+"/"+d+'/'+item
                  f.write(fullname+'\n')

                  if nFile == nFilesPerJob:
                     nFile = 0
                     nJob = nJob + 1
                  if nFile == 0:
                     forJobs[str(nJob)] = []
                  forJobs[str(nJob)].append(fullname+'\n')
                  nFile = nFile + 1
          # go to the next directory in the same sample
          break
    f.close()
    theSysts = systematics
    isData = 0
    if "Data" in sample.name:
      theSysts = "nominal"
      isData = 1

    for job in forJobs:
        jobName = sample.name+'_'+job
        infile = outputDir+"/input_"+jobName+'.txt'
        infullfile = outputDir+"/input_"+sample.name+'.txt'
        f = open(infile, 'w')
        for item in forJobs[job]:
            f.write(item)
        f.close()
        errfile = outputDir+"/stderr_"+jobName+'.txt'
        logfile = outputDir+"/stdout_"+jobName+'.txt'
        runfile = outputDir+"/run_"+jobName+'.sh'
        fr = open(runfile, 'w')
        fr.write('#!/bin/sh\n')
        fr.write('cd '+rundir+'\n')
        fr.write('source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n')
        fr.write('lsetup rcsetup\n')
        fr.write('cd TopNtupleAnalysis\n')
        fr.write('ls\n')
        fr.write('./read --applyEWK 0 --removeOverlapHighMtt 1 --data '+str(isData)+' --btags '+str(btags)+' --loose '+str(loose)+' --files '+infile+' --fullFiles '+infullfile+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+outfile+'_'+job+'.root,'+outputDir+'/resolved_mu_'+outfile+'_'+job+'.root,'+outputDir+'/boosted_e_'+outfile+'_'+job+'.root,'+outputDir+'/boosted_mu_'+outfile+'_'+job+'.root --systs '+theSysts+'\n')
        fr.close()
        os.system('chmod a+x '+runfile)
        subcmd = 'bsub -e '+errfile+' -o '+logfile+' -q '+queue+' -N -u '+email+' -J tna_'+jobName+' '+runfile
        os.system(subcmd)
        #print(subcmd)
        #import sys
        #sys.exit(0)

