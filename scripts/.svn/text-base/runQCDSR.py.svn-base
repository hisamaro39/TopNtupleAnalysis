
import HQTTtResonancesTools.DC15MC13TeV_25ns_EXOT4
import HQTTtResonancesTools.DC15Data13TeV_25ns_EXOT4

# input directory
ntuplesDir = '/afs/cern.ch/user/d/dferreir/work/eos/atlas/user/d/dferreir/topana/01122015v1'
#ntuplesDir = '/data/atlas/danilo/qcdcr_2341'
#ntuplesDir = '/data/atlas/danilo/qcdsr_2341'
ntuplesDir = '/data/atlas/danilo/qcdsr_ttva_2341'

# output directory
#outputDir = 'ttrescr'
outputDir = 'ttres'

# the default is AnaTtresSL, which produces many control pltos for tt res.
# The Mtt version produces a TTree to do the limit setting
# the QCD version aims at plots for QCD studies using the matrix method
# look into read.cxx to see what is available
# create yours, if you wish
analysisType='AnaTtresSL'
#analysisType='AnaTtresQCDfake'

# leave it for nominal to run only the nominal
systematics = 'nominal'
#systematics = 'nominal,EG_RESOLUTION_ALL__1down,EG_RESOLUTION_ALL__1up,EG_SCALE_ALL__1down,EG_SCALE_ALL__1up,JET_JER_SINGLE_NP__1up,JET_NPScenario1_JET_GroupedNP_1__1down,JET_NPScenario1_JET_GroupedNP_1__1up,JET_NPScenario1_JET_GroupedNP_2__1down,JET_NPScenario1_JET_GroupedNP_2__1up,JET_NPScenario1_JET_GroupedNP_3__1down,JET_NPScenario1_JET_GroupedNP_3__1up,MET_SoftTrk_ResoPara,MET_SoftTrk_ResoPerp,MET_SoftTrk_ScaleDown,MET_SoftTrk_ScaleUp,MUONS_ID__1down,MUONS_ID__1up,MUONS_MS__1down,MUONS_MS__1up,MUONS_SCALE__1down,MUONS_SCALE__1up'

btags = -1
#btags = 1

names   = []
#names  += ["Data15_13TeV_25ns_FS_EXOT4_3_3fb"]
names  += ["QCD15_13TeV_25ns_FS_EXOT4_3_3fb"]
# 25 ns datasets
#names  += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia']
#names  += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_mttsliced']
#names  += ['MC15_13TeV_25ns_FS_EXOT4_singletop']
#names  += ['MC15_13TeV_25ns_FS_EXOT4_Wjets']
#names  += ['MC15_13TeV_25ns_FS_EXOT4_Zjets']
#names  += ['MC15_13TeV_25ns_FS_EXOT4_VV']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime400']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime500']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime750']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1000']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1250']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1500']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1750']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2000']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2250']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2500']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2750']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime3000']
###names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime4000']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime5000']
#names += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarMCAtNLOHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythiaAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadLo', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadHi']

#import TopExamples.grid
import HQTTtResonancesTools.grid

import glob
files = glob.glob(ntuplesDir+'/*.root')
#samples = TopExamples.grid.Samples(names)
samples = HQTTtResonancesTools.grid.Samples(names)

import glob
import os
# get list of processed datasets
dirs = glob.glob(ntuplesDir+'/*')

# each "sample" below means an item in the list names above
# there may contain multiple datasets
# for each sample we want to read
for sample in samples:
    # write list of files to be read when processing this sample
    f = open(outputDir+"/input_"+sample.name+'.txt', 'w')
    # output file after running read
    #outfile = outputDir+"/"+sample.name
    outfile = sample.name
    
    # go over all directories in the ntuplesDir
    for d in dirs:
      # remove path and get only dir name in justfile
      justfile = d.split('/')[-1]
      if "group" in justfile.split('.')[0]:
          dsid_dir = justfile.split('.')[3] # get the DSID of the directory
      else:
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
          # go to the next directory in the same sample
          break
    f.close()
    theSysts = systematics
    isData = '0'
    qcdPar = ' --runMM 0 --loose 0'
    if "Data" in sample.name:
        theSysts = "nominal"
        isData = '1'
    elif "QCD" in sample.name:
        theSysts = "nominal"
        isData = '1'
        qcdPar = ' --runMM 1 --loose 1'
    if not "QCD" in sample.name:
        os.system('./read --removeOverlapHighMtt 0 --data '+isData+' '+qcdPar+' --btags '+str(btags)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+outfile+'.root,'+outputDir+'/resolved_mu_'+outfile+'.root,'+outputDir+'/boosted_e_'+outfile+'.root,'+outputDir+'/boosted_mu_'+outfile+'.root --systs '+theSysts)
    else:
        #os.system('./read --removeOverlapHighMtt 0 --data '+isData+' '+qcdPar+' --btags '+str(btags)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+outfile+'.root,'+outputDir+'/resolved_mu_'+outfile+'.root,'+outputDir+'/boosted_e_'+outfile+'.root,'+outputDir+'/boosted_mu_'+outfile+'.root --systs '+theSysts)
        onlyE = outfile.replace("QCD", "QCDe")
        onlyMu = outfile.replace("QCD", "QCDmu")
        os.system('./read --removeOverlapHighMtt 0 --onlyChannel 0,2 --data '+isData+' '+qcdPar+' --btags '+str(btags)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+onlyE+'.root,'+outputDir+'/resolved_mu_'+onlyE+'.root,'+outputDir+'/boosted_e_'+onlyE+'.root,'+outputDir+'/boosted_mu_'+onlyE+'.root --systs '+theSysts)
        os.system('./read --removeOverlapHighMtt 0 --onlyChannel 1,3 --data '+isData+' '+qcdPar+' --btags '+str(btags)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+onlyMu+'.root,'+outputDir+'/resolved_mu_'+onlyMu+'.root,'+outputDir+'/boosted_e_'+onlyMu+'.root,'+outputDir+'/boosted_mu_'+onlyMu+'.root --systs '+theSysts)
    #print('./read --removeOverlapHighMtt 0 --data '+isData+' '+qcdPar+' --btags '+str(btags)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+outfile+'.root,'+outputDir+'/resolved_mu_'+outfile+'.root,'+outputDir+'/boosted_e_'+outfile+'.root,'+outputDir+'/boosted_mu_'+outfile+'.root --systs '+theSysts)

