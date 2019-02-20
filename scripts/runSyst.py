
import HQTTtResonancesTools.DC15MC13TeV_25ns_EXOT4
import HQTTtResonancesTools.DC15Data13TeV_25ns_EXOT4

# input directory
#ntuplesDir = '/afs/cern.ch/user/d/dferreir/work/cernbox/user/d/dferreir/topana/20102015'
#ntuplesDir = '/afs/cern.ch/user/d/dferreir/work/eos/atlas/user/d/dferreir/topana/24102015mass'
ntuplesDir = '/afs/cern.ch/user/d/dferreir/work/eos/atlas/user/d/dferreir/topana/23102015TTSyst'

# output directory
#outputDir = 'test'
#outputDir = 'test25ns'
#outputDir = 'test_mttold_calo'
#outputDir = 'test_mttnew_calo'
#outputDir = 'test_new_track'
outputDir = 'res_calo_batch'
#outputDir = 'massonly_calo'

# the default is AnaTtresSL, which produces many control pltos for tt res.
# The Mtt version produces a TTree to do the limit setting
# the QCD version aims at plots for QCD studies using the matrix method
# look into read.cxx to see what is available
# create yours, if you wish
#analysisType='AnaTtresSLMtt'
analysisType='AnaTtresSL'

btags = 1

# leave it for nominal to run only the nominal
systematics = 'nominal'
#systematics = 'nominal,EG_RESOLUTION_ALL__1down,EG_RESOLUTION_ALL__1up,EG_SCALE_ALL__1down,EG_SCALE_ALL__1up,JET_JER_SINGLE_NP__1up,JET_NPScenario1_JET_GroupedNP_1__1down,JET_NPScenario1_JET_GroupedNP_1__1up,JET_NPScenario1_JET_GroupedNP_2__1down,JET_NPScenario1_JET_GroupedNP_2__1up,JET_NPScenario1_JET_GroupedNP_3__1down,JET_NPScenario1_JET_GroupedNP_3__1up,MET_SoftTrk_ResoPara,MET_SoftTrk_ResoPerp,MET_SoftTrk_ScaleDown,MET_SoftTrk_ScaleUp,MUONS_ID__1down,MUONS_ID__1up,MUONS_MS__1down,MUONS_MS__1up,MUONS_SCALE__1down,MUONS_SCALE__1up'

# set to 1 to run the loose selection for QCD
loose = 0

names   = []
#names  += ["Data15_13TeV_25ns_FS_EXOT4_1_4fb"]
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
names += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarMCAtNLOHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythiaAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadLo', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadHi']

import TopExamples.grid

import glob
files = glob.glob(ntuplesDir+'/*.root')
samples = TopExamples.grid.Samples(names)

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
    isData = 0
    if "Data" in sample.name:
        theSysts = "nominal"
        isData = 1
    #os.system('./read --btags -1 --loose '+str(loose)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outfile+'_re.root,'+outfile+'_rmu.root,'+outfile+'_be.root,'+outfile+'_bmu.root --systs '+systematics)
    #os.system('./read --data '+isData+' --btags -1 --loose '+str(loose)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+outfile+'.root,'+outputDir+'/resolved_mu_'+outfile+'.root,'+outputDir+'/boosted_e_'+outfile+'.root,'+outputDir+'/boosted_mu_'+outfile+'.root --systs '+theSysts)
    os.system('./read --removeOverlapHighMtt 1 --data '+str(isData)+' --btags '+str(btags)+' --loose '+str(loose)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+outfile+'.root,'+outputDir+'/resolved_mu_'+outfile+'.root,'+outputDir+'/boosted_e_'+outfile+'.root,'+outputDir+'/boosted_mu_'+outfile+'.root --systs '+theSysts)
    #os.system('./read --btags 1 --loose '+str(loose)+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+outfile+'.root,'+outputDir+'/resolved_mu_'+outfile+'.root,'+outputDir+'/boosted_e_'+outfile+'.root,'+outputDir+'/boosted_mu_'+outfile+'.root --systs '+systematics)

