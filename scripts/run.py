import HQTTtResonancesTools.DC15MC13TeV_25ns_mc15c_EXOT4
import HQTTtResonancesTools.DC15Data13TeV_25ns_207_EXOT4

print "Initialize!!"
# Initialize
ntuplesDir = '/home/mtanaka/WorkSpace/data5/ttbar_trigger'
trthr='0'
bdtthr='100'
outputDir = 'result_trigger/default'
analysisType='AnaTtresSL'
drlj='15'

# only run boosted channel?
onlyBoosted = '1'
# Debug mode 1:on 0:off
Debug = '0'
# Number of events to run
NEvent=-1
# Analysis mode
mode= '7'
# consider systematics
considerSyst = '0'
# <0:tjet 
# >0:calo-jet
btags = -1
# dr cut for closest lepton and small-jet
# default 15
drlj='15'# *0.1 is applied
# b-tag Working Point
#possible option 70% 77% 85% 100%
bwp='100' # % 100 means no b-tag

####################################
# Analysis mode for validation
# -1: off
# 0: release20 validation default
# 1: release21 validation default
# 2: release21 use met trigger
# 3: release20 check btag & toptag
# 4: release21 check met trigger efficiency
# 5: release21 validation use met trigger
# 6: release21 validation use pufit met trigger
# 7: release21 check btag & toptag
# 8: release20 check VR track b-tag
# 9: release20 default
# 10: release20 check multi jets
# 11: release21 check multi jets
# 12: release20 validation use calo btag
# 13: release20 validation use calo hybrid btag
# 14: release21 use calo btag
# 15: release21 use calo DL1rnn btag
# 16: release21 use calo DL1 btag
# 17: release21 validation use calo btag
# 18: release21 validation use calo DL1 btag
# 19: release21 validation use calo DL1rnn btag
# 20: release20 use substructure top tagger
# 21: release20 use smoothed top tagger
# 22: release20 validation smoothed top tagger
# 23: release21 use smoothed top tagger
# 24: release21 use dnn top tagger
mode_validation = 1
###################################

####################################
# Analysis mode for met trigger and flat b-tag study
# -1: off
# 0: default
# 1: use met trigger
# 2: use flat b-tag
# 3: met trigger & flat b-tag
# 4: NO OR for track jets
# 5: NO OR for track jets & flat b-tag
# 6: met trigger & NO OR for track jets & flat b-tag
# 7: validation
# 8: use calo b-tag
# 9: use hyb calo b-tag
mode_metbtag = -1
###################################

####################################
# Analysis mode for checking trigger
# -1: off
# 0: default
mode_check_trigger = -1
###################################

####################################
# Analysis mode for checking data driven QCD
# -1: off
# 0: default
mode_multi_jets = -1
###################################

# output directory
if mode_metbtag == 0:
  print 'Mode: default'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/default'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '0' # analysis type
elif mode_metbtag == 1:
  print 'Mode: use met trigger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/met_trigger'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '1' # analysis type
elif mode_metbtag == 2:
  print 'Mode: use flat b-tag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/flat_btag'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '2' # analysis type
elif mode_metbtag == 3:
  print 'Mode: use met trigger & flat b-tag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/met_trigger_flat_btag'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '3' # analysis type
elif mode_metbtag == 4:
  print 'Mode: NO OR for track jets'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/noOR_trackjet'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '4' # analysis type
elif mode_metbtag == 5:
  print 'Mode: NO OR for track jets & flat b-tag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/noOR_trackjet_flat_btag'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '5' # analysis type
elif mode_metbtag == 6:
  print 'Mode: MET trigger & NO OR for track jets & flat b-tag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/met_trigger_noOR_trackjet_flat_btag'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '6' # analysis type
elif mode_metbtag == 7:
  print 'Mode: validation'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/validation'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '7' # analysis type
elif mode_metbtag == 8:
  print 'Mode: use calo b-tag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/calo_btag'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '8' # analysis type
elif mode_metbtag == 9:
  print 'Mode: use hyblid calo b-tag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data7/notrigger-nobtag'
  outputDir = 'result_trigger_btag/hyb_calo_btag'
  bwp='100'
  analysisType='AnaTtresMetBtag'
  mode = '9' # analysis type
elif mode_validation == 0:
  print 'Mode: release20 validation'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel20'
  outputDir = 'result_validation/rel20_validation'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '0' 
elif mode_validation == 1:
  print 'Mode: release21 validation'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel21'
  outputDir = 'result_validation/rel21_validation'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '1' 
elif mode_validation == 2:
  print 'Mode: release21 use met trigger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/notrigger_nobtag_rel21'
  outputDir = 'result_validation/rel21_use_met_trigger'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '2' 
elif mode_validation == 3:
  print 'Mode: release20 check btag & toptag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/check_tagging_rel20'
  outputDir = 'result_validation/rel20_check_tagging'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '3' 
elif mode_validation == 4:
  print 'Mode: release21 check met trigger efficiency'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel21'
  outputDir = 'result_validation/rel21_met_trigeff'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '4' 
elif mode_validation == 5:
  print 'Mode: release21 validation use met trigger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel21'
  outputDir = 'result_validation/rel21_validation_use_met_trigger'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '5' 
elif mode_validation == 6:
  print 'Mode: release21 validation use pufit met trigger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel21'
  outputDir = 'result_validation/rel21_validation_use_pufit_met_trigger'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '6' 
elif mode_validation == 7:
  print 'Mode: release21 check btag & toptag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/check_tagging_rel21'
  outputDir = 'result_validation/rel21_check_tagging'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '3' 
elif mode_validation == 8:
  print 'Mode: release20 check VR track b-tag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data8/VRbtag'
  outputDir = 'result_validation/rel20_check_vrbtag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '3' 
elif mode_validation == 9:
  print 'Mode: release20 default'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/notrigger_nobtag_rel20'
  outputDir = 'result_validation/rel20_default'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '9' 
elif mode_validation == 10:
  print 'Mode: release20 check multi jet'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/notrigger_nobtag_rel20'
  outputDir = 'result_validation/rel20_check_multijets'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '10' 
elif mode_validation == 11:
  print 'Mode: release21 check multi jet'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/notrigger_nobtag_rel21'
  outputDir = 'result_validation/rel21_check_multijets'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '10' 
elif mode_validation == 12:
  print 'Mode: release20 validation use calo btag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel20'
  outputDir = 'result_validation/rel20_validation_calo_btag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '12' 
elif mode_validation == 13:
  print 'Mode: release20 validation use hyblid calo btag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel20'
  outputDir = 'result_validation/rel20_validation_hyb_calo_btag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '13' 
elif mode_validation == 14:
  print 'Mode: release21 use calo btag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/rel21'
  outputDir = 'result_validation/rel21_calo_btag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '14' 
elif mode_validation == 15:
  print 'Mode: release21 use calo DL1rnn btag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/rel21'
  outputDir = 'result_validation/rel21_calo_dl1rnn_btag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '15' 
elif mode_validation == 16:
  print 'Mode: release21 use calo DL1 btag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/rel21'
  outputDir = 'result_validation/rel21_calo_dl1_btag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '16' 
elif mode_validation == 17:
  print 'Mode: release21 validation use calo btag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel21'
  outputDir = 'result_validation/rel21_validation_calo_btag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '12' 
elif mode_validation == 18:
  print 'Mode: release21 validation use DL1 calo btag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel21'
  outputDir = 'result_validation/rel21_validation_calo_dl1_btag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '18' 
elif mode_validation == 19:
  print 'Mode: release21 validation use DL1rnn calo btag'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel21'
  outputDir = 'result_validation/rel21_validation_calo_dl1rnn_btag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '19' 
elif mode_validation == 20:
  print 'Mode: release20 use substructure top tagger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/check_tagging_rel20'
  outputDir = 'result_validation/rel20_substructure_toptag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '20' 
elif mode_validation == 21:
  print 'Mode: release20 use smoothed top tagger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/check_tagging_rel20'
  outputDir = 'result_validation/rel20_smoothed_toptag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '21' 
elif mode_validation == 22:
  print 'Mode: release20 validation'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/validation_rel20'
  outputDir = 'result_validation/rel20_validation_smoothed_toptag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '22' 
elif mode_validation == 23:
  print 'Mode: release21 use smoothed top tagger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/rel21'
  outputDir = 'result_validation/rel21_smoothed_toptag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '23' 
elif mode_validation == 24:
  print 'Mode: release21 use DNN top tagger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/rel21'
  outputDir = 'result_validation/rel21_dnn_toptag'
  bwp='100'
  analysisType='AnaTtresValidation'
  mode = '24' 
elif mode_check_trigger == 0:
  print 'Mode: check trigger'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/check_trigger_rel21'
  outputDir = 'result_validation/rel21_check_trigger'
  bwp='100'
  analysisType='AnaTtresCheckTrigger'
  mode = '0' 
elif mode_multi_jets == 0:
  print 'Mode: check data driven QCD'
  ntuplesDir = '/home/mtanaka/WorkSpace/data5/check_multi_jets_rel20'
  outputDir = 'result_multi_jets/default'
  bwp='70'
  analysisType='AnaTtresMultiJets'
  mode = '0' 

if Debug == '1': 
  outputDir = 'test'

# the default is AnaTtresSL, which produces many control pltos for tt res.
# The Mtt version produces a TTree to do the limit setting
# the QCD version aims at plots for QCD studies using the matrix method
# look into read.cxx to see what is available
# create yours, if you wish

# leave it for nominal to run only the nominal
systematics = 'nominal'
#systematics = 'nominal,JET_21NP_JET_EffectiveNP_2__1down,JET_21NP_JET_EffectiveNP_2__1up,JET_21NP_JET_BJES_Response__1down,JET_21NP_JET_BJES_Response__1up,JET_21NP_JET_Pileup_PtTerm__1down,JET_21NP_JET_Pileup_PtTerm__1up,JET_21NP_JET_Pileup_OffsetNPV__1down,JET_21NP_JET_Pileup_OffsetNPV__1up,JET_21NP_JET_Pileup_RhoTopology__1up,JET_21NP_JET_Pileup_RhoTopology__1down,JET_21NP_JET_Pileup_OffsetMu__1down,JET_21NP_JET_Pileup_OffsetMu__1up,JET_21NP_JET_Flavor_Response__1down,JET_21NP_JET_Flavor_Response__1up,JET_21NP_JET_Flavor_Composition__1down,JET_21NP_JET_Flavor_Composition__1up,MET_SoftTrk_ScaleDown,MET_SoftTrk_ScaleUp,MUON_SCALE__1down,MUON_SCALE__1up'

names   = []
#names  += ["Data15_13TeV_25ns_207_EXOT4"]
#names += ["Data16_13TeV_25ns_33257ipb_EXOT4"]
#names += ["Data16_13TeV_25ns_EXOT4"]
#names += ["Data17_13TeV_25ns_EXOT4"]
#names  += ["QCD15_13TeV_25ns_FS_EXOT4_3_3fb"]
#names  += ["QCD16_13TeV_25ns_EXOT4"]
## 25 ns datasets
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia']
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_mttsliced']
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_singletop']
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_Wjets221']
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_Zjets221']
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_VV']
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_ttbarV']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime400']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime500']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime750']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime1000']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime1250']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime1500']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime1750']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime2000']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime2250']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime2500']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime2750']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime3000']
#names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime4000']
names += ['MC15c_13TeV_25ns_FS_EXOT4_Zprime5000']
#names += ["MC15_13TeV_25ns_FS_EXOT4_dijets"]
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_singletop_tchannel']
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_singletop_schannel']
#names  += ['MC15c_13TeV_25ns_FS_EXOT4_singletop_wt']
#names += ['htaumu']
##names += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarMCAtNLOHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythiaAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadLo', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadHi']
#
import TopExamples.grid

import glob
files = glob.glob(ntuplesDir+'/*.root')
samples = TopExamples.grid.Samples(names)

import glob
import os
os.system('mkdir -p result_validation')
os.system('mkdir -p result_multi_jets')
os.system('mkdir -p result_trigger_btag')
os.system('mkdir -p result_check_trigger')
os.system('mkdir -p '+outputDir)
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
    isData = '0'
    qcdPar = ''
    if "Data" in sample.name:
        theSysts = "nominal"
        isData = '1'
    elif "QCD" in sample.name:
        theSysts = "nominal"
        isData = '1'
        qcdPar = ' --runMM 1 --loose 1 '
    os.system('./read --removeOverlapHighMtt 1 --data '+isData+' '+qcdPar+' --bwp '+bwp+' --drlj '+drlj+' --bdtthr '+bdtthr+' --debug '+Debug+' --trthr '+trthr+' --boost '+onlyBoosted+' --btags '+str(btags)+ ' --nentries '+str(NEvent)+' --mode '+mode+' --files '+outputDir+"/input_"+sample.name+'.txt'+' --analysis '+analysisType+' --output '+outputDir+'/resolved_e_'+outfile+'.root,'+outputDir+'/resolved_mu_'+outfile+'.root,'+outputDir+'/boosted_e_'+outfile+'.root,'+outputDir+'/boosted_mu_'+outfile+'.root --systs '+theSysts)
