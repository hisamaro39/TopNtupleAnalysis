
import HQTTtResonancesTools.DC15MC13TeV_25ns_EXOT4
import HQTTtResonancesTools.DC15Data13TeV_25ns_EXOT4
from subprocess import Popen,PIPE
import os

# output directory
outdir = 'ttres32'

names   = []
# 25 ns datasets
names  += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia']
names  += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_mttsliced']
names  += ['MC15_13TeV_25ns_FS_EXOT4_singletop']
names  += ['MC15_13TeV_25ns_FS_EXOT4_Wjets']
names  += ['MC15_13TeV_25ns_FS_EXOT4_Zjets']
names  += ['MC15_13TeV_25ns_FS_EXOT4_VV']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime400']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime500']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime750']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1000']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1250']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1500']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1750']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2000']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2250']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2500']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2750']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime3000']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime4000']
names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime5000']

names = []
#names += ['MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarMCAtNLOHerwigAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythiaAF2', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadLo', 'MC15_13TeV_25ns_FS_EXOT4_ttbarRadHi']
names += ['MC15_13TeV_25ns_FS_EXOT4_ttbaraMcAtNlo_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime400_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime500_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime750_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1000_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1250_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1500_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime1750_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2000_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2250_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2500_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime2750_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime3000_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime4000_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_Zprime5000_PDF']
#names += ['MC15_13TeV_25ns_FS_EXOT4_ttbarSherpaAF2']

channels = ['resolved_e', 'resolved_mu', 'boosted_e', 'boosted_mu']
#channels = ['boosted_e', 'boosted_mu']

import glob

for sample in names:
    for ch in channels:
        cmd = 'hadd -f -k '+outdir+"/"+ch+"_"+sample+".root"+"  "
        files = glob.glob(outdir+"/"+ch+"_"+sample+"_[0-9]*.root")
        if len(files) == 0:
           print "Failed to find any result from chanel "+ch+" for sample "+sample
           continue
        for i in files:
           cmd = cmd + "  " + i

        print(cmd)
        os.system(cmd)


#for ch in channels:
#    os.system('hadd -f -k '+outdir+"/"+ch+'_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_all.root '+outdir+"/"+ch+'_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia.root '+outdir+"/"+ch+'_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_mttsliced.root')

