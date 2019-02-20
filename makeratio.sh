
S=1
LUMI=3.20905

#hadd -f boosted_e_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_all.root boosted_e_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia.root boosted_e_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_mttsliced.root
#hadd -f boosted_mu_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_all.root boosted_mu_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia.root boosted_mu_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_mttsliced.root

#hadd -f resolved_e_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_all.root resolved_e_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia.root resolved_e_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_mttsliced.root
#hadd -f resolved_mu_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_all.root resolved_mu_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia.root resolved_mu_MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythia_mttsliced.root


#hist=lepPt

for hist in largeJetM largeJetEta largeJetPhi largeJetSd12 largeJet_tau32 largeJet_tau32_wta largeJet_tau21 largeJet_tau21_wta chi2 mwhad_res mthad_res mtlep_res jet0_m jet1_pt jet1_m jet2_pt jet2_m closejl_pt  closejl_minDeltaR MET_phi nBtagJets leadbJetPt lepEta lepPhi mwt  mtt nTrkBtagJets ; do
  for ch in boosted ; do
    ../plotting/plotChannelRatio -c e/mu -p $ch -h $hist -l $LUMI -T $ch --smoothen 1 --yTitle "e/#mu" -C config.txt
  done
done

for hist in lepPt MET largeJetPt mtlep_boo jet0_pt closeJetPt ; do
for ch in boosted ; do
    ../plotting/plotChannelRatio -c e/mu -p $ch -h $hist -l $LUMI -T $ch --smoothen 1 --yTitle "e/#mu" -C config.txt
done
done

for ch in boosted ; do
    ../plotting/plotChannelRatio -c e/mu -p $ch -h yields -l $LUMI -T $ch --smoothen 1 -C config.txt >yields_ratio_${ch}.txt
done

for ch in boosted ; do
  ../plotting/plotChannelRatio -c e/mu -p $ch -h mtt -l $LUMI -T $ch  --smoothen 1 --yTitle "e/#mu" -C config.txt
done

