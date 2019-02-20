#!/bin/sh
rm hist*
for ch in boosted 
do
  for lep in mu e 
  do
    #../../plotting/plot --mcOnly 1 --logY 0 -p boosted -c $lep -l 15 --normBinWidth 20 --yTitle "Events / 20 GeV" --saveTH1 boosted_$lep -h largeJetPt -C ./config_bkg.txt 
    #for hist in mtt largeJetM largeJet_tau32_wta mtlep_boo largeJetPt closeil_minDeltaR jet0_pt lepPt lepEta lepPhi mwt MET
    for hist in mtt 
    do
      ../../plotting/plot --normBinWidth 80 --mcOnly 1 --logY 0 -p $ch -c $lep -l 15 --saveTH1 $ch\_$lep -h $hist -C ./config_limit.txt 
    done
  done
done
