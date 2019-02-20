
S=1
#LUMI=3.31668
LUMI=3.20905

for ch in resolved boosted ; do
  for lep in e mu ; do
    ../plotting/plot -c $lep -p $ch -h largeJetM -l $LUMI --smoothen 0 --yTitle "Events / 10 GeV" --xMin 50 --yMax 850 --xTitle "Large-R jet mass [GeV]" --stamp 1
    ../plotting/plot -c $lep -p $ch -h mtlep_boo -l $LUMI --smoothen 0 --normBinWidth 20 --xMin 0 --yTitle "Events / 20 GeV" --xTitle "Mass of the leptonic top candidate [GeV]" --stamp 1
    ../plotting/plot -c $lep -p $ch -h MET -l $LUMI --smoothen 0 --normBinWidth 10 --xMin 20 --yTitle "Events / 10 GeV" --yMax 500 --xTitle "E_{T}^{miss} [GeV]" --stamp 1
    ../plotting/plot -c $lep -p $ch -h lepPt -l $LUMI --smoothen 0 --rebin 4 --xTitle "Lepton p_{T} [GeV]" --yTitle "Events / 20 GeV" --stamp 1
    ../plotting/plot -c $lep -p $ch -h closeJetPt -l $LUMI --smoothen 0 --normBinWidth 1 --yTitle "Events / GeV" --xTitle "Selected jet p_{T} [GeV]" --stamp 1
    ../plotting/plot -c $lep -p $ch -h largeJetPt -l $LUMI --smoothen 0 --normBinWidth 20 --yTitle "Events / 20 GeV" --xTitle "Large-R jet p_{T} [GeV]" --stamp 1
done
done



for hist in largeJetPt largeJetM largeJet_tau32_wta ; do
for ch in boosted resolved ; do
  for lep in e mu ; do
    #../plotting/plot -c $lep -p $ch -h $hist -l $LUMI --smoothen 1 --normBinWidth 0 -o ${ch}_${hist}_${lep}_smooth.pdf

    ../plotting/plotCompareNominal --syst LARGERJET_JET_Top_CrossCalib__1up,LARGERJET_JET_Top_CrossCalib__1down --systTitles "akt10 cross calib. up,akt10 cross calib. down" -l $LUMI -c $lep -p $ch -h $hist -o syst_${ch}_${hist}_${lep}_akt10xcalib.pdf --smoothen 1
    ../plotting/plotCompareNominal --syst LARGERJET_JET_Top_Run1__1up,LARGERJET_JET_Top_Run1__1down --systTitles "akt10 run 1 up,akt10 run 1 down" -l $LUMI -c $lep -p $ch -h $hist -o syst_${ch}_${hist}_${lep}_akt10run1.pdf --smoothen 1

    ../plotting/plotCompareNominal --syst CORRELATED_LARGERJET_CORR1__1up,CORRELATED_LARGERJET_CORR1__1down --systTitles "akt10 run 1 and akt4 JES up,akt10 run 1 and akt4 JES down" -l $LUMI -c $lep -p $ch -h $hist -o syst_${ch}_${hist}_${lep}_corr1.pdf --smoothen 1
    ../plotting/plotCompareNominal --syst CORRELATED_LARGERJET_CORR2A__1up,CORRELATED_LARGERJET_CORR2A__1down --systTitles "akt10 run 1 JES and akt4 JES up,akt10 run 1 JES and akt4 JES down" -l $LUMI -c $lep -p $ch -h $hist -o syst_${ch}_${hist}_${lep}_corr2a.pdf --smoothen 1
    ../plotting/plotCompareNominal --syst CORRELATED_LARGERJET_CORR2B__1up,CORRELATED_LARGERJET_CORR2B__1down --systTitles "akt10 run 1 JMS+Tau32 up,akt10 run 1 JMS+Tau32 down" -l $LUMI -c $lep -p $ch -h $hist -o syst_${ch}_${hist}_${lep}_corr2b.pdf --smoothen 1

  done
done
done

for ch in boosted resolved ; do
  for lep in e mu ; do

    ../plotting/plot -c $lep -p $ch -h yields -l $LUMI --smoothen 0 >yields_${ch}_${lep}.txt

    ../plotting/plotCompareNominal --syst btagbSF_0__1up,btagbSF_0__1down --systTitles "btag eff E0 up,btag eff E0 down" -l $LUMI -c $lep -p $ch -h leadTrkbJetPt -o syst_${ch}_leadTrkbJetPt_${lep}_btagb0.pdf --smoothen 0
    ../plotting/plotCompareNominal --syst btagcSF_0__1up,btagcSF_0__1down --systTitles "btag c mistag E0 up,btag c mistag E0 down" -l $LUMI -c $lep -p $ch -h leadTrkbJetPt -o syst_${ch}_leadTrkbJetPt_${lep}_btagc0.pdf --smoothen 0
    ../plotting/plotCompareNominal --syst btaglSF_0__1up,btaglSF_0__1down --systTitles "btag l mistag E0 up,btag l mistag E0 down" -l $LUMI -c $lep -p $ch -h leadTrkbJetPt -o syst_${ch}_leadTrkbJetPt_${lep}_btagl0.pdf --smoothen 0

    ../plotting/plotCompareNominal --syst btageSF_0__1up,btageSF_0__1down --systTitles "btag extrap up,btag extrap down" -l $LUMI -c $lep -p $ch -h leadTrkbJetPt -o syst_${ch}_leadTrkbJetPt_${lep}_btage0.pdf --smoothen 0

    # systematics in mtt
    #../plotting/plotCompareNominal --mcOnly 1 -c $lep -p $ch -h mtt -l $LUMI --other MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2,MC15_13TeV_25ns_FS_EXOT4_ttbarMCAtNLOHerwigAF2 --titles "Powheg+Herwig(AF2),MC@NLO+Herwig(AF2)" -o systmodel_${ch}_mtt_${lep}_mcgen.pdf --smoothen $S
    #../plotting/plotCompareNominal --mcOnly 1 -c $lep -p $ch -h mtt -l $LUMI --other MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythiaAF2,MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2 --titles "Powheg+Pythia(AF2),Powheg+Herwig(AF2)" -o systmodel_${ch}_mtt_${lep}_pshower.pdf --smoothen $S
    #../plotting/plotCompareNominal --mcOnly 1 -c $lep -p $ch -h mtt -l $LUMI --other MC15_13TeV_25ns_FS_EXOT4_ttbarRadLo,MC15_13TeV_25ns_FS_EXOT4_ttbarRadHi --titles "ISR/FSR(low),ISR/FSR(high)" -o systmodel_${ch}_mtt_${lep}_isrfsr.pdf --smoothen $S

    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 -c $lep -p $ch -h mtt -l $LUMI --other MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2,MC15_13TeV_25ns_FS_EXOT4_ttbarMCAtNLOHerwigAF2 --titles "MC gen up,MC gen dw" -o systmodel_${ch}_mtt_${lep}_mcgen.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 -c $lep -p $ch -h mtt -l $LUMI --other MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegPythiaAF2,MC15_13TeV_25ns_FS_EXOT4_ttbarPowhegHerwigAF2 --titles "p. shower up,p. shower dw" -o systmodel_${ch}_mtt_${lep}_pshower.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 -c $lep -p $ch -h mtt -l $LUMI --other MC15_13TeV_25ns_FS_EXOT4_ttbarRadLo,MC15_13TeV_25ns_FS_EXOT4_ttbarRadHi --titles "ISR/FSR(up),ISR/FSR(dw)" -o systmodel_${ch}_mtt_${lep}_isrfsr.pdf --smoothen $S

    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst JET_JER_SINGLE_NP__1up --systTitles "akt4 JER" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_akt4jer.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst JET_NPScenario1_JET_GroupedNP_1__1up,JET_NPScenario1_JET_GroupedNP_1__1down --systTitles "akt4 JES 1 up,akt4 JES 1 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_akt4jes1.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst JET_NPScenario1_JET_GroupedNP_2__1up,JET_NPScenario1_JET_GroupedNP_2__1down --systTitles "akt4 JES 2 up,akt4 JES 2 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_akt4jes2.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst JET_NPScenario1_JET_GroupedNP_3__1up,JET_NPScenario1_JET_GroupedNP_3__1down --systTitles "akt4 JES 3 up,akt4 JES 3 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_akt4jes3.pdf --smoothen $S

    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst MUONS_MS__1up,MUONS_MS__1down --systTitles "muon res. MS up,muon res. MS down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_muresms.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst MUONS_ID__1up,MUONS_ID__1down --systTitles "muon res. ID up,muon res. ID down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_muresid.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst MUONS_SCALE__1up,MUONS_SCALE__1down --systTitles "muon scale up,muon scale down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_muscale.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst EG_RESOLUTION_ALL__1up,EG_RESOLUTION_ALL__1down --systTitles "e res. up,e res. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_eres.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst EG_SCALE_ALL__1up,EG_SCALE_ALL__1down --systTitles "e scale up,e scale down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_escale.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst MET_SoftTrk_ResoPara --systTitles "MET res. para." -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_metrespara.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst MET_SoftTrk_ResoPerp --systTitles "MET res. perp." -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_metresperp.pdf --smoothen $S
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst MET_SoftTrk_ScaleUp,MET_SoftTrk_ScaleDown --systTitles "MET scale up,MET scale down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_metscale.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagbSF_0__1up,btagbSF_0__1down --systTitles "btag eff E0 up,btag eff E0 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagb0.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagbSF_0_pt1__1up,btagbSF_0_pt1__1down --systTitles "btag eff E0 pt1 up,btag eff E0 pt1 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagb0_pt1.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagbSF_0_pt2__1up,btagbSF_0_pt2__1down --systTitles "btag eff E0 pt2 up,btag eff E0 pt2 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagb0_pt2.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagbSF_0_pt3__1up,btagbSF_0_pt3__1down --systTitles "btag eff E0 pt3 up,btag eff E0 pt3 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagb0_pt3.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagbSF_1__1up,btagbSF_1__1down --systTitles "btag eff E1 up,btag eff E1 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagb1.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagbSF_2__1up,btagbSF_2__1down --systTitles "btag eff E2 up,btag eff E2 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagb2.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagbSF_3__1up,btagbSF_3__1down --systTitles "btag eff E3 up,btag eff E3 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagb3.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagbSF_4__1up,btagbSF_4__1down --systTitles "btag eff E4 up,btag eff E4 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagb4.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagcSF_0__1up,btagcSF_0__1down --systTitles "btag c mistag E0 up,btag c mistag E0 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagc0.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagcSF_0_pt1__1up,btagcSF_0_pt1__1down --systTitles "btag c mistag E0 pt1 up,btag c mistag E0 pt1 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagc0_pt1.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagcSF_0_pt2__1up,btagcSF_0_pt2__1down --systTitles "btag c mistag E0 pt2 up,btag c mistag E0 pt2 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagc0_pt2.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagcSF_0_pt3__1up,btagcSF_0_pt3__1down --systTitles "btag c mistag E0 pt3 up,btag c mistag E0 pt3 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagc0_pt3.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagcSF_1__1up,btagcSF_1__1down --systTitles "btag c mistag E1 up,btag c mistag E1 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagc1.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagcSF_2__1up,btagcSF_2__1down --systTitles "btag c mistag E2 up,btag c mistag E2 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagc2.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btagcSF_3__1up,btagcSF_3__1down --systTitles "btag c mistag E3 up,btag c mistag E3 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagc3.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_0__1up,btaglSF_0__1down --systTitles "btag l mistag E0 up,btag l mistag E0 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl0.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_0_pt1__1up,btaglSF_0_pt1__1down --systTitles "btag l mistag E0 pt1 up,btag l mistag E0 pt1 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl0_pt1.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_0_pt2__1up,btaglSF_0_pt2__1down --systTitles "btag l mistag E0 pt2 up,btag l mistag E0 pt2 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl0_pt2.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_0_pt3__1up,btaglSF_0_pt3__1down --systTitles "btag l mistag E0 pt3 up,btag l mistag E0 pt3 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl0_pt3.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_1__1up,btaglSF_1__1down --systTitles "btag l mistag E1 up,btag l mistag E1 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl1.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_2__1up,btaglSF_2__1down --systTitles "btag l mistag E2 up,btag l mistag E2 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl2.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_3__1up,btaglSF_3__1down --systTitles "btag l mistag E3 up,btag l mistag E3 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl3.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_4__1up,btaglSF_4__1down --systTitles "btag l mistag E4 up,btag l mistag E4 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl4.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_5__1up,btaglSF_5__1down --systTitles "btag l mistag E5 up,btag l mistag E5 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl5.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_6__1up,btaglSF_6__1down --systTitles "btag l mistag E6 up,btag l mistag E6 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl6.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_7__1up,btaglSF_7__1down --systTitles "btag l mistag E7 up,btag l mistag E7 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl7.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_8__1up,btaglSF_8__1down --systTitles "btag l mistag E8 up,btag l mistag E8 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl8.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_9__1up,btaglSF_9__1down --systTitles "btag l mistag E9 up,btag l mistag E9 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl9.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_10__1up,btaglSF_10__1down --systTitles "btag l mistag E10 up,btag l mistag E10 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl10.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btaglSF_11__1up,btaglSF_11__1down --systTitles "btag l mistag E11 up,btag l mistag E11 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btagl11.pdf --smoothen 1

    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btageSF_0__1up,btageSF_0__1down --systTitles "btag extrap up,btag extrap down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btage0.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst btageSF_1__1up,btageSF_1__1down --systTitles "btag extrap (c) up,btag extrap (c) down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_btage1.pdf --smoothen 1

    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst eTrigSF__1up,eTrigSF__1down --systTitles "e trig. up,e trig. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_etrig.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst eRecoSF__1up,eRecoSF__1down --systTitles "e rec. up,e rec. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_erec.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst eIDSF__1up,eIDSF__1down --systTitles "e ID up,e ID down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_eid.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst eIsolSF__1up,eIsolSF__1down --systTitles "e isol. up,e isol. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_eisol.pdf --smoothen 1

    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst muTrigStatSF__1up,muTrigStatSF__1down --systTitles "mu trig. stat. up,mu trig. stat. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_mutrigstat.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst muTrigSystSF__1up,muTrigSystSF__1down --systTitles "mu trig. syst. up,mu trig. syst. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_mutrigsyst.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst muIDStatSF__1up,muIDStatSF__1down --systTitles "mu ID stat. up,mu ID stat. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_muidstat.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst muIDSystSF__1up,muIDSystSF__1down --systTitles "mu ID syst. up,mu ID syst. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_muidsyst.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst muIsolStatSF__1up,muIsolStatSF__1down --systTitles "mu isol. stat. up,mu isol. stat. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_muisolstat.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst muIsolSystSF__1up,muIsolSystSF__1down --systTitles "mu isol syst. up,mu isol syst. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_muisolsyst.pdf --smoothen 1

    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst boostedWSF__1up,boostedWSF__1down --systTitles "boosted W C/A SF up,boosted W C/A SF down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_wmodel.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst ttEWK__1up,ttEWK__1down --systTitles "tt EWK corr up,tt EWK corr down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_ttewk.pdf --smoothen 1

    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst LARGERJET_JET_Top_CrossCalib__1up,LARGERJET_JET_Top_CrossCalib__1down --systTitles "akt10 cross calib. up,akt10 cross calib. down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_akt10xcalib.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst LARGERJET_JET_Top_Run1__1up,LARGERJET_JET_Top_Run1__1down --systTitles "akt10 run 1 up,akt10 run 1 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_akt10run1.pdf --smoothen 1

    ../plotting/plotCompareNominal --logY 1 --syst CORRELATED_LARGERJET_CORR1__1up,CORRELATED_LARGERJET_CORR1__1down --systTitles "akt10 run 1 and akt4 JES up,akt10 run 1 and akt4 JES down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_corr1.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --syst CORRELATED_LARGERJET_CORR2A__1up,CORRELATED_LARGERJET_CORR2A__1down --systTitles "akt10 run 1 JES and akt4 JES up,akt10 run 1 JES and akt4 JES down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_corr2a.pdf --smoothen 1
    ../plotting/plotCompareNominal --logY 1 --syst CORRELATED_LARGERJET_CORR2B__1up,CORRELATED_LARGERJET_CORR2B__1down --systTitles "akt10 run 1 JMS+Tau32 up,akt10 run 1 JMS+Tau32 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_${lep}_corr2b.pdf --smoothen 1

    ../plotting/plotCompareNominal --mcOnly 1 --syst btageSF_0__1up,btageSF_0__1down --systTitles "btag extrap up,btag extrap down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_signal_${lep}_btage0.pdf -C ../plotting/config_sig.txt --smoothen 0

    ../plotting/plotCompareNominal --mcOnly 1 --syst btagbSF_0__1up,btagbSF_0__1down --systTitles "btag eff. E0 up,btag eff. E0 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_signal_${lep}_btagb0.pdf -C ../plotting/config_sig.txt --smoothen 0
    ../plotting/plotCompareNominal --mcOnly 1 --syst btagcSF_0__1up,btagcSF_0__1down --systTitles "btag mistag c E0 up,btag mistag c E0 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_signal_${lep}_btagc0.pdf -C ../plotting/config_sig.txt --smoothen 0
    ../plotting/plotCompareNominal --mcOnly 1 --syst btaglSF_0__1up,btaglSF_0__1down --systTitles "btag mistag l E0 up,btag mistag l E0 down" -l $LUMI -c $lep -p $ch -h mtt -o syst_${ch}_mtt_signal_${lep}_btagl0.pdf -C ../plotting/config_sig.txt --smoothen 0

    for eig in $(seq 1 30) ; do
        ../plotting/plotCompareNominal --logY 1 --mcOnly 1 --syst pdf_PDF4LHC15_nlo_30_${eig} --systTitles "PDF4LHC15_nlo_30 E${eig}" -l $LUMI -c $lep -p $ch -h mtt -o systpdf_${ch}_mtt_${lep}_pdf${eig}.pdf --smoothen 1
    done

  done
done

for ch in boosted ; do
  for lep in e mu ; do
    ../plotting/plot --logY 1 --mcOnly 0 -c $lep -p $ch -h mtt -l $LUMI  --smoothen $S --xMax 3500 --stamp 1
  done
done

for ch in resolved ; do
  for lep in e mu ; do
    ../plotting/plot --logY 1 --mcOnly 1 -c $lep -p $ch -h mtt -l $LUMI  --smoothen $S --stamp 1
  done
done


for hist in largeJetEta largeJetPhi largeJetSd12 largeJet_tau32 largeJet_tau32_wta largeJet_tau21 largeJet_tau21_wta chi2 mwhad_res mthad_res mtlep_res jet0_m jet1_pt jet1_m jet2_pt jet2_m closejl_pt  closejl_minDeltaR MET_phi nBtagJets leadbJetPt lepEta lepPhi mwt  mtt nTrkBtagJets tjet1_pt tjet2_pt tjet3_pt ; do
for ch in resolved boosted ; do
  for lep in e mu ; do
    ../plotting/plot -c $lep -p $ch -h $hist -l $LUMI --smoothen 0 --stamp 1
  done
done
done

for hist in jet0_pt tjet0_pt ; do
for ch in resolved boosted ; do
  for lep in e mu ; do
    ../plotting/plot -c $lep -p $ch -h $hist -l $LUMI --smoothen 0 --normBinWidth 1 --stamp 1
  done
done
done

ch=boosted
hist=largeJetPt
for lep in e mu ; do
  ../plotting/plot -c $lep -p $ch -h $hist -l $LUMI --smoothen 0 --normBinWidth 1 --logY 1 --yMax 1e4 -o ${ch}_${hist}_${lep}_log.pdf --stamp 1
done


