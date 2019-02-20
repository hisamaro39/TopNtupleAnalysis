
LUMI=14.76518

PLOTTING=/afs/desy.de/user/d/danilo/xxl/af-atlas/Top2418/TopNtupleAnalysis/plotting/plot
PLOTCOMP=/afs/desy.de/user/d/danilo/xxl/af-atlas/Top2418/TopNtupleAnalysis/plotting/plotCompareNominal

for ch in be bmu re rmu be2015 bmu2015 re2015 rmu2015 be2016 bmu2016 re2016 rmu2016 ; do
    $PLOTTING -c $ch -h yields -l $LUMI -C config.txt >yields_${ch}.txt
    $PLOTTING -c $ch -h yieldsPos -l $LUMI -C config.txt >yieldsPos_${ch}.txt
    $PLOTTING -c $ch -h yieldsNeg -l $LUMI -C config.txt >yieldsNeg_${ch}.txt
done

for ch in be bmu re rmu be2015 bmu2015 re2015 rmu2015 be2016 bmu2016 re2016 rmu2016 ; do
    $PLOTTING -c $ch -h closestJetDr  -l $LUMI --xTitle "min #Delta R (jet, lepton)" --yTitle "Events / 0.1" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h lepPt --normBinWidth 20 -l $LUMI --rebin 4 --xTitle "Lepton p_{T} [GeV]" --yTitle "Events / 20 GeV" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h lepEta  -l $LUMI  --xTitle "Lepton #eta" --yTitle "Events / 0.25" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h nJets -l $LUMI --xTitle "Number of jets" --yTitle "Events" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h nTrkBtagJets -l $LUMI --xTitle "Number of b-tagged track jets" --yTitle "Events" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h MET --normBinWidth 20 -l $LUMI --yTitle "Events / 20 GeV" --xTitle "E_{T}^{miss} [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h MET_phi -l $LUMI --yTitle "Events / 0.2" --xTitle "E_{T}^{miss} #phi [rd]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mwt -l $LUMI --yTitle "Events / 10 GeV" --xTitle "m_{WT} [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h closeJetPt -l $LUMI --normBinWidth 10 --yTitle "Events / 10 GeV" --xTitle "Selected jet p_{T} [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h largeJetPt -l $LUMI --normBinWidth 20 --yTitle "Events / 20 GeV" --xTitle "Large-R jet p_{T} [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h largeJetM -l $LUMI  --yTitle "Events / 10 GeV" --xTitle "Large-R jet mass [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h largeJetPtMtt -l $LUMI  --yTitle "Events / 0.02" --xTitle "Large-R jet p_{T} / m_{t#bar{t}} " --stamp 0 -C config.txt
    $PLOTTING -c $ch -h largeJetEta -l $LUMI  --yTitle "Events / 0.1" --xTitle "Large-R jet #eta" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h largeJetPhi -l $LUMI  --yTitle "Events / 0.2" --xTitle "Large-R jet #phi [rd]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mtlep_boo -l $LUMI --normBinWidth 20  --yTitle "Events / 20 GeV" --xTitle "Mass of the leptonic top [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mtlep_res -l $LUMI --normBinWidth 20  --yTitle "Events / 20 GeV" --xTitle "Mass of the leptonic top [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mthad_res -l $LUMI --normBinWidth 20  --yTitle "Events / 20 GeV" --xTitle "Mass of the hadronic top [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mwhad_res -l $LUMI --normBinWidth 20  --yTitle "Events / 20 GeV" --xTitle "Mass of the hadronic W [GeV]" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h chi2 -l $LUMI  --yTitle "Events / 0.2" --xTitle "log(#chi^{2})" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h largeJet_tau32_wta -l $LUMI  --yTitle "Events / 0.05" --xTitle "Large-R jet #tau_{32} wta" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h largeJet_tau21_wta -l $LUMI  --yTitle "Events / 0.05" --xTitle "Large-R jet #tau_{21} wta" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h trueMtt --mcOnly 1 --normBinWidth 100 --logY 1 -l $LUMI --xTitle "true m_{t#bar{t}}" --yTitle "Events / 100 GeV" --stamp 0 -C config.txt
done

for ch in be bmu re rmu be2015 bmu2015 re2015 rmu2015 be2016 bmu2016 re2016 rmu2016 ; do

    $PLOTCOMP --logY 1 --mcOnly 1 --syst qcd__1up,qcd__1down --systTitles "QCD up,QCD down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_qcd.pdf
    $PLOTCOMP --mcOnly 0 --syst qcd__1up,qcd__1down --systTitles "QCD up,QCD down" -l $LUMI -c $ch -h MET -o syst_${ch}_MET_qcd.pdf
    $PLOTCOMP --mcOnly 0 --syst qcd__1up,qcd__1down --systTitles "QCD up,QCD down" -l $LUMI -c $ch -h mwt -o syst_${ch}_mwt_qcd.pdf
    $PLOTCOMP --mcOnly 0 --syst WJETS__1up,WJETS__1down --systTitles "Wjets 100% up,Wjets 100% down" -l $LUMI -c $ch -h mwt -o syst_${ch}_mwt_wjets100p.pdf

    # systematics in mtt
    $PLOTCOMP --logY 1 --mcOnly 1 -c $ch -h mtt -l $LUMI --other ttpowhegherwig,ttmcatnloherwig --titles "Powheg+Herwig,MC@NLO+Herwig" -o systmodel_${ch}_mtt_mcgen.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 -c $ch -h mtt -l $LUMI --other ttsyst,ttpowhegherwig --titles "Powheg+Pythia,Powheg+Herwig" -o systmodel_${ch}_mtt_pshower.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 -c $ch -h mtt -l $LUMI --other ttradlo,ttradhi --titles "ISR/FSR(low),ISR/FSR(high)" -o systmodel_${ch}_mtt_isrfsr.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_JER_SINGLE_NP__1up --systTitles "akt4 JER" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jer.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EffectiveNP_1__1up,JET_19NP_JET_EffectiveNP_1__1down --systTitles "akt4 JES 1 up,akt4 JES 1 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes1.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EffectiveNP_2__1up,JET_19NP_JET_EffectiveNP_2__1down --systTitles "akt4 JES 2 up,akt4 JES 2 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes2.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EffectiveNP_3__1up,JET_19NP_JET_EffectiveNP_3__1down --systTitles "akt4 JES 3 up,akt4 JES 3 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes3.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EffectiveNP_4__1up,JET_19NP_JET_EffectiveNP_4__1down --systTitles "akt4 JES 4 up,akt4 JES 4 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes4.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EffectiveNP_5__1up,JET_19NP_JET_EffectiveNP_5__1down --systTitles "akt4 JES 5 up,akt4 JES 5 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes5.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EffectiveNP_6restTerm__1up,JET_19NP_JET_EffectiveNP_6restTerm__1down --systTitles "akt4 JES 6 up,akt4 JES 6 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes6.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_Pileup_RhoTopology__1up,JET_19NP_JET_Pileup_RhoTopology__1down --systTitles "akt4 JES pile up rho topo. up,akt4 JES pile up rho topo. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jespurho.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_Pileup_OffsetNPV__1up,JET_19NP_JET_Pileup_OffsetNPV__1down --systTitles "akt4 JES pile up offset NPV up,akt4 JES pile up offset NPV down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jespuoffnpv.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_Pileup_OffsetMu__1up,JET_19NP_JET_Pileup_OffsetMu__1down --systTitles "akt4 JES pile up offset mu up,akt4 JES pile up offset mu down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jespuoffmu.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_Pileup_PtTerm__1up,JET_19NP_JET_Pileup_PtTerm__1down --systTitles "akt4 JES pile up pt term up,akt4 JES pile up pt term down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jespuoffmu.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_BJES_Response__1up,JET_19NP_JET_BJES_Response__1down --systTitles "akt4 bJES resp. up,akt4 bJES resp. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4bjes.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EtaIntercalibration_TotalStat__1up,JET_19NP_JET_EtaIntercalibration_TotalStat__1down --systTitles "akt4 eta intercalib. stat. up,akt4 eta intercalib. stat. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jesetastat.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EtaIntercalibration_NonClosure__1up,JET_19NP_JET_EtaIntercalibration_NonClosure__1down --systTitles "akt4 eta intercalib. non-closure up,akt4 eta intercalib. non-closure down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jesetanc.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_EtaIntercalibration_Modelling__1up,JET_19NP_JET_EtaIntercalibration_Modelling__1down --systTitles "akt4 eta intercalib. model. up,akt4 eta intercalib. model. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jesetamod.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_Flavor_Response__1up,JET_19NP_JET_Flavor_Response__1down --systTitles "akt4 JES flav. resp. up,akt4 JES flav. resp. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jesfr.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_Flavor_Composition__1up,JET_19NP_JET_Flavor_Composition__1down --systTitles "akt4 JES flav. comp. up,akt4 JES flav. comp. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jesfc.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_PunchThrough_MC15__1up,JET_19NP_JET_PunchThrough_MC15__1down --systTitles "akt4 JES punchthrough up,akt4 JES punchthrough down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jespt.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst JET_19NP_JET_SingleParticle_HighPt__1up,JET_19NP_JET_SingleParticle_HighPt__1down --systTitles "akt4 JES single part. up,akt4 JES single part. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jessp.pdf

    #$PLOTCOMP --logY 1 --mcOnly 1 --syst JET_NPScenario1_JET_GroupedNP_1__1up,JET_NPScenario1_JET_GroupedNP_1__1down --systTitles "akt4 JES 1 up,akt4 JES 1 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes1.pdf
    #$PLOTCOMP --logY 1 --mcOnly 1 --syst JET_NPScenario1_JET_GroupedNP_2__1up,JET_NPScenario1_JET_GroupedNP_2__1down --systTitles "akt4 JES 2 up,akt4 JES 2 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes2.pdf
    #$PLOTCOMP --logY 1 --mcOnly 1 --syst JET_NPScenario1_JET_GroupedNP_3__1up,JET_NPScenario1_JET_GroupedNP_3__1down --systTitles "akt4 JES 3 up,akt4 JES 3 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jes3.pdf
    #$PLOTCOMP --logY 1 --mcOnly 1 --syst JET_NPScenario1_JET_EtaIntercalibration_NonClosure__1up,JET_NPScenario1_JET_EtaIntercalibration_NonClosure__1down --systTitles "akt4 JES eta interc. up,akt4 JES eta interc. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt4jesetainter.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst MUONS_MS__1up,MUONS_MS__1down --systTitles "muon res. MS up,muon res. MS down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_muresms.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst MUONS_ID__1up,MUONS_ID__1down --systTitles "muon res. ID up,muon res. ID down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_muresid.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst MUONS_SCALE__1up,MUONS_SCALE__1down --systTitles "muon scale up,muon scale down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_muscale.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst EG_RESOLUTION_ALL__1up,EG_RESOLUTION_ALL__1down --systTitles "e res. up,e res. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_eres.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst EG_SCALE_ALL__1up,EG_SCALE_ALL__1down --systTitles "e scale up,e scale down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_escale.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst MET_SoftTrk_ResoPara --systTitles "MET res. para." -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_metrespara.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst MET_SoftTrk_ResoPerp --systTitles "MET res. perp." -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_metresperp.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst MET_SoftTrk_ScaleUp,MET_SoftTrk_ScaleDown --systTitles "MET scale up,MET scale down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_metscale.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagbSF_0__1up,btagbSF_0__1down --systTitles "btag eff E0 up,btag eff E0 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagb0.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagbSF_0_pt1__1up,btagbSF_0_pt1__1down --systTitles "btag eff E0 pt1 up,btag eff E0 pt1 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagb0_pt1.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagbSF_0_pt2__1up,btagbSF_0_pt2__1down --systTitles "btag eff E0 pt2 up,btag eff E0 pt2 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagb0_pt2.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagbSF_0_pt3__1up,btagbSF_0_pt3__1down --systTitles "btag eff E0 pt3 up,btag eff E0 pt3 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagb0_pt3.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagbSF_1__1up,btagbSF_1__1down --systTitles "btag eff E1 up,btag eff E1 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagb1.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagbSF_2__1up,btagbSF_2__1down --systTitles "btag eff E2 up,btag eff E2 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagb2.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagbSF_3__1up,btagbSF_3__1down --systTitles "btag eff E3 up,btag eff E3 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagb3.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagbSF_4__1up,btagbSF_4__1down --systTitles "btag eff E4 up,btag eff E4 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagb4.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagcSF_0__1up,btagcSF_0__1down --systTitles "btag c mistag E0 up,btag c mistag E0 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagc0.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagcSF_0_pt1__1up,btagcSF_0_pt1__1down --systTitles "btag c mistag E0 pt1 up,btag c mistag E0 pt1 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagc0_pt1.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagcSF_0_pt2__1up,btagcSF_0_pt2__1down --systTitles "btag c mistag E0 pt2 up,btag c mistag E0 pt2 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagc0_pt2.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagcSF_0_pt3__1up,btagcSF_0_pt3__1down --systTitles "btag c mistag E0 pt3 up,btag c mistag E0 pt3 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagc0_pt3.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagcSF_1__1up,btagcSF_1__1down --systTitles "btag c mistag E1 up,btag c mistag E1 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagc1.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagcSF_2__1up,btagcSF_2__1down --systTitles "btag c mistag E2 up,btag c mistag E2 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagc2.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btagcSF_3__1up,btagcSF_3__1down --systTitles "btag c mistag E3 up,btag c mistag E3 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagc3.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_0__1up,btaglSF_0__1down --systTitles "btag l mistag E0 up,btag l mistag E0 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl0.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_0_pt1__1up,btaglSF_0_pt1__1down --systTitles "btag l mistag E0 pt1 up,btag l mistag E0 pt1 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl0_pt1.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_0_pt2__1up,btaglSF_0_pt2__1down --systTitles "btag l mistag E0 pt2 up,btag l mistag E0 pt2 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl0_pt2.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_0_pt3__1up,btaglSF_0_pt3__1down --systTitles "btag l mistag E0 pt3 up,btag l mistag E0 pt3 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl0_pt3.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_1__1up,btaglSF_1__1down --systTitles "btag l mistag E1 up,btag l mistag E1 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl1.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_2__1up,btaglSF_2__1down --systTitles "btag l mistag E2 up,btag l mistag E2 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl2.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_3__1up,btaglSF_3__1down --systTitles "btag l mistag E3 up,btag l mistag E3 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl3.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_4__1up,btaglSF_4__1down --systTitles "btag l mistag E4 up,btag l mistag E4 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl4.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_5__1up,btaglSF_5__1down --systTitles "btag l mistag E5 up,btag l mistag E5 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl5.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_6__1up,btaglSF_6__1down --systTitles "btag l mistag E6 up,btag l mistag E6 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl6.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_7__1up,btaglSF_7__1down --systTitles "btag l mistag E7 up,btag l mistag E7 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl7.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_8__1up,btaglSF_8__1down --systTitles "btag l mistag E8 up,btag l mistag E8 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl8.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_9__1up,btaglSF_9__1down --systTitles "btag l mistag E9 up,btag l mistag E9 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl9.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btaglSF_10__1up,btaglSF_10__1down --systTitles "btag l mistag E10 up,btag l mistag E10 down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btagl10.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst btageSF_0__1up,btageSF_0__1down --systTitles "btag extrap up,btag extrap down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btage0.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst btageSF_1__1up,btageSF_1__1down --systTitles "btag extrap (c) up,btag extrap (c) down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_btage1.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst eTrigSF__1up,eTrigSF__1down --systTitles "e trig. up,e trig. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_etrig.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst eRecoSF__1up,eRecoSF__1down --systTitles "e rec. up,e rec. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_erec.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst eIDSF__1up,eIDSF__1down --systTitles "e ID up,e ID down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_eid.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst eIsolSF__1up,eIsolSF__1down --systTitles "e isol. up,e isol. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_eisol.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst muTrigStatSF__1up,muTrigStatSF__1down --systTitles "mu trig. stat. up,mu trig. stat. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_mutrigstat.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst muTrigSystSF__1up,muTrigSystSF__1down --systTitles "mu trig. syst. up,mu trig. syst. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_mutrigsyst.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst muIDStatSF__1up,muIDStatSF__1down --systTitles "mu ID stat. up,mu ID stat. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_muidstat.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst muIDSystSF__1up,muIDSystSF__1down --systTitles "mu ID syst. up,mu ID syst. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_muidsyst.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst muIsolStatSF__1up,muIsolStatSF__1down --systTitles "mu isol. stat. up,mu isol. stat. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_muisolstat.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst muIsolSystSF__1up,muIsolSystSF__1down --systTitles "mu isol syst. up,mu isol syst. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_muisolsyst.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst wnorm__1up,wnorm__1down --systTitles "W C/A SF up,W C/A SF down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_wnorm.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst wbb__1up,wbb__1down --systTitles "W+bb HF SF up,W+bb HF SF down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_wbb.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst wcc__1up,wcc__1down --systTitles "W+cc HF SF up,W+cc HF SF down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_wcc.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst wc__1up,wc__1down --systTitles "W+c HF SF up,W+c HF SF down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_wc.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst wl__1up,wl__1down --systTitles "W+l HF SF up,W+l HF SF down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_wl.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst ttEWK__1up,ttEWK__1down --systTitles "tt EWK corr up,tt EWK corr down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_ttewk.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst jvtSF__1up,jvtSF__1down --systTitles "JVT SF up,JVT SF down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_jvtsf.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst LARGERJET_Strong_JET_Rtrk_Modelling_All__1up,LARGERJET_Strong_JET_Rtrk_Modelling_All__1down --systTitles "akt10 rtrk mod. up,akt10 rtrk mod. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt10rtrkmodel.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst LARGERJET_Strong_JET_Rtrk_Baseline_All__1up,LARGERJET_Strong_JET_Rtrk_Baseline_All__1down --systTitles "akt10 rtrk bl. up,akt10 rtrk bl. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt10rtrkbl.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst LARGERJET_Strong_JET_Rtrk_Tracking_All__1up,LARGERJET_Strong_JET_Rtrk_Tracking_All__1down --systTitles "akt10 rtrk trk. up,akt10 rtrk trk. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt10rtrktrk.pdf
    $PLOTCOMP --logY 1 --mcOnly 1 --syst LARGERJET_Strong_JET_Rtrk_TotalStat_All__1up,LARGERJET_Strong_JET_Rtrk_TotalStat_All__1down --systTitles "akt10 rtrk stat. up,akt10 rtrk stat. down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_akt10rtrkstat.pdf

    $PLOTCOMP --logY 1 --mcOnly 1 --syst pileupSF__1up,pileupSF__1down --systTitles "pile up up,pile up down" -l $LUMI -c $ch -h mtt -o syst_${ch}_mtt_pile.pdf

    for eig in $(seq 0 30) ; do
        $PLOTCOMP --logY 1 --mcOnly 1 --syst pdf_PDF4LHC15_nlo_30_${eig} --systTitles "PDF4LHC15_nlo_30 E${eig}" -l $LUMI -c $ch -h mtt -o systpdf_${ch}_mtt_pdf${eig}.pdf
    done
done

for ch in be bmu re rmu re2015 rmu2015 be2016 bmu2016 re2016 rmu2016 ; do
    $PLOTTING -c $ch -h mtt --mcOnly 1 --normBinWidth 100 --logY 1 -l $LUMI --xTitle "m_{t#bar{t}}" --yTitle "Events / 100 GeV" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mttPos --mcOnly 1 --normBinWidth 100 --logY 1 -l $LUMI --xTitle "m_{t#bar{t}}^{+}" --yTitle "Events / 100 GeV" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mttNeg --mcOnly 1 --normBinWidth 100 --logY 1 -l $LUMI --xTitle "m_{t#bar{t}}^{-}" --yTitle "Events / 100 GeV" --stamp 0 -C config.txt
done

for ch in bmu2015 be2015 ; do
    $PLOTTING -c $ch -h mtt --mcOnly 0 --normBinWidth 100 --logY 1 -l $LUMI --xTitle "m_{t#bar{t}}" --yTitle "Events / 100 GeV" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mttPos --mcOnly 0 --normBinWidth 100 --logY 1 -l $LUMI --xTitle "m_{t#bar{t}}^{+}" --yTitle "Events / 100 GeV" --stamp 0 -C config.txt
    $PLOTTING -c $ch -h mttNeg --mcOnly 0 --normBinWidth 100 --logY 1 -l $LUMI --xTitle "m_{t#bar{t}}^{-}" --yTitle "Events / 100 GeV" --stamp 0 -C config.txt
done

