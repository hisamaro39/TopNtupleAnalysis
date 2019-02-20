
S=1
#LUMI=3.31668
LUMI=3.20905

for ch in resolved boosted ; do
  for lep in e mu ; do
    ../plotting/plot -c $lep -p $ch -h largeJetM -l $LUMI --smoothen 0 --yTitle "Events / 10 GeV" --xMin 50 --yMax 850 --xTitle "Large-R jet mass [GeV]" --stamp 0
    ../plotting/plot -c $lep -p $ch -h mtlep_boo -l $LUMI --smoothen 0 --normBinWidth 20 --xMin 0 --yTitle "Events / 20 GeV" --xTitle "Mass of the leptonic top candidate [GeV]" --stamp 0
    ../plotting/plot -c $lep -p $ch -h MET -l $LUMI --smoothen 0 --normBinWidth 10 --xMin 20 --yTitle "Events / 10 GeV" --yMax 500 --xTitle "E_{T}^{miss} [GeV]" --stamp 0
    ../plotting/plot -c $lep -p $ch -h lepPt -l $LUMI --smoothen 0 --rebin 4 --xTitle "Lepton p_{T} [GeV]" --yTitle "Events / 20 GeV" --stamp 0
    ../plotting/plot -c $lep -p $ch -h closeJetPt -l $LUMI --smoothen 0 --normBinWidth 1 --yTitle "Events / GeV" --xTitle "Selected jet p_{T} [GeV]" --stamp 0
    ../plotting/plot -c $lep -p $ch -h largeJetPt -l $LUMI --smoothen 0 --normBinWidth 20 --yTitle "Events / 20 GeV" --xTitle "Large-R jet p_{T} [GeV]" --stamp 0
done
done
