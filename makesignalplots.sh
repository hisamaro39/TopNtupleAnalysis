
LUMI=3.20905

for sig in Zprime400 Zprime500 Zprime750 Zprime1000 Zprime1250 Zprime1500 Zprime1750 Zprime2000 Zprime2250 Zprime2500 Zprime2750 Zprime3000 Zprime4000 Zprime5000 ; do
  echo "sample    ${sig}     MC15_13TeV_25ns_FS_EXOT4_${sig}                    "${sig}"         "${sig}"         0" >config_sig_tmp.txt
  for ch in resolved boosted ; do
    for lep in e mu ; do
      ../plotting/plotCompareNominal --mcOnly 1 -l $LUMI -c $lep -p $ch -h mtt -o ${sig}_${ch}_mtt_${lep}.pdf -C config_sig_tmp.txt --smoothen 0
    done
  done
done

