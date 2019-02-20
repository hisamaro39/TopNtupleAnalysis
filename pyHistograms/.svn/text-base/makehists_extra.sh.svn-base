
LUMI=14.76518

PLOTTING=/afs/desy.de/user/d/danilo/xxl/af-atlas/Top2418/TopNtupleAnalysis/plotting/plot

rm -f hist_*.root
for histogram in mtt mttPos mttNeg ; do
    for ch in be bmu re rmu ; do
        $PLOTTING  -c $ch -h $histogram -l $LUMI --saveTH1 ${histogram}${ch}  --smoothen 1 -C config_limit.txt
    done
done

#rm spectrum_all.tar.gz
#tar cvfz spectrum_all.tar.gz hist*.root

