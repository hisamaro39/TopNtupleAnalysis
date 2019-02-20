#!/bin/bash

for CHAN in be bmu re rmu
do

  root -l -b -q '../scripts/plot.C+("largeJetPt","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("largeJetM","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("largeJetSd12","./","_'${CHAN}'.root",false)'

  root -l -b -q '../scripts/plot.C+("lepPt","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("lepEta","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("lepPhi","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("leadJetPt","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("leadbJetPt","./","_'${CHAN}'.root",false)'

  root -l -b -q '../scripts/plot.C+("met","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("met_phi","./","_'${CHAN}'.root",false)'

  root -l -b -q '../scripts/plot.C+("closeJetPt","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("mtlep_boo","./","_'${CHAN}'.root",false)'

  root -l -b -q '../scripts/plot.C+("mtlep_res","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("mthad_res","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("mwhad_res","./","_'${CHAN}'.root",false)'
  root -l -b -q '../scripts/plot.C+("chi2","./","_'${CHAN}'.root",false)'

  root -l -b -q '../scripts/plot.C+("mtt","./","_'${CHAN}'.root",false)'

done

