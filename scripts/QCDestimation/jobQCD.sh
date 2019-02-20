#!/bin/bash
#setup 2.3.40
args=("$@")
echo 'name of the output folder: ' ${args[0]}
echo 'name of the input_QCD folder:  ' ${args[1]}

#to execute: 
# mv jobQCD.sh ../../.
#cd ../../
#sh jobQCD.sh outputDir

#create the lists
#ls /AtlasDisk/users/romano/25ns_2.3.30/EXOT4_data_${args[1]}2/25ns/DT*root > input_QCD_DT.txt

make

#for the data/MC comparison in the fake CR
./read --files fakeListSamples/input_fake_ttbar.txt --analysis AnaTtresSL --data 0 --btags -1 --loose 0 --output resolved_e_ttbar.root,resolved_mu_ttbar.root,boosted_e_ttbar.root,boosted_mu_ttbar.root
./read --files fakeListSamples/input_fake_Wev.txt   --analysis AnaTtresSL --data 0 --btags -1 --loose 0 --output resolved_e_Wev.root,resolved_mu_Wev.root,boosted_e_Wev.root,boosted_mu_Wev.root        
./read --files fakeListSamples/input_fake_Wmv.txt   --analysis AnaTtresSL --data 0 --btags -1 --loose 0 --output resolved_e_Wmv.root,resolved_mu_Wmv.root,boosted_e_Wmv.root,boosted_mu_Wmv.root        
./read --files fakeListSamples/input_fake_Wtv.txt   --analysis AnaTtresSL --data 0 --btags -1 --loose 0 --output resolved_e_Wtv.root,resolved_mu_Wtv.root,boosted_e_Wtv.root,boosted_mu_Wtv.root        
./read --files fakeListSamples/input_fake_Zee.txt   --analysis AnaTtresSL --data 0 --btags -1 --loose 0 --output resolved_e_Zee.root,resolved_mu_Zee.root,boosted_e_Zee.root,boosted_mu_Zee.root        
./read --files fakeListSamples/input_fake_Zmm.txt   --analysis AnaTtresSL --data 0 --btags -1 --loose 0 --output resolved_e_Zmm.root,resolved_mu_Zmm.root,boosted_e_Zmm.root,boosted_mu_Zmm.root        
./read --files fakeListSamples/input_fake_Ztt.txt   --analysis AnaTtresSL --data 0 --btags -1 --loose 0 --output resolved_e_Ztt.root,resolved_mu_Ztt.root,boosted_e_Ztt.root,boosted_mu_Ztt.root        
./read --files fakeListSamples/input_fake_st.txt    --analysis AnaTtresSL --data 0 --btags -1 --loose 0 --output resolved_e_st.root,resolved_mu_st.root,boosted_e_st.root,boosted_mu_st.root	      
./read --files fakeListSamples/input_fake_DT.txt    --analysis AnaTtresSL --data 1 --btags -1 --loose 0 --output resolved_e_DT.root,resolved_mu_DT.root,boosted_e_DT.root,boosted_mu_DT.root	      

#generate the multijet bkg
./read --files fakeListSamples/input_fake_DT.txt    --analysis AnaTtresSL --data 1 --btags -1 --loose 1 --runMM 1 --output resolved_e_QCD.root,resolved_mu_QCD.root,boosted_e_QCD.root,boosted_mu_QCD.root 

rm -f *mini.root

#resolved_QCD channel
mkdir -p ${args[0]}_resolved_QCD
mv resolved_*root ${args[0]}_resolved_QCD/.

#boosted_QCD
mkdir -p ${args[0]}_boosted_QCD
mv boosted_*root ${args[0]}_boosted_QCD/.

