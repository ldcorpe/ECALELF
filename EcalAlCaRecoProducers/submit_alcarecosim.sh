#!/bin/bash 
### prompt AOD reco'd with CMSSW_7_4_6_patch6
#./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep DYJetsToLL_M-50_13TeV-Asympt50ns-pythia8_21Jul15_v4`  -s ZSkim --type EcalCal  --tag 74X_dataRun2_Prompt_v0 --isMC 
#./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep DYJetsToLL_M-50_13TeV-Asympt50ns-pythia8_21Jul15_v5`  -s ZSkim --type EcalCal  --tag 74X_dataRun2_Prompt_v0 --isMC 
./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep DYJetsToLL_M-50_13TeV-Asympt50ns-pythia8_21Jul15_v5`  -s ZSkim --type ALCARECOSIM  --isMC --json_name v2 --check 
