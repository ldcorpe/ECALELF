#!/bin/bash 
### prompt AOD reco'd with CMSSW_7_4_6_patch6
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251022-251883`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0 
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep SingleEle | grep WSkim | grep 2015B | grep 251022-251883`  -s WSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep SingleEle | grep ZSkim | grep 2015B | grep 251022-251883`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0 
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep -v ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251022-251883`  -s ZSkim --type EcalCal  --tag 74X_dataRun2_Prompt_v0 
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep -v ALCARAW | grep SingleEle | grep WSkim | grep 2015B | grep 251022-251883`  -s WSkim --type EcalCal  --tag 74X_dataRun2_Prompt_v0
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep -v  ALCARAW | grep SingleEle | grep ZSkim | grep 2015B | grep 251022-251883`  -s ZSkim --type EcalCal  --tag 74X_dataRun2_Prompt_v0
