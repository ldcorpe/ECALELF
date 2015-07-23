#!/bin/bash 
### prompt AOD reco'd with CMSSW_7_4_6_patch6
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251022-251495`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251496-251496`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251497-251558`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251559-251559`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251560-251602`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check 
#./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251560-251781`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251563-251781`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251882-251883`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check



#./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep SingleEle | grep WSkim | grep 2015B | grep 251022-251883`  -s WSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check
#./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep ALCARAW | grep SingleEle | grep ZSkim | grep 2015B | grep 251022-251883`  -s ZSkim --type EcalUncal  --tag 74X_dataRun2_Prompt_v0  --check
#./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep -v ALCARAW | grep DoubleEG | grep ZSkim | grep 2015B | grep 251022-251883`  -s ZSkim --type EcalCal  --tag 74X_dataRun2_Prompt_v0  --check
#./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep -v ALCARAW | grep SingleEle | grep WSkim | grep 2015B | grep 251022-251883`  -s WSkim --type EcalCal  --tag 74X_dataRun2_Prompt_v0   --check
#./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep -v  ALCARAW | grep SingleEle | grep ZSkim | grep 2015B | grep 251022-251883`  -s ZSkim --type EcalCal  --tag 74X_dataRun2_Prompt_v0  --check
