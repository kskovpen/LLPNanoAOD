#!/bin/bash


# File with sample list
INPUT=$1

if [ -z "$INPUT" ]; then
    echo "Missing input argument."
    exit 1
fi
filename=$(basename "$INPUT" .conf)

# Paths
crabWorkspace=$CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test/crab
configWorkspace=$CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test

runFile=LLPminiAOD_cfg.py

source /cvmfs/cms.cern.ch/common/crab-setup.sh

nCores=4
maxMemory=$((1000 * $nCores))
maxRuntime=2750
# Setup for FileBased splitting, 0 to automatically set lowest number of files per job for max total 10000 jobs
filePerJob=0
# LLPminiAOD version
VERSION=1

whitelist="['T2_US_*', 'T2_US_*', 'T2_CH_*', 'T2_IT_*', 'T1_IT_*']"

# Year options for MC: 2016, 2017, 2018, 2022PreEE, 2022PostEE, 2023PreBPix, 2023PostBPix
# Year options for data: 2016HIPM, 2016 (no HIPM), 2017, 2018, 2022ReReco, 2022Prompt, 2023
year=2017

python $crabWorkspace/crab.py \
-p $configWorkspace/$runFile \
--site T2_DE_DESY \
-o /store/user/$USER/ttalps \
-t LLPminiAODv$VERSION \
-i $INPUT \
-s FileBased \
-n $filePerJob \
--num-cores $nCores \
--max-memory $maxMemory \
--max-runtime-min $maxRuntime \
--work-area $crabWorkspace/crab_projects/crab_${filename}_v$VERSION \
--publication \
--year $year \
--dryrun \
# --ignore_locality \
# --whitelist "$whitelist" \
# --runOnData \
# --test \
# --input-DBS 'phys03' \
# --set-input-dataset \
#--send-external \
