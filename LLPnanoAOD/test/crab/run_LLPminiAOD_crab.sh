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
filePerJob=1
VERSION=1

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
# --test \
# --dryrun \
# --input-DBS 'phys03' \
# --set-input-dataset \
#--send-external \
