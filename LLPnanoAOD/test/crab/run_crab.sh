#!/bin/bash



# File with sample list
INPUT=$1

if [ -z "$INPUT" ]; then
    echo "Missing input argument."
    exit 1
fi
filename=$(basename "$INPUT" .conf)

VERSION=1

# Paths
crabWorkspace=/afs/desy.de/user/l/lrygaard/TTALP/LLPNanoAOD/crab
LLPnanoAODWorkspace=/afs/desy.de/user/l/lrygaard/TTALP/LLPNanoAOD

runFile=run_LLPminiAOD.py

source /cvmfs/cms.cern.ch/common/crab-setup.sh

python $crabWorkspace/crab.py \
-p $LLPnanoAODWorkspace/$runFile \
--site T2_DE_DESY \
-o /store/user/$USER/TTALPs \
-t LLPnanoAOD-v$VERSION \
-i $INPUT \
-s FileBased \
-n 10 \
--includeDSAMuon \
--includeBS \
--includeGenPart \
--max-memory 3000 \
--max-runtime-min 2750 \
--work-area $crabWorkspace/crab_projects/crab_${filename}_$VERSION \
--dryrun
# --set-input-dataset \
#--send-external \
