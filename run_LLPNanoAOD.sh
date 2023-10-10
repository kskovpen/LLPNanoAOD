#!/bin/bash

#source /osg/osg3.2/osg-wn-client/setup.sh
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git


### INPUTS ###

process_id=$(( $1+0 ))
dataset_type=$2
dataset_name=$3
n_events=$4

### PATHS ###

userparams_input="UserParameters/ttalps_userparams.json"

python_path=$(jq -r ".python_path" "$userparams_input")
CMSSW_base_path=$(jq -r ".CMSSW_base_path" "$userparams_input")
LLP_working_path=$(jq -r ".LLP_working_path" "$userparams_input")
output_base_path=$(jq -r ".output_base_path" "$userparams_input")
output_base_path="$output_base_path/$dataset_type/$dataset_name/LLPNanoAOD"
mkdir -p "$output_base_path"

echo
echo -- Running $n_events number of events for dataset $dataset_name in $dataset_type
echo -- Using input parameter given in $userparams_input
echo 

### DATASET ###

dataset=$(jq -r ".datasets[\"$dataset_type\"][\"$dataset_name\"][\"dataset\"]" "$userparams_input")
input_format=$(jq -r ".datasets[\"$dataset_type\"][\"$dataset_name\"][\"input_format\"]" "$userparams_input")
if [ -z "$dataset" ] || [ -z "$input_format" ]; then
  echo "Error: Dataset information not found in the configuration file for dataset_type: $dataset_type, dataset_name: $dataset_name."
  exit 1
fi

### INPUT FILENAME ###

if [ "$input_format" = "das" ]; then
  filename="$(dasgoclient --query="file dataset=$dataset" -limit=0 | sed -n "$((process_id+2))p")"
  input_path="root://xrootd-cms.infn.it//$filename"
  output_path="$output_base_path/$(basename "$filename")"

elif [ "$input_format" = "local" ]; then
  input_path="file:$dataset/$(ls "$dataset" | sed -n "$((process_id+1))p")"
  output_path="$output_base_path/$(basename "$input_path")"

else
  echo "Error: input_format $input_format not found - exiting."
  exit 1
fi

if [ -z "$filename" ]; then
  echo "Error: filename from dataset $dataset_name in $dataset_type could not be found - exiting."
  exit 1
fi

### CMSSW ENVIRONMENT ###
cd "${CMSSW_base_path}/CMSSW_10_6_29/src"
eval `scramv1 runtime -sh`
echo $CMSSW_BASE
ls PhysicsTools/NanoAOD/plugins/

### RUN LLPNanoAOD ###
cmsRun "${LLP_working_path}/run_LLPNanoAOD.py" "$input_path" "$output_path" "$n_events"

echo "LLPNanoAOD file saved in: $output_path"
echo "Done"
