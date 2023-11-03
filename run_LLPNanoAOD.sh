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

home_path=$(jq -r ".home_path" "$userparams_input")
CMSSW_base_path=$(jq -r ".CMSSW_base_path" "$userparams_input")
LLP_working_path=$(jq -r ".LLP_working_path" "$userparams_input")
output_base_path=$(jq -r ".output_base_path" "$userparams_input")
output_base_path="$output_base_path/$dataset_type/$dataset_name/LLPNanoAOD"
mkdir -p "$output_base_path"
dasgoclient_path=$(jq -r ".dasgoclient_path" "$userparams_input")
proxy_path=$(jq -r ".proxy_path" "$userparams_input")

source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
export X509_USER_PROXY=$proxy_path
voms-proxy-info

export HOME=$home_path

echo
echo -- Running $n_events number of events for dataset $dataset_name in $dataset_type
echo -- Using input parameter given in $userparams_input
 

### DATASET ###

dataset=$(jq -r ".datasets[\"$dataset_type\"][\"$dataset_name\"][\"dataset\"]" "$userparams_input")
input_format=$(jq -r ".datasets[\"$dataset_type\"][\"$dataset_name\"][\"input_format\"]" "$userparams_input")
if [ -z "$dataset" ] || [ -z "$input_format" ]; then
  echo "Error: Dataset information not found in the configuration file for dataset_type: $dataset_type, dataset_name: $dataset_name."
  exit 1
fi

echo -- Using input format $input_format
echo 

### INPUT FILENAME ###

if [ "$input_format" = "das" ]; then
  echo $dasgoclient_path -query='dataset=$dataset'
  $dasgoclient_path -query='dataset=$dataset"'
  echo filename="$($dasgoclient_path --query="file dataset=$dataset" -limit=0 | sed -n "$((process_id+2))p")"
  filename="$($dasgoclient_path --query="file dataset=$dataset" -limit=0 | sed -n "$((process_id+2))p")"
  if [ -z "$filename" ]; then
    echo "Error: filename from dataset $dataset_name in $dataset_type could not be found - exiting."
    exit 1
  fi
  input_path="root://cms-xrd-global.cern.ch//$filename"
  output_path="$output_base_path/$(basename "$filename")"

elif [ "$input_format" = "local" ]; then
  filename=$(ls "$dataset" | sed -n "$((process_id+1))p")
  if [ -z "$filename" ]; then
    echo "Error: filename from dataset $dataset_name in $dataset_type could not be found - exiting."
    exit 1
  fi
  input_path="file:$dataset/$filename"
  output_path="$output_base_path/$(basename "$input_path")"

else
  echo "Error: input_format $input_format not found - exiting."
  exit 1
fi

### CMSSW ENVIRONMENT ###
cd "${CMSSW_base_path}/CMSSW_10_6_29/src"
eval `scramv1 runtime -sh`

### RUN LLPNanoAOD ###
cmsRun "${LLP_working_path}/run_LLPNanoAOD.py" "$input_path" "$output_path" "$n_events"

echo "LLPNanoAOD file saved in: $output_path"
echo "Done"
