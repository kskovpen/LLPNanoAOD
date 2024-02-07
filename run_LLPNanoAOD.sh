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
files_per_job=$4
n_events=$5
run_mini_and_nano=$6

### PATHS ###

userparams_input="UserParameters/ttalps_userparams_private.json"

home_path=$(jq -r ".home_path" "$userparams_input")
CMSSW_base_path=$(jq -r ".CMSSW_base_path" "$userparams_input")
LLP_working_path=$(jq -r ".LLP_working_path" "$userparams_input")
output_base_path=$(jq -r ".output_base_path" "$userparams_input")
tmp_miniAOD_path="$output_base_path/LLPtest/$dataset_type/$dataset_name/LLPminiAOD"
output_base_path="$output_base_path/LLPtest/$dataset_type/$dataset_name/LLPNanoAOD"
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

if [[ $run_mini_and_nano == 1 ]]; then
  echo -- Running both miniAOD and nanoAOD, miniAODs temporarly stored in: $tmp_miniAOD_path
fi

### DATASET ###

dataset=$(jq -r ".datasets[\"$dataset_type\"][\"$dataset_name\"][\"dataset\"]" "$userparams_input")
input_format=$(jq -r ".datasets[\"$dataset_type\"][\"$dataset_name\"][\"input_format\"]" "$userparams_input")
if [ -z "$dataset" ] || [ -z "$input_format" ]; then
  echo "Error: Dataset information not found in the configuration file for dataset_type: $dataset_type, dataset_name: $dataset_name."
  exit 1
fi

is_data=false
if [[ $dataset_type == "data"* ]]; then
  is_data=true
fi

echo -- Using input format $input_format
echo 

### INPUT FILENAME ###

filename_list=""

if [ "$input_format" = "das" ]; then

  total_files="$($dasgoclient_path --query="file dataset=$dataset" -limit=0 | wc -l)"

  for ((i=0; i<$files_per_job; i++)); do

    file_number=$((process_id*files_per_job+i+1))

    # Check that the file number doesn't exceed total number of files
    if [ $file_number -le $total_files ]; then 

      # Reading filename 
      filename="$($dasgoclient_path --query="file dataset=$dataset" -limit=0 | sed -n "${file_number}p")"
      # Check if filename was found for the first file
      if [ $i -eq 0 ]; then
        if [ -z "$filename" ]; then
          echo "Error: filename from dataset $dataset_name in $dataset_type could not be found - exiting."
          exit 1
        fi
      fi
      # Adding redirector for DAS root files
      # input_path="root://cms-xrd-global.cern.ch/$filename"
      input_path="root://xrootd-cms.infn.it/$filename"
      # input_path="root://cmsxrootd.fnal.gov/$filename"
      # input_path="root://cmseos.fnal.gov/$filename"
      # input_path="root://eoscms.cern.ch/$filename"
      # Adding path to list of filenames, separated by a comma
      filename_list+="$input_path,"

    fi
  done

elif [ "$input_format" = "local" ]; then

  total_files=$(ls "$dataset" | wc -l)

  for ((i=0; i<$files_per_job; i++)); do

    file_number=$((process_id*files_per_job+i+1))

    # Check that the file number doesn't exceed total number of files
    if [ $file_number -le $total_files ]; then

      # Reading filename 
      filename=$(ls "$dataset" | sed -n "${file_number}p")
      # Check if filename was found for the first file
      if [ $i -eq 0 ]; then
        if [ -z "$filename" ]; then
          echo "Error: filename from dataset $dataset_name in $dataset_type could not be found - exiting."
          exit 1
        fi
      fi
      # Adding redirector for local root file
      input_path="file:$dataset/$filename"
      # Adding path to list of filenames, separated by a comma
      filename_list+="$input_path,"
    fi
  done

else
  echo "Error: input_format $input_format not found - exiting."
  exit 1
fi
# Remove the trailing comma
filename_list=${filename_list%,}
echo -- "Running over input files: $filename_list"

# Setting output path
output_name="${dataset_name}_part-${process_id}.root"
output_path="$output_base_path/$output_name"
miniAOD_output_path="$tmp_miniAOD_path/$output_name"

### CMSSW ENVIRONMENT ###
cd "${CMSSW_base_path}/CMSSW_10_6_29/src"
eval `scramv1 runtime -sh`

if [[ $run_mini_and_nano == 1 ]]; then
  if [ ! -d "$miniAOD_output_path" ]; then
    mkdir -p "$miniAOD_output_path"
  fi
  ### RUN LLPminiAOD ###
  cmsRun "${LLP_working_path}/run_LLPminiAOD.py" "inputFiles=$filename_list" "outputFile=$miniAOD_output_path" "nEvents=$n_events" "runOnData=$is_data"
  n_events=0
  filename_list="file:$miniAOD_output_path"
  echo -- Updating input LLPminiAOD file to $filename_list
fi

### RUN LLPNanoAOD ###
cmsRun "${LLP_working_path}/run_LLPNanoAOD.py" "inputFiles=$filename_list" "outputFile=$output_path" "nEvents=$n_events" "runOnData=$is_data"

echo "LLPNanoAOD file saved in: $output_path"
echo "Removing temporary LLPminiAOD file: $miniAOD_output_path"
rm $miniAOD_output_path
echo "Done"
