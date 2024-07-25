#!/bin/bash

# Debugging output
echo "Hostname: $(hostname)"
echo "Current directory: $(pwd)"
echo "Initial SCRAM_ARCH: $SCRAM_ARCH"
echo "Initial CMSSW_BASE: $CMSSW_BASE"
echo "Initial LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

export X509_USER_PROXY=`pwd`/voms_proxy

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
# export SCRAM_ARCH=slc7_amd64_gcc700
export SCRAM_ARCH=slc7_amd64_gcc11
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git

echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

cd $CMSSW_BASE
eval `scramv1 runtime -sh`

# Check if /cvmfs is mounted
if [ ! -d /cvmfs/cms.cern.ch ]; then
  echo "/cvmfs/cms.cern.ch is not available"
  exit 1
fi

export LD_LIBRARY_PATH=/afs/desy.de/user/l/lrygaard/TTALP/CMSSW_13_0_13/biglib/slc7_amd64_gcc11:/afs/desy.de/user/l/lrygaard/TTALP/CMSSW_13_0_13/lib/slc7_amd64_gcc11:/afs/desy.de/user/l/lrygaard/TTALP/CMSSW_13_0_13/external/slc7_amd64_gcc11/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc11/cms/cmssw/CMSSW_13_0_13/biglib/slc7_amd64_gcc11:/cvmfs/cms.cern.ch/slc7_amd64_gcc11/cms/cmssw/CMSSW_13_0_13/lib/slc7_amd64_gcc11:/cvmfs/cms.cern.ch/slc7_amd64_gcc11/cms/cmssw/CMSSW_13_0_13/external/slc7_amd64_gcc11/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc11/external/llvm/12.0.1-476112d9475f69ecba350f77e3ec4975/lib64:/cvmfs/cms.cern.ch/slc7_amd64_gcc11/external/gcc/11.2.1-f9b9dfdd886f71cd63f5538223d8f161/lib64:/cvmfs/cms.cern.ch/slc7_amd64_gcc11/external/gcc/11.2.1-f9b9dfdd886f71cd63f5538223d8f161/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc11/external/cuda/11.5.2-66a9473808e7d5863d5bbec0824e2c4a/lib64/stubs:/cvmfs/grid.cern.ch/centos7-umd4-ui-211021/lib64:/cvmfs/grid.cern.ch/centos7-umd4-ui-211021/lib:/cvmfs/grid.cern.ch/centos7-umd4-ui-211021/usr/lib64:/cvmfs/grid.cern.ch/centos7-umd4-ui-211021/usr/lib:/cvmfs/grid.cern.ch/centos7-umd4-ui-211021/usr/lib64/dcap:/afs/desy.de/user/l/lrygaard/tools/MG5_aMC_v3_4_2/HEPTools/lhapdf6_py3/lib

# Debugging output after setting up the CMSSW environment
echo "After setup:"
echo "SCRAM_ARCH: $SCRAM_ARCH"
echo "CMSSW_BASE: $CMSSW_BASE"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "PATH: $PATH"
echo "PYTHONPATH: $PYTHONPATH"

job_number=$1
echo "Executing job number $job_number"

home_path="<home_path>"
export HOME=$home_path

dataset=<dataset>

inputFiles_path="${CMSSW_BASE}/src/<input_files_path>"
all_files=()
while IFS= read -r filename; do
  all_files+=("$filename")
done < "$inputFiles_path"

output="<output_path>"
outputPath="${output}_part-${job_number}.root"

LLPminiAOD="<LLPminiAOD_path>"
LLPminiAOD_path="${LLPminiAOD}_part-${job_number}.root"

nEvents=<n_events>
nFiles=<n_files>
runOnData=<is_data>
runOnDAS=<run_on_das>

includeDSAMuon=<include_DSAMuon>
includeBS=<include_BS>
includeGenPart=<include_GenPart>
includeDGLMuon=<include_DGLMuon>
includeRefittedTracks=<include_refittedTracks>

saveLLPminiAOD=<save_LLPminiAOD>

year="<year>"

miniAODrunFile=<miniAOD_runfile>
nanoAODrunFile=<nanoAOD_runfile>

total_files=${#all_files[@]}

for ((i=0; i<$nFiles; i++)); do
  file_number=$((job_number*nFiles+i))
  # Check that the file number doesn't exceed total number of files
    if [ $file_number -le $total_files ]; then 

      filename=${all_files[file_number]}
      if [ $i -eq 0 ]; then
        if [ -z "$filename" ]; then
          echo "Error: filename from dataset $dataset could not be found - exiting."
          exit 1
        fi
      fi
      if [[ $runOnDAS == "True" ]]; then
        # input_path="root://xrootd-cms.infn.it/$filename"
        input_path="root://cms-xrd-global.cern.ch/$filename"
      else
        input_path="file:$dataset/$filename"
      fi
      filename_list+="$input_path,"
    fi
done
# Remove the trailing comma
filename_list=${filename_list%,}

cd <work_dir>

echo cmsRun $CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test/$miniAODrunFile "inputFiles=$filename_list" "outputFile=$LLPminiAOD_path" "nEvents=$nEvents" "runOnData=$runOnData" "year=$year"
cmsRun $CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test/$miniAODrunFile "inputFiles=$filename_list" "outputFile=$LLPminiAOD_path" "nEvents=$nEvents" "runOnData=$runOnData" "year=$year"

echo cmsRun $CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test/$nanoAODrunFile "inputFiles=file:$LLPminiAOD_path" "outputFile=$outputPath" "nEvents=$nEvents" "runOnData=$runOnData" "includeDSAMuon=$includeDSAMuon" "includeBS=$includeBS" "includeGenPart=$includeGenPart" "includeDGLMuon=$includeDGLMuon" "includeRefittedTracks=$includeRefittedTracks" "year=$year"
cmsRun $CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test/$nanoAODrunFile "inputFiles=file:$LLPminiAOD_path" "outputFile=$outputPath" "nEvents=$nEvents" "runOnData=$runOnData" "includeDSAMuon=$includeDSAMuon" "includeBS=$includeBS" "includeGenPart=$includeGenPart" "includeDGLMuon=$includeDGLMuon" "includeRefittedTracks=$includeRefittedTracks" "year=$year"

echo "LLPNanoAOD file saved in: $outputPath"
if [[ $saveLLPminiAOD == "False" ]]; then
  echo rm $LLPminiAOD_path
  rm $LLPminiAOD_path
fi
echo "Done"