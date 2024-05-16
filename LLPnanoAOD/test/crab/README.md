# Run LLPminiAOD and LLPnanoAOD with crab #

Note: all functionality in crab.py has not yet been implemented - only use crab submit feature as in bash scripts

## File structure ##

### Logical script ###
Main functionalities of the crab submission and creation of crab configs are done in `crab.py` DNNTuples:
https://github.com/SohamBhattacharya/DNNTuples/blob/dca9f695499aa9e56932cac073e7bca2536c6deb/Ntupler/run/crab.py

### Run scripts ###
To run the script use the help bash scripts `run_LLPminiAOD_crab.sh` and `run_LLPnanoAOD_crab.sh`.

### Input files ###
Input files are stored in `inputs` as `.conf`-files. The already given example input files are named as `filetype_year_datatier.conf`.
The input files should be structured as:
- One dataset per line
- In each line given as: datasetname datasetpath
    example:
    `TTToSemiLeptonic /TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM`

### Crab projects ###
Crab projects will be created in the directory: `$CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test/crab/crab_projects`.
Project directory will be named after input file name and version, example: `crab_bkg_2018_AOD_v1`, and inside will be a directory for each dataset given in the input file.

## Run LLPminiAOD ##
Run LLPminiAOD crab submission with `run_LLPminiAOD_crab.sh` with an input file as argument:
Example:
```
cd $CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test/crab
source /cvmfs/cms.cern.ch/common/crab-setup.sh
./run_LLPminiAOD_crab.sh inputs/bkg_2018_AOD.conf
```
### Settings ####
Make sure to update the settings in `run_LLPminiAOD_crab.sh`:
- `nCores`: number of cores to use
- `maxMemory`: max number of memory calculated as: memory per core * `nCores`
- `maxRuntime`: maximum run time
- `filePerJob`: number of files to run per job for FileBased splitting - make sure to update for large datasets, max jobs is 10000 or set to 0 to automatically set smallest number of files for max 10000 jobs
- `VERSION`: LLPminiAOD version - used for output tag
- `--site`: Site where you have write permission
- `-o`: output directory on site
- `-t`: output tag
- `-i`: input filepath
- `-s`: splitting
- `-n` = `filePerJob`
- `--num-cores` = `nCores`
- `--max-memory` = `maxMemory`
- `--max-runtime-min` = `maxRuntime`
- `--work-area`: location of crab projects files
- `--publication`: include this to published the dataset
- `--input-DBS`: input dbs instance, default is 'global'
- `--test`: to test run on a small file - see run tests section
- `--dry-run`: for dry run

## Run LLPnanoAOD ##
Run LLPnanoAOD crab submission with `run_LLPnanoAOD_crab.sh` with an input file as argument:
Example:
```
cd $CMSSW_BASE/src/LLPNanoAOD/LLPnanoAOD/test/crab
source /cvmfs/cms.cern.ch/common/crab-setup.sh
./run_LLPnanoAOD_crab.sh inputs/bkg_2018_AOD.conf
```
### Settings ####
Make sure to update the settings in `run_LLPnanoAOD_crab.sh`:
- `nCores`: number of cores to use
- `maxMemory`: max number of memory calculated as: memory per core * `nCores`
- `maxRuntime`: maximum run time
- `filePerJob`: number of files to run per job for FileBased splitting - make sure to update for large datasets, max jobs is 10000
- `VERSION`: LLPnanoAOD version - used for output tag
- `--site`: Site where you have write permission
- `-o`: output directory on site
- `-t`: output tag
- `-i`: input filepath
- `-s`: splitting
- `-n` = `filePerJob`
- `--num-cores` = `nCores`
- `--max-memory` = `maxMemory`
- `--max-runtime-min` = `maxRuntime`
- `--work-area`: location of crab projects files
- `--publication`: include this to published the dataset
- `--includeDSAMuon`: LLPnanoAOD flag - to include DSAMuon Collections
- `--includeBS`: LLPnanoAOD flag - to include Beam Spot Collection
- `--includeGenPart`: LLPnanoAOD flag - to include GenPart collection
- `--input-DBS`: input dbs instance, default is 'global'
- `--test`: to test run on a small file - see run tests section
- `--dry-run`: for dry run

## Run test ##
If `--test` parameter is given it will only do a test run with 1 job with 1 file running over 10 events. Furthermore:
- the output tag will be the tagname+'_test'
- the output project dir is named work_area+'_test'
- publication is set to false
