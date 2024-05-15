# LLPNanoAOD #

LLPNanoAOD is an extension of NanoAOD with parameters useful for analyses with Long-Lived Particles (LLP).

LLPNanoAOD includes variables for:
* #### DisplacedStandAloneMuons (DSAMuon) ####
  from reco::Track displacedStandAloneMuon AOD/miniAOD collections
* #### BeamSpot (BS) ####
  including vertex-type variables
* #### Extended Muon variables (Muon collection) ####
  including impact parameters, indexing, matching to DSAMuon collection
* #### Muon vertices ####
  Muon vertex fits for two muons in combinations:
  * PatMuonVetex = between two Pat muons
  * PatDSAMuonVertex = between one Pat muon and one DSA muon
  * DSAMuonVertex = between two DSA muons
* #### DisplacedGlobalMuons (DGLMuon) ####
  from reco::Track displacedGlobalMuon AOD/miniAOD collections
* #### Muon vertices for DGL muons ####
  Muon vertex fits for two muons, with at least one DGL muon, in combinations:
  * PatDGLMuonVetex = between one Pat muon and one DGL muon
  * DGLDSAMuonVertex = between one DGL muon and one DSA muon
  * DGLMuonVertex = between two DGL muons
* #### Extended GenPart variables (GenPart collection) ####
  including displacent variables vx, vy, vz

A complete list of all branches can be found in `LLPnanoAOD_branches.txt`.

## Setup ##

LLPNanoAOD has been setup for `CMSSW_10_6_29`.

```
cmsrel CMSSW_10_6_29
cd CMSSW_10_6_29/src
cmsenv

git clone git@github.com:kerstinlovisa/LLPNanoAOD.git
scram b -j
```

Due to a bug to determine the charge of high pT tracks, corrections are made to RecoVertex scripts based on later (corrected) CMMSW releases.
There is a an updated version of CheckHitPattern for later CMMSW releases that we want to use:
```
git cms-addpkg RecoVertex/KalmanVertexFit
git cms-addpkg RecoVertex/VertexTools
git cms-addpkg RecoVertex/KinematicFitPrimitives
git cms-addpkg PhysicsTools/RecoUtils

cp LLPNanoAOD/RecoVertex_corrections/VertexTools/src/* RecoVertex/VertexTools/src/
cp LLPNanoAOD/RecoVertex_corrections/VertexTools/interface/* RecoVertex/VertexTools/interface/
cp LLPNanoAOD/RecoVertex_corrections/KalmanVertexFit/src/* RecoVertex/KalmanVertexFit/src/
cp LLPNanoAOD/RecoVertex_corrections/KinematicFitPrimitives/src/* RecoVertex/KinematicFitPrimitives/src/
cp LLPNanoAOD/PhysicsTools_corrections/RecoUtils/src/* PhysicsTools/RecoUtils/src/
cp LLPNanoAOD/PhysicsTools_corrections/RecoUtils/interface/* PhysicsTools/RecoUtils/interface/

scram b -j
```

This will also extend the values of the `trackerBoundsRadius` and `trackerBoundsHalfLength` to:
* `trackerBoundsRadius = 740` 
* `trackerBoundsHalfLength = 960`
  
In `RecoVertex/VertexTools/interface/SequencialVertexFitter.h` and `RecoVertex/KalmanVertexFit/src/SingleTrackVertexConstraint.cc`.

## Settings for running ##

Main settings while running is set in `LLPNanoAOD/LLPnanoAOD/test/run_LLPnanoAOD_config.py` with parameters:

* `output_base_path`: base path to output files
* `nEvents`: number of files to run on, if `nEvents = 0` it will run over all events in the file
* `filesPerJob`: number of files to run on per job
* `maxJobs`: max number of jobs, if `maxJobs = 0` it will create as many jobs as needed for number of dataset files
* `run_mini_and_nano`: if `run_mini_and_nano = True` it will first run `LLPminiAOD_cfg.py` (with necessary keep-statements) and then `LLPnanoAOD_cfg.py`, input should be AOD-files. if `run_mini_and_nano = False` it will only run `LLPnanoAOD_cfg.py`, input should be miniAOD-files.
* `save_mini`: if `run_mini_and_nano = True` you can decide if you want to save the LLPminiAOD-files after finishing running `LLPnanoAOD_cfg.py`
* `LLPcollections`: list of string with all extra LLPnanoAOD collections you would like to include. The options are:
  * `DSAMuon`: include DSAMuon collection, extended Muon collection and vertices collections: PatVertex, PatDSAVertex, DSAVertex
  * `BS`: include BeamSpot collection
  * `GenPart`: include extended GenPart collection
  * `DGLMuon`: include DGLMuon collection and vertices collections: PatDGLVertex, DGLDSAVertex, DGLVertex

* `datasets`: dictionary of input datasets in the format: `( output_dataset_path : dataset_name )`
  * `output_dataset_path`: this means that the complete output will be `output_base_path/output_dataset_path/`
  * `dataset_name`:
    * for DAS datasets: `dataset_name` = complete DAS dataset name
    * for local Datasets: `dataset_name` = path to directory with root-files
  Example for 2018 SingleMuon dataset A:
  ```
  datasets = (
    ("data2018/SingleMuonA", "/SingleMuon/Run2018A-12Nov2019_UL2018-v5/AOD"),
  )
  ```

## Run locally ##

To run locally use the run script `run_LLPnanoAOD.py` and the config file `run_LLPnanoAOD_config.py`:

Run from `$CMMSW_BASE/src`:
```
python LLPNanoAOD/LLPnanoAOD/test/run_LLPnanoAOD.py --config LLPNanoAOD/LLPnanoAOD/test/run_LLPnanoAOD_config.py --local
```

Make sure to update the config `run_LLPnanoAOD_config.py` with all settings needed to run on your computer.

## Run on condor ##

Run from `$CMMSW_BASE/src`:
```
python LLPNanoAOD/LLPnanoAOD/test/run_LLPnanoAOD.py --config LLPNanoAOD/LLPnanoAOD/test/run_LLPnanoAOD_config.py --condor
```

### Parameter options: ###
* `--job_flavour`: job flavour for condor job. Use help to see all flavour options. Default: "longlunch" = 2h
* `--run_time`: RequestRunTime for condor. If run_time is set it will overwrite the job_flavour input. Default is 0 - uses job_flavour.
* `--dry`: for condor dry run without submission

Make sure to update the config `run_LLPnanoAOD_config.py` with all settings needed to run on your computer.

## Crab submission ##
Crab submissions are also setup to store output on a tier2 site. See README in `LLPNanoAOD/LLPnanoAOD/test/crab` for more instructions.
