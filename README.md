# LLPNanoAOD #

LLPNanoAOD is an extension of NanoAOD with parameters useful for analyses with Long-Lived Particles (LLP).

LLPNanoAOD includes variables for:
* DisplacedStandAloneMuons (DSAMuon)
* BeamSpot (BS)
* Extended Muon viarables (Muon)
* Muon vertices (MuonVertex, MuonCombVertex, DSAMuonVertex)

## run_LLPNanoAOD ##

Several files are provided to run LLPNanoAOD, which will be run in order from top to bottom:
```
run_LLPNanoAOD.sub (for condor)
  run_LLPNanoAOD.sh
    LLPNanoAOD_cfg.py
```

## Setup CMSSW environment for LLPNanoAOD ##

Using CMSSW_10_6_29 release:

```
cmsrel CMSSW_10_6_29
cd CMSSW_10_6_29/src
cmsenv
git cms-addpkg PhysicsTools/NanoAOD
```

### Include LLP table producers: ###

#### Include producers: ####

Make sure to include the following plugins in PhysicsTools/NanoAOD/plugins:
* `DSAMuonTableProcuder.cc`
* `BeamSpotTableProcuder.cc`
* `MuonVertexTableProducer.cc`
And an updated build file:
* `BuildFile.xml` in the same plugin directory

Make sure to include the following python scripts in PhysicsTools/NanoAOD/python:
*  llpnano_cff.py

In your setup CMSSW_10_6_29 directory:
```
cd CMSSW_10_6_29/src
cmsenv
cp {LLPNanoAOD base path}/PhysicsTools/NanoAOD/plugins/*.cc PhysicsTools/NanoAOD/plugins/.
cp {LLPNanoAOD base path}/PhysicsTools/NanoAOD/plugins/BuildFile.xml PhysicsTools/NanoAOD/plugins/.
scram b -j
```

## Setup user parameters for LLPNanoAOD ##

Create your own `userparams.json` config file in `UserParameters`, following the example in `UserParameters/example_userparams.json`.

Settings that has to be set:
* `output_base_path`: desired output path
* `home_path`: path to local home directory, this is to export the home path when running on condor to run root modules
* `CMSSW_base_path`: path to the directory where your CMSSW_10_6_29 release (see above) is stored
* `LLP_working_path`: path to your local LLP directory
* `dasgoclient_path`: path to dasgoclient
* `proxy_path`: path to private proxy for connecting to das

### Datasets ###
The dataset you want to use as input has to be given in your userparams in the nested dictionary `datasets`:

```
"datasets" = {
  dataset_type: {
    dataset_name: {
      "dataset": dataset path,
      "input_format": "das" for input from DAS or "local" for locally stored file,
      "files": 0,
    },
  },
}
```
The `dataset_type` and `dataset_name` will also be used as the name of the output directory. The full output path will be:
```
output_base_path/dataset_type/dataset_name/LLPNanoAOD/filename.root
```
The variable `files` is not used but given for your own use to store the information if wanted. 

Example for TTToSemiLeptonic Summer20UL18 file on DAS:
```
"datasets" = {
  "backgrounds18": {
    "TTToSemiLeptonic": {
        "dataset": "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
        "input_format": "das",
        "files": 10010,
    },
  },
}
```
