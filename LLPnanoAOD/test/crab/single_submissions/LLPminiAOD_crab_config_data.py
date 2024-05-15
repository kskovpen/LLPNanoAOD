from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SingleMuonD_LLPminiAOD'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'LLPminiAOD_cfg.py'
# config.JobType.inputFiles = ['root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18RECO/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/10014/A7C63886-4306-4B40-95EA-4DA41FAA1C8E.root']
config.JobType.pyCfgParams = ['outputFile=test.root','nEvents=10','runOnData=True', 'nThreads=8']
NCORES = 8
config.JobType.maxMemoryMB = 2000 * NCORES
config.JobType.numCores = NCORES

config.Data.inputDataset = '/SingleMuon/Run2018D-12Nov2019_UL2018-v8/AOD'
config.Data.outLFNDirBase = '/store/user/lrygaard/ttalps/'
# config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
# NJOBS = 1  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
# config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True
# config.Data.outputDatasetTag = 'SingleMuonD_LLPminiAOD'

config.Site.storageSite = 'T2_DE_DESY'