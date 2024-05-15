from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DYJetsToMuMu_M-10to50_LLPminiAOD'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'LLPminiAOD_cfg.py'
config.JobType.pyCfgParams = ['outputFile=test.root','nEvents=0','runOnData=False', 'nThreads=8']
NCORES = 8
config.JobType.maxMemoryMB = 2000 * NCORES
config.JobType.numCores = NCORES

config.Data.inputDataset = '/DYJetsToMuMu_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM'
config.Data.outLFNDirBase = '/store/user/lrygaard/ttalps/'
config.Data.splitting = 'FileBased'
# config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 1
# NJOBS = 10  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
# config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True

config.Site.storageSite = 'T2_DE_DESY'