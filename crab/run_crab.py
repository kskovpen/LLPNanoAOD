from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'LLPminiAOD_SingleMuonD'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_MiniAOD.py'
# config.JobType.inputFiles = ['root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18RECO/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/10014/A7C63886-4306-4B40-95EA-4DA41FAA1C8E.root']
config.JobType.pyCfgParams = ['outputFile=test.root','nEvents=10','runOnData=1']

config.Data.inputDataset = '/SingleMuon/Run2018D-12Nov2019_UL2018-v8/AOD'
config.Data.outLFNDirBase = '/store/user/lrygaard/TTALPs/data2018/SingleMuonA/LLPminiAOD/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 1  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_test'

config.Site.storageSite = 'T2_DE_DESY'