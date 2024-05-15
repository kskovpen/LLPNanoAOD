from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SingleMuonC_LLPnanoAOD'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'LLPnanoAOD_cfg.py'
config.JobType.pyCfgParams = ['outputFile=test.root','nEvents=10','runOnData=True', 'nThreads=8']
NCORES = 8
config.JobType.maxMemoryMB = 2000 * NCORES
config.JobType.numCores = NCORES

config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/SingleMuon/lrygaard-crab_SingleMuonC_LLPminiAOD-9cdbfc999b77f606d32dabc67655eebd/USER'
config.Data.outLFNDirBase = '/store/user/lrygaard/ttalps/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
config.Data.ignoreLocality = True

config.Site.whitelist = ['T2_DE_DESY', 'T2_IT_Bari', 'T2_CH_CSCS', 'T2_EE_Estonia', 'T2_UK_London_IC', 'T2_BE_IIHE', 'T2_CH_CERN', 'T2_US_Wisconsin', 'T2_US_Vanderbilt', 'T2_US_MIT', 'T2_DE_RWTH', 'T2_US_Nebraska', 'T3_KR_KNU']
config.Site.storageSite = 'T2_DE_DESY'
