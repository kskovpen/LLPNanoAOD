from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ttalps_m-1GeV_ctau-1e2mm_LLPminiAOD'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'LLPminiAOD_cfg.py'
config.JobType.pyCfgParams = ['outputFile=test.root','nEvents=0','runOnData=False', 'nThreads=8']
NCORES = 8
config.JobType.maxMemoryMB = 2000 * NCORES
config.JobType.numCores = NCORES

config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/ttalps/lrygaard-crab_ttalps_m-1GeV_ctau-1e2mm_RECO-a5d501e738bc46974ac8d371aaff19e9/USER'
config.Data.outLFNDirBase = '/store/user/lrygaard/ttalps/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.outputDatasetTag = 'ttalps_m-1GeV_ctau-1e2mm_LLPminiAOD'
config.Data.publication = True

config.Site.storageSite = 'T2_DE_DESY'