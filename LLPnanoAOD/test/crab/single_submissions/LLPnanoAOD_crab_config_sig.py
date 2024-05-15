from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ttalps_m-0p35GeV_ctau-1e-5mm_LLPnanoAODv1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'LLPnanoAOD_cfg.py'
config.JobType.pyCfgParams = ['nEvents=0','runOnData=False', 'nThreads=4', 'includeDSAMuon=True', 'includeBS=True', 'includeGenPart=True', 'includeDGLMuon=False']
NCORES = 4
config.JobType.maxMemoryMB = 2000 * NCORES
config.JobType.numCores = NCORES

config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/ttalps/lrygaard-ttalps_m-0p35GeV_ctau-1e-5mm_LLPminiAOD-c15273f0b6812ff053a850f456209388/USER'
config.Data.outLFNDirBase = '/store/user/lrygaard/ttalps/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outputDatasetTag = 'ttalps_m-0p35GeV_ctau-1e-5mm_LLPnanoAODv1'
config.Data.publication = True

config.Site.storageSite = 'T2_DE_DESY'