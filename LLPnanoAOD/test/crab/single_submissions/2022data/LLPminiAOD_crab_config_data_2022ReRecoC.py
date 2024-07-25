from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = '/afs/desy.de/user/l/lrygaard/TTALP/CMSSW_13_0_13/src/LLPNanoAOD/LLPnanoAOD/test/crab/crab_projects/crab_data_2022ReReco_AOD_v1'
config.General.requestName = 'SingleMuonC_LLPminiAODv1'
config.section_('JobType')
config.JobType.numCores = 4
config.JobType.sendExternalFolder = False
config.JobType.pyCfgParams = ['nEvents=0', 'runOnData=True', 'nThreads=4', 'year=2022ReReco']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '/afs/desy.de/user/l/lrygaard/TTALP/CMSSW_13_0_13/src/LLPNanoAOD/LLPnanoAOD/test/LLPminiAOD_Run3_cfg.py'
config.JobType.maxJobRuntimeMin = 2750
config.JobType.maxMemoryMB = 4000
config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2022C-27Jun2023-v1/AOD'
config.Data.outputDatasetTag = 'LLPminiAODv1_Run2022C-27Jun2023_v1'
config.Data.publication = True
config.Data.unitsPerJob = 1
# config.Data.ignoreLocality = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '/store/user/lrygaard/ttalps'
config.section_('Site')
# config.Site.whitelist = ['T2_DE_*', 'T1_DE_*', 'T2_CH_*', 'T2_IT_*', 'T1_IT_*', 'T2_ES_*', 'T1_ES_*', 'T2_UK_*', 'T1_UK_*', 'T2_BR_*', 'T1_BR_*', 'T2_RU_*', 'T1_RU_*']
config.Site.storageSite = 'T2_DE_DESY'
config.section_('User')
config.section_('Debug')