from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = '/afs/desy.de/user/l/lrygaard/TTALP/CMSSW_13_0_13/src/LLPNanoAOD/LLPnanoAOD/test/crab/crab_projects/crab_bkg_2023PreBPix_AOD_v1'
config.General.requestName = 'TTtoLNu2Q_LLPminiAODv1'
config.section_('JobType')
config.JobType.numCores = 4
config.JobType.sendExternalFolder = False
config.JobType.pyCfgParams = ['nEvents=0', 'runOnData=False', 'nThreads=4', 'year=2023PreBPix']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '/afs/desy.de/user/l/lrygaard/TTALP/CMSSW_13_0_13/src/LLPNanoAOD/LLPnanoAOD/test/LLPminiAOD_Run3_cfg.py'
config.JobType.maxJobRuntimeMin = 2750
config.JobType.maxMemoryMB = 4000
config.section_('Data')
config.Data.inputDataset = '/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23DRPremix-130X_mcRun3_2023_realistic_v14-v2/AODSIM'
config.Data.outputDatasetTag = 'LLPminiAODv1_Run3Summer23DRPremix-130X_mcRun3_2023_realistic_v14-v2'
config.Data.publication = True
config.Data.unitsPerJob = 3
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