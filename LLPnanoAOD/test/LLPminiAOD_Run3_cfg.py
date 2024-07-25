import os
import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
from Configuration.Eras.Modifier_run3_miniAOD_12X_cff import run3_miniAOD_12X

# Input arguments
options = VarParsing('analysis')
options.outputFile = 'output.root'
options.inputFiles = 'root://xrootd-cms.infn.it//store/data/Run2022C/SingleMuon/AOD/27Jun2023-v1/2820000/ae0a7cb5-8034-4035-8fb9-07bd54a1bd73.root'
options.register('nEvents',
                    1000,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.int,
                    "Number of events to process"
                )
options.register('runOnData',
                    True,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.bool,
                    "If running on data"
                )
options.register('nThreads',
                    1,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.int,
                    "Number of threads to use"
                )
options.register('year',
                    "2022ReReco",
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    "Year of the dataset"
                )
options.parseArguments()

nevents = options.nEvents
if nevents == 0:
    nevents=-1

print('Running LLPminiAOD_Run3_cfg.py')
print('-- Running on '+str(nevents)+' number of events')
if options.runOnData:
    print('-- Running on data')
else:
    print('-- Running on mc')
print('-- Year: ' + options.year)

# load process for correct year and data/MC
if options.runOnData:
    if options.year == "2022ReReco":
        process = cms.Process('PAT',Run3,run3_miniAOD_12X)
    if options.year == "2022Prompt":
        process = cms.Process('PAT',Run3,run3_miniAOD_12X)
    if options.year == "2023":
        process = cms.Process('PAT',Run3_2023)
        
else:
    if options.year == "2022PreEE":
        process = cms.Process('PAT',Run3)
    if options.year == "2022PostEE":
        process = cms.Process('PAT',Run3)
    if options.year == "2023PreBPix":
        process = cms.Process('PAT',Run3_2023) 
    if options.year == "2023PostBPix":
        process = cms.Process('PAT',Run3_2023) 

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 100000  # Set reportEvery to control the frequency of report messages
# process.MessageLogger.threshold = cms.untracked.string('ERROR')  # Set the output threshold to ERROR

if options.runOnData:
    process.load('Configuration.StandardSequences.PAT_cff')
else:
    process.load('Configuration.StandardSequences.PATMC_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
datatier = 'MINIAODSIM'
outputcommands = process.MINIAODSIMEventContent.outputCommands
if options.runOnData:
    datatier = 'MINIAOD'
    outputcommands = process.MINIAODEventContent.outputCommands

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(datatier),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:'+options.outputFile),
    outputCommands = outputcommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.MINIAODSIMoutput.outputCommands.append('keep *_generalTracks_*_*')
process.MINIAODSIMoutput.outputCommands.append('keep recoTrackExtras_muonReducedTrackExtras__RECO')
process.MINIAODSIMoutput.outputCommands.append('keep recoTracks_displacedStandAloneMuons__RECO')
process.MINIAODSIMoutput.outputCommands.append('keep recoTrackExtras_displacedStandAloneMuons__RECO') 
process.MINIAODSIMoutput.outputCommands.append('keep TrackingRecHitsOwned_displacedStandAloneMuons__RECO')
process.MINIAODSIMoutput.outputCommands.append('keep recoTrackExtras_globalMuons__RECO')
process.MINIAODSIMoutput.outputCommands.append('keep recoTrackExtras_standAloneMuons__RECO')


# Path and EndPath definitions
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_hfNoisyHitsFilter = cms.Path(process.hfNoisyHitsFilter)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.Flag_BadPFMuonDzFilter = cms.Path(process.BadPFMuonDzFilter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
globalTag = ""
if options.runOnData:
    if options.year == "2022ReReco":
        globalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_v2', '')
    if options.year == "2022Prompt":
        globalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '')
    if options.year == "2023":
        globalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '')
else:
    if options.year == "2022PreEE":
        globalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v5', '')
    if options.year == "2022PostEE":
        globalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_postEE_v6', '')
    if options.year == "2023PreBPix":
        globalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_v14', '')
    if options.year == "2023PostBPix":
        globalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_postBPix_v2', '')
print("year: "+options.year)
process.GlobalTag = globalTag

# Schedule definition
process.schedule = cms.Schedule(process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadPFMuonDzFilter,process.Flag_hfNoisyHitsFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,process.MINIAODSIMoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.nThreads)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

if options.runOnData:
    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 
    process = miniAOD_customizeAllData(process)
else:
    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 
    process = miniAOD_customizeAllMC(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
