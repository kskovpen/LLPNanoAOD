# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_upgrade2018_realistic_v16_L1v1 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_106Xv2 --filein file:step-1.root --fileout file:step0.root
import os
import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Modifier_run2_nanoAOD_106Xv2_cff import run2_nanoAOD_106Xv2

# Input arguments
options = VarParsing('analysis')
options.outputFile = 'output.root'
options.inputFiles = 'root://cms-xrd-global.cern.ch///store/data/Run2018A/SingleMuon/AOD/12Nov2019_UL2018-v5/270000/D437BC65-B40B-1944-8DD0-9710C088B916.root'
options.register('nEvents',
                    1000,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.int,
                    "Number of events to process"
                )
options.register('runOnData',
                    False,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.bool,
                    "If running on data"
                )
options.register('includeDSAMuon',
                    False,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.bool,
                    "Flag to include Displaced StandAlone muon information"
                )
options.register('includeBS',
                    False,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.bool,
                    "Flag to include BeamSpot information"
                )
options.register('includeGenPart',
                    False,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.bool,
                    "Flag to include extended GenPart information"
                )
options.register('includeDGLMuon',
                    False,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.bool,
                    "Flag to include Displaced GLobal Muon information"
                )
options.register('includeRefittedTracks',
                    True,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.bool,
                    "Flag to include Refitted Tracks information of Muon Vertex fits"
                )
options.register('nThreads',
                    1,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.int,
                    "Number of threads to use"
                )
options.register('year',
                    2018,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    "Year of the dataset"
                )
options.parseArguments()

runRefittedTracks_=options.includeRefittedTracks
nevents = options.nEvents
if nevents == 0:
    nevents=-1
    
print('Running run_LLPNanoAOD.py')
# print('-- Input MiniAOD files: ')
# for file in options.inputFiles:
#     print(file)
print('-- Output NanoAOD file: '+options.outputFile)
print('-- Running on '+str(nevents)+' number of events')
if options.runOnData:
    print('-- Running on data')
else:
    print('-- Running on mc')

process = cms.Process('NANO',Run2_2018,run2_nanoAOD_106Xv2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 100000  # Set reportEvery to control the frequency of report messages
process.MessageLogger.threshold = cms.untracked.string('ERROR')  # Set the output threshold to ERROR

if not options.runOnData:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevents)
)

# Input source
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

datatier = 'NANOAOD'
outputcommands = process.NANOAODSIMEventContent.outputCommands
if options.runOnData:
    datatier = 'NANOAODSIM'
    outputcommands = process.NANOAODEventContent.outputCommands

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(datatier),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:{0}'.format(options.outputFile)),
    outputCommands = outputcommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
if options.runOnData:
    if options.year == "2016":
        globalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v36', '')
    if options.year == "2017":
        globalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')
    if options.year == "2018":
        globalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')
    if options.year == "2022ReReco":
        globalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_v2', '')
    if options.year == "2022Prompt":
        globalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '')
    if options.year == "2023":
        globalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '')
else:
    if options.year == "2016":
        globalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v17', '')
    if options.year == "2017":
        globalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v9', '')
    if options.year == "2018":
        globalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')
    if options.year == "2022PreEE":
        globalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v5', '')
    if options.year == "2022PostEE":
        globalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_postEE_v6', '')
    if options.year == "2023PreBPix":
        globalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_v14', '')
    if options.year == "2023PostBPix":
        globalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_postBPix_v2', '')
process.GlobalTag = globalTag


# LLPnanoAOD custom producers
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi import *

# DisplacedStandAlone (DSA) Muon table
process.dSAMuonsTable = cms.EDProducer("DSAMuonTableProducer",
    dsaMuons=cms.InputTag("displacedStandAloneMuons"),
    muons=cms.InputTag("linkedObjects","muons"),
    primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamspot = cms.InputTag("offlineBeamSpot")
)
# DisplacedGLobal (DGL) Muon table
process.dGlMuonsTable = cms.EDProducer("DGLMuonTableProducer",
    dglMuons=cms.InputTag("displacedGlobalMuons"),
    muons=cms.InputTag("linkedObjects","muons"),
    primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamspot = cms.InputTag("offlineBeamSpot")
)
# BeamSpot table
process.beamSpotTable = cms.EDProducer("BeamSpotTableProducer",
    beamSpot = cms.InputTag("offlineBeamSpot")
)
# GenPart table
process.genPartExtendedTable = cms.EDProducer("GenParticlesExtendedTableProducer",
    genparticles = cms.InputTag("finalGenParticles")
)
# Vertex between two muons (pat-pat, pat-dsa or dsa-dsa)
process.muonVertexTable = cms.EDProducer("MuonVertexTableProducer",
    dsaMuons=cms.InputTag("displacedStandAloneMuons"),
    patMuons=cms.InputTag("linkedObjects","muons"),
    beamspot=cms.InputTag("offlineBeamSpot"),
    generalTracks=cms.InputTag("generalTracks"),
    primaryVertex=cms.InputTag("offlineSlimmedPrimaryVertices")
    # runRefittedTracks=cms.bool(runRefittedTracks_)
)
# Vertex between two muons (pat-pat, pat-dgl or dgl-dgl)
process.dGlMuonVertexTable = cms.EDProducer("MuonVertexTableProducer",
    dsaMuons=cms.InputTag("displacedGlobalMuons"),
    patMuons=cms.InputTag("linkedObjects","muons"),
    beamspot=cms.InputTag("offlineBeamSpot"),
    generalTracks=cms.InputTag("generalTracks"),
    primaryVertex=cms.InputTag("offlineSlimmedPrimaryVertices")
    # runRefittedTracks=cms.bool(runRefittedTracks_)
)
# LLPnanoAOD Muon extended table
process.muonExtendedTable = cms.EDProducer("MuonExtendedTableProducer",
    name=cms.string("Muon"),
    muons=cms.InputTag("linkedObjects","muons"),
    dsaMuons=cms.InputTag("displacedStandAloneMuons"),
    primaryVertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamspot=cms.InputTag("offlineBeamSpot"),
    generalTracks=cms.InputTag("generalTracks")
)

# Path and EndPath definitions
if options.runOnData:
    process.nanoAOD_step = cms.Path(process.nanoSequence)
else:
    process.nanoAOD_step = cms.Path(process.nanoSequenceMC)

if options.includeDSAMuon:
    process.nanoAOD_step += process.dSAMuonsTable
    process.nanoAOD_step += process.muonExtendedTable
    process.nanoAOD_step += process.muonVertexTable
if options.includeDGLMuon:
    process.nanoAOD_step += process.dGlMuonsTable
    process.nanoAOD_step += process.dGlMuonVertexTable
if options.includeBS:
    process.nanoAOD_step += process.beamSpotTable
if options.includeGenPart and not options.runOnData:
    process.nanoAOD_step += process.genPartExtendedTable

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.nThreads)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# customisation of the process.

# customisation of the process.
if options.runOnData:
    # Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData 
    #call to customisation function nanoAOD_customizeData imported from PhysicsTools.NanoAOD.nano_cff
    process = nanoAOD_customizeData(process)
else:
    # Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 
    #call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
    process = nanoAOD_customizeMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
