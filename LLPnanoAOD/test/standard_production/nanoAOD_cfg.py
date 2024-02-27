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
options.register('nEvents',
                    '',
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.int,
                    "Number of events to process"
                )
options.register('runOnData',
                    '',
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.bool,
                    "If running on data"
                )
options.parseArguments()

nevents = options.nEvents
if nevents == 0:
    nevents=-1
    
print('Running run_nanoAOD.py')
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
globalTag = ""
if options.runOnData:
    globalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')
else:
    globalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')
process.GlobalTag = globalTag


# # LLPnanoAOD custom producers

# # DisplacedStandAloneMuon table
# process.dSAMuonsTable = cms.EDProducer("DSAMuonTableProducer",
#     dsaMuons=cms.InputTag("displacedStandAloneMuons"),
#     muons=cms.InputTag("linkedObjects","muons"),
#     primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
#     beamspot = cms.InputTag("offlineBeamSpot")
# )
# # BeamSpot table
# process.beamSpotTable = cms.EDProducer("BeamSpotTableProducer",
#     beamSpot = cms.InputTag("offlineBeamSpot")
# )
# # GenPart table
# process.beamSpotTable = cms.EDProducer("GenParticlesExtendedTableProducer",
#     genparticles = cms.InputTag("finalGenParticles")
# )
# # Vertex between two muons (pat-pat, pat-dsa or dsa-dsa)
# process.muonVertexTable = cms.EDProducer("MuonVertexTableProducer",
#     dsaMuons=cms.InputTag("displacedStandAloneMuons"),
#     patMuons=cms.InputTag("linkedObjects","muons"),
#     beamspot=cms.InputTag("offlineBeamSpot"),
#     generalTracks=cms.InputTag("generalTracks")
# )
# # LLPnanoAOD Muon extended table
# process.muonExtendedTable = cms.EDProducer("MuonExtendedTableProducer",
#     name=cms.string("Muon"),
#     muons=cms.InputTag("linkedObjects","muons"),
#     dsaMuons=cms.InputTag("displacedStandAloneMuons"),
#     primaryVertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
#     beamspot=cms.InputTag("offlineBeamSpot"),
#     generalTracks=cms.InputTag("generalTracks")
# )

# # Path and EndPath definitions
# if options.runOnData:
#     process.nanoAOD_step = cms.Path(process.nanoSequence
#                                     +process.dSAMuonsTable
#                                     +process.muonExtendedTable
#                                     +process.beamSpotTable
#                                     +process.muonVertexTable
#                                     +process.muonVertexTable
#                                     )
# else:
#     process.nanoAOD_step = cms.Path(process.nanoSequenceMC
#                                 +process.dSAMuonsTable
#                                 +process.muonExtendedTable
#                                 +process.beamSpotTable
#                                 +process.muonVertexTable
#                                 )
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
# process.options.numberOfThreads=cms.untracked.uint32(2)
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
