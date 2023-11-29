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
options.parseArguments()

nevents = options.nEvents
if nevents == 0:
    nevents=-1
    
print('Running run_LLPNanoAOD.py')
# print('-- Input MiniAOD files: ')
# for file in options.inputFiles:
#     print(file)
print('-- Output NanoAOD file: '+options.outputFile)
print('-- Running on '+str(nevents)+' number of events')

process = cms.Process('NANO',Run2_2018,run2_nanoAOD_106Xv2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    duplicateCheckMode = cms.untracked.string('checkEachRealDataFile')
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

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:{0}'.format(options.outputFile)),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')

# DisplacedStandAloneMuon table
process.dSAMuonsTable = cms.EDProducer("DSAMuonTableProducer",
    displacedMuons=cms.InputTag("displacedStandAloneMuons"),
    primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot")
)
# BeamSpot table
process.beamSpotTable = cms.EDProducer("BeamSpotTableProducer",
    beamSpot = cms.InputTag("offlineBeamSpot")
)
# Vertex between two (default) muons
process.muonVertexTable = cms.EDProducer("MuonVertexTableProducer",
    displacedMuons=cms.InputTag("displacedStandAloneMuons"),
    muons=cms.InputTag("linkedObjects","muons"),
    muonCombination = cms.string("muon")
)
# Vertex between one (default) muon and one displaced standalone muon
process.muonCombVertexTable = cms.EDProducer("MuonVertexTableProducer",
    displacedMuons=cms.InputTag("displacedStandAloneMuons"),
    muons=cms.InputTag("linkedObjects","muons"),
    muonCombination = cms.string("comb")
)
# Vertex between two displaced standalone muon
process.dsamuonVertexTable = cms.EDProducer("MuonVertexTableProducer",
    displacedMuons=cms.InputTag("displacedStandAloneMuons"),
    muons=cms.InputTag("linkedObjects","muons"),
    muonCombination = cms.string("dsa")
)

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC
                                +process.dSAMuonsTable
                                +process.beamSpotTable
                                +process.muonVertexTable
                                +process.muonCombVertexTable
                                +process.dsamuonVertexTable
                                )

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(2)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 
#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeMC(process)

from PhysicsTools.NanoAOD.llpnano_cff import LLPNanoAOD_customize
process = LLPNanoAOD_customize(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
