import FWCore.ParameterSet.Config as cms

unpackedTracksAndVertices = cms.EDProducer('PATTrackAndVertexUnpacker',
                                           slimmedVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                           slimmedSecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
                                           additionalTracks= cms.InputTag("lostTracks"),
                                           packedCandidates = cms.InputTag("packedPFCandidates"))
    
displacedInclusiveVertexFinder  = cms.EDProducer("InclusiveVertexFinder",
                                                 beamSpot = cms.InputTag("offlineBeamSpot"),
                                                 primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                 tracks = cms.InputTag("unpackedTracksAndVertices"),
                                                 minHits = cms.uint32(6), #old 8 -> 0 AOD produciton has problems with nhits
                                                 maximumLongitudinalImpactParameter = cms.double(20), #old  .3 -> infty
                                                 minPt = cms.double(0.8), #old .8 -> 1
                                                 maxNTracks = cms.uint32(100), #old 30 -> 100
                                                 
                                                 clusterizer = cms.PSet(
                                                     seedMax3DIPSignificance = cms.double(9999.),
                                                     seedMax3DIPValue = cms.double(9999.),
                                                     seedMin3DIPSignificance = cms.double(1.2),
                                                     seedMin3DIPValue = cms.double(0.005),
                                                     clusterMaxDistance = cms.double(0.4), #500um #old .05 -> 1
                                                     clusterMaxSignificance = cms.double(4.5), #4.5 sigma  #old  4.5 ---> infty
                                                     distanceRatio = cms.double(20), # was cluster scale = 1 / density factor =0.05
                                                     clusterMinAngleCosine = cms.double(0.5), # only forward decays   #old accept backward decays (unboosted topologies) .5 -> -9999
                                                 ),
                                                     
                                                 vertexMinAngleCosine = cms.double(0.95), # scalar prod direction of tracks and flight dir  #old accept backward decays  .95 -> .6 me 0
                                                 vertexMinDLen2DSig = cms.double(2.5), #2.5 sigma
                                                 vertexMinDLenSig = cms.double(0.5), #0.5 sigma
                                                 fitterSigmacut =  cms.double(3),
                                                 fitterTini = cms.double(256),
                                                 fitterRatio = cms.double(0.25),
                                                 useDirectVertexFitter = cms.bool(True),
                                                 useVertexReco  = cms.bool(True),
                                                 vertexReco = cms.PSet(
                                                     finder = cms.string('avr'),
                                                     primcut = cms.double(1.0),
                                                     seccut = cms.double(3),
                                                     smoothing = cms.bool(True))
                                            )
                                            
displacedVertexMerger = cms.EDProducer("VertexMerger",
                                       secondaryVertices = cms.InputTag("displacedInclusiveVertexFinder"),
                                       maxFraction = cms.double(0.7), #old .7 -> .5
                                       minSignificance = cms.double(2))
                                           
displacedTrackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
                                                beamSpot = cms.InputTag("offlineBeamSpot"),
                                                primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                tracks = cms.InputTag("unpackedTracksAndVertices"),
                                                secondaryVertices = cms.InputTag("displacedVertexMerger"),
                                                dLenFraction = cms.double(0.333), #old .333 -> .2
                                                dRCut = cms.double(1), # old .4 -> 1   me 3
                                                distCut = cms.double(0.1), #old .04 -> .1 
                                                sigCut = cms.double(5), #old 5->10
                                                fitterSigmacut =  cms.double(3),
                                                fitterTini = cms.double(256),
                                                fitterRatio = cms.double(0.25),
                                                trackMinLayers = cms.int32(4), #old 4-> 0
                                                trackMinPt = cms.double(.4),
                                                trackMinPixels = cms.int32(0) #old 1 -> 0
)
    
displacedInclusiveSecondaryVertices = displacedVertexMerger.clone()
displacedInclusiveSecondaryVertices.secondaryVertices = cms.InputTag("displacedTrackVertexArbitrator")
displacedInclusiveSecondaryVertices.maxFraction = 0.2 #0.05 #old .2 -> .05
displacedInclusiveSecondaryVertices.minSignificance = 10
    
displacedInclusiveVertexing = cms.Sequence(unpackedTracksAndVertices * displacedInclusiveVertexFinder  * displacedVertexMerger * displacedTrackVertexArbitrator * displacedInclusiveSecondaryVertices)
