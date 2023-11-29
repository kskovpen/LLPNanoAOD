import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.electrons_cff import *
from PhysicsTools.NanoAOD.lowPtElectrons_cff import *

def LLPNanoAOD_customize(process):

  myMuonTable = muonTable.clone()
  myMuonTable.variables.vx = Var("bestTrack().vx()", float)
  myMuonTable.variables.vy = Var("bestTrack().vy()", float)
  myMuonTable.variables.vz = Var("bestTrack().vz()", float)

  myMuonTable.variables.trkNumPlanes = Var("bestTrack().hitPattern().muonStationsWithValidHits()", float)
  myMuonTable.variables.trkNumHits = Var("bestTrack().hitPattern().numberOfValidMuonHits()", float)
  myMuonTable.variables.trkNumDTHits = Var("bestTrack().hitPattern().numberOfValidMuonDTHits()", float)
  myMuonTable.variables.trkNumCSCHits = Var("bestTrack().hitPattern().numberOfValidMuonCSCHits()", float)
  myMuonTable.variables.normChi2 = Var("bestTrack().normalizedChi2()", float)

  process.globalReplace("muonTable", myMuonTable)

  myGenParticleTable = genParticleTable.clone()
  myGenParticleTable.variables.vx = Var("vertex.X", float)
  myGenParticleTable.variables.vy = Var("vertex.Y", float)
  myGenParticleTable.variables.vz = Var("vertex.Z", float)
  myGenParticleTable.variables.Rho = Var("vertex.Rho", float)
  myGenParticleTable.variables.R = Var("vertex.R", float)

  process.globalReplace("genParticleTable", myGenParticleTable)

  return process