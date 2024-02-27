#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TLorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include <vector>
#include <iostream>

class DGLMuonVertexTableProducer : public edm::global::EDProducer<> {

  public:
    explicit DGLMuonVertexTableProducer(const edm::ParameterSet &iConfig)
      :
      dglMuonTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dglMuons"))),
      dsaMuonTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dsaMuons"))),
      patMuonTag_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("patMuons"))),
      bsTag_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
      generalTrackTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("generalTracks")))
      {
      produces<nanoaod::FlatTable>("PatDGLMuonVertex");
      produces<nanoaod::FlatTable>("DGLDSAMuonVertex");
      produces<nanoaod::FlatTable>("DGLMuonVertex");
    }

    ~DGLMuonVertexTableProducer() override {}

    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("patMuons")->setComment("input pat muon collection");
      desc.add<edm::InputTag>("dsaMuons")->setComment("input displaced standalone muon collection");
      desc.add<edm::InputTag>("dglMuons")->setComment("input displaced global muon collection");
      desc.add<edm::InputTag>("beamspot")->setComment("input beamspot collection");
      desc.add<edm::InputTag>("generalTracks")->setComment("input generalTracks collection");
      descriptions.add("dglMuonVertexTables", desc);
    }

  private:
    void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

    std::pair<float,float> getVxy(const reco::Vertex muonVertex) const;

    template <typename MuonType1 = pat::Muon, typename MuonType2>
    float getDisplacedTrackerIsolation(const std::vector<reco::Track>& generalTracks, const MuonType1& muon_1,
                                      const reco::Vertex muonVertex, const reco::BeamSpot& beamspot, 
                                      const MuonType2* muon_2 = nullptr, float maxDR = 0.3, float minDR = 0.01,
                                      float maxDz = 0.5, float maxDxy = 0.2) const;

    const edm::EDGetTokenT<std::vector<reco::Track>> dglMuonTag_;
    const edm::EDGetTokenT<std::vector<reco::Track>> dsaMuonTag_;
    const edm::EDGetTokenT<std::vector<pat::Muon>> patMuonTag_;
    const edm::EDGetTokenT<reco::BeamSpot> bsTag_;
    const edm::EDGetTokenT<std::vector<reco::Track>> generalTrackTag_;

};

void DGLMuonVertexTableProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const 
{
  
  edm::Handle<std::vector<reco::Track>> dglMuonHandle;
  iEvent.getByToken(dglMuonTag_, dglMuonHandle);
  edm::Handle<std::vector<reco::Track>> dsaMuonHandle;
  iEvent.getByToken(dsaMuonTag_, dsaMuonHandle);
  edm::Handle<std::vector<pat::Muon>> patMuonHandle;
  iEvent.getByToken(patMuonTag_, patMuonHandle);
  edm::Handle<reco::BeamSpot> beamspots;
  iEvent.getByToken(bsTag_, beamspots);
  edm::Handle<std::vector<reco::Track>> generalTracks;
  iEvent.getByToken(generalTrackTag_, generalTracks);

  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  KalmanVertexFitter vertexFitter(true);

  std::vector<float> patdglVxy,patdglVx,patdglVy,patdglVz,patdglChi2,patdglVxySigma,patdglDR,patdglIdx1,patdglIdx2,patdglIsDSAMuon1,patdglIsDSAMuon2,patdglIsDGLMuon1,patdglIsDGLMuon2;
  std::vector<float> dgldsaVxy,dgldsaVx,dgldsaVy,dgldsaVz,dgldsaChi2,dgldsaVxySigma,dgldsaDR,dgldsaIdx1,dgldsaIdx2,dgldsaIsDSAMuon1,dgldsaIsDSAMuon2,dgldsaIsDGLMuon1,dgldsaIsDGLMuon2;
  std::vector<float> dgldglVxy,dgldglVx,dgldglVy,dgldglVz,dgldglChi2,dgldglVxySigma,dgldglDR,dgldglIdx1,dgldglIdx2,dgldglIsDSAMuon1,dgldglIsDSAMuon2,dgldglIsDGLMuon1,dgldglIsDGLMuon2;
  std::vector<float> patdglDisplacedTrackIso03Dimuon1,patdglDisplacedTrackIso03Dimuon2,patdglDisplacedTrackIso04Dimuon1,patdglDisplacedTrackIso04Dimuon2;
  std::vector<float> dgldsaDisplacedTrackIso03Dimuon1,dgldsaDisplacedTrackIso03Dimuon2,dgldsaDisplacedTrackIso04Dimuon1,dgldsaDisplacedTrackIso04Dimuon2;
  std::vector<float> dgldglDisplacedTrackIso03Dimuon1,dgldglDisplacedTrackIso03Dimuon2,dgldglDisplacedTrackIso04Dimuon1,dgldglDisplacedTrackIso04Dimuon2;
  std::vector<float> patdglDisplacedTrackIso03Muon1,patdglDisplacedTrackIso03Muon2,patdglDisplacedTrackIso04Muon1,patdglDisplacedTrackIso04Muon2;
  std::vector<float> dgldsaDisplacedTrackIso03Muon1,dgldsaDisplacedTrackIso03Muon2,dgldsaDisplacedTrackIso04Muon1,dgldsaDisplacedTrackIso04Muon2;
  std::vector<float> dgldglDisplacedTrackIso03Muon1,dgldglDisplacedTrackIso03Muon2,dgldglDisplacedTrackIso04Muon1,dgldglDisplacedTrackIso04Muon2;

  // pat muons
  for(size_t i = 0; i < patMuonHandle->size(); i++){

    const pat::Muon& muon_i = patMuonHandle->at(i);
    const auto& muonTrack_i = muon_i.bestTrack();
    reco::TransientTrack muonTransientTrack_i = builder->build(muonTrack_i);

    // pat-dgl muon vertex
    for(size_t j = 0; j < dglMuonHandle->size(); j++){
      const auto& muonTrack_j = dglMuonHandle->at(j);
      reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

      std::vector<reco::TransientTrack> muonTransientTracks{};
      muonTransientTracks.push_back(muonTransientTrack_i);
      muonTransientTracks.push_back(muonTransientTrack_j);

      TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
      reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

      std::pair<float,float> vxy = getVxy(muonVertex);
      patdglVxy.push_back(vxy.first);
      patdglVxySigma.push_back(vxy.second);
      patdglChi2.push_back(muonVertex.normalizedChi2());
      patdglVx.push_back(muonVertex.x());
      patdglVy.push_back(muonVertex.y());
      patdglVz.push_back(muonVertex.z());
      patdglDR.push_back(reco::deltaR(muon_i.eta(), muon_i.phi(), muonTrack_j.eta(), muonTrack_j.phi()));
      patdglIdx1.push_back(i);
      patdglIdx2.push_back(j);
      patdglIsDSAMuon1.push_back(0);
      patdglIsDSAMuon2.push_back(0);
      patdglIsDGLMuon1.push_back(0);
      patdglIsDGLMuon2.push_back(1);

      patdglDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.3));
      patdglDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.4));
      patdglDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muon_i, 0.3));
      patdglDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muon_i, 0.4));
      patdglDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), nullptr, 0.3));
      patdglDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), nullptr, 0.4));
      patdglDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.3));
      patdglDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.4));
    }
  }

  // dgl muons
  for(size_t i = 0; i < dglMuonHandle->size(); i++){

    const auto& muonTrack_i = dglMuonHandle->at(i);
    reco::TransientTrack muonTransientTrack_i = builder->build(muonTrack_i);

    //dgl-dgl muon vertex
    for(size_t j = i + 1; j < dglMuonHandle->size(); j++){
      
      const auto& muonTrack_j = dglMuonHandle->at(j);
      reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

      std::vector<reco::TransientTrack> muonTransientTracks{};
      muonTransientTracks.push_back(muonTransientTrack_i);
      muonTransientTracks.push_back(muonTransientTrack_j);

      TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
      reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

      std::pair<float,float> vxy = getVxy(muonVertex);
      dgldglVxy.push_back(vxy.first);
      dgldglVxySigma.push_back(vxy.second);
      dgldglChi2.push_back(muonVertex.normalizedChi2());
      dgldglVx.push_back(muonVertex.x());
      dgldglVy.push_back(muonVertex.y());
      dgldglVz.push_back(muonVertex.z());
      dgldglDR.push_back(reco::deltaR(muonTrack_i.eta(), muonTrack_i.phi(), muonTrack_j.eta(), muonTrack_j.phi()));
      dgldglIdx1.push_back(i);
      dgldglIdx2.push_back(j);
      dgldglIsDSAMuon1.push_back(0);
      dgldglIsDSAMuon2.push_back(0);
      dgldglIsDGLMuon1.push_back(1);
      dgldglIsDGLMuon2.push_back(1);

      dgldglDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.3));
      dgldglDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.4));
      dgldglDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muonTrack_i, 0.3));
      dgldglDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muonTrack_i, 0.4));
      dgldglDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), nullptr, 0.3));
      dgldglDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), nullptr, 0.4));
      dgldglDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.3));
      dgldglDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.4));
    }
    
    //dgl-dsa muon vertex
    for(size_t j = 0; j < dsaMuonHandle->size(); j++){
      
      const auto& muonTrack_j = dsaMuonHandle->at(j);
      reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

      std::vector<reco::TransientTrack> muonTransientTracks{};
      muonTransientTracks.push_back(muonTransientTrack_i);
      muonTransientTracks.push_back(muonTransientTrack_j);

      TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
      reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

      std::pair<float,float> vxy = getVxy(muonVertex);
      dgldsaVxy.push_back(vxy.first);
      dgldsaVxySigma.push_back(vxy.second);
      dgldsaChi2.push_back(muonVertex.normalizedChi2());
      dgldsaVx.push_back(muonVertex.x());
      dgldsaVy.push_back(muonVertex.y());
      dgldsaVz.push_back(muonVertex.z());
      dgldsaDR.push_back(reco::deltaR(muonTrack_i.eta(), muonTrack_i.phi(), muonTrack_j.eta(), muonTrack_j.phi()));
      dgldsaIdx1.push_back(i);
      dgldsaIdx2.push_back(j);
      dgldsaIsDSAMuon1.push_back(0);
      dgldsaIsDSAMuon2.push_back(1);
      dgldsaIsDGLMuon1.push_back(1);
      dgldsaIsDGLMuon2.push_back(0);

      dgldsaDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.3));
      dgldsaDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.4));
      dgldsaDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muonTrack_i, 0.3));
      dgldsaDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muonTrack_i, 0.4));
      dgldsaDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), nullptr, 0.3));
      dgldsaDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), nullptr, 0.4));
      dgldsaDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.3));
      dgldsaDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.4));
    }
  }
  
  auto patdglVertexTab = std::make_unique<nanoaod::FlatTable>(patMuonHandle->size() * dglMuonHandle->size(), "PatDGLMuonVertex", false, false);
  auto dgldsaVertexTab = std::make_unique<nanoaod::FlatTable>(dglMuonHandle->size() * dsaMuonHandle->size(), "DGLDSAMuonVertex", false, false);
  auto dglVertexTab = std::make_unique<nanoaod::FlatTable>(dglMuonHandle->size() * (dglMuonHandle->size() - 1) / 2, "DGLMuonVertex", false, false);

  patdglVertexTab->addColumn<float>("vxy", patdglVxy, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("vxySigma", patdglVxySigma, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("chi2", patdglChi2, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("vx", patdglVx, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("vy", patdglVy, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("vz", patdglVz, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("dR", patdglDR, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("idx1", patdglIdx1, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("idx2", patdglIdx2, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("isDSAMuon1", patdglIsDSAMuon1, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("isDSAMuon2", patdglIsDSAMuon2, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("isDGLMuon1", patdglIsDGLMuon1, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("isDGLMuon2", patdglIsDGLMuon2, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", patdglDisplacedTrackIso03Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", patdglDisplacedTrackIso04Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", patdglDisplacedTrackIso03Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", patdglDisplacedTrackIso04Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("displacedTrackIso03Muon1", patdglDisplacedTrackIso03Muon1, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("displacedTrackIso04Muon1", patdglDisplacedTrackIso04Muon1, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("displacedTrackIso03Muon2", patdglDisplacedTrackIso03Muon2, "",  nanoaod::FlatTable::FloatColumn);
  patdglVertexTab->addColumn<float>("displacedTrackIso04Muon2", patdglDisplacedTrackIso04Muon2, "",  nanoaod::FlatTable::FloatColumn);

  dgldsaVertexTab->addColumn<float>("vxy", dgldsaVxy, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("vxySigma", dgldsaVxySigma, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("chi2", dgldsaChi2, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("vx", dgldsaVx, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("vy", dgldsaVy, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("vz", dgldsaVz, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("dR", dgldsaDR, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("idx1", dgldsaIdx1, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("idx2", dgldsaIdx2, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("isDSAMuon1", dgldsaIsDSAMuon1, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("isDSAMuon2", dgldsaIsDSAMuon2, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("isDGLMuon1", dgldsaIsDGLMuon1, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("isDGLMuon2", dgldsaIsDGLMuon2, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", dgldsaDisplacedTrackIso03Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", dgldsaDisplacedTrackIso04Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", dgldsaDisplacedTrackIso03Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", dgldsaDisplacedTrackIso04Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("displacedTrackIso03Muon1", dgldsaDisplacedTrackIso03Muon1, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("displacedTrackIso04Muon1", dgldsaDisplacedTrackIso04Muon1, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("displacedTrackIso03Muon2", dgldsaDisplacedTrackIso03Muon2, "",  nanoaod::FlatTable::FloatColumn);
  dgldsaVertexTab->addColumn<float>("displacedTrackIso04Muon2", dgldsaDisplacedTrackIso04Muon2, "",  nanoaod::FlatTable::FloatColumn);

  dglVertexTab->addColumn<float>("vxy", dgldglVxy, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("vxySigma", dgldglVxySigma, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("chi2", dgldglChi2, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("vx", dgldglVx, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("vy", dgldglVy, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("vz", dgldglVz, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("dR", dgldglDR, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("idx1", dgldglIdx1, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("idx2", dgldglIdx2, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("isDSAMuon1", dgldglIsDSAMuon1, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("isDSAMuon2", dgldglIsDSAMuon2, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("isDGLMuon1", dgldglIsDGLMuon1, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("isDGLMuon2", dgldglIsDGLMuon2, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", dgldglDisplacedTrackIso03Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", dgldglDisplacedTrackIso04Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", dgldglDisplacedTrackIso03Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", dgldglDisplacedTrackIso04Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("displacedTrackIso03Muon1", dgldglDisplacedTrackIso03Muon1, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("displacedTrackIso04Muon1", dgldglDisplacedTrackIso04Muon1, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("displacedTrackIso03Muon2", dgldglDisplacedTrackIso03Muon2, "",  nanoaod::FlatTable::FloatColumn);
  dglVertexTab->addColumn<float>("displacedTrackIso04Muon2", dgldglDisplacedTrackIso04Muon2, "",  nanoaod::FlatTable::FloatColumn);
  
  iEvent.put(std::move(patdglVertexTab), "PatDGLMuonVertex");
  iEvent.put(std::move(dgldsaVertexTab), "DGLDSAMuonVertex");
  iEvent.put(std::move(dglVertexTab), "DGLMuonVertex");
}

std::pair<float,float> DGLMuonVertexTableProducer::getVxy(const reco::Vertex muonVertex) const {
  float vxy = sqrt(muonVertex.x()*muonVertex.x() + muonVertex.y()*muonVertex.y());
  float vxySigma = (1/vxy)*sqrt(muonVertex.x()*muonVertex.x()*muonVertex.xError()*muonVertex.xError() + muonVertex.y()*muonVertex.y()*muonVertex.yError()*muonVertex.yError());
  return std::make_pair(vxy,vxySigma);
}

template <typename MuonType1, typename MuonType2>
float DGLMuonVertexTableProducer::getDisplacedTrackerIsolation(const std::vector<reco::Track>& generalTracks, 
                                                            const MuonType1& muon_1, const reco::Vertex muonVertex, 
                                                            const reco::BeamSpot& beamspot, const MuonType2* muon_2, 
                                                            float maxDR, float minDR, float maxDz, float maxDxy) const 
{
  float trackPtSum = 0;

  int nGeneralTracks = generalTracks.size();
  float muonTrack2_pt = 0;
  float muonTrack2_minDR = 9999;

  for (int i = 0; i < nGeneralTracks; i++) {
    const reco::Track & generalTrack = (generalTracks)[i];

    // Muon POG Tracker Isolation recommendation
    float dR = deltaR(muon_1.eta(), muon_1.phi(), generalTrack.eta(), generalTrack.phi());
    if (dR > maxDR) continue;
    if (abs(generalTrack.vz() - muonVertex.z()) > maxDz) continue;
    if (generalTrack.dxy(beamspot) > maxDxy) continue;
    if (dR < minDR) continue;

    // Determine if track belongs to other muon and get pt of the track
    // Only if muon is given as input
    if(muon_2 != nullptr) {
      float dR_2 = deltaR(muon_2->eta(), muon_2->phi(), generalTrack.eta(), generalTrack.phi());
      if (dR_2 < minDR && dR_2 < muonTrack2_minDR) {
        muonTrack2_pt = generalTrack.pt();
        muonTrack2_minDR = dR_2;
      }
    }

    trackPtSum += generalTrack.pt();
  }

  // Remove pt of track that belongs to other muon
  trackPtSum -= muonTrack2_pt;

  float ptRatio = trackPtSum / muon_1.pt();
  return ptRatio;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DGLMuonVertexTableProducer);
