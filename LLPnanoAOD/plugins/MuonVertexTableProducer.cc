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

class MuonVertexTableProducer : public edm::global::EDProducer<> {

  public:
    explicit MuonVertexTableProducer(const edm::ParameterSet &iConfig)
      :
      dsaMuonTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dsaMuons"))),
      patMuonTag_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("patMuons"))),
      bsTag_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
      generalTrackTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("generalTracks")))
      {
      produces<nanoaod::FlatTable>("PatMuonVertex");
      produces<nanoaod::FlatTable>("PatDSAMuonVertex");
      produces<nanoaod::FlatTable>("DSAMuonVertex");
    }

    ~MuonVertexTableProducer() override {}

    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("patMuons")->setComment("input pat muon collection");
      desc.add<edm::InputTag>("dsaMuons")->setComment("input displaced standalone muon collection");
      desc.add<edm::InputTag>("beamspot")->setComment("input beamspot collection");
      desc.add<edm::InputTag>("generalTracks")->setComment("input generalTracks collection");
      descriptions.add("muonVertexTables", desc);
    }

  private:
    void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

    std::pair<float,float> getVxy(const reco::Vertex muonVertex) const;

    template <typename MuonType1 = pat::Muon, typename MuonType2>
    float getDisplacedTrackerIsolation(const std::vector<reco::Track>& generalTracks, const MuonType1& muon_1,
                                      const reco::Vertex muonVertex, const reco::BeamSpot& beamspot, 
                                      const MuonType2* muon_2 = nullptr, float maxDR = 0.3, float minDR = 0.01,
                                      float maxDz = 0.5, float maxDxy = 0.2) const;

    const edm::EDGetTokenT<std::vector<reco::Track>> dsaMuonTag_;
    const edm::EDGetTokenT<std::vector<pat::Muon>> patMuonTag_;
    const edm::EDGetTokenT<reco::BeamSpot> bsTag_;
    const edm::EDGetTokenT<std::vector<reco::Track>> generalTrackTag_;

};

void MuonVertexTableProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const 
{
  

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

  std::vector<float> ppVxy,ppVx,ppVy,ppVz,ppChi2,ppVxySigma,ppDR,ppIdx1,ppIdx2,ppIsDSAMuon1,ppIsDSAMuon2;
  std::vector<float> pdVxy,pdVx,pdVy,pdVz,pdChi2,pdVxySigma,pdDR,pdIdx1,pdIdx2,pdIsDSAMuon1,pdIsDSAMuon2;
  std::vector<float> ddVxy,ddVx,ddVy,ddVz,ddChi2,ddVxySigma,ddDR,ddIdx1,ddIdx2,ddIsDSAMuon1,ddIsDSAMuon2;
  std::vector<float> ppDisplacedTrackIso03Dimuon1,ppDisplacedTrackIso03Dimuon2,ppDisplacedTrackIso04Dimuon1,ppDisplacedTrackIso04Dimuon2;
  std::vector<float> pdDisplacedTrackIso03Dimuon1,pdDisplacedTrackIso03Dimuon2,pdDisplacedTrackIso04Dimuon1,pdDisplacedTrackIso04Dimuon2;
  std::vector<float> ddDisplacedTrackIso03Dimuon1,ddDisplacedTrackIso03Dimuon2,ddDisplacedTrackIso04Dimuon1,ddDisplacedTrackIso04Dimuon2;
  std::vector<float> ppDisplacedTrackIso03Muon1,ppDisplacedTrackIso03Muon2,ppDisplacedTrackIso04Muon1,ppDisplacedTrackIso04Muon2;
  std::vector<float> pdDisplacedTrackIso03Muon1,pdDisplacedTrackIso03Muon2,pdDisplacedTrackIso04Muon1,pdDisplacedTrackIso04Muon2;
  std::vector<float> ddDisplacedTrackIso03Muon1,ddDisplacedTrackIso03Muon2,ddDisplacedTrackIso04Muon1,ddDisplacedTrackIso04Muon2;

  // pat muons
  for(size_t i = 0; i < patMuonHandle->size(); i++){

    const pat::Muon& muon_i = patMuonHandle->at(i);
    const auto& muonTrack_i = muon_i.bestTrack();
    reco::TransientTrack muonTransientTrack_i = builder->build(muonTrack_i);

    // pat-pat muon vertex
    for(size_t j = i+1; j < patMuonHandle->size(); j++){
    
      const pat::Muon& muon_j = patMuonHandle->at(j);
      const auto& muonTrack_j = muon_j.bestTrack();
      reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

      std::vector<reco::TransientTrack> muonTransientTracks{};
      muonTransientTracks.push_back(muonTransientTrack_i);
      muonTransientTracks.push_back(muonTransientTrack_j);

      TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
      reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

      std::pair<float,float> vxy = getVxy(muonVertex);
      ppVxy.push_back(vxy.first);
      ppVxySigma.push_back(vxy.second);
      ppChi2.push_back(muonVertex.normalizedChi2());
      ppVx.push_back(muonVertex.x());
      ppVy.push_back(muonVertex.y());
      ppVz.push_back(muonVertex.z());
      ppDR.push_back(reco::deltaR(muon_i, muon_j));
      ppIdx1.push_back(i);
      ppIdx2.push_back(j);
      ppIsDSAMuon1.push_back(0);
      ppIsDSAMuon2.push_back(0);

      ppDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), &muon_j, 0.3));
      ppDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), &muon_j, 0.4));
      ppDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_j, muonVertex,*beamspots.product(), &muon_i, 0.3));
      ppDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_j, muonVertex,*beamspots.product(), &muon_i, 0.4));
      ppDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), nullptr, 0.3));
      ppDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), nullptr, 0.4));
      ppDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_j, muonVertex,*beamspots.product(), nullptr, 0.3));
      ppDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_j, muonVertex,*beamspots.product(), nullptr, 0.4));
    }

    // pat-dsa muon vertex
    for(size_t j = 0; j < dsaMuonHandle->size(); j++){
      const auto& muonTrack_j = dsaMuonHandle->at(j);
      reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

      std::vector<reco::TransientTrack> muonTransientTracks{};
      muonTransientTracks.push_back(muonTransientTrack_i);
      muonTransientTracks.push_back(muonTransientTrack_j);

      TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
      reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

      std::pair<float,float> vxy = getVxy(muonVertex);
      pdVxy.push_back(vxy.first);
      pdVxySigma.push_back(vxy.second);
      pdChi2.push_back(muonVertex.normalizedChi2());
      pdVx.push_back(muonVertex.x());
      pdVy.push_back(muonVertex.y());
      pdVz.push_back(muonVertex.z());
      pdDR.push_back(reco::deltaR(muon_i.eta(), muon_i.phi(), muonTrack_j.eta(), muonTrack_j.phi()));
      pdIdx1.push_back(i);
      pdIdx2.push_back(j);
      pdIsDSAMuon1.push_back(0);
      pdIsDSAMuon2.push_back(1);

      pdDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.3));
      pdDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.4));
      pdDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muon_i, 0.3));
      pdDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muon_i, 0.4));
      pdDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), nullptr, 0.3));
      pdDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamspots.product(), nullptr, 0.4));
      pdDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.3));
      pdDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.4));
    }
  }

  // dsa muons
  for(size_t i = 0; i < dsaMuonHandle->size(); i++){

    const auto& muonTrack_i = dsaMuonHandle->at(i);
    reco::TransientTrack muonTransientTrack_i = builder->build(muonTrack_i);

    //dsa-dsa muon vertex
    for(size_t j = i + 1; j < dsaMuonHandle->size(); j++){
      
      const auto& muonTrack_j = dsaMuonHandle->at(j);
      reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

      std::vector<reco::TransientTrack> muonTransientTracks{};
      muonTransientTracks.push_back(muonTransientTrack_i);
      muonTransientTracks.push_back(muonTransientTrack_j);

      TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
      reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

      std::pair<float,float> vxy = getVxy(muonVertex);
      ddVxy.push_back(vxy.first);
      ddVxySigma.push_back(vxy.second);
      ddChi2.push_back(muonVertex.normalizedChi2());
      ddVx.push_back(muonVertex.x());
      ddVy.push_back(muonVertex.y());
      ddVz.push_back(muonVertex.z());
      ddDR.push_back(reco::deltaR(muonTrack_i.eta(), muonTrack_i.phi(), muonTrack_j.eta(), muonTrack_j.phi()));
      ddIdx1.push_back(i);
      ddIdx2.push_back(j);
      ddIsDSAMuon1.push_back(1);
      ddIsDSAMuon2.push_back(1);

      ddDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.3));
      ddDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), &muonTrack_j, 0.4));
      ddDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muonTrack_i, 0.3));
      ddDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), &muonTrack_i, 0.4));
      ddDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), nullptr, 0.3));
      ddDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamspots.product(), nullptr, 0.4));
      ddDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.3));
      ddDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamspots.product(), nullptr, 0.4));
    }
  }
  
  auto patVertexTab = std::make_unique<nanoaod::FlatTable>(patMuonHandle->size() * (patMuonHandle->size() - 1) / 2, "PatMuonVertex", false, false);
  auto patdsaVertexTab = std::make_unique<nanoaod::FlatTable>(patMuonHandle->size() * dsaMuonHandle->size(), "PatDSAMuonVertex", false, false);
  auto dsaVertexTab = std::make_unique<nanoaod::FlatTable>(dsaMuonHandle->size() * (dsaMuonHandle->size() - 1) / 2, "DSAMuonVertex", false, false);

  patVertexTab->addColumn<float>("vxy", ppVxy, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("vxySigma", ppVxySigma, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("chi2", ppChi2, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("vx", ppVx, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("vy", ppVy, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("vz", ppVz, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("dR", ppDR, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("idx1", ppIdx1, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("idx2", ppIdx2, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("isDSAMuon1", ppIsDSAMuon1, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("isDSAMuon2", ppIsDSAMuon2, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", ppDisplacedTrackIso03Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", ppDisplacedTrackIso04Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", ppDisplacedTrackIso03Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", ppDisplacedTrackIso04Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("displacedTrackIso03Muon1", ppDisplacedTrackIso03Muon1, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("displacedTrackIso04Muon1", ppDisplacedTrackIso04Muon1, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("displacedTrackIso03Muon2", ppDisplacedTrackIso03Muon2, "",  nanoaod::FlatTable::FloatColumn);
  patVertexTab->addColumn<float>("displacedTrackIso04Muon2", ppDisplacedTrackIso04Muon2, "",  nanoaod::FlatTable::FloatColumn);

  patdsaVertexTab->addColumn<float>("vxy", pdVxy, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("vxySigma", pdVxySigma, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("chi2", pdChi2, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("vx", pdVx, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("vy", pdVy, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("vz", pdVz, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("dR", pdDR, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("idx1", pdIdx1, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("idx2", pdIdx2, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("isDSAMuon1", pdIsDSAMuon1, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("isDSAMuon2", pdIsDSAMuon2, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", pdDisplacedTrackIso03Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", pdDisplacedTrackIso04Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", pdDisplacedTrackIso03Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", pdDisplacedTrackIso04Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("displacedTrackIso03Muon1", pdDisplacedTrackIso03Muon1, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("displacedTrackIso04Muon1", pdDisplacedTrackIso04Muon1, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("displacedTrackIso03Muon2", pdDisplacedTrackIso03Muon2, "",  nanoaod::FlatTable::FloatColumn);
  patdsaVertexTab->addColumn<float>("displacedTrackIso04Muon2", pdDisplacedTrackIso04Muon2, "",  nanoaod::FlatTable::FloatColumn);

  dsaVertexTab->addColumn<float>("vxy", ddVxy, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("vxySigma", ddVxySigma, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("chi2", ddChi2, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("vx", ddVx, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("vy", ddVy, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("vz", ddVz, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("dR", ddDR, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("idx1", ddIdx1, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("idx2", ddIdx2, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("isDSAMuon1", ddIsDSAMuon1, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("isDSAMuon2", ddIsDSAMuon2, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", ddDisplacedTrackIso03Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", ddDisplacedTrackIso04Dimuon1, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", ddDisplacedTrackIso03Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", ddDisplacedTrackIso04Dimuon2, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("displacedTrackIso03Muon1", ddDisplacedTrackIso03Muon1, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("displacedTrackIso04Muon1", ddDisplacedTrackIso04Muon1, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("displacedTrackIso03Muon2", ddDisplacedTrackIso03Muon2, "",  nanoaod::FlatTable::FloatColumn);
  dsaVertexTab->addColumn<float>("displacedTrackIso04Muon2", ddDisplacedTrackIso04Muon2, "",  nanoaod::FlatTable::FloatColumn);

  
  iEvent.put(std::move(patVertexTab), "PatMuonVertex");
  iEvent.put(std::move(patdsaVertexTab), "PatDSAMuonVertex");
  iEvent.put(std::move(dsaVertexTab), "DSAMuonVertex");
}

std::pair<float,float> MuonVertexTableProducer::getVxy(const reco::Vertex muonVertex) const {
  float vxy = sqrt(muonVertex.x()*muonVertex.x() + muonVertex.y()*muonVertex.y());
  float vxySigma = (1/vxy)*sqrt(muonVertex.x()*muonVertex.x()*muonVertex.xError()*muonVertex.xError() + muonVertex.y()*muonVertex.y()*muonVertex.yError()*muonVertex.yError());
  return std::make_pair(vxy,vxySigma);
}

template <typename MuonType1, typename MuonType2>
float MuonVertexTableProducer::getDisplacedTrackerIsolation(const std::vector<reco::Track>& generalTracks, 
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
DEFINE_FWK_MODULE(MuonVertexTableProducer);
