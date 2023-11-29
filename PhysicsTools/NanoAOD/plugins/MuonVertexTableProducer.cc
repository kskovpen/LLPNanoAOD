#include "FWCore/Framework/interface/stream/EDProducer.h"
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

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include <vector>
#include <iostream>

class MuonVertexTableProducer : public edm::stream::EDProducer<> {
protected:
  const edm::InputTag displacedMuonLabel;
  const edm::EDGetTokenT<std::vector<reco::Track>> displacedMuonTag_;
  const edm::InputTag muonLabel;
  const edm::EDGetTokenT<std::vector<pat::Muon>> muonTag_;
  const std::string type; 
  
public:
  MuonVertexTableProducer(edm::ParameterSet const& params)
    :
    displacedMuonLabel(params.getParameter<edm::InputTag>("displacedMuons")),
    displacedMuonTag_(consumes<std::vector<reco::Track>>(displacedMuonLabel)),
    muonLabel(params.getParameter<edm::InputTag>("muons")),
    muonTag_(consumes<std::vector<pat::Muon>>(muonLabel)),
    type(params.getParameter<std::string>("muonCombination"))
    {
    std::string tabName;
    if(type == "muon"){
      tabName = "MuonVertex";
    } else if(type == "comb"){
      tabName = "MuonCombVertex";
    } else if(type == "dsa"){
      tabName = "DSAMuonVertex";
    }
      
    produces<nanoaod::FlatTable>(tabName);
  }

  ~MuonVertexTableProducer() override {}

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    

    edm::Handle<std::vector<reco::Track>> displacedMuonHandle;
    iEvent.getByToken(displacedMuonTag_, displacedMuonHandle);

    edm::Handle<std::vector<pat::Muon>> muonHandle;
    iEvent.getByToken(muonTag_, muonHandle);
    size_t muon_size = muonHandle->size();
    size_t dsamuon_size = displacedMuonHandle->size();

    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

    KalmanVertexFitter vertexFitter(true);

    std::vector<float> Vxy,Vz,Chi2,VxySigma,DR;
    float vxy_;

    size_t tabSize;
    std::string tabName;
    if(type == "muon"){
      tabSize = muon_size * (muon_size - 1) / 2;
      tabName = "MuonVertex";
    } else if(type == "comb"){
      tabSize = muon_size * dsamuon_size;
      tabName = "MuonCombVertex";
    } else if(type == "dsa"){
      tabSize = dsamuon_size * (dsamuon_size - 1) / 2;
      tabName = "DSAMuonVertex";
    } else {
      std::cout << "Error incorrect input type: " << type; 
      std::cout << ". Correct types are muon, comb or dsa." << std::endl;
      return;
    }

    auto vertexTab = std::make_unique<nanoaod::FlatTable>(tabSize, tabName, false, false);

    if(type == "muon"){

      for(size_t i = 0; i < muon_size; i++){

        const pat::Muon& muon_i = muonHandle->at(i);
        const auto& muonTrack_i = muon_i.bestTrack();
        reco::TransientTrack muonTransientTrack_i = builder->build(muonTrack_i);

        for(size_t j = i+1; j < muon_size; j++){
        
          const pat::Muon& muon_j = muonHandle->at(j);
          const auto& muonTrack_j = muon_j.bestTrack();
          reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

          std::vector<reco::TransientTrack> muonTransientTracks{};
          muonTransientTracks.push_back(muonTransientTrack_i);
          muonTransientTracks.push_back(muonTransientTrack_j);

          TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
          reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

          vxy_ = sqrt(muonVertex.x()*muonVertex.x() + muonVertex.y()*muonVertex.y());
          Vxy.push_back(vxy_);
          VxySigma.push_back((1/vxy_)*sqrt(muonVertex.x()*muonVertex.x()*muonVertex.xError()*muonVertex.xError() + muonVertex.y()*muonVertex.y()*muonVertex.yError()*muonVertex.yError()));
          Chi2.push_back(muonVertex.normalizedChi2());
          Vz.push_back(muonVertex.z());
          DR.push_back(reco::deltaR(muon_i, muon_j));
        }
      }
    } else if(type == "comb"){

      for(size_t i = 0; i < muon_size; i++){

        const pat::Muon& muon_i = muonHandle->at(i);
        const auto& muonTrack_i = muon_i.bestTrack();
        reco::TransientTrack muonTransientTrack_i = builder->build(muonTrack_i);

        for(size_t j = 0; j < dsamuon_size; j++){
        
          const auto& muonTrack_j = displacedMuonHandle->at(j);
          reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

          std::vector<reco::TransientTrack> muonTransientTracks{};
          muonTransientTracks.push_back(muonTransientTrack_i);
          muonTransientTracks.push_back(muonTransientTrack_j);

          TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
          reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

          vxy_ = sqrt(muonVertex.x()*muonVertex.x() + muonVertex.y()*muonVertex.y());
          Vxy.push_back(vxy_);
          VxySigma.push_back((1/vxy_)*sqrt(muonVertex.x()*muonVertex.x()*muonVertex.xError()*muonVertex.xError() + muonVertex.y()*muonVertex.y()*muonVertex.yError()*muonVertex.yError()));
          Chi2.push_back(muonVertex.normalizedChi2());
          Vz.push_back(muonVertex.z());
          DR.push_back(reco::deltaR(muon_i.eta(), muon_i.phi(), muonTrack_j.eta(), muonTrack_j.phi()));
        }
      }
    } else {
      for(size_t i = 0; i < dsamuon_size; i++){

        const auto& dsamuonTrack_i = displacedMuonHandle->at(i);
        reco::TransientTrack dsamuonTransientTrack_i = builder->build(dsamuonTrack_i);

        for(size_t j = i + 1; j < dsamuon_size; j++){
          
          const auto& dsamuonTrack_j = displacedMuonHandle->at(j);
          reco::TransientTrack dsamuonTransientTrack_j = builder->build(dsamuonTrack_j);

          std::vector<reco::TransientTrack> dsamuonTransientTracks{};
          dsamuonTransientTracks.push_back(dsamuonTransientTrack_i);
          dsamuonTransientTracks.push_back(dsamuonTransientTrack_j);

          TransientVertex transientDsamuonVertex = vertexFitter.vertex(dsamuonTransientTracks);
          reco::Vertex dsamuonVertex = reco::Vertex(transientDsamuonVertex);

          vxy_ = sqrt(dsamuonVertex.x()*dsamuonVertex.x() + dsamuonVertex.y()*dsamuonVertex.y());
          Vxy.push_back(vxy_);
          VxySigma.push_back((1/vxy_)*sqrt(dsamuonVertex.x()*dsamuonVertex.x()*dsamuonVertex.xError()*dsamuonVertex.xError() + dsamuonVertex.y()*dsamuonVertex.y()*dsamuonVertex.yError()*dsamuonVertex.yError()));
          Chi2.push_back(dsamuonVertex.normalizedChi2());
          Vz.push_back(dsamuonVertex.z());
          DR.push_back(reco::deltaR(dsamuonTrack_i.eta(), dsamuonTrack_i.phi(), dsamuonTrack_j.eta(), dsamuonTrack_j.phi()));

        }
      }
    }

    vertexTab->addColumn<float>("vxy", Vxy, "",  nanoaod::FlatTable::FloatColumn);
    vertexTab->addColumn<float>("vxySigma", VxySigma, "",  nanoaod::FlatTable::FloatColumn);
    vertexTab->addColumn<float>("chi2", Chi2, "",  nanoaod::FlatTable::FloatColumn);
    vertexTab->addColumn<float>("vz", Vz, "",  nanoaod::FlatTable::FloatColumn);
    vertexTab->addColumn<float>("dR", DR, "",  nanoaod::FlatTable::FloatColumn);    

    iEvent.put(std::move(vertexTab), tabName);
  }

};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonVertexTableProducer);

