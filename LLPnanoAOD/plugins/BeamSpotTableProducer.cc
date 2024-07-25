#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "TLorentzVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


#include <vector>
#include <iostream>

class BeamSpotTableProducer : public edm::stream::EDProducer<> {
protected:
  const edm::InputTag beamSpotLabel;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotTag_;
  
public:
  BeamSpotTableProducer(edm::ParameterSet const& params)
    :
    beamSpotLabel(params.getParameter<edm::InputTag>("beamSpot")),
    beamSpotTag_(consumes<reco::BeamSpot>(beamSpotLabel))
    {
    produces<nanoaod::FlatTable>("BS");
    
  }

  ~BeamSpotTableProducer() override {}

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotTag_, beamSpotHandle);
    auto beamSpotTab = std::make_unique<nanoaod::FlatTable>(1, "BS", false, false);

    const auto& beamSpot = *beamSpotHandle;
    reco::Vertex beamSpotVertex(beamSpot.position(), beamSpot.covariance3D());

    std::vector<float> x, y, z, ndof, chi2, ntracks;
    x.push_back(beamSpot.position().x());
    y.push_back(beamSpot.position().y());
    z.push_back(beamSpot.position().z());
    ndof.push_back(beamSpotVertex.ndof());
    chi2.push_back(beamSpotVertex.normalizedChi2());
    ntracks.push_back(beamSpotVertex.nTracks());

    beamSpotTab->addColumn<float>("x", x, "");
    beamSpotTab->addColumn<float>("y", y, "");
    beamSpotTab->addColumn<float>("z", z, "");
    beamSpotTab->addColumn<float>("ndof", ndof, "");
    beamSpotTab->addColumn<float>("chi2", chi2, "");
    beamSpotTab->addColumn<uint8_t>("ntracks", ntracks, "");

    iEvent.put(std::move(beamSpotTab), "BS");
  }

};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BeamSpotTableProducer);
