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

    float x, y, z, ndof, chi2;
    uint8_t ntracks;
    x = beamSpot.position().x();
    y = beamSpot.position().y();
    z = beamSpot.position().z();
    ndof = beamSpotVertex.ndof();
    chi2 = beamSpotVertex.normalizedChi2();
    ntracks = beamSpotVertex.nTracks();

    beamSpotTab->addColumn<float>("x", {x}, "",  nanoaod::FlatTable::FloatColumn);
    beamSpotTab->addColumn<float>("y", {y}, "",  nanoaod::FlatTable::FloatColumn);
    beamSpotTab->addColumn<float>("z", {z}, "",  nanoaod::FlatTable::FloatColumn);
    beamSpotTab->addColumn<float>("ndof", {ndof}, "",  nanoaod::FlatTable::FloatColumn);
    beamSpotTab->addColumn<float>("chi2", {chi2}, "",  nanoaod::FlatTable::FloatColumn);
    beamSpotTab->addColumn<uint8_t>("ntracks", {ntracks}, "",  nanoaod::FlatTable::UInt8Column);

    iEvent.put(std::move(beamSpotTab), "BS");
  }

};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BeamSpotTableProducer);
