#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#include <vector>
#include <iostream>

class DSAMuonTableProducer : public edm::stream::EDProducer<> {
protected:
  const edm::InputTag displacedMuonLabel;
  const edm::EDGetTokenT<std::vector<reco::Track>> displacedMuonTag_;
  
public:
  DSAMuonTableProducer(edm::ParameterSet const& params)
    :
    displacedMuonLabel(params.getParameter<edm::InputTag>("displacedMuons")),
    displacedMuonTag_(consumes<std::vector<reco::Track>>(displacedMuonLabel))
    {
    produces<nanoaod::FlatTable>("DSAMuon");
    
  }

  ~DSAMuonTableProducer() override {}

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    

    edm::Handle<std::vector<reco::Track>> displacedMuonHandle;
    iEvent.getByToken(displacedMuonTag_, displacedMuonHandle);
    
    auto displacedMuonTab = std::make_unique<nanoaod::FlatTable>(displacedMuonHandle->size(), "DSAMuon", false, false);

    std::vector<float> muonpt;
    for(auto muon : *displacedMuonHandle) {
      muonpt.push_back(muon.pt());
    }
    displacedMuonTab->addColumn<float>("pt", muonpt, "");

    iEvent.put(std::move(displacedMuonTab), "DSAMuon");
  }

};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DSAMuonTableProducer);
