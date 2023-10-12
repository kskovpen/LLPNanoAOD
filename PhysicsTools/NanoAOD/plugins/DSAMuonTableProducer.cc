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
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <vector>
#include <iostream>

class DSAMuonTableProducer : public edm::stream::EDProducer<> {
protected:
  const edm::InputTag displacedMuonLabel;
  const edm::EDGetTokenT<std::vector<reco::Track>> displacedMuonTag_;
  const edm::InputTag primaryVertexLabel;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVertexTag_;
  const edm::InputTag beamSpotLabel;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotTag_;
  
public:
  DSAMuonTableProducer(edm::ParameterSet const& params)
    :
    displacedMuonLabel(params.getParameter<edm::InputTag>("displacedMuons")),
    displacedMuonTag_(consumes<std::vector<reco::Track>>(displacedMuonLabel)),
    primaryVertexLabel(params.getParameter<edm::InputTag>("primaryVertex")),
    primaryVertexTag_(consumes<reco::VertexCollection>(primaryVertexLabel)),
    beamSpotLabel(params.getParameter<edm::InputTag>("beamSpot")),
    beamSpotTag_(consumes<reco::BeamSpot>(beamSpotLabel))
    {
    produces<nanoaod::FlatTable>("DSAMuon");
    
  }

  ~DSAMuonTableProducer() override {}

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    

    edm::Handle<std::vector<reco::Track>> displacedMuonHandle;
    iEvent.getByToken(displacedMuonTag_, displacedMuonHandle);
    auto displacedMuonTab = std::make_unique<nanoaod::FlatTable>(displacedMuonHandle->size(), "DSAMuon", false, false);

    edm::Handle<reco::VertexCollection> primaryVertexHandle;
    iEvent.getByToken(primaryVertexTag_, primaryVertexHandle);
    const auto& pv = primaryVertexHandle->at(0);
    GlobalPoint primaryVertex(pv.x(), pv.y(), pv.z());

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotTag_, beamSpotHandle);
    const auto& bs = beamSpotHandle->position();
    GlobalPoint beamSpot(bs.x(), bs.y(), bs.z());
    reco::Vertex beamSpotVertex(beamSpotHandle->position(), beamSpotHandle->covariance3D());

    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

    std::vector<float> pt,ptErr,eta,etaErr,phi,phiErr,charge,dxy,dz,vx,vy,vz,chi2,ndof;
    std::vector<float> dxyPV,dzPV,dxyTrajPV,dxyTrajPVErr,dxyTrajPVAbs,dxyTrajPVAbsErr,dxyTrajPVSigned,dxyTrajPVSignedErr;
    std::vector<float> ip3DPVAbs,ip3DPVAbsErr,ip3DPVSigned,ip3DPVSignedErr;
    std::vector<float> dxyBS,dzBS,dxyTrajBS,dxyTrajBSErr,dxyTrajBSAbs,dxyTrajBSAbsErr,dxyTrajBSSigned,dxyTrajBSSignedErr;
    std::vector<float> ip3DBSAbs,ip3DBSAbsErr,ip3DBSSigned,ip3DBSSignedErr;

    for(auto muon : *displacedMuonHandle) {
      pt.push_back(muon.pt());
      ptErr.push_back(muon.ptError());
      eta.push_back(muon.eta());
      etaErr.push_back(muon.etaError());
      phi.push_back(muon.phi());
      phiErr.push_back(muon.phiError());
      charge.push_back(muon.charge());
      dxy.push_back(muon.dxy());
      dz.push_back(muon.dz());
      vx.push_back(muon.vx());
      vy.push_back(muon.vy());
      vz.push_back(muon.vz());
      chi2.push_back(muon.chi2());
      ndof.push_back(muon.ndof());

      dxyPV.push_back(muon.dxy(pv.position()));
      dzPV.push_back(muon.dz(pv.position()));

      reco::TransientTrack transientTrack = builder->build(muon);
      TrajectoryStateClosestToPoint trajectoryPV = transientTrack.trajectoryStateClosestToPoint(primaryVertex);
      dxyTrajPV.push_back(trajectoryPV.perigeeParameters().transverseImpactParameter());
      dxyTrajPVErr.push_back(trajectoryPV.perigeeError().transverseImpactParameterError());

      dxyTrajPVAbs.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, pv).second.value());
      dxyTrajPVAbsErr.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, pv).second.error());
      GlobalVector muonRefTrackDir(muon.px(),muon.py(),muon.pz());
      dxyTrajPVSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, pv).second.value());
      dxyTrajPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, pv).second.error());

      ip3DPVAbs.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.value());
      ip3DPVAbsErr.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.error());
      ip3DPVSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.value());
      ip3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.error());  

      dxyBS.push_back(muon.dxy(bs));
      dzBS.push_back(muon.dz(bs));

      TrajectoryStateClosestToBeamLine trajectoryBS = transientTrack.stateAtBeamLine();
      dxyTrajBS.push_back(trajectoryBS.transverseImpactParameter().value());
      dxyTrajBSErr.push_back(trajectoryBS.transverseImpactParameter().error());

      dxyTrajBSAbs.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, beamSpotVertex).second.value());
      dxyTrajBSAbsErr.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, beamSpotVertex).second.error());
      dxyTrajBSSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, beamSpotVertex).second.value());
      dxyTrajBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, beamSpotVertex).second.error());  

      ip3DBSAbs.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.value());
      ip3DBSAbsErr.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.error());
      ip3DBSSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.value());
      ip3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.error());      

    }

    displacedMuonTab->addColumn<float>("pt", pt, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("ptErr", ptErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("eta", eta, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("etaErr", etaErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("phi", phi, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("phiErr", phiErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("charge", charge, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxy", dxy, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dz", dz, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("vx", vx, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("vy", vy, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("vz", vz, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("chi2", chi2, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("ndof", ndof, "",  nanoaod::FlatTable::FloatColumn);

    displacedMuonTab->addColumn<float>("dxyPV", dxyPV, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dzPV", dzPV, "",  nanoaod::FlatTable::FloatColumn);

    displacedMuonTab->addColumn<float>("dxyTrajPV", dxyTrajPV, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajPVErr", dxyTrajPVErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajPVAbs", dxyTrajPVAbs, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajPVAbsErr", dxyTrajPVAbsErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajPVSigned", dxyTrajPVSigned, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajPVSignedErr", dxyTrajPVSignedErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("ip3DPVSigned", ip3DPVSigned, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("ip3DPVSignedErr", ip3DPVSignedErr, "",  nanoaod::FlatTable::FloatColumn);

    displacedMuonTab->addColumn<float>("dxyBS", dxyBS, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dzBS", dzBS, "",  nanoaod::FlatTable::FloatColumn);

    displacedMuonTab->addColumn<float>("dxyTrajBS", dxyTrajBS, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajBSErr", dxyTrajBSErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajBSAbs", dxyTrajBSAbs, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajBSAbsErr", dxyTrajBSAbsErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajBSSigned", dxyTrajBSSigned, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("dxyTrajBSSignedErr", dxyTrajBSSignedErr, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("ip3DBSSigned", ip3DBSSigned, "",  nanoaod::FlatTable::FloatColumn);
    displacedMuonTab->addColumn<float>("ip3DBSSignedErr", ip3DBSSignedErr, "",  nanoaod::FlatTable::FloatColumn);

    iEvent.put(std::move(displacedMuonTab), "DSAMuon");
  }

};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DSAMuonTableProducer);
