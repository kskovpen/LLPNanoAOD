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

#include "DataFormats/PatCandidates/interface/Electron.h"

#include <vector>
#include <iostream>

class LowPtElectronExtendedTableProducer : public edm::stream::EDProducer<> {
protected:
  const edm::InputTag electronLabel;
  const edm::EDGetTokenT<std::vector<pat::Electron>> electronTag_;
  const edm::InputTag primaryVertexLabel;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVertexTag_;
  const edm::InputTag beamSpotLabel;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotTag_;
  
public:
  LowPtElectronExtendedTableProducer(edm::ParameterSet const& params)
    :
    electronLabel(params.getParameter<edm::InputTag>("lowPtElectrons")),
    electronTag_(consumes<std::vector<pat::Electron>>(electronLabel)),
    primaryVertexLabel(params.getParameter<edm::InputTag>("primaryVertex")),
    primaryVertexTag_(consumes<reco::VertexCollection>(primaryVertexLabel)),
    beamSpotLabel(params.getParameter<edm::InputTag>("beamSpot")),
    beamSpotTag_(consumes<reco::BeamSpot>(beamSpotLabel))
    {
    produces<nanoaod::FlatTable>("LowPtElectronExtended");
    
  }

  ~LowPtElectronExtendedTableProducer() override {}

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    

    edm::Handle<std::vector<pat::Electron>> electronHandle;
    iEvent.getByToken(electronTag_, electronHandle);
    auto electronTab = std::make_unique<nanoaod::FlatTable>(electronHandle->size(), "LowPtElectronExtended", false, false);

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

    std::vector<float> vx,vy,vz;
    std::vector<float> normChi2;
    // std::vector<float> outerEta,outerPhi;
    std::vector<float> dxyPV,dxyPVErr,dzPV,dzPVErr,dxyPVTraj,dxyPVTrajErr,dxyPVAbs,dxyPVAbsErr,dxyPVSigned,dxyPVSignedErr;
    std::vector<float> ip3DPVAbs,ip3DPVAbsErr,ip3DPVSigned,ip3DPVSignedErr;
    std::vector<float> dxyBS,dxyBSErr,dzBS,dzBSErr,dxyBSTraj,dxyBSTrajErr,dxyBSAbs,dxyBSAbsErr,dxyBSSigned,dxyBSSignedErr;
    std::vector<float> ip3DBSAbs,ip3DBSAbsErr,ip3DBSSigned,ip3DBSSignedErr;
    std::vector<float> dxyBS_test;

    for(auto electron : *electronHandle) {

      const auto& track = electron.bestTrack();
      reco::TransientTrack transientTrack = builder->build(track);

      vx.push_back(track->vx());
      vy.push_back(track->vy());
      vz.push_back(track->vz());

      normChi2.push_back(track->normalizedChi2());

      // outerEta.push_back(track->outerEta());
      // outerPhi.push_back(track->outerPhi());

      dxyPV.push_back(track->dxy(pv.position()));
      dxyPVErr.push_back(track->dxyError(pv.position(), pv.covariance()));
      dzPV.push_back(track->dz(pv.position()));
      dzPVErr.push_back(std::hypot(track->dzError(), pv.zError()));

      TrajectoryStateClosestToPoint trajectoryPV = transientTrack.trajectoryStateClosestToPoint(primaryVertex);
      dxyPVTraj.push_back(trajectoryPV.perigeeParameters().transverseImpactParameter());
      dxyPVTrajErr.push_back(trajectoryPV.perigeeError().transverseImpactParameterError());

      dxyPVAbs.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, pv).second.value());
      dxyPVAbsErr.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, pv).second.error());
      GlobalVector electronRefTrackDir(electron.px(),electron.py(),electron.pz());
      dxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, electronRefTrackDir, pv).second.value());
      dxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, electronRefTrackDir, pv).second.error());

      ip3DPVAbs.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.value());
      ip3DPVAbsErr.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.error());
      ip3DPVSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, electronRefTrackDir, beamSpotVertex).second.value());
      ip3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, electronRefTrackDir, beamSpotVertex).second.error());  

      dxyBS_test.push_back(electron.dB(pat::Electron::BS2D));

      dxyBS.push_back(track->dxy(bs));
      dxyBSErr.push_back(track->dxyError(bs, beamSpotVertex.covariance()));
      dzBS.push_back(track->dz(bs));
      dzBSErr.push_back(std::hypot(track->dzError(), beamSpotVertex.zError()));

      TrajectoryStateClosestToBeamLine trajectoryBS = transientTrack.stateAtBeamLine();
      dxyBSTraj.push_back(trajectoryBS.transverseImpactParameter().value());
      dxyBSTrajErr.push_back(trajectoryBS.transverseImpactParameter().error());

      dxyBSAbs.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, beamSpotVertex).second.value());
      dxyBSAbsErr.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, beamSpotVertex).second.error());
      dxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, electronRefTrackDir, beamSpotVertex).second.value());
      dxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, electronRefTrackDir, beamSpotVertex).second.error());  

      ip3DBSAbs.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.value());
      ip3DBSAbsErr.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.error());
      ip3DBSSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, electronRefTrackDir, beamSpotVertex).second.value());
      ip3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, electronRefTrackDir, beamSpotVertex).second.error());      

    }

    electronTab->addColumn<float>("vx", vx, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("vy", vy, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("vz", vz, "",  nanoaod::FlatTable::FloatColumn);

    electronTab->addColumn<float>("normChi2", normChi2, "",  nanoaod::FlatTable::FloatColumn);

    // electronTab->addColumn<float>("outerEta", outerEta, "",  nanoaod::FlatTable::FloatColumn);
    // electronTab->addColumn<float>("outerPhi", outerPhi, "",  nanoaod::FlatTable::FloatColumn);

    electronTab->addColumn<float>("dxyPV", dxyPV, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyPVErr", dxyPVErr, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dzPV", dzPV, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dzPVErr", dzPVErr, "",  nanoaod::FlatTable::FloatColumn);

    electronTab->addColumn<float>("dxyPVTraj", dxyPVTraj, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyPVTrajErr", dxyPVTrajErr, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyPVAbs", dxyPVAbs, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyPVAbsErr", dxyPVAbsErr, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyPVSigned", dxyPVSigned, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyPVSignedErr", dxyPVSignedErr, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("ip3DPVSigned", ip3DPVSigned, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("ip3DPVSignedErr", ip3DPVSignedErr, "",  nanoaod::FlatTable::FloatColumn);

    electronTab->addColumn<float>("dxyBS_test", dxyBS_test, "",  nanoaod::FlatTable::FloatColumn);

    electronTab->addColumn<float>("dxyBS", dxyBS, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyBSErr", dxyBSErr, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dzBS", dzBS, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dzBSErr", dzBSErr, "",  nanoaod::FlatTable::FloatColumn);

    electronTab->addColumn<float>("dxyBSTraj", dxyBSTraj, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyBSTrajErr", dxyBSTrajErr, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyBSAbs", dxyBSAbs, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyBSAbsErr", dxyBSAbsErr, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyBSSigned", dxyBSSigned, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("dxyBSSignedErr", dxyBSSignedErr, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("ip3DBSSigned", ip3DBSSigned, "",  nanoaod::FlatTable::FloatColumn);
    electronTab->addColumn<float>("ip3DBSSignedErr", ip3DBSSignedErr, "",  nanoaod::FlatTable::FloatColumn);

    iEvent.put(std::move(electronTab), "LowPtElectronExtended");
  }

};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LowPtElectronExtendedTableProducer);
