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

#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

#include <vector>
#include <iostream>

class IsoTrackExtendedTableProducer : public edm::stream::EDProducer<> {
protected:
  const edm::InputTag isoTrackLabel;
  const edm::EDGetTokenT<std::vector<pat::IsolatedTrack>> isoTrackTag_;
  const edm::InputTag primaryVertexLabel;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVertexTag_;
  const edm::InputTag beamSpotLabel;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotTag_;
  
public:
  IsoTrackExtendedTableProducer(edm::ParameterSet const& params)
    :
    isoTrackLabel(params.getParameter<edm::InputTag>("isoTracks")),
    isoTrackTag_(consumes<std::vector<pat::IsolatedTrack>>(isoTrackLabel)),
    primaryVertexLabel(params.getParameter<edm::InputTag>("primaryVertex")),
    primaryVertexTag_(consumes<reco::VertexCollection>(primaryVertexLabel)),
    beamSpotLabel(params.getParameter<edm::InputTag>("beamSpot")),
    beamSpotTag_(consumes<reco::BeamSpot>(beamSpotLabel))
    {
    produces<nanoaod::FlatTable>("IsoTrackExtended");
    
  }

  ~IsoTrackExtendedTableProducer() override {}

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    

    edm::Handle<std::vector<pat::IsolatedTrack>> isoTrackHandle;
    iEvent.getByToken(isoTrackTag_, isoTrackHandle);
    auto isoTrackTab = std::make_unique<nanoaod::FlatTable>(isoTrackHandle->size(), "IsoTrackExtended", false, false);

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

    // std::vector<float> dxyPV,dzPV,dxyPVTraj,dxyPVTrajErr,dxyPVAbs,dxyPVAbsErr,dxyPVSigned,dxyPVSignedErr;
    // std::vector<float> ip3DPVAbs,ip3DPVAbsErr,ip3DPVSigned,ip3DPVSignedErr;
    // std::vector<float> dxyBS,dzBS,dxyBSTraj,dxyBSTrajErr,dxyBSAbs,dxyBSAbsErr,dxyBSSigned,dxyBSSignedErr;
    // std::vector<float> ip3DBSAbs,ip3DBSAbsErr,ip3DBSSigned,ip3DBSSignedErr;
    std::vector<float> dxyBS;

    for(auto isoTrack : *isoTrackHandle) {

      // const auto track = isoTrack.bestTrack();
      // reco::TransientTrack transientTrack = builder->build(track);

      // dxyPV.push_back(track->dxy(pv.position()));
      // dzPV.push_back(track->dz(pv.position()));

      // TrajectoryStateClosestToPoint trajectoryPV = transientTrack.trajectoryStateClosestToPoint(primaryVertex);
      // dxyPVTraj.push_back(trajectoryPV.perigeeParameters().transverseImpactParameter());
      // dxyPVTrajErr.push_back(trajectoryPV.perigeeError().transverseImpactParameterError());

      // dxyPVAbs.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, pv).second.value());
      // dxyPVAbsErr.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, pv).second.error());
      // GlobalVector isoTrackRefTrackDir(isoTrack.px(),isoTrack.py(),isoTrack.pz());
      // dxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, isoTrackRefTrackDir, pv).second.value());
      // dxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, isoTrackRefTrackDir, pv).second.error());

      // ip3DPVAbs.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.value());
      // ip3DPVAbsErr.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.error());
      // ip3DPVSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, isoTrackRefTrackDir, beamSpotVertex).second.value());
      // ip3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, isoTrackRefTrackDir, beamSpotVertex).second.error());  

      dxyBS.push_back(isoTrack.dxy());
      // dxyBS.push_back(track->dxy(bs));
      // dzBS.push_back(track->dz(bs));

      // TrajectoryStateClosestToBeamLine trajectoryBS = transientTrack.stateAtBeamLine();
      // dxyBSTraj.push_back(trajectoryBS.transverseImpactParameter().value());
      // dxyBSTrajErr.push_back(trajectoryBS.transverseImpactParameter().error());

      // dxyBSAbs.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, beamSpotVertex).second.value());
      // dxyBSAbsErr.push_back(IPTools::absoluteTransverseImpactParameter(transientTrack, beamSpotVertex).second.error());
      // dxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, isoTrackRefTrackDir, beamSpotVertex).second.value());
      // dxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, isoTrackRefTrackDir, beamSpotVertex).second.error());  

      // ip3DBSAbs.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.value());
      // ip3DBSAbsErr.push_back(IPTools::absoluteImpactParameter3D(transientTrack, beamSpotVertex).second.error());
      // ip3DBSSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, isoTrackRefTrackDir, beamSpotVertex).second.value());
      // ip3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, isoTrackRefTrackDir, beamSpotVertex).second.error());      

    }

    // isoTrackTab->addColumn<float>("dxyPV", dxyPV, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dzPV", dzPV, "",  nanoaod::FlatTable::FloatColumn);

    // isoTrackTab->addColumn<float>("dxyPVTraj", dxyPVTraj, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyPVTrajErr", dxyPVTrajErr, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyPVAbs", dxyPVAbs, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyPVAbsErr", dxyPVAbsErr, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyPVSigned", dxyPVSigned, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyPVSignedErr", dxyPVSignedErr, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("ip3DPVSigned", ip3DPVSigned, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("ip3DPVSignedErr", ip3DPVSignedErr, "",  nanoaod::FlatTable::FloatColumn);

    isoTrackTab->addColumn<float>("dxyBS", dxyBS, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dzBS", dzBS, "",  nanoaod::FlatTable::FloatColumn);

    // isoTrackTab->addColumn<float>("dxyBSTraj", dxyBSTraj, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyBSTrajErr", dxyBSTrajErr, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyBSAbs", dxyBSAbs, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyBSAbsErr", dxyBSAbsErr, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyBSSigned", dxyBSSigned, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("dxyBSSignedErr", dxyBSSignedErr, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("ip3DBSSigned", ip3DBSSigned, "",  nanoaod::FlatTable::FloatColumn);
    // isoTrackTab->addColumn<float>("ip3DBSSignedErr", ip3DBSSignedErr, "",  nanoaod::FlatTable::FloatColumn);

    iEvent.put(std::move(isoTrackTab), "IsoTrackExtended");
  }

};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IsoTrackExtendedTableProducer);
