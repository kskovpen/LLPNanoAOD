#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

#include <vector>
#include <iostream>

class MuonVertexTableProducer : public edm::global::EDProducer<> {

  public:
    explicit MuonVertexTableProducer(const edm::ParameterSet &iConfig)
      :
      dsaMuonTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dsaMuons"))),
      patMuonTag_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("patMuons"))),
      bsTag_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
      pvTag_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertex"))),
      generalTrackTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("generalTracks"))),
      propagatorToken_(esConsumes(edm::ESInputTag("", "SteppingHelixPropagatorAny"))),
      magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      tkerGeomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
      tkerTopoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder")))
      {
      produces<nanoaod::FlatTable>("PatMuonVertex");
      produces<nanoaod::FlatTable>("PatDSAMuonVertex");
      produces<nanoaod::FlatTable>("DSAMuonVertex");
      if(runRefittedTracks_) {
        produces<nanoaod::FlatTable>("PatMuonVertexRefittedTracks");
        produces<nanoaod::FlatTable>("PatDSAMuonVertexRefittedTracks");
        produces<nanoaod::FlatTable>("DSAMuonVertexRefittedTracks");
      }
    }

    ~MuonVertexTableProducer() override {}

    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("patMuons")->setComment("input pat muon collection");
      desc.add<edm::InputTag>("dsaMuons")->setComment("input displaced standalone muon collection");
      desc.add<edm::InputTag>("beamspot")->setComment("input beamspot collection");
      desc.add<edm::InputTag>("primaryVertex")->setComment("input primaryVertex collection");
      desc.add<edm::InputTag>("generalTracks")->setComment("input generalTracks collection");
      descriptions.add("muonVertexTables", desc);
    }

  private:
    void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

    std::pair<float,float> getVxy(const reco::Vertex muonVertex) const;
    std::pair<float,float> getVxyz(const reco::Vertex muonVertex) const;

    template <typename MuonType1 = pat::Muon, typename MuonType2>
    float getDisplacedTrackerIsolation(const std::vector<reco::Track>& generalTracks, const MuonType1& muon_1,
                                      const reco::Vertex muonVertex, const reco::BeamSpot& beamspot, 
                                      const MuonType2* muon_2 = nullptr, float maxDR = 0.3, float minDR = 0.01,
                                      float maxDz = 0.5, float maxDxy = 0.2) const;
    
    template <typename MuonType=reco::Track>
    float getProximityDeltaR(const MuonType& track, 
                            const MuonType& trackRef,
                            const edm::ESHandle<MagneticField>& magneticField,
                            const edm::ESHandle<Propagator>& propagator) const;

    std::tuple<float, float, GlobalPoint> getDistanceBetweenMuonTracks(const reco::Track& track1, 
                                          const reco::Track& track2,
                                          const edm::ESHandle<MagneticField>& magneticField) const;

    const edm::EDGetTokenT<std::vector<reco::Track>> dsaMuonTag_;
    const edm::EDGetTokenT<std::vector<pat::Muon>> patMuonTag_;
    const edm::EDGetTokenT<reco::BeamSpot> bsTag_;
    const edm::EDGetTokenT<reco::VertexCollection> pvTag_;
    const edm::EDGetTokenT<std::vector<reco::Track>> generalTrackTag_;
    // FIX ME - make this a configurable parameter
    const bool runRefittedTracks_ = true;
    const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorToken_;
    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkerGeomToken_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tkerTopoToken_;
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;

};

void MuonVertexTableProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const 
{
  
  edm::Handle<std::vector<reco::Track>> dsaMuonHandle;
  iEvent.getByToken(dsaMuonTag_, dsaMuonHandle);
  edm::Handle<std::vector<pat::Muon>> patMuonHandle;
  iEvent.getByToken(patMuonTag_, patMuonHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(bsTag_, beamSpotHandle);
  const auto& bs = beamSpotHandle->position();
  GlobalPoint beamSpot(bs.x(), bs.y(), bs.z());
  reco::Vertex beamSpotVertex(beamSpotHandle->position(), beamSpotHandle->covariance3D());

  edm::Handle<reco::VertexCollection> primaryVertexHandle;
  iEvent.getByToken(pvTag_, primaryVertexHandle);
  const auto& pv = primaryVertexHandle->at(0);
  GlobalPoint primaryVertex(pv.x(), pv.y(), pv.z());

  edm::Handle<std::vector<reco::Track>> generalTracks;
  iEvent.getByToken(generalTrackTag_, generalTracks);

  // edm::ESHandle<Propagator> propagator;
  // iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", propagator);
  
  // edm::ESHandle<MagneticField> magneticField;
  // iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  // edm::ESHandle<TransientTrackBuilder> builder;
  // iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  auto const& propagator = &iSetup.getData(propagatorToken_);
  auto const& magneticField = &iSetup.getData(magneticFieldToken_);

  auto const& tkerGeom = &iSetup.getData(tkerGeomToken_);
  auto const& tkerTopo = &iSetup.getData(tkerTopoToken_);
  edm::ESHandle<TransientTrackBuilder> builder = iSetup.getHandle(transientTrackBuilderToken_);

  KalmanVertexFitter vertexFitter(true);

  int nPatPatVertices = 0;
  int nPatDSAVertices = 0;
  int nDSADSAVertices = 0;
  int ppRefittedTrackIdx_counter = 0;
  int pdRefittedTrackIdx_counter = 0;
  int ddRefittedTrackIdx_counter = 0;

  std::vector<bool> ppVertexIsValid, pdVertexIsValid, ddVertexIsValid;
  std::vector<float> ppVxy,ppVxyz,ppVx,ppVy,ppVz,ppT,ppVxySigma,ppVxyzSigma,ppVxErr,ppVyErr,ppVzErr,ppTErr;
  std::vector<float> pdVxy,pdVxyz,pdVx,pdVy,pdVz,pdT,pdVxySigma,pdVxyzSigma,pdVxErr,pdVyErr,pdVzErr,pdTErr;
  std::vector<float> ddVxy,ddVxyz,ddVx,ddVy,ddVz,ddT,ddVxySigma,ddVxyzSigma,ddVxErr,ddVyErr,ddVzErr,ddTErr;
  std::vector<float> ppChi2,ppNdof,ppNormChi2,ppDR,ppOriginalMuonIdx1,ppOriginalMuonIdx2,ppRefittedTrackIdx1,ppRefittedTrackIdx2,ppIsDSAMuon1,ppIsDSAMuon2;
  std::vector<float> pdChi2,pdNdof,pdNormChi2,pdDR,pdOriginalMuonIdx1,pdOriginalMuonIdx2,pdRefittedTrackIdx1,pdRefittedTrackIdx2,pdIsDSAMuon1,pdIsDSAMuon2;
  std::vector<float> ddChi2,ddNdof,ddNormChi2,ddDR,ddOriginalMuonIdx1,ddOriginalMuonIdx2,ddRefittedTrackIdx1,ddRefittedTrackIdx2,ddIsDSAMuon1,ddIsDSAMuon2;
  std::vector<float> ppDisplacedTrackIso03Dimuon1,ppDisplacedTrackIso03Dimuon2,ppDisplacedTrackIso04Dimuon1,ppDisplacedTrackIso04Dimuon2;
  std::vector<float> pdDisplacedTrackIso03Dimuon1,pdDisplacedTrackIso03Dimuon2,pdDisplacedTrackIso04Dimuon1,pdDisplacedTrackIso04Dimuon2;
  std::vector<float> ddDisplacedTrackIso03Dimuon1,ddDisplacedTrackIso03Dimuon2,ddDisplacedTrackIso04Dimuon1,ddDisplacedTrackIso04Dimuon2;
  std::vector<float> ppDisplacedTrackIso03Muon1,ppDisplacedTrackIso03Muon2,ppDisplacedTrackIso04Muon1,ppDisplacedTrackIso04Muon2,ppProxDeltaR;
  std::vector<float> pdDisplacedTrackIso03Muon1,pdDisplacedTrackIso03Muon2,pdDisplacedTrackIso04Muon1,pdDisplacedTrackIso04Muon2,pdProxDeltaR;
  std::vector<float> ddDisplacedTrackIso03Muon1,ddDisplacedTrackIso03Muon2,ddDisplacedTrackIso04Muon1,ddDisplacedTrackIso04Muon2,ddProxDeltaR;
  std::vector<float> ppDCA,ppDCAstatus,ppDCAx,ppDCAy,ppDCAz;
  std::vector<float> pdDCA,pdDCAstatus,pdDCAx,pdDCAy,pdDCAz;
  std::vector<float> ddDCA,ddDCAstatus,ddDCAx,ddDCAy,ddDCAz;
  std::vector<float> ppHitsInFrontOfVert1,ppMissHitsAfterVert1,ppHitsInFrontOfVert2,ppMissHitsAfterVert2;
  std::vector<float> pdHitsInFrontOfVert1,pdMissHitsAfterVert1,pdHitsInFrontOfVert2,pdMissHitsAfterVert2;
  std::vector<float> ddHitsInFrontOfVert1,ddMissHitsAfterVert1,ddHitsInFrontOfVert2,ddMissHitsAfterVert2;

  std::vector<float> ppRefittedTrackIdx,ppRefittedTrackIsDSAMuon,ppRefittedTrackOriginalIdx,ppRefittedTrackPt,ppRefittedTrackPtErr,ppRefittedTrackPx,ppRefittedTrackPy,ppRefittedTrackPz;
  std::vector<float> pdRefittedTrackIdx,pdRefittedTrackIsDSAMuon,pdRefittedTrackOriginalIdx,pdRefittedTrackPt,pdRefittedTrackPtErr,pdRefittedTrackPx,pdRefittedTrackPy,pdRefittedTrackPz;
  std::vector<float> ddRefittedTrackIdx,ddRefittedTrackIsDSAMuon,ddRefittedTrackOriginalIdx,ddRefittedTrackPt,ddRefittedTrackPtErr,ddRefittedTrackPx,ddRefittedTrackPy,ddRefittedTrackPz;
  std::vector<float> ppRefittedTrackEta,ppRefittedTrackEtaErr,ppRefittedTrackPhi,ppRefittedTrackPhiErr,ppRefittedTrackCharge,ppRefittedTrackNormChi2,ppRefittedTrackNdof,ppRefittedTrackChi2;
  std::vector<float> pdRefittedTrackEta,pdRefittedTrackEtaErr,pdRefittedTrackPhi,pdRefittedTrackPhiErr,pdRefittedTrackCharge,pdRefittedTrackNormChi2,pdRefittedTrackNdof,pdRefittedTrackChi2;
  std::vector<float> ddRefittedTrackEta,ddRefittedTrackEtaErr,ddRefittedTrackPhi,ddRefittedTrackPhiErr,ddRefittedTrackCharge,ddRefittedTrackNormChi2,ddRefittedTrackNdof,ddRefittedTrackChi2;
  std::vector<float> ppRefittedTrackDzPV,ppRefittedTrackDzPVErr,ppRefittedTrackDxyPVTraj,ppRefittedTrackDxyPVTrajErr,ppRefittedTrackDxyPVSigned,ppRefittedTrackDxyPVSignedErr,ppRefittedTrackIp3DPVSigned,ppRefittedTrackIp3DPVSignedErr;
  std::vector<float> pdRefittedTrackDzPV,pdRefittedTrackDzPVErr,pdRefittedTrackDxyPVTraj,pdRefittedTrackDxyPVTrajErr,pdRefittedTrackDxyPVSigned,pdRefittedTrackDxyPVSignedErr,pdRefittedTrackIp3DPVSigned,pdRefittedTrackIp3DPVSignedErr;
  std::vector<float> ddRefittedTrackDzPV,ddRefittedTrackDzPVErr,ddRefittedTrackDxyPVTraj,ddRefittedTrackDxyPVTrajErr,ddRefittedTrackDxyPVSigned,ddRefittedTrackDxyPVSignedErr,ddRefittedTrackIp3DPVSigned,ddRefittedTrackIp3DPVSignedErr;
  std::vector<float> ppRefittedTrackDxyBS,ppRefittedTrackDxyBSErr,ppRefittedTrackDzBS,ppRefittedTrackDzBSErr,ppRefittedTrackDxyBSTraj,ppRefittedTrackDxyBSTrajErr,ppRefittedTrackDxyBSSigned,ppRefittedTrackDxyBSSignedErr,ppRefittedTrackIp3DBSSigned,ppRefittedTrackIp3DBSSignedErr;
  std::vector<float> pdRefittedTrackDxyBS,pdRefittedTrackDxyBSErr,pdRefittedTrackDzBS,pdRefittedTrackDzBSErr,pdRefittedTrackDxyBSTraj,pdRefittedTrackDxyBSTrajErr,pdRefittedTrackDxyBSSigned,pdRefittedTrackDxyBSSignedErr,pdRefittedTrackIp3DBSSigned,pdRefittedTrackIp3DBSSignedErr;
  std::vector<float> ddRefittedTrackDxyBS,ddRefittedTrackDxyBSErr,ddRefittedTrackDzBS,ddRefittedTrackDzBSErr,ddRefittedTrackDxyBSTraj,ddRefittedTrackDxyBSTrajErr,ddRefittedTrackDxyBSSigned,ddRefittedTrackDxyBSSignedErr,ddRefittedTrackIp3DBSSigned,ddRefittedTrackIp3DBSSignedErr;

  std::vector<float> ppRefittedTrackIso03Dimuon1,ppRefittedTrackIso03Dimuon2,ppRefittedTrackIso04Dimuon1,ppRefittedTrackIso04Dimuon2,ppRefittedTrackIso03Muon1,ppRefittedTrackIso03Muon2,ppRefittedTrackIso04Muon1,ppRefittedTrackIso04Muon2;
  std::vector<float> pdRefittedTrackIso03Dimuon1,pdRefittedTrackIso03Dimuon2,pdRefittedTrackIso04Dimuon1,pdRefittedTrackIso04Dimuon2,pdRefittedTrackIso03Muon1,pdRefittedTrackIso03Muon2,pdRefittedTrackIso04Muon1,pdRefittedTrackIso04Muon2;
  std::vector<float> ddRefittedTrackIso03Dimuon1,ddRefittedTrackIso03Dimuon2,ddRefittedTrackIso04Dimuon1,ddRefittedTrackIso04Dimuon2,ddRefittedTrackIso03Muon1,ddRefittedTrackIso03Muon2,ddRefittedTrackIso04Muon1,ddRefittedTrackIso04Muon2;

  // pat muons
  for(size_t i = 0; i < patMuonHandle->size(); i++){

    const pat::Muon& muon_i = patMuonHandle->at(i);
    
    reco::TrackRef trackRef_i;
    if(muon_i.isGlobalMuon()) trackRef_i = muon_i.combinedMuon();
    else if (muon_i.isStandAloneMuon()) trackRef_i = muon_i.standAloneMuon();
    else trackRef_i = muon_i.tunePMuonBestTrack();

    const auto& muonTrack_i = trackRef_i.get();
    reco::TransientTrack muonTransientTrack_i = builder->build(muonTrack_i);

    // pat-pat muon vertex
    for(size_t j = i+1; j < patMuonHandle->size(); j++){
    
      const pat::Muon& muon_j = patMuonHandle->at(j);

      reco::TrackRef trackRef_j;
      if(muon_j.isGlobalMuon()) trackRef_j = muon_j.combinedMuon();
      else if (muon_j.isStandAloneMuon()) trackRef_j = muon_j.standAloneMuon();
      else trackRef_j = muon_j.tunePMuonBestTrack();

      const auto& muonTrack_j = trackRef_j.get();
      reco::TransientTrack muonTransientTrack_j = builder->build(muonTrack_j);

      std::vector<reco::TransientTrack> muonTransientTracks{};
      muonTransientTracks.push_back(muonTransientTrack_i);
      muonTransientTracks.push_back(muonTransientTrack_j);

      TransientVertex transientMuonVertex = vertexFitter.vertex(muonTransientTracks);
      reco::Vertex muonVertex = reco::Vertex(transientMuonVertex);

      if (!transientMuonVertex.isValid()) continue;
      std::tuple<float,float,GlobalPoint> distanceTuple = getDistanceBetweenMuonTracks(*muonTrack_i, *muonTrack_j, magneticField);
      // if dca status is good but dca is more than 15 cm
      if(std::get<1>(distanceTuple) && std::get<0>(distanceTuple) > 15) continue;

      nPatPatVertices++;

      ppVertexIsValid.push_back(transientMuonVertex.isValid());
      std::pair<float,float> vxy = getVxy(muonVertex);
      ppVxy.push_back(vxy.first);
      ppVxySigma.push_back(vxy.second);
      std::pair<float,float> vxyz = getVxyz(muonVertex);
      ppVxyz.push_back(vxyz.first);
      ppVxyzSigma.push_back(vxyz.second);
      ppChi2.push_back(muonVertex.chi2());
      ppNdof.push_back(muonVertex.ndof());
      ppNormChi2.push_back(muonVertex.normalizedChi2());
      ppVx.push_back(muonVertex.x());
      ppVy.push_back(muonVertex.y());
      ppVz.push_back(muonVertex.z());
      ppT.push_back(muonVertex.t());
      ppVxErr.push_back(muonVertex.xError());
      ppVyErr.push_back(muonVertex.yError());
      ppVzErr.push_back(muonVertex.zError());
      ppTErr.push_back(muonVertex.tError());

      ppDR.push_back(reco::deltaR(muon_i, muon_j));
      ppOriginalMuonIdx1.push_back(i);
      ppOriginalMuonIdx2.push_back(j);
      ppIsDSAMuon1.push_back(0);
      ppIsDSAMuon2.push_back(0);

      ppDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_i, muonVertex,*beamSpotHandle.product(), &muon_j, 0.3));
      ppDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_i, muonVertex,*beamSpotHandle.product(), &muon_j, 0.4));
      ppDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_j, muonVertex,*beamSpotHandle.product(), &muon_i, 0.3));
      ppDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_j, muonVertex,*beamSpotHandle.product(), &muon_i, 0.4));
      ppDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_i, muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
      ppDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_i, muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
      ppDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_j, muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
      ppDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<pat::Muon,pat::Muon>(*generalTracks.product(), muon_j, muonVertex,*beamSpotHandle.product(), nullptr, 0.4));

      // PAT muons cannot be only tracker muons for proximity deltaR calculations
      if((muon_i.isGlobalMuon() || muon_i.isStandAloneMuon()) && (muon_j.isGlobalMuon() || muon_j.isStandAloneMuon())) {
        ppProxDeltaR.push_back(getProximityDeltaR(*muonTrack_i, *muonTrack_j, magneticField, propagator));
      }
      else ppProxDeltaR.push_back(-1);

      ppDCA.push_back(std::get<0>(distanceTuple));
      ppDCAstatus.push_back(std::get<1>(distanceTuple));
      ppDCAx.push_back(std::get<2>(distanceTuple).x());
      ppDCAy.push_back(std::get<2>(distanceTuple).y());
      ppDCAz.push_back(std::get<2>(distanceTuple).z());

      CheckHitPattern checkHitPattern;
      checkHitPattern.init(tkerTopo,*tkerGeom,*builder);
      // checkHitPattern.init(iSetup);
      if(muon_i.isTrackerMuon()) {
        CheckHitPattern::Result hitPattern_i = checkHitPattern(*muonTrack_i, transientMuonVertex.vertexState());
        ppHitsInFrontOfVert1.push_back(hitPattern_i.hitsInFrontOfVert);
        ppMissHitsAfterVert1.push_back(hitPattern_i.missHitsAfterVert);
      }
      else {
        ppHitsInFrontOfVert1.push_back(-1);
        ppMissHitsAfterVert1.push_back(-1);
      }
      if(muon_j.isTrackerMuon()) {
        CheckHitPattern::Result hitPattern_j = checkHitPattern(*muonTrack_j, transientMuonVertex.vertexState());
        ppHitsInFrontOfVert2.push_back(hitPattern_j.hitsInFrontOfVert);
        ppMissHitsAfterVert2.push_back(hitPattern_j.missHitsAfterVert);
      }
      else {
        ppHitsInFrontOfVert2.push_back(-1);
        ppMissHitsAfterVert2.push_back(-1);
      }

      if(runRefittedTracks_) {
        reco::TransientTrack refittedTrack_i = transientMuonVertex.refittedTrack(muonTransientTrack_i);
        reco::TransientTrack refittedTrack_j = transientMuonVertex.refittedTrack(muonTransientTrack_j);
        ppRefittedTrackIsDSAMuon.push_back(0);
        ppRefittedTrackOriginalIdx.push_back(i);
        ppRefittedTrackPt.push_back(refittedTrack_i.track().pt());
        ppRefittedTrackPtErr.push_back(refittedTrack_i.track().ptError());
        ppRefittedTrackPx.push_back(refittedTrack_i.track().px());
        ppRefittedTrackPy.push_back(refittedTrack_i.track().py());
        ppRefittedTrackPz.push_back(refittedTrack_i.track().pz());
        ppRefittedTrackEta.push_back(refittedTrack_i.track().eta());
        ppRefittedTrackEtaErr.push_back(refittedTrack_i.track().etaError());
        ppRefittedTrackPhi.push_back(refittedTrack_i.track().phi());
        ppRefittedTrackPhiErr.push_back(refittedTrack_i.track().phiError());
        ppRefittedTrackCharge.push_back(refittedTrack_i.track().charge());
        ppRefittedTrackNormChi2.push_back(refittedTrack_i.normalizedChi2());
        ppRefittedTrackNdof.push_back(refittedTrack_i.ndof());
        ppRefittedTrackChi2.push_back(refittedTrack_i.chi2());

        ppRefittedTrackDzPV.push_back(refittedTrack_i.track().dz(pv.position()));
        ppRefittedTrackDzPVErr.push_back(std::hypot(refittedTrack_i.track().dzError(), pv.zError()));
        TrajectoryStateClosestToPoint trajectoryPV_i = refittedTrack_i.trajectoryStateClosestToPoint(primaryVertex);
        ppRefittedTrackDxyPVTraj.push_back(trajectoryPV_i.perigeeParameters().transverseImpactParameter());
        ppRefittedTrackDxyPVTrajErr.push_back(trajectoryPV_i.perigeeError().transverseImpactParameterError());
        GlobalVector muonRefTrackDir_i(refittedTrack_i.track().px(),refittedTrack_i.track().py(),refittedTrack_i.track().pz());
        ppRefittedTrackDxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, pv).second.value());
        ppRefittedTrackDxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, pv).second.error());
        ppRefittedTrackIp3DPVSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, pv).second.value());
        ppRefittedTrackIp3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, pv).second.error());  
        ppRefittedTrackDxyBS.push_back(refittedTrack_i.track().dxy(bs));
        ppRefittedTrackDxyBSErr.push_back(std::hypot(refittedTrack_i.track().dxyError(), beamSpotVertex.zError()));
        ppRefittedTrackDzBS.push_back(refittedTrack_i.track().dz(bs));
        ppRefittedTrackDzBSErr.push_back(std::hypot(refittedTrack_i.track().dzError(), beamSpotVertex.zError()));
        TrajectoryStateClosestToBeamLine trajectoryBS_i = refittedTrack_i.stateAtBeamLine();
        ppRefittedTrackDxyBSTraj.push_back(trajectoryBS_i.transverseImpactParameter().value());
        ppRefittedTrackDxyBSTrajErr.push_back(trajectoryBS_i.transverseImpactParameter().error());
        ppRefittedTrackDxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.value());
        ppRefittedTrackDxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.error());  
        ppRefittedTrackIp3DBSSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.value());
        ppRefittedTrackIp3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.error());
        ppRefittedTrackIdx.push_back(ppRefittedTrackIdx_counter);
        ppRefittedTrackIdx1.push_back(ppRefittedTrackIdx_counter);
        ppRefittedTrackIdx_counter++;

        ppRefittedTrackIsDSAMuon.push_back(0);
        ppRefittedTrackOriginalIdx.push_back(j);
        ppRefittedTrackPt.push_back(refittedTrack_j.track().pt());
        ppRefittedTrackPtErr.push_back(refittedTrack_j.track().ptError());
        ppRefittedTrackPx.push_back(refittedTrack_j.track().px());
        ppRefittedTrackPy.push_back(refittedTrack_j.track().py());
        ppRefittedTrackPz.push_back(refittedTrack_j.track().pz());
        ppRefittedTrackEta.push_back(refittedTrack_j.track().eta());
        ppRefittedTrackEtaErr.push_back(refittedTrack_j.track().etaError());
        ppRefittedTrackPhi.push_back(refittedTrack_j.track().phi());
        ppRefittedTrackPhiErr.push_back(refittedTrack_j.track().phiError());
        ppRefittedTrackCharge.push_back(refittedTrack_j.track().charge());
        ppRefittedTrackNormChi2.push_back(refittedTrack_j.normalizedChi2());
        ppRefittedTrackNdof.push_back(refittedTrack_j.ndof());
        ppRefittedTrackChi2.push_back(refittedTrack_j.chi2());

        ppRefittedTrackDzPV.push_back(refittedTrack_j.track().dz(pv.position()));
        ppRefittedTrackDzPVErr.push_back(std::hypot(refittedTrack_j.track().dzError(), pv.zError()));
        TrajectoryStateClosestToPoint trajectoryPV_j = refittedTrack_j.trajectoryStateClosestToPoint(primaryVertex);
        ppRefittedTrackDxyPVTraj.push_back(trajectoryPV_j.perigeeParameters().transverseImpactParameter());
        ppRefittedTrackDxyPVTrajErr.push_back(trajectoryPV_j.perigeeError().transverseImpactParameterError());
        GlobalVector muonRefTrackDir_j(refittedTrack_j.track().px(),refittedTrack_j.track().py(),refittedTrack_j.track().pz());
        ppRefittedTrackDxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, pv).second.value());
        ppRefittedTrackDxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, pv).second.error());
        ppRefittedTrackIp3DPVSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, pv).second.value());
        ppRefittedTrackIp3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, pv).second.error());  
        ppRefittedTrackDxyBS.push_back(refittedTrack_j.track().dxy(bs));
        ppRefittedTrackDxyBSErr.push_back(std::hypot(refittedTrack_j.track().dxyError(), beamSpotVertex.zError()));
        ppRefittedTrackDzBS.push_back(refittedTrack_j.track().dz(bs));
        ppRefittedTrackDzBSErr.push_back(std::hypot(refittedTrack_j.track().dzError(), beamSpotVertex.zError()));
        TrajectoryStateClosestToBeamLine trajectoryBS_j = refittedTrack_j.stateAtBeamLine();
        ppRefittedTrackDxyBSTraj.push_back(trajectoryBS_j.transverseImpactParameter().value());
        ppRefittedTrackDxyBSTrajErr.push_back(trajectoryBS_j.transverseImpactParameter().error());
        ppRefittedTrackDxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.value());
        ppRefittedTrackDxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.error());  
        ppRefittedTrackIp3DBSSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.value());
        ppRefittedTrackIp3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.error());
        ppRefittedTrackIdx.push_back(ppRefittedTrackIdx_counter);
        ppRefittedTrackIdx2.push_back(ppRefittedTrackIdx_counter);
        ppRefittedTrackIdx_counter++;

        ppRefittedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_j.track(), 0.3));
        ppRefittedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_j.track(), 0.4));
        ppRefittedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_i.track(), 0.3));
        ppRefittedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_i.track(), 0.4));
        ppRefittedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
        ppRefittedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
        ppRefittedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
        ppRefittedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
      }
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

      if (!transientMuonVertex.isValid()) continue;
      std::tuple<float,float,GlobalPoint> distanceTuple = getDistanceBetweenMuonTracks(*muonTrack_i, muonTrack_j, magneticField);
      // if dca status is good but dca is more than 15 cm
      if(std::get<1>(distanceTuple) && std::get<0>(distanceTuple) > 15) continue;

      nPatDSAVertices++;

      pdVertexIsValid.push_back(transientMuonVertex.isValid());
      std::pair<float,float> vxy = getVxy(muonVertex);
      pdVxy.push_back(vxy.first);
      pdVxySigma.push_back(vxy.second);
      std::pair<float,float> vxyz = getVxyz(muonVertex);
      pdVxyz.push_back(vxyz.first);
      pdVxyzSigma.push_back(vxyz.second);
      pdChi2.push_back(muonVertex.chi2());
      pdNdof.push_back(muonVertex.ndof());
      pdNormChi2.push_back(muonVertex.normalizedChi2());
      pdVx.push_back(muonVertex.x());
      pdVy.push_back(muonVertex.y());
      pdVz.push_back(muonVertex.z());
      pdT.push_back(muonVertex.t());
      pdVxErr.push_back(muonVertex.xError());
      pdVyErr.push_back(muonVertex.yError());
      pdVzErr.push_back(muonVertex.zError());
      pdTErr.push_back(muonVertex.tError());

      pdDR.push_back(reco::deltaR(muon_i.eta(), muonTrack_j.eta(), muon_i.phi(), muonTrack_j.phi()));
      pdOriginalMuonIdx1.push_back(i);
      pdOriginalMuonIdx2.push_back(j);
      pdIsDSAMuon1.push_back(0);
      pdIsDSAMuon2.push_back(1);

      pdDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamSpotHandle.product(), &muonTrack_j, 0.3));
      pdDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamSpotHandle.product(), &muonTrack_j, 0.4));
      pdDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamSpotHandle.product(), &muon_i, 0.3));
      pdDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamSpotHandle.product(), &muon_i, 0.4));
      pdDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
      pdDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<pat::Muon,reco::Track>(*generalTracks.product(), muon_i, muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
      pdDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
      pdDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,pat::Muon>(*generalTracks.product(), muonTrack_j, muonVertex,*beamSpotHandle.product(), nullptr, 0.4));

      // PAT muons cannot be only tracker muons for proximity deltaR calculations
      if(muon_i.isGlobalMuon() || muon_i.isStandAloneMuon()) {
        pdProxDeltaR.push_back(getProximityDeltaR(*muonTrack_i, muonTrack_j, magneticField, propagator));
      }
      else pdProxDeltaR.push_back(-1);

      pdDCA.push_back(std::get<0>(distanceTuple));
      pdDCAstatus.push_back(std::get<1>(distanceTuple));
      pdDCAx.push_back(std::get<2>(distanceTuple).x());
      pdDCAy.push_back(std::get<2>(distanceTuple).y());
      pdDCAz.push_back(std::get<2>(distanceTuple).z());

      CheckHitPattern checkHitPattern;
      checkHitPattern.init(tkerTopo,*tkerGeom,*builder);
      // checkHitPattern.init(iSetup);
      if(muon_i.isTrackerMuon()) {
        CheckHitPattern::Result hitPattern_i = checkHitPattern(*muonTrack_i, transientMuonVertex.vertexState());
        pdHitsInFrontOfVert1.push_back(hitPattern_i.hitsInFrontOfVert);
        pdMissHitsAfterVert1.push_back(hitPattern_i.missHitsAfterVert);
      }
      else {
        pdHitsInFrontOfVert1.push_back(-1);
        pdMissHitsAfterVert1.push_back(-1);
      }
      pdHitsInFrontOfVert2.push_back(-1);
      pdMissHitsAfterVert2.push_back(-1);

      if(runRefittedTracks_) {
        reco::TransientTrack refittedTrack_i = transientMuonVertex.refittedTrack(muonTransientTrack_i);
        reco::TransientTrack refittedTrack_j = transientMuonVertex.refittedTrack(muonTransientTrack_j);
        pdRefittedTrackIsDSAMuon.push_back(0);
        pdRefittedTrackOriginalIdx.push_back(i);
        pdRefittedTrackPt.push_back(refittedTrack_i.track().pt());
        pdRefittedTrackPtErr.push_back(refittedTrack_i.track().ptError());
        pdRefittedTrackPx.push_back(refittedTrack_i.track().px());
        pdRefittedTrackPy.push_back(refittedTrack_i.track().py());
        pdRefittedTrackPz.push_back(refittedTrack_i.track().pz());
        pdRefittedTrackEta.push_back(refittedTrack_i.track().eta());
        pdRefittedTrackEtaErr.push_back(refittedTrack_i.track().etaError());
        pdRefittedTrackPhi.push_back(refittedTrack_i.track().phi());
        pdRefittedTrackPhiErr.push_back(refittedTrack_i.track().phiError());
        pdRefittedTrackCharge.push_back(refittedTrack_i.track().charge());
        pdRefittedTrackNormChi2.push_back(refittedTrack_i.normalizedChi2());
        pdRefittedTrackNdof.push_back(refittedTrack_i.ndof());
        pdRefittedTrackChi2.push_back(refittedTrack_i.chi2());

        pdRefittedTrackDzPV.push_back(refittedTrack_i.track().dz(pv.position()));
        pdRefittedTrackDzPVErr.push_back(std::hypot(refittedTrack_i.track().dzError(), pv.zError()));
        TrajectoryStateClosestToPoint trajectoryPV_i = refittedTrack_i.trajectoryStateClosestToPoint(primaryVertex);
        pdRefittedTrackDxyPVTraj.push_back(trajectoryPV_i.perigeeParameters().transverseImpactParameter());
        pdRefittedTrackDxyPVTrajErr.push_back(trajectoryPV_i.perigeeError().transverseImpactParameterError());
        GlobalVector muonRefTrackDir_i(refittedTrack_i.track().px(),refittedTrack_i.track().py(),refittedTrack_i.track().pz());
        pdRefittedTrackDxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, pv).second.value());
        pdRefittedTrackDxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, pv).second.error());
        pdRefittedTrackIp3DPVSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, pv).second.value());
        pdRefittedTrackIp3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, pv).second.error());  
        pdRefittedTrackDxyBS.push_back(refittedTrack_i.track().dxy(bs));
        pdRefittedTrackDxyBSErr.push_back(std::hypot(refittedTrack_i.track().dxyError(), beamSpotVertex.zError()));
        pdRefittedTrackDzBS.push_back(refittedTrack_i.track().dz(bs));
        pdRefittedTrackDzBSErr.push_back(std::hypot(refittedTrack_i.track().dzError(), beamSpotVertex.zError()));
        TrajectoryStateClosestToBeamLine trajectoryBS_i = refittedTrack_i.stateAtBeamLine();
        pdRefittedTrackDxyBSTraj.push_back(trajectoryBS_i.transverseImpactParameter().value());
        pdRefittedTrackDxyBSTrajErr.push_back(trajectoryBS_i.transverseImpactParameter().error());
        pdRefittedTrackDxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.value());
        pdRefittedTrackDxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.error());  
        pdRefittedTrackIp3DBSSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.value());
        pdRefittedTrackIp3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.error());
        pdRefittedTrackIdx.push_back(pdRefittedTrackIdx_counter);
        pdRefittedTrackIdx1.push_back(pdRefittedTrackIdx_counter);
        pdRefittedTrackIdx_counter++;

        pdRefittedTrackIsDSAMuon.push_back(1);
        pdRefittedTrackOriginalIdx.push_back(j);
        pdRefittedTrackPt.push_back(refittedTrack_j.track().pt());
        pdRefittedTrackPtErr.push_back(refittedTrack_j.track().ptError());
        pdRefittedTrackPx.push_back(refittedTrack_j.track().px());
        pdRefittedTrackPy.push_back(refittedTrack_j.track().py());
        pdRefittedTrackPz.push_back(refittedTrack_j.track().pz());
        pdRefittedTrackEta.push_back(refittedTrack_j.track().eta());
        pdRefittedTrackEtaErr.push_back(refittedTrack_j.track().etaError());
        pdRefittedTrackPhi.push_back(refittedTrack_j.track().phi());
        pdRefittedTrackPhiErr.push_back(refittedTrack_j.track().phiError());
        pdRefittedTrackCharge.push_back(refittedTrack_j.track().charge());
        pdRefittedTrackNormChi2.push_back(refittedTrack_j.normalizedChi2());
        pdRefittedTrackNdof.push_back(refittedTrack_j.ndof());
        pdRefittedTrackChi2.push_back(refittedTrack_j.chi2());

        pdRefittedTrackDzPV.push_back(refittedTrack_j.track().dz(pv.position()));
        pdRefittedTrackDzPVErr.push_back(std::hypot(refittedTrack_j.track().dzError(), pv.zError()));
        TrajectoryStateClosestToPoint trajectoryPV_j = refittedTrack_j.trajectoryStateClosestToPoint(primaryVertex);
        pdRefittedTrackDxyPVTraj.push_back(trajectoryPV_j.perigeeParameters().transverseImpactParameter());
        pdRefittedTrackDxyPVTrajErr.push_back(trajectoryPV_j.perigeeError().transverseImpactParameterError());
        GlobalVector muonRefTrackDir_j(refittedTrack_j.track().px(),refittedTrack_j.track().py(),refittedTrack_j.track().pz());
        pdRefittedTrackDxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, pv).second.value());
        pdRefittedTrackDxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, pv).second.error());
        pdRefittedTrackIp3DPVSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, pv).second.value());
        pdRefittedTrackIp3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, pv).second.error());  
        pdRefittedTrackDxyBS.push_back(refittedTrack_j.track().dxy(bs));
        pdRefittedTrackDxyBSErr.push_back(std::hypot(refittedTrack_j.track().dxyError(), beamSpotVertex.zError()));
        pdRefittedTrackDzBS.push_back(refittedTrack_j.track().dz(bs));
        pdRefittedTrackDzBSErr.push_back(std::hypot(refittedTrack_j.track().dzError(), beamSpotVertex.zError()));
        TrajectoryStateClosestToBeamLine trajectoryBS_j = refittedTrack_j.stateAtBeamLine();
        pdRefittedTrackDxyBSTraj.push_back(trajectoryBS_j.transverseImpactParameter().value());
        pdRefittedTrackDxyBSTrajErr.push_back(trajectoryBS_j.transverseImpactParameter().error());
        pdRefittedTrackDxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.value());
        pdRefittedTrackDxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.error());  
        pdRefittedTrackIp3DBSSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.value());
        pdRefittedTrackIp3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.error());
        pdRefittedTrackIdx.push_back(pdRefittedTrackIdx_counter);
        pdRefittedTrackIdx2.push_back(pdRefittedTrackIdx_counter);
        pdRefittedTrackIdx_counter++;

        pdRefittedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_j.track(), 0.3));
        pdRefittedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_j.track(), 0.4));
        pdRefittedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_i.track(), 0.3));
        pdRefittedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_i.track(), 0.4));
        pdRefittedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
        pdRefittedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
        pdRefittedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
        pdRefittedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
      }
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

      if (!transientMuonVertex.isValid()) continue;
      std::tuple<float,float,GlobalPoint> distanceTuple = getDistanceBetweenMuonTracks(muonTrack_i, muonTrack_j, magneticField);
      // if dca status is good but dca is more than 15 cm
      if(std::get<1>(distanceTuple) && std::get<0>(distanceTuple) > 15) continue;

      nDSADSAVertices++;
      
      ddVertexIsValid.push_back(transientMuonVertex.isValid());

      std::pair<float,float> vxy = getVxy(muonVertex);
      ddVxy.push_back(vxy.first);
      ddVxySigma.push_back(vxy.second);
      std::pair<float,float> vxyz = getVxyz(muonVertex);
      ddVxyz.push_back(vxyz.first);
      ddVxyzSigma.push_back(vxyz.second);
      ddChi2.push_back(muonVertex.chi2());
      ddNdof.push_back(muonVertex.ndof());
      ddNormChi2.push_back(muonVertex.normalizedChi2());
      ddVx.push_back(muonVertex.x());
      ddVy.push_back(muonVertex.y());
      ddVz.push_back(muonVertex.z());
      ddT.push_back(muonVertex.t());
      ddVxErr.push_back(muonVertex.xError());
      ddVyErr.push_back(muonVertex.yError());
      ddVzErr.push_back(muonVertex.zError());
      ddTErr.push_back(muonVertex.tError());
      ddDR.push_back(reco::deltaR(muonTrack_i.eta(), muonTrack_j.eta(), muonTrack_i.phi(), muonTrack_j.phi()));
      ddOriginalMuonIdx1.push_back(i);
      ddOriginalMuonIdx2.push_back(j);
      ddIsDSAMuon1.push_back(1);
      ddIsDSAMuon2.push_back(1);

      ddDisplacedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamSpotHandle.product(), &muonTrack_j, 0.3));
      ddDisplacedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamSpotHandle.product(), &muonTrack_j, 0.4));
      ddDisplacedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamSpotHandle.product(), &muonTrack_i, 0.3));
      ddDisplacedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamSpotHandle.product(), &muonTrack_i, 0.4));
      ddDisplacedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
      ddDisplacedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_i, muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
      ddDisplacedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
      ddDisplacedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), muonTrack_j, muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
    
      ddProxDeltaR.push_back(getProximityDeltaR(muonTrack_i, muonTrack_j, magneticField, propagator));

      ddDCA.push_back(std::get<0>(distanceTuple));
      ddDCAstatus.push_back(std::get<1>(distanceTuple));
      ddDCAx.push_back(std::get<2>(distanceTuple).x());
      ddDCAy.push_back(std::get<2>(distanceTuple).y());
      ddDCAz.push_back(std::get<2>(distanceTuple).z());

      ddHitsInFrontOfVert1.push_back(-1);
      ddMissHitsAfterVert1.push_back(-1);
      ddHitsInFrontOfVert2.push_back(-1);
      ddMissHitsAfterVert2.push_back(-1);

      if(runRefittedTracks_) {
        reco::TransientTrack refittedTrack_i = transientMuonVertex.refittedTrack(muonTransientTrack_i);
        reco::TransientTrack refittedTrack_j = transientMuonVertex.refittedTrack(muonTransientTrack_j);
        ddRefittedTrackIsDSAMuon.push_back(1);
        ddRefittedTrackOriginalIdx.push_back(i);
        ddRefittedTrackPt.push_back(refittedTrack_i.track().pt());
        ddRefittedTrackPtErr.push_back(refittedTrack_i.track().ptError());
        ddRefittedTrackPx.push_back(refittedTrack_i.track().px());
        ddRefittedTrackPy.push_back(refittedTrack_i.track().py());
        ddRefittedTrackPz.push_back(refittedTrack_i.track().pz());
        ddRefittedTrackEta.push_back(refittedTrack_i.track().eta());
        ddRefittedTrackEtaErr.push_back(refittedTrack_i.track().etaError());
        ddRefittedTrackPhi.push_back(refittedTrack_i.track().phi());
        ddRefittedTrackPhiErr.push_back(refittedTrack_i.track().phiError());
        ddRefittedTrackCharge.push_back(refittedTrack_i.track().charge());
        ddRefittedTrackNormChi2.push_back(refittedTrack_i.normalizedChi2());
        ddRefittedTrackNdof.push_back(refittedTrack_i.ndof());
        ddRefittedTrackChi2.push_back(refittedTrack_i.chi2());
        
        ddRefittedTrackDzPV.push_back(refittedTrack_i.track().dz(pv.position()));
        ddRefittedTrackDzPVErr.push_back(std::hypot(refittedTrack_i.track().dzError(), pv.zError()));
        TrajectoryStateClosestToPoint trajectoryPV_i = refittedTrack_i.trajectoryStateClosestToPoint(primaryVertex);
        ddRefittedTrackDxyPVTraj.push_back(trajectoryPV_i.perigeeParameters().transverseImpactParameter());
        ddRefittedTrackDxyPVTrajErr.push_back(trajectoryPV_i.perigeeError().transverseImpactParameterError());
        GlobalVector muonRefTrackDir_i(refittedTrack_i.track().px(),refittedTrack_i.track().py(),refittedTrack_i.track().pz());
        ddRefittedTrackDxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, pv).second.value());
        ddRefittedTrackDxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, pv).second.error());
        ddRefittedTrackIp3DPVSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, pv).second.value());
        ddRefittedTrackIp3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, pv).second.error());  
        ddRefittedTrackDxyBS.push_back(refittedTrack_i.track().dxy(bs));
        ddRefittedTrackDxyBSErr.push_back(std::hypot(refittedTrack_i.track().dxyError(), beamSpotVertex.zError()));
        ddRefittedTrackDzBS.push_back(refittedTrack_i.track().dz(bs));
        ddRefittedTrackDzBSErr.push_back(std::hypot(refittedTrack_i.track().dzError(), beamSpotVertex.zError()));
        TrajectoryStateClosestToBeamLine trajectoryBS_i = refittedTrack_i.stateAtBeamLine();
        ddRefittedTrackDxyBSTraj.push_back(trajectoryBS_i.transverseImpactParameter().value());
        ddRefittedTrackDxyBSTrajErr.push_back(trajectoryBS_i.transverseImpactParameter().error());
        ddRefittedTrackDxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.value());
        ddRefittedTrackDxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.error());  
        ddRefittedTrackIp3DBSSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.value());
        ddRefittedTrackIp3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_i, muonRefTrackDir_i, beamSpotVertex).second.error());
        ddRefittedTrackIdx.push_back(ddRefittedTrackIdx_counter);
        ddRefittedTrackIdx1.push_back(ddRefittedTrackIdx_counter);
        ddRefittedTrackIdx_counter++;

        ddRefittedTrackIsDSAMuon.push_back(1);
        ddRefittedTrackOriginalIdx.push_back(j);
        ddRefittedTrackPt.push_back(refittedTrack_j.track().pt());
        ddRefittedTrackPtErr.push_back(refittedTrack_j.track().ptError());
        ddRefittedTrackPx.push_back(refittedTrack_j.track().px());
        ddRefittedTrackPy.push_back(refittedTrack_j.track().py());
        ddRefittedTrackPz.push_back(refittedTrack_j.track().pz());
        ddRefittedTrackEta.push_back(refittedTrack_j.track().eta());
        ddRefittedTrackEtaErr.push_back(refittedTrack_j.track().etaError());
        ddRefittedTrackPhi.push_back(refittedTrack_j.track().phi());
        ddRefittedTrackPhiErr.push_back(refittedTrack_j.track().phiError());
        ddRefittedTrackCharge.push_back(refittedTrack_j.track().charge());
        ddRefittedTrackNormChi2.push_back(refittedTrack_j.normalizedChi2());
        ddRefittedTrackNdof.push_back(refittedTrack_j.ndof());
        ddRefittedTrackChi2.push_back(refittedTrack_j.chi2());
        
        ddRefittedTrackDzPV.push_back(refittedTrack_j.track().dz(pv.position()));
        ddRefittedTrackDzPVErr.push_back(std::hypot(refittedTrack_j.track().dzError(), pv.zError()));
        TrajectoryStateClosestToPoint trajectoryPV_j = refittedTrack_j.trajectoryStateClosestToPoint(primaryVertex);
        ddRefittedTrackDxyPVTraj.push_back(trajectoryPV_j.perigeeParameters().transverseImpactParameter());
        ddRefittedTrackDxyPVTrajErr.push_back(trajectoryPV_j.perigeeError().transverseImpactParameterError());
        GlobalVector muonRefTrackDir_j(refittedTrack_j.track().px(),refittedTrack_j.track().py(),refittedTrack_j.track().pz());
        ddRefittedTrackDxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, pv).second.value());
        ddRefittedTrackDxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, pv).second.error());
        ddRefittedTrackIp3DPVSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, pv).second.value());
        ddRefittedTrackIp3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, pv).second.error());  
        ddRefittedTrackDxyBS.push_back(refittedTrack_j.track().dxy(bs));
        ddRefittedTrackDxyBSErr.push_back(std::hypot(refittedTrack_j.track().dxyError(), beamSpotVertex.zError()));
        ddRefittedTrackDzBS.push_back(refittedTrack_j.track().dz(bs));
        ddRefittedTrackDzBSErr.push_back(std::hypot(refittedTrack_j.track().dzError(), beamSpotVertex.zError()));
        TrajectoryStateClosestToBeamLine trajectoryBS_j = refittedTrack_j.stateAtBeamLine();
        ddRefittedTrackDxyBSTraj.push_back(trajectoryBS_j.transverseImpactParameter().value());
        ddRefittedTrackDxyBSTrajErr.push_back(trajectoryBS_j.transverseImpactParameter().error());
        ddRefittedTrackDxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.value());
        ddRefittedTrackDxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.error());  
        ddRefittedTrackIp3DBSSigned.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.value());
        ddRefittedTrackIp3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(refittedTrack_j, muonRefTrackDir_j, beamSpotVertex).second.error());
        ddRefittedTrackIdx.push_back(ddRefittedTrackIdx_counter);
        ddRefittedTrackIdx2.push_back(ddRefittedTrackIdx_counter);
        ddRefittedTrackIdx_counter++;

        ddRefittedTrackIso03Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_j.track(), 0.3));
        ddRefittedTrackIso04Dimuon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_j.track(), 0.4));
        ddRefittedTrackIso03Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_i.track(), 0.3));
        ddRefittedTrackIso04Dimuon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), &refittedTrack_i.track(), 0.4));
        ddRefittedTrackIso03Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
        ddRefittedTrackIso04Muon1.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_i.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
        ddRefittedTrackIso03Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.3));
        ddRefittedTrackIso04Muon2.push_back(getDisplacedTrackerIsolation<reco::Track,reco::Track>(*generalTracks.product(), refittedTrack_j.track(), muonVertex,*beamSpotHandle.product(), nullptr, 0.4));
      }
    }
  }
  
  // auto patVertexTab = std::make_unique<nanoaod::FlatTable>(patMuonHandle->size() * (patMuonHandle->size() - 1) / 2, "PatMuonVertex", false, false);
  // auto patdsaVertexTab = std::make_unique<nanoaod::FlatTable>(patMuonHandle->size() * dsaMuonHandle->size(), "PatDSAMuonVertex", false, false);
  // auto dsaVertexTab = std::make_unique<nanoaod::FlatTable>(dsaMuonHandle->size() * (dsaMuonHandle->size() - 1) / 2, "DSAMuonVertex", false, false);
  auto patVertexTab = std::make_unique<nanoaod::FlatTable>(nPatPatVertices, "PatMuonVertex", false, false);
  auto patdsaVertexTab = std::make_unique<nanoaod::FlatTable>(nPatDSAVertices, "PatDSAMuonVertex", false, false);
  auto dsaVertexTab = std::make_unique<nanoaod::FlatTable>(nDSADSAVertices, "DSAMuonVertex", false, false);

  auto patRefittedTracksTab = std::make_unique<nanoaod::FlatTable>(nPatPatVertices*2, "PatMuonVertexRefittedTracks", false, false);
  auto patdsaRefittedTracksTab = std::make_unique<nanoaod::FlatTable>(nPatDSAVertices*2, "PatDSAMuonVertexRefittedTracks", false, false);
  auto dsaRefittedTracksTab = std::make_unique<nanoaod::FlatTable>(nDSADSAVertices*2, "DSAMuonVertexRefittedTracks", false, false);

  patVertexTab->addColumn<float>("isValid", ppVertexIsValid, "");
  patVertexTab->addColumn<float>("vxy", ppVxy, "");
  patVertexTab->addColumn<float>("vxySigma", ppVxySigma, "");
  patVertexTab->addColumn<float>("vxyz", ppVxyz, "");
  patVertexTab->addColumn<float>("vxyzSigma", ppVxyzSigma, "");
  patVertexTab->addColumn<float>("chi2", ppChi2, "");
  patVertexTab->addColumn<float>("ndof", ppNdof, "");
  patVertexTab->addColumn<float>("normChi2", ppNormChi2, "");
  patVertexTab->addColumn<float>("vx", ppVx, "");
  patVertexTab->addColumn<float>("vy", ppVy, "");
  patVertexTab->addColumn<float>("vz", ppVz, "");
  patVertexTab->addColumn<float>("t", ppT, "");
  patVertexTab->addColumn<float>("vxErr", ppVxErr, "");
  patVertexTab->addColumn<float>("vyErr", ppVyErr, "");
  patVertexTab->addColumn<float>("vzErr", ppVzErr, "");
  patVertexTab->addColumn<float>("tErr", ppTErr, "");
  patVertexTab->addColumn<float>("dR", ppDR, "");
  patVertexTab->addColumn<float>("originalMuonIdx1", ppOriginalMuonIdx1, "");
  patVertexTab->addColumn<float>("originalMuonIdx2", ppOriginalMuonIdx2, "");
  patVertexTab->addColumn<float>("isDSAMuon1", ppIsDSAMuon1, "");
  patVertexTab->addColumn<float>("isDSAMuon2", ppIsDSAMuon2, "");
  patVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", ppDisplacedTrackIso03Dimuon1, "");
  patVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", ppDisplacedTrackIso04Dimuon1, "");
  patVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", ppDisplacedTrackIso03Dimuon2, "");
  patVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", ppDisplacedTrackIso04Dimuon2, "");
  patVertexTab->addColumn<float>("displacedTrackIso03Muon1", ppDisplacedTrackIso03Muon1, "");
  patVertexTab->addColumn<float>("displacedTrackIso04Muon1", ppDisplacedTrackIso04Muon1, "");
  patVertexTab->addColumn<float>("displacedTrackIso03Muon2", ppDisplacedTrackIso03Muon2, "");
  patVertexTab->addColumn<float>("displacedTrackIso04Muon2", ppDisplacedTrackIso04Muon2, "");
  patVertexTab->addColumn<float>("dRprox", ppProxDeltaR, "");
  patVertexTab->addColumn<float>("dca", ppDCA, "");
  patVertexTab->addColumn<float>("dcaStatus", ppDCAstatus, "");
  patVertexTab->addColumn<float>("dcax", ppDCAx, "");
  patVertexTab->addColumn<float>("dcay", ppDCAy, "");
  patVertexTab->addColumn<float>("dcaz", ppDCAz, "");
  patVertexTab->addColumn<float>("hitsInFrontOfVert1", ppHitsInFrontOfVert1, "");
  patVertexTab->addColumn<float>("hitsInFrontOfVert2", ppHitsInFrontOfVert2, "");
  patVertexTab->addColumn<float>("missHitsAfterVert1", ppMissHitsAfterVert1, "");
  patVertexTab->addColumn<float>("missHitsAfterVert2", ppMissHitsAfterVert2, "");
  if(runRefittedTracks_) {
    patVertexTab->addColumn<float>("refittedTrackIdx1", ppRefittedTrackIdx1, "");
    patVertexTab->addColumn<float>("refittedTrackIdx2", ppRefittedTrackIdx2, "");
    patVertexTab->addColumn<float>("refittedTrackIso03Dimuon1", ppRefittedTrackIso03Dimuon1, "");
    patVertexTab->addColumn<float>("refittedTrackIso04Dimuon1", ppRefittedTrackIso04Dimuon1, "");
    patVertexTab->addColumn<float>("refittedTrackIso03Dimuon2", ppRefittedTrackIso03Dimuon2, "");
    patVertexTab->addColumn<float>("refittedTrackIso04Dimuon2", ppRefittedTrackIso04Dimuon2, "");
    patVertexTab->addColumn<float>("refittedTrackIso03Muon1", ppRefittedTrackIso03Muon1, "");
    patVertexTab->addColumn<float>("refittedTrackIso04Muon1", ppRefittedTrackIso04Muon1, "");
    patVertexTab->addColumn<float>("refittedTrackIso03Muon2", ppRefittedTrackIso03Muon2, "");
    patVertexTab->addColumn<float>("refittedTrackIso04Muon2", ppRefittedTrackIso04Muon2, "");
  }

  patdsaVertexTab->addColumn<float>("isValid", pdVertexIsValid, "");
  patdsaVertexTab->addColumn<float>("vxy", pdVxy, "");
  patdsaVertexTab->addColumn<float>("vxySigma", pdVxySigma, "");
  patdsaVertexTab->addColumn<float>("vxyz", pdVxyz, "");
  patdsaVertexTab->addColumn<float>("vxyzSigma", pdVxyzSigma, "");
  patdsaVertexTab->addColumn<float>("chi2", pdChi2, "");
  patdsaVertexTab->addColumn<float>("ndof", pdNdof, "");
  patdsaVertexTab->addColumn<float>("normChi2", pdNormChi2, "");
  patdsaVertexTab->addColumn<float>("vx", pdVx, "");
  patdsaVertexTab->addColumn<float>("vy", pdVy, "");
  patdsaVertexTab->addColumn<float>("vz", pdVz, "");
  patdsaVertexTab->addColumn<float>("t", pdT, "");
  patdsaVertexTab->addColumn<float>("vxErr", pdVxErr, "");
  patdsaVertexTab->addColumn<float>("vyErr", pdVyErr, "");
  patdsaVertexTab->addColumn<float>("vzErr", pdVzErr, "");
  patdsaVertexTab->addColumn<float>("tErr", pdTErr, "");
  patdsaVertexTab->addColumn<float>("dR", pdDR, "");
  patdsaVertexTab->addColumn<float>("originalMuonIdx1", pdOriginalMuonIdx1, "");
  patdsaVertexTab->addColumn<float>("originalMuonIdx2", pdOriginalMuonIdx2, "");
  patdsaVertexTab->addColumn<float>("isDSAMuon1", pdIsDSAMuon1, "");
  patdsaVertexTab->addColumn<float>("isDSAMuon2", pdIsDSAMuon2, "");
  patdsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", pdDisplacedTrackIso03Dimuon1, "");
  patdsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", pdDisplacedTrackIso04Dimuon1, "");
  patdsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", pdDisplacedTrackIso03Dimuon2, "");
  patdsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", pdDisplacedTrackIso04Dimuon2, "");
  patdsaVertexTab->addColumn<float>("displacedTrackIso03Muon1", pdDisplacedTrackIso03Muon1, "");
  patdsaVertexTab->addColumn<float>("displacedTrackIso04Muon1", pdDisplacedTrackIso04Muon1, "");
  patdsaVertexTab->addColumn<float>("displacedTrackIso03Muon2", pdDisplacedTrackIso03Muon2, "");
  patdsaVertexTab->addColumn<float>("displacedTrackIso04Muon2", pdDisplacedTrackIso04Muon2, "");
  patdsaVertexTab->addColumn<float>("dRprox", pdProxDeltaR, "");
  patdsaVertexTab->addColumn<float>("dca", pdDCA, "");
  patdsaVertexTab->addColumn<float>("dcaStatus", pdDCAstatus, "");
  patdsaVertexTab->addColumn<float>("dcax", pdDCAx, "");
  patdsaVertexTab->addColumn<float>("dcay", pdDCAy, "");
  patdsaVertexTab->addColumn<float>("dcaz", pdDCAz, "");
  patdsaVertexTab->addColumn<float>("hitsInFrontOfVert1", pdHitsInFrontOfVert1, "");
  patdsaVertexTab->addColumn<float>("hitsInFrontOfVert2", pdHitsInFrontOfVert2, "");
  patdsaVertexTab->addColumn<float>("missHitsAfterVert1", pdMissHitsAfterVert1, "");
  patdsaVertexTab->addColumn<float>("missHitsAfterVert2", pdMissHitsAfterVert2, "");
  if(runRefittedTracks_) {
    patdsaVertexTab->addColumn<float>("refittedTrackIdx1", pdRefittedTrackIdx1, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIdx2", pdRefittedTrackIdx2, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIso03Dimuon1", pdRefittedTrackIso03Dimuon1, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIso04Dimuon1", pdRefittedTrackIso04Dimuon1, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIso03Dimuon2", pdRefittedTrackIso03Dimuon2, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIso04Dimuon2", pdRefittedTrackIso04Dimuon2, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIso03Muon1", pdRefittedTrackIso03Muon1, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIso04Muon1", pdRefittedTrackIso04Muon1, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIso03Muon2", pdRefittedTrackIso03Muon2, "");
    patdsaVertexTab->addColumn<float>("refittedTrackIso04Muon2", pdRefittedTrackIso04Muon2, "");
  }

  dsaVertexTab->addColumn<float>("isValid", ddVertexIsValid, "");
  dsaVertexTab->addColumn<float>("vxy", ddVxy, "");
  dsaVertexTab->addColumn<float>("vxySigma", ddVxySigma, "");
  dsaVertexTab->addColumn<float>("vxyz", ddVxyz, "");
  dsaVertexTab->addColumn<float>("vxyzSigma", ddVxyzSigma, "");
  dsaVertexTab->addColumn<float>("chi2", ddChi2, "");
  dsaVertexTab->addColumn<float>("ndof", ddNdof, "");
  dsaVertexTab->addColumn<float>("normChi2", ddNormChi2, "");
  dsaVertexTab->addColumn<float>("vx", ddVx, "");
  dsaVertexTab->addColumn<float>("vy", ddVy, "");
  dsaVertexTab->addColumn<float>("vz", ddVz, "");
  dsaVertexTab->addColumn<float>("t", ddT, "");
  dsaVertexTab->addColumn<float>("vxErr", ddVxErr, "");
  dsaVertexTab->addColumn<float>("vyErr", ddVyErr, "");
  dsaVertexTab->addColumn<float>("vzErr", ddVzErr, "");
  dsaVertexTab->addColumn<float>("tErr", ddTErr, "");
  dsaVertexTab->addColumn<float>("dR", ddDR, "");
  dsaVertexTab->addColumn<float>("originalMuonIdx1", ddOriginalMuonIdx1, "");
  dsaVertexTab->addColumn<float>("originalMuonIdx2", ddOriginalMuonIdx2, "");
  dsaVertexTab->addColumn<float>("isDSAMuon1", ddIsDSAMuon1, "");
  dsaVertexTab->addColumn<float>("isDSAMuon2", ddIsDSAMuon2, "");
  dsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon1", ddDisplacedTrackIso03Dimuon1, "");
  dsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon1", ddDisplacedTrackIso04Dimuon1, "");
  dsaVertexTab->addColumn<float>("displacedTrackIso03Dimuon2", ddDisplacedTrackIso03Dimuon2, "");
  dsaVertexTab->addColumn<float>("displacedTrackIso04Dimuon2", ddDisplacedTrackIso04Dimuon2, "");
  dsaVertexTab->addColumn<float>("displacedTrackIso03Muon1", ddDisplacedTrackIso03Muon1, "");
  dsaVertexTab->addColumn<float>("displacedTrackIso04Muon1", ddDisplacedTrackIso04Muon1, "");
  dsaVertexTab->addColumn<float>("displacedTrackIso03Muon2", ddDisplacedTrackIso03Muon2, "");
  dsaVertexTab->addColumn<float>("displacedTrackIso04Muon2", ddDisplacedTrackIso04Muon2, "");
  dsaVertexTab->addColumn<float>("dRprox", ddProxDeltaR, "");
  dsaVertexTab->addColumn<float>("dca", ddDCA, "");
  dsaVertexTab->addColumn<float>("dcaStatus", ddDCAstatus, "");
  dsaVertexTab->addColumn<float>("dcax", ddDCAx, "");
  dsaVertexTab->addColumn<float>("dcay", ddDCAy, "");
  dsaVertexTab->addColumn<float>("dcaz", ddDCAz, "");
  dsaVertexTab->addColumn<float>("hitsInFrontOfVert1", ddHitsInFrontOfVert1, "");
  dsaVertexTab->addColumn<float>("hitsInFrontOfVert2", ddHitsInFrontOfVert2, "");
  dsaVertexTab->addColumn<float>("missHitsAfterVert1", ddMissHitsAfterVert1, "");
  dsaVertexTab->addColumn<float>("missHitsAfterVert2", ddMissHitsAfterVert2, "");
  if(runRefittedTracks_) {
    dsaVertexTab->addColumn<float>("refittedTrackIdx1", ddRefittedTrackIdx1, "");
    dsaVertexTab->addColumn<float>("refittedTrackIdx2", ddRefittedTrackIdx2, "");
    dsaVertexTab->addColumn<float>("refittedTrackIso03Dimuon1", ddRefittedTrackIso03Dimuon1, "");
    dsaVertexTab->addColumn<float>("refittedTrackIso04Dimuon1", ddRefittedTrackIso04Dimuon1, "");
    dsaVertexTab->addColumn<float>("refittedTrackIso03Dimuon2", ddRefittedTrackIso03Dimuon2, "");
    dsaVertexTab->addColumn<float>("refittedTrackIso04Dimuon2", ddRefittedTrackIso04Dimuon2, "");
    dsaVertexTab->addColumn<float>("refittedTrackIso03Muon1", ddRefittedTrackIso03Muon1, "");
    dsaVertexTab->addColumn<float>("refittedTrackIso04Muon1", ddRefittedTrackIso04Muon1, "");
    dsaVertexTab->addColumn<float>("refittedTrackIso03Muon2", ddRefittedTrackIso03Muon2, "");
    dsaVertexTab->addColumn<float>("refittedTrackIso04Muon2", ddRefittedTrackIso04Muon2, "");
  }
  iEvent.put(std::move(patVertexTab), "PatMuonVertex");
  iEvent.put(std::move(patdsaVertexTab), "PatDSAMuonVertex");
  iEvent.put(std::move(dsaVertexTab), "DSAMuonVertex");

  if(runRefittedTracks_) {

    patRefittedTracksTab->addColumn<float>("idx", ppRefittedTrackIdx, "");
    patRefittedTracksTab->addColumn<float>("isDSAMuon", ppRefittedTrackIsDSAMuon, "");
    patRefittedTracksTab->addColumn<float>("originalMuonIdx", ppRefittedTrackOriginalIdx, "");
    patRefittedTracksTab->addColumn<float>("pt", ppRefittedTrackPt, "");
    patRefittedTracksTab->addColumn<float>("ptErr", ppRefittedTrackPtErr, "");
    patRefittedTracksTab->addColumn<float>("px", ppRefittedTrackPx, "");
    patRefittedTracksTab->addColumn<float>("py", ppRefittedTrackPy, "");
    patRefittedTracksTab->addColumn<float>("pz", ppRefittedTrackPz, "");
    patRefittedTracksTab->addColumn<float>("eta", ppRefittedTrackEta, "");
    patRefittedTracksTab->addColumn<float>("etaErr", ppRefittedTrackEtaErr, "");
    patRefittedTracksTab->addColumn<float>("phi", ppRefittedTrackPhi, "");
    patRefittedTracksTab->addColumn<float>("phiErr", ppRefittedTrackPhiErr, "");
    patRefittedTracksTab->addColumn<float>("charge", ppRefittedTrackCharge, "");
    patRefittedTracksTab->addColumn<float>("normChi2", ppRefittedTrackNormChi2, "");
    patRefittedTracksTab->addColumn<float>("ndof", ppRefittedTrackNdof, "");
    patRefittedTracksTab->addColumn<float>("chi2", ppRefittedTrackChi2, "");
    patRefittedTracksTab->addColumn<float>("dzPV", ppRefittedTrackDzPV, "");
    patRefittedTracksTab->addColumn<float>("dzPVErr", ppRefittedTrackDzPVErr, "");
    patRefittedTracksTab->addColumn<float>("dxyPVTraj", ppRefittedTrackDxyPVTraj, "");
    patRefittedTracksTab->addColumn<float>("dxyPVTrajErr", ppRefittedTrackDxyPVTrajErr, "");
    patRefittedTracksTab->addColumn<float>("dxyPVSigned", ppRefittedTrackDxyPVSigned, "");
    patRefittedTracksTab->addColumn<float>("dxyPVSignedErr", ppRefittedTrackDxyPVSignedErr, "");
    patRefittedTracksTab->addColumn<float>("ip3DPVSigned", ppRefittedTrackIp3DPVSigned, "");
    patRefittedTracksTab->addColumn<float>("ip3DPVSignedErr", ppRefittedTrackIp3DPVSignedErr, "");
    patRefittedTracksTab->addColumn<float>("dxyBS", ppRefittedTrackDxyBS, "");
    patRefittedTracksTab->addColumn<float>("dxyBSErr", ppRefittedTrackDxyBSErr, "");
    patRefittedTracksTab->addColumn<float>("dzBS", ppRefittedTrackDzBS, "");
    patRefittedTracksTab->addColumn<float>("dzBSErr", ppRefittedTrackDzBSErr, "");
    patRefittedTracksTab->addColumn<float>("dxyBSTraj", ppRefittedTrackDxyBSTraj, "");
    patRefittedTracksTab->addColumn<float>("dxyBSTrajErr", ppRefittedTrackDxyBSTrajErr, "");
    patRefittedTracksTab->addColumn<float>("dxyBSSigned", ppRefittedTrackDxyBSSigned, "");
    patRefittedTracksTab->addColumn<float>("dxyBSSignedErr", ppRefittedTrackDxyBSSignedErr, "");
    patRefittedTracksTab->addColumn<float>("ip3DBSSigned", ppRefittedTrackIp3DBSSigned, "");
    patRefittedTracksTab->addColumn<float>("ip3DBSSignedErr", ppRefittedTrackIp3DBSSignedErr, "");

    patdsaRefittedTracksTab->addColumn<float>("idx", pdRefittedTrackIdx, "");
    patdsaRefittedTracksTab->addColumn<float>("isDSAMuon", pdRefittedTrackIsDSAMuon, "");
    patdsaRefittedTracksTab->addColumn<float>("originalMuonIdx", pdRefittedTrackOriginalIdx, "");
    patdsaRefittedTracksTab->addColumn<float>("pt", pdRefittedTrackPt, "");
    patdsaRefittedTracksTab->addColumn<float>("ptErr", pdRefittedTrackPtErr, "");
    patdsaRefittedTracksTab->addColumn<float>("px", pdRefittedTrackPx, "");
    patdsaRefittedTracksTab->addColumn<float>("py", pdRefittedTrackPy, "");
    patdsaRefittedTracksTab->addColumn<float>("pz", pdRefittedTrackPz, "");
    patdsaRefittedTracksTab->addColumn<float>("eta", pdRefittedTrackEta, "");
    patdsaRefittedTracksTab->addColumn<float>("etaErr", pdRefittedTrackEtaErr, "");
    patdsaRefittedTracksTab->addColumn<float>("phi", pdRefittedTrackPhi, "");
    patdsaRefittedTracksTab->addColumn<float>("phiErr", pdRefittedTrackPhiErr, "");
    patdsaRefittedTracksTab->addColumn<float>("charge", pdRefittedTrackCharge, "");
    patdsaRefittedTracksTab->addColumn<float>("normChi2", pdRefittedTrackNormChi2, "");
    patdsaRefittedTracksTab->addColumn<float>("ndof", pdRefittedTrackNdof, "");
    patdsaRefittedTracksTab->addColumn<float>("chi2", pdRefittedTrackChi2, "");
    patdsaRefittedTracksTab->addColumn<float>("dzPV", pdRefittedTrackDzPV, "");
    patdsaRefittedTracksTab->addColumn<float>("dzPVErr", pdRefittedTrackDzPVErr, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyPVTraj", pdRefittedTrackDxyPVTraj, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyPVTrajErr", pdRefittedTrackDxyPVTrajErr, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyPVSigned", pdRefittedTrackDxyPVSigned, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyPVSignedErr", pdRefittedTrackDxyPVSignedErr, "");
    patdsaRefittedTracksTab->addColumn<float>("ip3DPVSigned", pdRefittedTrackIp3DPVSigned, "");
    patdsaRefittedTracksTab->addColumn<float>("ip3DPVSignedErr", pdRefittedTrackIp3DPVSignedErr, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyBS", pdRefittedTrackDxyBS, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyBSErr", pdRefittedTrackDxyBSErr, "");
    patdsaRefittedTracksTab->addColumn<float>("dzBS", pdRefittedTrackDzBS, "");
    patdsaRefittedTracksTab->addColumn<float>("dzBSErr", pdRefittedTrackDzBSErr, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyBSTraj", pdRefittedTrackDxyBSTraj, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyBSTrajErr", pdRefittedTrackDxyBSTrajErr, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyBSSigned", pdRefittedTrackDxyBSSigned, "");
    patdsaRefittedTracksTab->addColumn<float>("dxyBSSignedErr", pdRefittedTrackDxyBSSignedErr, "");
    patdsaRefittedTracksTab->addColumn<float>("ip3DBSSigned", pdRefittedTrackIp3DBSSigned, "");
    patdsaRefittedTracksTab->addColumn<float>("ip3DBSSignedErr", pdRefittedTrackIp3DBSSignedErr, "");

    dsaRefittedTracksTab->addColumn<float>("idx", ddRefittedTrackIdx, "");
    dsaRefittedTracksTab->addColumn<float>("isDSAMuon", ddRefittedTrackIsDSAMuon, "");
    dsaRefittedTracksTab->addColumn<float>("originalMuonIdx", ddRefittedTrackOriginalIdx, "");
    dsaRefittedTracksTab->addColumn<float>("pt", ddRefittedTrackPt, "");
    dsaRefittedTracksTab->addColumn<float>("ptErr", ddRefittedTrackPtErr, "");
    dsaRefittedTracksTab->addColumn<float>("px", ddRefittedTrackPx, "");
    dsaRefittedTracksTab->addColumn<float>("py", ddRefittedTrackPy, "");
    dsaRefittedTracksTab->addColumn<float>("pz", ddRefittedTrackPz, "");
    dsaRefittedTracksTab->addColumn<float>("eta", ddRefittedTrackEta, "");
    dsaRefittedTracksTab->addColumn<float>("etaErr", ddRefittedTrackEtaErr, "");
    dsaRefittedTracksTab->addColumn<float>("phi", ddRefittedTrackPhi, "");
    dsaRefittedTracksTab->addColumn<float>("phiErr", ddRefittedTrackPhiErr, "");
    dsaRefittedTracksTab->addColumn<float>("charge", ddRefittedTrackCharge, "");
    dsaRefittedTracksTab->addColumn<float>("normChi2", ddRefittedTrackNormChi2, "");
    dsaRefittedTracksTab->addColumn<float>("ndof", ddRefittedTrackNdof, "");
    dsaRefittedTracksTab->addColumn<float>("chi2", ddRefittedTrackChi2, "");
    dsaRefittedTracksTab->addColumn<float>("dzPV", ddRefittedTrackDzPV, "");
    dsaRefittedTracksTab->addColumn<float>("dzPVErr", ddRefittedTrackDzPVErr, "");
    dsaRefittedTracksTab->addColumn<float>("dxyPVTraj", ddRefittedTrackDxyPVTraj, "");
    dsaRefittedTracksTab->addColumn<float>("dxyPVTrajErr", ddRefittedTrackDxyPVTrajErr, "");
    dsaRefittedTracksTab->addColumn<float>("dxyPVSigned", ddRefittedTrackDxyPVSigned, "");
    dsaRefittedTracksTab->addColumn<float>("dxyPVSignedErr", ddRefittedTrackDxyPVSignedErr, "");
    dsaRefittedTracksTab->addColumn<float>("ip3DPVSigned", ddRefittedTrackIp3DPVSigned, "");
    dsaRefittedTracksTab->addColumn<float>("ip3DPVSignedErr", ddRefittedTrackIp3DPVSignedErr, "");
    dsaRefittedTracksTab->addColumn<float>("dxyBS", ddRefittedTrackDxyBS, "");
    dsaRefittedTracksTab->addColumn<float>("dxyBSErr", ddRefittedTrackDxyBSErr, "");
    dsaRefittedTracksTab->addColumn<float>("dzBS", ddRefittedTrackDzBS, "");
    dsaRefittedTracksTab->addColumn<float>("dzBSErr", ddRefittedTrackDzBSErr, "");
    dsaRefittedTracksTab->addColumn<float>("dxyBSTraj", ddRefittedTrackDxyBSTraj, "");
    dsaRefittedTracksTab->addColumn<float>("dxyBSTrajErr", ddRefittedTrackDxyBSTrajErr, "");
    dsaRefittedTracksTab->addColumn<float>("dxyBSSigned", ddRefittedTrackDxyBSSigned, "");
    dsaRefittedTracksTab->addColumn<float>("dxyBSSignedErr", ddRefittedTrackDxyBSSignedErr, "");
    dsaRefittedTracksTab->addColumn<float>("ip3DBSSigned", ddRefittedTrackIp3DBSSigned, "");
    dsaRefittedTracksTab->addColumn<float>("ip3DBSSignedErr", ddRefittedTrackIp3DBSSignedErr, "");

    iEvent.put(std::move(patRefittedTracksTab), "PatMuonVertexRefittedTracks");
    iEvent.put(std::move(patdsaRefittedTracksTab), "PatDSAMuonVertexRefittedTracks");
    iEvent.put(std::move(dsaRefittedTracksTab), "DSAMuonVertexRefittedTracks");
  }
  
}

std::pair<float,float> MuonVertexTableProducer::getVxy(const reco::Vertex muonVertex) const {
  float vxy = sqrt(muonVertex.x()*muonVertex.x() + muonVertex.y()*muonVertex.y());
  float vxySigma = (1/vxy)*sqrt(muonVertex.x()*muonVertex.x()*muonVertex.xError()*muonVertex.xError() + muonVertex.y()*muonVertex.y()*muonVertex.yError()*muonVertex.yError());
  return std::make_pair(vxy,vxySigma);
}

std::pair<float,float> MuonVertexTableProducer::getVxyz(const reco::Vertex muonVertex) const {
  float vxyz = sqrt(muonVertex.x()*muonVertex.x() + muonVertex.y()*muonVertex.y() + muonVertex.z()*muonVertex.z());
  float vxyzSigma = (1/vxyz)*sqrt(muonVertex.x()*muonVertex.x()*muonVertex.xError()*muonVertex.xError() + muonVertex.y()*muonVertex.y()*muonVertex.yError()*muonVertex.yError() + muonVertex.z()*muonVertex.z()*muonVertex.zError()*muonVertex.zError());
  return std::make_pair(vxyz,vxyzSigma);
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

/**
*  Proximity match based on EXO-23-010
*  Calculating deltaR between the innermost hit of the DSAMuon trackRef
*  and the extracted closest position of PATMuon track
**/ 
template <typename MuonType>
float MuonVertexTableProducer::getProximityDeltaR(const MuonType& track, 
                        const MuonType& trackRef,
                        const edm::ESHandle<MagneticField>& magneticField,
                        const edm::ESHandle<Propagator>& propagator) const
{
  FreeTrajectoryState trajState(GlobalPoint(track.vx(), track.vy(), track.vz()),
      GlobalVector(track.px(), track.py(), track.pz()),
      track.charge(), magneticField.product());

  GlobalPoint refPos(trackRef.innerPosition().x(),
      trackRef.innerPosition().y(),
      trackRef.innerPosition().z());
  FreeTrajectoryState trajStatePCA(propagator->propagate(trajState, refPos));

  float dR = deltaR(trajStatePCA.position().eta(), trajStatePCA.position().phi(),
                    trackRef.innerPosition().eta(), trackRef.innerPosition().phi());
  return dR;
}

/**
*  Proximity between the muons based on EXO-23-010
*  Getting Distance of Closest Approach between muon tracks using TwoTrackMinimumDistance
*  Returns tuple of distance (float), error of distance (float) and crossing point (GlobalPoint)
**/ 
std::tuple<float, float, GlobalPoint> MuonVertexTableProducer::getDistanceBetweenMuonTracks(const reco::Track& track1, 
                                                            const reco::Track& track2,
                                                            const edm::ESHandle<MagneticField>& magneticField) const 
{
  TwoTrackMinimumDistance ttmd;
  FreeTrajectoryState fts1(GlobalPoint(track1.vx(), track1.vy(), track1.vz()),
                            GlobalVector(track1.px(), track1.py(), track1.pz()),
                            track1.charge(), magneticField.product());
  FreeTrajectoryState fts2(GlobalPoint(track2.vx(), track2.vy(), track2.vz()),
                            GlobalVector(track2.px(), track2.py(), track2.pz()),
                            track2.charge(), magneticField.product());
  bool status = ttmd.calculate(fts1, fts2);
  if (!status) return std::tuple(-999.f, status, GlobalPoint(-999.f,-999.f,-999.f));
  return std::make_tuple(ttmd.distance(), status, ttmd.crossingPoint());
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonVertexTableProducer);
