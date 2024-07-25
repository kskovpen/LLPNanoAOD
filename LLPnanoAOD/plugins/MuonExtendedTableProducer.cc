// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


class MuonExtendedTableProducer : public edm::global::EDProducer<> {
  public:
    explicit MuonExtendedTableProducer(const edm::ParameterSet &iConfig) :
      name_(iConfig.getParameter<std::string>("name")),
      muonTag_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      dsaMuonTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dsaMuons"))),
      vtxTag_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertex"))),
      bsTag_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
      generalTrackTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("generalTracks"))),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder")))
    {
      produces<nanoaod::FlatTable>();
    }

    ~MuonExtendedTableProducer() override {};

    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("muons")->setComment("input muon collection");
      desc.add<edm::InputTag>("dsaMuons")->setComment("input displaced standalone muon collection");
      desc.add<edm::InputTag>("primaryVertex")->setComment("input primary vertex collection");
      desc.add<edm::InputTag>("beamspot")->setComment("input beamspot collection");
      desc.add<edm::InputTag>("generalTracks")->setComment("input generalTracks collection");
      desc.add<std::string>("name")->setComment("name of the muon nanoaod::FlatTable we are extending");
      descriptions.add("muonTable", desc);
    }

  private:
    void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

    int getMatches(const pat::Muon& muon, const reco::Track& dsaMuon, const float minPositionDiff) const;

    std::string name_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonTag_;
    edm::EDGetTokenT<std::vector<reco::Track>> dsaMuonTag_;
    edm::EDGetTokenT<reco::VertexCollection> vtxTag_;
    edm::EDGetTokenT<reco::BeamSpot> bsTag_;
    edm::EDGetTokenT<std::vector<reco::Track>> generalTrackTag_;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;

};

void MuonExtendedTableProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const 
{

  float minPositionDiffForMatching = 1e-6;

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonTag_, muons);
  edm::Handle<std::vector<reco::Track>> dsaMuons;
  iEvent.getByToken(dsaMuonTag_, dsaMuons);
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(vtxTag_, primaryVertices);
  edm::Handle<reco::BeamSpot> beamspots;
  iEvent.getByToken(bsTag_, beamspots);
  edm::Handle<std::vector<reco::Track>> generalTracks;
  iEvent.getByToken(generalTrackTag_, generalTracks);

  const auto& pv = primaryVertices->at(0);
  GlobalPoint primaryVertex(pv.x(), pv.y(), pv.z());

  const auto& bs = beamspots->position();
  GlobalPoint beamSpot(bs.x(), bs.y(), bs.z());
  reco::Vertex beamSpotVertex(beamspots->position(), beamspots->covariance3D());

  // edm::ESHandle<TransientTrackBuilder> builder;
  // iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  // auto const& builder = &iSetup.getData(transientTrackBuilderToken_);
  edm::ESHandle<TransientTrackBuilder> builder = iSetup.getHandle(transientTrackBuilderToken_);

  unsigned int nMuons = muons->size();
  unsigned int nDsaMuons = dsaMuons->size();
  

  std::vector<float> idx, charge, trkPt, trkPtErr;

  std::vector<float> dzPV,dzPVErr,dxyPVTraj,dxyPVTrajErr,dxyPVSigned,dxyPVSignedErr,ip3DPVSigned,ip3DPVSignedErr;
  std::vector<float> dxyBS,dxyBSErr,dzBS,dzBSErr,dxyBSTraj,dxyBSTrajErr,dxyBSSigned,dxyBSSignedErr,ip3DBSSigned,ip3DBSSignedErr;

  std::vector<float> trkNumPlanes,trkNumHits,trkNumDTHits,trkNumCSCHits,trkNumPixelHits(nMuons,-1),trkNumTrkLayers(nMuons,-1),normChi2;
  std::vector<float> outerEta(nMuons,-5),outerPhi(nMuons,-5);
  std::vector<float> innerVx(nMuons,-1),innerVy(nMuons,-1),innerVz(nMuons,-1),innerPt(nMuons,-1),innerEta(nMuons,-5),innerPhi(nMuons,-5);

  std::vector<std::vector<float>> nMatchesPerDSA;
  std::vector<float> dsaMatch1,dsaMatch1idx,dsaMatch2,dsaMatch2idx,dsaMatch3,dsaMatch3idx,dsaMatch4,dsaMatch4idx,dsaMatch5,dsaMatch5idx;

  // std::vector<float> displacedTrackIso03, displacedTrackIso04;

  for (unsigned int i = 0; i < nMuons; i++) {
    const pat::Muon & muon = (*muons)[i];
    const pat::MuonRef muonRef(muons,i);

    // const auto& track = muon.bestTrack();
    // const reco::Track* track; // = muon.tunePMuonBestTrack();
    reco::TrackRef trackRef; // = muonRef->combinedMuon();

    if(muon.isGlobalMuon()) {
      trackRef = muon.combinedMuon();
    }
    else if (muon.isStandAloneMuon()) {
      trackRef = muon.standAloneMuon();
    }
    else {
      trackRef = muon.tunePMuonBestTrack();
    }

    // if (!trackRef.isNonnull()) continue;
    idx.push_back(i);

    const auto& track = trackRef.get();
    reco::TransientTrack transientTrack = builder->build(track);

    charge.push_back(muon.charge());

    trkPt.push_back(track->pt());
    trkPtErr.push_back(track->ptError());

    dzPV.push_back(track->dz(pv.position()));
    dzPVErr.push_back(std::hypot(track->dzError(), pv.zError()));
    TrajectoryStateClosestToPoint trajectoryPV = transientTrack.trajectoryStateClosestToPoint(primaryVertex);
    dxyPVTraj.push_back(trajectoryPV.perigeeParameters().transverseImpactParameter());
    dxyPVTrajErr.push_back(trajectoryPV.perigeeError().transverseImpactParameterError());
    GlobalVector muonRefTrackDir(muon.px(),muon.py(),muon.pz());
    dxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, pv).second.value());
    dxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, pv).second.error());
    
    ip3DPVSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.value());
    ip3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.error());  

    dxyBS.push_back(track->dxy(bs));
    dxyBSErr.push_back(std::hypot(track->dxyError(), beamSpotVertex.zError()));
    dzBS.push_back(track->dz(bs));
    dzBSErr.push_back(std::hypot(track->dzError(), beamSpotVertex.zError()));
    TrajectoryStateClosestToBeamLine trajectoryBS = transientTrack.stateAtBeamLine();
    dxyBSTraj.push_back(trajectoryBS.transverseImpactParameter().value());
    dxyBSTrajErr.push_back(trajectoryBS.transverseImpactParameter().error());
    dxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, beamSpotVertex).second.value());
    dxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, beamSpotVertex).second.error());  

    ip3DBSSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.value());
    ip3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.error()); 

    // if (!trackRef.isNonnull()) continue;
    trkNumPlanes.push_back(track->hitPattern().muonStationsWithValidHits());
    trkNumHits.push_back(track->hitPattern().numberOfValidMuonHits());
    trkNumDTHits.push_back(track->hitPattern().numberOfValidMuonDTHits());
    trkNumCSCHits.push_back(track->hitPattern().numberOfValidMuonCSCHits());

    normChi2.push_back(track->normalizedChi2());

    // Cannot get outer track for tracker muons
    if (track->extra().isNonnull() && track->extra().isAvailable() && track->outerOk()) {
      outerEta[i] = track->outerEta();
      outerPhi[i] = track->outerPhi();
    }

    if(muon.innerTrack().isNonnull() && muon.innerTrack().isAvailable()){
      innerVx[i] = muon.innerTrack()->vx();
      innerVy[i] = muon.innerTrack()->vy();
      innerVz[i] = muon.innerTrack()->vz();
      innerPt[i] = muon.innerTrack()->pt();
      innerEta[i] = muon.innerTrack()->eta();
      innerPhi[i] = muon.innerTrack()->phi();
      trkNumPixelHits[i] = muon.innerTrack()->hitPattern().numberOfValidPixelHits();
      trkNumTrkLayers[i] = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement();
    }

    std::vector<std::pair<float, float>> dsaMatches(5, std::make_pair(-1.0,-1.0));
    std::vector<float> nDsaMatches;
    for (unsigned int j = 0; j < nDsaMuons; j++){
      if (j > 4) break;
      const reco::Track & dsaMuon = (*dsaMuons)[j];
      // Muon-DSA Matches Table
      int nMatches = getMatches(muon, dsaMuon, minPositionDiffForMatching);
      dsaMatches[j] = std::make_pair(nMatches, j);
      nDsaMatches.push_back(nMatches);
    }
    nMatchesPerDSA.push_back(nDsaMatches);
    std::sort(dsaMatches.rbegin(), dsaMatches.rend());
    dsaMatch1.push_back(dsaMatches[0].first);
    dsaMatch1idx.push_back(dsaMatches[0].second);
    dsaMatch2.push_back(dsaMatches[1].first);
    dsaMatch2idx.push_back(dsaMatches[1].second);
    dsaMatch3.push_back(dsaMatches[2].first);
    dsaMatch3idx.push_back(dsaMatches[2].second);
    dsaMatch4.push_back(dsaMatches[3].first);
    dsaMatch4idx.push_back(dsaMatches[3].second);
    dsaMatch5.push_back(dsaMatches[4].first);
    dsaMatch5idx.push_back(dsaMatches[4].second);
  }

  auto tab  = std::make_unique<nanoaod::FlatTable>(nMuons, name_, false, true);
  tab->addColumn<float>("idx", idx, "LLPnanoAOD muon index");

  tab->addColumn<float>("trkPt", trkPt, "");
  tab->addColumn<float>("trkPtErr", trkPtErr, "");

  tab->addColumn<float>("dzPV", dzPV, "");
  tab->addColumn<float>("dzPVErr", dzPVErr, "");
  tab->addColumn<float>("dxyPVTraj", dxyPVTraj, "");
  tab->addColumn<float>("dxyPVTrajErr", dxyPVTrajErr, "");
  tab->addColumn<float>("dxyPVSigned", dxyPVSigned, "");
  tab->addColumn<float>("dxyPVSignedErr", dxyPVSignedErr, "");
  tab->addColumn<float>("ip3DPVSigned", ip3DPVSigned, "");
  tab->addColumn<float>("ip3DPVSignedErr", ip3DPVSignedErr, "");

  tab->addColumn<float>("dxyBS", dxyBS, "");
  tab->addColumn<float>("dxyBSErr", dxyBSErr, "");
  tab->addColumn<float>("dzBS", dzBS, "");
  tab->addColumn<float>("dzBSErr", dzBSErr, "");
  tab->addColumn<float>("dxyBSTraj", dxyBSTraj, "");
  tab->addColumn<float>("dxyBSTrajErr", dxyBSTrajErr, "");
  tab->addColumn<float>("dxyBSSigned", dxyBSSigned, "");
  tab->addColumn<float>("dxyBSSignedErr", dxyBSSignedErr, "");
  tab->addColumn<float>("ip3DBSSigned", ip3DBSSigned, "");
  tab->addColumn<float>("ip3DBSSignedErr", ip3DBSSignedErr, "");

  tab->addColumn<float>("trkNumPlanes", trkNumPlanes, "");
  tab->addColumn<float>("trkNumHits", trkNumHits, "");
  tab->addColumn<float>("trkNumDTHits", trkNumDTHits, "");
  tab->addColumn<float>("trkNumCSCHits", trkNumCSCHits, "");
  tab->addColumn<float>("normChi2", normChi2, "");
  tab->addColumn<float>("trkNumPixelHits", trkNumPixelHits, "");
  tab->addColumn<float>("trkNumTrkLayers", trkNumTrkLayers, "");

  tab->addColumn<float>("outerEta", outerEta, "");
  tab->addColumn<float>("outerPhi", outerPhi, "");

  tab->addColumn<float>("innerVx", innerVx, "");
  tab->addColumn<float>("innerVy", innerVy, "");
  tab->addColumn<float>("innerVz", innerVz, "");
  tab->addColumn<float>("innerPt", innerPt, "");
  tab->addColumn<float>("innerEta", innerEta, "");
  tab->addColumn<float>("innerPhi", innerPhi, "");

  tab->addColumn<float>("dsaMatch1", dsaMatch1, "");
  tab->addColumn<float>("dsaMatch1idx", dsaMatch1idx, "");
  tab->addColumn<float>("dsaMatch2", dsaMatch2, "");
  tab->addColumn<float>("dsaMatch2idx", dsaMatch2idx, "");
  tab->addColumn<float>("dsaMatch3", dsaMatch3, "");
  tab->addColumn<float>("dsaMatch3idx", dsaMatch3idx, "");
  tab->addColumn<float>("dsaMatch4", dsaMatch4, "");
  tab->addColumn<float>("dsaMatch4idx", dsaMatch4idx, "");
  tab->addColumn<float>("dsaMatch5", dsaMatch5, "");
  tab->addColumn<float>("dsaMatch5idx", dsaMatch5idx, "");

  iEvent.put(std::move(tab));
}

int MuonExtendedTableProducer::getMatches(const pat::Muon& muon, const reco::Track& dsaMuon, const float minPositionDiff=1e-6) const {

  int nMatches = 0;

  for (auto& hit : dsaMuon.recHits()){

    if (!hit->isValid()) continue;
    DetId id = hit->geographicalId();
    if (id.det() != DetId::Muon) continue;

    if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC){

      for (auto& chamber : muon.matches()) {

        if (chamber.id.rawId() != id.rawId()) continue;

        for (auto& segment : chamber.segmentMatches) {

          if (fabs(segment.x - hit->localPosition().x()) < minPositionDiff &&
              fabs(segment.y - hit->localPosition().y()) < minPositionDiff) {
              nMatches++;
              break;
          }
        }
      }
    }
  }
  return nMatches;
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(MuonExtendedTableProducer);
