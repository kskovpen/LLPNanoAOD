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

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


class DGLMuonTableProducer : public edm::global::EDProducer<> {
  public:
    DGLMuonTableProducer(const edm::ParameterSet &iConfig)
      :
      dglMuonTag_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dglMuons"))),
      muonTag_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      vtxTag_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertex"))),
      bsTag_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot")))
      {
      produces<nanoaod::FlatTable>("DGLMuon");
    }

    ~DGLMuonTableProducer() override {}

    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("dglMuons")->setComment("input displaced global muon collection");
      desc.add<edm::InputTag>("muons")->setComment("input muon collection");
      desc.add<edm::InputTag>("primaryVertex")->setComment("input primary vertex collection");
      desc.add<edm::InputTag>("beamspot")->setComment("input beamspot collection");
      descriptions.add("DGLMuonTable", desc);
    }

  private:
    void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

    bool passesDisplacedID(const reco::Track& DGLMuon) const;
    int getMatches(const pat::Muon& muon, const reco::Track& DGLMuon, const float minPositionDiff) const;

    edm::EDGetTokenT<std::vector<reco::Track>> dglMuonTag_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonTag_;
    edm::EDGetTokenT<reco::VertexCollection> vtxTag_;
    edm::EDGetTokenT<reco::BeamSpot> bsTag_;

};


void DGLMuonTableProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  
  float minPositionDiffForMatching = 1e-6;

  edm::Handle<std::vector<reco::Track>> DGLMuonHandle;
  iEvent.getByToken(dglMuonTag_, DGLMuonHandle);
  edm::Handle<std::vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonTag_, muonHandle);

  edm::Handle<reco::VertexCollection> primaryVertexHandle;
  iEvent.getByToken(vtxTag_, primaryVertexHandle);
  const auto& pv = primaryVertexHandle->at(0);
  GlobalPoint primaryVertex(pv.x(), pv.y(), pv.z());

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(bsTag_, beamSpotHandle);
  const auto& bs = beamSpotHandle->position();
  GlobalPoint beamSpot(bs.x(), bs.y(), bs.z());
  reco::Vertex beamSpotVertex(beamSpotHandle->position(), beamSpotHandle->covariance3D());

  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  unsigned int nDGLMuons = DGLMuonHandle->size();
  unsigned int nMuons = muonHandle->size();

  std::vector<float> idx,pt,ptErr,eta,etaErr,phi,phiErr,charge,dxy,dz,vx,vy,vz,chi2,ndof;

  std::vector<float> trkNumPlanes,trkNumHits,trkNumDTHits,trkNumCSCHits,normChi2,outerEta, outerPhi;

  std::vector<float> dzPV,dzPVErr,dxyPVTraj,dxyPVTrajErr,dxyPVSigned,dxyPVSignedErr,ip3DPVSigned,ip3DPVSignedErr;
  std::vector<float> dxyBS,dxyBSErr,dzBS,dzBSErr,dxyBSTraj,dxyBSTrajErr,dxyBSSigned,dxyBSSignedErr,ip3DBSSigned,ip3DBSSignedErr;

  std::vector<float> displacedId;
  std::vector<std::vector<float>> nMatchesPerMuon;
  std::vector<float> muonMatch1,muonMatch1idx,muonMatch2,muonMatch2idx,muonMatch3,muonMatch3idx,muonMatch4,muonMatch4idx,muonMatch5,muonMatch5idx;

  for(unsigned int i = 0; i < nDGLMuons; i++) {

    const reco::Track & DGLMuon = (*DGLMuonHandle)[i];
    idx.push_back(i);

    pt.push_back(DGLMuon.pt());
    ptErr.push_back(DGLMuon.ptError());
    eta.push_back(DGLMuon.eta());
    etaErr.push_back(DGLMuon.etaError());
    phi.push_back(DGLMuon.phi());
    phiErr.push_back(DGLMuon.phiError());
    charge.push_back(DGLMuon.charge());
    dxy.push_back(DGLMuon.dxy());
    dz.push_back(DGLMuon.dz());
    vx.push_back(DGLMuon.vx());
    vy.push_back(DGLMuon.vy());
    vz.push_back(DGLMuon.vz());
    chi2.push_back(DGLMuon.chi2());
    ndof.push_back(DGLMuon.ndof());

    trkNumPlanes.push_back(DGLMuon.hitPattern().muonStationsWithValidHits());
    trkNumHits.push_back(DGLMuon.hitPattern().numberOfValidMuonHits());
    trkNumDTHits.push_back(DGLMuon.hitPattern().numberOfValidMuonDTHits());
    trkNumCSCHits.push_back(DGLMuon.hitPattern().numberOfValidMuonCSCHits());
    normChi2.push_back(DGLMuon.normalizedChi2());

    outerEta.push_back(DGLMuon.outerEta());
    outerPhi.push_back(DGLMuon.outerPhi());

    dzPV.push_back(DGLMuon.dz(pv.position()));
    dzPVErr.push_back(std::hypot(DGLMuon.dzError(), pv.zError()));
    reco::TransientTrack transientTrack = builder->build(DGLMuon);
    TrajectoryStateClosestToPoint trajectoryPV = transientTrack.trajectoryStateClosestToPoint(primaryVertex);
    dxyPVTraj.push_back(trajectoryPV.perigeeParameters().transverseImpactParameter());
    dxyPVTrajErr.push_back(trajectoryPV.perigeeError().transverseImpactParameterError());
    GlobalVector muonRefTrackDir(DGLMuon.px(),DGLMuon.py(),DGLMuon.pz());
    dxyPVSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, pv).second.value());
    dxyPVSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, pv).second.error());

    ip3DPVSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, pv).second.value());
    ip3DPVSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, pv).second.error());  

    dxyBS.push_back(DGLMuon.dxy(bs));
    dxyBSErr.push_back(std::hypot(DGLMuon.dxyError(), beamSpotVertex.zError()));
    dzBS.push_back(DGLMuon.dz(bs));
    dzBSErr.push_back(std::hypot(DGLMuon.dzError(), beamSpotVertex.zError()));
    TrajectoryStateClosestToBeamLine trajectoryBS = transientTrack.stateAtBeamLine();
    dxyBSTraj.push_back(trajectoryBS.transverseImpactParameter().value());
    dxyBSTrajErr.push_back(trajectoryBS.transverseImpactParameter().error());
    dxyBSSigned.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, beamSpotVertex).second.value());
    dxyBSSignedErr.push_back(IPTools::signedTransverseImpactParameter(transientTrack, muonRefTrackDir, beamSpotVertex).second.error());  

    ip3DBSSigned.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.value());
    ip3DBSSignedErr.push_back(IPTools::signedImpactParameter3D(transientTrack, muonRefTrackDir, beamSpotVertex).second.error());

    float passesDisplacedId = 0;
    if(passesDisplacedID(DGLMuon)) passesDisplacedId=1;
    displacedId.push_back(passesDisplacedId);

    // Assigning 5 best matches and corresponding muon indices
    std::vector<std::pair<float, float>> muonMatches(5, std::make_pair(-1.0,-1.0));
    std::vector<float> nMuonMatches;
    for (unsigned int j = 0; j < nMuons; j++){
      if (j > 4) break;
      const pat::Muon & muon = (*muonHandle)[j];
      // Muon-DSA Matches Table
      int nMatches = getMatches(muon, DGLMuon, minPositionDiffForMatching);
      muonMatches[j] = std::make_pair(nMatches, j);
      nMuonMatches.push_back(nMatches);
    }
    nMatchesPerMuon.push_back(nMuonMatches);
    std::sort(muonMatches.rbegin(), muonMatches.rend());
    muonMatch1.push_back(muonMatches[0].first);
    muonMatch1idx.push_back(muonMatches[0].second);
    muonMatch2.push_back(muonMatches[1].first);
    muonMatch2idx.push_back(muonMatches[1].second);
    muonMatch3.push_back(muonMatches[2].first);
    muonMatch3idx.push_back(muonMatches[2].second);
    muonMatch4.push_back(muonMatches[3].first);
    muonMatch4idx.push_back(muonMatches[3].second);
    muonMatch5.push_back(muonMatches[4].first);
    muonMatch5idx.push_back(muonMatches[4].second);

  }

  auto DGLMuonTab = std::make_unique<nanoaod::FlatTable>(DGLMuonHandle->size(), "DGLMuon", false, false);

  DGLMuonTab->addColumn<float>("idx", idx, "",  nanoaod::FlatTable::FloatColumn);

  DGLMuonTab->addColumn<float>("pt", pt, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("ptErr", ptErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("eta", eta, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("etaErr", etaErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("phi", phi, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("phiErr", phiErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("charge", charge, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxy", dxy, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dz", dz, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("vx", vx, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("vy", vy, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("vz", vz, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("chi2", chi2, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("ndof", ndof, "",  nanoaod::FlatTable::FloatColumn);

  DGLMuonTab->addColumn<float>("trkNumPlanes", trkNumPlanes, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("trkNumHits", trkNumHits, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("trkNumDTHits", trkNumDTHits, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("trkNumCSCHits", trkNumCSCHits, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("normChi2", normChi2, "",  nanoaod::FlatTable::FloatColumn);

  DGLMuonTab->addColumn<float>("outerEta", outerEta, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("outerPhi", outerPhi, "",  nanoaod::FlatTable::FloatColumn);
  
  DGLMuonTab->addColumn<float>("dzPV", dzPV, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dzPVErr", dzPVErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyPVTraj", dxyPVTraj, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyPVTrajErr", dxyPVTrajErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyPVSigned", dxyPVSigned, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyPVSignedErr", dxyPVSignedErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("ip3DPVSigned", ip3DPVSigned, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("ip3DPVSignedErr", ip3DPVSignedErr, "",  nanoaod::FlatTable::FloatColumn);

  DGLMuonTab->addColumn<float>("dxyBS", dxyBS, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyBSErr", dxyBSErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dzBS", dzBS, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dzBSErr", dzBSErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyBSTraj", dxyBSTraj, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyBSTrajErr", dxyBSTrajErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyBSSigned", dxyBSSigned, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("dxyBSSignedErr", dxyBSSignedErr, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("ip3DBSSigned", ip3DBSSigned, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("ip3DBSSignedErr", ip3DBSSignedErr, "",  nanoaod::FlatTable::FloatColumn);

  DGLMuonTab->addColumn<float>("displacedID", displacedId, "",  nanoaod::FlatTable::FloatColumn);

  DGLMuonTab->addColumn<float>("muonMatch1", muonMatch1, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch1idx", muonMatch1idx, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch2", muonMatch2, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch2idx", muonMatch2idx, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch3", muonMatch3, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch3idx", muonMatch3idx, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch4", muonMatch4, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch4idx", muonMatch4idx, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch5", muonMatch5, "",  nanoaod::FlatTable::FloatColumn);
  DGLMuonTab->addColumn<float>("muonMatch5idx", muonMatch5idx, "",  nanoaod::FlatTable::FloatColumn);

  iEvent.put(std::move(DGLMuonTab), "DGLMuon");
}

bool DGLMuonTableProducer::passesDisplacedID(const reco::Track& DGLMuon) const {
  // displaced muon Id as recommended by Muon POG
  float validHits =  DGLMuon.hitPattern().numberOfValidMuonCSCHits() + DGLMuon.hitPattern().numberOfValidMuonDTHits();
  if(validHits > 12){
    if(DGLMuon.hitPattern().numberOfValidMuonCSCHits() != 0 || (DGLMuon.hitPattern().numberOfValidMuonCSCHits() == 0 && DGLMuon.hitPattern().numberOfValidMuonDTHits() > 18)){
      if(DGLMuon.normalizedChi2() < 2.5) {
        if(DGLMuon.ptError()/DGLMuon.pt() < 1){
          return true;
        }
      }
    }
  }
  return false;
}

int DGLMuonTableProducer::getMatches(const pat::Muon& muon, const reco::Track& DGLMuon, const float minPositionDiff=1e-6) const {

  int nMatches = 0;
  if (!(muon.isTrackerMuon() && muon::isGoodMuon(muon, muon::TrackerMuonArbitrated))) return -1;

  for (auto& hit : DGLMuon.recHits()){

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
DEFINE_FWK_MODULE(DGLMuonTableProducer);
