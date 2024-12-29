// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class DispJetTableProducer : public edm::stream::EDProducer<> {
   
  public:

   DispJetTableProducer(const edm::ParameterSet &iConfig)
     :
     rhoTag_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
     electronTag_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
     muonTag_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
     jetTag_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
     vtxTag_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertex"))),
     secVtxTag_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("secondaryVertex"))) {	
	produces<nanoaod::FlatTable>("DispJetElectron");
	produces<nanoaod::FlatTable>("DispJetElectronVtx");
	produces<nanoaod::FlatTable>("DispJetElectronTrk");
	produces<nanoaod::FlatTable>("DispJetMuon");
	produces<nanoaod::FlatTable>("DispJetMuonVtx");
	produces<nanoaod::FlatTable>("DispJetMuonTrk");
     }

   ~DispJetTableProducer() override {}

   static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("rho")->setComment("input rho collection");
      desc.add<edm::InputTag>("electrons")->setComment("input electron collection");
      desc.add<edm::InputTag>("muons")->setComment("input muon collection");
      desc.add<edm::InputTag>("jets")->setComment("input jet collection");
      desc.add<edm::InputTag>("primaryVertex")->setComment("input primary vertex collection");
      desc.add<edm::InputTag>("secondaryVertex")->setComment("input secondary vertex collection");
      descriptions.add("DispJetTable", desc);
   }

 private:

   virtual void produce(edm::Event&, edm::EventSetup const&) override;
   
   float dEtaInSeed(const pat::Electron* el) const;
   float getPFIso(const pat::Muon& muon) const;
   float getPFIso(const pat::Electron& electron) const;
   void fillLeptonJetVariables(const reco::GsfElectron *el, const reco::Muon *mu, edm::Handle< std::vector< pat::Jet > >& jets, const reco::Vertex& vertex, const double rho, const bool oldMatching);
   
   edm::EDGetTokenT<double> rhoTag_;
   edm::EDGetTokenT<std::vector<pat::Electron> > electronTag_;
   edm::EDGetTokenT<std::vector<pat::Muon> > muonTag_;
   edm::EDGetTokenT<std::vector<pat::Jet> > jetTag_;
   edm::EDGetTokenT<reco::VertexCollection> vtxTag_;
   edm::EDGetTokenT<reco::VertexCollection> secVtxTag_;

   std::vector<int> el_idx;
   std::vector<bool> el_lIVF_match;
   std::vector<bool> el_isEB, el_isEE;
   std::vector<float> el_superClusterOverP, el_ecalEnergy, el_dEtaInSeed;
   std::vector<int> el_numberInnerHitsMissing, el_numberOfValidPixelHits, el_numberOfValidTrackerHits;
   std::vector<int> el_IVF_df, el_IVF_ntracks, el_IVF_elid;
   std::vector<float> el_IVF_x, el_IVF_y, el_IVF_z, el_IVF_cx, el_IVF_cy, el_IVF_cz, el_IVF_chi2, el_IVF_pt, el_IVF_eta, el_IVF_phi, el_IVF_E, el_IVF_mass;
   std::vector<int> el_IVF_trackcharge, el_IVF_trackelid, el_IVF_trackvtxid;
   std::vector<float> el_IVF_trackpt, el_IVF_tracketa, el_IVF_trackphi, el_IVF_trackE, el_IVF_trackdxy, el_IVF_trackdz;

   std::vector<float> el_ptRatio, el_ptRel, el_closestJetTag_b, el_closestJetTag_bb, el_closestJetTag_lepb, el_closestJetTag;
   std::vector<int> el_selectedTrackMult;
   std::vector<float> el_relIso0p4;
   std::vector<float> el_sigmaIetaIeta, el_deltaPhiSuperClusterTrack, el_deltaEtaSuperClusterTrack, el_eInvMinusPInv, el_hOverE;
   
   std::vector<float> el_dxy, el_dz, el_3dIP, el_3dIPSig;
   
   std::vector<int> mu_idx;
   std::vector<bool> mu_lIVF_match;
   std::vector<float> mu_innerTrackValidFraction, mu_globalTrackNormalizedChi2, mu_CQChi2Position, mu_CQTrackKink;
   std::vector<int> mu_numberOfMatchedStation, mu_numberOfValidPixelHits, mu_numberOfValidTrackerHits, mu_numberInnerHitsMissing, mu_trackerLayersWithMeasurement, mu_numberInnerHits;
   std::vector<int> mu_IVF_df, mu_IVF_ntracks, mu_IVF_muid;
   std::vector<float> mu_IVF_x, mu_IVF_y, mu_IVF_z, mu_IVF_cx, mu_IVF_cy, mu_IVF_cz, mu_IVF_chi2, mu_IVF_pt, mu_IVF_eta, mu_IVF_phi, mu_IVF_E, mu_IVF_mass;
   std::vector<int> mu_IVF_trackcharge, mu_IVF_trackmuid, mu_IVF_trackvtxid;
   std::vector<float> mu_IVF_trackpt, mu_IVF_tracketa, mu_IVF_trackphi, mu_IVF_trackE, mu_IVF_trackdxy, mu_IVF_trackdz;
   
   std::vector<float> mu_ptRatio, mu_ptRel, mu_closestJetTag_b, mu_closestJetTag_bb, mu_closestJetTag_lepb, mu_closestJetTag;
   std::vector<int> mu_selectedTrackMult;
   std::vector<float> mu_relIso0p4;
   std::vector<float> mu_trackPt, mu_trackPtErr;
   
   std::vector<float> mu_dxy, mu_dz, mu_3dIP, mu_3dIPSig;
};

void DispJetTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   
   float minPositionDiffForMatching = 1e-6;

   edm::Handle<double> rhoHandle;
   iEvent.getByToken(rhoTag_, rhoHandle);
   
   edm::Handle<std::vector<pat::Electron> > electronHandle;
   iEvent.getByToken(electronTag_, electronHandle);
   
   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByToken(muonTag_, muonHandle);
   
   edm::Handle<std::vector<pat::Jet> > jetHandle;
   iEvent.getByToken(jetTag_, jetHandle);
   
   edm::Handle<reco::VertexCollection> primaryVertexHandle;
   iEvent.getByToken(vtxTag_, primaryVertexHandle);
   
   edm::Handle<reco::VertexCollection> secondaryVertexHandle;
   iEvent.getByToken(secVtxTag_, secondaryVertexHandle);
   
   const auto& pv = primaryVertexHandle->at(0);
   GlobalPoint primaryVertex(pv.x(), pv.y(), pv.z());
   
   const reco::VertexCollection & secVertices = (*secondaryVertexHandle);
   
   unsigned int nElectrons = electronHandle->size();
   unsigned int nMuons = muonHandle->size();
   unsigned int nJets = jetHandle->size();
  
   el_idx.clear();
   el_lIVF_match.clear();
   el_isEB.clear(); el_isEE.clear();
   el_superClusterOverP.clear(); el_ecalEnergy.clear(); el_dEtaInSeed.clear();
   el_numberInnerHitsMissing.clear(); el_numberOfValidPixelHits.clear(); el_numberOfValidTrackerHits.clear();
   el_IVF_df.clear(); el_IVF_ntracks.clear(); el_IVF_elid.clear();
   el_IVF_x.clear(); el_IVF_y.clear(); el_IVF_z.clear(); el_IVF_cx.clear(); el_IVF_cy.clear(); el_IVF_cz.clear(); el_IVF_chi2.clear(); el_IVF_pt.clear(); el_IVF_eta.clear(); el_IVF_phi.clear(); el_IVF_E.clear(); el_IVF_mass.clear();
   el_IVF_trackcharge.clear(); el_IVF_trackelid.clear(); el_IVF_trackvtxid.clear();
   el_IVF_trackpt.clear(); el_IVF_tracketa.clear(); el_IVF_trackphi.clear(); el_IVF_trackE.clear(); el_IVF_trackdxy.clear(); el_IVF_trackdz.clear();
   
   el_ptRatio.clear(); el_ptRel.clear(); el_closestJetTag_b.clear(); el_closestJetTag_bb.clear(); el_closestJetTag_lepb.clear(); el_closestJetTag.clear();
   el_selectedTrackMult.clear();
   el_relIso0p4.clear();
   el_sigmaIetaIeta.clear(); el_deltaPhiSuperClusterTrack.clear(); el_deltaEtaSuperClusterTrack.clear(); el_eInvMinusPInv.clear(); el_hOverE.clear();

   el_dxy.clear(); el_dz.clear(); el_3dIP.clear(); el_3dIPSig.clear();
   
   mu_idx.clear();
   mu_lIVF_match.clear();
   mu_innerTrackValidFraction.clear(); mu_globalTrackNormalizedChi2.clear(); mu_CQChi2Position.clear(); mu_CQTrackKink.clear();
   mu_numberOfMatchedStation.clear(); mu_numberOfValidPixelHits.clear(); mu_numberOfValidTrackerHits.clear(); mu_numberInnerHitsMissing.clear(); mu_trackerLayersWithMeasurement.clear(); mu_numberInnerHits.clear();
   mu_IVF_df.clear(); mu_IVF_ntracks.clear(); mu_IVF_muid.clear();
   mu_IVF_x.clear(); mu_IVF_y.clear(); mu_IVF_z.clear(); mu_IVF_cx.clear(); mu_IVF_cy.clear(); mu_IVF_cz.clear(); mu_IVF_chi2.clear(); mu_IVF_pt.clear(); mu_IVF_eta.clear(); mu_IVF_phi.clear(); mu_IVF_E.clear(); mu_IVF_mass.clear();
   mu_IVF_trackcharge.clear(); mu_IVF_trackmuid.clear(); mu_IVF_trackvtxid.clear();
   mu_IVF_trackpt.clear(); mu_IVF_tracketa.clear(); mu_IVF_trackphi.clear(); mu_IVF_trackE.clear(); mu_IVF_trackdxy.clear(); mu_IVF_trackdz.clear();
   
   mu_ptRatio.clear(); mu_ptRel.clear(); mu_closestJetTag_b.clear(); mu_closestJetTag_bb.clear(); mu_closestJetTag_lepb.clear(); mu_closestJetTag.clear();
   mu_selectedTrackMult.clear();
   mu_relIso0p4.clear();
   mu_trackPt.clear(); mu_trackPtErr.clear();
   
   mu_dxy.clear(); mu_dz.clear(); mu_3dIP.clear(); mu_3dIPSig.clear();
   
   int ntrack_max = 100;
   int nElectronsSel = 0;
   int nMuonsSel = 0;
   int nJetsSel = 0;

   for(unsigned int i = 0; i < nJets; i++) {
      
      const pat::Jet & jet = (*jetHandle)[i];
   }   
   
   for(unsigned int i = 0; i < nElectrons; i++) {
      
      const pat::Electron & el = (*electronHandle)[i];
      
      if(el.gsfTrack().isNull()) continue;
      if(el.pt() < 7) continue;     
      if(fabs(el.eta()) > 2.5) continue;
      
      int ielCand = 0;
      for(edm::Ref<pat::PackedCandidateCollection> cand : el.associatedPackedPFCandidates()){	  
	 if( cand.isNonnull() and cand.isAvailable() ) {
	    break;
	 }       
	 ielCand++;
      }
      
      el_idx.push_back(i);
      el_lIVF_match.push_back(false);
      
      el_isEB.push_back(el.isEB());
      el_isEE.push_back(el.isEE());
      el_superClusterOverP.push_back(el.eSuperClusterOverP());
      el_ecalEnergy.push_back(el.ecalEnergy());
      el_dEtaInSeed.push_back(std::abs(dEtaInSeed(&el)));
      el_numberInnerHitsMissing.push_back(el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      el_numberOfValidPixelHits.push_back((!el.gsfTrack().isNull())? el.gsfTrack()->hitPattern().numberOfValidPixelHits() : 0);
      el_numberOfValidTrackerHits.push_back((!el.gsfTrack().isNull())? el.gsfTrack()->hitPattern().numberOfValidTrackerHits() : 0);

      el_relIso0p4.push_back(getPFIso(el));
      el_sigmaIetaIeta.push_back(el.full5x5_sigmaIetaIeta());
      el_deltaPhiSuperClusterTrack.push_back(fabs(el.deltaPhiSuperClusterTrackAtVtx()));
      el_deltaEtaSuperClusterTrack.push_back(fabs(el.deltaEtaSuperClusterTrackAtVtx()));
      el_eInvMinusPInv.push_back((1.0 - el.eSuperClusterOverP())/el.correctedEcalEnergy());
      el_hOverE.push_back(el.hadronicOverEm());
      
      el_dxy.push_back(el.dB(pat::Electron::PV2D));
      el_dz.push_back(el.dB(pat::Electron::PVDZ));
      el_3dIP.push_back(el.dB(pat::Electron::PV3D));
      el_3dIPSig.push_back(fabs(el.dB(pat::Electron::PV3D)/el.edB(pat::Electron::PV3D)));
      
      fillLeptonJetVariables(&el, NULL, jetHandle, pv, *rhoHandle, false);
      
      bool new_vtx = false;
      double dR, deta, normchi2;
      double mindR = 20, minnormchi2 = 10000;
      int nVtx = 0;
      for(const reco::Vertex& vtx : secVertices) {
	 for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++) {
	    reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
	    for(edm::Ref<pat::PackedCandidateCollection> cand : el.associatedPackedPFCandidates()) {	       
	       dR       = reco::deltaR(cand->eta(), cand->phi(), vtxTrack->eta(), vtxTrack->phi());
	       deta     = fabs(cand->eta() - vtxTrack->eta());
	       normchi2 = fabs(vtx.chi2()/vtx.ndof());
	       
	       if((dR < 0.05 or (dR < 0.1 and deta < 0.03)) and (dR < mindR or (dR == mindR and normchi2 < minnormchi2))) {	       
		  new_vtx = true;
		  el_lIVF_match[nElectronsSel] = true;
		  mindR       = dR;
		  minnormchi2 = normchi2;
	       }
	    }
	 }
	 
	 if(new_vtx) {
	    el_IVF_x.push_back(vtx.x());
	    el_IVF_y.push_back(vtx.y());
	    el_IVF_z.push_back(vtx.z());
	    el_IVF_cx.push_back(vtx.xError());
	    el_IVF_cy.push_back(vtx.yError());
	    el_IVF_cz.push_back(vtx.zError());
	    el_IVF_df.push_back(vtx.ndof());
	    el_IVF_chi2.push_back(vtx.chi2());
	    el_IVF_pt.push_back(vtx.p4().pt());
	    el_IVF_eta.push_back(vtx.p4().eta());
	    el_IVF_phi.push_back(vtx.p4().phi());
	    el_IVF_E.push_back(vtx.p4().energy());
	    el_IVF_mass.push_back(vtx.p4().mass());
	    el_IVF_elid.push_back(nElectronsSel);
	    
	    el_IVF_ntracks.push_back(0);
	    for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++) {
	       if(el_IVF_ntracks.back() == ntrack_max) break;
	       reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
	       el_IVF_trackpt.push_back(vtxTrack->pt());
	       el_IVF_tracketa.push_back(vtxTrack->eta());
	       el_IVF_trackphi.push_back(vtxTrack->phi());
	       el_IVF_trackE.push_back(vtxTrack->p());
	       el_IVF_trackcharge.push_back(vtxTrack->charge());
	       el_IVF_trackdxy.push_back(std::abs(vtxTrack->dxy(pv.position())));
	       el_IVF_trackdz.push_back(std::abs(vtxTrack->dz(pv.position())));
	       el_IVF_trackelid.push_back(nElectronsSel);
	       el_IVF_trackvtxid.push_back(nVtx);
	       el_IVF_ntracks.back()++;
	    }
	    nVtx++;
	    new_vtx = false;
	 }
      }
      nElectronsSel += 1;
   }
   
   for(unsigned int i = 0; i < nMuons; i++) {
      
      const pat::Muon & mu = (*muonHandle)[i];
      
      if(mu.innerTrack().isNull()) continue;
      if(mu.pt() < 5) continue;
      if(fabs(mu.eta()) > 2.4) continue;
      if(!mu.isPFMuon()) continue;
      if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue;
      
      mu_idx.push_back(i);
      mu_lIVF_match.push_back(false);
      
      mu_innerTrackValidFraction.push_back((!mu.innerTrack().isNull()) ? mu.innerTrack()->validFraction() : -1);
      mu_globalTrackNormalizedChi2.push_back((!mu.globalTrack().isNull()) ? mu.globalTrack()->normalizedChi2() : -1);
      mu_CQChi2Position.push_back(mu.combinedQuality().chi2LocalPosition);
      mu_CQTrackKink.push_back(mu.combinedQuality().trkKink);
      mu_numberOfMatchedStation.push_back(mu.numberOfMatchedStations());
      mu_numberOfValidPixelHits.push_back((!mu.innerTrack().isNull()) ? mu.innerTrack()->hitPattern().numberOfValidPixelHits() : 0);
      mu_numberOfValidTrackerHits.push_back((!mu.innerTrack().isNull()) ? mu.innerTrack()->hitPattern().numberOfValidTrackerHits() : 0);
      mu_numberInnerHitsMissing.push_back(mu.innerTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      mu_trackerLayersWithMeasurement.push_back((!mu.innerTrack().isNull()) ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0);
      mu_numberInnerHits.push_back((!mu.globalTrack().isNull()) ? mu.globalTrack()->hitPattern().numberOfValidMuonHits() : (!mu.outerTrack().isNull() ? mu.outerTrack()->hitPattern().numberOfValidMuonHits() : 0));

      mu_relIso0p4.push_back(getPFIso(mu));
      mu_trackPt.push_back(mu.innerTrack()->pt());
      mu_trackPtErr.push_back(mu.innerTrack()->ptError());
      
      mu_dxy.push_back(mu.dB(pat::Muon::PV2D));
      mu_dz.push_back(mu.dB(pat::Muon::PVDZ));
      mu_3dIP.push_back(mu.dB(pat::Muon::PV3D));
      mu_3dIPSig.push_back(fabs(mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D)));
      
      fillLeptonJetVariables(NULL, &mu, jetHandle, pv, *rhoHandle, false);
      
      bool new_vtx = false;
      double ptdiff, normchi2;
      double minptdiff = 10, minnormchi2 = 10000;
      int nVtx = 0;
      for(const reco::Vertex& vtx : secVertices) {
	 for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++) {
	    reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
	    for( size_t iCand=0;iCand<mu.numberOfSourceCandidatePtrs();++iCand ) {
	       if( !(mu.sourceCandidatePtr(iCand).isNonnull() and mu.sourceCandidatePtr(iCand).isAvailable()) ) continue;
	       ptdiff   = fabs(mu.sourceCandidatePtr(iCand)->pt() - vtxTrack->pt());
	       normchi2 = fabs(vtx.chi2()/vtx.ndof());
	       
	       if(ptdiff < 0.001 and (ptdiff < minptdiff or (ptdiff == minptdiff and normchi2 < minnormchi2))) {
		  new_vtx = true;
		  mu_lIVF_match[nMuonsSel] = true;
		  minptdiff   = ptdiff;
		  minnormchi2 = normchi2;
	       }
	    }
	 }       
	 if(new_vtx) {
	    mu_IVF_x.push_back(vtx.x());
	    mu_IVF_y.push_back(vtx.y());
	    mu_IVF_z.push_back(vtx.z());
	    mu_IVF_cx.push_back(vtx.xError());
	    mu_IVF_cy.push_back(vtx.yError());
	    mu_IVF_cz.push_back(vtx.zError());
	    mu_IVF_df.push_back(vtx.ndof());
	    mu_IVF_chi2.push_back(vtx.chi2());
	    mu_IVF_pt.push_back(vtx.p4().pt());
	    mu_IVF_eta.push_back(vtx.p4().eta());
	    mu_IVF_phi.push_back(vtx.p4().phi());
	    mu_IVF_E.push_back(vtx.p4().energy());
	    mu_IVF_mass.push_back(vtx.p4().mass());
	    mu_IVF_muid.push_back(nMuonsSel);
	    
	    mu_IVF_ntracks.push_back(0);
	    for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++) {
	       if(mu_IVF_ntracks.back() == ntrack_max) break;
	       reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
	       mu_IVF_trackpt.push_back(vtxTrack->pt());
	       mu_IVF_tracketa.push_back(vtxTrack->eta());
	       mu_IVF_trackphi.push_back(vtxTrack->phi());
	       mu_IVF_trackE.push_back(vtxTrack->p());
	       mu_IVF_trackcharge.push_back(vtxTrack->charge());
	       mu_IVF_trackdxy.push_back(std::abs(vtxTrack->dxy(pv.position())));
	       mu_IVF_trackdz.push_back(std::abs(vtxTrack->dz(pv.position())));
	       mu_IVF_trackmuid.push_back(nMuonsSel);
	       mu_IVF_trackvtxid.push_back(nVtx);
	       mu_IVF_ntracks.back()++;
	    }
	    nVtx++;
	    new_vtx = false;
	 }	 
      }      
      nMuonsSel += 1;
   }
   
   auto dispJetElectronTab = std::make_unique<nanoaod::FlatTable>(nElectronsSel, "DispJetElectron", false, false);
   auto dispJetMuonTab = std::make_unique<nanoaod::FlatTable>(nMuonsSel, "DispJetMuon", false, false);

   dispJetElectronTab->addColumn<int>("idx", el_idx, "");
   dispJetElectronTab->addColumn<bool>("lIVF_match", el_lIVF_match, "");
   dispJetElectronTab->addColumn<bool>("isEB", el_isEB, "");
   dispJetElectronTab->addColumn<bool>("isEE", el_isEE, "");
   dispJetElectronTab->addColumn<float>("superClusterOverP", el_superClusterOverP, "");
   dispJetElectronTab->addColumn<float>("ecalEnergy", el_ecalEnergy, "");
   dispJetElectronTab->addColumn<float>("dEtaInSeed", el_dEtaInSeed, "");
   dispJetElectronTab->addColumn<int>("numberInnerHitsMissing", el_numberInnerHitsMissing, "");
   dispJetElectronTab->addColumn<int>("numberOfValidPixelHits", el_numberOfValidPixelHits, "");
   dispJetElectronTab->addColumn<int>("numberOfValidTrackerHits", el_numberOfValidTrackerHits, "");
   dispJetElectronTab->addColumn<float>("relIso0p4", el_relIso0p4, "");
   dispJetElectronTab->addColumn<float>("ptRatio", el_ptRatio, "");
   dispJetElectronTab->addColumn<float>("ptRel", el_ptRel, "");
   dispJetElectronTab->addColumn<float>("closestJetTag_b", el_closestJetTag_b, "");
   dispJetElectronTab->addColumn<float>("closestJetTag_bb", el_closestJetTag_bb, "");
   dispJetElectronTab->addColumn<float>("closestJetTag_lepb", el_closestJetTag_lepb, "");
   dispJetElectronTab->addColumn<float>("closestJetTag", el_closestJetTag, "");
   dispJetElectronTab->addColumn<int>("selectedTrackMult", el_selectedTrackMult, "");
   dispJetElectronTab->addColumn<float>("sigmaIetaIeta", el_sigmaIetaIeta, "");
   dispJetElectronTab->addColumn<float>("deltaPhiSuperClusterTrack", el_deltaPhiSuperClusterTrack, "");
   dispJetElectronTab->addColumn<float>("deltaEtaSuperClusterTrack", el_deltaEtaSuperClusterTrack, "");
   dispJetElectronTab->addColumn<float>("eInvMinusPInv", el_eInvMinusPInv, "");
   dispJetElectronTab->addColumn<float>("hOverE", el_hOverE, "");
   dispJetElectronTab->addColumn<float>("dxy", el_dxy, "");
   dispJetElectronTab->addColumn<float>("dz", el_dz, "");
   dispJetElectronTab->addColumn<float>("3dIP", el_3dIP, "");
   dispJetElectronTab->addColumn<float>("3dIPSig", el_3dIPSig, "");
   
   auto dispJetElectronVtxTab = std::make_unique<nanoaod::FlatTable>(el_IVF_x.size(), "DispJetElectronVtx", false, false);
   dispJetElectronVtxTab->addColumn<int>("IVF_df", el_IVF_df, "");
   dispJetElectronVtxTab->addColumn<int>("IVF_ntracks", el_IVF_ntracks, "");
   dispJetElectronVtxTab->addColumn<int>("IVF_elid", el_IVF_elid, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_x", el_IVF_x, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_y", el_IVF_y, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_z", el_IVF_z, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_cx", el_IVF_cx, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_cy", el_IVF_cy, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_cz", el_IVF_cz, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_chi2", el_IVF_chi2, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_pt", el_IVF_pt, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_eta", el_IVF_eta, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_phi", el_IVF_phi, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_E", el_IVF_E, "");
   dispJetElectronVtxTab->addColumn<float>("IVF_mass", el_IVF_mass, "");   
   
   int nTracksElectron = 0;
   for( unsigned int iv=0;iv<el_IVF_ntracks.size();iv++ ) {
      nTracksElectron += std::min(el_IVF_ntracks[iv], ntrack_max);
   }
   auto dispJetElectronTrkTab = std::make_unique<nanoaod::FlatTable>(nTracksElectron, "DispJetElectronTrk", false, false);
   dispJetElectronTrkTab->addColumn<int>("IVF_trackcharge", el_IVF_trackcharge, "");
   dispJetElectronTrkTab->addColumn<float>("IVF_trackpt", el_IVF_trackpt, "");
   dispJetElectronTrkTab->addColumn<float>("IVF_tracketa", el_IVF_tracketa, "");
   dispJetElectronTrkTab->addColumn<float>("IVF_trackphi", el_IVF_trackphi, "");
   dispJetElectronTrkTab->addColumn<float>("IVF_trackE", el_IVF_trackE, "");
   dispJetElectronTrkTab->addColumn<float>("IVF_trackdxy", el_IVF_trackdxy, "");
   dispJetElectronTrkTab->addColumn<float>("IVF_trackdz", el_IVF_trackdz, "");
   dispJetElectronTrkTab->addColumn<int>("IVF_trackelid", el_IVF_trackelid, "");
   dispJetElectronTrkTab->addColumn<int>("IVF_trackvtxid", el_IVF_trackvtxid, "");

   dispJetMuonTab->addColumn<int>("idx", mu_idx, "");
   dispJetMuonTab->addColumn<bool>("lIVF_match", mu_lIVF_match, "");
   dispJetMuonTab->addColumn<float>("innerTrackValidFraction", mu_innerTrackValidFraction, "");
   dispJetMuonTab->addColumn<float>("globalTrackNormalizedChi2", mu_globalTrackNormalizedChi2, "");
   dispJetMuonTab->addColumn<float>("CQChi2Position", mu_CQChi2Position, "");
   dispJetMuonTab->addColumn<float>("CQTrackKink", mu_CQTrackKink, "");
   dispJetMuonTab->addColumn<int>("numberOfMatchedStation", mu_numberOfMatchedStation, "");
   dispJetMuonTab->addColumn<int>("numberOfValidPixelHits", mu_numberOfValidPixelHits, "");
   dispJetMuonTab->addColumn<int>("numberOfValidTrackerHits", mu_numberOfValidTrackerHits, "");
   dispJetMuonTab->addColumn<int>("numberInnerHitsMissing", mu_numberInnerHitsMissing, "");
   dispJetMuonTab->addColumn<int>("trackerLayersWithMeasurement", mu_trackerLayersWithMeasurement, "");
   dispJetMuonTab->addColumn<int>("numberInnerHits", mu_numberInnerHits, "");
   dispJetMuonTab->addColumn<float>("relIso0p4", mu_relIso0p4, "");
   dispJetMuonTab->addColumn<float>("ptRatio", mu_ptRatio, "");
   dispJetMuonTab->addColumn<float>("ptRel", mu_ptRel, "");
   dispJetMuonTab->addColumn<float>("closestJetTag_b", mu_closestJetTag_b, "");
   dispJetMuonTab->addColumn<float>("closestJetTag_bb", mu_closestJetTag_bb, "");
   dispJetMuonTab->addColumn<float>("closestJetTag_lepb", mu_closestJetTag_lepb, "");
   dispJetMuonTab->addColumn<float>("closestJetTag", mu_closestJetTag, "");
   dispJetMuonTab->addColumn<int>("selectedTrackMult", mu_selectedTrackMult, "");
   dispJetMuonTab->addColumn<float>("trackPt", mu_trackPt, "");
   dispJetMuonTab->addColumn<float>("trackPtErr", mu_trackPtErr, "");
   dispJetMuonTab->addColumn<float>("dxy", mu_dxy, "");
   dispJetMuonTab->addColumn<float>("dz", mu_dz, "");
   dispJetMuonTab->addColumn<float>("3dIP", mu_3dIP, "");
   dispJetMuonTab->addColumn<float>("3dIPSig", mu_3dIPSig, "");

   auto dispJetMuonVtxTab = std::make_unique<nanoaod::FlatTable>(mu_IVF_x.size(), "DispJetMuonVtx", false, false);
   dispJetMuonVtxTab->addColumn<int>("IVF_df", mu_IVF_df, "");
   dispJetMuonVtxTab->addColumn<int>("IVF_ntracks", mu_IVF_ntracks, "");
   dispJetMuonVtxTab->addColumn<int>("IVF_muid", mu_IVF_muid, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_x", mu_IVF_x, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_y", mu_IVF_y, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_z", mu_IVF_z, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_cx", mu_IVF_cx, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_cy", mu_IVF_cy, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_cz", mu_IVF_cz, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_chi2", mu_IVF_chi2, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_pt", mu_IVF_pt, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_eta", mu_IVF_eta, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_phi", mu_IVF_phi, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_E", mu_IVF_E, "");
   dispJetMuonVtxTab->addColumn<float>("IVF_mass", mu_IVF_mass, "");
   
   int nTracksMuon = 0;
   for( unsigned int iv=0;iv<mu_IVF_ntracks.size();iv++ ) {
      nTracksMuon += std::min(mu_IVF_ntracks[iv], ntrack_max);
   }
   auto dispJetMuonTrkTab = std::make_unique<nanoaod::FlatTable>(nTracksMuon, "DispJetMuonTrk", false, false);
   dispJetMuonTrkTab->addColumn<int>("IVF_trackcharge", mu_IVF_trackcharge, "");
   dispJetMuonTrkTab->addColumn<float>("IVF_trackpt", mu_IVF_trackpt, "");
   dispJetMuonTrkTab->addColumn<float>("IVF_tracketa", mu_IVF_tracketa, "");
   dispJetMuonTrkTab->addColumn<float>("IVF_trackphi", mu_IVF_trackphi, "");
   dispJetMuonTrkTab->addColumn<float>("IVF_trackE", mu_IVF_trackE, "");
   dispJetMuonTrkTab->addColumn<float>("IVF_trackdxy", mu_IVF_trackdxy, "");
   dispJetMuonTrkTab->addColumn<float>("IVF_trackdz", mu_IVF_trackdz, "");
   dispJetMuonTrkTab->addColumn<int>("IVF_trackmuid", mu_IVF_trackmuid, "");
   dispJetMuonTrkTab->addColumn<int>("IVF_trackvtxid", mu_IVF_trackvtxid, "");

   iEvent.put(std::move(dispJetElectronTab), "DispJetElectron");
   iEvent.put(std::move(dispJetElectronVtxTab), "DispJetElectronVtx");
   iEvent.put(std::move(dispJetElectronTrkTab), "DispJetElectronTrk");
   iEvent.put(std::move(dispJetMuonTab), "DispJetMuon");
   iEvent.put(std::move(dispJetMuonVtxTab), "DispJetMuonVtx");
   iEvent.put(std::move(dispJetMuonTrkTab), "DispJetMuonTrk");
}

float DispJetTableProducer::dEtaInSeed(const pat::Electron* el) const
{   
   if( el->superCluster().isNonnull() and el->superCluster()->seed().isNonnull()) return el->deltaEtaSuperClusterTrackAtVtx() - el->superCluster()->eta() + el->superCluster()->seed()->eta();
   else return std::numeric_limits<float>::max();
}

template< typename T1, typename T2 > bool isSourceCandidatePtrMatch( const T1& lhs, const T2& rhs ) {
   
   for( size_t lhsIndex = 0; lhsIndex < lhs.numberOfSourceCandidatePtrs(); ++lhsIndex ) {
      auto lhsSourcePtr = lhs.sourceCandidatePtr( lhsIndex );
      for( size_t rhsIndex = 0; rhsIndex < rhs.numberOfSourceCandidatePtrs(); ++rhsIndex ) {
	 auto rhsSourcePtr = rhs.sourceCandidatePtr( rhsIndex );
	 if( lhsSourcePtr == rhsSourcePtr ) {
	    return true;
	 }
      }
   }
   
   return false;
}

const pat::Jet* findMatchedJet(const reco::Candidate& lepton, const edm::Handle< std::vector< pat::Jet > >& jets, const bool oldMatching) {

   const pat::Jet* matchedJetPtr = nullptr;
   
   if( oldMatching ) {
      for( auto& jet : *jets ) {
	 if( jet.pt() <= 5 || fabs( jet.eta() ) >= 3 ) continue;
	 if( ( matchedJetPtr == nullptr) || reco::deltaR( jet, lepton ) < reco::deltaR( *matchedJetPtr, lepton ) ) {
	    matchedJetPtr = &jet;
	 }	 
      }
      if( matchedJetPtr != nullptr && reco::deltaR( lepton, *matchedJetPtr ) > 0.4 ) {
	 matchedJetPtr = nullptr;
      }
   } else {
      for( auto& jet : *jets ) {
	 if( isSourceCandidatePtrMatch( lepton, jet ) ) {
	    return &jet;
	 }
      }
   }
   
   return matchedJetPtr;
}

void DispJetTableProducer::fillLeptonJetVariables(const reco::GsfElectron *el, const reco::Muon *mu, edm::Handle< std::vector< pat::Jet > >& jets, const reco::Vertex& vertex, const double rho, const bool oldMatching) {
   
   const reco::Candidate *cand = (el) ? dynamic_cast<const reco::Candidate*>(el) : dynamic_cast<const reco::Candidate*>(mu);
   const pat::Jet* matchedJetPtr = findMatchedJet( *cand, jets, oldMatching );
   
   bool isElectron = (el);
   
   if( matchedJetPtr == nullptr ) {
      if( isElectron ) {
	 float ptRatio = ( oldMatching ? 1. : 1. / ( 1. + el_relIso0p4.back() ) );
	 el_ptRatio.push_back(ptRatio);
	 el_ptRel.push_back(0);
	 el_selectedTrackMult.push_back(0);
	 el_closestJetTag_b.push_back(0);
	 el_closestJetTag_bb.push_back(0);
	 el_closestJetTag_lepb.push_back(0);
	 el_closestJetTag.push_back(0);
      } else {	   
	 float ptRatio = ( oldMatching ? 1. : 1. / ( 1. + mu_relIso0p4.back() ) );
	 mu_ptRatio.push_back(ptRatio);	 
	 mu_ptRel.push_back(0);
	 mu_selectedTrackMult.push_back(0);
	 mu_closestJetTag_b.push_back(0);
	 mu_closestJetTag_bb.push_back(0);
	 mu_closestJetTag_lepb.push_back(0);
	 mu_closestJetTag.push_back(0);
      }
   } else {
      const pat::Jet& jet = *matchedJetPtr;
      
      auto rawJetP4 = jet.correctedP4("Uncorrected");
      auto leptonP4 = cand->p4();

      bool leptonEqualsJet = ( ( rawJetP4 - leptonP4 ).P() < 1e-4 );
      
      if( leptonEqualsJet && !oldMatching ) {
	 if( isElectron ) {
	    el_ptRatio.push_back(1);
	    el_ptRel.push_back(0);
	    el_selectedTrackMult.push_back(0);
	 } else {
	    mu_ptRatio.push_back(1);
	    mu_ptRel.push_back(0);
	    mu_selectedTrackMult.push_back(0);	    
	 }	 
      } else {
	 auto L1JetP4 = jet.correctedP4("L1FastJet");
	 double L2L3JEC = jet.pt()/L1JetP4.pt();
	 auto lepAwareJetP4 = ( L1JetP4 - leptonP4 )*L2L3JEC + leptonP4;

	 float ptRatio = cand->pt() / lepAwareJetP4.pt();
	 float ptRel = leptonP4.Vect().Cross( (lepAwareJetP4 - leptonP4 ).Vect().Unit() ).R();
	 if( isElectron ) {
	    el_ptRatio.push_back(ptRatio);
	    el_ptRel.push_back(ptRel);
	    el_selectedTrackMult.push_back(0);
	 } else {
	    mu_ptRatio.push_back(ptRatio);
	    mu_ptRel.push_back(ptRel);
	    mu_selectedTrackMult.push_back(0);
	 }
	 	 
	 for( const auto daughterPtr : jet.daughterPtrVector() ) {
	    const pat::PackedCandidate& daughter = *( (const pat::PackedCandidate*) daughterPtr.get() );
	    
	    if( daughter.charge() == 0 ) continue;
	    if( daughter.fromPV() < 2 ) continue;
	    if( reco::deltaR( daughter, *cand ) > 0.4 ) continue;
	    if( !daughter.hasTrackDetails() ) continue;
	    
	    auto daughterTrack = daughter.pseudoTrack();
	    
	    if( daughterTrack.pt() <= 1 ) continue;
	    if( daughterTrack.hitPattern().numberOfValidHits() < 8 ) continue;
	    if( daughterTrack.hitPattern().numberOfValidPixelHits() < 2 ) continue;
	    if( daughterTrack.normalizedChi2() >= 5 ) continue;
	    if( std::abs( daughterTrack.dz( vertex.position() ) ) >= 17 ) continue;
	    if( std::abs( daughterTrack.dxy( vertex.position() ) ) >= 0.2 ) continue;
	    if( isElectron ) ++el_selectedTrackMult.back();
	    else ++mu_selectedTrackMult.back();
	 }
      }

      float probb = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
      float probbb = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
      float problepb = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
      if( isElectron ) {
	 el_closestJetTag_b.push_back(probb);
	 el_closestJetTag_bb.push_back(probbb);
	 el_closestJetTag_lepb.push_back(problepb);
	 el_closestJetTag.push_back(probb + probbb + problepb);
	 if( std::isnan( el_closestJetTag.back() ) ) el_closestJetTag.back() = 0.;
      }
      else {
	 mu_closestJetTag_b.push_back(probb);
	 mu_closestJetTag_bb.push_back(probbb);
	 mu_closestJetTag_lepb.push_back(problepb);
	 mu_closestJetTag.push_back(probb + probbb + problepb);
	 if( std::isnan( mu_closestJetTag.back() ) ) mu_closestJetTag.back() = 0.;	 
      }
   }
}

float DispJetTableProducer::getPFIso(const pat::Muon& muon) const {
   return (muon.pfIsolationR04().sumChargedHadronPt +
	   std::max(0.,
		    muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt -
		    0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();
}

float DispJetTableProducer::getPFIso(const pat::Electron& electron) const {
   return electron.userFloat("PFIsoAll04") / electron.pt();
}
   
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DispJetTableProducer);
