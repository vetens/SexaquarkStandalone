#include "../interface/FlatTreeProducerTracking.h"
#include <typeinfo>

FlatTreeProducerTracking::FlatTreeProducerTracking(edm::ParameterSet const& pset):
  m_bsTag(pset.getParameter<edm::InputTag>("beamspot")),
  m_offlinePVTag(pset.getParameter<edm::InputTag>("offlinePV")),
  m_genParticlesTag_GEN(pset.getParameter<edm::InputTag>("genCollection_GEN")),
  m_genParticlesTag_SIM_GEANT(pset.getParameter<edm::InputTag>("genCollection_SIM_GEANT")),
  m_generalTracksTag(pset.getParameter<edm::InputTag>("generalTracksCollection")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_V0KsTag(pset.getParameter<edm::InputTag>("V0KsCollection")),
  m_V0LTag(pset.getParameter<edm::InputTag>("V0LCollection")),
  m_trackAssociatorTag(pset.getParameter<edm::InputTag>("trackAssociators")),
  m_TPTag(pset.getParameter<edm::InputTag>("TrackingParticles")),
//  m_PileupInfoTag(pset.getParameter<edm::InputTag>("PileupInfo")),

  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_offlinePVToken    (consumes<vector<reco::Vertex>>(m_offlinePVTag)),
  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM_GEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM_GEANT)),
  //m_generalTracksToken(consumes<vector<reco::Track> >(m_generalTracksTag)),
  m_generalTracksToken(consumes<View<reco::Track> >(m_generalTracksTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_V0KsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0KsTag)),
  m_V0LToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0LTag)),
  m_trackAssociatorToken(consumes<reco::TrackToTrackingParticleAssociator> (m_trackAssociatorTag)),
  m_TPToken(consumes<vector<TrackingParticle> >(m_TPTag))
//  m_PileupInfoToken(consumes<vector<PileupSummaryInfo> >(m_PileupInfoTag))
  


{

}


void FlatTreeProducerTracking::beginJob() {

	// Initialize when class is created                                                              |  m_TPTag(pset.getParameter<edm::InputTag>("TrackingParticles")),
        edm::Service<TFileService> fs ; 
	
	_tree_tracks = fs->make <TTree>("FlatTreeTracks","treeTracks");
	
	_tree_tracks->Branch("_tp_pt",&_tp_pt);
	_tree_tracks->Branch("_tp_eta",&_tp_eta);
	_tree_tracks->Branch("_tp_phi",&_tp_phi);
	_tree_tracks->Branch("_tp_pz",&_tp_pz);

	_tree_tracks->Branch("_tp_Lxy_beamspot",&_tp_Lxy_beamspot);
	_tree_tracks->Branch("_tp_vz_beamspot",&_tp_vz_beamspot);
	_tree_tracks->Branch("_tp_dxy_beamspot",&_tp_dxy_beamspot);
	_tree_tracks->Branch("_tp_dz_beamspot",&_tp_dz_beamspot);

	_tree_tracks->Branch("_tp_numberOfTrackerLayers",&_tp_numberOfTrackerLayers);
	_tree_tracks->Branch("_tp_charge",&_tp_charge);

	_tree_tracks->Branch("_tp_reconstructed",&_tp_reconstructed);
	_tree_tracks->Branch("_tp_isAntiSTrack",&_tp_isAntiSTrack);

	_tree_tracks->Branch("_tp_etaOfGrandMotherAntiS",&_tp_etaOfGrandMotherAntiS);

	_tree_tracks->Branch("_matchedTrack_pt",&_matchedTrack_pt);
	_tree_tracks->Branch("_matchedTrack_eta",&_matchedTrack_eta);
	_tree_tracks->Branch("_matchedTrack_phi",&_matchedTrack_phi);
	_tree_tracks->Branch("_matchedTrack_pz",&_matchedTrack_pz);

	_tree_tracks->Branch("_matchedTrack_chi2",&_matchedTrack_chi2);
	_tree_tracks->Branch("_matchedTrack_ndof",&_matchedTrack_ndof);
	_tree_tracks->Branch("_matchedTrack_charge",&_matchedTrack_charge);
	_tree_tracks->Branch("_matchedTrack_dxy_beamspot",&_matchedTrack_dxy_beamspot);
	_tree_tracks->Branch("_matchedTrack_dz_beamspot",&_matchedTrack_dz_beamspot);

	_tree_tracks->Branch("_matchedTrack_trackQuality",&_matchedTrack_trackQuality);
	_tree_tracks->Branch("_matchedTrack_isLooper",&_matchedTrack_isLooper);



}

void FlatTreeProducerTracking::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

 
  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

  //primary vertex
  edm::Handle<vector<reco::Vertex>> h_offlinePV;
  iEvent.getByToken(m_offlinePVToken, h_offlinePV);
  int nPVs = h_offlinePV->size(); 

  //SIM particles: normal Gen particles or PlusGEANT
  edm::Handle<vector<reco::GenParticle>> h_genParticles;
  //iEvent.getByToken(m_genParticlesToken_GEN, h_genParticles);
  iEvent.getByToken(m_genParticlesToken_SIM_GEANT, h_genParticles);

  //General tracks particles
  //edm::Handle<vector<reco::Track>> h_generalTracks;
  edm::Handle<View<reco::Track>> h_generalTracks;
  iEvent.getByToken(m_generalTracksToken, h_generalTracks);

  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);

  //V0 Kshorts
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks;
  iEvent.getByToken(m_V0KsToken, h_V0Ks);

  //V0 Lambdas
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L;
  iEvent.getByToken(m_V0LToken, h_V0L);

  //trackingparticle collection
  edm::Handle<TrackingParticleCollection>  h_TP ;
  iEvent.getByToken(m_TPToken,h_TP);
  TrackingParticleCollection const & TPColl = *(h_TP.product());
  //track associator on hits
  edm::Handle< reco::TrackToTrackingParticleAssociator>  h_trackAssociator;
  iEvent.getByToken(m_trackAssociatorToken, h_trackAssociator);

  //Pileup information from GEN level, so this is the true pileup
//  edm::Handle<vector<PileupSummaryInfo> > h_PileupInfo;
//  iEvent.getByToken(m_PileupInfoToken, h_PileupInfo);
//  int nPU = h_PileupInfo->size(); 



  //beamspot
  TVector3 beamspot(0,0,0);
  TVector3 beamspotVariance(0,0,0);
  if(h_bs.isValid()){  
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));			
  }
  else{
	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!beamspot collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  }

  if(!h_generalTracks.isValid()){ std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!generalTracks collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;}
  if(!h_TP.isValid()){ std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!trackingParticles collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;}
  if(!h_trackAssociator.isValid()){ std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!trackAssociator collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;}

//evaluate tracking performance
  if(h_generalTracks.isValid() && h_TP.isValid() && h_trackAssociator.isValid()){
	for(size_t i=0; i<TPColl.size(); ++i) {
	  //first investigate whether this tp is a dauhgter of an AntiS. You need this because these granddaughter tracks you want to investigate separately from the rest of the tracks. 
	  const TrackingParticle& tp = TPColl[i];
	  std::vector<double> tpIsGrandDaughterAntiS = AnalyzerAllSteps::isTpGrandDaughterAntiS(TPColl, tp);
	  //std::cout << "---------------------------" << std::endl;
	  //std::cout << "tp pdgId: " << tp.pdgId() << std::endl;
	  //std::cout << "tp is a granddaughter of an antiS " << tpIsGrandDaughterAntiS << std::endl;
	  bool matchingTrackFound = false;
	  const reco::Track *matchedTrackPointer = nullptr;
	  TrackingParticleRef tpr(h_TP,i);
	  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
	  reco::SimToRecoCollection simRecCollL;
	  reco::SimToRecoCollection const * simRecCollP=nullptr;
	  simRecCollL = std::move(h_trackAssociator->associateSimToReco(h_generalTracks,h_TP));
	  simRecCollP = &simRecCollL;
          reco::SimToRecoCollection const & simRecColl = *simRecCollP;

	  if(simRecColl.find(tpr) != simRecColl.end()){
		  auto const & rt = simRecColl[tpr];
		  if (rt.size()!=0) {
		    // isRecoMatched = true; // UNUSED
		    matchedTrackPointer = rt.begin()->first.get();
		    matchingTrackFound  = true;
		    //std::cout << "found a matching track with quality: " << matchedTrackQuality << std::endl;
	  	    //std::cout << "Track matching result: " << matchingTrackFound << " ,for a GEN track with a pt of " << tpr->pt() << " and a RECO matched track with a pt of " << matchedTrackPointer->pt() << std::endl;	
		  }
	  }else{
		    //std::cout << "no matching track found" << std::endl;
		    matchingTrackFound = false;
	  }
        
	
	if(tp.charge() != 0 )FillTreesTracks(tp, beamspot, nPVs, matchedTrackPointer, matchingTrackFound,tpIsGrandDaughterAntiS);
       	  
	  //now in the above two catagories you are missing all the tracks from the antiS granddaughters, so when you encounter an antiS you should still have to fill four it's 4 potential granddaughters the histograms and also one where you count how many antiS have four granddaughters reconstructed
	  //if(tp.pdgId() == AnalyzerAllSteps::pdgIdAntiS)FillHistosAntiSTracks(tp, beamspot, TPColl,  h_TP, h_trackAssociator, h_generalTracks, h_V0Ks, h_V0L);
	}
  }
  
/*
  int nAntiSThisEvent = 0;
  int nAntiSInteractThisEvent = 0;
  int nKsThisEvent = 0;
  int nAntiLambdaThisEvent = 0;
  if(h_genParticles.isValid()){
      for(unsigned int i = 0; i < h_genParticles->size(); ++i){//loop all genparticlesPlusGEANT
		const reco::Candidate * genParticle = &h_genParticles->at(i);
		//all antiS, so you will also be counting the loopers here, so there are duplicates, but I have the flattrees where I make these histos so it's ok
		bool genParticleIsAntiS = false;
		if(genParticle->pdgId() == AnalyzerAllSteps::pdgIdAntiS) genParticleIsAntiS = true;
		if(genParticleIsAntiS){histos_th1f["h_GEN_nAntiSTotal"]->Fill(0); nAntiSThisEvent++; FillHistosGENAntiS(genParticle, beamspot);}

		//all Ks
		bool genParticleIsKs = false;	
		bool genParticleMotherIsAntiS = false;//Ks as daughter of an antiS	
		if(abs(genParticle->pdgId()) == AnalyzerAllSteps::pdgIdKs)genParticleIsKs =true;
		if(genParticle->numberOfMothers()==1){if(genParticle->mother(0)->pdgId()==AnalyzerAllSteps::pdgIdAntiS){genParticleMotherIsAntiS=true;}};
		if(genParticleIsKs && !genParticleMotherIsAntiS){nKsThisEvent++; FillHistosGENKsNonAntiS(genParticle,beamspot);};
		if(genParticleIsKs && genParticleMotherIsAntiS)FillHistosGENKsAntiS(genParticle,beamspot);
		//all antiLambda
		bool genParticleIsAntiLambda = false;
		if(genParticle->pdgId() == AnalyzerAllSteps::pdgIdAntiLambda)genParticleIsAntiLambda=true;
		if(genParticleIsAntiLambda && !genParticleMotherIsAntiS){nAntiLambdaThisEvent++;FillHistosGENAntiLambdaNonAntiS(genParticle,beamspot);};
		if(genParticleIsAntiLambda && genParticleMotherIsAntiS)FillHistosGENAntiLambdaAntiS(genParticle,beamspot);


		//now start looking at the RECO efficiencies of daughters and antiS, you can do this based on deltaR, because these particles are not charged. For all of these you want however that at least they have the correct daughters and for the antiS the correct granddaughters:
		if(genParticle->numberOfDaughters()==2){
			if(genParticleIsAntiS){histos_th1f["h_GEN_nAntiSTotal"]->Fill(1);nAntiSInteractThisEvent++;}			
			int daughterParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle);
			//for Ks
			if(genParticleIsKs && !genParticleMotherIsAntiS && daughterParticlesTypes == 1)RecoEvaluationKsNonAntiS(genParticle,h_V0Ks,beamspot);
			if(genParticleIsKs && genParticleMotherIsAntiS && daughterParticlesTypes == 1)RecoEvaluationKsAntiS(genParticle,h_V0Ks,beamspot);
			//for AntiL
			if(genParticleIsAntiLambda && !genParticleMotherIsAntiS && daughterParticlesTypes == 2)RecoEvaluationAntiLambdaNonAntiS(genParticle,h_V0L,beamspot);
			if(genParticleIsAntiLambda && genParticleMotherIsAntiS && daughterParticlesTypes == 2)RecoEvaluationAntiLambdaAntiS(genParticle,h_V0L,beamspot);
			//for antiS
			if(genParticle->daughter(0)->numberOfDaughters()==2 && genParticle->daughter(1)->numberOfDaughters()==2 && daughterParticlesTypes == 3){
				FillHistosGENInteractingAntiS(genParticle, beamspot);
				int graddaughters0ParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle->daughter(0));
				int graddaughters1ParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle->daughter(1));
				if(genParticleIsAntiS && ((graddaughters0ParticlesTypes == 1 && graddaughters1ParticlesTypes == 2) || (graddaughters1ParticlesTypes == 1 && graddaughters0ParticlesTypes == 2))){histos_th1f["h_GEN_nAntiSTotal"]->Fill(2); RecoEvaluationAntiS(genParticle,h_sCands,beamspot, beamspotVariance, h_offlinePV);}
			}
		}
      }//for(unsigned int i = 0; i < h_genParticles->size(); ++i)
  }//if(h_genParticles.isValid())
  histos_th1f["h_GEN_nAntiS"]->Fill(nAntiSThisEvent);
  histos_th1f["h_GEN_nAntiSInteract"]->Fill(nAntiSInteractThisEvent);
  histos_th1f["h_GEN_nKs_NonAntiS"]->Fill(nKsThisEvent);
  histos_th1f["h_GEN_nAntiLambda_NonAntiS"]->Fill(nAntiLambdaThisEvent); 
*/

} //end of analyzer


void FlatTreeProducerTracking::FillTreesTracks(const TrackingParticle& tp, TVector3 beamspot, int nPVs, const reco::Track *matchedTrackPointer, bool matchingTrackFound, std::vector<double> tpIsGrandDaughterAntiS){

	math::XYZPoint beamspotPoint(beamspot.X(),beamspot.Y(),beamspot.Z());

	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	double dz = AnalyzerAllSteps::dz_line_point(tpCreationVertex,tpMomentum,beamspot);

	InitTracking();	

	_tp_pt.push_back(tp.pt());
	_tp_eta.push_back(tp.eta());
	_tp_phi.push_back(tp.phi());
	_tp_pz.push_back(tp.pz());

	_tp_Lxy_beamspot.push_back(Lxy);
	_tp_vz_beamspot.push_back(tp.vz());
	_tp_dxy_beamspot.push_back(dxy);
	_tp_dz_beamspot.push_back(dz);

	_tp_numberOfTrackerLayers.push_back(tp.numberOfTrackerLayers());
	_tp_charge.push_back(tp.charge());

	_tp_reconstructed.push_back(matchingTrackFound);
	_tp_isAntiSTrack.push_back(tpIsGrandDaughterAntiS[0]);

	_tp_etaOfGrandMotherAntiS.push_back(tpIsGrandDaughterAntiS[1]);

	//save also the information of the actual reconstructed track if the tp was matched to a track
	if(matchingTrackFound){
		_matchedTrack_pt.push_back(matchedTrackPointer->pt());
		_matchedTrack_eta.push_back(matchedTrackPointer->eta());
		_matchedTrack_phi.push_back(matchedTrackPointer->phi());
		_matchedTrack_pz.push_back(matchedTrackPointer->pz());

		_matchedTrack_chi2.push_back(matchedTrackPointer->chi2());
		_matchedTrack_ndof.push_back(matchedTrackPointer->ndof());
		_matchedTrack_charge.push_back(matchedTrackPointer->charge());
		_matchedTrack_dxy_beamspot.push_back(matchedTrackPointer->dxy(beamspotPoint));
		_matchedTrack_dz_beamspot.push_back(matchedTrackPointer->dsz(beamspotPoint));
		
		int matchedTrackQuality = AnalyzerAllSteps::trackQualityAsInt(matchedTrackPointer);
		_matchedTrack_trackQuality.push_back(matchedTrackQuality);

		_matchedTrack_isLooper.push_back(matchedTrackPointer->isLooper());
	}
	else{
		_matchedTrack_pt.push_back(-999);
		_matchedTrack_eta.push_back(-999);
		_matchedTrack_phi.push_back(-999);
		_matchedTrack_pz.push_back(-999);

		_matchedTrack_chi2.push_back(-999);
		_matchedTrack_ndof.push_back(-999);
		_matchedTrack_charge.push_back(-999);
		_matchedTrack_dxy_beamspot.push_back(-999);
		_matchedTrack_dz_beamspot.push_back(-999);
		
		int matchedTrackQuality = -999;
		_matchedTrack_trackQuality.push_back(matchedTrackQuality);

		_matchedTrack_isLooper.push_back(-999);

	}

	_tree_tracks->Fill();
	


}

/*
void FlatTreeProducerTracking::FillHistosNonAntiSTracksRECO(const TrackingParticle& tp, TVector3 beamspot, int nPVs, int matchedTrackQuality){
	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	double dz = AnalyzerAllSteps::dz_line_point(tpCreationVertex,tpMomentum,beamspot);
	histos_th1f["h_NonAntiSTrack_RECO_pt"]->Fill(tp.pt());	
	histos_th1f["h_NonAntiSTrack_RECO_eta"]->Fill(tp.eta());	
	histos_th1f["h_NonAntiSTrack_RECO_phi"]->Fill(tp.phi());	
	histos_th2f["h2_NonAntiSTrack_RECO_vx_vy"]->Fill(tp.vx(),tp.vy());
	histos_th2f["h2_NonAntiSTrack_RECO_lxy_vz"]->Fill(Lxy,tp.vz());
	histos_th1f["h_NonAntiSTrack_RECO_lxy"]->Fill(Lxy);	
	histos_th1f["h_NonAntiSTrack_RECO_vz"]->Fill(tp.vz());	
	histos_th1f["h_NonAntiSTrack_RECO_dxy"]->Fill(dxy);
	histos_th1f["h_NonAntiSTrack_RECO_dz"]->Fill(dz);

	if(abs(tp.eta()) < 2.5 && Lxy < 3.5) histos_th1f["h_NonAntiSTrack_RECO_pt_cut_eta_Lxy"]->Fill(tp.pt());
	if(tp.pt() > 0.9 && abs(tp.eta()) < 2.5){

		 histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta"]->Fill(Lxy);
		 if(matchedTrackQuality == 0 || matchedTrackQuality == 1 || matchedTrackQuality == 2)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_purity_loose_tight_high"]->Fill(Lxy);
		 if(matchedTrackQuality == 1 || matchedTrackQuality == 2)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_purity_tight_high"]->Fill(Lxy);
		 if(abs(dxy) < 3.5 && abs(dz) < 30) histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_dxy_dz"]->Fill(Lxy);

		 if(nPVs == 1)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_dxy_nPVs_1"]->Fill(Lxy);
		 if(nPVs < 5)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_dxy_nPVs_smaller5"]->Fill(Lxy);
		 if(nPVs >= 5 && nPVs < 10)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_dxy_nPVs_from5to10"]->Fill(Lxy);
		 if(nPVs >= 10 && nPVs < 20)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_dxy_nPVs_from10to20"]->Fill(Lxy);
		 if(nPVs >= 20 && nPVs < 30)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_dxy_nPVs_from20to30"]->Fill(Lxy);
		 if(nPVs >= 30)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_pt_eta_dxy_nPVs_larger30"]->Fill(Lxy);

	}
	if(tp.pt() > 2. && abs(tp.eta()) < 1)histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_tight_pt_eta"]->Fill(Lxy);
	if(tp.pt() > 0.9 && Lxy < 3.5) histos_th1f["h_NonAntiSTrack_RECO_eta_cut_pt_lxy"]->Fill(tp.eta());

	//investigate inneficiency wrt vz
	if(tp.pt() > 0.9) histos_th1f["h_NonAntiSTrack_RECO_vz_cut_pt"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && abs(tp.eta()) < 1) histos_th1f["h_NonAntiSTrack_RECO_vz_cut_pt_eta"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && Lxy < 3.5) histos_th1f["h_NonAntiSTrack_RECO_vz_cut_pt_lxy"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && Lxy < 3.5 && abs(tp.pz()) > 2) histos_th1f["h_NonAntiSTrack_RECO_vz_cut_pt_lxy_pz"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && abs(dxy) < 0.1) histos_th1f["h_NonAntiSTrack_RECO_vz_cut_pt_dxy"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && Lxy < 3.5      && abs(tp.eta()) < 1 ) histos_th1f["h_NonAntiSTrack_RECO_vz_cut_pt_lxy_tight_eta"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && abs(dxy) < 0.1 && abs(tp.eta()) < 1 ) histos_th1f["h_NonAntiSTrack_RECO_vz_cut_pt_dxy_tight_eta"]->Fill(tp.vz());

	if(tp.pt() > 0.9 && abs(dxy) < 0.1 && Lxy <3.5) histos_th1f["h_NonAntiSTrack_RECO_vz_cut_pt_dxy_lxy"]->Fill(tp.vz());

	//to investigate inneficiency wrt vz check the difference between the trackcollections for |vz| < 20cm and |vz| > 20cm:
	if(abs(tp.vz()) < 20 && tp.pt() > 0.9){
		histos_th1f["h_NonAntiSTrack_RECO_pt_cut_vzSmall"]->Fill(tp.pt());
		histos_th1f["h_NonAntiSTrack_RECO_pz_cut_vzSmall"]->Fill(tp.pz());
		histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_vzSmall"]->Fill(Lxy);
		histos_th1f["h_NonAntiSTrack_RECO_dxy_cut_vzSmall"]->Fill(dxy);
		histos_th1f["h_NonAntiSTrack_RECO_dz_cut_vzSmall"]->Fill(dz);
		histos_th1f["h_NonAntiSTrack_RECO_vz_cut_vzSmall"]->Fill(tp.vz());
	}
	else if(abs(tp.vz()) > 20 && tp.pt() > 2. && abs(tp.pz()) > 5. && Lxy < 20){
		histos_th1f["h_NonAntiSTrack_RECO_pt_cut_vzLarge"]->Fill(tp.pt());
		histos_th1f["h_NonAntiSTrack_RECO_pz_cut_vzLarge"]->Fill(tp.pz());
		histos_th1f["h_NonAntiSTrack_RECO_lxy_cut_vzLarge"]->Fill(Lxy);
		histos_th1f["h_NonAntiSTrack_RECO_dxy_cut_vzLarge"]->Fill(dxy);
		histos_th1f["h_NonAntiSTrack_RECO_dz_cut_vzLarge"]->Fill(dz);
		histos_th1f["h_NonAntiSTrack_RECO_vz_cut_vzLarge"]->Fill(tp.vz());
	}

}

void FlatTreeProducerTracking::FillHistosNonAntiSTracksAll(const TrackingParticle& tp, TVector3 beamspot, int nPVs, int matchedTrackQuality){
	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	double dz = AnalyzerAllSteps::dz_line_point(tpCreationVertex,tpMomentum,beamspot);
	histos_th1f["h_NonAntiSTrack_All_pt"]->Fill(tp.pt());	
	histos_th1f["h_NonAntiSTrack_All_eta"]->Fill(tp.eta());	
	histos_th1f["h_NonAntiSTrack_All_phi"]->Fill(tp.phi());	
	histos_th2f["h2_NonAntiSTrack_All_vx_vy"]->Fill(tp.vx(),tp.vy());
	histos_th2f["h2_NonAntiSTrack_All_lxy_vz"]->Fill(Lxy,tp.vz());
	histos_th1f["h_NonAntiSTrack_All_lxy"]->Fill(Lxy);	
	histos_th1f["h_NonAntiSTrack_All_vz"]->Fill(tp.vz());	
	histos_th1f["h_NonAntiSTrack_All_dxy"]->Fill(dxy);
	histos_th1f["h_NonAntiSTrack_All_dz"]->Fill(dz);

	if(abs(tp.eta()) < 2.5 && Lxy < 3.5) histos_th1f["h_NonAntiSTrack_All_pt_cut_eta_Lxy"]->Fill(tp.pt());
	if(tp.pt() > 0.9 && abs(tp.eta()) < 2.5) {
  	
		 histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta"]->Fill(Lxy);
		 if(matchedTrackQuality == 0 || matchedTrackQuality == 1 || matchedTrackQuality == 2)histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_purity_loose_tight_high"]->Fill(Lxy);
		 if(matchedTrackQuality == 1 || matchedTrackQuality == 2)histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_purity_tight_high"]->Fill(Lxy);
		 if(abs(dxy) < 3.5 && abs(dz) < 30) histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_dxy_dz"]->Fill(Lxy);
       
		 if(nPVs == 1)histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_dxy_nPVs_1"]->Fill(Lxy);
		 if(nPVs < 5)histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_dxy_nPVs_smaller5"]->Fill(Lxy);
		 if(nPVs >= 5 && nPVs < 10)histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_dxy_nPVs_from5to10"]->Fill(Lxy);
		 if(nPVs >= 10 && nPVs < 20)histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_dxy_nPVs_from10to20"]->Fill(Lxy);
		 if(nPVs >= 20 && nPVs < 30)histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_dxy_nPVs_from20to30"]->Fill(Lxy);
		 if(nPVs >= 30)histos_th1f["h_NonAntiSTrack_All_lxy_cut_pt_eta_dxy_nPVs_larger30"]->Fill(Lxy);
 	}
	if(tp.pt() > 2. && abs(tp.eta()) < 1)histos_th1f["h_NonAntiSTrack_All_lxy_cut_tight_pt_eta"]->Fill(Lxy);
	if(tp.pt() > 0.9 && Lxy < 3.5) histos_th1f["h_NonAntiSTrack_All_eta_cut_pt_lxy"]->Fill(tp.eta());

	//investigate inneficiency wrt vz (to check the cut off at 20cm)
	if(tp.pt() > 0.9) histos_th1f["h_NonAntiSTrack_All_vz_cut_pt"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && abs(tp.eta()) < 1) histos_th1f["h_NonAntiSTrack_All_vz_cut_pt_eta"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && Lxy < 3.5) histos_th1f["h_NonAntiSTrack_All_vz_cut_pt_lxy"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && Lxy < 3.5 && abs(tp.pz()) > 2) histos_th1f["h_NonAntiSTrack_All_vz_cut_pt_lxy_pz"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && abs(dxy) < 0.1) histos_th1f["h_NonAntiSTrack_All_vz_cut_pt_dxy"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && Lxy < 3.5      && abs(tp.eta()) < 1 ) histos_th1f["h_NonAntiSTrack_All_vz_cut_pt_lxy_tight_eta"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && abs(dxy) < 0.1 && abs(tp.eta()) < 1 ) histos_th1f["h_NonAntiSTrack_All_vz_cut_pt_dxy_tight_eta"]->Fill(tp.vz());
	if(tp.pt() > 0.9 && abs(dxy) < 0.1 && Lxy < 3.5) histos_th1f["h_NonAntiSTrack_All_vz_cut_pt_dxy_lxy"]->Fill(tp.vz());
	
	//to investigate inneficiency wrt vz check the difference between the trackcollections for |vz| < 20cm and |vz| > 20cm:
	if(abs(tp.vz()) < 20 && tp.pt() > 0.9){
		histos_th1f["h_NonAntiSTrack_All_pt_cut_vzSmall"]->Fill(tp.pt());
		histos_th1f["h_NonAntiSTrack_All_pz_cut_vzSmall"]->Fill(tp.pz());
		histos_th1f["h_NonAntiSTrack_All_lxy_cut_vzSmall"]->Fill(Lxy);
		histos_th1f["h_NonAntiSTrack_All_dxy_cut_vzSmall"]->Fill(dxy);
		histos_th1f["h_NonAntiSTrack_All_dz_cut_vzSmall"]->Fill(dz);
		histos_th1f["h_NonAntiSTrack_All_vz_cut_vzSmall"]->Fill(tp.vz());
	}
	else if(abs(tp.vz()) > 20 && tp.pt() > 2. && abs(tp.pz()) > 5. && Lxy < 20){
		histos_th1f["h_NonAntiSTrack_All_pt_cut_vzLarge"]->Fill(tp.pt());
		histos_th1f["h_NonAntiSTrack_All_pz_cut_vzLarge"]->Fill(tp.pz());
		histos_th1f["h_NonAntiSTrack_All_lxy_cut_vzLarge"]->Fill(Lxy);
		histos_th1f["h_NonAntiSTrack_All_dxy_cut_vzLarge"]->Fill(dxy);
		histos_th1f["h_NonAntiSTrack_All_dz_cut_vzLarge"]->Fill(dz);
		histos_th1f["h_NonAntiSTrack_All_vz_cut_vzLarge"]->Fill(tp.vz());
	}
	

}


void FlatTreeProducerTracking::FillHistosAntiSTracks(const TrackingParticle& tp, TVector3 beamspot, TrackingParticleCollection const & TPColl, edm::Handle<TrackingParticleCollection> h_TP, edm::Handle< reco::TrackToTrackingParticleAssociator> h_trackAssociator, edm::Handle<View<reco::Track>> h_generalTracks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L){
	//now start from this tp and go down in gen particles to check if you find daughter and granddaughters
	vector<bool> granddaughterTrackMatched;
	granddaughterTrackMatched.push_back(false);granddaughterTrackMatched.push_back(false);granddaughterTrackMatched.push_back(false);granddaughterTrackMatched.push_back(false);
	const reco::Track *matchedTrackPointerKsPosPion = nullptr;
	const reco::Track *matchedTrackPointerKsNegPion = nullptr;
	const reco::Track *matchedTrackPointerAntiLPosPion = nullptr;
	const reco::Track *matchedTrackPointerAntiLNegProton = nullptr;
	
	int nCorrectGrandDaughters = 0;
	
	//get antiS decay vertex
	tv_iterator tp_firstDecayVertex = tp.decayVertices_begin();
	double tpdecayVx = (**tp_firstDecayVertex).position().X(); double tpdecayVy = (**tp_firstDecayVertex).position().Y(); double tpdecayVz = (**tp_firstDecayVertex).position().Z();
	//don't consider antiS which have an interaction vertex equal to their decay vertex because these are loopers, the return here is not strictly nessesary because loopers do not have daughters, so the below loops would not work, but with the return it is much faster of course 
	if(tp.vertex().X() == tpdecayVx && tp.vertex().Y() == tpdecayVy && tp.vertex().Z() == tpdecayVz){return;}
	//std::cout << "-----------------" << std::endl;
	//now loop over the TP and try to find daughters
	for(size_t j=0; j<TPColl.size(); ++j) {
		const TrackingParticle& tp_daughter = TPColl[j];
		if(abs(tp_daughter.pdgId()) == 310 || tp_daughter.pdgId() == -3122){//daughter has to be a Ks or Lambda
			
			double daughterVx = tp_daughter.vx(); double daughterVy = tp_daughter.vy(); double daughterVz = tp_daughter.vz();	
			tv_iterator tp_daughter_firstDecayVertex = tp_daughter.decayVertices_begin();
			double tp_daughter_decayVx = (**tp_daughter_firstDecayVertex).position().X(); double tp_daughter_decayVy = (**tp_daughter_firstDecayVertex).position().Y(); double tp_daughter_decayVz = (**tp_daughter_firstDecayVertex).position().Z();

			if( daughterVx == tpdecayVx && daughterVy == tpdecayVy && daughterVz == tpdecayVz){//daughter vertex has to match the mother decay vertex
				for(size_t k=0; k<TPColl.size(); ++k) {//loop to find the granddaughters
					const TrackingParticle& tp_granddaughter = TPColl[k];

					if(tp_granddaughter.vx() == tp_daughter_decayVx && tp_granddaughter.vy() == tp_daughter_decayVy && tp_granddaughter.vz() == tp_daughter_decayVz){
						if(tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdAntiProton || tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdPosPion ||tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdNegPion ){//found a tp that is a granddaughter of the antiS
							  //now check if you have matching track to this granddaughter
													
							  const reco::Track *matchedTrackPointer = nullptr;
							  TrackingParticleRef tpr(h_TP,k);
							  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
							  reco::SimToRecoCollection simRecCollL;
							  reco::SimToRecoCollection const * simRecCollP=nullptr;
							  simRecCollL = std::move(h_trackAssociator->associateSimToReco(h_generalTracks,h_TP));
							  simRecCollP = &simRecCollL;
							  reco::SimToRecoCollection const & simRecColl = *simRecCollP;

							  if(simRecColl.find(tpr) != simRecColl.end()){
								  auto const & rt = simRecColl[tpr];
								  if (rt.size()!=0) {
									    
									    matchedTrackPointer = rt.begin()->first.get();
									    //matchingTrackFound = true;
									    if(abs(tp_daughter.pdgId()) == 310 && tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdPosPion){granddaughterTrackMatched[0]=true;matchedTrackPointerKsPosPion= matchedTrackPointer; FillHistosAntiSKsDaughterTracksRECO(tp_granddaughter,beamspot);}
									    else if(abs(tp_daughter.pdgId()) == 310 && tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdNegPion){granddaughterTrackMatched[1]=true;matchedTrackPointerKsNegPion=matchedTrackPointer;FillHistosAntiSKsDaughterTracksRECO(tp_granddaughter,beamspot);}
									    else if(tp_daughter.pdgId() == -3122 && tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdAntiProton){granddaughterTrackMatched[2]=true;matchedTrackPointerAntiLPosPion=matchedTrackPointer; FillHistosAntiSAntiLAntiProtonDaughterTracksRECO(tp_granddaughter,beamspot);}
									    else if(tp_daughter.pdgId() == -3122 && tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdPosPion){granddaughterTrackMatched[3]=true;matchedTrackPointerAntiLNegProton=matchedTrackPointer;FillHistosAntiSAntiLPosPionDaughterTracksRECO(tp_granddaughter,beamspot);}
								  }
							  }
							  if(abs(tp_daughter.pdgId()) == 310 && tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdPosPion){nCorrectGrandDaughters++;FillHistosAntiSKsDaughterTracksAll(tp_granddaughter,beamspot);}
							  else if(abs(tp_daughter.pdgId()) == 310 && tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdNegPion){nCorrectGrandDaughters++;FillHistosAntiSKsDaughterTracksAll(tp_granddaughter,beamspot);}
							  else if(tp_daughter.pdgId() == -3122 && tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdAntiProton){nCorrectGrandDaughters++;FillHistosAntiSAntiLAntiProtonDaughterTracksAll(tp_granddaughter,beamspot);}
							  else if(tp_daughter.pdgId() == -3122 && tp_granddaughter.pdgId()==AnalyzerAllSteps::pdgIdPosPion){nCorrectGrandDaughters++;FillHistosAntiSAntiLPosPionDaughterTracksAll(tp_granddaughter,beamspot);}
							  //}else{
							  //	    	    matchingTrackFound = false;
							  //}
						}
					}
				}
			}	
		}	
    	}
  //std::cout << "Ks, pi+ " << granddaughterTrackMatched[0] << std::endl; 
  //std::cout << "Ks, pi- " << granddaughterTrackMatched[1] << std::endl; 
  //std::cout << " L, p-  " << granddaughterTrackMatched[2] << std::endl; 
  //std::cout << " L, pi+ " << granddaughterTrackMatched[3] << std::endl; 
  
  //std::cout << "number of correct (=charged) granddaughters of this antiS: " << nCorrectGrandDaughters << std::endl;
  if(nCorrectGrandDaughters == 4)FillMajorEfficiencyPlot(granddaughterTrackMatched, matchedTrackPointerKsPosPion,matchedTrackPointerKsNegPion,matchedTrackPointerAntiLPosPion,matchedTrackPointerAntiLNegProton,h_V0Ks,h_V0L);

}

void FlatTreeProducerTracking::FillHistosAntiSKsDaughterTracksRECO(const TrackingParticle& tp, TVector3 beamspot){
	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	histos_th1f["h_AntiSKsDaughterTracks_RECO_pt"]->Fill(tp.pt());	
	histos_th1f["h_AntiSKsDaughterTracks_RECO_eta"]->Fill(tp.eta());	
	histos_th1f["h_AntiSKsDaughterTracks_RECO_phi"]->Fill(tp.phi());	

	histos_th2f["h2_AntiSKsDaughterTracks_RECO_vx_vy"]->Fill(tp.vx(),tp.vy());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSKsDaughterTracks_RECO_vx_vy_cut_pt_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSKsDaughterTracks_RECO_vx_vy_cut_pt_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>1)histos_th2f["h2_AntiSKsDaughterTracks_RECO_vx_vy_cut_pt_1"]->Fill(tp.vx(),tp.vy());

	if(tp.pz()<0.5)histos_th2f["h2_AntiSKsDaughterTracks_RECO_vx_vy_cut_pz_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSKsDaughterTracks_RECO_vx_vy_cut_pz_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>1)histos_th2f["h2_AntiSKsDaughterTracks_RECO_vx_vy_cut_pz_1"]->Fill(tp.vx(),tp.vy());

	histos_th2f["h2_AntiSKsDaughterTracks_RECO_lxy_vz"]->Fill(Lxy,tp.vz());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSKsDaughterTracks_RECO_lxy_vz_cut_pt_0p5"]->Fill(Lxy,tp.vz());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSKsDaughterTracks_RECO_lxy_vz_cut_pt_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(tp.pt()>1)histos_th2f["h2_AntiSKsDaughterTracks_RECO_lxy_vz_cut_pt_1"]->Fill(Lxy,tp.vz());

	if(abs(tp.pz())<0.5)histos_th2f["h2_AntiSKsDaughterTracks_RECO_lxy_vz_cut_pz_0p5"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>0.5 && abs(tp.pz())<1)histos_th2f["h2_AntiSKsDaughterTracks_RECO_lxy_vz_cut_pz_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>1)histos_th2f["h2_AntiSKsDaughterTracks_RECO_lxy_vz_cut_pz_1"]->Fill(Lxy,tp.vz());

	histos_th1f["h_AntiSKsDaughterTracks_RECO_lxy"]->Fill(Lxy);	
	histos_th1f["h_AntiSKsDaughterTracks_RECO_vz"]->Fill(tp.vz());	
	histos_th1f["h_AntiSKsDaughterTracks_RECO_dxy"]->Fill(dxy);
	
}
void FlatTreeProducerTracking::FillHistosAntiSKsDaughterTracksAll(const TrackingParticle& tp, TVector3 beamspot){
	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	histos_th1f["h_AntiSKsDaughterTracks_All_pt"]->Fill(tp.pt());	
	histos_th1f["h_AntiSKsDaughterTracks_All_eta"]->Fill(tp.eta());	
	histos_th1f["h_AntiSKsDaughterTracks_All_phi"]->Fill(tp.phi());	

	histos_th2f["h2_AntiSKsDaughterTracks_All_vx_vy"]->Fill(tp.vx(),tp.vy());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSKsDaughterTracks_All_vx_vy_cut_pt_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSKsDaughterTracks_All_vx_vy_cut_pt_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>1)histos_th2f["h2_AntiSKsDaughterTracks_All_vx_vy_cut_pt_1"]->Fill(tp.vx(),tp.vy());

	if(tp.pz()<0.5)histos_th2f["h2_AntiSKsDaughterTracks_All_vx_vy_cut_pz_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSKsDaughterTracks_All_vx_vy_cut_pz_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>1)histos_th2f["h2_AntiSKsDaughterTracks_All_vx_vy_cut_pz_1"]->Fill(tp.vx(),tp.vy());

	histos_th2f["h2_AntiSKsDaughterTracks_All_lxy_vz"]->Fill(Lxy,tp.vz());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSKsDaughterTracks_All_lxy_vz_cut_pt_0p5"]->Fill(Lxy,tp.vz());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSKsDaughterTracks_All_lxy_vz_cut_pt_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(tp.pt()>1)histos_th2f["h2_AntiSKsDaughterTracks_All_lxy_vz_cut_pt_1"]->Fill(Lxy,tp.vz());

	if(abs(tp.pz())<0.5)histos_th2f["h2_AntiSKsDaughterTracks_All_lxy_vz_cut_pz_0p5"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>0.5 && abs(tp.pz())<1)histos_th2f["h2_AntiSKsDaughterTracks_All_lxy_vz_cut_pz_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>1)histos_th2f["h2_AntiSKsDaughterTracks_All_lxy_vz_cut_pz_1"]->Fill(Lxy,tp.vz());

	histos_th1f["h_AntiSKsDaughterTracks_All_lxy"]->Fill(Lxy);	
	histos_th1f["h_AntiSKsDaughterTracks_All_vz"]->Fill(tp.vz());	
	histos_th1f["h_AntiSKsDaughterTracks_All_dxy"]->Fill(dxy);
	
}

void FlatTreeProducerTracking::FillHistosAntiSAntiLAntiProtonDaughterTracksRECO(const TrackingParticle& tp, TVector3 beamspot){
	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_RECO_pt"]->Fill(tp.pt());	
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_RECO_eta"]->Fill(tp.eta());	
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_RECO_phi"]->Fill(tp.phi());	

	histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_vx_vy"]->Fill(tp.vx(),tp.vy());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_vx_vy_cut_pt_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_vx_vy_cut_pt_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_vx_vy_cut_pt_1"]->Fill(tp.vx(),tp.vy());

	if(tp.pz()<0.5)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_vx_vy_cut_pz_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_vx_vy_cut_pz_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_vx_vy_cut_pz_1"]->Fill(tp.vx(),tp.vy());

	histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_lxy_vz"]->Fill(Lxy,tp.vz());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_lxy_vz_cut_pt_0p5"]->Fill(Lxy,tp.vz());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_lxy_vz_cut_pt_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(tp.pt()>1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_lxy_vz_cut_pt_1"]->Fill(Lxy,tp.vz());

	if(abs(tp.pz())<0.5)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_lxy_vz_cut_pz_0p5"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>0.5 && abs(tp.pz())<1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_lxy_vz_cut_pz_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_RECO_lxy_vz_cut_pz_1"]->Fill(Lxy,tp.vz());


	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_RECO_lxy"]->Fill(Lxy);	
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_RECO_vz"]->Fill(tp.vz());	
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_RECO_dxy"]->Fill(dxy);
	
}
void FlatTreeProducerTracking::FillHistosAntiSAntiLAntiProtonDaughterTracksAll(const TrackingParticle& tp, TVector3 beamspot){
	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_All_pt"]->Fill(tp.pt());	
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_All_eta"]->Fill(tp.eta());	
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_All_phi"]->Fill(tp.phi());	

	histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_vx_vy"]->Fill(tp.vx(),tp.vy());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_vx_vy_cut_pt_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_vx_vy_cut_pt_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_vx_vy_cut_pt_1"]->Fill(tp.vx(),tp.vy());

	if(tp.pz()<0.5)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_vx_vy_cut_pz_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_vx_vy_cut_pz_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_vx_vy_cut_pz_1"]->Fill(tp.vx(),tp.vy());

	histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_lxy_vz"]->Fill(Lxy,tp.vz());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_lxy_vz_cut_pt_0p5"]->Fill(Lxy,tp.vz());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_lxy_vz_cut_pt_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(tp.pt()>1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_lxy_vz_cut_pt_1"]->Fill(Lxy,tp.vz());

	if(abs(tp.pz())<0.5)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_lxy_vz_cut_pz_0p5"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>0.5 && abs(tp.pz())<1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_lxy_vz_cut_pz_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>1)histos_th2f["h2_AntiSAntiLAntiProtonDaughterTracks_All_lxy_vz_cut_pz_1"]->Fill(Lxy,tp.vz());



	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_All_lxy"]->Fill(Lxy);	
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_All_vz"]->Fill(tp.vz());	
	histos_th1f["h_AntiSAntiLAntiProtonDaughterTracks_All_dxy"]->Fill(dxy);
	
}

void FlatTreeProducerTracking::FillHistosAntiSAntiLPosPionDaughterTracksRECO(const TrackingParticle& tp, TVector3 beamspot){
	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_RECO_pt"]->Fill(tp.pt());	
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_RECO_eta"]->Fill(tp.eta());	
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_RECO_phi"]->Fill(tp.phi());	

	histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_vx_vy"]->Fill(tp.vx(),tp.vy());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_vx_vy_cut_pt_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_vx_vy_cut_pt_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_vx_vy_cut_pt_1"]->Fill(tp.vx(),tp.vy());

	if(tp.pz()<0.5)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_vx_vy_cut_pz_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_vx_vy_cut_pz_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_vx_vy_cut_pz_1"]->Fill(tp.vx(),tp.vy());

	histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_lxy_vz"]->Fill(Lxy,tp.vz());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_lxy_vz_cut_pt_0p5"]->Fill(Lxy,tp.vz());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_lxy_vz_cut_pt_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(tp.pt()>1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_lxy_vz_cut_pt_1"]->Fill(Lxy,tp.vz());

	if(abs(tp.pz())<0.5)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_lxy_vz_cut_pz_0p5"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>0.5 && abs(tp.pz())<1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_lxy_vz_cut_pz_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_RECO_lxy_vz_cut_pz_1"]->Fill(Lxy,tp.vz());

	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_RECO_lxy"]->Fill(Lxy);	
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_RECO_vz"]->Fill(tp.vz());	
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_RECO_dxy"]->Fill(dxy);
	
}
void FlatTreeProducerTracking::FillHistosAntiSAntiLPosPionDaughterTracksAll(const TrackingParticle& tp, TVector3 beamspot){
	TVector3 tpCreationVertex(tp.vx(),tp.vy(),tp.vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(tp.px(),tp.py(),tp.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_All_pt"]->Fill(tp.pt());	
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_All_eta"]->Fill(tp.eta());	
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_All_phi"]->Fill(tp.phi());	

	histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_vx_vy"]->Fill(tp.vx(),tp.vy());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_vx_vy_cut_pt_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_vx_vy_cut_pt_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pt()>1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_vx_vy_cut_pt_1"]->Fill(tp.vx(),tp.vy());

	if(tp.pz()<0.5)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_vx_vy_cut_pz_0p5"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_vx_vy_cut_pz_0p5_and_1"]->Fill(tp.vx(),tp.vy());
	if(tp.pz()>1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_vx_vy_cut_pz_1"]->Fill(tp.vx(),tp.vy());

	histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_lxy_vz"]->Fill(Lxy,tp.vz());

	if(tp.pt()<0.5)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_lxy_vz_cut_pt_0p5"]->Fill(Lxy,tp.vz());
	if(tp.pt()>0.5 && tp.pt()<1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_lxy_vz_cut_pt_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(tp.pt()>1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_lxy_vz_cut_pt_1"]->Fill(Lxy,tp.vz());

	if(abs(tp.pz())<0.5)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_lxy_vz_cut_pz_0p5"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>0.5 && abs(tp.pz())<1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_lxy_vz_cut_pz_0p5_and_1"]->Fill(Lxy,tp.vz());
	if(abs(tp.pz())>1)histos_th2f["h2_AntiSAntiLPosPionDaughterTracks_All_lxy_vz_cut_pz_1"]->Fill(Lxy,tp.vz());


	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_All_lxy"]->Fill(Lxy);	
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_All_vz"]->Fill(tp.vz());	
	histos_th1f["h_AntiSAntiLPosPionDaughterTracks_All_dxy"]->Fill(dxy);
	
}


void FlatTreeProducerTracking::FillMajorEfficiencyPlot(std::vector<bool>granddaughterTrackMatched, const reco::Track *matchedTrackPointerKsPosPion,const reco::Track *matchedTrackPointerKsNegPion,const reco::Track *matchedTrackPointerAntiLPosPion,const reco::Track *matchedTrackPointerAntiLNegProton, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L){
	//if the Ks daughters were found calculate the sum of the daughters momenta and try to find a reco Ks which matches this in deltaR
	bool matchingV0KsFound = false;
	if(matchedTrackPointerKsPosPion && matchedTrackPointerKsNegPion){
		  math::XYZVector sumV0Momenta = matchedTrackPointerKsPosPion->momentum() + matchedTrackPointerKsNegPion->momentum();
		  double deltaRmin = 999;
		  if(h_V0Ks.isValid()){
		      for(unsigned int i = 0; i < h_V0Ks->size(); ++i){//loop all Ks candidates
			const reco::VertexCompositeCandidate * V0 = &h_V0Ks->at(i);	
			double deltaPhi = reco::deltaPhi(sumV0Momenta.phi(),V0->phi());
			double deltaEta = sumV0Momenta.eta()-V0->eta();
			double deltaR = sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
			if(deltaR < deltaRmin) deltaRmin = deltaR;
		      }
		      if(deltaRmin<AnalyzerAllSteps::deltaRCutV0RECOKs)matchingV0KsFound = true;
		  }
		  histos_th1f["h_deltaRmin_V0Ks_momentumSumKsDaughterTracks"]->Fill(deltaRmin);
	}

	//if the AntiL daughters were found calculate the sum of the daughters momenta and try to find a reco AntiL which matches this in deltaR
	bool matchingV0AntiLFound = false;
	if(matchedTrackPointerAntiLPosPion && matchedTrackPointerAntiLNegProton){
		  math::XYZVector sumV0Momenta = matchedTrackPointerAntiLPosPion->momentum() + matchedTrackPointerAntiLNegProton->momentum();
		  double deltaRmin = 999;
		  if(h_V0L.isValid()){
		      for(unsigned int i = 0; i < h_V0L->size(); ++i){//loop all Ks candidates
			const reco::VertexCompositeCandidate * V0 = &h_V0L->at(i);	
			double deltaPhi = reco::deltaPhi(sumV0Momenta.phi(),V0->phi());
			double deltaEta = sumV0Momenta.eta()-V0->eta();
			double deltaR = sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
			if(deltaR < deltaRmin) deltaRmin = deltaR;
		      }
		      if(deltaRmin<AnalyzerAllSteps::deltaRCutV0RECOLambda)matchingV0AntiLFound = true;
		  }
		histos_th1f["h_deltaRmin_V0AntiL_momentumSumAntiLDaughterTracks"]->Fill(deltaRmin);
	}

	if(!granddaughterTrackMatched[0] && !granddaughterTrackMatched[1] && !granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(0.,0.);
	if(granddaughterTrackMatched[0] && !granddaughterTrackMatched[1] && !granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(0.,1.);
	if(!granddaughterTrackMatched[0] && granddaughterTrackMatched[1] && !granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(0.,2.);
	if(granddaughterTrackMatched[0] && granddaughterTrackMatched[1] && !granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(0.,3.);

	if(!granddaughterTrackMatched[0] && !granddaughterTrackMatched[1] && granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(1.,0.);
	if(granddaughterTrackMatched[0]  && !granddaughterTrackMatched[1] && granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(1.,1.);
	if(!granddaughterTrackMatched[0] && granddaughterTrackMatched[1] && granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(1.,2.);
	if(granddaughterTrackMatched[0] && granddaughterTrackMatched[1] && granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(1.,3.);

	if(!granddaughterTrackMatched[0] && !granddaughterTrackMatched[1] && !granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(2.,0.);
	if(granddaughterTrackMatched[0] && !granddaughterTrackMatched[1] && !granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(2.,1.);
	if(!granddaughterTrackMatched[0] && granddaughterTrackMatched[1] && !granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(2.,2.);
	if(granddaughterTrackMatched[0] && granddaughterTrackMatched[1] && !granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(2.,3.);

	if(!granddaughterTrackMatched[0] && !granddaughterTrackMatched[1] && granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(3.,0.);
	if(granddaughterTrackMatched[0] && !granddaughterTrackMatched[1] && granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(3.,1.);
	if(!granddaughterTrackMatched[0] && granddaughterTrackMatched[1] && granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(3.,2.);
	if(granddaughterTrackMatched[0] && granddaughterTrackMatched[1] && granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(3.,3.);

	if(matchingV0KsFound && !matchingV0AntiLFound){//fill the top row of the matrix
		if( !granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(0.,4.);
		if(granddaughterTrackMatched[2] && !granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(1.,4.);
		if(!granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(2.,4.);
		if(granddaughterTrackMatched[2] && granddaughterTrackMatched[3])histos_th2i["h2_GlobalEfficiencies"]->Fill(3.,4.);
	}
	if(!matchingV0KsFound && matchingV0AntiLFound){//fill the right column of the matrix
		if( !granddaughterTrackMatched[0] && !granddaughterTrackMatched[1])histos_th2i["h2_GlobalEfficiencies"]->Fill(4.,0.);
		if(granddaughterTrackMatched[0] && !granddaughterTrackMatched[1])histos_th2i["h2_GlobalEfficiencies"]->Fill(4.,1.);
		if(!granddaughterTrackMatched[0] && granddaughterTrackMatched[1])histos_th2i["h2_GlobalEfficiencies"]->Fill(4.,2.);
		if(granddaughterTrackMatched[0] && granddaughterTrackMatched[1])histos_th2i["h2_GlobalEfficiencies"]->Fill(4.,3.);

	}
	if(matchingV0KsFound && matchingV0AntiLFound){//this means you found both a RECO V0-Ks and RECO V0-AntiL
		histos_th2i["h2_GlobalEfficiencies"]->Fill(4.,4.);
	}

	

}
*/

/*
void FlatTreeProducerTracking::FillHistosGENAntiS(const reco::Candidate  * GENantiS, TVector3 beamspot){

	TVector3 momentumGENantiS(GENantiS->px(),GENantiS->py(),GENantiS->pz());
	histos_th1f["h_GEN_AntiS_pt"]->Fill(GENantiS->pt());	
	histos_th1f["h_GEN_AntiS_p"]->Fill(GENantiS->p());	
	histos_th1f["h_GEN_AntiS_energy"]->Fill(GENantiS->energy());	
	histos_th1f["h_GEN_AntiS_eta"]->Fill(GENantiS->eta());	
	histos_th1f["h_GEN_AntiS_phi"]->Fill(GENantiS->phi());
		
	if(GENantiS->numberOfDaughters()==2){
		TVector3 momentumKs(GENantiS->daughter(0)->px(),GENantiS->daughter(0)->py(),GENantiS->daughter(0)->pz());
		TVector3 momentumAntiLambda(GENantiS->daughter(1)->px(),GENantiS->daughter(1)->py(),GENantiS->daughter(1)->pz());
		if(abs(GENantiS->daughter(1)->pdgId())==AnalyzerAllSteps::pdgIdKs){momentumKs.SetX(GENantiS->daughter(1)->px());momentumKs.SetY(GENantiS->daughter(1)->py());momentumKs.SetZ(GENantiS->daughter(1)->pz());}
		if(GENantiS->daughter(0)->pdgId()==AnalyzerAllSteps::pdgIdAntiLambda){momentumAntiLambda.SetX(GENantiS->daughter(0)->px());momentumAntiLambda.SetY(GENantiS->daughter(0)->py());momentumAntiLambda.SetZ(GENantiS->daughter(0)->pz());}

		double deltaPhiDaughters = reco::deltaPhi(GENantiS->daughter(0)->phi(),GENantiS->daughter(1)->phi());
		double deltaEtaDaughters = GENantiS->daughter(0)->eta()-GENantiS->daughter(1)->eta();
		double deltaRDaughters = sqrt(deltaPhiDaughters*deltaPhiDaughters+deltaEtaDaughters*deltaEtaDaughters);
		TVector3 interactionVertexAntiS(GENantiS->daughter(0)->vx(),GENantiS->daughter(0)->vy(),GENantiS->daughter(0)->vz());
		double lxyInteractionVertexAntiS = AnalyzerAllSteps::lxy(beamspot, interactionVertexAntiS);
		double lxyzInteractionVertexAntiS = AnalyzerAllSteps::lxyz(beamspot, interactionVertexAntiS);
		histos_th1f["h_GEN_AntiS_deltaR_daughters"]->Fill(deltaRDaughters);
		histos_th1f["h_GEN_AntiS_deltaPhi_daughters"]->Fill(deltaPhiDaughters);
		histos_th1f["h_GEN_AntiS_lxy_interactioVertex"]->Fill(lxyInteractionVertexAntiS);
		histos_th1f["h_GEN_AntiS_lxyz_interactioVertex"]->Fill(lxyzInteractionVertexAntiS);
		histos_th2f["h2_GEN_AntiS_deltaPhi_daughters_lxy_interactionVertex"]->Fill(deltaPhiDaughters,lxyInteractionVertexAntiS);
		histos_th1f["h_GEN_AntiS_AnglePAntiSPKs"]->Fill(TMath::ACos(AnalyzerAllSteps::CosOpeningsAngle(momentumKs,momentumGENantiS)));
		histos_th2f["h2_GEN_AntiS_PAntiS_AnglePAntiSPKs"]->Fill(GENantiS->p(),TMath::ACos(AnalyzerAllSteps::CosOpeningsAngle(momentumKs,momentumGENantiS)));
		histos_th1f["h_GEN_AntiS_AnglePAntiSPAntiLambda"]->Fill(TMath::ACos(AnalyzerAllSteps::CosOpeningsAngle(momentumAntiLambda,momentumGENantiS)));
		histos_th2f["h2_GEN_AntiS_PAntiS_AnglePAntiSPAntiLambda"]->Fill(GENantiS->p(),TMath::ACos(AnalyzerAllSteps::CosOpeningsAngle(momentumAntiLambda,momentumGENantiS)));
		histos_th2f["h2_GEN_AntiS_PointingAngleKs_PointingAngleAntiLambda"]->Fill(AnalyzerAllSteps::XYpointingAngle(GENantiS->daughter(0),beamspot),AnalyzerAllSteps::XYpointingAngle(GENantiS->daughter(1),beamspot));
	   	
	}	

}
*/

/*
void FlatTreeProducerTracking::FillHistosGENInteractingAntiS(const reco::Candidate  * GENantiS, TVector3 beamspot){

	TVector3 momentumGENantiS(GENantiS->px(),GENantiS->py(),GENantiS->pz());
	histos_th1f["h_GEN_AntiSInteract_pt"]->Fill(GENantiS->pt());	
	histos_th1f["h_GEN_AntiSInteract_p"]->Fill(GENantiS->p());	
	histos_th1f["h_GEN_AntiSInteract_energy"]->Fill(GENantiS->energy());	
	histos_th1f["h_GEN_AntiSInteract_eta"]->Fill(GENantiS->eta());	
	histos_th1f["h_GEN_AntiSInteract_phi"]->Fill(GENantiS->phi());
		
	if(GENantiS->numberOfDaughters()==2){
		TVector3 momentumKs(GENantiS->daughter(0)->px(),GENantiS->daughter(0)->py(),GENantiS->daughter(0)->pz());
		TVector3 momentumAntiLambda(GENantiS->daughter(1)->px(),GENantiS->daughter(1)->py(),GENantiS->daughter(1)->pz());
		if(abs(GENantiS->daughter(1)->pdgId())==AnalyzerAllSteps::pdgIdKs){momentumKs.SetX(GENantiS->daughter(1)->px());momentumKs.SetY(GENantiS->daughter(1)->py());momentumKs.SetZ(GENantiS->daughter(1)->pz());}
		if(GENantiS->daughter(0)->pdgId()==AnalyzerAllSteps::pdgIdAntiLambda){momentumAntiLambda.SetX(GENantiS->daughter(0)->px());momentumAntiLambda.SetY(GENantiS->daughter(0)->py());momentumAntiLambda.SetZ(GENantiS->daughter(0)->pz());}

		double deltaPhiDaughters = reco::deltaPhi(GENantiS->daughter(0)->phi(),GENantiS->daughter(1)->phi());
		double deltaEtaDaughters = GENantiS->daughter(0)->eta()-GENantiS->daughter(1)->eta();
		double deltaRDaughters = sqrt(deltaPhiDaughters*deltaPhiDaughters+deltaEtaDaughters*deltaEtaDaughters);
		TVector3 interactionVertexAntiS(GENantiS->daughter(0)->vx(),GENantiS->daughter(0)->vy(),GENantiS->daughter(0)->vz());
		double lxyInteractionVertexAntiS = AnalyzerAllSteps::lxy(beamspot, interactionVertexAntiS);
		double lxyzInteractionVertexAntiS = AnalyzerAllSteps::lxyz(beamspot, interactionVertexAntiS);



		histos_th1f["h_GEN_AntiSInteract_deltaR_daughters"]->Fill(deltaRDaughters);
		histos_th1f["h_GEN_AntiSInteract_deltaPhi_daughters"]->Fill(deltaPhiDaughters);
		histos_th1f["h_GEN_AntiSInteract_lxy_interactioVertex"]->Fill(lxyInteractionVertexAntiS);
		histos_th1f["h_GEN_AntiSInteract_lxyz_interactioVertex"]->Fill(lxyzInteractionVertexAntiS);
		histos_th2f["h2_GEN_AntiSInteract_deltaPhi_daughters_lxy_interactionVertex"]->Fill(deltaPhiDaughters,lxyInteractionVertexAntiS);
		histos_th1f["h_GEN_AntiSInteract_AnglePAntiSPKs"]->Fill(TMath::ACos(AnalyzerAllSteps::CosOpeningsAngle(momentumKs,momentumGENantiS)));
		histos_th2f["h2_GEN_AntiSInteract_PAntiS_AnglePAntiSPKs"]->Fill(GENantiS->p(),TMath::ACos(AnalyzerAllSteps::CosOpeningsAngle(momentumKs,momentumGENantiS)));
		histos_th1f["h_GEN_AntiSInteract_AnglePAntiSPAntiLambda"]->Fill(TMath::ACos(AnalyzerAllSteps::CosOpeningsAngle(momentumAntiLambda,momentumGENantiS)));
		histos_th2f["h2_GEN_AntiSInteract_PAntiS_AnglePAntiSPAntiLambda"]->Fill(GENantiS->p(),TMath::ACos(AnalyzerAllSteps::CosOpeningsAngle(momentumAntiLambda,momentumGENantiS)));
		histos_th2f["h2_GEN_AntiSInteract_PointingAngleKs_PointingAngleAntiLambda"]->Fill(AnalyzerAllSteps::XYpointingAngle(GENantiS->daughter(0),beamspot),AnalyzerAllSteps::XYpointingAngle(GENantiS->daughter(1),beamspot));
	   	
	}	

}
*/


/*
void FlatTreeProducerTracking::FillHistosGENKsNonAntiS(const reco::Candidate  * GENKsNonAntiS, TVector3 beamspot){
	TVector3 KsCreationVertex(GENKsNonAntiS->vx(),GENKsNonAntiS->vy(),GENKsNonAntiS->vz());
	
	double Lxy = AnalyzerAllSteps::lxy(beamspot,KsCreationVertex);
	TVector3 KsMomentum(GENKsNonAntiS->px(),GENKsNonAntiS->py(),GENKsNonAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,beamspot);
	double xypointingAngle = AnalyzerAllSteps::XYpointingAngle(GENKsNonAntiS,beamspot);

	histos_th1f["h_GEN_KsNonAntiS_pt"]->Fill(GENKsNonAntiS->pt());	
	histos_th1f["h_GEN_KsNonAntiS_eta"]->Fill(GENKsNonAntiS->eta());	
	histos_th1f["h_GEN_KsNonAntiS_phi"]->Fill(GENKsNonAntiS->phi());	
	histos_th2f["h2_GEN_KsNonAntiS_vx_vy"]->Fill(GENKsNonAntiS->vx(),GENKsNonAntiS->vy());
	histos_th2f["h2_GEN_KsNonAntiS_lxy_vz"]->Fill(Lxy,GENKsNonAntiS->vz());
	histos_th1f["h_GEN_KsNonAntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_GEN_KsNonAntiS_vz"]->Fill(GENKsNonAntiS->vz());	
	histos_th1f["h_GEN_KsNonAntiS_dxy"]->Fill(dxy);	
	if(GENKsNonAntiS->numberOfDaughters() == 2)histos_th1f["h_GEN_KsNonAntiS_XYpointingAngle"]->Fill(xypointingAngle);	
}
*/


/*
void FlatTreeProducerTracking::FillHistosGENAntiLambdaNonAntiS(const reco::Candidate  * GENAntiLambdaNonAntiS, TVector3 beamspot){
	TVector3 AntiLambdaCreationVertex(GENAntiLambdaNonAntiS->vx(),GENAntiLambdaNonAntiS->vy(),GENAntiLambdaNonAntiS->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,AntiLambdaCreationVertex);
	TVector3 AntiLambdaMomentum(GENAntiLambdaNonAntiS->px(),GENAntiLambdaNonAntiS->py(),GENAntiLambdaNonAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(AntiLambdaCreationVertex,AntiLambdaMomentum,beamspot);
	double xypointingAngle = AnalyzerAllSteps::XYpointingAngle(GENAntiLambdaNonAntiS,beamspot);

	histos_th1f["h_GEN_AntiLambdaNonAntiS_pt"]->Fill(GENAntiLambdaNonAntiS->pt());	
	histos_th1f["h_GEN_AntiLambdaNonAntiS_eta"]->Fill(GENAntiLambdaNonAntiS->eta());	
	histos_th1f["h_GEN_AntiLambdaNonAntiS_phi"]->Fill(GENAntiLambdaNonAntiS->phi());	
	histos_th2f["h2_GEN_AntiLambdaNonAntiS_vx_vy"]->Fill(GENAntiLambdaNonAntiS->vx(),GENAntiLambdaNonAntiS->vy());
	histos_th2f["h2_GEN_AntiLambdaNonAntiS_lxy_vz"]->Fill(Lxy,GENAntiLambdaNonAntiS->vz());
	histos_th1f["h_GEN_AntiLambdaNonAntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_GEN_AntiLambdaNonAntiS_vz"]->Fill(GENAntiLambdaNonAntiS->vz());	
	histos_th1f["h_GEN_AntiLambdaNonAntiS_dxy"]->Fill(dxy);	
	if(GENAntiLambdaNonAntiS->numberOfDaughters() == 2)histos_th1f["h_GEN_AntiLambdaNonAntiS_XYpointingAngle"]->Fill(xypointingAngle);	
}
*/

/*
void FlatTreeProducerTracking::FillHistosGENKsAntiS(const reco::Candidate  * GENKsAntiS, TVector3 beamspot){
	TVector3 KsAntiSCreationVertex(GENKsAntiS->vx(),GENKsAntiS->vy(),GENKsAntiS->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,KsAntiSCreationVertex);
	double Lxyz = AnalyzerAllSteps::lxyz(beamspot,KsAntiSCreationVertex);
	TVector3 KsAntiSMomentum(GENKsAntiS->px(),GENKsAntiS->py(),GENKsAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(KsAntiSCreationVertex,KsAntiSMomentum,beamspot);
	double dz = AnalyzerAllSteps::dz_line_point(KsAntiSCreationVertex,KsAntiSMomentum,beamspot);
	double xypointingAngle = AnalyzerAllSteps::XYpointingAngle(GENKsAntiS,beamspot);

	histos_th1f["h_GEN_KsAntiS_pdgId"]->Fill(GENKsAntiS->pdgId());	
	histos_th1f["h_GEN_KsAntiS_pt"]->Fill(GENKsAntiS->pt());	
	histos_th1f["h_GEN_KsAntiS_pz"]->Fill(GENKsAntiS->pz());	
	histos_th2f["h2_GEN_KsAntiS_pt_pz"]->Fill(GENKsAntiS->pt(),GENKsAntiS->pz());	
	histos_th1f["h_GEN_KsAntiS_eta"]->Fill(GENKsAntiS->eta());	
	histos_th1f["h_GEN_KsAntiS_phi"]->Fill(GENKsAntiS->phi());	
	histos_th2f["h2_GEN_KsAntiS_vx_vy"]->Fill(GENKsAntiS->vx(),GENKsAntiS->vy());
	histos_th2f["h2_GEN_KsAntiS_lxy_vz"]->Fill(Lxy,GENKsAntiS->vz());
	histos_th1f["h_GEN_KsAntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_GEN_KsAntiS_vz"]->Fill(GENKsAntiS->vz());	
	histos_th1f["h_GEN_KsAntiS_lxyz"]->Fill(Lxyz);	
	histos_th1f["h_GEN_KsAntiS_dxy"]->Fill(dxy);
	histos_th1f["h_GEN_KsAntiS_dz"]->Fill(dz);
	histos_th2f["h2_GEN_KsAntiS_dxy_dz"]->Fill(dxy,dz);
	if(GENKsAntiS->numberOfDaughters() == 2)histos_th1f["h_GEN_KsAntiS_XYpointingAngle"]->Fill(xypointingAngle);
		
}
*/

/*
void FlatTreeProducerTracking::FillHistosGENAntiLambdaAntiS(const reco::Candidate  * GENAntiLambdaAntiS, TVector3 beamspot){
	TVector3 AntiLambdaAntiSCreationVertex(GENAntiLambdaAntiS->vx(),GENAntiLambdaAntiS->vy(),GENAntiLambdaAntiS->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,AntiLambdaAntiSCreationVertex);
	double Lxyz = AnalyzerAllSteps::lxyz(beamspot,AntiLambdaAntiSCreationVertex);
	TVector3 AntiLambdaMomentum(GENAntiLambdaAntiS->px(),GENAntiLambdaAntiS->py(),GENAntiLambdaAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(AntiLambdaAntiSCreationVertex,AntiLambdaMomentum,beamspot);
	double dz = AnalyzerAllSteps::dz_line_point(AntiLambdaAntiSCreationVertex,AntiLambdaMomentum,beamspot);
	double xypointingAngle = AnalyzerAllSteps::XYpointingAngle(GENAntiLambdaAntiS,beamspot);

	histos_th1f["h_GEN_AntiLambdaAntiS_pdgId"]->Fill(GENAntiLambdaAntiS->pdgId());	
	histos_th1f["h_GEN_AntiLambdaAntiS_pt"]->Fill(GENAntiLambdaAntiS->pt());	
	histos_th1f["h_GEN_AntiLambdaAntiS_pz"]->Fill(GENAntiLambdaAntiS->pz());	
	histos_th2f["h2_GEN_AntiLambdaAntiS_pt_pz"]->Fill(GENAntiLambdaAntiS->pt(),GENAntiLambdaAntiS->pz());	
	histos_th1f["h_GEN_AntiLambdaAntiS_eta"]->Fill(GENAntiLambdaAntiS->eta());	
	histos_th1f["h_GEN_AntiLambdaAntiS_phi"]->Fill(GENAntiLambdaAntiS->phi());	
	histos_th2f["h2_GEN_AntiLambdaAntiS_vx_vy"]->Fill(GENAntiLambdaAntiS->vx(),GENAntiLambdaAntiS->vy());
	histos_th2f["h2_GEN_AntiLambdaAntiS_lxy_vz"]->Fill(Lxy,GENAntiLambdaAntiS->vz());
	histos_th1f["h_GEN_AntiLambdaAntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_GEN_AntiLambdaAntiS_vz"]->Fill(GENAntiLambdaAntiS->vz());	
	histos_th1f["h_GEN_AntiLambdaAntiS_lxyz"]->Fill(Lxyz);	
	histos_th1f["h_GEN_AntiLambdaAntiS_dxy"]->Fill(dxy);
	histos_th1f["h_GEN_AntiLambdaAntiS_dz"]->Fill(dz);
	histos_th2f["h2_GEN_AntiLambdaAntiS_dxy_dz"]->Fill(dxy,dz);
	if(GENAntiLambdaAntiS->numberOfDaughters() == 2)histos_th1f["h_GEN_AntiLambdaAntiS_XYpointingAngle"]->Fill(xypointingAngle);
	
}
*/



/*
void FlatTreeProducerTracking::RecoEvaluationKsNonAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, TVector3 beamspot){
  const reco::Candidate  * GENKsNonAntiS = genParticle; 
  if(h_V0Ks.isValid()){
      double deltaRmin = 999.;
      const reco::VertexCompositeCandidate * RECOKs = nullptr;
      for(unsigned int i = 0; i < h_V0Ks->size(); ++i){//loop all Ks candidates
	const reco::VertexCompositeCandidate * Ks = &h_V0Ks->at(i);	
	double deltaPhi = reco::deltaPhi(Ks->phi(),genParticle->phi());
	double deltaEta = Ks->eta() - genParticle->eta();
	double deltaR = pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
        if(deltaR < deltaRmin){ deltaRmin = deltaR; RECOKs = Ks;}
       }	
	histos_th1f["h_GENRECO_KsNonAntiS_deltaRmin"]->Fill(deltaRmin);
	TVector3 KsCreationVertex(GENKsNonAntiS->vx(),GENKsNonAntiS->vy(),GENKsNonAntiS->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,KsCreationVertex);
	TVector3 KsMomentum(GENKsNonAntiS->px(),GENKsNonAntiS->py(),GENKsNonAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,beamspot);

	double RECOKsLxy = -999;
	double RECOKsdxy = -999;
	double RECOKsdz  = -999;
	if(RECOKs){
		TVector3 RECOKsCreationVertex(RECOKs->vx(),RECOKs->vy(),RECOKs->vz());
		RECOKsLxy = AnalyzerAllSteps::lxy(beamspot,RECOKsCreationVertex);
		TVector3 RECOKsMomentum(RECOKs->px(),RECOKs->py(),RECOKs->pz());
		RECOKsdxy = AnalyzerAllSteps::dxy_signed_line_point(RECOKsCreationVertex,RECOKsMomentum,beamspot);
		RECOKsdz = AnalyzerAllSteps::dz_line_point(RECOKsCreationVertex,RECOKsMomentum,beamspot);
	}

	if(deltaRmin<AnalyzerAllSteps::deltaRCutV0RECOKs){//matched
		
		histos_th1f["h_GENRECO_RECO_KsNonAntiS_pt"]->Fill(GENKsNonAntiS->pt());	
		histos_th1f["h_GENRECO_RECO_KsNonAntiS_eta"]->Fill(GENKsNonAntiS->eta());	
		histos_th1f["h_GENRECO_RECO_KsNonAntiS_phi"]->Fill(GENKsNonAntiS->phi());	
		histos_th2f["h2_GENRECO_RECO_KsNonAntiS_vx_vy"]->Fill(GENKsNonAntiS->vx(),GENKsNonAntiS->vy());
		histos_th2f["h2_GENRECO_RECO_KsNonAntiS_lxy_vz"]->Fill(Lxy,GENKsNonAntiS->vz());
		histos_th1f["h_GENRECO_RECO_KsNonAntiS_lxy"]->Fill(Lxy);	
		histos_th1f["h_GENRECO_RECO_KsNonAntiS_vz"]->Fill(GENKsNonAntiS->vz());	
		histos_th1f["h_GENRECO_RECO_KsNonAntiS_dxy"]->Fill(dxy);

		if(RECOKs->mass() > 0.48 && RECOKs->mass() < 0.52 && RECOKs){
			histos_th2f["h2_GENRECO_RECO_RECO_Ks_lxy_vz"]->Fill(abs(RECOKs->vz()),RECOKsLxy);
                        histos_th2f["h2_GENRECO_RECO_RECO_Ks_lxy_dz"]->Fill(abs(RECOKsdz),RECOKsLxy);
                        histos_th1f["h_GENRECO_RECO_RECO_Ks_dxy"]->Fill(RECOKsdxy);
                        histos_th1f["h_GENRECO_RECO_RECO_Ks_dz"]->Fill(RECOKsdz);
                        histos_th2f["h2_GENRECO_RECO_RECO_Ks_mass_dz"]->Fill(abs(RECOKsdz),RECOKs->mass());
		}

	}


	histos_th1f["h_GENRECO_All_KsNonAntiS_pt"]->Fill(GENKsNonAntiS->pt());	
	histos_th1f["h_GENRECO_All_KsNonAntiS_eta"]->Fill(GENKsNonAntiS->eta());	
	histos_th1f["h_GENRECO_All_KsNonAntiS_phi"]->Fill(GENKsNonAntiS->phi());	
	histos_th2f["h2_GENRECO_All_KsNonAntiS_vx_vy"]->Fill(GENKsNonAntiS->vx(),GENKsNonAntiS->vy());
	histos_th2f["h2_GENRECO_All_KsNonAntiS_lxy_vz"]->Fill(Lxy,GENKsNonAntiS->vz());
	histos_th1f["h_GENRECO_All_KsNonAntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_GENRECO_All_KsNonAntiS_vz"]->Fill(GENKsNonAntiS->vz());	
	histos_th1f["h_GENRECO_All_KsNonAntiS_dxy"]->Fill(dxy);	


      }	
  
}
*/

/*
void FlatTreeProducerTracking::RecoEvaluationKsAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, TVector3 beamspot){
  const reco::Candidate  * GENKsAntiS = genParticle; 
  if(h_V0Ks.isValid()){
      double deltaRmin = 999.;
      for(unsigned int i = 0; i < h_V0Ks->size(); ++i){//loop all Ks candidates
	const reco::VertexCompositeCandidate * Ks = &h_V0Ks->at(i);	
	double deltaPhi = reco::deltaPhi(Ks->phi(),genParticle->phi());
	double deltaEta = Ks->eta() - genParticle->eta();
	double deltaR = pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
        if(deltaR < deltaRmin) deltaRmin = deltaR;
      }
      histos_th1f["h_GENRECO_KsAntiS_deltaRmin"]->Fill(deltaRmin);

	TVector3 KsCreationVertex(GENKsAntiS->vx(),GENKsAntiS->vy(),GENKsAntiS->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,KsCreationVertex);
	TVector3 KsMomentum(GENKsAntiS->px(),GENKsAntiS->py(),GENKsAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,beamspot);


	if(deltaRmin<AnalyzerAllSteps::deltaRCutV0RECOKs){//matched
		
		histos_th1f["h_GENRECO_RECO_KsAntiS_pt"]->Fill(GENKsAntiS->pt());	
		histos_th1f["h_GENRECO_RECO_KsAntiS_eta"]->Fill(GENKsAntiS->eta());	
		histos_th1f["h_GENRECO_RECO_KsAntiS_phi"]->Fill(GENKsAntiS->phi());	
		histos_th2f["h2_GENRECO_RECO_KsAntiS_vx_vy"]->Fill(GENKsAntiS->vx(),GENKsAntiS->vy());
		histos_th2f["h2_GENRECO_RECO_KsAntiS_lxy_vz"]->Fill(Lxy,GENKsAntiS->vz());
		histos_th1f["h_GENRECO_RECO_KsAntiS_lxy"]->Fill(Lxy);	
		histos_th1f["h_GENRECO_RECO_KsAntiS_vz"]->Fill(GENKsAntiS->vz());	
		histos_th1f["h_GENRECO_RECO_KsAntiS_dxy"]->Fill(dxy);	

	}


	histos_th1f["h_GENRECO_All_KsAntiS_pt"]->Fill(GENKsAntiS->pt());	
	histos_th1f["h_GENRECO_All_KsAntiS_eta"]->Fill(GENKsAntiS->eta());	
	histos_th1f["h_GENRECO_All_KsAntiS_phi"]->Fill(GENKsAntiS->phi());	
	histos_th2f["h2_GENRECO_All_KsAntiS_vx_vy"]->Fill(GENKsAntiS->vx(),GENKsAntiS->vy());
	histos_th2f["h2_GENRECO_All_KsAntiS_lxy_vz"]->Fill(Lxy,GENKsAntiS->vz());
	histos_th1f["h_GENRECO_All_KsAntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_GENRECO_All_KsAntiS_vz"]->Fill(GENKsAntiS->vz());	
	histos_th1f["h_GENRECO_All_KsAntiS_dxy"]->Fill(dxy);	


      	
  }
}
*/

/*
void FlatTreeProducerTracking::RecoEvaluationAntiLambdaNonAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L, TVector3 beamspot){
  const reco::Candidate  * GENAntiLambdaNonAntiS = genParticle; 
  if(h_V0L.isValid()){
      double deltaRmin = 999.;
      const reco::VertexCompositeCandidate * RECOLambda = nullptr;
      for(unsigned int i = 0; i < h_V0L->size(); ++i){//loop all AntiLambda candidates
	const reco::VertexCompositeCandidate * AntiLambda = &h_V0L->at(i);	
	double deltaPhi = reco::deltaPhi(AntiLambda->phi(),genParticle->phi());
	double deltaEta = AntiLambda->eta() - genParticle->eta();
	double deltaR = pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
	if(deltaR < deltaRmin){ deltaRmin = deltaR;RECOLambda = AntiLambda;}
      }	
	histos_th1f["h_GENRECO_AntiLambdaNonAntiS_deltaRmin"]->Fill(deltaRmin);

	TVector3 AntiLambdaCreationVertex(GENAntiLambdaNonAntiS->vx(),GENAntiLambdaNonAntiS->vy(),GENAntiLambdaNonAntiS->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,AntiLambdaCreationVertex);
	TVector3 AntiLambdaMomentum(GENAntiLambdaNonAntiS->px(),GENAntiLambdaNonAntiS->py(),GENAntiLambdaNonAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(AntiLambdaCreationVertex,AntiLambdaMomentum,beamspot);

	double RECOLambdaLxy = -999;
	double RECOLambdadxy = -999;
	double RECOLambdadz = -999;
	if(RECOLambda){
		TVector3 RECOLambdaCreationVertex(RECOLambda->vx(),RECOLambda->vy(),RECOLambda->vz());
		RECOLambdaLxy = AnalyzerAllSteps::lxy(beamspot,RECOLambdaCreationVertex);
		TVector3 RECOLambdaMomentum(RECOLambda->px(),RECOLambda->py(),RECOLambda->pz());
		RECOLambdadxy = AnalyzerAllSteps::dxy_signed_line_point(RECOLambdaCreationVertex,RECOLambdaMomentum,beamspot);
		RECOLambdadz = AnalyzerAllSteps::dz_line_point(RECOLambdaCreationVertex,RECOLambdaMomentum,beamspot);
	}

	if(deltaRmin<AnalyzerAllSteps::deltaRCutV0RECOLambda){//matched
		
		histos_th1f["h_GENRECO_RECO_AntiLambdaNonAntiS_pt"]->Fill(GENAntiLambdaNonAntiS->pt());	
		histos_th1f["h_GENRECO_RECO_AntiLambdaNonAntiS_eta"]->Fill(GENAntiLambdaNonAntiS->eta());	
		histos_th1f["h_GENRECO_RECO_AntiLambdaNonAntiS_phi"]->Fill(GENAntiLambdaNonAntiS->phi());	
		histos_th2f["h2_GENRECO_RECO_AntiLambdaNonAntiS_vx_vy"]->Fill(GENAntiLambdaNonAntiS->vx(),GENAntiLambdaNonAntiS->vy());
		histos_th2f["h2_GENRECO_RECO_AntiLambdaNonAntiS_lxy_vz"]->Fill(Lxy,GENAntiLambdaNonAntiS->vz());
		histos_th1f["h_GENRECO_RECO_AntiLambdaNonAntiS_lxy"]->Fill(Lxy);	
		histos_th1f["h_GENRECO_RECO_AntiLambdaNonAntiS_vz"]->Fill(GENAntiLambdaNonAntiS->vz());	
		histos_th1f["h_GENRECO_RECO_AntiLambdaNonAntiS_dxy"]->Fill(dxy);	

		if(RECOLambda->mass() > 1.105 && RECOLambda->mass() < 1.125 && RECOLambda){
                        histos_th2f["h2_GENRECO_RECO_RECO_Lambda_lxy_vz"]->Fill(abs(RECOLambda->vz()),RECOLambdaLxy);
                        histos_th2f["h2_GENRECO_RECO_RECO_Lambda_lxy_dz"]->Fill(abs(RECOLambdadz),RECOLambdaLxy);
                        histos_th1f["h_GENRECO_RECO_RECO_Lambda_dxy"]->Fill(RECOLambdadxy);
                        histos_th1f["h_GENRECO_RECO_RECO_Lambda_dz"]->Fill(RECOLambdadz);
                        histos_th2f["h2_GENRECO_RECO_RECO_Lambda_mass_dz"]->Fill(abs(RECOLambdadz),RECOLambda->mass());
                }

	}


	histos_th1f["h_GENRECO_All_AntiLambdaNonAntiS_pt"]->Fill(GENAntiLambdaNonAntiS->pt());	
	histos_th1f["h_GENRECO_All_AntiLambdaNonAntiS_eta"]->Fill(GENAntiLambdaNonAntiS->eta());	
	histos_th1f["h_GENRECO_All_AntiLambdaNonAntiS_phi"]->Fill(GENAntiLambdaNonAntiS->phi());	
	histos_th2f["h2_GENRECO_All_AntiLambdaNonAntiS_vx_vy"]->Fill(GENAntiLambdaNonAntiS->vx(),GENAntiLambdaNonAntiS->vy());
	histos_th2f["h2_GENRECO_All_AntiLambdaNonAntiS_lxy_vz"]->Fill(Lxy,GENAntiLambdaNonAntiS->vz());
	histos_th1f["h_GENRECO_All_AntiLambdaNonAntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_GENRECO_All_AntiLambdaNonAntiS_vz"]->Fill(GENAntiLambdaNonAntiS->vz());	
	histos_th1f["h_GENRECO_All_AntiLambdaNonAntiS_dxy"]->Fill(dxy);	


      	
  }
}
*/

/*
void FlatTreeProducerTracking::RecoEvaluationAntiLambdaAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L, TVector3 beamspot){
  const reco::Candidate  * GENAntiLambdaAntiS = genParticle; 
  if(h_V0L.isValid()){
	double deltaRmin = 999.;
      for(unsigned int i = 0; i < h_V0L->size(); ++i){//loop all AntiLambda candidates
	const reco::VertexCompositeCandidate * AntiLambda = &h_V0L->at(i);	
	double deltaPhi = reco::deltaPhi(AntiLambda->phi(),genParticle->phi());
	double deltaEta = AntiLambda->eta() - genParticle->eta();
	double deltaR = pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
	if(deltaR < deltaRmin) deltaRmin = deltaR;
      }
	histos_th1f["h_GENRECO_AntiLambdaAntiS_deltaRmin"]->Fill(deltaRmin);

	
	TVector3 AntiLambdaCreationVertex(GENAntiLambdaAntiS->vx(),GENAntiLambdaAntiS->vy(),GENAntiLambdaAntiS->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,AntiLambdaCreationVertex);
	TVector3 AntiLambdaMomentum(GENAntiLambdaAntiS->px(),GENAntiLambdaAntiS->py(),GENAntiLambdaAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(AntiLambdaCreationVertex,AntiLambdaMomentum,beamspot);

	if(deltaRmin<AnalyzerAllSteps::deltaRCutV0RECOLambda){//matched
		
		histos_th1f["h_GENRECO_RECO_AntiLambdaAntiS_pt"]->Fill(GENAntiLambdaAntiS->pt());	
		histos_th1f["h_GENRECO_RECO_AntiLambdaAntiS_eta"]->Fill(GENAntiLambdaAntiS->eta());	
		histos_th1f["h_GENRECO_RECO_AntiLambdaAntiS_phi"]->Fill(GENAntiLambdaAntiS->phi());	
		histos_th2f["h2_GENRECO_RECO_AntiLambdaAntiS_vx_vy"]->Fill(GENAntiLambdaAntiS->vx(),GENAntiLambdaAntiS->vy());
		histos_th2f["h2_GENRECO_RECO_AntiLambdaAntiS_lxy_vz"]->Fill(Lxy,GENAntiLambdaAntiS->vz());
		histos_th1f["h_GENRECO_RECO_AntiLambdaAntiS_lxy"]->Fill(Lxy);	
		histos_th1f["h_GENRECO_RECO_AntiLambdaAntiS_vz"]->Fill(GENAntiLambdaAntiS->vz());	
		histos_th1f["h_GENRECO_RECO_AntiLambdaAntiS_dxy"]->Fill(dxy);	

	}


	histos_th1f["h_GENRECO_All_AntiLambdaAntiS_pt"]->Fill(GENAntiLambdaAntiS->pt());	
	histos_th1f["h_GENRECO_All_AntiLambdaAntiS_eta"]->Fill(GENAntiLambdaAntiS->eta());	
	histos_th1f["h_GENRECO_All_AntiLambdaAntiS_phi"]->Fill(GENAntiLambdaAntiS->phi());	
	histos_th2f["h2_GENRECO_All_AntiLambdaAntiS_vx_vy"]->Fill(GENAntiLambdaAntiS->vx(),GENAntiLambdaAntiS->vy());
	histos_th2f["h2_GENRECO_All_AntiLambdaAntiS_lxy_vz"]->Fill(Lxy,GENAntiLambdaAntiS->vz());
	histos_th1f["h_GENRECO_All_AntiLambdaAntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_GENRECO_All_AntiLambdaAntiS_vz"]->Fill(GENAntiLambdaAntiS->vz());	
	histos_th1f["h_GENRECO_All_AntiLambdaAntiS_dxy"]->Fill(dxy);	


  }
}
*/

/*
void FlatTreeProducerTracking::RecoEvaluationAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV){
  const reco::Candidate  * GENAntiS = genParticle;   
  
  if(h_sCands.isValid()){
        double deltaRmin = 999.;
        const reco::VertexCompositeCandidate * bestRECOAntiS = nullptr;
        for(unsigned int i = 0; i < h_sCands->size(); ++i){//loop all AntiS RECO candidates
		const reco::VertexCompositeCandidate * AntiS = &h_sCands->at(i);	
		double deltaPhi = reco::deltaPhi(AntiS->phi(),genParticle->phi());
		double deltaEta = AntiS->eta() - genParticle->eta();
		double deltaR = pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
		if(deltaR < deltaRmin){ deltaRmin = deltaR; bestRECOAntiS = AntiS;}
        }
	histos_th1f["h_GENRECO_AntiS_deltaRmin"]->Fill(deltaRmin);

	

	//calculate some kinmatic variables for the GEN antiS	
	TVector3 GENAntiSCreationVertex(GENAntiS->vx(),GENAntiS->vy(),GENAntiS->vz());
	double GENLxy = AnalyzerAllSteps::lxy(beamspot,GENAntiSCreationVertex);
	TVector3 GENAntiSMomentum(GENAntiS->px(),GENAntiS->py(),GENAntiS->pz());
	double GENdxy = AnalyzerAllSteps::dxy_signed_line_point(GENAntiSCreationVertex,GENAntiSMomentum,beamspot);
	double GENdz = AnalyzerAllSteps::dz_line_point(GENAntiSCreationVertex,GENAntiSMomentum,beamspot);
	TVector3 GENAntiSInteractionVertex(GENAntiS->daughter(0)->vx(),GENAntiS->daughter(0)->vy(),GENAntiS->daughter(0)->vz());
	double GENLxy_interactionVertex = AnalyzerAllSteps::lxy(beamspot,GENAntiSInteractionVertex);
	
	//calculate some kinematic variables for the RECO AntiS
	TVector3 RECOAntiSInteractionVertex(bestRECOAntiS->vx(),bestRECOAntiS->vy(),bestRECOAntiS->vz());//this is the interaction vertex of the antiS and the neutron. Check in the skimming code if you want to check.
	TVector3 RECOAntiSMomentumVertex(bestRECOAntiS->px(),bestRECOAntiS->py(),bestRECOAntiS->pz());
	double RECOLxy_interactionVertex = AnalyzerAllSteps::lxy(beamspot,RECOAntiSInteractionVertex);
	double RECOErrorLxy_interactionVertex = AnalyzerAllSteps::std_dev_lxy(bestRECOAntiS->vx(), bestRECOAntiS->vy(), bestRECOAntiS->vertexCovariance(0,0), bestRECOAntiS->vertexCovariance(1,1), beamspot.X(), beamspot.Y(), beamspotVariance.X(), beamspotVariance.Y());
	double RECODeltaPhiDaughters = reco::deltaPhi(bestRECOAntiS->daughter(0)->phi(),bestRECOAntiS->daughter(1)->phi());
	double RECODeltaEtaDaughters = bestRECOAntiS->daughter(0)->eta()-bestRECOAntiS->daughter(1)->eta();
	double RECODeltaRDaughters = pow(RECODeltaPhiDaughters*RECODeltaPhiDaughters+RECODeltaEtaDaughters*RECODeltaEtaDaughters,0.5);

	reco::LeafCandidate::LorentzVector n_(0,0,0,0.939565);
	double bestRECOAntiSmass = (bestRECOAntiS->p4()-n_).mass();

	//the dxy of the Ks and Lambda
	TVector3 RECOAntiSDaug0Momentum(bestRECOAntiS->daughter(0)->px(),bestRECOAntiS->daughter(0)->py(),bestRECOAntiS->daughter(0)->pz());
        TVector3 RECOAntiSDaug1Momentum(bestRECOAntiS->daughter(1)->px(),bestRECOAntiS->daughter(1)->py(),bestRECOAntiS->daughter(1)->pz());
	reco::Candidate::Vector vRECOAntiSDaug0Momentum(bestRECOAntiS->daughter(0)->px(),bestRECOAntiS->daughter(0)->py(),bestRECOAntiS->daughter(0)->pz());
        reco::Candidate::Vector vRECOAntiSDaug1Momentum(bestRECOAntiS->daughter(1)->px(),bestRECOAntiS->daughter(1)->py(),bestRECOAntiS->daughter(1)->pz());
	double RECOOpeningsAngleDaughters = AnalyzerAllSteps::openings_angle(vRECOAntiSDaug0Momentum,vRECOAntiSDaug1Momentum);
        double RECO_dxy_daughter0 = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot);
        double RECO_dxy_daughter1 = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot);
	//the dz of the Ks and Lambda
	double RECO_dz_daughter0 = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot);
	double RECO_dz_daughter1 = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot);
	//the dxyz of the Ks and Lambda
	double RECO_dxyz_daughter0 = AnalyzerAllSteps::dxyz_signed_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot);
	double RECO_dxyz_daughter1 = AnalyzerAllSteps::dxyz_signed_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot);
	//collapse the RECO_dxy, RECO_dz, RECO_dxyz variables on one variable
	double relDiff_RECO_dxy_daughters = (abs(RECO_dxy_daughter0)-abs(RECO_dxy_daughter1))/(abs(RECO_dxy_daughter0)+abs(RECO_dxy_daughter1));
	double relDiff_RECO_dz_daughters = (abs(RECO_dz_daughter0)-abs(RECO_dz_daughter1))/(abs(RECO_dz_daughter0)+abs(RECO_dz_daughter1));
	double relDiff_RECO_dxyz_daughters = (abs(RECO_dxyz_daughter0)-abs(RECO_dxyz_daughter1))/(abs(RECO_dxyz_daughter0)+abs(RECO_dxyz_daughter1));
	//dxy and dz of the AntiS itself
	double RECO_dxy_antiS = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,beamspot);
	double RECO_dz_antiS = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,beamspot);
	

	//calculate the sign of the dot product between the displacement vector and the dxy vector for both the Ks and the Lambda
	double signLxyDotdxy_daughter0 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot)*RECOAntiSInteractionVertex);
	double signLxyDotdxy_daughter1 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot)*RECOAntiSInteractionVertex);
	//calculate the sign of the dot product between the displacement vector and the pt of both the Ks and the Lambda
	double signPtDotdxy_daughter0 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot)*RECOAntiSDaug0Momentum);
	double signPtDotdxy_daughter1 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot)*RECOAntiSDaug1Momentum);

	//loop over all PVs and find the one which minimises the dxy of the antiS
	double RECOdzAntiSPVmin = 999.;
	TVector3  bestPVdzAntiS;
	for(unsigned int i = 0; i < h_offlinePV->size(); ++i){
		TVector3 PV(h_offlinePV->at(i).x(),h_offlinePV->at(i).y(),h_offlinePV->at(i).z());
		double dzAntiSPV = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,PV);
		if(abs(dzAntiSPV) < abs(RECOdzAntiSPVmin)) {RECOdzAntiSPVmin = dzAntiSPV; bestPVdzAntiS = PV;}
	}
	
	double dxyAntiSPVmin = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,bestPVdzAntiS);
	
	if(deltaRmin<AnalyzerAllSteps::deltaRCutRECOAntiS){//matched
		//kinematics of the matched GEN particles	
		histos_th1f["h_GENRECO_RECO_AntiS_m"]->Fill(bestRECOAntiSmass);	
		histos_th1f["h_GENRECO_RECO_AntiS_pt"]->Fill(bestRECOAntiS->pt());	
		histos_th1f["h_GENRECO_RECO_AntiS_pz"]->Fill(bestRECOAntiS->pz());	
		histos_th1f["h_GENRECO_RECO_AntiS_eta"]->Fill(bestRECOAntiS->eta());	
		histos_th1f["h_GENRECO_RECO_AntiS_phi"]->Fill(bestRECOAntiS->phi());	
		histos_th2f["h2_GENRECO_RECO_AntiS_vx_vy"]->Fill(bestRECOAntiS->vx(),bestRECOAntiS->vy());
		histos_th2f["h2_GENRECO_RECO_AntiS_lxy_vz"]->Fill(RECOLxy_interactionVertex,bestRECOAntiS->vz());
		histos_th1f["h_GENRECO_RECO_AntiS_lxy_interactionVertex"]->Fill(RECOLxy_interactionVertex);	
		histos_th1f["h_GENRECO_RECO_AntiS_vz"]->Fill(bestRECOAntiS->vz());	
		histos_th1f["h_GENRECO_RECO_AntiS_vz_interactionVertex"]->Fill(RECOAntiSInteractionVertex.Z());	
		histos_th1f["h_GENRECO_RECO_AntiS_dxy"]->Fill(RECO_dxy_antiS);	
		histos_th1f["h_GENRECO_RECO_AntiS_dz"]->Fill(RECO_dz_antiS);
		histos_th1f["h_GENRECO_RECO_AntiS_deltaRDaughters"]->Fill(RECODeltaRDaughters);	
		histos_th1f["h_GENRECO_RECO_AntiS_openingsAngleDaughters"]->Fill(RECOOpeningsAngleDaughters);	
		histos_th1f["h_GENRECO_RECO_AntiS_pt_Ks"]->Fill(bestRECOAntiS->daughter(1)->pt());	
		histos_th1f["h_GENRECO_RECO_AntiS_pt_AntiL"]->Fill(bestRECOAntiS->daughter(0)->pt());	
		
		
		//so now if you found a RECO antiS compare it to the GEN antiS. Compare parameters such as GEN vs RECO pt, GEN vs RECO interaction vertex, ...
		if(RECOLxy_interactionVertex > AnalyzerAllSteps::MinLxyCut && RECOErrorLxy_interactionVertex < AnalyzerAllSteps::MaxErrorLxyCut && bestRECOAntiSmass > 0){
			histos_th1f["h_RECOAcc_AntiS_pt"]->Fill(GENAntiS->pt()-bestRECOAntiS->pt());
			histos_th1f["h_RECOAcc_AntiS_phi"]->Fill(reco::deltaPhi(GENAntiS->phi(),bestRECOAntiS->phi()));
			histos_th1f["h_RECOAcc_AntiS_eta"]->Fill(GENAntiS->eta()-bestRECOAntiS->eta());
			histos_th1f["h_RECOAcc_AntiS_mass"]->Fill(GENAntiS->mass()-bestRECOAntiSmass);
			histos_th1f["h_RECOAcc_AntiS_lxy_interactionVertex"]->Fill(GENLxy_interactionVertex-RECOLxy_interactionVertex);
			histos_th1f["h_RECOAcc_AntiS_vz_interactionVertex"]->Fill(GENAntiSInteractionVertex.Z()-RECOAntiSInteractionVertex.Z());
		}	
	
		//so now you know that these RECO AntiS are really signal, so what you can do is plot the variables which you will cut on to get rid of the background for the signal
		histos_th1f["h1_RECOSignal_AntiS_corr_dxy_Ks"]->Fill(RECO_dxy_daughter1);
		histos_th1f["h1_RECOSignal_AntiS_corr_dxy_over_lxy_Ks"]->Fill(RECO_dxy_daughter1/RECOLxy_interactionVertex);
		histos_th1f["h1_RECOSignal_AntiS_corr_dxy_AntiL"]->Fill(RECO_dxy_daughter0);
		histos_th1f["h1_RECOSignal_AntiS_corr_dxy_over_lxy_AntiL"]->Fill(RECO_dxy_daughter0/RECOLxy_interactionVertex);
		histos_th1f["h1_RECOSignal_AntiS_corr_dz_Ks"]->Fill(RECO_dz_daughter1);
		histos_th1f["h1_RECOSignal_AntiS_corr_dz_AntiL"]->Fill(RECO_dz_daughter0);


		histos_th2f["h2_RECOSignal_AntiS_dxy_Ks_vs_lxy_Ks"]->Fill(RECOLxy_interactionVertex,RECO_dxy_daughter1);
		histos_th2f["h2_RECOSignal_AntiS_dxy_AntiL_vs_lxy_AntiL"]->Fill(RECOLxy_interactionVertex,RECO_dxy_daughter0);
		

		histos_th2f["h2_RECOSignal_AntiS_corr_dxy_daughters"]->Fill(RECO_dxy_daughter0, RECO_dxy_daughter1);
		histos_th2f["h2_RECOSignal_AntiS_corr_dxy_sign_dotProdLxydxy_daughters"]->Fill(abs(RECO_dxy_daughter0)*signLxyDotdxy_daughter0,abs(RECO_dxy_daughter1)*signLxyDotdxy_daughter1);
		histos_th2f["h2_RECOSignal_AntiS_corr_dxy_sign_dotProdptdxy_daughters"]->Fill(abs(RECO_dxy_daughter0)*signPtDotdxy_daughter0,abs(RECO_dxy_daughter1)*signPtDotdxy_daughter1);
		histos_th2f["h2_RECOSignal_AntiS_corr_dxy_over_lxy_daughters"]->Fill(RECO_dxy_daughter0/RECOLxy_interactionVertex, RECO_dxy_daughter1/RECOLxy_interactionVertex);
		histos_th2f["h2_RECOSignal_AntiS_corr_dz_daughters"]->Fill(RECO_dz_daughter0, RECO_dz_daughter1);
		histos_th2f["h2_RECOSignal_AntiS_corr_dxyz_daughters"]->Fill(RECO_dxyz_daughter0, RECO_dxyz_daughter1);
		histos_th1f["h_RECOSignal_AntiS_relDiff_RECO_dxy_daughters"]->Fill(relDiff_RECO_dxy_daughters);
		histos_th1f["h_RECOSignal_AntiS_relDiff_RECO_dz_daughters"]->Fill(relDiff_RECO_dz_daughters);
		histos_th1f["h_RECOSignal_AntiS_relDiff_RECO_dxyz_daughters"]->Fill(relDiff_RECO_dxyz_daughters);
		histos_th2f["h2_RECOSignal_AntiS_lxy_error_lxy"]->Fill(RECOLxy_interactionVertex,  RECOErrorLxy_interactionVertex);
		histos_th1f["h_RECOSignal_AntiS_dxy"]->Fill(RECO_dxy_antiS);	
		histos_th1f["h_RECOSignal_AntiS_dxy_over_lxy"]->Fill(RECO_dxy_antiS/RECOLxy_interactionVertex);	
		histos_th2f["h2_RECOSignal_AntiS_dxy_vs_dxy_over_lxy"]->Fill(RECO_dxy_antiS/RECOLxy_interactionVertex, RECO_dxy_antiS);
		histos_th1f["h_RECOSignal_AntiS_dz"]->Fill(RECO_dz_antiS);	
		histos_th1f["h_RECOSignal_AntiS_vertexNormalizedChi2"]->Fill(bestRECOAntiS->vertexNormalizedChi2());
			
		histos_th2f["h2_RECOSignal_AntiS_deltaPhi_deltaEta"]->Fill(RECODeltaPhiDaughters, RECODeltaEtaDaughters);
		histos_th2f["h2_RECOSignal_AntiS_openings_angle_vs_deltaEta"]->Fill(RECODeltaEtaDaughters,RECOOpeningsAngleDaughters);

		histos_th1f["h_RECOSignal_AntiS_dxyAntiSPVmin"]->Fill(dxyAntiSPVmin);
		histos_th1f["h_RECOSignal_AntiS_RECOdzAntiSPVmin"]->Fill(RECOdzAntiSPVmin);
	
		//plot the important background cuts:
		histos_th1f["h_RECO_AntiS_lxy_LCuts"]->Fill(RECOLxy_interactionVertex);
		histos_th1f["h_RECO_AntiS_error_lxy_LCuts"]->Fill(RECOErrorLxy_interactionVertex);
		histos_th1f["h_RECO_AntiS_mass_LCuts"]->Fill(bestRECOAntiSmass);
		histos_th1f["h_RECO_AntiS_vertex_chi2_ndof_LCuts"]->Fill(bestRECOAntiS->vertexNormalizedChi2());
		histos_th1f["h_RECO_AntiS_deltaPhiDaughters_LCuts"]->Fill(RECODeltaPhiDaughters);
		histos_th1f["h_RECO_AntiS_openingsAngleDaughters_LCuts"]->Fill(RECOOpeningsAngleDaughters);
		histos_th1f["h_RECO_AntiS_deltaEtaDaughters_LCuts"]->Fill(RECODeltaEtaDaughters);
		histos_th1f["h_RECO_AntiS_dxy_over_lxy_LCuts"]->Fill(RECO_dxy_antiS/RECOLxy_interactionVertex);
		histos_th1f["h_RECO_AntiS_dxy_over_lxy_Ks_LCuts"]->Fill(RECO_dxy_daughter0/RECOLxy_interactionVertex);
		histos_th1f["h_RECO_AntiS_dxy_over_lxy_AntiL_LCuts"]->Fill(RECO_dxy_daughter1/RECOLxy_interactionVertex);
		histos_th1f["h_RECO_AntiS_pt_AntiL_LCuts"]->Fill(bestRECOAntiS->daughter(0)->pt());
		histos_th1f["h_RECO_AntiS_dzPVmin_LCuts"]->Fill(RECOdzAntiSPVmin);
		histos_th1f["h_RECO_AntiS_vz_LCuts"]->Fill(abs(bestRECOAntiS->vz()));
		histos_th1f["h_RECO_AntiS_eta_LCuts"]->Fill(bestRECOAntiS->eta());


		//for these RECO antiS which for sure are signal plot how much yo will cut on the signal using the bkg cuts
		histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(0);
		if(RECOLxy_interactionVertex>AnalyzerAllSteps::MinLxyCut)				histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(1);
		if(RECOErrorLxy_interactionVertex < AnalyzerAllSteps::MaxErrorLxyCut) 			histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(2);
		if(bestRECOAntiSmass > AnalyzerAllSteps::MinErrorMassCut)				histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(3);
		if(bestRECOAntiS->vertexNormalizedChi2() < AnalyzerAllSteps::MaxNormChi2Cut)		histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(4);

		if(abs(RECODeltaPhiDaughters) > AnalyzerAllSteps::MinAbsDeltaPhiDaughtersCut)		histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(5);
		if(abs(RECODeltaPhiDaughters) < AnalyzerAllSteps::MaxAbsDeltaPhiDaughtersCut)		histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(6);
	   	if(RECOOpeningsAngleDaughters<AnalyzerAllSteps::MaxOpeningsAngleCut)			histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(7);		
	   	if(abs(RECODeltaEtaDaughters)<AnalyzerAllSteps::MaxAbsDeltaEtaDaughCut)			histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(8);
		if(RECO_dxy_antiS/RECOLxy_interactionVertex > AnalyzerAllSteps::MinDxyOverLxyCut)	histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(9);
		if(RECO_dxy_antiS/RECOLxy_interactionVertex < AnalyzerAllSteps::MaxDxyOverLxyCut)	histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(10);

		if(RECO_dxy_daughter0/RECOLxy_interactionVertex>AnalyzerAllSteps::DxyKsExclusionRangeMaxCut || RECO_dxy_daughter0/RECOLxy_interactionVertex<AnalyzerAllSteps::DxyKsExclusionRangeMinCut)histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(11);
		if(RECO_dxy_daughter1/RECOLxy_interactionVertex>AnalyzerAllSteps::DxyAntiLExclusionRangeMaxCut || RECO_dxy_daughter1/RECOLxy_interactionVertex<AnalyzerAllSteps::DxyAntiLExclusionRangeMinCut) histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(12);

		if(bestRECOAntiS->daughter(0)->pt() > AnalyzerAllSteps::MinLambdaPtCut)			histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(13);
		if(abs(RECOdzAntiSPVmin)<AnalyzerAllSteps::dzAntiSPVminCut)				histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(14);
		if(abs(bestRECOAntiS->vz())>AnalyzerAllSteps::vzAntiSInteractionVertexCut)		histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(15);
		if(abs(bestRECOAntiS->eta())>AnalyzerAllSteps::antiSEtaCut)				histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(16);

		if(	RECOLxy_interactionVertex>AnalyzerAllSteps::MinLxyCut && 
			RECOErrorLxy_interactionVertex < AnalyzerAllSteps::MaxErrorLxyCut &&  
			bestRECOAntiSmass > AnalyzerAllSteps::MinErrorMassCut && 
			bestRECOAntiS->vertexNormalizedChi2() < AnalyzerAllSteps::MaxNormChi2Cut && 

			abs(RECODeltaPhiDaughters) > AnalyzerAllSteps::MinAbsDeltaPhiDaughtersCut && 
			abs(RECODeltaPhiDaughters) < AnalyzerAllSteps::MaxAbsDeltaPhiDaughtersCut && 
			RECOOpeningsAngleDaughters < AnalyzerAllSteps::MaxOpeningsAngleCut && 
			abs(RECODeltaEtaDaughters) < AnalyzerAllSteps::MaxAbsDeltaEtaDaughCut && 
			RECO_dxy_antiS/RECOLxy_interactionVertex > AnalyzerAllSteps::MinDxyOverLxyCut && 
			RECO_dxy_antiS/RECOLxy_interactionVertex < AnalyzerAllSteps::MaxDxyOverLxyCut && 
		
			(RECO_dxy_daughter0/RECOLxy_interactionVertex > AnalyzerAllSteps::DxyKsExclusionRangeMaxCut || RECO_dxy_daughter0/RECOLxy_interactionVertex < AnalyzerAllSteps::DxyKsExclusionRangeMinCut) && 
			(RECO_dxy_daughter1/RECOLxy_interactionVertex > AnalyzerAllSteps::DxyAntiLExclusionRangeMaxCut || RECO_dxy_daughter1/RECOLxy_interactionVertex<AnalyzerAllSteps::DxyAntiLExclusionRangeMinCut) &&
		
			bestRECOAntiS->daughter(0)->pt() > AnalyzerAllSteps::MinLambdaPtCut &&
			abs(RECOdzAntiSPVmin) < AnalyzerAllSteps::dzAntiSPVminCut &&
			abs(bestRECOAntiS->vz())>AnalyzerAllSteps::vzAntiSInteractionVertexCut && 
			abs(bestRECOAntiS->eta())>AnalyzerAllSteps::antiSEtaCut) histos_th1f["h_RECOSignal_AntiS_BackgroundCuts_cutFlow"]->Fill(17);	



	}

	histos_th1f["h_GENRECO_All_AntiS_pt"]->Fill(GENAntiS->pt());	
	histos_th1f["h_GENRECO_All_AntiS_eta"]->Fill(GENAntiS->eta());	
	histos_th1f["h_GENRECO_All_AntiS_phi"]->Fill(GENAntiS->phi());	
	histos_th2f["h2_GENRECO_All_AntiS_vx_vy"]->Fill(GENAntiS->vx(),GENAntiS->vy());
	histos_th2f["h2_GENRECO_All_AntiS_lxy_vz"]->Fill(GENLxy,GENAntiS->vz());
	histos_th1f["h_GENRECO_All_AntiS_lxy"]->Fill(GENLxy);	
	histos_th1f["h_GENRECO_All_AntiS_lxy_interactionVertex"]->Fill(GENLxy_interactionVertex);	
	histos_th1f["h_GENRECO_All_AntiS_vz"]->Fill(GENAntiS->vz());	
	histos_th1f["h_GENRECO_All_AntiS_vz_interactionVertex"]->Fill(GENAntiSInteractionVertex.Z());	
	histos_th1f["h_GENRECO_All_AntiS_dxy"]->Fill(GENdxy);	
  }
}
*/


void FlatTreeProducerTracking::InitTracking()
{

	_tp_pt.clear();
	_tp_eta.clear();
	_tp_phi.clear();
	_tp_pz.clear();

	_tp_Lxy_beamspot.clear();
	_tp_vz_beamspot.clear();
	_tp_dxy_beamspot.clear();
	_tp_dz_beamspot.clear();
	
	_tp_numberOfTrackerLayers.clear();
	_tp_charge.clear();

	_tp_reconstructed.clear();
	_tp_isAntiSTrack.clear();

	_tp_etaOfGrandMotherAntiS.clear();

	_matchedTrack_pt.clear();
	_matchedTrack_eta.clear();
	_matchedTrack_phi.clear();
	_matchedTrack_pz.clear();

	_matchedTrack_chi2.clear();
	_matchedTrack_ndof.clear();
	_matchedTrack_charge.clear();
	_matchedTrack_dxy_beamspot.clear();
	_matchedTrack_dz_beamspot.clear();
	
	_matchedTrack_trackQuality.clear();
	_matchedTrack_isLooper.clear();	



}


void FlatTreeProducerTracking::endJob()
{
}

void
FlatTreeProducerTracking::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void
FlatTreeProducerTracking::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
FlatTreeProducerTracking::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
FlatTreeProducerTracking::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatTreeProducerTracking::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

FlatTreeProducerTracking::~FlatTreeProducerTracking()
{
}


DEFINE_FWK_MODULE(FlatTreeProducerTracking);
