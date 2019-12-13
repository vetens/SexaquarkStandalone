#include "../interface/FlatTreeProducerTracking.h"
#include <typeinfo>

// pdg mass constants, I copies (sorry, I know its ugly) a part of the V0Fitter code below to check where the V0 fitting fails for Sbar candidates. The constants defined
// here are used in this piece of code
namespace {
   const double piMass = 0.13957018;
   const double piMassSquared = piMass*piMass;
   const double protonMass = 0.938272046;
   const double protonMassSquared = protonMass*protonMass;
   const double kShortMass = 0.497614;
   const double lambdaMass = 1.115683;
}

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

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
   //the below configurations are also for the V0Fitting

   useVertex_ = pset.getParameter<bool>("useVertex");
   // whether to reconstruct KShorts
   doKShorts_ = pset.getParameter<bool>("doKShorts");
   // whether to reconstruct Lambdas
   doLambdas_ = pset.getParameter<bool>("doLambdas");

   // cuts on initial track selection
   tkChi2Cut_ = pset.getParameter<double>("tkChi2Cut");
   tkNHitsCut_ = pset.getParameter<int>("tkNHitsCut");
   tkPtCut_ = pset.getParameter<double>("tkPtCut");
   tkIPSigXYCut_ = pset.getParameter<double>("tkIPSigXYCut");
   tkIPSigZCut_ = pset.getParameter<double>("tkIPSigZCut");
   
   // cuts on vertex
   vtxChi2Cut_ = pset.getParameter<double>("vtxChi2Cut");
   vtxDecaySigXYZCut_ = pset.getParameter<double>("vtxDecaySigXYZCut");
   vtxDecaySigXYCut_ = pset.getParameter<double>("vtxDecaySigXYCut");
   // miscellaneous cuts
   tkDCACut_ = pset.getParameter<double>("tkDCACut");
   mPiPiCut_ = pset.getParameter<double>("mPiPiCut");
   innerHitPosCut_ = pset.getParameter<double>("innerHitPosCut");
   cosThetaXYCut_ = pset.getParameter<double>("cosThetaXYCut");
   cosThetaXYZCut_ = pset.getParameter<double>("cosThetaXYZCut");
   // cuts on the V0 candidate mass
   kShortMassCut_ = pset.getParameter<double>("kShortMassCut");
   lambdaMassCut_ = pset.getParameter<double>("lambdaMassCut");
}


void FlatTreeProducerTracking::beginJob() {

	// Initialize when class is created
        edm::Service<TFileService> fs ; 

	//PV info
	_tree_PV = fs->make <TTree>("FlatTreePV","treePV");
	_tree_PV->Branch("_goodPVxPOG",&_goodPVxPOG);
	_tree_PV->Branch("_goodPVyPOG",&_goodPVyPOG);
	_tree_PV->Branch("_goodPVzPOG",&_goodPVzPOG);
	_tree_PV->Branch("_goodPV_weightPU",&_goodPV_weightPU);

	//counting the number of reco antiS and the total number of GEN antiS
	_tree_counter = fs->make <TTree>("FlatTreeCounter","treeCounter");
	_tree_counter->Branch("_nGENAntiS",&_nGENAntiS);
	_tree_counter->Branch("_nRECOAntiS",&_nRECOAntiS);

	//tree for all the tracks, normally I don't use this as it way too heavy (there are a looooot of tracks)	
	_tree_tracks = fs->make <TTree>("FlatTreeTracks","treeTracks");
	//GEN (trackingparticle) level info
	_tree_tracks->Branch("_tp_pt",&_tp_pt);
	_tree_tracks->Branch("_tp_eta",&_tp_eta);
	_tree_tracks->Branch("_tp_phi",&_tp_phi);
	_tree_tracks->Branch("_tp_pz",&_tp_pz);
	_tree_tracks->Branch("_tp_Lxy_beamspot",&_tp_Lxy_beamspot);
	_tree_tracks->Branch("_tp_vz_beamspot",&_tp_vz_beamspot);
	_tree_tracks->Branch("_tp_dxy_beamspot",&_tp_dxy_beamspot);
	_tree_tracks->Branch("_tp_dz_beamspot",&_tp_dz_beamspot);
	_tree_tracks->Branch("_tp_numberOfTrackerHits",&_tp_numberOfTrackerHits);
	_tree_tracks->Branch("_tp_charge",&_tp_charge);
	_tree_tracks->Branch("_tp_reconstructed",&_tp_reconstructed);
	_tree_tracks->Branch("_tp_isAntiSTrack",&_tp_isAntiSTrack);
	_tree_tracks->Branch("_tp_etaOfGrandMotherAntiS",&_tp_etaOfGrandMotherAntiS);
	//RECO (matched to trackingparticle) level info
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

	//tree to store the tps in an Sbar event, so for each branch there will be 7 entries in the vector: 0th is the Sbar, 1st is the Ks, 2nd is the Lambda, 
	//3rd pi+ from Ks, 4th pi- from Ks, 5th pi+ from antiLambda, 6th pi- from antiproton
	_tree_tpsAntiS = fs->make <TTree>("FlatTreeTpsAntiS","tree_tpsAntiS");
	//GEN (trackingparticle) level info
	_tree_tpsAntiS->Branch("_tpsAntiS_type",&_tpsAntiS_type);
	_tree_tpsAntiS->Branch("_tpsAntiS_pdgId",&_tpsAntiS_pdgId);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestDeltaRWithRECO",&_tpsAntiS_bestDeltaRWithRECO);
	_tree_tpsAntiS->Branch("_tpsAntiS_deltaLInteractionVertexAntiSmin",&_tpsAntiS_deltaLInteractionVertexAntiSmin);
	_tree_tpsAntiS->Branch("_tpsAntiS_mass",&_tpsAntiS_mass);
	_tree_tpsAntiS->Branch("_tpsAntiS_pt",&_tpsAntiS_pt);
	_tree_tpsAntiS->Branch("_tpsAntiS_eta",&_tpsAntiS_eta);
	_tree_tpsAntiS->Branch("_tpsAntiS_phi",&_tpsAntiS_phi);
	_tree_tpsAntiS->Branch("_tpsAntiS_pz",&_tpsAntiS_pz);
	_tree_tpsAntiS->Branch("_tpsAntiS_Lxy_beampipeCenter",&_tpsAntiS_Lxy_beampipeCenter);
	_tree_tpsAntiS->Branch("_tpsAntiS_Lxy_beamspot",&_tpsAntiS_Lxy_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_vz",&_tpsAntiS_vz);
	_tree_tpsAntiS->Branch("_tpsAntiS_vz_beamspot",&_tpsAntiS_vz_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_dxy_beamspot",&_tpsAntiS_dxy_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_dz_beamspot",&_tpsAntiS_dz_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_dz_AntiSCreationVertex",&_tpsAntiS_dz_AntiSCreationVertex);
	_tree_tpsAntiS->Branch("_tpsAntiS_dxyTrack_beamspot",&_tpsAntiS_dxyTrack_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_dzTrack_beamspot",&_tpsAntiS_dzTrack_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_numberOfTrackerHits",&_tpsAntiS_numberOfTrackerHits);
	_tree_tpsAntiS->Branch("_tpsAntiS_charge",&_tpsAntiS_charge);
	_tree_tpsAntiS->Branch("_tpsAntiS_reconstructed",&_tpsAntiS_reconstructed);
	//RECO (matched to trackingparticle) level info
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_mass",&_tpsAntiS_bestRECO_mass);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_massMinusNeutron",&_tpsAntiS_bestRECO_massMinusNeutron);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_pt",&_tpsAntiS_bestRECO_pt);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_eta",&_tpsAntiS_bestRECO_eta);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_phi",&_tpsAntiS_bestRECO_phi);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_pz",&_tpsAntiS_bestRECO_pz);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_Lxy_beampipeCenter",&_tpsAntiS_bestRECO_Lxy_beampipeCenter);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_error_Lxy_beampipeCenter",&_tpsAntiS_bestRECO_error_Lxy_beampipeCenter);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_Lxy_beamspot",&_tpsAntiS_bestRECO_Lxy_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_error_Lxy_beamspot",&_tpsAntiS_bestRECO_error_Lxy_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_vz",&_tpsAntiS_bestRECO_vz);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_vz_beamspot",&_tpsAntiS_bestRECO_vz_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_dxy_beamspot",&_tpsAntiS_bestRECO_dxy_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_dz_beamspot",&_tpsAntiS_bestRECO_dz_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_dxyTrack_beamspot",&_tpsAntiS_bestRECO_dxyTrack_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_dzTrack_beamspot",&_tpsAntiS_bestRECO_dzTrack_beamspot);
	_tree_tpsAntiS->Branch("_tpsAntiS_bestRECO_charge",&_tpsAntiS_bestRECO_charge);
	_tree_tpsAntiS->Branch("_tpsAntiS_returnCodeV0Fitter",&_tpsAntiS_returnCodeV0Fitter);
	_tree_tpsAntiS->Branch("_tpsAntiS_event_weighting_factor",&_tpsAntiS_event_weighting_factor);
	_tree_tpsAntiS->Branch("_tpsAntiS_event_weighting_factorPU",&_tpsAntiS_event_weighting_factorPU);

}

void FlatTreeProducerTracking::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

 
  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);
  const reco::BeamSpot* theBeamSpot = h_bs.product();

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

  edm::ESHandle<MagneticField> theMagneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagneticFieldHandle);
  const MagneticField* theMagneticField = theMagneticFieldHandle.product();


  //save some info on the PVs
  InitPV();
  unsigned int nGoodPV = 0;
  if(h_offlinePV.isValid()){

	//first count the number of good vertices
	for(unsigned int i = 0; i < h_offlinePV->size(); i++){
		double r = sqrt(h_offlinePV->at(i).x()*h_offlinePV->at(i).x()+h_offlinePV->at(i).y()*h_offlinePV->at(i).y());
                if(h_offlinePV->at(i).ndof() > 4 && abs(h_offlinePV->at(i).z()) < 24 && r < 2)nGoodPV++;
	}

	//now that you know the good number of vertices store the location of the vertex and the reweighing factor (you need the nGoodPV to calculate the weighing factor)
	for(unsigned int i = 0; i < h_offlinePV->size(); i++){
		double r = sqrt(h_offlinePV->at(i).x()*h_offlinePV->at(i).x()+h_offlinePV->at(i).y()*h_offlinePV->at(i).y());
                if(h_offlinePV->at(i).ndof() > 4 && abs(h_offlinePV->at(i).z()) < 24 && r < 2){
			_goodPVxPOG.push_back(h_offlinePV->at(i).x());
                        _goodPVyPOG.push_back(h_offlinePV->at(i).y());
                        _goodPVzPOG.push_back(h_offlinePV->at(i).z());
			double weightPU = 0.;
			if(nGoodPV < AnalyzerAllSteps::v_mapPU.size()) weightPU = AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[nGoodPV], h_offlinePV->at(i).z());
        		_goodPV_weightPU.push_back(weightPU);
		}
	}
	
  }
  _tree_PV->Fill();

  //beamspot
  TVector3 beamspot(0,0,0);
  TVector3 beamspotVariance(0,0,0);
  reco::BeamSpot::Point beamspotPoint;
  if(h_bs.isValid()){  
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));		
	reco::BeamSpot::Point beamspotPoint(h_bs->position());
	std::cout << "beamspot location: " << h_bs->x0() << ", " << h_bs->y0() << ", "<< h_bs->z0() << std::endl;
	std::cout << "beamspot errors: " << h_bs->x0Error() << ", "<< h_bs->y0Error() << ", "<< h_bs->z0Error() << std::endl;
  }
  else{
	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!beamspot collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  }

  if(!h_generalTracks.isValid()){ std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!generalTracks collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;}
  if(!h_TP.isValid()){ std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!trackingParticles collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;}
  if(!h_trackAssociator.isValid()){ std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!trackAssociator collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;}

  //loop over the gen particles, check for this antiS if there are any antiS with the same eta, so duplicates
  //save the number of duplicates in a vector of vectors. Each vector has  as a first element the eta of the antiS and 2nd element the # of antiS with this eta. 
  vector<vector<float>> v_antiS_momenta_and_itt; 
  for(unsigned int i = 0; i < h_genParticles->size(); ++i){
	const reco::Candidate * genParticle = &h_genParticles->at(i);

  	if(genParticle->pdgId() != AnalyzerAllSteps::pdgIdAntiS) continue;

	int duplicateIt = -1; 
	for(unsigned int j = 0; j<v_antiS_momenta_and_itt.size();j++){//check if this antiS is already in the list.
		if(v_antiS_momenta_and_itt[j][0] == genParticle->eta()){//you found a duplicate
			duplicateIt = j;
		}
	}

	if(duplicateIt>-1){	
			v_antiS_momenta_and_itt[duplicateIt][1]++;
	}
	else{//this is a new antiS
		double weight_PU = 0.;
		if(nGoodPV < AnalyzerAllSteps::v_mapPU.size()) weight_PU = AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[nGoodPV],genParticle->vz());
		nTotalUniqueGenS_weighted = nTotalUniqueGenS_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;
		nTotalUniqueGenS_Nonweighted = nTotalUniqueGenS_Nonweighted + 1;
		vector<float> dummyVec; 
		dummyVec.push_back(genParticle->eta());
		dummyVec.push_back(1.);
		v_antiS_momenta_and_itt.push_back(dummyVec);


	}
  }
  std::cout << "In this event found " << v_antiS_momenta_and_itt.size() << " unique AntiS, with following #duplicates: " << std::endl;
  for(unsigned int j = 0; j<v_antiS_momenta_and_itt.size();j++){
	std::cout << v_antiS_momenta_and_itt[j][1] << " with eta " << v_antiS_momenta_and_itt[j][0] << std::endl;
  }
//evaluate tracking performance, the below part works, but normally I do not use it, because this tree is very heavy
/*  if(h_generalTracks.isValid() && h_TP.isValid() && h_trackAssociator.isValid()){
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
        
	//fill the tree for all charged particles, so this will also be charged particles from the antiS, but you save in the trees which ones are from the antiS, but in this tree one entry is just one track, so you loose the knowledge of which tp belong to the same antiS, that is why I fill also a second tree as below.
	if(tp.charge() != 0 )FillTreesTracks(tp, beamspot, nPVs, matchedTrackPointer, matchingTrackFound,tpIsGrandDaughterAntiS);

	  //now in the above two catagories you are missing all the tracks from the antiS granddaughters, so when you encounter an antiS you should still have to fill four it's 4 potential granddaughters the histograms and also one where you count how many antiS have four granddaughters reconstructed
	  //if(tp.pdgId() == AnalyzerAllSteps::pdgIdAntiS)FillHistosAntiSTracks(tp, beamspot, TPColl,  h_TP, h_trackAssociator, h_generalTracks, h_V0Ks, h_V0L);
	}
  }
*/

//fill a second tree: one entry in this tree will have the parameters of the antiS, the parameters of the daughters and the parameters of the granddaughters. And will check if they were reconstructed or not. 
  if(h_generalTracks.isValid() && h_TP.isValid() && h_trackAssociator.isValid() && h_V0Ks.isValid() && h_V0L.isValid()){

	//to do the trackmatching on hits
	reco::SimToRecoCollection simRecCollL;
	reco::SimToRecoCollection const * simRecCollP=nullptr;
	simRecCollL= std::move(h_trackAssociator->associateSimToReco(h_generalTracks,h_TP));
	simRecCollP= &simRecCollL;
	reco::SimToRecoCollection const & simRecColl= *simRecCollP;

	//some counters
	int nUniqueAntiSInThisEvent = 0;
	int nUniqueAntiSWithCorrectGranddaughtersThisEvent = 0;
	int nUniqueAntiSWithCorrectGranddaughtersRECONSTRUCTEDThisEvent= 0;
	for(size_t i=0; i<TPColl.size(); ++i) {
	
		const TrackingParticle& tp = TPColl[i];
		if(tp.pdgId() != AnalyzerAllSteps::pdgIdAntiS)continue;

		//first check if the antiS is a duplicate: 
		tv_iterator antiS_firstDecayVertex = tp.decayVertices_begin();
		double antisDecayVx = (**antiS_firstDecayVertex).position().X(); double antisDecayVy = (**antiS_firstDecayVertex).position().Y(); double antisDecayVz = (**antiS_firstDecayVertex).position().Z();
		//the duplicate antiS from the looping have their creation vertex at the same location as the production vertex, so exclude these
		if(tp.vx()==antisDecayVx && tp.vy() == antisDecayVy && tp.vz() == antisDecayVz) continue;
		nUniqueAntiSInThisEvent++;
		totalNumberOfUniqueAntiS_FromTrackingParticles++;		

		int returnCodeFillTreesAntiSAndDaughters = FillTreesAntiSAndDaughters(tp, beamspot, beamspotPoint, beamspotVariance, nPVs, h_generalTracks, h_TP, h_trackAssociator, h_V0Ks, h_V0L, h_sCands, TPColl, simRecColl, theBeamSpot, theMagneticField,nGoodPV);
 		if(returnCodeFillTreesAntiSAndDaughters == 0 || returnCodeFillTreesAntiSAndDaughters == 1) nUniqueAntiSWithCorrectGranddaughtersThisEvent++;
 		if(returnCodeFillTreesAntiSAndDaughters == 1) nUniqueAntiSWithCorrectGranddaughtersRECONSTRUCTEDThisEvent++;

	}
	std::cout << "***************************************************************************************" << std::endl;
	std::cout << "---->number of unique AntiS in this event: " << nUniqueAntiSInThisEvent << " <------" << std::endl;
	std::cout << "---->number of unique AntiS with correct granddaughters in this event: " << nUniqueAntiSWithCorrectGranddaughtersThisEvent << " <------" << std::endl;
	if(nUniqueAntiSInThisEvent == 1 && nUniqueAntiSWithCorrectGranddaughtersThisEvent == 1 && nUniqueAntiSWithCorrectGranddaughtersRECONSTRUCTEDThisEvent ==0) std::cout << "---->Found an event with only one unique AntiS and this AntiS goes to the correct granddaughters and is NOT reconstructed <------" << std::endl;
	if(nUniqueAntiSInThisEvent == 1 && nUniqueAntiSWithCorrectGranddaughtersThisEvent == 1 && nUniqueAntiSWithCorrectGranddaughtersRECONSTRUCTEDThisEvent ==1) std::cout << "---->Found an event with only one unique AntiS and this AntiS goes to the correct granddaughters and is reconstructed <------" << std::endl;
	std::cout << "***************************************************************************************" << std::endl;
  }
  else{
	std::cout << "one of the collections for filling the tree for the antiS related particles is not valid:" << std::endl;
	std::cout << "h_generalTracks.isValid " << h_generalTracks.isValid() << std::endl;
	std::cout << "h_TP.isValid " << h_TP.isValid() << std::endl;
	std::cout << "h_trackAssociator.isValid " << h_trackAssociator.isValid() << std::endl;
	std::cout << "h_V0Ks.isValid " << h_V0Ks.isValid() << std::endl;
	std::cout << "h_V0L.isValid " << h_V0L.isValid() << std::endl;
  }
 

} //end of analyzer

//for all possible tracks
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

	_tp_numberOfTrackerHits.push_back(tp.numberOfTrackerHits());
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

//Fill tree with info on the Sbar, its daughters and granddaughters
int FlatTreeProducerTracking::FillTreesAntiSAndDaughters(const TrackingParticle& tp, TVector3 beamspot, reco::BeamSpot::Point beamspotPoint,  TVector3 beamspotVariance, int nPVs, edm::Handle<View<reco::Track>> h_generalTracks, edm::Handle<TrackingParticleCollection> h_TP, edm::Handle< reco::TrackToTrackingParticleAssociator> h_trackAssociator, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L, edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands, TrackingParticleCollection const & TPColl, reco::SimToRecoCollection const & simRecColl, const reco::BeamSpot* theBeamSpot, const MagneticField* theMagneticField, unsigned int nGoodPV ){

	//loop through the trackingparticles to find the antiS daughters and granddaughters and save them so that you can compare later the dauhgters 
	//to the V0 collections using deltaR matching and the granddaughters to the track collection using track matching on hits

	//to save the trackingparticles
	TrackingParticle tp_Ks;
	TrackingParticle tp_AntiLambda;
	TrackingParticle tp_Ks_posPion;
	TrackingParticle tp_Ks_negPion;
	TrackingParticle tp_AntiLambda_posPion;
	TrackingParticle tp_AntiLambda_AntiProton;
	
	int tp_it_Ks_posPion = -1;
	int tp_it_Ks_negPion = -1;
	int tp_it_AntiLambda_posPion = -1;
	int tp_it_AntiLambda_AntiProton = -1;

	int numberOfGranddaughtersFoundThisEvent = 0;
	
	tv_iterator antiS_firstDecayVertex = tp.decayVertices_begin();
	double antisDecayVx = (**antiS_firstDecayVertex).position().X(); double antisDecayVy = (**antiS_firstDecayVertex).position().Y(); double antisDecayVz = (**antiS_firstDecayVertex).position().Z();

	//the trackingparticles do not have mother-daughter relations stored, so have to go through all tps and try to find onderlapping decay and creation vertices.
        for(size_t j=0; j<TPColl.size(); ++j) {//now find the daughters which have a production vertex = decay vertex of the antiS

                const TrackingParticle& tp_daughter = TPColl[j];

                if(abs(tp_daughter.pdgId()) == AnalyzerAllSteps::pdgIdKs || tp_daughter.pdgId() == AnalyzerAllSteps::pdgIdAntiLambda){//daughter has to be a Ks or Lambda

			double daughterProductionVx = tp_daughter.vx(); double daughterProductionVy = tp_daughter.vy(); double daughterProductionVz = tp_daughter.vz();

                        if(antisDecayVx == daughterProductionVx && antisDecayVy == daughterProductionVy &&  antisDecayVz == daughterProductionVz){//daughter prod vertex has to match the mother decay

				//save the daughter:
				if( abs(tp_daughter.pdgId()) == AnalyzerAllSteps::pdgIdKs ) tp_Ks = tp_daughter;
				if( tp_daughter.pdgId() == AnalyzerAllSteps::pdgIdAntiLambda ) tp_AntiLambda = tp_daughter;

				tv_iterator tp_daughter_firstDecayVertex = tp_daughter.decayVertices_begin();
				double daughterdecayVx = (**tp_daughter_firstDecayVertex).position().X(); double daughterdecayVy = (**tp_daughter_firstDecayVertex).position().Y(); double daughterdecayVz = (**tp_daughter_firstDecayVertex).position().Z();

                                for(unsigned int k=0; k<TPColl.size(); ++k) {

                                        const TrackingParticle& tp_granddaughter = TPColl[k];

                                        if(tp_granddaughter.pdgId() == AnalyzerAllSteps::pdgIdPosPion || tp_granddaughter.pdgId() == AnalyzerAllSteps::pdgIdNegPion || tp_granddaughter.pdgId() == AnalyzerAllSteps::pdgIdAntiProton){
                                                double granddaughterProductionVx = tp_granddaughter.vx(); double granddaughterProductionVy = tp_granddaughter.vy();double granddaughterProductionVz= tp_granddaughter.vz();

                                                if(daughterdecayVx == granddaughterProductionVx && daughterdecayVy == granddaughterProductionVy && daughterdecayVz == granddaughterProductionVz){
							//save the granddaughters
							if(abs(tp_daughter.pdgId()) == AnalyzerAllSteps::pdgIdKs && tp_granddaughter.pdgId() == AnalyzerAllSteps::pdgIdPosPion){ tp_Ks_posPion = tp_granddaughter;numberOfGranddaughtersFoundThisEvent++;tp_it_Ks_posPion = k;}
							if(abs(tp_daughter.pdgId()) == AnalyzerAllSteps::pdgIdKs && tp_granddaughter.pdgId() == AnalyzerAllSteps::pdgIdNegPion){ tp_Ks_negPion = tp_granddaughter;numberOfGranddaughtersFoundThisEvent++;tp_it_Ks_negPion = k;}
							if( tp_daughter.pdgId() == AnalyzerAllSteps::pdgIdAntiLambda && tp_granddaughter.pdgId() == AnalyzerAllSteps::pdgIdPosPion){ tp_AntiLambda_posPion = tp_granddaughter;numberOfGranddaughtersFoundThisEvent++;tp_it_AntiLambda_posPion = k;}
							if( tp_daughter.pdgId() == AnalyzerAllSteps::pdgIdAntiLambda && tp_granddaughter.pdgId() == AnalyzerAllSteps::pdgIdAntiProton){ tp_AntiLambda_AntiProton = tp_granddaughter;numberOfGranddaughtersFoundThisEvent++;tp_it_AntiLambda_AntiProton = k;}
                                                }
                                        }//end if good granddaughter
                                }//end loop over the tp to find granddaughter
                        }// end check if antiS decay vertex matches daughter production vertex
                }//end check for pdgId daughter
        }//end loop over tp to find daughters

	std::cout << "number trackingParticle granddaughters found in this event: " << numberOfGranddaughtersFoundThisEvent << std::endl;

	//now only when there are 4 correct tp granddaughters found you have the chance to actually reconstruct them, so save only those events in the tree
	if(tp_it_Ks_posPion == -1 || tp_it_Ks_negPion == -1 || tp_it_AntiLambda_posPion == -1 || tp_it_AntiLambda_AntiProton == -1) return -1;
	numberOfAntiSWithCorrectGranddaughters++;	

	//but first: check if they were reconstructed.
	//for the Ks and AntiLambda you can do this based on deltaR, however for the AntiS now have to adopt a bit different strategy 
	//because the RECO antiS will not be pointing necessarily to the GEN antiS due to the Fermi momentum of the neutron. 
	//So what we can compare is the interaction vertex of the antiS in RECO and in GEN.
	
	double deltaRminAntiS = 999.;
	double deltaLInteractionVertexAntiSmin = 999.;
	int bestMatchingAntiS = -1;
	bool RECOAntiSFound = false;
	if(h_sCands.isValid()){
		for(size_t j=0; j<h_sCands->size(); ++j) {

			if(h_sCands->at(j).charge() != -1 )continue; //only save antiS
			TVector3 RECOAnitSCreationVertex(h_sCands->at(j).vx(),h_sCands->at(j).vy(),h_sCands->at(j).vz());
        		double Lxy = AnalyzerAllSteps::lxy(beamspot,RECOAnitSCreationVertex);			
			if(Lxy <= AnalyzerAllSteps::MinLxyCut) continue;
			

			double deltaR = AnalyzerAllSteps::deltaR( h_sCands->at(j).phi(), h_sCands->at(j).eta(), tp.phi(), tp.eta() );
			//have to use the vertex of the daughter Ks (at GEN level) as the interaction vertex of the AntiS and compare it to the vertex of the RECO antiS which is the annihilation vertex
			double deltaLInteractionVertexAntiS = sqrt( pow(tp_Ks.vx() - h_sCands->at(j).vx(),2) + pow(tp_Ks.vy() - h_sCands->at(j).vy(),2) +  pow(tp_Ks.vz() - h_sCands->at(j).vz(),2)  ); 
			if(deltaR < deltaRminAntiS){
				 deltaRminAntiS = deltaR; 
			}
			if(deltaLInteractionVertexAntiS < deltaLInteractionVertexAntiSmin){
				deltaLInteractionVertexAntiSmin = deltaLInteractionVertexAntiS;
				bestMatchingAntiS = j;	
			}

		}
	}

	//For the Ks, we can use the deltaR
	double deltaRminKs = 999.;
	int bestMatchingKs = -1;
	bool RECOKsFound = false;
	if(h_V0Ks.isValid()){
		for(size_t j=0; j<h_V0Ks->size(); ++j) {
			double deltaR = AnalyzerAllSteps::deltaR( h_V0Ks->at(j).phi(), h_V0Ks->at(j).eta(), tp_Ks.phi(), tp_Ks.eta() ); 
			if(deltaR < deltaRminKs){
				deltaRminKs = deltaR;
				bestMatchingKs = j;
			}
		}
	}

	//also calculate the 3D distance betwen GEN and RECO decay vertex  as an extra check
	double deltaRminKs_deltaL = 999.;
	if(bestMatchingKs > -1) deltaRminKs_deltaL = sqrt( pow(tp_Ks_posPion.vx() - h_V0Ks->at(bestMatchingKs).vx(),2) + pow(tp_Ks_posPion.vy() - h_V0Ks->at(bestMatchingKs).vy(),2) + pow(tp_Ks_posPion.vz() - h_V0Ks->at(bestMatchingKs).vz(),2) );
	//now use both the deltaR and deltaL cut to define a matching criterium
	if( deltaRminKs <  AnalyzerAllSteps::deltaRCutV0RECOKs && deltaRminKs_deltaL < AnalyzerAllSteps::deltaLCutV0RECOKs ) RECOKsFound = true;

	//For the Lambdabar, use deltaR
	double deltaRminAntiL = 999.;
	int bestMatchingAntiL = -1;
	bool RECOAntiLambdaFound = false;
	if(h_V0L.isValid()){
		for(size_t j=0; j<h_V0L->size(); ++j) {
			double deltaR = AnalyzerAllSteps::deltaR( h_V0L->at(j).phi(), h_V0L->at(j).eta(), tp_AntiLambda.phi(), tp_AntiLambda.eta() ); 
			if(deltaR < deltaRminAntiL){
				deltaRminAntiL = deltaR;
				bestMatchingAntiL = j;
			}
		}
	}
	//also calculate the 3D distance betwen GEN and RECO decay vertex  as an extra check
	double deltaRminAntiL_deltaL = 999.;
	if(bestMatchingAntiL > -1) deltaRminAntiL_deltaL = sqrt( pow(tp_AntiLambda_posPion.vx() - h_V0L->at(bestMatchingAntiL).vx(),2) + pow(tp_AntiLambda_posPion.vy() - h_V0L->at(bestMatchingAntiL).vy(),2) + pow(tp_AntiLambda_posPion.vz() - h_V0L->at(bestMatchingAntiL).vz(),2) );
	//now use both the deltaR and deltaL cut to define a matching criterium
	if( deltaRminAntiL <  AnalyzerAllSteps::deltaRCutV0RECOLambda && deltaRminAntiL_deltaL < AnalyzerAllSteps::deltaLCutV0RECOLambda) RECOAntiLambdaFound = true;


	//for the granddaughters you have to do the matching based on hits

	//For the Ks PosPion:
	bool RECOKsPosPionFound = false;
	const reco::Track *matchedTrackPointer_KsPosPion = nullptr;

	TrackingParticleRef tpr_KsPosPion(h_TP,tp_it_Ks_posPion);

	if(simRecColl.find(tpr_KsPosPion) != simRecColl.end()){
	        auto const & rt = simRecColl[tpr_KsPosPion];
	        if (rt.size()!=0) {
	          matchedTrackPointer_KsPosPion = rt.begin()->first.get();
	          RECOKsPosPionFound  = true;
	        }
	}
	else{
	          RECOKsPosPionFound = false;
	}


	//For the Ks NegPion:
	bool RECOKsNegPionFound = false;
	const reco::Track *matchedTrackPointer_KsNegPion = nullptr;

	TrackingParticleRef tpr_KsNegPion(h_TP,tp_it_Ks_negPion);

	if(simRecColl.find(tpr_KsNegPion) != simRecColl.end()){
	        auto const & rt = simRecColl[tpr_KsNegPion];
	        if (rt.size()!=0) {
	          matchedTrackPointer_KsNegPion = rt.begin()->first.get();
	          RECOKsNegPionFound  = true;
	        }
	}
	else{
	          RECOKsNegPionFound = false;
	}


	//For the AntiLambda pos pion:
	bool RECOAntiLambda_posPionFound = false;
	const reco::Track *matchedTrackPointer_AntiLambda_posPion = nullptr;

	TrackingParticleRef tpr_AntiLambda_posPion(h_TP,tp_it_AntiLambda_posPion);

	if(simRecColl.find(tpr_AntiLambda_posPion) != simRecColl.end()){
	        auto const & rt = simRecColl[tpr_AntiLambda_posPion];
	        if (rt.size()!=0) {
	          matchedTrackPointer_AntiLambda_posPion = rt.begin()->first.get();
	          RECOAntiLambda_posPionFound  = true;
	        }
	}
	else{
	          RECOAntiLambda_posPionFound = false;
	}


	//For the AntiLambda anti proton:
	bool RECOAntiLambda_AntiProtonFound = false;
	const reco::Track *matchedTrackPointer_AntiLambda_AntiProton = nullptr;

	TrackingParticleRef tpr_AntiLambda_AntiProton(h_TP,tp_it_AntiLambda_AntiProton);

	if(simRecColl.find(tpr_AntiLambda_AntiProton) != simRecColl.end()){
	        auto const & rt = simRecColl[tpr_AntiLambda_AntiProton];
	        if (rt.size()!=0) {
	          matchedTrackPointer_AntiLambda_AntiProton = rt.begin()->first.get();
	          RECOAntiLambda_AntiProtonFound  = true;
	        }
	}
	else{
	          RECOAntiLambda_AntiProtonFound = false;
	}

	TVector3 AntiSCreationVertex(tp.vx(),tp.vy(),tp.vz());

	//now define the requirement for reconstructed antiS. This counting also uses the fact if matches between tracks and V0s were found
	if(RECOKsPosPionFound && RECOKsNegPionFound && RECOAntiLambda_posPionFound && RECOAntiLambda_AntiProtonFound &&
		RECOKsFound && RECOAntiLambdaFound &&
		  deltaLInteractionVertexAntiSmin < AnalyzerAllSteps::deltaLCutInteractionVertexAntiSmin && deltaRminAntiS < AnalyzerAllSteps::deltaRCutRECOAntiS){ 
			RECOAntiSFound = true;
	}


	//now if you have both tracks of a V0 in theory you can reconstruct this V0. Pass the tracks to the V0Fitter and check if the V0 gets reconstructed, 
	//but more importantly if it does not get reconstructed, check where things go wrong
	//for the Ks, first check for both daughter tracks if they were reconstructed if they pass the cuts in the V0Fitter
	int returnCodeV0Fitter_KsPosPion_track = -1;
	if(matchedTrackPointer_KsPosPion) returnCodeV0Fitter_KsPosPion_track =  V0Fitter_trackSelection(matchedTrackPointer_KsPosPion,theBeamSpot);
	
	int returnCodeV0Fitter_KsNegPion_track = -1;
	if(matchedTrackPointer_KsNegPion) returnCodeV0Fitter_KsNegPion_track =  V0Fitter_trackSelection(matchedTrackPointer_KsNegPion,theBeamSpot);


	//now if two tracks from a KS got reconstructed check where the V0 fitting failed: is it both tracks which did not pass the 
	//preselection in the V0Fitter? Is it the trackpair that does not pass a cut? Did the V0 get reconstructed?
	int returnCodeV0Fitter_Ks = -1;
	if(matchedTrackPointer_KsNegPion && matchedTrackPointer_KsPosPion){
		if(returnCodeV0Fitter_KsPosPion_track != 0 && returnCodeV0Fitter_KsNegPion_track != 0) returnCodeV0Fitter_Ks = 50;
		else if(returnCodeV0Fitter_KsPosPion_track != 0) returnCodeV0Fitter_Ks = 51;
		else if(returnCodeV0Fitter_KsNegPion_track != 0) returnCodeV0Fitter_Ks = 52;
		else if(returnCodeV0Fitter_KsPosPion_track == 0 && returnCodeV0Fitter_KsNegPion_track == 0) returnCodeV0Fitter_Ks = V0Fitter(matchedTrackPointer_KsPosPion,matchedTrackPointer_KsNegPion,theBeamSpot, theMagneticField, true, false, false);
	}

	//for the AntiLambda, first check for both daughter tracks if they were reconstructed if they pass the cuts in the V0Fitter
	int returnCodeV0Fitter_AntiLambdaPosPion_track = -1;
	if(matchedTrackPointer_AntiLambda_posPion) returnCodeV0Fitter_AntiLambdaPosPion_track =  V0Fitter_trackSelection(matchedTrackPointer_AntiLambda_posPion,theBeamSpot);

	int returnCodeV0Fitter_AntiLambda_AntiProton = -1;
	if(matchedTrackPointer_AntiLambda_AntiProton) returnCodeV0Fitter_AntiLambda_AntiProton =  V0Fitter_trackSelection(matchedTrackPointer_AntiLambda_AntiProton,theBeamSpot);


	//now if the tracks passed the cuts in the V0 Fitter you can go and look at the V0 itself which will be created from this track
	int returnCodeV0Fitter_AntiLambda = -1;
	if(matchedTrackPointer_AntiLambda_posPion && matchedTrackPointer_AntiLambda_AntiProton){
		if(returnCodeV0Fitter_AntiLambdaPosPion_track != 0 && returnCodeV0Fitter_AntiLambda_AntiProton != 0) returnCodeV0Fitter_AntiLambda = 50;
		else if(returnCodeV0Fitter_AntiLambdaPosPion_track != 0) returnCodeV0Fitter_AntiLambda = 51;
		else if(returnCodeV0Fitter_AntiLambda_AntiProton   != 0) returnCodeV0Fitter_AntiLambda = 52;
		else if(returnCodeV0Fitter_AntiLambdaPosPion_track == 0 && returnCodeV0Fitter_AntiLambda_AntiProton == 0) returnCodeV0Fitter_AntiLambda = V0Fitter(matchedTrackPointer_AntiLambda_posPion,matchedTrackPointer_AntiLambda_AntiProton,theBeamSpot, theMagneticField, false, false, true);
	}


	//now fill the trees for each of the trackingparticles	
	InitTrackingAntiS();	

	//for each of the 7 particles in the game fill first some branches (FillFlatTreeTpsAntiS) which contains the kinematics of the tp, so this is on GEN level. 
	//Then fill a some branches (FillFlatTreeTpsAntiSRECO) containing info on the best matching RECO object. 
	//For the first 3 particles these are VertexCompositeCandidate, for the other 4 these matching RECO objects are tracks. I can only fill this second part of the 
	//tree if indeed a RECO object was found. If a RECO object was not found I should still make a dummy entry in the vector to keep the order in the vector correct, 
	//as entry 0 is the the antiS, 1 is the Ks, 2 is the antiLambda, 3 is the pi plus of the ks, 4 the pi minus of the Ks, 5 the pos pion of the antiLambda, 6 the anti proton of the antiLambda. 
	FillFlatTreeTpsAntiS(beamspot,AntiSCreationVertex,tp,RECOAntiSFound,0,deltaRminAntiS,999,deltaLInteractionVertexAntiSmin, theMagneticField);

	//calculate the weight parameter for the pathlength through the beampipe and the reweighing for the PV
	double weightBeampipe = AnalyzerAllSteps::EventWeightingFactor(tp.theta());
	_tpsAntiS_event_weighting_factor.push_back(weightBeampipe);
	double weightPV = 0.;
	if(nGoodPV < AnalyzerAllSteps::v_mapPU.size()) weightPV = AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[nGoodPV],tp.vz());
	_tpsAntiS_event_weighting_factorPU.push_back(weightPV);

	if(bestMatchingAntiS>-1){//for the antiS just save a few extras, which you do not save for the other particles:

		FillFlatTreeTpsAntiSRECO(beamspot,beamspotPoint,RECOAntiSFound,0,h_sCands->at(bestMatchingAntiS));

		//for the antiS just save a few extras:

		//the invariant mass minus the neutron mass, for any of the other 7 objects this is not relevant
		reco::LeafCandidate::LorentzVector n_(0,0,0,0.939565);
        	double RECO_Smass = (h_sCands->at(bestMatchingAntiS).p4()-n_).mass();
		//the error on the interaction vertex
		double RECOErrorLxy_interactionVertex = AnalyzerAllSteps::std_dev_lxy(h_sCands->at(bestMatchingAntiS).vx(), h_sCands->at(bestMatchingAntiS).vy(), h_sCands->at(bestMatchingAntiS).vertexCovariance(0,0), h_sCands->at(bestMatchingAntiS).vertexCovariance(1,1), beamspot.X(), beamspot.Y(), beamspotVariance.X(), beamspotVariance.Y());
		double RECOErrorLxy_interactionVertex_beampipeCenter = AnalyzerAllSteps::std_dev_lxy(h_sCands->at(bestMatchingAntiS).vx(), h_sCands->at(bestMatchingAntiS).vy(), h_sCands->at(bestMatchingAntiS).vertexCovariance(0,0), h_sCands->at(bestMatchingAntiS).vertexCovariance(1,1), 0., 0., 0., 0.);

		_tpsAntiS_bestRECO_massMinusNeutron.push_back(RECO_Smass);	
		_tpsAntiS_bestRECO_error_Lxy_beamspot.push_back(RECOErrorLxy_interactionVertex);
		_tpsAntiS_bestRECO_error_Lxy_beampipeCenter.push_back(RECOErrorLxy_interactionVertex_beampipeCenter);

	}
	else{
		 _tpsAntiS_bestRECO_massMinusNeutron.push_back(999.);
		 _tpsAntiS_bestRECO_error_Lxy_beamspot.push_back(999.);
		 _tpsAntiS_bestRECO_error_Lxy_beampipeCenter.push_back(999.);
		 FillFlatTreeTpsAntiSRECODummy();
	}

	//Ks
	FillFlatTreeTpsAntiS(beamspot,AntiSCreationVertex,tp_Ks,RECOKsFound,1,deltaRminKs,returnCodeV0Fitter_Ks,deltaRminKs_deltaL,theMagneticField);
	if(bestMatchingKs > -1)FillFlatTreeTpsAntiSRECO(beamspot,beamspotPoint,RECOKsFound,1,h_V0Ks->at(bestMatchingKs));
	else FillFlatTreeTpsAntiSRECODummy();
	//AntiLambda
	FillFlatTreeTpsAntiS(beamspot,AntiSCreationVertex,tp_AntiLambda,RECOAntiLambdaFound,2,deltaRminAntiL,returnCodeV0Fitter_AntiLambda,deltaRminAntiL_deltaL,theMagneticField);
	if(bestMatchingAntiL > -1)FillFlatTreeTpsAntiSRECO(beamspot,beamspotPoint,RECOAntiLambdaFound,2,h_V0L->at(bestMatchingAntiL));
	else FillFlatTreeTpsAntiSRECODummy();
	//Ks_posPion
	FillFlatTreeTpsAntiS(beamspot,AntiSCreationVertex,tp_Ks_posPion,RECOKsPosPionFound,3,-999,returnCodeV0Fitter_KsPosPion_track,999.,theMagneticField);
	if(matchedTrackPointer_KsPosPion)FillFlatTreeTpsAntiSRECO(beamspot,beamspotPoint,RECOKsPosPionFound,3,matchedTrackPointer_KsPosPion);
	else FillFlatTreeTpsAntiSRECODummy();
	//Ks_negPion
	FillFlatTreeTpsAntiS(beamspot,AntiSCreationVertex,tp_Ks_negPion,RECOKsNegPionFound,4,-999,returnCodeV0Fitter_KsNegPion_track,999.,theMagneticField);
	if(matchedTrackPointer_KsNegPion)FillFlatTreeTpsAntiSRECO(beamspot,beamspotPoint,RECOKsNegPionFound,4,matchedTrackPointer_KsNegPion);
	else FillFlatTreeTpsAntiSRECODummy();
	//AntiLambda_posPion
	FillFlatTreeTpsAntiS(beamspot,AntiSCreationVertex,tp_AntiLambda_posPion,RECOAntiLambda_posPionFound,5,-999,returnCodeV0Fitter_AntiLambdaPosPion_track,999.,theMagneticField);
	if(matchedTrackPointer_AntiLambda_posPion)FillFlatTreeTpsAntiSRECO(beamspot,beamspotPoint,RECOAntiLambda_posPionFound,5,matchedTrackPointer_AntiLambda_posPion);
	else FillFlatTreeTpsAntiSRECODummy();
	//AntiLambda_AntiProton
	FillFlatTreeTpsAntiS(beamspot,AntiSCreationVertex,tp_AntiLambda_AntiProton,RECOAntiLambda_AntiProtonFound,6,-999,returnCodeV0Fitter_AntiLambda_AntiProton,999.,theMagneticField);
	if(matchedTrackPointer_AntiLambda_AntiProton)FillFlatTreeTpsAntiSRECO(beamspot,beamspotPoint,RECOAntiLambda_AntiProtonFound,6,matchedTrackPointer_AntiLambda_AntiProton);
	else FillFlatTreeTpsAntiSRECODummy();

	_tree_tpsAntiS->Fill();


	if(RECOAntiSFound) weighedRecoAntiS  += weightBeampipe*weightPV; 
	if(RECOAntiSFound) nonweighedRecoAntiS  = nonweighedRecoAntiS + 1; 
	
	return RECOAntiSFound;
}



void FlatTreeProducerTracking::FillFlatTreeTpsAntiS(TVector3 beamspot, TVector3 AntiSCreationVertex, TrackingParticle trackingParticle, bool RECOFound, int type, double besteDeltaR, int returnCodeV0Fitter, double besteDeltaL, const MagneticField* theMagneticField){

	//some kinematic variables
	TVector3 tpCreationVertex(trackingParticle.vx(),trackingParticle.vy(),trackingParticle.vz());
	TVector3 ZeroZeroZero(0.,0.,0.);
	double Lxy_beampipeCenter = AnalyzerAllSteps::lxy(ZeroZeroZero,tpCreationVertex);
	double Lxy = AnalyzerAllSteps::lxy(beamspot,tpCreationVertex);
	TVector3 tpMomentum(trackingParticle.px(),trackingParticle.py(),trackingParticle.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(tpCreationVertex,tpMomentum,beamspot);
	double dz_beamspot = AnalyzerAllSteps::dz_line_point(tpCreationVertex,tpMomentum,beamspot);
	double dz_AntiSCreationVertex = AnalyzerAllSteps::dz_line_point(tpCreationVertex,tpMomentum,AntiSCreationVertex);

	//calculate the dxy and dz of the tp, but like you would for a real track, that is: extrapolating the track vertex 
	//and momentum to the point of closest approach to the beamspot (as in the example starting at 
	//https://github.com/cms-sw/cmssw/blob/9a33fb13bee1a546877a4b581fa63876043f38f0/SimTracker/TrackHistory/src/TrackClassifier.cc#L161)

	FreeTrajectoryState ftsAtProduction(
	      GlobalPoint(trackingParticle.vx(), trackingParticle.vy(), trackingParticle.vz()),
	      GlobalVector(trackingParticle.px(), trackingParticle.py(), trackingParticle.pz()),
	      TrackCharge(trackingParticle.charge()),
	      theMagneticField);

    	TSCPBuilderNoMaterial tscpBuilder;
    	TrajectoryStateClosestToPoint tsAtClosestApproach = tscpBuilder(ftsAtProduction, GlobalPoint(beamspot.X(),beamspot.Y(),beamspot.Z()));

    	GlobalVector v = tsAtClosestApproach.theState().position() - GlobalPoint(beamspot.X(),beamspot.Y(),beamspot.Z());
    	GlobalVector p = tsAtClosestApproach.theState().momentum();

    	// Simulated dxy
    	double dxySim = -v.x() * sin(p.phi()) + v.y() * cos(p.phi());

    	// Simulated dz
    	double dzSim = v.z() - (v.x() * p.x() + v.y() * p.y()) * p.z() / p.perp2();	

	//the type is just to use in the ROOT viewer to be able to select which particle I want to look at
	_tpsAntiS_type.push_back(type);
	_tpsAntiS_pdgId.push_back(trackingParticle.pdgId());
	_tpsAntiS_bestDeltaRWithRECO.push_back(besteDeltaR);
	_tpsAntiS_deltaLInteractionVertexAntiSmin.push_back(besteDeltaL);
	_tpsAntiS_mass.push_back(trackingParticle.mass());
	_tpsAntiS_pt.push_back(trackingParticle.pt());
	_tpsAntiS_eta.push_back(trackingParticle.eta());
	_tpsAntiS_phi.push_back(trackingParticle.phi());
	_tpsAntiS_pz.push_back(trackingParticle.pz());
	_tpsAntiS_Lxy_beampipeCenter.push_back(Lxy_beampipeCenter);
	_tpsAntiS_Lxy_beamspot.push_back(Lxy);
	_tpsAntiS_vz.push_back(trackingParticle.vz());
	_tpsAntiS_vz_beamspot.push_back(trackingParticle.vz()-beamspot.Z());
	_tpsAntiS_dxy_beamspot.push_back(dxy);
	_tpsAntiS_dz_beamspot.push_back(dz_beamspot);
	_tpsAntiS_dz_AntiSCreationVertex.push_back(dz_AntiSCreationVertex);
	_tpsAntiS_dxyTrack_beamspot.push_back(dxySim);
	_tpsAntiS_dzTrack_beamspot.push_back(dzSim);
	_tpsAntiS_numberOfTrackerHits.push_back(trackingParticle.numberOfTrackerHits());
	_tpsAntiS_charge.push_back(trackingParticle.charge());
	_tpsAntiS_reconstructed.push_back(RECOFound);
	_tpsAntiS_returnCodeV0Fitter.push_back(returnCodeV0Fitter);

}

void FlatTreeProducerTracking::FillFlatTreeTpsAntiSRECO(TVector3 beamspot, reco::BeamSpot::Point beamspotPoint, bool RECOFound,int type, reco::VertexCompositeCandidate bestRECOCompositeCandidate){

	//calculate some kinematic variables
	TVector3 bestRECOCompositeCandidateCreationVertex(bestRECOCompositeCandidate.vx(),bestRECOCompositeCandidate.vy(),bestRECOCompositeCandidate.vz());
	TVector3 ZeroZeroZero(0.,0.,0.);
	double Lxy_beampipeCenter = AnalyzerAllSteps::lxy(ZeroZeroZero,bestRECOCompositeCandidateCreationVertex);
	double Lxy = AnalyzerAllSteps::lxy(beamspot,bestRECOCompositeCandidateCreationVertex);
	TVector3 bestRECOCompositeCandidateMomentum(bestRECOCompositeCandidate.px(),bestRECOCompositeCandidate.py(),bestRECOCompositeCandidate.pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(bestRECOCompositeCandidateCreationVertex,bestRECOCompositeCandidateMomentum,beamspot);
	double dz_beamspot = AnalyzerAllSteps::dz_line_point(bestRECOCompositeCandidateCreationVertex,bestRECOCompositeCandidateMomentum,beamspot);

	//the type is just to use in the ROOT viewer to be able to select which particle I want to look at
	_tpsAntiS_bestRECO_mass.push_back(bestRECOCompositeCandidate.mass());
	_tpsAntiS_bestRECO_pt.push_back(bestRECOCompositeCandidate.pt());
	_tpsAntiS_bestRECO_eta.push_back(bestRECOCompositeCandidate.eta());
	_tpsAntiS_bestRECO_phi.push_back(bestRECOCompositeCandidate.phi());
	_tpsAntiS_bestRECO_pz.push_back(bestRECOCompositeCandidate.pz());
	_tpsAntiS_bestRECO_Lxy_beampipeCenter.push_back(Lxy_beampipeCenter);
	_tpsAntiS_bestRECO_Lxy_beamspot.push_back(Lxy);
	_tpsAntiS_bestRECO_vz.push_back(bestRECOCompositeCandidate.vz());
	_tpsAntiS_bestRECO_vz_beamspot.push_back(bestRECOCompositeCandidate.vz()-beamspot.Z());
	_tpsAntiS_bestRECO_dxy_beamspot.push_back(dxy);
	_tpsAntiS_bestRECO_dz_beamspot.push_back(dz_beamspot);
	//VertexCompositeCandidate does not have a dxy or dz
	_tpsAntiS_bestRECO_dxyTrack_beamspot.push_back(999.);
        _tpsAntiS_bestRECO_dzTrack_beamspot.push_back(999.);
	_tpsAntiS_bestRECO_charge.push_back(bestRECOCompositeCandidate.charge());

}

void FlatTreeProducerTracking::FillFlatTreeTpsAntiSRECO(TVector3 beamspot, reco::BeamSpot::Point beamspotPoint, bool RECOFound,int type, const reco::Track *matchedTrackPointer){
	
	//some kinematic variables
	TVector3 bestRECOCompositeCandidateCreationVertex(matchedTrackPointer->vx(),matchedTrackPointer->vy(),matchedTrackPointer->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,bestRECOCompositeCandidateCreationVertex);
        TVector3 ZeroZeroZero(0.,0.,0.);
        double Lxy_beampipeCenter = AnalyzerAllSteps::lxy(ZeroZeroZero,bestRECOCompositeCandidateCreationVertex);
	TVector3 bestRECOCompositeCandidateMomentum(matchedTrackPointer->px(),matchedTrackPointer->py(),matchedTrackPointer->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(bestRECOCompositeCandidateCreationVertex,bestRECOCompositeCandidateMomentum,beamspot);
	double dz_beamspot = AnalyzerAllSteps::dz_line_point(bestRECOCompositeCandidateCreationVertex,bestRECOCompositeCandidateMomentum,beamspot);
	double dxyTrackbeamspot = matchedTrackPointer->dxy(beamspotPoint);
	double dzTrackbeamspot = matchedTrackPointer->dz(beamspotPoint);

	_tpsAntiS_bestRECO_mass.push_back(999.); //tracks dont have a mass so save a dummy value, I need to save dummy values to keep the order in the vectors
	_tpsAntiS_bestRECO_pt.push_back(matchedTrackPointer->pt());
	_tpsAntiS_bestRECO_eta.push_back(matchedTrackPointer->eta());
	_tpsAntiS_bestRECO_phi.push_back(matchedTrackPointer->phi());
	_tpsAntiS_bestRECO_pz.push_back(matchedTrackPointer->pz());
	_tpsAntiS_bestRECO_Lxy_beampipeCenter.push_back(Lxy_beampipeCenter);
	_tpsAntiS_bestRECO_Lxy_beamspot.push_back(Lxy);
	_tpsAntiS_bestRECO_vz.push_back(matchedTrackPointer->vz());
	_tpsAntiS_bestRECO_vz_beamspot.push_back(matchedTrackPointer->vz()-beamspot.Z());
	_tpsAntiS_bestRECO_dxy_beamspot.push_back(dxy);
	_tpsAntiS_bestRECO_dz_beamspot.push_back(dz_beamspot);
	_tpsAntiS_bestRECO_dxyTrack_beamspot.push_back(dxyTrackbeamspot);
	_tpsAntiS_bestRECO_dzTrack_beamspot.push_back(dzTrackbeamspot);
	_tpsAntiS_bestRECO_charge.push_back(matchedTrackPointer->charge());

}

void FlatTreeProducerTracking::FillFlatTreeTpsAntiSRECODummy(){
	

	_tpsAntiS_bestRECO_mass.push_back(999.); //tracks dont have a mass so save a dummy value, I need to save dummy values to keep the order in the vectors
	_tpsAntiS_bestRECO_pt.push_back(999.);
	_tpsAntiS_bestRECO_eta.push_back(999.);
	_tpsAntiS_bestRECO_phi.push_back(999.);
	_tpsAntiS_bestRECO_pz.push_back(999.);
	_tpsAntiS_bestRECO_Lxy_beampipeCenter.push_back(999.);
	_tpsAntiS_bestRECO_Lxy_beamspot.push_back(999.);
	_tpsAntiS_bestRECO_vz.push_back(999.);
	_tpsAntiS_bestRECO_vz_beamspot.push_back(999.);
	_tpsAntiS_bestRECO_dxy_beamspot.push_back(999.);
	_tpsAntiS_bestRECO_dz_beamspot.push_back(999.);
	_tpsAntiS_bestRECO_dxyTrack_beamspot.push_back(999.);
        _tpsAntiS_bestRECO_dzTrack_beamspot.push_back(999.);
	_tpsAntiS_bestRECO_charge.push_back(999);

}

//this is a stupid copy paste from the V0Fitter code in CMSSSW
int FlatTreeProducerTracking::V0Fitter_trackSelection(const reco::Track *matchedTrackPointer1, const reco::BeamSpot* theBeamSpot){

      const reco::Track* tmpTrack = matchedTrackPointer1;
      double ipsigXY = std::abs(tmpTrack->dxy(*theBeamSpot)/tmpTrack->dxyError());
      double ipsigZ = std::abs(tmpTrack->dz(theBeamSpot->position())/tmpTrack->dzError());

      //save the result of the cut in a string
      std::string cut1,cut2,cut3,cut4,cut5 = "0"; 
      if (tmpTrack->normalizedChi2() >=  tkChi2Cut_) cut1 = "1";
      if (tmpTrack->numberOfValidHits() < tkNHitsCut_) cut2 = "1";
      if (tmpTrack->pt() <= tkPtCut_ ) cut3 = "1";
      if (ipsigXY <= tkIPSigXYCut_) cut4 = "1";
      if (ipsigZ <= tkIPSigZCut_) cut5 = "1";

      std::string cutSummary = cut1+cut2+cut3+cut4+cut5;
      int intCutSummary = std::stoi(cutSummary, nullptr, 2);
      return intCutSummary;

}

//this is a stupid copy paste from the V0Fitter code in CMSSSW
int FlatTreeProducerTracking::V0Fitter(const reco::Track *matchedTrackPointer1, const reco::Track *matchedTrackPointer2, const reco::BeamSpot* theBeamSpot, const MagneticField* theMagneticField, bool isGENKs, bool isGENLambda, bool isGENAntiLambda){
      
        math::XYZPoint referencePos(theBeamSpot->position());

	if(matchedTrackPointer1->charge() == matchedTrackPointer2->charge()) return 1;
	if(abs(matchedTrackPointer1->charge()) != 1) return 2;
	if(abs(matchedTrackPointer2->charge()) != 1) return 3;


	reco::Track positiveTrack;
        reco::Track negativeTrack;
	if(matchedTrackPointer1->charge() == 1) positiveTrack = *matchedTrackPointer1;
	else if(matchedTrackPointer1->charge() == -1) negativeTrack = *matchedTrackPointer1;

	if(matchedTrackPointer2->charge() == 1) positiveTrack = *matchedTrackPointer2;
	else if(matchedTrackPointer2->charge() == -1) negativeTrack = *matchedTrackPointer2;

	reco::TransientTrack posTransTk(positiveTrack, theMagneticField);
	reco::TransientTrack negTransTk(negativeTrack, theMagneticField);

	reco::TransientTrack* posTransTkPtr = &posTransTk; 
	reco::TransientTrack* negTransTkPtr = &negTransTk; 

      // measure distance between tracks at their closest approach
      if (!posTransTkPtr->impactPointTSCP().isValid() || !negTransTkPtr->impactPointTSCP().isValid()) return 4;
      FreeTrajectoryState const & posState = posTransTkPtr->impactPointTSCP().theState();
      FreeTrajectoryState const & negState = negTransTkPtr->impactPointTSCP().theState();
      ClosestApproachInRPhi cApp;
      cApp.calculate(posState, negState);
      if (!cApp.status()) return 5;
      float dca = std::abs(cApp.distance());
      std::cout << "the distance of closest approach of the tracks: " << dca << std::endl;
      if (dca > tkDCACut_) return 6;

      // the POCA should at least be in the sensitive volume
      GlobalPoint cxPt = cApp.crossingPoint();
      if (sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) return 7;

      // the tracks should at least point in the same quadrant
      TrajectoryStateClosestToPoint const & posTSCP = posTransTkPtr->trajectoryStateClosestToPoint(cxPt);
      TrajectoryStateClosestToPoint const & negTSCP = negTransTkPtr->trajectoryStateClosestToPoint(cxPt);
      if (!posTSCP.isValid() || !negTSCP.isValid()) return 8;
      if (posTSCP.momentum().dot(negTSCP.momentum())  < 0) return 9;
     
      // calculate mPiPi
      double totalE = sqrt(posTSCP.momentum().mag2() + piMassSquared) + sqrt(negTSCP.momentum().mag2() + piMassSquared);
      double totalESq = totalE*totalE;
      double totalPSq = (posTSCP.momentum() + negTSCP.momentum()).mag2();
      double mass = sqrt(totalESq - totalPSq);
      if (mass > mPiPiCut_) return 10;

      // Fill the vector of TransientTracks to send to KVF
      std::vector<reco::TransientTrack> transTracks;
      transTracks.reserve(2);
      transTracks.push_back(*posTransTkPtr);
      transTracks.push_back(*negTransTkPtr);

      // create the vertex fitter object and vertex the tracks
      TransientVertex theRecoVertex;
      if (vertexFitter_) {
         KalmanVertexFitter theKalmanFitter(useRefTracks_ == 0 ? false : true);
         theRecoVertex = theKalmanFitter.vertex(transTracks);
      } else if (!vertexFitter_) {
         useRefTracks_ = false;
         AdaptiveVertexFitter theAdaptiveFitter;
         theRecoVertex = theAdaptiveFitter.vertex(transTracks);
      }
      if (!theRecoVertex.isValid()) return 11;
     
      reco::Vertex theVtx = theRecoVertex;
      if (theVtx.normalizedChi2() > vtxChi2Cut_) return 12;
      GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

      // 2D decay significance
      SMatrixSym3D totalCov = theBeamSpot->rotatedCovariance3D() + theVtx.covariance();
      //comment below because always using the beamspot
      //if (useVertex_) totalCov = *theBeamSpot.covariance() + theVtx.covariance();
      
      SVector3 distVecXY(vtxPos.x()-referencePos.x(), vtxPos.y()-referencePos.y(), 0.);
      double distMagXY = ROOT::Math::Mag(distVecXY);
      double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
      if (distMagXY/sigmaDistMagXY < vtxDecaySigXYCut_) return 13;

      // 3D decay significance
      SVector3 distVecXYZ(vtxPos.x()-referencePos.x(), vtxPos.y()-referencePos.y(), vtxPos.z()-referencePos.z());
      double distMagXYZ = ROOT::Math::Mag(distVecXYZ);
      double sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;
      if (distMagXYZ/sigmaDistMagXYZ < vtxDecaySigXYZCut_) return 14;

      //comment the below part because also in my reconstruction I donnot use it
      /*
      // make sure the vertex radius is within the inner track hit radius
      if (innerHitPosCut_ > 0. && positiveTrackRef->innerOk()) {
         reco::Vertex::Point posTkHitPos = positiveTrackRef->innerPosition();
         double posTkHitPosD2 =  (posTkHitPos.x()-referencePos.x())*(posTkHitPos.x()-referencePos.x()) +
            (posTkHitPos.y()-referencePos.y())*(posTkHitPos.y()-referencePos.y());
         if (sqrt(posTkHitPosD2) < (distMagXY - sigmaDistMagXY*innerHitPosCut_)) continue;
      }
      if (innerHitPosCut_ > 0. && negativeTrackRef->innerOk()) {
         reco::Vertex::Point negTkHitPos = negativeTrackRef->innerPosition();
         double negTkHitPosD2 = (negTkHitPos.x()-referencePos.x())*(negTkHitPos.x()-referencePos.x()) +
            (negTkHitPos.y()-referencePos.y())*(negTkHitPos.y()-referencePos.y());
         if (sqrt(negTkHitPosD2) < (distMagXY - sigmaDistMagXY*innerHitPosCut_)) continue;
      }
      */
      std::auto_ptr<TrajectoryStateClosestToPoint> trajPlus;
      std::auto_ptr<TrajectoryStateClosestToPoint> trajMins;
      std::vector<reco::TransientTrack> theRefTracks;
      if (theRecoVertex.hasRefittedTracks()) {
         theRefTracks = theRecoVertex.refittedTracks();
      }

      if (useRefTracks_ && theRefTracks.size() > 1) {
         reco::TransientTrack* thePositiveRefTrack = 0;
         reco::TransientTrack* theNegativeRefTrack = 0;
         for (std::vector<reco::TransientTrack>::iterator iTrack = theRefTracks.begin(); iTrack != theRefTracks.end(); ++iTrack) {
            if (iTrack->track().charge() > 0.) {
               thePositiveRefTrack = &*iTrack;
            } else if (iTrack->track().charge() < 0.) {
               theNegativeRefTrack = &*iTrack;
            }
         }
         if (thePositiveRefTrack == 0 || theNegativeRefTrack == 0) return 15;
         trajPlus.reset(new TrajectoryStateClosestToPoint(thePositiveRefTrack->trajectoryStateClosestToPoint(vtxPos)));
         trajMins.reset(new TrajectoryStateClosestToPoint(theNegativeRefTrack->trajectoryStateClosestToPoint(vtxPos)));
      } else {
         trajPlus.reset(new TrajectoryStateClosestToPoint(posTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
         trajMins.reset(new TrajectoryStateClosestToPoint(negTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
      }

      if (trajPlus.get() == 0 || trajMins.get() == 0 || !trajPlus->isValid() || !trajMins->isValid()) return 16;

      GlobalVector positiveP(trajPlus->momentum());
      GlobalVector negativeP(trajMins->momentum());
      GlobalVector totalP(positiveP + negativeP);

      // 2D pointing angle
      double dx = theVtx.x()-referencePos.x();
      double dy = theVtx.y()-referencePos.y();
      double px = totalP.x();
      double py = totalP.y();
      double angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
      if (angleXY < cosThetaXYCut_) return 17;

      // 3D pointing angle
      double dz = theVtx.z()-referencePos.z();
      double pz = totalP.z();
      double angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
      if (angleXYZ < cosThetaXYZCut_) return 18;

      // calculate total energy of V0 3 ways: assume it's a kShort, a Lambda, or a LambdaBar.
      double piPlusE = sqrt(positiveP.mag2() + piMassSquared);
      double piMinusE = sqrt(negativeP.mag2() + piMassSquared);
      double protonE = sqrt(positiveP.mag2() + protonMassSquared);
      double antiProtonE = sqrt(negativeP.mag2() + protonMassSquared);
      double kShortETot = piPlusE + piMinusE;
      double lambdaEtot = protonE + piMinusE;
      double lambdaBarEtot = antiProtonE + piPlusE;

      // Create momentum 4-vectors for the 3 candidate types
      const reco::Particle::LorentzVector kShortP4(totalP.x(), totalP.y(), totalP.z(), kShortETot);
      const reco::Particle::LorentzVector lambdaP4(totalP.x(), totalP.y(), totalP.z(), lambdaEtot);
      const reco::Particle::LorentzVector lambdaBarP4(totalP.x(), totalP.y(), totalP.z(), lambdaBarEtot);

      reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
      const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
      double vtxChi2(theVtx.chi2());
      double vtxNdof(theVtx.ndof());

      // Create the VertexCompositeCandidate object that will be stored in the Event
      reco::VertexCompositeCandidate* theKshort = nullptr;
      reco::VertexCompositeCandidate* theLambda = nullptr;
      reco::VertexCompositeCandidate* theLambdaBar = nullptr;

      if (doKShorts_) {
         theKshort = new reco::VertexCompositeCandidate(0, kShortP4, vtx, vtxCov, vtxChi2, vtxNdof);
      }
      if (doLambdas_) {
         if (positiveP.mag2() > negativeP.mag2()) {
            theLambda = new reco::VertexCompositeCandidate(0, lambdaP4, vtx, vtxCov, vtxChi2, vtxNdof);
         } else {
            theLambdaBar = new reco::VertexCompositeCandidate(0, lambdaBarP4, vtx, vtxCov, vtxChi2, vtxNdof);
         }
      }

      // Create daughter candidates for the VertexCompositeCandidates
      //comment the below because I donnot have access to the trackref and anyway it is not useful here
/*      reco::RecoChargedCandidate thePiPlusCand(
         1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), piPlusE), vtx);
      thePiPlusCand.setTrack(positiveTrackRef);
      
      reco::RecoChargedCandidate thePiMinusCand(
         -1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), piMinusE), vtx);
      thePiMinusCand.setTrack(negativeTrackRef);
      
      reco::RecoChargedCandidate theProtonCand(
         1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), protonE), vtx);
      theProtonCand.setTrack(positiveTrackRef);

      reco::RecoChargedCandidate theAntiProtonCand(
         -1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), antiProtonE), vtx);
      theAntiProtonCand.setTrack(negativeTrackRef);
*/
      //AddFourMomenta addp4;
      // Store the daughter Candidates in the VertexCompositeCandidates if they pass mass cuts
      if (doKShorts_) {
  //       theKshort->addDaughter(thePiPlusCand);
  //       theKshort->addDaughter(thePiMinusCand);
         theKshort->setPdgId(310);
         //addp4.set(*theKshort);
         if (theKshort->mass() < kShortMass + kShortMassCut_ && theKshort->mass() > kShortMass - kShortMassCut_) {
            //theKshorts.push_back(std::move(*theKshort));
            if(isGENKs)return 0;
         }
         else if(isGENKs)return 19;
         
      }

      if (doLambdas_ && theLambda) {
   //      theLambda->addDaughter(theProtonCand);
   //      theLambda->addDaughter(thePiMinusCand);
         theLambda->setPdgId(3122);
         //addp4.set( *theLambda );
         if (theLambda->mass() < lambdaMass + lambdaMassCut_ && theLambda->mass() > lambdaMass - lambdaMassCut_) {
            //theLambdas.push_back(std::move(*theLambda));
            if(isGENLambda){return 0;}
         }
         else if(isGENLambda){return 20;}
      } 
      else if (doLambdas_ && theLambdaBar) {
   //      theLambdaBar->addDaughter(theAntiProtonCand);
   //      theLambdaBar->addDaughter(thePiPlusCand);
         theLambdaBar->setPdgId(-3122);
         //addp4.set(*theLambdaBar);
         if (theLambdaBar->mass() < lambdaMass + lambdaMassCut_ && theLambdaBar->mass() > lambdaMass - lambdaMassCut_) {
            //theLambdas.push_back(std::move(*theLambdaBar));
            if(isGENAntiLambda)return 0;
         }
         else if(isGENAntiLambda){ return 21;}
      }

      delete theKshort;
      delete theLambda;
      delete theLambdaBar;
      theKshort = theLambda = theLambdaBar = nullptr;

      return 22;
}


void FlatTreeProducerTracking::InitPV()
{
	_goodPVxPOG.clear();
	_goodPVyPOG.clear();
	_goodPVzPOG.clear();
	_goodPV_weightPU.clear();
}

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
	
	_tp_numberOfTrackerHits.clear();
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

void FlatTreeProducerTracking::InitTrackingAntiS(){

        _tpsAntiS_type.clear();
        _tpsAntiS_pdgId.clear();
        _tpsAntiS_bestDeltaRWithRECO.clear();
        _tpsAntiS_deltaLInteractionVertexAntiSmin.clear();
        _tpsAntiS_mass.clear();

	_tpsAntiS_pt.clear();
	_tpsAntiS_eta.clear();
	_tpsAntiS_phi.clear();
	_tpsAntiS_pz.clear();

	_tpsAntiS_Lxy_beampipeCenter.clear();
	_tpsAntiS_Lxy_beamspot.clear();
	_tpsAntiS_vz.clear();
	_tpsAntiS_vz_beamspot.clear();
	_tpsAntiS_dxy_beamspot.clear();
	_tpsAntiS_dz_beamspot.clear();
	_tpsAntiS_dz_AntiSCreationVertex.clear();
	_tpsAntiS_dxyTrack_beamspot.clear();
	_tpsAntiS_dzTrack_beamspot.clear();

	_tpsAntiS_numberOfTrackerHits.clear();
	_tpsAntiS_charge.clear();

	_tpsAntiS_reconstructed.clear();

	_tpsAntiS_bestRECO_mass.clear();
	_tpsAntiS_bestRECO_massMinusNeutron.clear();
	
	_tpsAntiS_bestRECO_pt.clear();
	_tpsAntiS_bestRECO_eta.clear();
	_tpsAntiS_bestRECO_phi.clear();
	_tpsAntiS_bestRECO_pz.clear();

	_tpsAntiS_bestRECO_Lxy_beampipeCenter.clear();
	_tpsAntiS_bestRECO_Lxy_beamspot.clear();
	_tpsAntiS_bestRECO_error_Lxy_beamspot.clear();
	_tpsAntiS_bestRECO_error_Lxy_beampipeCenter.clear();
	_tpsAntiS_bestRECO_vz.clear();
	_tpsAntiS_bestRECO_vz_beamspot.clear();
	_tpsAntiS_bestRECO_dxy_beamspot.clear();
	_tpsAntiS_bestRECO_dz_beamspot.clear();

	_tpsAntiS_bestRECO_dxyTrack_beamspot.clear();
	_tpsAntiS_bestRECO_dzTrack_beamspot.clear();

	_tpsAntiS_bestRECO_charge.clear();
	_tpsAntiS_returnCodeV0Fitter.clear();

	_tpsAntiS_event_weighting_factor.clear();
	_tpsAntiS_event_weighting_factorPU.clear();

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

	std::cout << "number of antiS found which interact (based on the number of trackingParticles which should have their decay vertex different from their creation vertex): " << totalNumberOfUniqueAntiS_FromTrackingParticles<< std::endl;
	std::cout << "number of antiS found with the correct trackingparticle granddaughters: " << numberOfAntiSWithCorrectGranddaughters << std::endl;

	std::cout << "weighed number of generated antiS (unique): " << nTotalUniqueGenS_weighted << std::endl;
	std::cout << "weighed number of reconstructed antiS: " << weighedRecoAntiS << std::endl;

	_nGENAntiS.clear();
	_nRECOAntiS.clear();
	_nGENAntiS.push_back(nTotalUniqueGenS_weighted);
	_nRECOAntiS.push_back(weighedRecoAntiS);
	_tree_counter->Fill();

	std::cout << "non weighed number of generated antiS (unique): " << nTotalUniqueGenS_Nonweighted << std::endl;
	std::cout << "non weighed number of reconstructed antiS: " << nonweighedRecoAntiS << std::endl;
}


DEFINE_FWK_MODULE(FlatTreeProducerTracking);
