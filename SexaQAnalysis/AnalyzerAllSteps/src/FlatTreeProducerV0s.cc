#include "../interface/FlatTreeProducerV0s.h"
#include <typeinfo>

typedef math::XYZTLorentzVector LorentzVector;

FlatTreeProducerV0s::FlatTreeProducerV0s(edm::ParameterSet const& pset):
  m_lookAtAntiS(pset.getUntrackedParameter<bool>("lookAtAntiS")),
  m_runningOnData(pset.getUntrackedParameter<bool>("runningOnData")),
  m_bsTag(pset.getParameter<edm::InputTag>("beamspot")),
  m_offlinePVTag(pset.getParameter<edm::InputTag>("offlinePV")),
  m_genParticlesTag_GEN(pset.getParameter<edm::InputTag>("genCollection_GEN")),
  m_genParticlesTag_SIM_GEANT(pset.getParameter<edm::InputTag>("genCollection_SIM_GEANT")),
  //m_generalTracksTag(pset.getParameter<edm::InputTag>("generalTracksCollection")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_V0KsTag(pset.getParameter<edm::InputTag>("V0KsCollection")),
  m_V0LTag(pset.getParameter<edm::InputTag>("V0LCollection")),
  m_muonsTag(pset.getParameter<edm::InputTag>("muonsCollection")),
  m_jetsTag(pset.getParameter<edm::InputTag>("jetsCollection")),


  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_offlinePVToken    (consumes<vector<reco::Vertex>>(m_offlinePVTag)),
  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM_GEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM_GEANT)),
  //m_generalTracksToken(consumes<vector<reco::Track> >(m_generalTracksTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_V0KsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0KsTag)),
  m_V0LToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0LTag)),
  m_muonsToken(consumes<vector<reco::Muon>  >(m_muonsTag)),
  m_jetsToken(consumes<vector<reco::PFJet>  >(m_jetsTag))



{
	HLTTagToken_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
	triggerPrescalesToken_ = consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
}


void FlatTreeProducerV0s::beginJob() {

    
        // Initialize when class is created
        edm::Service<TFileService> fs ;

	//for the GEN Ks: just fill tree with few variables so you have an idea if the GEN Ks are really correctly modeuled
	_tree_GEN_Ks = fs->make <TTree>("FlatTreeGENKs","treeGENKs");
	_tree_GEN_Ks->Branch("_GEN_Ks_mass",&_GEN_Ks_mass);
	_tree_GEN_Ks->Branch("_GEN_Ks_pt",&_GEN_Ks_pt);
        
	//for the Ks
	_tree_Ks = fs->make <TTree>("FlatTreeKs","treeKs");

	_tree_Ks->Branch("_Ks_mass",&_Ks_mass);
	_tree_Ks->Branch("_Ks_pt",&_Ks_pt);
	_tree_Ks->Branch("_Ks_pz",&_Ks_pz);
	_tree_Ks->Branch("_Ks_Lxy",&_Ks_Lxy);
	_tree_Ks->Branch("_Ks_vz",&_Ks_vz);
	_tree_Ks->Branch("_Ks_eta",&_Ks_eta);
	_tree_Ks->Branch("_Ks_phi",&_Ks_phi);
	_tree_Ks->Branch("_Ks_dxy_beamspot",&_Ks_dxy_beamspot);
	_tree_Ks->Branch("_Ks_dxy_min_PV",&_Ks_dxy_min_PV);
	_tree_Ks->Branch("_Ks_dxy_PV0",&_Ks_dxy_PV0);
	_tree_Ks->Branch("_Ks_dxy_000",&_Ks_dxy_000);


	_tree_Ks->Branch("_Ks_dz_beamspot",&_Ks_dz_beamspot);
	_tree_Ks->Branch("_Ks_dz_min_PV",&_Ks_dz_min_PV);
	_tree_Ks->Branch("_Ks_dz_PV0",&_Ks_dz_PV0);
	_tree_Ks->Branch("_Ks_dz_000",&_Ks_dz_000);
	_tree_Ks->Branch("_Ks_vz_dz_min_PV",&_Ks_vz_dz_min_PV);
	_tree_Ks->Branch("_Ks_deltaRBestMatchingGENParticle",&_Ks_deltaRBestMatchingGENParticle);

	_tree_Ks->Branch("_Ks_trackPair_mindeltaR",&_Ks_trackPair_mindeltaR);
	_tree_Ks->Branch("_Ks_trackPair_mass",&_Ks_trackPair_mass);

	_tree_Ks->Branch("_Ks_Track1Track2_openingsAngle",&_Ks_Track1Track2_openingsAngle);
	_tree_Ks->Branch("_Ks_Track1Track2_deltaR",&_Ks_Track1Track2_deltaR);

	_tree_Ks->Branch("_Ks_Track1_openingsAngle",&_Ks_Track1_openingsAngle);
	_tree_Ks->Branch("_Ks_Track2_openingsAngle",&_Ks_Track2_openingsAngle);
	_tree_Ks->Branch("_Ks_Track1_deltaR",&_Ks_Track1_deltaR);
	_tree_Ks->Branch("_Ks_Track2_deltaR",&_Ks_Track2_deltaR);


	_tree_Ks->Branch("_Ks_daughterTrack1_charge",&_Ks_daughterTrack1_charge);
	_tree_Ks->Branch("_Ks_daughterTrack1_chi2",&_Ks_daughterTrack1_chi2);
	_tree_Ks->Branch("_Ks_daughterTrack1_ndof",&_Ks_daughterTrack1_ndof);
	_tree_Ks->Branch("_Ks_daughterTrack1_eta",&_Ks_daughterTrack1_eta);
	_tree_Ks->Branch("_Ks_daughterTrack1_phi",&_Ks_daughterTrack1_phi);
	_tree_Ks->Branch("_Ks_daughterTrack1_pt",&_Ks_daughterTrack1_pt);
	_tree_Ks->Branch("_Ks_daughterTrack1_pz",&_Ks_daughterTrack1_pz);
	_tree_Ks->Branch("_Ks_daughterTrack1_dxy_beamspot",&_Ks_daughterTrack1_dxy_beamspot);
	_tree_Ks->Branch("_Ks_daughterTrack1_dz_beamspot",&_Ks_daughterTrack1_dz_beamspot);
	_tree_Ks->Branch("_Ks_daughterTrack1_dz_min_PV",&_Ks_daughterTrack1_dz_min_PV);
	_tree_Ks->Branch("_Ks_daughterTrack1_dz_PV0",&_Ks_daughterTrack1_dz_PV0);
	_tree_Ks->Branch("_Ks_daughterTrack1_dz_000",&_Ks_daughterTrack1_dz_000);

	_tree_Ks->Branch("_Ks_daughterTrack2_charge",&_Ks_daughterTrack2_charge);
	_tree_Ks->Branch("_Ks_daughterTrack2_chi2",&_Ks_daughterTrack2_chi2);
	_tree_Ks->Branch("_Ks_daughterTrack2_ndof",&_Ks_daughterTrack2_ndof);
	_tree_Ks->Branch("_Ks_daughterTrack2_eta",&_Ks_daughterTrack2_eta);
	_tree_Ks->Branch("_Ks_daughterTrack2_phi",&_Ks_daughterTrack2_phi);
	_tree_Ks->Branch("_Ks_daughterTrack2_pt",&_Ks_daughterTrack2_pt);
	_tree_Ks->Branch("_Ks_daughterTrack2_pz",&_Ks_daughterTrack2_pz);
	_tree_Ks->Branch("_Ks_daughterTrack2_dxy_beamspot",&_Ks_daughterTrack2_dxy_beamspot);
	_tree_Ks->Branch("_Ks_daughterTrack2_dz_beamspot",&_Ks_daughterTrack2_dz_beamspot);
	_tree_Ks->Branch("_Ks_daughterTrack2_dz_min_PV",&_Ks_daughterTrack2_dz_min_PV);
	_tree_Ks->Branch("_Ks_daughterTrack2_dz_PV0",&_Ks_daughterTrack2_dz_PV0);
	_tree_Ks->Branch("_Ks_daughterTrack2_dz_000",&_Ks_daughterTrack2_dz_000);

	//for the Lambda
        _tree_Lambda = fs->make <TTree>("FlatTreeLambda","treeLambda");

	_tree_Lambda->Branch("_Lambda_mass",&_Lambda_mass);
	_tree_Lambda->Branch("_Lambda_pt",&_Lambda_pt);
	_tree_Lambda->Branch("_Lambda_pz",&_Lambda_pz);
	_tree_Lambda->Branch("_Lambda_Lxy",&_Lambda_Lxy);
	_tree_Lambda->Branch("_Lambda_vz",&_Lambda_vz);
	_tree_Lambda->Branch("_Lambda_eta",&_Lambda_eta);
	_tree_Lambda->Branch("_Lambda_phi",&_Lambda_phi);
	_tree_Lambda->Branch("_Lambda_dxy_beamspot",&_Lambda_dxy_beamspot);
	_tree_Lambda->Branch("_Lambda_dxy_min_PV",&_Lambda_dxy_min_PV);
	_tree_Lambda->Branch("_Lambda_dxy_PV0",&_Lambda_dxy_PV0);
	_tree_Lambda->Branch("_Lambda_dxy_000",&_Lambda_dxy_000);


	_tree_Lambda->Branch("_Lambda_dz_beamspot",&_Lambda_dz_beamspot);
	_tree_Lambda->Branch("_Lambda_dz_min_PV",&_Lambda_dz_min_PV);
	_tree_Lambda->Branch("_Lambda_dz_PV0",&_Lambda_dz_PV0);
	_tree_Lambda->Branch("_Lambda_dz_000",&_Lambda_dz_000);
	_tree_Lambda->Branch("_Lambda_vz_dz_min_PV",&_Lambda_vz_dz_min_PV);
	_tree_Lambda->Branch("_Lambda_deltaRBestMatchingGENParticle",&_Lambda_deltaRBestMatchingGENParticle);

	_tree_Lambda->Branch("_Lambda_trackPair_mindeltaR",&_Lambda_trackPair_mindeltaR);
	_tree_Lambda->Branch("_Lambda_trackPair_mass",&_Lambda_trackPair_mass);

	_tree_Lambda->Branch("_Lambda_Track1Track2_openingsAngle",&_Lambda_Track1Track2_openingsAngle);
	_tree_Lambda->Branch("_Lambda_Track1Track2_deltaR",&_Lambda_Track1Track2_deltaR);

	_tree_Lambda->Branch("_Lambda_Track1_openingsAngle",&_Lambda_Track1_openingsAngle);
	_tree_Lambda->Branch("_Lambda_Track2_openingsAngle",&_Lambda_Track2_openingsAngle);
	_tree_Lambda->Branch("_Lambda_Track1_deltaR",&_Lambda_Track1_deltaR);
	_tree_Lambda->Branch("_Lambda_Track2_deltaR",&_Lambda_Track2_deltaR);


	_tree_Lambda->Branch("_Lambda_daughterTrack1_charge",&_Lambda_daughterTrack1_charge);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_chi2",&_Lambda_daughterTrack1_chi2);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_ndof",&_Lambda_daughterTrack1_ndof);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_eta",&_Lambda_daughterTrack1_eta);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_phi",&_Lambda_daughterTrack1_phi);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_pt",&_Lambda_daughterTrack1_pt);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_pz",&_Lambda_daughterTrack1_pz);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_dxy_beamspot",&_Lambda_daughterTrack1_dxy_beamspot);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_dz_beamspot",&_Lambda_daughterTrack1_dz_beamspot);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_dz_min_PV",&_Lambda_daughterTrack1_dz_min_PV);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_dz_PV0",&_Lambda_daughterTrack1_dz_PV0);
	_tree_Lambda->Branch("_Lambda_daughterTrack1_dz_000",&_Lambda_daughterTrack1_dz_000);

	_tree_Lambda->Branch("_Lambda_daughterTrack2_charge",&_Lambda_daughterTrack2_charge);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_chi2",&_Lambda_daughterTrack2_chi2);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_ndof",&_Lambda_daughterTrack2_ndof);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_eta",&_Lambda_daughterTrack2_eta);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_phi",&_Lambda_daughterTrack2_phi);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_pt",&_Lambda_daughterTrack2_pt);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_pz",&_Lambda_daughterTrack2_pz);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_dxy_beamspot",&_Lambda_daughterTrack2_dxy_beamspot);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_dz_beamspot",&_Lambda_daughterTrack2_dz_beamspot);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_dz_min_PV",&_Lambda_daughterTrack2_dz_min_PV);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_dz_PV0",&_Lambda_daughterTrack2_dz_PV0);
	_tree_Lambda->Branch("_Lambda_daughterTrack2_dz_000",&_Lambda_daughterTrack2_dz_000);


	//for the Z
        _tree_Z = fs->make <TTree>("FlatTreeZ","treeZ");

        _tree_Z->Branch("_Z_mass",&_Z_mass);
        _tree_Z->Branch("_Z_dz_PV_muon1",&_Z_dz_PV_muon1);
        _tree_Z->Branch("_Z_dz_PV_muon2",&_Z_dz_PV_muon2);
        _tree_Z->Branch("_Z_ptMuMu",&_Z_ptMuMu);


	//for the PV 
        _tree_PV = fs->make <TTree>("FlatTreePV","treePV");

        _tree_PV->Branch("_PV_n",&_PV_n);
        _tree_PV->Branch("_PV0_lxy",&_PV0_lxy);
        _tree_PV->Branch("_PV0_vz",&_PV0_vz);

        //for the beamspot
        _tree_beamspot = fs->make <TTree>("FlatTreeBeamspot","treeBeamspot");
        _tree_beamspot->Branch("_beampot_lxy",&_beampot_lxy);
        _tree_beamspot->Branch("_beampot_vz",&_beampot_vz);

	//some generalities:
	_tree_general = fs->make <TTree>("FlatTreeGeneral","treeGeneral");
	_tree_general->Branch("_general_triggerFired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&_general_triggerFired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
	_tree_general->Branch("_general_triggerFired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&_general_triggerFired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
	_tree_general->Branch("_general_eventTrackMultiplicity",&_general_eventTrackMultiplicity);
	_tree_general->Branch("_general_eventTrackMultiplicity_highPurity",&_general_eventTrackMultiplicity_highPurity);
	
}

void FlatTreeProducerV0s::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

  //primary vertex
  edm::Handle<vector<reco::Vertex>> h_offlinePV;
  iEvent.getByToken(m_offlinePVToken, h_offlinePV);

  //SIM particles: normal Gen particles or PlusGEANT
  edm::Handle<vector<reco::GenParticle>> h_genParticles;
  iEvent.getByToken(m_genParticlesToken_GEN, h_genParticles);
  //iEvent.getByToken(m_genParticlesToken_SIM_GEANT, h_genParticles);

  //General tracks particles
  //edm::Handle<vector<reco::Track>> h_generalTracks;
  //edm::Handle<View<reco::Track>> h_generalTracks;
  //iEvent.getByToken(m_generalTracksToken, h_generalTracks);

  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);

  //V0 Kshorts
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks;
  iEvent.getByToken(m_V0KsToken, h_V0Ks);

  //V0 Lambdas
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L;
  iEvent.getByToken(m_V0LToken, h_V0L);

  //muons
  edm::Handle<vector<reco::Muon> > h_muons;
  iEvent.getByToken(m_muonsToken, h_muons);

  //jets
  edm::Handle<vector<reco::PFJet> > h_jets;
  iEvent.getByToken(m_jetsToken, h_jets);

  //trigger information
  edm::Handle< pat::PackedTriggerPrescales > triggerPrescales;
  edm::Handle< edm::TriggerResults > HLTResHandle;
  iEvent.getByToken(triggerPrescalesToken_, triggerPrescales);
  iEvent.getByToken(HLTTagToken_, HLTResHandle);  


  //beamspot
  TVector3 beamspot(0,0,0);
  TVector3 beamspotVariance(0,0,0);
  if(h_bs.isValid()){ 
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));			
  }

  if(h_genParticles.isValid()){//loop over the gen particles, find the Ks and save some of the kinematic variables at GEN level
	for(unsigned int i = 0; i < h_genParticles->size(); ++i){
		if(h_genParticles->at(i).pdgId() == AnalyzerAllSteps::pdgIdKs){
			InitGENKs();
			_GEN_Ks_mass.push_back(h_genParticles->at(i).mass());	
			_GEN_Ks_pt.push_back(h_genParticles->at(i).pt());	
			_tree_GEN_Ks->Fill();
		}
	}
  }

  bool Fired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = false;
  bool Fired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = false;
  //std::cout << "-----------------------------------------------" << std::endl;
  if ( HLTResHandle.isValid() && !HLTResHandle.failedToGet() ) {

	edm::RefProd<edm::TriggerNames> trigNames( &(iEvent.triggerNames( *HLTResHandle )) );
	int ntrigs = (int)trigNames->size();
	for (int i = 0; i < ntrigs; i++) {

	      //std::cout << trigNames->triggerName(i) << std::endl;
	      //if(HLTResHandle->accept(i) == 1)std::cout << trigNames->triggerName(i) << " : " << HLTResHandle->accept(i) << std::endl;

	      if(trigNames->triggerName(i).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ") != std::string::npos){
			 Fired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = HLTResHandle->accept(i);			
	      }
	      if(trigNames->triggerName(i).find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ") != std::string::npos){
			 Fired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = HLTResHandle->accept(i);			
	      }
	}

  }
  else std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!HLTResHandle collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;


  //before you start filling the trees with reconstructed Ks and Lambdas in this event you have to check wether this event contains a Z boson candidate.
  bool ZCandidatePresent = false;
  double ZCandidateMass = 0;
  double ZCandidatePhi = 0;
  double dz_PV_muon1 = 999.;
  double dz_PV_muon2 = 999.;
  double pTMuMu = 999.;
  double nLooseMuons = 0; 
//  double ptMuon1 = 0; 
//  double ptMuon2 = 0; 
//  double phiMuon1 = 0; 
//  double phiMuon2 = 0; 
//  double etaMuon1 = 0; 
//  double etaMuon2 = 0; 
  if(h_muons.isValid()){
	for(unsigned int i = 0; i < h_muons->size(); ++i){

		//now check if this muon i is tight and isolated, has a minmal pt and is in the acceptance, if not do not consider this muon

		if(h_muons->at(i).pt() < 20) continue;	
		if(fabs(h_muons->at(i).eta()) > 2.4) continue;	
		if(! IsolationCriterium(h_muons->at(i))) continue;
		//number of loose muons will be used later to veto events with more than 2 loose muons.	
		bool muon1Loose = muon::isLooseMuon(h_muons->at(i));
                if(muon1Loose) nLooseMuons++;

		bool muon1Tight = muon::isTightMuon(h_muons->at(i), h_offlinePV->at(0));
		if(!muon1Tight) continue;

		for(unsigned int j = i + 1; j < h_muons->size(); ++j){
			//check if the muons have opposite charge:
			if(h_muons->at(j).charge() == h_muons->at(i).charge()) continue;
			
			//now check if this muon j is tight and isolated, has a minmal pt and is in the acceptance, if not do not consider this muon
			if(h_muons->at(j).pt() < 20) continue;	
			if(fabs(h_muons->at(j).eta()) > 2.4) continue;	
			if(! IsolationCriterium(h_muons->at(j))) continue;

			bool muon2Tight = muon::isTightMuon(h_muons->at(j), h_offlinePV->at(0));
			if(!muon2Tight) continue;

			LorentzVector p4ZCandidate = h_muons->at(i).p4() + h_muons->at(j).p4();

			//now check if the ONLY jet in the event is going back to back with the Z, jets with momenta lower than 30 GeV can be ignored
			bool ContaminatingHighPtJetFound = false;

			for(unsigned int k = 0; k < h_jets->size(); ++k){

				if(h_jets->at(k).pt() < 30) continue;//should be .energy()  ????
				double deltaPhiJetZ = reco::deltaPhi(p4ZCandidate.phi(),h_jets->at(k).phi());
				//if the jet is outside of the phi-cone with opening pi/4 around the backToBack of the Z then this jet can affect the cleanlyness of the the transverse region
				bool jetInBackToBackRegion = false;
				bool jetInForwardRegion    = false;
				if(abs(reco::deltaPhi(deltaPhiJetZ,TMath::Pi())) < TMath::Pi()/4 ) jetInBackToBackRegion = true;
				if(abs(reco::deltaPhi(deltaPhiJetZ, 0.)) < TMath::Pi()/4 ) jetInForwardRegion = true;
				if(jetInBackToBackRegion || jetInForwardRegion) continue;
				ContaminatingHighPtJetFound = true;

			}

			if(ContaminatingHighPtJetFound) continue;

			//event survived all the cuts :-).
			ZCandidatePresent = true; 
			
			//some variables to be saved later
			ZCandidateMass =  p4ZCandidate.mass();
			ZCandidatePhi = p4ZCandidate.phi();

			TVector3 PV0(h_offlinePV->at(0).x(),h_offlinePV->at(0).y(),h_offlinePV->at(0).z());

			TVector3 muon1Vertex(h_muons->at(i).vx(),h_muons->at(i).vy(),h_muons->at(i).vz());
			TVector3 muon1Momentum(h_muons->at(i).px(),h_muons->at(i).py(),h_muons->at(i).pz());
			//dz_PV_muon1  = AnalyzerAllSteps::dz_line_point(muon1Vertex,muon1Momentum,PV0);
			dz_PV_muon1  =  h_muons->at(i).muonBestTrack()->dz( h_offlinePV->at(0).position()) ;

			TVector3 muon2Vertex(h_muons->at(j).vx(),h_muons->at(j).vy(),h_muons->at(j).vz());
			TVector3 muon2Momentum(h_muons->at(j).px(),h_muons->at(j).py(),h_muons->at(j).pz());
			//dz_PV_muon2  = AnalyzerAllSteps::dz_line_point(muon2Vertex,muon2Momentum,PV0);
			dz_PV_muon2  =  h_muons->at(j).muonBestTrack()->dz( h_offlinePV->at(0).position()) ;

			pTMuMu  = sqrt( pow( h_muons->at(i).px() + h_muons->at(j).px() , 2 ) +  pow( h_muons->at(i).py() + h_muons->at(j).py() , 2 ) );

			//ptMuon1 = h_muons->at(i).pt();		
			//ptMuon2 = h_muons->at(j).pt();		
			//etaMuon1 = h_muons->at(i).eta();		
			//etaMuon2 = h_muons->at(j).eta();		
			//phiMuon1 = h_muons->at(i).phi();		
			//phiMuon2 = h_muons->at(j).phi();		
		}


	}
  }


  //now that you have the phi of the Z candidate you know in which direction the hard event goes, so now look at the V0s in the underlying event, which are the V0s which are not in the deltaPhi cone of the hard event, this cone is defined as 60° around the Z and 60° around the Z in the back to back

  if(ZCandidatePresent  &&  nLooseMuons <= 2 && abs(dz_PV_muon1)  < 0.5 && abs(dz_PV_muon2) < 0.5 && abs(ZCandidateMass - 91.1876) < 15/2){

	std::cout << "ZCandMass = " << ZCandidateMass << " ZCandidatePhi " << ZCandidatePhi  << std::endl;

//	std::cout << "ptMuon1 = " << ptMuon1 << std::endl;
//	std::cout << "ptMuon2 = " << ptMuon2 << std::endl;
//	std::cout << "etaMuon1 = " << etaMuon1 << std::endl;
//	std::cout << "etaMuon2 = " << etaMuon2 << std::endl;
//	std::cout << "phiMuon1 = " << phiMuon1 << std::endl;
//	std::cout << "phiMuon2 = " << phiMuon2 << std::endl;
	
	//save some variables to the Z tree
	InitZ();
	_Z_mass.push_back(ZCandidateMass);
	_Z_dz_PV_muon1.push_back(dz_PV_muon1);
	_Z_dz_PV_muon2.push_back(dz_PV_muon2);
	_Z_ptMuMu.push_back(pTMuMu);
	_tree_Z->Fill();

	InitPV();
	_PV_n.push_back(h_offlinePV->size());
	_PV0_lxy.push_back( sqrt( pow( h_offlinePV->at(0).x()- h_bs->x0() , 2) + pow( h_offlinePV->at(0).y()- h_bs->y0() , 2)  ) );
	_PV0_vz.push_back(h_offlinePV->at(0).z());
	_tree_PV->Fill();

	InitBeamspot();
	_beampot_lxy.push_back( sqrt( pow( h_bs->x0() , 2) + pow( h_bs->y0() , 2) ));
	_beampot_vz.push_back(h_bs->z0());
	_tree_beamspot->Fill();       


	//loop over the track collection and count the number of high purity tracks to have an idea 
	//int nHighPurityTracks = 0;
	//if(h_generalTracks.isValid()){
	//	for(unsigned int i = 0; i < h_generalTracks->size(); ++i){
	//		const reco::Track *ptr = &h_generalTracks->at(i);
	//		if( AnalyzerAllSteps::trackQualityAsInt(ptr) == 2) nHighPurityTracks++;
	//	}
	//}

	InitGeneral();
	_general_triggerFired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ.push_back(Fired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
	_general_triggerFired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ.push_back(Fired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
	//_general_eventTrackMultiplicity.push_back(h_generalTracks->size());
	//_general_eventTrackMultiplicity_highPurity.push_back(nHighPurityTracks);
	_tree_general->Fill();
	

	TVector3 PV0(h_offlinePV->at(0).x(),h_offlinePV->at(0).y(),h_offlinePV->at(0).z());

	InitKs();
	  if(h_V0Ks.isValid()){
	      for(unsigned int i = 0; i < h_V0Ks->size(); ++i){//loop all RECO Ks
		const reco::VertexCompositeCandidate * Ks = &h_V0Ks->at(i);
		double deltaPhiKsHardCone = reco::deltaPhi( Ks->phi() , ZCandidatePhi ); 
		double deltaPhiKsBackToBackHardCone = reco::deltaPhi( Ks->phi() , -TMath::Pi() + ZCandidatePhi ); 
		//only select Ks in the transverse region OR if the Ks is pointing in dz far enough from the hard PV
		
		TVector3 V0CreationVertex(Ks->vx(),Ks->vy(),Ks->vz());
        	TVector3 V0Momentum(Ks->px(),Ks->py(),Ks->pz());
		double dz_PV0 = AnalyzerAllSteps::dz_line_point(V0CreationVertex,V0Momentum,PV0);

		if(  (abs(deltaPhiKsHardCone) >  TMath::Pi()/3 &&  abs(deltaPhiKsBackToBackHardCone) >  TMath::Pi()/3)  || abs(dz_PV0)>1  )  FillBranchesV0(Ks, beamspot, beamspotVariance, h_offlinePV, h_genParticles,  "Ks");
	      }
	  }
	_tree_Ks->Fill();
	InitLambda();
	  if(h_V0L.isValid()){
	      for(unsigned int i = 0; i < h_V0L->size(); ++i){//loop all RECO Lambdas
		const reco::VertexCompositeCandidate * L = &h_V0L->at(i);
		double deltaPhiLHardCone = reco::deltaPhi( L->phi() , ZCandidatePhi );
		double deltaPhiLBackToBackHardCone = reco::deltaPhi( L->phi() , -TMath::Pi() + ZCandidatePhi );

		TVector3 V0CreationVertex(L->vx(),L->vy(),L->vz());
        	TVector3 V0Momentum(L->px(),L->py(),L->pz());
		double dz_PV0 = AnalyzerAllSteps::dz_line_point(V0CreationVertex,V0Momentum,PV0);

		if(  (abs(deltaPhiLHardCone) >  TMath::Pi()/3 &&  abs(deltaPhiLBackToBackHardCone) >  TMath::Pi()/3)  || abs(dz_PV0)>1  )  FillBranchesV0(L, beamspot, beamspotVariance, h_offlinePV,h_genParticles, "Lambda");
	      }
	  }
	_tree_Lambda->Fill();

 }


 } //end of analyzer




void FlatTreeProducerV0s::FillBranchesV0(const reco::VertexCompositeCandidate * RECOV0, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV, edm::Handle<vector<reco::GenParticle>> h_genParticles, std::string V0Type){


	math::XYZPoint beamspotPoint(beamspot.X(),beamspot.Y(),beamspot.Z());
	
	//loop over the GEN particles and try to find a Kshort which mathches this RECO Kshort. Then save the status of this particle so you know from where it comes: PV, material or maybe a fake?
	int bestMatchingGENParticle = -1;
	double deltaRBestMatchingGENParticle = 99;
	if(h_genParticles.isValid() && V0Type == "Ks"){
		for(unsigned int i = 0; i < h_genParticles->size(); ++i){

			if(h_genParticles->at(i).pdgId() == AnalyzerAllSteps::pdgIdKs){
				double deltaPhi = reco::deltaPhi(h_genParticles->at(i).phi(), RECOV0->phi()); 
				double deltaEta = abs(h_genParticles->at(i).eta() - RECOV0->eta()); 
				double deltaR = sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
				if(deltaR < deltaRBestMatchingGENParticle){ 
					bestMatchingGENParticle = i;
					deltaRBestMatchingGENParticle = deltaR;
				}

			}

		}		
	}

	if(h_genParticles.isValid() && V0Type == "Lambda"){
		for(unsigned int i = 0; i < h_genParticles->size(); ++i){

			if(abs(h_genParticles->at(i).pdgId()) == abs(AnalyzerAllSteps::pdgIdAntiLambda) ){
				double deltaPhi = reco::deltaPhi(h_genParticles->at(i).phi(), RECOV0->phi()); 
				double deltaEta = abs(h_genParticles->at(i).eta() - RECOV0->eta()); 
				double deltaR = sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
				if(deltaR < deltaRBestMatchingGENParticle){ 
					bestMatchingGENParticle = i;
					deltaRBestMatchingGENParticle = deltaR;
				}

			}

		}		
	}


        TVector3 V0CreationVertex(RECOV0->vx(),RECOV0->vy(),RECOV0->vz());
        double Lxy = AnalyzerAllSteps::lxy(beamspot,V0CreationVertex);
        TVector3 V0Momentum(RECOV0->px(),RECOV0->py(),RECOV0->pz());

	//different reference points
	TVector3 PVmin = AnalyzerAllSteps::dz_line_point_min(V0CreationVertex,V0Momentum,h_offlinePV);
	TVector3 PV0(h_offlinePV->at(0).x(),h_offlinePV->at(0).y(),h_offlinePV->at(0).z());
	TVector3 ZeroZeroZero(0.,0.,0.);

	//different dxy definitions
	double dxy_beamspot = AnalyzerAllSteps::dxy_signed_line_point(V0CreationVertex,V0Momentum,beamspot);
	double dxy_min_PV = AnalyzerAllSteps::dxy_signed_line_point(V0CreationVertex,V0Momentum,PVmin);
	double dxy_PV0 = AnalyzerAllSteps::dxy_signed_line_point(V0CreationVertex,V0Momentum,PV0);
	double dxy_000 = AnalyzerAllSteps::dxy_signed_line_point(V0CreationVertex,V0Momentum,ZeroZeroZero);
	//different dz definitions
        double dz_beamspot = AnalyzerAllSteps::dz_line_point(V0CreationVertex,V0Momentum,beamspot);
        double dz_min_PV = AnalyzerAllSteps::dz_line_point(V0CreationVertex,V0Momentum,PVmin);
        double dz_PV0 = AnalyzerAllSteps::dz_line_point(V0CreationVertex,V0Momentum,PV0);
        double dz_000 = AnalyzerAllSteps::dz_line_point(V0CreationVertex,V0Momentum,ZeroZeroZero);

	//save things related to the track daughters of the Ks

	//calculate the invariant mass of the trackpair as if the two tracks were pions (which is the case for Ks)
	TLorentzVector p4bestTrack1(RECOV0->daughter(0)->px(),RECOV0->daughter(0)->py(),RECOV0->daughter(0)->pz(), sqrt( pow(RECOV0->daughter(0)->p(),2) + pow(AnalyzerAllSteps::pdgMassChargedPion,2)  ));
	TLorentzVector p4bestTrack2(RECOV0->daughter(1)->px(),RECOV0->daughter(1)->py(),RECOV0->daughter(1)->pz(), sqrt( pow(RECOV0->daughter(1)->p(),2) + pow(AnalyzerAllSteps::pdgMassChargedPion,2)  ));
	TLorentzVector p4bestTrackPair = p4bestTrack1 + p4bestTrack2;	

	//calculate some angles between the tracks and between the Ks and the tracks, the RECOV0->daughter(0) and RECOV0->daughter(1) objects should have the angles calculated at the vertex of the two tracks, which is a good thing as this is where you want to evaluate the angles. See: https://github.com/jarnedc/cmssw/blob/from-CMSSW_8_0_30/RecoVertex/V0Producer/src/V0Fitter.cc#L435-L440
	//between the two tracks
	double openingsAngleTrack1Track2 = AnalyzerAllSteps::openings_angle(RECOV0->daughter(0)->momentum(),RECOV0->daughter(1)->momentum());
	double deltaRTrack1Track2 = AnalyzerAllSteps::deltaR(RECOV0->daughter(0)->phi(),RECOV0->daughter(0)->eta(),RECOV0->daughter(1)->phi(),RECOV0->daughter(1)->eta());

	//between the tracks and the Ks
	double openingsAngleTrack1V0 = AnalyzerAllSteps::openings_angle(RECOV0->daughter(0)->momentum(),RECOV0->momentum());
	double openingsAngleTrack2V0 = AnalyzerAllSteps::openings_angle(RECOV0->daughter(1)->momentum(),RECOV0->momentum());

	double deltaRTrack1V0 = AnalyzerAllSteps::deltaR(RECOV0->daughter(0)->phi(),RECOV0->daughter(0)->eta(),RECOV0->phi(),RECOV0->eta());
	double deltaRTrack2V0 = AnalyzerAllSteps::deltaR(RECOV0->daughter(1)->phi(),RECOV0->daughter(1)->eta(),RECOV0->phi(),RECOV0->eta());

	TVector3 V0Track1creationVertex(RECOV0->daughter(0)->vx(),RECOV0->daughter(0)->vy(),RECOV0->daughter(0)->vz());
        TVector3 V0Track1Momentum(RECOV0->daughter(0)->px(),RECOV0->daughter(0)->py(),RECOV0->daughter(0)->pz());
	double dxy_Track1V0_beamspot = AnalyzerAllSteps::dxy_signed_line_point(V0Track1creationVertex,V0Track1Momentum,beamspot);	
        double dz_Track1V0_beamspot = AnalyzerAllSteps::dz_line_point(V0Track1creationVertex,V0Track1Momentum,beamspot);
        double dz_Track1V0_min_PV = AnalyzerAllSteps::dz_line_point(V0Track1creationVertex,V0Track1Momentum,PVmin);
        double dz_Track1V0_PV0 = AnalyzerAllSteps::dz_line_point(V0Track1creationVertex,V0Track1Momentum,PV0);
        double dz_Track1V0_000 = AnalyzerAllSteps::dz_line_point(V0Track1creationVertex,V0Track1Momentum,ZeroZeroZero);

	TVector3 V0Track2creationVertex(RECOV0->daughter(1)->vx(),RECOV0->daughter(1)->vy(),RECOV0->daughter(1)->vz());
        TVector3 V0Track2Momentum(RECOV0->daughter(1)->px(),RECOV0->daughter(1)->py(),RECOV0->daughter(1)->pz());
	double dxy_Track2V0_beamspot = AnalyzerAllSteps::dxy_signed_line_point(V0Track2creationVertex,V0Track2Momentum,beamspot);	
        double dz_Track2V0_beamspot = AnalyzerAllSteps::dz_line_point(V0Track2creationVertex,V0Track2Momentum,beamspot);
        double dz_Track2V0_min_PV = AnalyzerAllSteps::dz_line_point(V0Track2creationVertex,V0Track2Momentum,PVmin);
        double dz_Track2V0_PV0 = AnalyzerAllSteps::dz_line_point(V0Track2creationVertex,V0Track2Momentum,PV0);
        double dz_Track2V0_000 = AnalyzerAllSteps::dz_line_point(V0Track2creationVertex,V0Track2Momentum,ZeroZeroZero);


	if(V0Type == "Ks"){
		_Ks_mass.push_back(RECOV0->mass());	

		_Ks_pt.push_back(RECOV0->pt());	
		_Ks_pz.push_back(RECOV0->pz());	
		_Ks_Lxy.push_back(Lxy);	
		_Ks_vz.push_back(RECOV0->vz());	

		_Ks_eta.push_back(RECOV0->eta());
		_Ks_phi.push_back(RECOV0->phi());

		_Ks_dxy_beamspot.push_back(dxy_beamspot);
		_Ks_dxy_min_PV.push_back(dxy_min_PV);
		_Ks_dxy_PV0.push_back(dxy_PV0);
		_Ks_dxy_000.push_back(dxy_000);

		_Ks_dz_beamspot.push_back(dz_beamspot);	
		_Ks_dz_min_PV.push_back(dz_min_PV);	
		_Ks_dz_PV0.push_back(dz_PV0);	
		_Ks_dz_000.push_back(dz_000);

		_Ks_vz_dz_min_PV.push_back(PVmin.Z());

		_Ks_deltaRBestMatchingGENParticle.push_back(deltaRBestMatchingGENParticle);


		//_Ks_trackPair_mindeltaR.push_back(deltaRRecoKsTrackPair);
		_Ks_trackPair_mass.push_back(p4bestTrackPair.M());

		//between the two tracks
		_Ks_Track1Track2_openingsAngle.push_back(openingsAngleTrack1Track2);
		_Ks_Track1Track2_deltaR.push_back(deltaRTrack1Track2);

		//between the two tracks and the Ks
		_Ks_Track1_openingsAngle.push_back(openingsAngleTrack1V0);
		_Ks_Track2_openingsAngle.push_back(openingsAngleTrack2V0);
		_Ks_Track1_deltaR.push_back(deltaRTrack1V0);
		_Ks_Track2_deltaR.push_back(deltaRTrack2V0);
		
		//somehow the RECOV0->daughter(0) should be containing the reference to the track (see e.g. https://github.com/jarnedc/cmssw/blob/from-CMSSW_8_0_30/RecoVertex/V0Producer/src/V0Fitter.cc#L515), but I am not able to get it out...
		//const reco::Candidate* Ks_daug0 = RECOV0->daughter(0);
		//reco::TrackRef Ks_daug0_trackref = Ks_daug0->track();

		
		_Ks_daughterTrack1_charge.push_back(RECOV0->daughter(0)->charge());
		//_Ks_daughterTrack1_chi2.push_back(RECOV0->daughter(0)->chi2());
		//_Ks_daughterTrack1_ndof.push_back(RECOV0->daughter(0)->ndof());
		_Ks_daughterTrack1_eta.push_back(RECOV0->daughter(0)->eta());
		_Ks_daughterTrack1_phi.push_back(RECOV0->daughter(0)->phi());
		_Ks_daughterTrack1_pt.push_back(RECOV0->daughter(0)->pt());
		_Ks_daughterTrack1_pz.push_back(RECOV0->daughter(0)->pz());
		_Ks_daughterTrack1_dxy_beamspot.push_back(dxy_Track1V0_beamspot);
		_Ks_daughterTrack1_dz_beamspot.push_back(dz_Track1V0_beamspot);
		_Ks_daughterTrack1_dz_min_PV.push_back(dz_Track1V0_min_PV);
		_Ks_daughterTrack1_dz_PV0.push_back(dz_Track1V0_PV0);
		_Ks_daughterTrack1_dz_000.push_back(dz_Track1V0_000);

		_Ks_daughterTrack2_charge.push_back(RECOV0->daughter(1)->charge());
		//_Ks_daughterTrack2_chi2.push_back(RECOV0->daughter(1)->chi2());
		//_Ks_daughterTrack2_ndof.push_back(RECOV0->daughter(1)->ndof());
		_Ks_daughterTrack2_eta.push_back(RECOV0->daughter(1)->eta());
		_Ks_daughterTrack2_phi.push_back(RECOV0->daughter(1)->phi());
		_Ks_daughterTrack2_pt.push_back(RECOV0->daughter(1)->pt());
		_Ks_daughterTrack2_pz.push_back(RECOV0->daughter(1)->pz());
		_Ks_daughterTrack2_dxy_beamspot.push_back(dxy_Track2V0_beamspot);
		_Ks_daughterTrack2_dz_beamspot.push_back(dz_Track2V0_beamspot);
		_Ks_daughterTrack2_dz_min_PV.push_back(dz_Track2V0_min_PV);
		_Ks_daughterTrack2_dz_PV0.push_back(dz_Track2V0_PV0);
		_Ks_daughterTrack2_dz_000.push_back(dz_Track2V0_000);
	}
	else if (V0Type == "Lambda"){
		_Lambda_mass.push_back(RECOV0->mass());	

		_Lambda_pt.push_back(RECOV0->pt());	
		_Lambda_pz.push_back(RECOV0->pz());	
		_Lambda_Lxy.push_back(Lxy);	
		_Lambda_vz.push_back(RECOV0->vz());	

		_Lambda_eta.push_back(RECOV0->eta());
		_Lambda_phi.push_back(RECOV0->phi());

		_Lambda_dxy_beamspot.push_back(dxy_beamspot);
		_Lambda_dxy_min_PV.push_back(dxy_min_PV);
		_Lambda_dxy_PV0.push_back(dxy_PV0);
		_Lambda_dxy_000.push_back(dxy_000);

		_Lambda_dz_beamspot.push_back(dz_beamspot);	
		_Lambda_dz_min_PV.push_back(dz_min_PV);	
		_Lambda_dz_PV0.push_back(dz_PV0);	
		_Lambda_dz_000.push_back(dz_000);

		_Lambda_vz_dz_min_PV.push_back(PVmin.Z());

		_Lambda_deltaRBestMatchingGENParticle.push_back(deltaRBestMatchingGENParticle);


		//_Ks_trackPair_mindeltaR.push_back(deltaRRecoKsTrackPair);
		_Lambda_trackPair_mass.push_back(p4bestTrackPair.M());

		//between the two tracks
		_Lambda_Track1Track2_openingsAngle.push_back(openingsAngleTrack1Track2);
		_Lambda_Track1Track2_deltaR.push_back(deltaRTrack1Track2);

		//between the two tracks and the Ks
		_Lambda_Track1_openingsAngle.push_back(openingsAngleTrack1V0);
		_Lambda_Track2_openingsAngle.push_back(openingsAngleTrack2V0);
		_Lambda_Track1_deltaR.push_back(deltaRTrack1V0);
		_Lambda_Track2_deltaR.push_back(deltaRTrack2V0);
		
		//somehow the RECOV0->daughter(0) should be containing the reference to the track (see e.g. https://github.com/jarnedc/cmssw/blob/from-CMSSW_8_0_30/RecoVertex/V0Producer/src/V0Fitter.cc#L515), but I am not able to get it out...
		//const reco::Candidate* Ks_daug0 = RECOV0->daughter(0);
		//reco::TrackRef Ks_daug0_trackref = Ks_daug0->track();

		
		_Lambda_daughterTrack1_charge.push_back(RECOV0->daughter(0)->charge());
		//_Ks_daughterTrack1_chi2.push_back(RECOV0->daughter(0)->chi2());
		//_Ks_daughterTrack1_ndof.push_back(RECOV0->daughter(0)->ndof());
		_Lambda_daughterTrack1_eta.push_back(RECOV0->daughter(0)->eta());
		_Lambda_daughterTrack1_phi.push_back(RECOV0->daughter(0)->phi());
		_Lambda_daughterTrack1_pt.push_back(RECOV0->daughter(0)->pt());
		_Lambda_daughterTrack1_pz.push_back(RECOV0->daughter(0)->pz());
		_Lambda_daughterTrack1_dxy_beamspot.push_back(dxy_Track1V0_beamspot);
		_Lambda_daughterTrack1_dz_beamspot.push_back(dz_Track1V0_beamspot);
		_Lambda_daughterTrack1_dz_min_PV.push_back(dz_Track1V0_min_PV);
		_Lambda_daughterTrack1_dz_PV0.push_back(dz_Track1V0_PV0);
		_Lambda_daughterTrack1_dz_000.push_back(dz_Track1V0_000);

		_Lambda_daughterTrack2_charge.push_back(RECOV0->daughter(1)->charge());
		//_Ks_daughterTrack2_chi2.push_back(RECOV0->daughter(1)->chi2());
		//_Ks_daughterTrack2_ndof.push_back(RECOV0->daughter(1)->ndof());
		_Lambda_daughterTrack2_eta.push_back(RECOV0->daughter(1)->eta());
		_Lambda_daughterTrack2_phi.push_back(RECOV0->daughter(1)->phi());
		_Lambda_daughterTrack2_pt.push_back(RECOV0->daughter(1)->pt());
		_Lambda_daughterTrack2_pz.push_back(RECOV0->daughter(1)->pz());
		_Lambda_daughterTrack2_dxy_beamspot.push_back(dxy_Track2V0_beamspot);
		_Lambda_daughterTrack2_dz_beamspot.push_back(dz_Track2V0_beamspot);
		_Lambda_daughterTrack2_dz_min_PV.push_back(dz_Track2V0_min_PV);
		_Lambda_daughterTrack2_dz_PV0.push_back(dz_Track2V0_PV0);
		_Lambda_daughterTrack2_dz_000.push_back(dz_Track2V0_000);


	}

	//just a consitency check: if deltaR is small then also the 3D distance between the GEN decay vertex and the RECO decay vertex should be small
	if(h_genParticles.isValid()){
		if( RECOV0->mass() > 0.48 && RECOV0->mass() < 0.52 &&  abs(RECOV0->eta()) < 2 && dxy_beamspot < 0.1 && dxy_beamspot > 0. && abs(dz_PV0) < 0.2 ){
			std::cout << "number of daughters: "<< h_genParticles->at(bestMatchingGENParticle).numberOfDaughters() << std::endl;
			if(deltaRBestMatchingGENParticle < 0.01 && h_genParticles->at(bestMatchingGENParticle).numberOfDaughters() > 0){
				double distance = sqrt(  pow( h_genParticles->at(bestMatchingGENParticle).daughter(0)->vx() - RECOV0->vx() ,2) + pow( h_genParticles->at(bestMatchingGENParticle).daughter(0)->vy() - RECOV0->vy(),2) + pow( h_genParticles->at(bestMatchingGENParticle).daughter(0)->vz() - RECOV0->vz(),2) );
				double GEN_lxyz = sqrt( pow(h_genParticles->at(bestMatchingGENParticle).daughter(0)->vx(),2) + pow(h_genParticles->at(bestMatchingGENParticle).daughter(0)->vy(),2) + pow(h_genParticles->at(bestMatchingGENParticle).daughter(0)->vz(),2) );
				std::cout << "Found a GEN Ks matching a RECO Ks in deltaR, the 3D distance between GEN and RECO decay vertex is: "<< distance << " the GEN lxyz is: " << GEN_lxyz << std::endl;
			}
		}
	}	

}


bool FlatTreeProducerV0s::IsolationCriterium(reco::Muon muon)
{
      bool passedIso = false;

      double chargedHadronIso = muon.pfIsolationR04().sumChargedHadronPt;
      double chargedHadronIsoPU = muon.pfIsolationR04().sumPUPt;
      double neutralHadronIso  = muon.pfIsolationR04().sumNeutralHadronEt;
      double photonIso  = muon.pfIsolationR04().sumPhotonEt;
      double a=0.5;
      // OPTION 1: DeltaBeta corrections for iosolation
      float RelativeIsolationDBetaCorr = (chargedHadronIso + std::max(photonIso+neutralHadronIso - 0.5*chargedHadronIsoPU,0.))/std::max(a, muon.pt());
      if(RelativeIsolationDBetaCorr < 0.15) passedIso = true;

      return passedIso;

}

void FlatTreeProducerV0s::endJob()
{
}

void
FlatTreeProducerV0s::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void
FlatTreeProducerV0s::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
FlatTreeProducerV0s::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
FlatTreeProducerV0s::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatTreeProducerV0s::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

FlatTreeProducerV0s::~FlatTreeProducerV0s()
{
}


void FlatTreeProducerV0s::InitKs()
{

	_Ks_mass.clear();

	_Ks_pt.clear();
	_Ks_pz.clear();
	_Ks_Lxy.clear();	
	_Ks_vz.clear();	

	_Ks_eta.clear();
	_Ks_phi.clear();

	_Ks_dxy_beamspot.clear();
	_Ks_dxy_min_PV.clear();
	_Ks_dxy_PV0.clear();
	_Ks_dxy_000.clear();

	_Ks_dz_beamspot.clear();
	_Ks_dz_min_PV.clear();	
	_Ks_dz_PV0.clear();	
	_Ks_dz_000.clear();	

	_Ks_vz_dz_min_PV.clear();

	_Ks_deltaRBestMatchingGENParticle.clear();	

	_Ks_trackPair_mindeltaR.clear();	
	_Ks_trackPair_mass.clear();	

	_Ks_Track1Track2_openingsAngle.clear();
	_Ks_Track1Track2_deltaR.clear();

	_Ks_Track1_openingsAngle.clear();
	_Ks_Track2_openingsAngle.clear();
	_Ks_Track1_deltaR.clear();
	_Ks_Track2_deltaR.clear();



	_Ks_daughterTrack1_charge.clear();
	_Ks_daughterTrack1_chi2.clear();
	_Ks_daughterTrack1_ndof.clear();
	_Ks_daughterTrack1_eta.clear();	
	_Ks_daughterTrack1_phi.clear();	
	_Ks_daughterTrack1_pt.clear();	
	_Ks_daughterTrack1_pz.clear();	
        _Ks_daughterTrack1_dxy_beamspot.clear();
        _Ks_daughterTrack1_dz_beamspot.clear();
        _Ks_daughterTrack1_dz_min_PV.clear();
        _Ks_daughterTrack1_dz_PV0.clear();
        _Ks_daughterTrack1_dz_000.clear();



	_Ks_daughterTrack2_charge.clear();
	_Ks_daughterTrack2_chi2.clear();
	_Ks_daughterTrack2_ndof.clear();
	_Ks_daughterTrack2_eta.clear();	
	_Ks_daughterTrack2_phi.clear();	
	_Ks_daughterTrack2_pt.clear();	
	_Ks_daughterTrack2_pz.clear();	
        _Ks_daughterTrack2_dxy_beamspot.clear();
        _Ks_daughterTrack2_dz_beamspot.clear();
        _Ks_daughterTrack2_dz_min_PV.clear();
        _Ks_daughterTrack2_dz_PV0.clear();
        _Ks_daughterTrack2_dz_000.clear();
}

void FlatTreeProducerV0s::InitLambda()
{

	_Lambda_mass.clear();

	_Lambda_pt.clear();
	_Lambda_pz.clear();
	_Lambda_Lxy.clear();	
	_Lambda_vz.clear();	

	_Lambda_eta.clear();
	_Lambda_phi.clear();

	_Lambda_dxy_beamspot.clear();
	_Lambda_dxy_min_PV.clear();
	_Lambda_dxy_PV0.clear();
	_Lambda_dxy_000.clear();

	_Lambda_dz_beamspot.clear();
	_Lambda_dz_min_PV.clear();	
	_Lambda_dz_PV0.clear();	
	_Lambda_dz_000.clear();	

	_Lambda_vz_dz_min_PV.clear();

	_Lambda_deltaRBestMatchingGENParticle.clear();	

	_Lambda_trackPair_mindeltaR.clear();	
	_Lambda_trackPair_mass.clear();	

	_Lambda_Track1Track2_openingsAngle.clear();
	_Lambda_Track1Track2_deltaR.clear();

	_Lambda_Track1_openingsAngle.clear();
	_Lambda_Track2_openingsAngle.clear();
	_Lambda_Track1_deltaR.clear();
	_Lambda_Track2_deltaR.clear();



	_Lambda_daughterTrack1_charge.clear();
	_Lambda_daughterTrack1_chi2.clear();
	_Lambda_daughterTrack1_ndof.clear();
	_Lambda_daughterTrack1_eta.clear();	
	_Lambda_daughterTrack1_phi.clear();	
	_Lambda_daughterTrack1_pt.clear();	
	_Lambda_daughterTrack1_pz.clear();	
        _Lambda_daughterTrack1_dxy_beamspot.clear();
        _Lambda_daughterTrack1_dz_beamspot.clear();
        _Lambda_daughterTrack1_dz_min_PV.clear();
        _Lambda_daughterTrack1_dz_PV0.clear();
        _Lambda_daughterTrack1_dz_000.clear();



	_Lambda_daughterTrack2_charge.clear();
	_Lambda_daughterTrack2_chi2.clear();
	_Lambda_daughterTrack2_ndof.clear();
	_Lambda_daughterTrack2_eta.clear();	
	_Lambda_daughterTrack2_phi.clear();	
	_Lambda_daughterTrack2_pt.clear();	
	_Lambda_daughterTrack2_pz.clear();	
        _Lambda_daughterTrack2_dxy_beamspot.clear();
        _Lambda_daughterTrack2_dz_beamspot.clear();
        _Lambda_daughterTrack2_dz_min_PV.clear();
        _Lambda_daughterTrack2_dz_PV0.clear();
        _Lambda_daughterTrack2_dz_000.clear();


}

void FlatTreeProducerV0s::InitGENKs()
{
	_GEN_Ks_mass.clear();
	_GEN_Ks_pt.clear();
}

void FlatTreeProducerV0s::InitZ()
{

	_Z_mass.clear();
	_Z_dz_PV_muon1.clear();
	_Z_dz_PV_muon2.clear();
	_Z_ptMuMu.clear();

}

void FlatTreeProducerV0s::InitPV()
{
	_PV_n.clear();
	_PV0_lxy.clear();
	_PV0_vz.clear();

}

void FlatTreeProducerV0s::InitBeamspot(){
	_beampot_lxy.clear();
	_beampot_vz.clear();
}

void FlatTreeProducerV0s::InitGeneral(){
	_general_triggerFired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ.clear();
	_general_triggerFired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ.clear();

	_general_eventTrackMultiplicity.clear();
	_general_eventTrackMultiplicity_highPurity.clear();
}

DEFINE_FWK_MODULE(FlatTreeProducerV0s);
