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
  m_generalTracksTag(pset.getParameter<edm::InputTag>("generalTracksCollection")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_V0KsTag(pset.getParameter<edm::InputTag>("V0KsCollection")),
  m_V0LTag(pset.getParameter<edm::InputTag>("V0LCollection")),
  m_muonsTag(pset.getParameter<edm::InputTag>("muonsCollection")),
  m_jetsTag(pset.getParameter<edm::InputTag>("jetsCollection")),


  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_offlinePVToken    (consumes<vector<reco::Vertex>>(m_offlinePVTag)),
  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM_GEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM_GEANT)),
  m_generalTracksToken(consumes<View<reco::Track> >(m_generalTracksTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_V0KsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0KsTag)),
  m_V0LToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0LTag)),
  m_muonsToken(consumes<vector<reco::Muon>  >(m_muonsTag)),
  m_jetsToken(consumes<vector<reco::PFJet>  >(m_jetsTag))



{

}


void FlatTreeProducerV0s::beginJob() {

    
        // Initialize when class is created
        edm::Service<TFileService> fs ;
        
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
	_tree_Lambda->Branch("_Lambda_dxy",&_Lambda_dxy);
	_tree_Lambda->Branch("_Lambda_dxy_min_PV",&_Lambda_dxy_min_PV);
	_tree_Lambda->Branch("_Lambda_dxy_PV0",&_Lambda_dxy_PV0);
	_tree_Lambda->Branch("_Lambda_dxy_000",&_Lambda_dxy_000);


	_tree_Lambda->Branch("_Lambda_dz_beamspot",&_Lambda_dz_beamspot);
	_tree_Lambda->Branch("_Lambda_dz_min_PV",&_Lambda_dz_min_PV);
	_tree_Lambda->Branch("_Lambda_dz_PV0",&_Lambda_dz_PV0);
	_tree_Lambda->Branch("_Lambda_dz_000",&_Lambda_dz_000);
	_tree_Lambda->Branch("_Lambda_vz_dz_min_PV",&_Lambda_vz_dz_min_PV);
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

  //muons
  edm::Handle<vector<reco::Muon> > h_muons;
  iEvent.getByToken(m_muonsToken, h_muons);

  //jets
  edm::Handle<vector<reco::PFJet> > h_jets;
  iEvent.getByToken(m_jetsToken, h_jets);

  //beamspot
  TVector3 beamspot(0,0,0);
  TVector3 beamspotVariance(0,0,0);
  if(h_bs.isValid()){ 
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));			
  }

  //before you start filling the trees with reconstructed Ks and Lambdas in this event you have to check wether this event contains a Z boson candidate.
  bool ZCandidatePresent = false;
  double ZCandidateMass = 0;
  double ZCandidatePhi = 0;
  double dz_PV_muon1 = 999.;
  double dz_PV_muon2 = 999.;
  double pTMuMu = 999.;
  double nLooseMuons = 0; 
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

				if(h_jets->at(k).pt() < 30) continue;
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
			dz_PV_muon1  = AnalyzerAllSteps::dz_line_point(muon1Vertex,muon1Momentum,PV0);

			TVector3 muon2Vertex(h_muons->at(j).vx(),h_muons->at(j).vy(),h_muons->at(j).vz());
			TVector3 muon2Momentum(h_muons->at(j).px(),h_muons->at(j).py(),h_muons->at(j).pz());
			dz_PV_muon2  = AnalyzerAllSteps::dz_line_point(muon2Vertex,muon2Momentum,PV0);

			pTMuMu  = sqrt( pow( h_muons->at(i).px() + h_muons->at(j).px() , 2 ) +  pow( h_muons->at(i).py() + h_muons->at(j).py() , 2 ) );
		}


	}
  }


  //now that you have the phi of the Z candidate you know in which direction the hard event goes, so now look at the V0s in the underlying event, which are the V0s which are not in the deltaPhi cone of the hard event, this cone is defined as 60° around the Z and 60° around the Z in the back to back

  if(ZCandidatePresent  &&  nLooseMuons <= 2 && abs(dz_PV_muon1)  < 0.01 && abs(dz_PV_muon2) < 0.01 && abs(ZCandidateMass - 91.1876) < 15/2){

	std::cout << "ZCandMass = " << ZCandidateMass << " ZCandidatePhi " << ZCandidatePhi  << std::endl;
	
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
 

	InitKs();
	  if(h_V0Ks.isValid()){
	      for(unsigned int i = 0; i < h_V0Ks->size(); ++i){//loop all RECO Ks
		const reco::VertexCompositeCandidate * Ks = &h_V0Ks->at(i);
		double deltaPhiKsHardCone = reco::deltaPhi( Ks->phi() , ZCandidatePhi ); 
		double deltaPhiKsBackToBackHardCone = reco::deltaPhi( Ks->phi() , -TMath::Pi() + ZCandidatePhi ); 
		if(abs(deltaPhiKsHardCone) >  TMath::Pi()/3 &&  abs(deltaPhiKsBackToBackHardCone) >  TMath::Pi()/3)FillBranchesKs(Ks, beamspot, beamspotVariance, h_offlinePV, h_genParticles,  h_generalTracks);
	      }
	  }
	_tree_Ks->Fill();

	InitLambda();
	  if(h_V0L.isValid()){
	      for(unsigned int i = 0; i < h_V0L->size(); ++i){//loop all RECO Lambdas
		const reco::VertexCompositeCandidate * L = &h_V0L->at(i);
		double deltaPhiLHardCone = reco::deltaPhi( L->phi() , ZCandidatePhi );
		double deltaPhiLBackToBackHardCone = reco::deltaPhi( L->phi() , -TMath::Pi() + ZCandidatePhi );
		if(abs(deltaPhiLHardCone) >  TMath::Pi()/3 &&  abs(deltaPhiLBackToBackHardCone) >  TMath::Pi()/3)FillBranchesLambda(L, beamspot, beamspotVariance, h_offlinePV,h_genParticles);
	      }
	  }
	_tree_Lambda->Fill();


 }


 } //end of analyzer




void FlatTreeProducerV0s::FillBranchesKs(const reco::VertexCompositeCandidate * RECOKs, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV, edm::Handle<vector<reco::GenParticle>> h_genParticles, edm::Handle<View<reco::Track>> h_generalTracks){


	math::XYZPoint beamspotPoint(beamspot.X(),beamspot.Y(),beamspot.Z());
	
	//loop over the GEN particles and try to find a Kshort which mathches this RECO Kshort. Then save the status of this particle so you know from where it comes: PV, material or maybe a fake?
	//int bestMatchingGENParticle = -1;
	double deltaRBestMatchingGENParticle = 10;
	if(h_genParticles.isValid()){
		for(unsigned int i = 0; i < h_genParticles->size(); ++i){

			if(h_genParticles->at(i).pdgId() == AnalyzerAllSteps::pdgIdKs){
				double deltaPhi = reco::deltaPhi(h_genParticles->at(i).phi(), RECOKs->phi()); 
				double deltaEta = abs(h_genParticles->at(i).eta() - RECOKs->eta()); 
				double deltaR = sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
				if(deltaR < deltaRBestMatchingGENParticle){ 
					//bestMatchingGENParticle = i;
					deltaRBestMatchingGENParticle = deltaR;
				}

			}

		}		
	}

        TVector3 KsCreationVertex(RECOKs->vx(),RECOKs->vy(),RECOKs->vz());
        double Lxy = AnalyzerAllSteps::lxy(beamspot,KsCreationVertex);
        TVector3 KsMomentum(RECOKs->px(),RECOKs->py(),RECOKs->pz());

	//different reference points
	TVector3 PVmin = AnalyzerAllSteps::dz_line_point_min(KsCreationVertex,KsMomentum,h_offlinePV);
	TVector3 PV0(h_offlinePV->at(0).x(),h_offlinePV->at(0).y(),h_offlinePV->at(0).z());
	TVector3 ZeroZeroZero(0.,0.,0.);

	//different dxy definitions
	double dxy_beamspot = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,beamspot);
	double dxy_min_PV = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,PVmin);
	double dxy_PV0 = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,PV0);
	double dxy_000 = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,ZeroZeroZero);
	//different dz definitions
        double dz_beamspot = AnalyzerAllSteps::dz_line_point(KsCreationVertex,KsMomentum,beamspot);
        double dz_min_PV = AnalyzerAllSteps::dz_line_point(KsCreationVertex,KsMomentum,PVmin);
        double dz_PV0 = AnalyzerAllSteps::dz_line_point(KsCreationVertex,KsMomentum,PV0);
        double dz_000 = AnalyzerAllSteps::dz_line_point(KsCreationVertex,KsMomentum,ZeroZeroZero);



	_Ks_mass.push_back(RECOKs->mass());	

	_Ks_pt.push_back(RECOKs->pt());	
	_Ks_pz.push_back(RECOKs->pz());	
	_Ks_Lxy.push_back(Lxy);	
	_Ks_vz.push_back(RECOKs->vz());	

	_Ks_eta.push_back(RECOKs->eta());
	_Ks_phi.push_back(RECOKs->phi());

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

	//save things related to the track daughters of the Ks

	//calculate the invariant mass of the trackpair as if the two tracks were pions (which is the case for Ks)
	TLorentzVector p4bestTrack1(RECOKs->daughter(0)->px(),RECOKs->daughter(0)->py(),RECOKs->daughter(0)->pz(), sqrt( pow(RECOKs->daughter(0)->p(),2) + pow(AnalyzerAllSteps::pdgMassChargedPion,2)  ));
	TLorentzVector p4bestTrack2(RECOKs->daughter(1)->px(),RECOKs->daughter(1)->py(),RECOKs->daughter(1)->pz(), sqrt( pow(RECOKs->daughter(1)->p(),2) + pow(AnalyzerAllSteps::pdgMassChargedPion,2)  ));
	TLorentzVector p4bestTrackPair = p4bestTrack1 + p4bestTrack2;	

	//calculate some angles between the tracks and between the Ks and the tracks, the RECOKs->daughter(0) and RECOKs->daughter(1) objects should have the angles calculated at the vertex of the two tracks, which is a good thing as this is where you want to evaluate the angles. See: https://github.com/jarnedc/cmssw/blob/from-CMSSW_8_0_30/RecoVertex/V0Producer/src/V0Fitter.cc#L435-L440
	//between the two tracks
	double openingsAngleTrack1Track2 = AnalyzerAllSteps::openings_angle(RECOKs->daughter(0)->momentum(),RECOKs->daughter(1)->momentum());
	double deltaRTrack1Track2 = AnalyzerAllSteps::deltaR(RECOKs->daughter(0)->phi(),RECOKs->daughter(0)->eta(),RECOKs->daughter(1)->phi(),RECOKs->daughter(1)->eta());

	//between the tracks and the Ks
	double openingsAngleTrack1Ks = AnalyzerAllSteps::openings_angle(RECOKs->daughter(0)->momentum(),RECOKs->momentum());
	double openingsAngleTrack2Ks = AnalyzerAllSteps::openings_angle(RECOKs->daughter(1)->momentum(),RECOKs->momentum());

	double deltaRTrack1Ks = AnalyzerAllSteps::deltaR(RECOKs->daughter(0)->phi(),RECOKs->daughter(0)->eta(),RECOKs->phi(),RECOKs->eta());
	double deltaRTrack2Ks = AnalyzerAllSteps::deltaR(RECOKs->daughter(1)->phi(),RECOKs->daughter(1)->eta(),RECOKs->phi(),RECOKs->eta());

	TVector3 KsTrack1creationVertex(RECOKs->daughter(0)->vx(),RECOKs->daughter(0)->vy(),RECOKs->daughter(0)->vz());
        TVector3 KsTrack1Momentum(RECOKs->daughter(0)->px(),RECOKs->daughter(0)->py(),RECOKs->daughter(0)->pz());
	double dxy_Track1Ks_beamspot = AnalyzerAllSteps::dxy_signed_line_point(KsTrack1creationVertex,KsTrack1Momentum,beamspot);	
        double dz_Track1Ks_beamspot = AnalyzerAllSteps::dz_line_point(KsTrack1creationVertex,KsTrack1Momentum,beamspot);
        double dz_Track1Ks_min_PV = AnalyzerAllSteps::dz_line_point(KsTrack1creationVertex,KsTrack1Momentum,PVmin);
        double dz_Track1Ks_PV0 = AnalyzerAllSteps::dz_line_point(KsTrack1creationVertex,KsTrack1Momentum,PV0);
        double dz_Track1Ks_000 = AnalyzerAllSteps::dz_line_point(KsTrack1creationVertex,KsTrack1Momentum,ZeroZeroZero);

	TVector3 KsTrack2creationVertex(RECOKs->daughter(1)->vx(),RECOKs->daughter(1)->vy(),RECOKs->daughter(1)->vz());
        TVector3 KsTrack2Momentum(RECOKs->daughter(1)->px(),RECOKs->daughter(1)->py(),RECOKs->daughter(1)->pz());
	double dxy_Track2Ks_beamspot = AnalyzerAllSteps::dxy_signed_line_point(KsTrack2creationVertex,KsTrack2Momentum,beamspot);	
        double dz_Track2Ks_beamspot = AnalyzerAllSteps::dz_line_point(KsTrack2creationVertex,KsTrack2Momentum,beamspot);
        double dz_Track2Ks_min_PV = AnalyzerAllSteps::dz_line_point(KsTrack2creationVertex,KsTrack2Momentum,PVmin);
        double dz_Track2Ks_PV0 = AnalyzerAllSteps::dz_line_point(KsTrack2creationVertex,KsTrack2Momentum,PV0);
        double dz_Track2Ks_000 = AnalyzerAllSteps::dz_line_point(KsTrack2creationVertex,KsTrack2Momentum,ZeroZeroZero);

	//_Ks_trackPair_mindeltaR.push_back(deltaRRecoKsTrackPair);
	_Ks_trackPair_mass.push_back(p4bestTrackPair.M());

	//between the two tracks
	_Ks_Track1Track2_openingsAngle.push_back(openingsAngleTrack1Track2);
	_Ks_Track1Track2_deltaR.push_back(deltaRTrack1Track2);

	//between the two tracks and the Ks
	_Ks_Track1_openingsAngle.push_back(openingsAngleTrack1Ks);
	_Ks_Track2_openingsAngle.push_back(openingsAngleTrack2Ks);
	_Ks_Track1_deltaR.push_back(deltaRTrack1Ks);
	_Ks_Track2_deltaR.push_back(deltaRTrack2Ks);
	
	//somehow the RECOKs->daughter(0) should be containing the reference to the track (see e.g. https://github.com/jarnedc/cmssw/blob/from-CMSSW_8_0_30/RecoVertex/V0Producer/src/V0Fitter.cc#L515), but I am not able to get it out...
	//const reco::Candidate* Ks_daug0 = RECOKs->daughter(0);
	//reco::TrackRef Ks_daug0_trackref = Ks_daug0->track();

	
	_Ks_daughterTrack1_charge.push_back(RECOKs->daughter(0)->charge());
	//_Ks_daughterTrack1_chi2.push_back(RECOKs->daughter(0)->chi2());
	//_Ks_daughterTrack1_ndof.push_back(RECOKs->daughter(0)->ndof());
	_Ks_daughterTrack1_eta.push_back(RECOKs->daughter(0)->eta());
	_Ks_daughterTrack1_phi.push_back(RECOKs->daughter(0)->phi());
	_Ks_daughterTrack1_pt.push_back(RECOKs->daughter(0)->pt());
	_Ks_daughterTrack1_pz.push_back(RECOKs->daughter(0)->pz());
	_Ks_daughterTrack1_dxy_beamspot.push_back(dxy_Track1Ks_beamspot);
	_Ks_daughterTrack1_dz_beamspot.push_back(dz_Track1Ks_beamspot);
	_Ks_daughterTrack1_dz_min_PV.push_back(dz_Track1Ks_min_PV);
	_Ks_daughterTrack1_dz_PV0.push_back(dz_Track1Ks_PV0);
	_Ks_daughterTrack1_dz_000.push_back(dz_Track1Ks_000);

	_Ks_daughterTrack2_charge.push_back(RECOKs->daughter(1)->charge());
	//_Ks_daughterTrack2_chi2.push_back(RECOKs->daughter(1)->chi2());
	//_Ks_daughterTrack2_ndof.push_back(RECOKs->daughter(1)->ndof());
	_Ks_daughterTrack2_eta.push_back(RECOKs->daughter(1)->eta());
	_Ks_daughterTrack2_phi.push_back(RECOKs->daughter(1)->phi());
	_Ks_daughterTrack2_pt.push_back(RECOKs->daughter(1)->pt());
	_Ks_daughterTrack2_pz.push_back(RECOKs->daughter(1)->pz());
	_Ks_daughterTrack2_dxy_beamspot.push_back(dxy_Track2Ks_beamspot);
	_Ks_daughterTrack2_dz_beamspot.push_back(dz_Track2Ks_beamspot);
	_Ks_daughterTrack2_dz_min_PV.push_back(dz_Track2Ks_min_PV);
	_Ks_daughterTrack2_dz_PV0.push_back(dz_Track2Ks_PV0);
	_Ks_daughterTrack2_dz_000.push_back(dz_Track2Ks_000);

	

}

void FlatTreeProducerV0s::FillBranchesLambda(const reco::VertexCompositeCandidate * RECOLambda, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV, edm::Handle<vector<reco::GenParticle>> h_genParticles){


	TVector3 LambdaCreationVertex(RECOLambda->vx(),RECOLambda->vy(),RECOLambda->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,LambdaCreationVertex);
	TVector3 LambdaMomentum(RECOLambda->px(),RECOLambda->py(),RECOLambda->pz());

	//different reference points
	TVector3 PVmin = AnalyzerAllSteps::dz_line_point_min(LambdaCreationVertex,LambdaMomentum,h_offlinePV);
	TVector3 PV0(h_offlinePV->at(0).x(),h_offlinePV->at(0).y(),h_offlinePV->at(0).z());
	TVector3 ZeroZeroZero(0.,0.,0.);
	//different dxy
	double dxy_beamspot = AnalyzerAllSteps::dxy_signed_line_point(LambdaCreationVertex,LambdaMomentum,beamspot);
	double dxy_min_PV = AnalyzerAllSteps::dxy_signed_line_point(LambdaCreationVertex,LambdaMomentum,PVmin);
	double dxy_PV0 = AnalyzerAllSteps::dxy_signed_line_point(LambdaCreationVertex,LambdaMomentum,PV0);
	double dxy_000 = AnalyzerAllSteps::dxy_signed_line_point(LambdaCreationVertex,LambdaMomentum,ZeroZeroZero);

	//different dz's
	double dz_beamspot = AnalyzerAllSteps::dz_line_point(LambdaCreationVertex,LambdaMomentum,beamspot);
	double dz_min_PV = AnalyzerAllSteps::dz_line_point(LambdaCreationVertex,LambdaMomentum,PVmin);
	double dz_PV0 = AnalyzerAllSteps::dz_line_point(LambdaCreationVertex,LambdaMomentum,PV0);
        double dz_000 = AnalyzerAllSteps::dz_line_point(LambdaCreationVertex,LambdaMomentum,ZeroZeroZero);


	_Lambda_mass.push_back(RECOLambda->mass());	

	_Lambda_pt.push_back(RECOLambda->pt());	
	_Lambda_pz.push_back(RECOLambda->pz());	
	_Lambda_Lxy.push_back(Lxy);	
	_Lambda_vz.push_back(RECOLambda->vz());	

	_Lambda_eta.push_back(RECOLambda->eta());
	_Lambda_phi.push_back(RECOLambda->phi());

	_Lambda_dxy.push_back(dxy_beamspot);
	_Lambda_dxy_min_PV.push_back(dxy_min_PV);
	_Lambda_dxy_PV0.push_back(dxy_PV0);
	_Lambda_dxy_000.push_back(dxy_000);
	
	_Lambda_dz_beamspot.push_back(dz_beamspot);	
	_Lambda_dz_min_PV.push_back(dz_min_PV);	
	_Lambda_dz_PV0.push_back(dz_PV0);	
	_Lambda_dz_000.push_back(dz_000);	

	_Lambda_vz_dz_min_PV.push_back(PVmin.Z());


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

	_Lambda_dxy.clear();
	_Lambda_dxy_min_PV.clear();
	_Lambda_dxy_PV0.clear();
	_Lambda_dxy_000.clear();

	_Lambda_dz_beamspot.clear();
	_Lambda_dz_min_PV.clear();	
	_Lambda_dz_PV0.clear();	
	_Lambda_dz_000.clear();	

	_Lambda_vz_dz_min_PV.clear();	


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

DEFINE_FWK_MODULE(FlatTreeProducerV0s);
