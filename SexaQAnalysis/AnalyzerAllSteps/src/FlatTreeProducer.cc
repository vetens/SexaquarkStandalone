#include "../interface/FlatTreeProducer.h"
#include <typeinfo>

FlatTreeProducer::FlatTreeProducer(edm::ParameterSet const& pset):
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



  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_offlinePVToken    (consumes<vector<reco::Vertex>>(m_offlinePVTag)),
  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM_GEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM_GEANT)),
  m_generalTracksToken(consumes<View<reco::Track> >(m_generalTracksTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_V0KsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0KsTag)),
  m_V0LToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0LTag))
  


{

}


void FlatTreeProducer::beginJob() {

    
        // Initialize when class is created
        edm::Service<TFileService> fs ;
        _tree = fs->make <TTree>("FlatTree","tree");

        // Declare tree's branches
        // Event
        //just for checking you are looking at the correct one: S or antiS
	_tree->Branch("_S_charge",&_S_charge);
	_tree->Branch("_S_deltaRmin_GEN_RECO",&_S_deltaRmin_GEN_RECO);
	_tree->Branch("_S_lxy_interaction_vertex",&_S_lxy_interaction_vertex);
	_tree->Branch("_S_error_lxy_interaction_vertex",&_S_error_lxy_interaction_vertex);
	_tree->Branch("_Ks_lxy_decay_vertex",&_Ks_lxy_decay_vertex);
	_tree->Branch("_Lambda_lxy_decay_vertex",&_Lambda_lxy_decay_vertex);
	_tree->Branch("_S_mass",&_S_mass);
	_tree->Branch("_S_chi2_ndof",&_S_chi2_ndof);

	_tree->Branch("_S_daughters_deltaphi",&_S_daughters_deltaphi);
	_tree->Branch("_S_daughters_deltaeta",&_S_daughters_deltaeta);
	_tree->Branch("_S_daughters_openingsangle",&_S_daughters_openingsangle);
	_tree->Branch("_S_Ks_openingsangle",&_S_Ks_openingsangle);
	_tree->Branch("_S_Lambda_openingsangle",&_S_Lambda_openingsangle);
	_tree->Branch("_S_daughters_DeltaR",&_S_daughters_DeltaR);
	_tree->Branch("_S_eta",&_S_eta);
	_tree->Branch("_Ks_eta",&_Ks_eta);
	_tree->Branch("_Lambda_eta",&_Lambda_eta);

	_tree->Branch("_S_dxy",&_S_dxy);
	_tree->Branch("_Ks_dxy",&_Ks_dxy);
	_tree->Branch("_Lambda_dxy",&_Lambda_dxy);
	_tree->Branch("_S_dxy_dzPVmin",&_S_dxy_dzPVmin);
	_tree->Branch("_Ks_dxy_dzPVmin",&_Ks_dxy_dzPVmin);
	_tree->Branch("_Lambda_dxy_dzPVmin",&_Lambda_dxy_dzPVmin);

	_tree->Branch("_S_dxy_over_lxy",&_S_dxy_over_lxy);
	_tree->Branch("_Ks_dxy_over_lxy",&_Ks_dxy_over_lxy);
	_tree->Branch("_Lambda_dxy_over_lxy",&_Lambda_dxy_over_lxy);

	_tree->Branch("_S_dz",&_S_dz);
	_tree->Branch("_Ks_dz",&_Ks_dz);
	_tree->Branch("_Lambda_dz",&_Lambda_dz);
	_tree->Branch("_S_dz_min",&_S_dz_min);
	_tree->Branch("_Ks_dz_min",&_Ks_dz_min);
	_tree->Branch("_Lambda_dz_min",&_Lambda_dz_min);

	_tree->Branch("_S_pt",&_S_pt);
	_tree->Branch("_Ks_pt",&_Ks_pt);
	_tree->Branch("_Lambda_pt",&_Lambda_pt);

	_tree->Branch("_S_pz",&_S_pz);
	_tree->Branch("_Ks_pz",&_Ks_pz);
	_tree->Branch("_Lambda_pz",&_Lambda_pz);

	_tree->Branch("_S_vz_interaction_vertex",&_S_vz_interaction_vertex);
	_tree->Branch("_Ks_vz_decay_vertex",&_Ks_vz_decay_vertex);
	_tree->Branch("_Lambda_vz_decay_vertex",&_Lambda_vz_decay_vertex);

	_tree->Branch("_S_vx",&_S_vx);
	_tree->Branch("_S_vy",&_S_vy);
	_tree->Branch("_S_vz",&_S_vz);



}

void FlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  //iEvent.getByToken(m_bsToken, h_bs);

  //primary vertex
  edm::Handle<vector<reco::Vertex>> h_offlinePV;
  iEvent.getByToken(m_offlinePVToken, h_offlinePV);

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



  //beamspot
  TVector3 beamspot(0,0,0);
  TVector3 beamspotVariance(0,0,0);
  if(h_bs.isValid()){  
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));			
  }


  //in case of data: just fill the flatTree with data from all the reconstructed S particles.
  if(m_runningOnData){
	  //loop over the RECO AntiS to plot the kinematics
	  if(h_sCands.isValid()){
	      for(unsigned int i = 0; i < h_sCands->size(); ++i){//loop all S candidates
		const reco::VertexCompositeCandidate * antiS = &h_sCands->at(i);
		//the if below is important, if it is -1 you are looking at signal (antiS). If it is +1 you are looking at background (S).
		int chargeAntiProton = 1; //by default look at the background
		if(m_lookAtAntiS) chargeAntiProton = -1;//only when the m_lookAtAntiS flag is enabled look at the antiS, which has a charge of -1 (representing the charge of the proton)
		if(antiS->charge()==chargeAntiProton){
			FillBranches(antiS, beamspot, beamspotVariance, h_offlinePV,-999);
		}
	      }
	  }
  }


  //in case of MC: first find the GEN particles which are proper antiS and then find the corresponding RECO particles
  else if(!m_runningOnData && m_lookAtAntiS){
	  if(h_genParticles.isValid()){
	      for(unsigned int i = 0; i < h_genParticles->size(); ++i){//loop all genparticlesPlusGEANT
			const reco::Candidate * genParticle = &h_genParticles->at(i);
			bool genParticleIsAntiS = false;
			if(genParticle->pdgId() == AnalyzerAllSteps::pdgIdAntiS) {genParticleIsAntiS = true;}

			//find the antiS:
			if(genParticle->numberOfDaughters()==2){
				int daughterParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle);
				//for antiS
				if(genParticle->daughter(0)->numberOfDaughters()==2 && genParticle->daughter(1)->numberOfDaughters()==2 && daughterParticlesTypes == 3){
					int graddaughters0ParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle->daughter(0));
					int graddaughters1ParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle->daughter(1));
					if(genParticleIsAntiS && ((graddaughters0ParticlesTypes == 1 && graddaughters1ParticlesTypes == 2) || (graddaughters1ParticlesTypes == 1 && graddaughters0ParticlesTypes == 2))){ FindRecoAntiS(genParticle,h_sCands,beamspot, beamspotVariance, h_offlinePV);}
				}
			}
	      }//for(unsigned int i = 0; i < h_genParticles->size(); ++i)
	  }//if(h_genParticles.isValid())
  }


 } //end of analyzer



//Find the best matching RECO AntiS which matches the GEN S best in deltaR
void FlatTreeProducer::FindRecoAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV){
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

		if(deltaRmin<AnalyzerAllSteps::deltaRCutRECOAntiS)FillBranches(bestRECOAntiS, beamspot, beamspotVariance, h_offlinePV, deltaRmin);
       }
}

void FlatTreeProducer::FillBranches(const reco::VertexCompositeCandidate * RECO_S, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV, double deltaRmin){
	nTotalRECOS++;
	//calculate some kinematic variables for the RECO AntiS
	TVector3 RECOAntiSInteractionVertex(RECO_S->vx(),RECO_S->vy(),RECO_S->vz());//this is the interaction vertex of the antiS and the neutron. Check in the skimming code if you want to check.
	TVector3 RECOAntiSMomentumVertex(RECO_S->px(),RECO_S->py(),RECO_S->pz());
	double RECOLxy_interactionVertex = AnalyzerAllSteps::lxy(beamspot,RECOAntiSInteractionVertex);
	double RECOErrorLxy_interactionVertex = AnalyzerAllSteps::std_dev_lxy(RECO_S->vx(), RECO_S->vy(), RECO_S->vertexCovariance(0,0), RECO_S->vertexCovariance(1,1), beamspot.X(), beamspot.Y(), beamspotVariance.X(), beamspotVariance.Y());
	double RECODeltaPhiDaughters = reco::deltaPhi(RECO_S->daughter(0)->phi(),RECO_S->daughter(1)->phi());
	double RECODeltaEtaDaughters = RECO_S->daughter(0)->eta()-RECO_S->daughter(1)->eta();
	double RECODeltaRDaughters = pow(RECODeltaPhiDaughters*RECODeltaPhiDaughters+RECODeltaEtaDaughters*RECODeltaEtaDaughters,0.5);

	reco::LeafCandidate::LorentzVector n_(0,0,0,0.939565);
	double RECO_Smass = (RECO_S->p4()-n_).mass();

	//the lxy of the Lambda and Ks decay vertex
	TVector3 RECOAntiSDaug0Vertex(RECO_S->daughter(0)->vx(),RECO_S->daughter(0)->vy(),RECO_S->daughter(0)->vz());
        TVector3 RECOAntiSDaug1Vertex(RECO_S->daughter(1)->vx(),RECO_S->daughter(1)->vy(),RECO_S->daughter(1)->vz());
	double RECOLxy_Lambda = AnalyzerAllSteps::lxy(beamspot,RECOAntiSDaug0Vertex);
	double RECOLxy_Ks = AnalyzerAllSteps::lxy(beamspot,RECOAntiSDaug1Vertex);

	//the dxy of the Lambda and Ks
	TVector3 RECOAntiSDaug0Momentum(RECO_S->daughter(0)->px(),RECO_S->daughter(0)->py(),RECO_S->daughter(0)->pz());
        TVector3 RECOAntiSDaug1Momentum(RECO_S->daughter(1)->px(),RECO_S->daughter(1)->py(),RECO_S->daughter(1)->pz());
	reco::Candidate::Vector vRECOAntiSMomentum(RECO_S->px(),RECO_S->py(),RECO_S->pz());
	reco::Candidate::Vector vRECOAntiSDaug0Momentum(RECO_S->daughter(0)->px(),RECO_S->daughter(0)->py(),RECO_S->daughter(0)->pz());
        reco::Candidate::Vector vRECOAntiSDaug1Momentum(RECO_S->daughter(1)->px(),RECO_S->daughter(1)->py(),RECO_S->daughter(1)->pz());

	double RECOOpeningsAngleAntiSLambda = AnalyzerAllSteps::openings_angle(vRECOAntiSDaug0Momentum,vRECOAntiSMomentum);
	double RECOOpeningsAngleAntiSKs = AnalyzerAllSteps::openings_angle(vRECOAntiSDaug1Momentum,vRECOAntiSMomentum);

	double RECOOpeningsAngleDaughters = AnalyzerAllSteps::openings_angle(vRECOAntiSDaug0Momentum,vRECOAntiSDaug1Momentum);
        double RECO_dxy_daughter0 = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot);
        double RECO_dxy_daughter1 = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot);
	//the dz of the Ks and Lambda
	double RECO_dz_daughter0 = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot);
	double RECO_dz_daughter1 = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot);
	//the dxyz of the Ks and Lambda
//	double RECO_dxyz_daughter0 = AnalyzerAllSteps::dxyz_signed_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot);
//	double RECO_dxyz_daughter1 = AnalyzerAllSteps::dxyz_signed_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot);
	//collapse the RECO_dxy, RECO_dz, RECO_dxyz variables on one variable
//	double relDiff_RECO_dxy_daughters = (abs(RECO_dxy_daughter0)-abs(RECO_dxy_daughter1))/(abs(RECO_dxy_daughter0)+abs(RECO_dxy_daughter1));
//	double relDiff_RECO_dz_daughters = (abs(RECO_dz_daughter0)-abs(RECO_dz_daughter1))/(abs(RECO_dz_daughter0)+abs(RECO_dz_daughter1));
//	double relDiff_RECO_dxyz_daughters = (abs(RECO_dxyz_daughter0)-abs(RECO_dxyz_daughter1))/(abs(RECO_dxyz_daughter0)+abs(RECO_dxyz_daughter1));
	//dxy and dz of the AntiS itself
	double RECO_dxy_antiS = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,beamspot);
	double RECO_dz_antiS = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,beamspot);
	

	//calculate the sign of the dot product between the displacement vector and the dxy vector for both the Ks and the Lambda
//	double signLxyDotdxy_daughter0 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot)*RECOAntiSInteractionVertex);
//	double signLxyDotdxy_daughter1 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot)*RECOAntiSInteractionVertex);
	//calculate the sign of the dot product between the displacement vector and the pt of both the Ks and the Lambda
//	double signPtDotdxy_daughter0 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug0Momentum,beamspot)*RECOAntiSDaug0Momentum);
//	double signPtDotdxy_daughter1 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(RECOAntiSInteractionVertex, RECOAntiSDaug1Momentum,beamspot)*RECOAntiSDaug1Momentum);

	//loop over all PVs and find the one which minimises the dz of the antiS
	double RECOdzAntiSPVmin = 999.;
	double dxyAntiSPVmin = 999.;
	TVector3  bestPVdzAntiS;
	if(h_offlinePV.isValid()){
		bestPVdzAntiS = AnalyzerAllSteps::dz_line_point_min(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,h_offlinePV);
		RECOdzAntiSPVmin = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,bestPVdzAntiS);
		dxyAntiSPVmin = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex,RECOAntiSMomentumVertex,bestPVdzAntiS);
	}

	double RECOdzLambdaPVmin = 999.;
	double dxyLambdaPVmin = 999.;
	TVector3  bestPVdzLambda;
	if(h_offlinePV.isValid()){
		bestPVdzLambda = AnalyzerAllSteps::dz_line_point_min(RECOAntiSInteractionVertex,RECOAntiSDaug0Momentum,h_offlinePV);
		RECOdzLambdaPVmin = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex,RECOAntiSDaug0Momentum,bestPVdzLambda);
		dxyLambdaPVmin = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex,RECOAntiSDaug0Momentum,bestPVdzLambda);
	}


	double RECOdzKsPVmin = 999.;
	double dxyKsPVmin = 999.;
	TVector3  bestPVdzKs;
	if(h_offlinePV.isValid()){
		bestPVdzKs = AnalyzerAllSteps::dz_line_point_min(RECOAntiSInteractionVertex,RECOAntiSDaug1Momentum,h_offlinePV);
		RECOdzKsPVmin = AnalyzerAllSteps::dz_line_point(RECOAntiSInteractionVertex,RECOAntiSDaug1Momentum,bestPVdzKs);
		dxyKsPVmin = AnalyzerAllSteps::dxy_signed_line_point(RECOAntiSInteractionVertex,RECOAntiSDaug1Momentum,bestPVdzKs);
	}



	//if the RECO S particle fails one of the below cuts than don't fill the tree. These already cut the majority of the background, so the background trees will be much smaller, which is nice.
	if(RECOLxy_interactionVertex < AnalyzerAllSteps::MinLxyCut || RECOErrorLxy_interactionVertex > AnalyzerAllSteps::MaxErrorLxyCut || RECO_Smass < 0.)return;


	nSavedRECOS++;

	Init();	

	_S_charge.push_back(RECO_S->charge());
	_S_deltaRmin_GEN_RECO.push_back(deltaRmin);

	_S_lxy_interaction_vertex.push_back(RECOLxy_interactionVertex);
	_S_error_lxy_interaction_vertex.push_back(RECOErrorLxy_interactionVertex);
	_Ks_lxy_decay_vertex.push_back(RECOLxy_Ks);
	_Lambda_lxy_decay_vertex.push_back(RECOLxy_Lambda);
	_S_mass.push_back(RECO_Smass);
	_S_chi2_ndof.push_back(RECO_S->vertexNormalizedChi2());

	_S_daughters_deltaphi.push_back(RECODeltaPhiDaughters);
	_S_daughters_deltaeta.push_back(RECODeltaEtaDaughters);
	_S_daughters_openingsangle.push_back(RECOOpeningsAngleDaughters);
	_S_Ks_openingsangle.push_back(RECOOpeningsAngleAntiSKs);
	_S_Lambda_openingsangle.push_back(RECOOpeningsAngleAntiSLambda);
	_S_daughters_DeltaR.push_back(RECODeltaRDaughters);
	_S_eta.push_back(RECO_S->eta());
	_Lambda_eta.push_back(RECO_S->daughter(0)->eta());
	_Ks_eta.push_back(RECO_S->daughter(1)->eta());

	_S_dxy.push_back(RECO_dxy_antiS);
	_Lambda_dxy.push_back(RECO_dxy_daughter0);
	_Ks_dxy.push_back(RECO_dxy_daughter1);
	_S_dxy_dzPVmin.push_back(dxyAntiSPVmin);
	_Ks_dxy_dzPVmin.push_back(dxyKsPVmin);
	_Lambda_dxy_dzPVmin.push_back(dxyLambdaPVmin);

	_S_dxy_over_lxy.push_back(RECO_dxy_antiS/RECOLxy_interactionVertex);
	_Lambda_dxy_over_lxy.push_back(RECO_dxy_daughter0/RECOLxy_interactionVertex);
	_Ks_dxy_over_lxy.push_back(RECO_dxy_daughter1/RECOLxy_interactionVertex);

	_S_dz.push_back(RECO_dz_antiS);
	_Lambda_dz.push_back(RECO_dz_daughter0);
	_Ks_dz.push_back(RECO_dz_daughter1);
	_S_dz_min.push_back(RECOdzAntiSPVmin);
	_Ks_dz_min.push_back(RECOdzKsPVmin);
	_Lambda_dz_min.push_back(RECOdzLambdaPVmin);

	_S_pt.push_back(RECO_S->pt());
	_Lambda_pt.push_back(RECO_S->daughter(0)->pt());
	_Ks_pt.push_back(RECO_S->daughter(1)->pt());
	
	_S_pz.push_back(RECO_S->pz());
	_Lambda_pz.push_back(RECO_S->daughter(0)->pz());
	_Ks_pz.push_back(RECO_S->daughter(1)->pz());

	_S_vz_interaction_vertex.push_back(RECO_S->vz());
        _Lambda_vz_decay_vertex.push_back(RECOAntiSDaug0Vertex.Z());
        _Ks_vz_decay_vertex.push_back(RECOAntiSDaug1Vertex.Z());

	_S_vx.push_back(RECO_S->vx());	
	_S_vy.push_back(RECO_S->vy());	
	_S_vz.push_back(RECO_S->vz());	

  	_tree->Fill();

}


void FlatTreeProducer::endJob()
{
}

void
FlatTreeProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void
FlatTreeProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
FlatTreeProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
FlatTreeProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

FlatTreeProducer::~FlatTreeProducer()
{
	if(m_lookAtAntiS){
			std::cout << "The total number RECO anti-S that were found is: " << nTotalRECOS << std::endl; 
			std::cout << "The total number RECO anti-S that were saved is: " << nSavedRECOS << std::endl; 
	}
	if(!m_lookAtAntiS){
			std::cout << "The total number RECO S that were found is: " << nTotalRECOS << std::endl; 
			std::cout << "The total number RECO S that were saved is: " << nSavedRECOS << std::endl; 
	}
	std::cout << "saved/found = " << (double)nSavedRECOS/(double)nTotalRECOS << std::endl;
}


void
FlatTreeProducer::Init()
{


    	_S_charge.clear();
    	_S_deltaRmin_GEN_RECO.clear();

    	_S_lxy_interaction_vertex.clear();
        _S_error_lxy_interaction_vertex.clear();
	_Ks_lxy_decay_vertex.clear();
	_Lambda_lxy_decay_vertex.clear();
        _S_mass.clear();
        _S_chi2_ndof.clear();

	_S_daughters_deltaphi.clear();
	_S_daughters_deltaeta.clear();
	_S_daughters_openingsangle.clear();
	_S_Ks_openingsangle.clear();
	_S_Lambda_openingsangle.clear();
	_S_daughters_DeltaR.clear();
	_S_eta.clear();
	_Ks_eta.clear();
	_Lambda_eta.clear();

	_S_dxy.clear();
	_Ks_dxy.clear();
	_Lambda_dxy.clear();
	_S_dxy_dzPVmin.clear();
	_Ks_dxy_dzPVmin.clear();
	_Lambda_dxy_dzPVmin.clear();

	_S_dxy_over_lxy.clear();
	_Ks_dxy_over_lxy.clear();
	_Lambda_dxy_over_lxy.clear();

	_S_dz.clear();
	_Ks_dz.clear();
	_Lambda_dz.clear();
	_S_dz_min.clear();
	_Ks_dz_min.clear();
	_Lambda_dz_min.clear();

	_S_pt.clear();
	_Ks_pt.clear();
	_Lambda_pt.clear();

	_S_pz.clear();
	_Ks_pz.clear();
	_Lambda_pz.clear();

	_S_vz_interaction_vertex.clear();
	_Ks_vz_decay_vertex.clear();
	_Lambda_vz_decay_vertex.clear();

	_S_vx.clear();
	_S_vy.clear();
	_S_vz.clear();


}

DEFINE_FWK_MODULE(FlatTreeProducer);
