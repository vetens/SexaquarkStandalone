#include "../interface/FlatTreeProducerBDT.h"
#include <typeinfo>

FlatTreeProducerBDT::FlatTreeProducerBDT(edm::ParameterSet const& pset):
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


void FlatTreeProducerBDT::beginJob() {

        // Initialize when class is created
        edm::Service<TFileService> fs ;

	//PV information
        _tree_PV = fs->make <TTree>("FlatTreePV","tree_PV");
	_tree_PV->Branch("_nPV",&_nPV);
	_tree_PV->Branch("_nGoodPV",&_nGoodPV);
	_tree_PV->Branch("_nGoodPVPOG",&_nGoodPVPOG);
	_tree_PV->Branch("_PVx",&_PVx);
	_tree_PV->Branch("_PVy",&_PVy);
	_tree_PV->Branch("_PVz",&_PVz);
	_tree_PV->Branch("_goodPVx",&_goodPVx);
	_tree_PV->Branch("_goodPVy",&_goodPVy);
	_tree_PV->Branch("_goodPVz",&_goodPVz);
	_tree_PV->Branch("_goodPVxPOG",&_goodPVxPOG);
	_tree_PV->Branch("_goodPVyPOG",&_goodPVyPOG);
	_tree_PV->Branch("_goodPVzPOG",&_goodPVzPOG);
    
        // Initialize when class is created
        _tree = fs->make <TTree>("FlatTree","tree");

        // Declare tree's branches
        // Event
        //just for checking you are looking at the correct one: S or antiS
	_tree->Branch("_S_charge",&_S_charge);
	_tree->Branch("_S_deltaLInteractionVertexAntiSmin",&_S_deltaLInteractionVertexAntiSmin);
	_tree->Branch("_S_deltaRAntiSmin",&_S_deltaRAntiSmin);
	_tree->Branch("_S_deltaRKsAntiSmin",&_S_deltaRKsAntiSmin);
	_tree->Branch("_S_deltaRLambdaAntiSmin",&_S_deltaRLambdaAntiSmin);

	_tree->Branch("_S_lxy_interaction_vertex",&_S_lxy_interaction_vertex);
	_tree->Branch("_S_lxy_interaction_vertex_beampipeCenter",&_S_lxy_interaction_vertex_beampipeCenter);
	_tree->Branch("_S_error_lxy_interaction_vertex",&_S_error_lxy_interaction_vertex);
	_tree->Branch("_S_error_lxy_interaction_vertex_beampipeCenter",&_S_error_lxy_interaction_vertex_beampipeCenter);
	_tree->Branch("_Ks_lxy_decay_vertex",&_Ks_lxy_decay_vertex);
	_tree->Branch("_Lambda_lxy_decay_vertex",&_Lambda_lxy_decay_vertex);
	_tree->Branch("_S_mass",&_S_mass);
	_tree->Branch("_S_chi2_ndof",&_S_chi2_ndof);
	_tree->Branch("_S_event_weighting_factor",&_S_event_weighting_factor);
	_tree->Branch("_S_event_weighting_factorPU",&_S_event_weighting_factorPU);
	_tree->Branch("_S_event_weighting_factorALL",&_S_event_weighting_factorALL);

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

	_tree->Branch("_Lambda_mass",&_Lambda_mass);
	_tree->Branch("_Ks_mass",&_Ks_mass);

	_tree->Branch("_RECO_Lambda_daughter0_charge",&_RECO_Lambda_daughter0_charge);
	_tree->Branch("_RECO_Lambda_daughter0_pt",&_RECO_Lambda_daughter0_pt);
	_tree->Branch("_RECO_Lambda_daughter0_pz",&_RECO_Lambda_daughter0_pz);
	_tree->Branch("_RECO_Lambda_daughter0_dxy_beamspot",&_RECO_Lambda_daughter0_dxy_beamspot);
	_tree->Branch("_RECO_Lambda_daughter0_dz_beamspot",&_RECO_Lambda_daughter0_dz_beamspot);

	_tree->Branch("_RECO_Lambda_daughter1_charge",&_RECO_Lambda_daughter1_charge);
	_tree->Branch("_RECO_Lambda_daughter1_pt",&_RECO_Lambda_daughter1_pt);
	_tree->Branch("_RECO_Lambda_daughter1_pz",&_RECO_Lambda_daughter1_pz);
	_tree->Branch("_RECO_Lambda_daughter1_dxy_beamspot",&_RECO_Lambda_daughter1_dxy_beamspot);
	_tree->Branch("_RECO_Lambda_daughter1_dz_beamspot",&_RECO_Lambda_daughter1_dz_beamspot);

	
	_tree->Branch("_RECO_Ks_daughter0_charge",&_RECO_Ks_daughter0_charge);
	_tree->Branch("_RECO_Ks_daughter0_pt",&_RECO_Ks_daughter0_pt);
	_tree->Branch("_RECO_Ks_daughter0_pz",&_RECO_Ks_daughter0_pz);
	_tree->Branch("_RECO_Ks_daughter0_dxy_beamspot",&_RECO_Ks_daughter0_dxy_beamspot);
	_tree->Branch("_RECO_Ks_daughter0_dz_beamspot",&_RECO_Ks_daughter0_dz_beamspot);

	_tree->Branch("_RECO_Ks_daughter1_charge",&_RECO_Ks_daughter1_charge);
	_tree->Branch("_RECO_Ks_daughter1_pt",&_RECO_Ks_daughter1_pt);
	_tree->Branch("_RECO_Ks_daughter1_pz",&_RECO_Ks_daughter1_pz);
	_tree->Branch("_RECO_Ks_daughter1_dxy_beamspot",&_RECO_Ks_daughter1_dxy_beamspot);
	_tree->Branch("_RECO_Ks_daughter1_dz_beamspot",&_RECO_Ks_daughter1_dz_beamspot);

        _tree_counter = fs->make <TTree>("FlatTreeCounter","tree_counter");
	_tree_counter->Branch("_nGENAntiSWithCorrectGranddaughters",&_nGENAntiSWithCorrectGranddaughters);


}

void FlatTreeProducerBDT::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

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

  //PV
  int nPVs = 0;
  int ngoodPVs = 0;
  unsigned int ngoodPVsPOG = 0;
  Init_PV();
  if(h_offlinePV.isValid()){
	for(unsigned int i = 0; i < h_offlinePV->size(); i++ ){
		if(h_offlinePV->at(i).isValid()){
			nPVs++;
			_PVx.push_back(h_offlinePV->at(i).x());
			_PVy.push_back(h_offlinePV->at(i).y());
			_PVz.push_back(h_offlinePV->at(i).z());
		}
		if(h_offlinePV->at(i).isValid() && h_offlinePV->at(i).tracksSize() >= 4){
			ngoodPVs++;
			_goodPVx.push_back(h_offlinePV->at(i).x());
			_goodPVy.push_back(h_offlinePV->at(i).y());
			_goodPVz.push_back(h_offlinePV->at(i).z());
		}

                double r = sqrt(h_offlinePV->at(i).x()*h_offlinePV->at(i).x()+h_offlinePV->at(i).y()*h_offlinePV->at(i).y());
                if(h_offlinePV->at(i).ndof() > 4 && abs(h_offlinePV->at(i).z()) < 24 && r < 2){
			ngoodPVsPOG++;
                        _goodPVxPOG.push_back(h_offlinePV->at(i).x());
                        _goodPVyPOG.push_back(h_offlinePV->at(i).y());
                        _goodPVzPOG.push_back(h_offlinePV->at(i).z());
		}


	}
  }
//  srand(time(NULL));
  int randomIndexPV = rand()%(ngoodPVsPOG);
  double randomPVz = _goodPVzPOG[randomIndexPV];
  

  _nPV.push_back(nPVs);
  _nGoodPV.push_back(ngoodPVs);
  _nGoodPVPOG.push_back(ngoodPVsPOG);
  //_tree_PV->Fill();

  //beamspot
  TVector3 beamspot(999999,999999,999999);
  TVector3 beamspotVariance(999999,999999,999999);
  if(h_bs.isValid()){  
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));			
  }

  int nGENAntiSWithCorrectGranddaughtersThisEvent = 0;
  
  //for both data and MC: loop over all entries in h_sCands, both the ones with positive and the ones with negative charge. For the MC ones with negative charge I will check if they have a matching GEN antiS 
  if(h_sCands.isValid()){
      for(unsigned int i = 0; i < h_sCands->size(); ++i){//loop all S candidates
	const reco::VertexCompositeCandidate * antiS = &h_sCands->at(i);
	FillBranches(antiS, beamspot, beamspotVariance, h_offlinePV,m_runningOnData,  h_genParticles, h_V0Ks,  h_V0L, ngoodPVsPOG, randomPVz);
      }
  }

  Init_Counter();
  _nGENAntiSWithCorrectGranddaughters.push_back(nGENAntiSWithCorrectGranddaughtersThisEvent);
  _tree_counter->Fill();

 } //end of analyzer




void FlatTreeProducerBDT::FillBranches(const reco::VertexCompositeCandidate * RECO_S, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV, bool m_runningOnData, edm::Handle<vector<reco::GenParticle>> h_genParticles, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L, unsigned int ngoodPVsPOG, double randomPVz){

	TVector3 RECOAntiSInteractionVertex(RECO_S->vx(),RECO_S->vy(),RECO_S->vz());//this is the interaction vertex of the antiS and the neutron. Check in the skimming code if you want to check.
	double RECOLxy_interactionVertex_beampipeCenter = sqrt(RECOAntiSInteractionVertex.X()*RECOAntiSInteractionVertex.X() + RECOAntiSInteractionVertex.Y()*RECOAntiSInteractionVertex.Y() );
	//indeed, if you are running on MC then the above calculation, which is with respect to (0,0,0) is correct for RECOLxy_interactionVertex_beampipeCenter, but if you run on data then you should actually calculate wrt the center of the beampipe, which is offsset wrt (0,0,0)
	if(m_runningOnData) RECOLxy_interactionVertex_beampipeCenter = sqrt( pow(RECOAntiSInteractionVertex.X()-AnalyzerAllSteps::center_beampipe_x , 2) + pow(RECOAntiSInteractionVertex.Y()-AnalyzerAllSteps::center_beampipe_y, 2) ) ;
	
	double RECOLxy_interactionVertex_beampipeCenter_error = 1/RECOLxy_interactionVertex_beampipeCenter*sqrt( pow(RECOAntiSInteractionVertex.X(),2)*pow(RECO_S->vertexCovariance(0,0),2 ) + pow(RECOAntiSInteractionVertex.Y(),2)*pow(RECO_S->vertexCovariance(1,1),2 ) );
	double RECOLxy_interactionVertex = AnalyzerAllSteps::lxy(beamspot,RECOAntiSInteractionVertex);

	//if running on MC and the RECO_S charge is negative and the RECO_S has an interaction vertex which is far enough in lxy, check if this is a real AntiS by looking at the difference in lxyz between the RECO and the GEN antiS. Save this deltaLInteractionVertexAntiSmin in the tree, like this later I can easily select in the tree on signal antiS and background antiS
        double deltaLInteractionVertexAntiSmin = 999.;
        double deltaRAntiSmin = 999.;
	int bestMatchingAntiS = -1;
	if(!m_runningOnData && RECO_S->charge() == -1  && RECOLxy_interactionVertex >= AnalyzerAllSteps::MinLxyCut){
		if(h_genParticles.isValid()){
			for(unsigned int i = 0; i < h_genParticles->size(); ++i){//loop all genparticlesPlusGEANT and only for the ones with the correct pdgId check if this RECO antiS is matching a GEN particle and is thus not a fake antiS
				if(h_genParticles->at(i).pdgId() != AnalyzerAllSteps::pdgIdAntiS) continue;
				if(h_genParticles->at(i).numberOfDaughters() != 2) continue;
				//have to use the vertex of the daughter Ks (at GEN level) as the interaction vertex of the AntiS and compare it to the vertex of the RECO antiS which is the annihilation vertex
				double deltaLInteractionVertexAntiS = sqrt( pow(h_genParticles->at(i).daughter(0)->vx() - RECO_S->vx(),2) + pow(h_genParticles->at(i).daughter(0)->vy() - RECO_S->vy(),2) +  pow(h_genParticles->at(i).daughter(0)->vz() - RECO_S->vz(),2)  );
				if(deltaLInteractionVertexAntiS < deltaLInteractionVertexAntiSmin){
					deltaLInteractionVertexAntiSmin = deltaLInteractionVertexAntiS;
					deltaRAntiSmin = AnalyzerAllSteps::deltaR(h_genParticles->at(i).phi(), h_genParticles->at(i).eta(),RECO_S->phi(),RECO_S->eta());
					bestMatchingAntiS = i;
				}


		       }
		}
	}
	//the weighting factor for events will depend on their pathlength through the beampipe
	double event_weighting_factor = AnalyzerAllSteps::EventWeightingFactor(RECO_S->theta()); 
	double event_weighting_factorPU = 1.; 
	//you only need to calculate a reweighing parameter for the PU and z location if you are running on MC
	std::cout << "just before event_weighting_factorPU: " << ngoodPVsPOG << ", " << AnalyzerAllSteps::v_mapPU.size() << " , " << bestMatchingAntiS << std::endl;
        if(ngoodPVsPOG < AnalyzerAllSteps::v_mapPU.size() && bestMatchingAntiS > -1) {
		std::cout << "event_weighting_factorPU for signal" << event_weighting_factorPU << std::endl;
		event_weighting_factorPU = AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[ngoodPVsPOG],h_genParticles->at(bestMatchingAntiS).vz());
	}
	else if(ngoodPVsPOG < AnalyzerAllSteps::v_mapPU.size()){ //but if the MC does not contain any antiS you have to reweigh on the 'event', so pick a random PVz location to reweigh on
		event_weighting_factorPU = AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[ngoodPVsPOG],randomPVz);
		event_weighting_factorPU = event_weighting_factorPU * ngoodPVsPOG / 18.479;
		std::cout << "event_weighting_factorPU for BKG" << event_weighting_factorPU << std::endl;
	}

	nTotalRECOS++;
	//calculate some kinematic variables for the RECO AntiS
	TVector3 RECOAntiSMomentumVertex(RECO_S->px(),RECO_S->py(),RECO_S->pz());
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


	//for the granddaughters: problem is the Ks and Lambda as daughter of the AntiS do not have daughters, so I need to go through the reconstructed Ks and Lambda collection and find the best matching ones
	//for the Lambda
	const reco::Candidate* Lambda_fromAntiS = RECO_S->daughter(0);
	const reco::Candidate* Ks_fromAntiS = RECO_S->daughter(1);

	double deltaRMinLambda = 999;
	double bestMatchingLambda = 0;
	for(unsigned int i_l = 0; i_l <  h_V0L->size(); i_l++){
		double deltaPhi = reco::deltaPhi(h_V0L->at(i_l).phi() , Lambda_fromAntiS->phi());
		double deltaEta = h_V0L->at(i_l).eta() - Lambda_fromAntiS->eta();
		double deltaR = sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta);
		if(deltaR < deltaRMinLambda){deltaRMinLambda=deltaR;bestMatchingLambda=i_l;}
	}

	const reco::VertexCompositeCandidate Lambda = h_V0L->at(bestMatchingLambda); 

	double deltaRMinKs = 999;
	double bestMatchingKs = 0;
	for(unsigned int i_k = 0; i_k <  h_V0Ks->size(); i_k++){
		double deltaPhi = reco::deltaPhi(h_V0Ks->at(i_k).phi() , Ks_fromAntiS->phi());
		double deltaEta = h_V0Ks->at(i_k).eta() - Ks_fromAntiS->eta();
		double deltaR = sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta);
		if(deltaR < deltaRMinKs){deltaRMinKs=deltaR;bestMatchingKs=i_k;}
	}
	const reco::VertexCompositeCandidate Ks = h_V0Ks->at(bestMatchingKs);
	 
	//track1
	double RECO_Lambda_daughter0_charge = Lambda.daughter(0)->charge();
	double RECO_Lambda_daughter0_pt = Lambda.daughter(0)->pt();
	double RECO_Lambda_daughter0_pz = Lambda.daughter(0)->pz();
	TVector3 RECO_Lambda_Daughter0Momentum( Lambda.daughter(0)->px(), Lambda.daughter(0)->py(), Lambda.daughter(0)->pz() );	
	TVector3 RECO_Lambda_Daughter0Vertex( Lambda.daughter(0)->vx(), Lambda.daughter(0)->vy(), Lambda.daughter(0)->vz() );	
	double RECO_Lambda_daughter0_dxy_beamspot = AnalyzerAllSteps::dxy_signed_line_point(RECO_Lambda_Daughter0Vertex, RECO_Lambda_Daughter0Momentum, beamspot);
	double RECO_Lambda_daughter0_dz_beamspot = AnalyzerAllSteps::dz_line_point(RECO_Lambda_Daughter0Vertex, RECO_Lambda_Daughter0Momentum, beamspot);
	//track2	
	double RECO_Lambda_daughter1_charge = Lambda.daughter(1)->charge();
	double RECO_Lambda_daughter1_pt = Lambda.daughter(1)->pt();
	double RECO_Lambda_daughter1_pz = Lambda.daughter(1)->pz();
	TVector3 RECO_Lambda_Daughter1Momentum( Lambda.daughter(1)->px(), Lambda.daughter(1)->py(), Lambda.daughter(1)->pz() );	
	TVector3 RECO_Lambda_Daughter1Vertex( Lambda.daughter(1)->vx(), Lambda.daughter(1)->vy(), Lambda.daughter(1)->vz() );	
	double RECO_Lambda_daughter1_dxy_beamspot = AnalyzerAllSteps::dxy_signed_line_point(RECO_Lambda_Daughter1Vertex, RECO_Lambda_Daughter1Momentum, beamspot);
	double RECO_Lambda_daughter1_dz_beamspot = AnalyzerAllSteps::dz_line_point(RECO_Lambda_Daughter1Vertex, RECO_Lambda_Daughter1Momentum, beamspot);

	//for the Ks
	//track1
	double RECO_Ks_daughter0_charge = Ks.daughter(0)->charge();
	double RECO_Ks_daughter0_pt = Ks.daughter(0)->pt();
	double RECO_Ks_daughter0_pz = Ks.daughter(0)->pz();
	TVector3 RECO_Ks_Daughter0Momentum( Ks.daughter(0)->px(), Ks.daughter(0)->py(), Ks.daughter(0)->pz() );	
	TVector3 RECO_Ks_Daughter0Vertex( Ks.daughter(0)->vx(), Ks.daughter(0)->vy(), Ks.daughter(0)->vz() );	
	double RECO_Ks_daughter0_dxy_beamspot = AnalyzerAllSteps::dxy_signed_line_point(RECO_Ks_Daughter0Vertex, RECO_Ks_Daughter0Momentum, beamspot);
	double RECO_Ks_daughter0_dz_beamspot = AnalyzerAllSteps::dz_line_point(RECO_Ks_Daughter0Vertex, RECO_Ks_Daughter0Momentum, beamspot);
	//track2	
	double RECO_Ks_daughter1_charge = Ks.daughter(1)->charge();
	double RECO_Ks_daughter1_pt = Ks.daughter(1)->pt();
	double RECO_Ks_daughter1_pz = Ks.daughter(1)->pz();
	TVector3 RECO_Ks_Daughter1Momentum( Ks.daughter(1)->px(), Ks.daughter(1)->py(), Ks.daughter(1)->pz() );	
	TVector3 RECO_Ks_Daughter1Vertex( Ks.daughter(1)->vx(), Ks.daughter(1)->vy(), Ks.daughter(1)->vz() );	
	double RECO_Ks_daughter1_dxy_beamspot = AnalyzerAllSteps::dxy_signed_line_point(RECO_Ks_Daughter1Vertex, RECO_Ks_Daughter1Momentum, beamspot);
	double RECO_Ks_daughter1_dz_beamspot = AnalyzerAllSteps::dz_line_point(RECO_Ks_Daughter1Vertex, RECO_Ks_Daughter1Momentum, beamspot);


	//if the RECO S particle fails one of the below cuts than don't fill the tree. These already cut the majority of the background, so the background trees will be much smaller, which is nice for computational reasons
	//if(RECOLxy_interactionVertex < AnalyzerAllSteps::MinLxyCut || RECOErrorLxy_interactionVertex > AnalyzerAllSteps::MaxErrorLxyCut)return;
	if(RECOLxy_interactionVertex_beampipeCenter < AnalyzerAllSteps::MinLxyCut  )return;


	nSavedRECOS++;

	Init();	

	_S_charge.push_back(RECO_S->charge());
	_S_deltaLInteractionVertexAntiSmin.push_back(deltaLInteractionVertexAntiSmin);
	_S_deltaRAntiSmin.push_back(deltaRAntiSmin);
	_S_deltaRKsAntiSmin.push_back(deltaRMinKs);
	_S_deltaRLambdaAntiSmin.push_back(deltaRMinLambda);

	_S_lxy_interaction_vertex.push_back(RECOLxy_interactionVertex);
	_S_lxy_interaction_vertex_beampipeCenter.push_back(RECOLxy_interactionVertex_beampipeCenter);
	_S_error_lxy_interaction_vertex.push_back(RECOErrorLxy_interactionVertex);
	_S_error_lxy_interaction_vertex_beampipeCenter.push_back(RECOLxy_interactionVertex_beampipeCenter_error);
	_Ks_lxy_decay_vertex.push_back(RECOLxy_Ks);
	_Lambda_lxy_decay_vertex.push_back(RECOLxy_Lambda);
	_S_mass.push_back(RECO_Smass);
	_S_chi2_ndof.push_back(RECO_S->vertexNormalizedChi2());
	_S_event_weighting_factor.push_back(event_weighting_factor);
	_S_event_weighting_factorPU.push_back(event_weighting_factorPU);
	_S_event_weighting_factorALL.push_back(event_weighting_factor*event_weighting_factorPU);

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

	_S_vz_interaction_vertex.push_back(RECO_S->vz()-beamspot.Z());
        _Lambda_vz_decay_vertex.push_back(RECOAntiSDaug0Vertex.Z()-beamspot.Z());
        _Ks_vz_decay_vertex.push_back(RECOAntiSDaug1Vertex.Z()-beamspot.Z());

	_S_vx.push_back(RECO_S->vx());	
	_S_vy.push_back(RECO_S->vy());	
	_S_vz.push_back(RECO_S->vz()-beamspot.Z());	

	_Lambda_mass.push_back(RECO_S->daughter(0)->mass());
	_Ks_mass.push_back(RECO_S->daughter(1)->mass());

	_RECO_Lambda_daughter0_charge.push_back(RECO_Lambda_daughter0_charge);
	_RECO_Lambda_daughter0_pt.push_back(RECO_Lambda_daughter0_pt);
	_RECO_Lambda_daughter0_pz.push_back(RECO_Lambda_daughter0_pz);
	_RECO_Lambda_daughter0_dxy_beamspot.push_back(RECO_Lambda_daughter0_dxy_beamspot);
	_RECO_Lambda_daughter0_dz_beamspot.push_back(RECO_Lambda_daughter0_dz_beamspot);

	_RECO_Lambda_daughter1_charge.push_back(RECO_Lambda_daughter1_charge);
	_RECO_Lambda_daughter1_pt.push_back(RECO_Lambda_daughter1_pt);
	_RECO_Lambda_daughter1_pz.push_back(RECO_Lambda_daughter1_pz);
	_RECO_Lambda_daughter1_dxy_beamspot.push_back(RECO_Lambda_daughter1_dxy_beamspot);
	_RECO_Lambda_daughter1_dz_beamspot.push_back(RECO_Lambda_daughter1_dz_beamspot);

	_RECO_Ks_daughter0_charge.push_back(RECO_Ks_daughter0_charge);
	_RECO_Ks_daughter0_pt.push_back(RECO_Ks_daughter0_pt);
	_RECO_Ks_daughter0_pz.push_back(RECO_Ks_daughter0_pz);
	_RECO_Ks_daughter0_dxy_beamspot.push_back(RECO_Ks_daughter0_dxy_beamspot);
	_RECO_Ks_daughter0_dz_beamspot.push_back(RECO_Ks_daughter0_dz_beamspot);

	_RECO_Ks_daughter1_charge.push_back(RECO_Ks_daughter1_charge);
	_RECO_Ks_daughter1_pt.push_back(RECO_Ks_daughter1_pt);
	_RECO_Ks_daughter1_pz.push_back(RECO_Ks_daughter1_pz);
	_RECO_Ks_daughter1_dxy_beamspot.push_back(RECO_Ks_daughter1_dxy_beamspot);
	_RECO_Ks_daughter1_dz_beamspot.push_back(RECO_Ks_daughter1_dz_beamspot);

  	_tree->Fill();


	//very temporary: for the SingleElectron2016 RunG I found some S which have a high BDT variable. Now try to find exactly the events where these S occur
	bool foundHighBDTS = false;
	if(_S_eta[0] == -2.99413752556) foundHighBDTS = true;
	if(_S_eta[0] == 1.80451142788) foundHighBDTS = true;
	if(_S_eta[0] == 2.85730171204) foundHighBDTS = true;
	if(_S_eta[0] == 2.88255834579) foundHighBDTS = true;
	if(_S_eta[0] == 2.62389588356)  foundHighBDTS = true;
	if(_S_eta[0] == 1.38479459286) foundHighBDTS = true;
	if(_S_eta[0] == 1.48003470898) foundHighBDTS = true;
	//if(_S_eta[0] == -2.99413752556 &&  _S_lxy_interaction_vertex[0] == 2.0698633194 && _S_vz_interaction_vertex[0] == -9.11111164093) foundHighBDTS = true;
	//if(_S_eta[0] == 1.80451142788 &&  _S_lxy_interaction_vertex[0] == 2.48056674004 && _S_vz_interaction_vertex[0] == 12.9399452209) foundHighBDTS = true;
	//if(_S_eta[0] == 2.85730171204 &&  _S_lxy_interaction_vertex[0] == 2.09732484818 && _S_vz_interaction_vertex[0] == 11.9825458527) foundHighBDTS = true;
	//if(_S_eta[0] == 2.88255834579 &&  _S_lxy_interaction_vertex[0] == 2.19579172134 && _S_vz_interaction_vertex[0] == 14.0933685303) foundHighBDTS = true;
	//if(_S_eta[0] == 2.62389588356 &&  _S_lxy_interaction_vertex[0] == 2.3096203804 && _S_vz_interaction_vertex[0] == 9.77902030945) foundHighBDTS = true;
	//if(_S_eta[0] == 1.38479459286 &&  _S_lxy_interaction_vertex[0] == 2.16217327118 && _S_vz_interaction_vertex[0] == 4.99546194077) foundHighBDTS = true;
	//if(_S_eta[0] == 1.48003470898 &&  _S_lxy_interaction_vertex[0] == 2.01532459259 && _S_vz_interaction_vertex[0] == 9.65248298645) foundHighBDTS = true;

	if(foundHighBDTS)std::cout << "found an S which matches the S with high BDT _S_eta, _S_lxy_interaction_vertex, _S_vz_interaction_vertex: " << _S_eta[0] << ", " << _S_lxy_interaction_vertex[0] << "," << _S_vz_interaction_vertex[0]  << std::endl; 


}


void FlatTreeProducerBDT::endJob()
{
}

void
FlatTreeProducerBDT::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  srand(time(NULL));
}

// ------------ method called when ending the processing of a run  ------------
void
FlatTreeProducerBDT::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
FlatTreeProducerBDT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
FlatTreeProducerBDT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatTreeProducerBDT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

FlatTreeProducerBDT::~FlatTreeProducerBDT()
{
	if(m_lookAtAntiS){
			std::cout << "The total number of GEN anti-S with all correct granddaughters: " << nGENAntiSWithCorrectGranddaughters << std::endl;
			std::cout << "The total number RECO anti-S that were found is: " << nTotalRECOS << std::endl; 
			std::cout << "The total number RECO anti-S that were saved is: " << nSavedRECOS << std::endl; 
	}
	if(!m_lookAtAntiS){
			std::cout << "The total number RECO S that were found is: " << nTotalRECOS << std::endl; 
			std::cout << "The total number RECO S that were saved is: " << nSavedRECOS << std::endl; 
	}
	std::cout << "saved/found = " << (double)nSavedRECOS/(double)nTotalRECOS << std::endl;
}

void FlatTreeProducerBDT::Init_PV()
{
	_nPV.clear();
        _nGoodPV.clear();
        _nGoodPVPOG.clear();
        _PVx.clear();
        _PVy.clear();
        _PVz.clear();
        _goodPVx.clear();
        _goodPVy.clear();
        _goodPVz.clear();
        _goodPVxPOG.clear();
        _goodPVyPOG.clear();
        _goodPVzPOG.clear();
}

void FlatTreeProducerBDT::Init_Counter()
{
	_nGENAntiSWithCorrectGranddaughters.clear();
}

void FlatTreeProducerBDT::Init()
{


    	_S_charge.clear();
	_S_deltaLInteractionVertexAntiSmin.clear();
	_S_deltaRAntiSmin.clear();
	_S_deltaRKsAntiSmin.clear();
	_S_deltaRLambdaAntiSmin.clear();

    	_S_lxy_interaction_vertex.clear();
    	_S_lxy_interaction_vertex_beampipeCenter.clear();
        _S_error_lxy_interaction_vertex.clear();
        _S_error_lxy_interaction_vertex_beampipeCenter.clear();
	_Ks_lxy_decay_vertex.clear();
	_Lambda_lxy_decay_vertex.clear();
        _S_mass.clear();
        _S_chi2_ndof.clear();
        _S_event_weighting_factor.clear();
        _S_event_weighting_factorPU.clear();
        _S_event_weighting_factorALL.clear();

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

	_Lambda_mass.clear();
	_Ks_mass.clear();

	_RECO_Lambda_daughter0_charge.clear();
	_RECO_Lambda_daughter0_pt.clear();
	_RECO_Lambda_daughter0_pz.clear();
	_RECO_Lambda_daughter0_dxy_beamspot.clear();
	_RECO_Lambda_daughter0_dz_beamspot.clear();

	_RECO_Lambda_daughter1_charge.clear();
	_RECO_Lambda_daughter1_pt.clear();
	_RECO_Lambda_daughter1_pz.clear();
	_RECO_Lambda_daughter1_dxy_beamspot.clear();
	_RECO_Lambda_daughter1_dz_beamspot.clear();

	_RECO_Ks_daughter0_charge.clear();
	_RECO_Ks_daughter0_pt.clear();
	_RECO_Ks_daughter0_pz.clear();
	_RECO_Ks_daughter0_dxy_beamspot.clear();
	_RECO_Ks_daughter0_dz_beamspot.clear();
	
	_RECO_Ks_daughter1_charge.clear();
	_RECO_Ks_daughter1_pt.clear();
	_RECO_Ks_daughter1_pz.clear();
	_RECO_Ks_daughter1_dxy_beamspot.clear();
	_RECO_Ks_daughter1_dz_beamspot.clear();

}

DEFINE_FWK_MODULE(FlatTreeProducerBDT);
