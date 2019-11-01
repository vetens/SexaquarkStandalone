#include "../interface/FlatTreeProducerGENSIM.h"
#include <typeinfo>

FlatTreeProducerGENSIM::FlatTreeProducerGENSIM(edm::ParameterSet const& pset):
  m_lookAtAntiS(pset.getUntrackedParameter<bool>("lookAtAntiS")),
  m_runningOnData(pset.getUntrackedParameter<bool>("runningOnData")),
  m_bsTag(pset.getParameter<edm::InputTag>("beamspot")),
  m_offlinePVTag(pset.getParameter<edm::InputTag>("offlinePV")),
  m_genParticlesTag_GEN(pset.getParameter<edm::InputTag>("genCollection_GEN")),
  m_genParticlesTag_SIM_GEANT(pset.getParameter<edm::InputTag>("genCollection_SIM_GEANT")),
  m_TPTag(pset.getParameter<edm::InputTag>("TrackingParticles")),

  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_offlinePVToken    (consumes<vector<reco::Vertex>>(m_offlinePVTag)),
  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM_GEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM_GEANT)),
  m_TPToken(consumes<vector<TrackingParticle> >(m_TPTag))
  
{

}


void FlatTreeProducerGENSIM::beginJob() {

    
        // Initialize when class is created
        edm::Service<TFileService> fs ;


	_treeAllAntiS = fs->make <TTree>("FlatTreeGENLevelAllAntiS","treeAllAntiS");
	_treeAllAntiS->Branch("_S_eta_all",&_S_eta_all);
	_treeAllAntiS->Branch("_S_reconstructable_all",&_S_reconstructable_all);
	_treeAllAntiS->Branch("_S_event_weighting_factor_all",&_S_event_weighting_factor_all);
	_treeAllAntiS->Branch("_S_event_weighting_factor_PU_all",&_S_event_weighting_factor_PU_all);
	_treeAllAntiS->Branch("_S_vz_creation_vertex_all",&_S_vz_creation_vertex_all);
	_treeAllAntiS->Branch("_S_nGoodPV_all",&_S_nGoodPV_all);
	_treeAllAntiS->Branch("_S_pt_all",&_S_pt_all);
	_treeAllAntiS->Branch("_S_pz_all",&_S_pz_all);

        _tree = fs->make <TTree>("FlatTreeGENLevel","tree");
        // Declare tree's branches
	_tree->Branch("_S_n_loops",&_S_n_loops);
	_tree->Branch("_S_charge",&_S_charge);
	_tree->Branch("_S_nGoodPV",&_S_nGoodPV);
	_tree->Branch("_S_event_weighting_factor",&_S_event_weighting_factor);
	_tree->Branch("_S_event_weighting_factor_PU",&_S_event_weighting_factor_PU);
	_tree->Branch("_S_lxy_interaction_vertex",&_S_lxy_interaction_vertex);
	_tree->Branch("_S_lxy_interaction_vertex_beamspot",&_S_lxy_interaction_vertex_beamspot);
	_tree->Branch("_S_lxy_interaction_vertex_beampipeCenterData",&_S_lxy_interaction_vertex_beampipeCenterData);
	_tree->Branch("_S_lxyz_interaction_vertex",&_S_lxyz_interaction_vertex);
	_tree->Branch("_S_error_lxy_interaction_vertex",&_S_error_lxy_interaction_vertex);
	_tree->Branch("_S_mass",&_S_mass);
	_tree->Branch("_S_Mt",&_S_Mt);
	_tree->Branch("_n_M",&_n_M);
	_tree->Branch("_n_p",&_n_p);
	_tree->Branch("_S_chi2_ndof",&_S_chi2_ndof);

	_tree->Branch("_S_daughters_deltaphi",&_S_daughters_deltaphi);
	_tree->Branch("_S_daughters_deltaeta",&_S_daughters_deltaeta);
	_tree->Branch("_S_daughters_openingsangle",&_S_daughters_openingsangle);
	_tree->Branch("_S_Ks_openingsangle",&_S_Ks_openingsangle);
	_tree->Branch("_S_Lambda_openingsangle",&_S_Lambda_openingsangle);
	_tree->Branch("_S_sumDaughters_openingsangle",&_S_sumDaughters_openingsangle);
	_tree->Branch("_S_sumDaughters_deltaPhi",&_S_sumDaughters_deltaPhi);
	_tree->Branch("_S_sumDaughters_deltaEta",&_S_sumDaughters_deltaEta);
	_tree->Branch("_S_sumDaughters_deltaR",&_S_sumDaughters_deltaR);
	_tree->Branch("_S_daughters_DeltaR",&_S_daughters_DeltaR);
	_tree->Branch("_S_eta",&_S_eta);
	_tree->Branch("_Ks_eta",&_Ks_eta);
	_tree->Branch("_Lambda_eta",&_Lambda_eta);

	_tree->Branch("_S_dxy",&_S_dxy);
	_tree->Branch("_Ks_dxy",&_Ks_dxy);
	_tree->Branch("_Lambda_dxy",&_Lambda_dxy);

	_tree->Branch("_S_dxy_over_lxy",&_S_dxy_over_lxy);
	_tree->Branch("_Ks_dxy_over_lxy",&_Ks_dxy_over_lxy);
	_tree->Branch("_Lambda_dxy_over_lxy",&_Lambda_dxy_over_lxy);

	_tree->Branch("_S_dz",&_S_dz);
	_tree->Branch("_Ks_dz",&_Ks_dz);
	_tree->Branch("_Lambda_dz",&_Lambda_dz);
	_tree->Branch("_S_dz_min",&_S_dz_min);
	_tree->Branch("_Ks_dz_min",&_Ks_dz_min);
	_tree->Branch("_Lambda_dz_min",&_Lambda_dz_min);

	_tree->Branch("_deltaR_sumDaughterMomenta_antiSMomentum",&_deltaR_sumDaughterMomenta_antiSMomentum);

	_tree->Branch("_Ks_openings_angle_displacement_momentum",&_Ks_openings_angle_displacement_momentum);
	_tree->Branch("_Lambda_openings_angle_displacement_momentum",&_Lambda_openings_angle_displacement_momentum);

	_tree->Branch("_S_pt",&_S_pt);
	_tree->Branch("_Ks_pt",&_Ks_pt);
	_tree->Branch("_Lambda_pt",&_Lambda_pt);

	_tree->Branch("_S_pz",&_S_pz);
	_tree->Branch("_Ks_pz",&_Ks_pz);
	_tree->Branch("_Lambda_pz",&_Lambda_pz);

	_tree->Branch("_S_vx_interaction_vertex",&_S_vx_interaction_vertex);
	_tree->Branch("_S_vy_interaction_vertex",&_S_vy_interaction_vertex);
	_tree->Branch("_S_vz_interaction_vertex",&_S_vz_interaction_vertex);

	_tree->Branch("_S_vx",&_S_vx);
	_tree->Branch("_S_vy",&_S_vy);
	_tree->Branch("_S_vz",&_S_vz);

	_tree->Branch("_GEN_Ks_daughter0_px",&_GEN_Ks_daughter0_px);
	_tree->Branch("_GEN_Ks_daughter0_py",&_GEN_Ks_daughter0_py);
	_tree->Branch("_GEN_Ks_daughter0_pz",&_GEN_Ks_daughter0_pz);
	_tree->Branch("_GEN_Ks_daughter0_pt",&_GEN_Ks_daughter0_pt);
	_tree->Branch("_GEN_Ks_daughter0_eta",&_GEN_Ks_daughter0_eta);
	_tree->Branch("_GEN_Ks_daughter0_phi",&_GEN_Ks_daughter0_phi);
	_tree->Branch("_GEN_Ks_daughter0_vx",&_GEN_Ks_daughter0_vx);
	_tree->Branch("_GEN_Ks_daughter0_vy",&_GEN_Ks_daughter0_vy);
	_tree->Branch("_GEN_Ks_daughter0_vz",&_GEN_Ks_daughter0_vz);
	_tree->Branch("_GEN_Ks_daughter0_lxy",&_GEN_Ks_daughter0_lxy);
	_tree->Branch("_GEN_Ks_daughter0_dxy",&_GEN_Ks_daughter0_dxy);
	_tree->Branch("_GEN_Ks_daughter0_dz",&_GEN_Ks_daughter0_dz);
	_tree->Branch("_GEN_Ks_daughter0_openings_angle_displacement_momentum",&_GEN_Ks_daughter0_openings_angle_displacement_momentum);

	_tree->Branch("_GEN_Ks_daughter1_px",&_GEN_Ks_daughter1_px);
	_tree->Branch("_GEN_Ks_daughter1_py",&_GEN_Ks_daughter1_py);
	_tree->Branch("_GEN_Ks_daughter1_pz",&_GEN_Ks_daughter1_pz);
	_tree->Branch("_GEN_Ks_daughter1_pt",&_GEN_Ks_daughter1_pt);
	_tree->Branch("_GEN_Ks_daughter1_eta",&_GEN_Ks_daughter1_eta);
	_tree->Branch("_GEN_Ks_daughter1_phi",&_GEN_Ks_daughter1_phi);
	_tree->Branch("_GEN_Ks_daughter1_vx",&_GEN_Ks_daughter1_vx);
	_tree->Branch("_GEN_Ks_daughter1_vy",&_GEN_Ks_daughter1_vy);
	_tree->Branch("_GEN_Ks_daughter1_vz",&_GEN_Ks_daughter1_vz);
	_tree->Branch("_GEN_Ks_daughter1_lxy",&_GEN_Ks_daughter1_lxy);
	_tree->Branch("_GEN_Ks_daughter1_dxy",&_GEN_Ks_daughter1_dxy);
	_tree->Branch("_GEN_Ks_daughter1_dz",&_GEN_Ks_daughter1_dz);
	_tree->Branch("_GEN_Ks_daughter1_openings_angle_displacement_momentum",&_GEN_Ks_daughter1_openings_angle_displacement_momentum);

	_tree->Branch("_GEN_AntiLambda_AntiProton_px",&_GEN_AntiLambda_AntiProton_px);
	_tree->Branch("_GEN_AntiLambda_AntiProton_py",&_GEN_AntiLambda_AntiProton_py);
	_tree->Branch("_GEN_AntiLambda_AntiProton_pz",&_GEN_AntiLambda_AntiProton_pz);
	_tree->Branch("_GEN_AntiLambda_AntiProton_pt",&_GEN_AntiLambda_AntiProton_pt);
	_tree->Branch("_GEN_AntiLambda_AntiProton_eta",&_GEN_AntiLambda_AntiProton_eta);
	_tree->Branch("_GEN_AntiLambda_AntiProton_phi",&_GEN_AntiLambda_AntiProton_phi);
	_tree->Branch("_GEN_AntiLambda_AntiProton_vx",&_GEN_AntiLambda_AntiProton_vx);
	_tree->Branch("_GEN_AntiLambda_AntiProton_vy",&_GEN_AntiLambda_AntiProton_vy);
	_tree->Branch("_GEN_AntiLambda_AntiProton_vz",&_GEN_AntiLambda_AntiProton_vz);
	_tree->Branch("_GEN_AntiLambda_AntiProton_lxy",&_GEN_AntiLambda_AntiProton_lxy);
	_tree->Branch("_GEN_AntiLambda_AntiProton_dxy",&_GEN_AntiLambda_AntiProton_dxy);
	_tree->Branch("_GEN_AntiLambda_AntiProton_dz",&_GEN_AntiLambda_AntiProton_dz);
	_tree->Branch("_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum",&_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum);

	_tree->Branch("_GEN_AntiLambda_Pion_px",&_GEN_AntiLambda_Pion_px);
	_tree->Branch("_GEN_AntiLambda_Pion_py",&_GEN_AntiLambda_Pion_py);
	_tree->Branch("_GEN_AntiLambda_Pion_pz",&_GEN_AntiLambda_Pion_pz);
	_tree->Branch("_GEN_AntiLambda_Pion_pt",&_GEN_AntiLambda_Pion_pt);
	_tree->Branch("_GEN_AntiLambda_Pion_eta",&_GEN_AntiLambda_Pion_eta);
	_tree->Branch("_GEN_AntiLambda_Pion_phi",&_GEN_AntiLambda_Pion_phi);
	_tree->Branch("_GEN_AntiLambda_Pion_vx",&_GEN_AntiLambda_Pion_vx);
	_tree->Branch("_GEN_AntiLambda_Pion_vy",&_GEN_AntiLambda_Pion_vy);
	_tree->Branch("_GEN_AntiLambda_Pion_vz",&_GEN_AntiLambda_Pion_vz);
	_tree->Branch("_GEN_AntiLambda_Pion_lxy",&_GEN_AntiLambda_Pion_lxy);
	_tree->Branch("_GEN_AntiLambda_Pion_dxy",&_GEN_AntiLambda_Pion_dxy);
	_tree->Branch("_GEN_AntiLambda_Pion_dz",&_GEN_AntiLambda_Pion_dz);
	_tree->Branch("_GEN_AntiLambda_Pion_openings_angle_displacement_momentum",&_GEN_AntiLambda_Pion_openings_angle_displacement_momentum);

  	_tree->Branch("_GEN_Ks_daughter0_numberOfTrackerLayers",&_GEN_Ks_daughter0_numberOfTrackerLayers);
        _tree->Branch("_GEN_Ks_daughter1_numberOfTrackerLayers",&_GEN_Ks_daughter1_numberOfTrackerLayers);
        _tree->Branch("_GEN_AntiLambda_AntiProton_numberOfTrackerLayers",&_GEN_AntiLambda_AntiProton_numberOfTrackerLayers);
        _tree->Branch("_GEN_AntiLambda_Pion_numberOfTrackerLayers",&_GEN_AntiLambda_Pion_numberOfTrackerLayers);

        _tree->Branch("_GEN_Ks_daughter0_numberOfTrackerHits",&_GEN_Ks_daughter0_numberOfTrackerHits);
        _tree->Branch("_GEN_Ks_daughter1_numberOfTrackerHits",&_GEN_Ks_daughter1_numberOfTrackerHits);
        _tree->Branch("_GEN_AntiLambda_AntiProton_numberOfTrackerHits",&_GEN_AntiLambda_AntiProton_numberOfTrackerHits);
        _tree->Branch("_GEN_AntiLambda_Pion_numberOfTrackerHits",&_GEN_AntiLambda_Pion_numberOfTrackerHits);



}

void FlatTreeProducerGENSIM::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

 

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

  //trackingparticle collection
  edm::Handle<TrackingParticleCollection>  h_TP ;
  iEvent.getByToken(m_TPToken,h_TP);


  //beamspot
  TVector3 beamspot(0,0,0);
  TVector3 beamspotVariance(0,0,0);
  if(h_bs.isValid()){  
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));
	std::cout << "beamspot location: "<< sqrt(h_bs->x0()*h_bs->x0()+h_bs->y0()*h_bs->y0()) << " , " << h_bs->z0() << std::endl;
  }
  else{
	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!beamspot collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  }

  unsigned int nGoodPV = 0;
  if(h_offlinePV.isValid()){
        for(unsigned int i = 0; i < h_offlinePV->size(); i++){
                if(h_offlinePV->at(i).isValid() && h_offlinePV->at(i).tracksSize() >= 4){
                        nGoodPV++;
                }
        }
  }

  //loop over the gen particles, check for this antiS if there are any antiS with the same eta, so duplicates
  //save the number of duplicates in a vector of vectors. Each vector has  as a first element the eta of the antiS and 2nd element the # of antiS with this eta. 
  vector<vector<float>> v_antiS_momenta_and_itt; 
  for(unsigned int i = 0; i < h_genParticles->size(); ++i){
	const reco::Candidate * genParticle = &h_genParticles->at(i);

  	if(genParticle->pdgId() != AnalyzerAllSteps::pdgIdAntiS) continue;
   	nTotalGENS++;	
	/*if(v_antiS_momenta_and_itt.size() == 0){//to get the thing started save the first antiS you encounter to the vector and continue
		vector<float> dummyVec;
		dummyVec.push_back(genParticle->eta());
		dummyVec.push_back(1.);
		v_antiS_momenta_and_itt.push_back(dummyVec);
		continue;
	}*/
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
		nTotalUniqueGenS++;
		double weight_PU = 0.;
		if(nGoodPV < AnalyzerAllSteps::v_mapPU.size()) weight_PU = AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[nGoodPV],genParticle->vz());
		nTotalUniqueGenS_weighted = nTotalUniqueGenS_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;
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

  //first find the GEN particles which are proper antiS
  if(!m_runningOnData && m_lookAtAntiS){
	  if(h_genParticles.isValid()){
	      //double etaPrevious = 999;
	      //double eventWeighingPrevious = 999;
	      //bool reconstructablePrevious = false;
	      //vector of vectors containing for each of the unique antiS a vector and this vector has as first element the eta of the antiS and as second entry if it is reconstructable
	      std::vector<std::vector<double>> v_antiS_eta_reconstructable; 
	      for(unsigned int i = 0; i < h_genParticles->size(); ++i){//loop all genparticlesPlusGEANT

			const reco::Candidate * genParticle = &h_genParticles->at(i);
			if(genParticle->pdgId() != AnalyzerAllSteps::pdgIdAntiS) continue;

			//check if you already have an entry for this antiS in v_antiS_eta_reconstructable
			bool newAntiS = true;
			int antiS_it = -1;
			for(unsigned int j = 0; j < v_antiS_eta_reconstructable.size(); j++){//check for all the entries in v_antiS_eta_reconstructable if there is an eta match, if there is an eta match this antiS has already been encountered
				if(v_antiS_eta_reconstructable[j][0] == genParticle->eta() ){ newAntiS = false;antiS_it = j;}
			}
			//if there was no eta match newAntiS is still true and I should make an new entry in the vector for this guy
			if(newAntiS){

				std::vector<double> newAntiSEntry; 
				newAntiSEntry.push_back(genParticle->eta()); 
				newAntiSEntry.push_back(0);//take it as not reconstructable at first 
				newAntiSEntry.push_back(AnalyzerAllSteps::EventWeightingFactor(genParticle->theta()));
				newAntiSEntry.push_back(genParticle->vz());//save where it got produced.
				newAntiSEntry.push_back(genParticle->pt());
				newAntiSEntry.push_back(genParticle->pz());

				v_antiS_eta_reconstructable.push_back(newAntiSEntry);
				antiS_it = v_antiS_eta_reconstructable.size()-1;
			}
			
			

			bool AntiSReconstructable = false;

			double weight_PU = 0.;
                	if(nGoodPV < AnalyzerAllSteps::v_mapPU.size()) weight_PU = AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[nGoodPV],genParticle->vz());
			//check if this is a reconstructable antiS, so should have 2 daughters of correct type, each daughter should have 2 daughters with the correct type
			//the below implicitely neglects the duplitcate antiS due to looping, because only 1 of the duplicates will interact and give daughters
			if(genParticle->numberOfDaughters()==2){

				nTotalGiving2DaughtersGENS++;
				nTotalGiving2DaughtersGENS_weighted = nTotalGiving2DaughtersGENS_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;

				int daughterParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle);//returns 3 if daughters from antiS are Ks and AntiLambda
 
				if(genParticle->daughter(0)->numberOfDaughters()==2 && genParticle->daughter(1)->numberOfDaughters()==2 && daughterParticlesTypes == 3){

					nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS++;
					nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS_weighted = nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;
					int graddaughters0ParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle->daughter(0));
					int graddaughters1ParticlesTypes = AnalyzerAllSteps::getDaughterParticlesTypes(genParticle->daughter(1));
					if((graddaughters0ParticlesTypes == 1 && graddaughters1ParticlesTypes == 2) || (graddaughters1ParticlesTypes == 1 && graddaughters0ParticlesTypes == 2)){
						std::cout << "AntiS with all CORRECT: correct types and numbers of daughters and granddaughters and eta " << genParticle->eta() << std::endl;
						nTotalCorrectGENS++;
						nTotalCorrectGENS_weighted = nTotalCorrectGENS_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;
						if(genParticle->eta()>0) {nTotalGENSPosEta++; nTotalGENSPosEta_weighted = nTotalGENSPosEta_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;}
						if(genParticle->eta()<0) {nTotalGENSNegEta++; nTotalGENSNegEta_weighted = nTotalGENSNegEta_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;}
						AntiSReconstructable = FillBranchesGENAntiS(genParticle,beamspot, beamspotVariance, v_antiS_momenta_and_itt,  h_TP,nGoodPV);
						if(AntiSReconstructable)nTotalCorrectGENS_Reconstructable_weighted = nTotalCorrectGENS_Reconstructable_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;
						
					}

				}

				
			}

			//now save for this antiS in v_antiS_eta_reconstructable the or of the reconstructability flag which was already in this vector and the the current AntiSReconstructable, like that you will check for any of the duplicates if it was reconstructable
			v_antiS_eta_reconstructable[antiS_it][1] = (int)v_antiS_eta_reconstructable[antiS_it][1] | AntiSReconstructable;
		
	      }//for(unsigned int i = 0; i < h_genParticles->size(); ++i)

	    for(unsigned int j = 0; j < v_antiS_eta_reconstructable.size(); j++){
		std::cout << "v_antiS_eta_reconstructable: " << v_antiS_eta_reconstructable[j][1] << ", " << v_antiS_eta_reconstructable[j][0] << std::endl;

		_S_eta_all.push_back(v_antiS_eta_reconstructable[j][0]);
		_S_reconstructable_all.push_back(v_antiS_eta_reconstructable[j][1]);
		_S_event_weighting_factor_all.push_back(v_antiS_eta_reconstructable[j][2]);
		if(nGoodPV < AnalyzerAllSteps::v_mapPU.size()) _S_event_weighting_factor_PU_all.push_back(AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[nGoodPV],v_antiS_eta_reconstructable[j][3]));
        	else _S_event_weighting_factor_PU_all.push_back(0.);
		_S_vz_creation_vertex_all.push_back(v_antiS_eta_reconstructable[j][3]);
		_S_pt_all.push_back(v_antiS_eta_reconstructable[j][4]);
		_S_pz_all.push_back(v_antiS_eta_reconstructable[j][5]);
		_S_nGoodPV_all.push_back(nGoodPV);

		_treeAllAntiS->Fill();

		_S_eta_all.clear();
		_S_reconstructable_all.clear();
		_S_event_weighting_factor_all.clear();
		_S_event_weighting_factor_PU_all.clear();
		_S_vz_creation_vertex_all.clear();
		_S_pt_all.clear();
		_S_pz_all.clear();
		_S_nGoodPV_all.clear();

		if(v_antiS_eta_reconstructable[j][1]){
			if(v_antiS_eta_reconstructable[j][0]>0)nTotalRecoconstructableGENS_posEta++;
			if(v_antiS_eta_reconstructable[j][0]<0)nTotalRecoconstructableGENS_negEta++;
		}

	    }

	  }//if(h_genParticles.isValid())
	else{
		std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!genparticle collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	}

  }


 } //end of analyzer



bool FlatTreeProducerGENSIM::FillBranchesGENAntiS(const reco::Candidate  * genParticle, TVector3 beamspot, TVector3 beamspotVariance, vector<vector<float>> v_antiS_momenta_and_itt, edm::Handle<TrackingParticleCollection>  h_TP, unsigned int nGoodPV){
  

	//Find how many duplicates there are for this AntiS
	unsigned int itDuplicateAntiS = -1;
	for(unsigned int j = 0; j < v_antiS_momenta_and_itt.size();j++){
		if(v_antiS_momenta_and_itt[j][0] == genParticle->eta()){
			itDuplicateAntiS = j;
		}
	} 

	std::cout << "number of daughters of the GEN antiS: " << genParticle->numberOfDaughters() << " pdgID daughter0 and daughter1: " <<  genParticle->daughter(0)->pdgId() << " " << genParticle->daughter(1)->pdgId()  << std::endl;
	
	//calculate some kinematic variables for the GEN AntiS
	TVector3 GENAntiSInteractionVertex(genParticle->daughter(0)->vx(),genParticle->daughter(0)->vy(),genParticle->daughter(0)->vz());//this is the interaction vertex of the antiS and the neutron.
	TVector3 GENAntiSMomentumVertex(genParticle->px(),genParticle->py(),genParticle->pz());
	TVector3 ZeroZeroZero(0.,0.,0.);
	double GENLxy_interactionVertex= AnalyzerAllSteps::lxy(ZeroZeroZero,GENAntiSInteractionVertex);
	double GENLxy_interactionVertex_beamspot = AnalyzerAllSteps::lxy(beamspot,GENAntiSInteractionVertex);
	TVector3 BeampipeCenterData(AnalyzerAllSteps::center_beampipe_x,AnalyzerAllSteps::center_beampipe_y,0);
	double GENLxy_interactionVertex_beampipeCenterData = AnalyzerAllSteps::lxy(BeampipeCenterData,GENAntiSInteractionVertex);
	double GENLxyz_interactionVertex = AnalyzerAllSteps::lxyz(beamspot,GENAntiSInteractionVertex);
	double GENDeltaPhiDaughters = reco::deltaPhi(genParticle->daughter(0)->phi(),genParticle->daughter(1)->phi());
	double GENDeltaEtaDaughters = genParticle->daughter(0)->eta()-genParticle->daughter(1)->eta();
	double GENDeltaRDaughters = pow(GENDeltaPhiDaughters*GENDeltaPhiDaughters+GENDeltaEtaDaughters*GENDeltaEtaDaughters,0.5);

	double weight_PU = 0.;
	if(nGoodPV < AnalyzerAllSteps::v_mapPU.size()) weight_PU = AnalyzerAllSteps::PUReweighingFactor(AnalyzerAllSteps::v_mapPU[nGoodPV],genParticle->vz());
	//count the number of AntiS going to fully correct granddaughters and which interact in the beampipe 
	if(GENLxy_interactionVertex < 2.26){ 
		nTotalCorrectGENSInteractingInBeampipe++; 
		nTotalCorrectGENSInteractingInBeampipe_weighted = nTotalCorrectGENSInteractingInBeampipe_weighted + AnalyzerAllSteps::EventWeightingFactor(genParticle->theta())*weight_PU;
	} 

	reco::LeafCandidate::LorentzVector n_(0,0,0,0.939565);
	double GEN_Smass = ((genParticle->daughter(0)->daughter(0)->p4()+genParticle->daughter(0)->daughter(1)->p4()+genParticle->daughter(1)->daughter(0)->p4()+genParticle->daughter(1)->daughter(1)->p4())-n_).mass();
	double GEN_Smass_trans = ((genParticle->daughter(0)->daughter(0)->p4()+genParticle->daughter(0)->daughter(1)->p4()+genParticle->daughter(1)->daughter(0)->p4()+genParticle->daughter(1)->daughter(1)->p4())-n_).Mt();

	//calculate properties of the neutron. This is interesting when running over simulation with neutrons with a Fermi momentum
	double GEN_n_invM =  ( genParticle->daughter(0)->p4() + genParticle->daughter(1)->p4() - genParticle->p4() ).mass();
	double GEN_n_p =  ( genParticle->daughter(0)->p4() + genParticle->daughter(1)->p4() - genParticle->p4() ).P();

	//the dxy of the Ks and Lambda
	TVector3 GENAntiSDaug0Momentum(genParticle->daughter(0)->px(),genParticle->daughter(0)->py(),genParticle->daughter(0)->pz());
	TVector3 GENAntiSDaug1Momentum(genParticle->daughter(1)->px(),genParticle->daughter(1)->py(),genParticle->daughter(1)->pz());
	reco::Candidate::Vector vGENAntiSMomentum(genParticle->px(),genParticle->py(),genParticle->pz());
	reco::Candidate::Vector vGENAntiSDaug0Momentum(genParticle->daughter(0)->px(),genParticle->daughter(0)->py(),genParticle->daughter(0)->pz());
	reco::Candidate::Vector vGENAntiSDaug1Momentum(genParticle->daughter(1)->px(),genParticle->daughter(1)->py(),genParticle->daughter(1)->pz());
	double GENOpeningsAngleAntiSKs = AnalyzerAllSteps::openings_angle(vGENAntiSDaug0Momentum,vGENAntiSMomentum);
	double GENOpeningsAngleAntiSLambda = AnalyzerAllSteps::openings_angle(vGENAntiSDaug1Momentum,vGENAntiSMomentum);
	//the openingsangle between the GEN AntiS and the sum of the momenta of the daughters. Normally these should be pointing exactly in the direction of the mother, but due to the momentum of the neutron this is no longer the case
	double S_sumDaughters_openingsangle = AnalyzerAllSteps::openings_angle(vGENAntiSDaug0Momentum+vGENAntiSDaug1Momentum,vGENAntiSMomentum);	
	double S_sumDaughters_deltaPhi = reco::deltaPhi((vGENAntiSDaug0Momentum+vGENAntiSDaug1Momentum).phi(),vGENAntiSMomentum.phi());	
	double S_sumDaughters_deltaEta = (vGENAntiSDaug0Momentum+vGENAntiSDaug1Momentum).eta()-vGENAntiSMomentum.eta();	
	double S_sumDaughters_deltaR =  sqrt(S_sumDaughters_deltaPhi*S_sumDaughters_deltaPhi+S_sumDaughters_deltaEta*S_sumDaughters_deltaEta);

	double GENOpeningsAngleDaughters = AnalyzerAllSteps::openings_angle(vGENAntiSDaug0Momentum,vGENAntiSDaug1Momentum);
	double GEN_dxy_daughter0 = AnalyzerAllSteps::dxy_signed_line_point(GENAntiSInteractionVertex, GENAntiSDaug0Momentum,beamspot);
	double GEN_dxy_daughter1 = AnalyzerAllSteps::dxy_signed_line_point(GENAntiSInteractionVertex, GENAntiSDaug1Momentum,beamspot);
	//the dz of the Ks and Lambda
	double GEN_dz_daughter0 = AnalyzerAllSteps::dz_line_point(GENAntiSInteractionVertex, GENAntiSDaug0Momentum,beamspot);
	double GEN_dz_daughter1 = AnalyzerAllSteps::dz_line_point(GENAntiSInteractionVertex, GENAntiSDaug1Momentum,beamspot);
	//dxy and dz of the AntiS itself
	double GEN_dxy_antiS = AnalyzerAllSteps::dxy_signed_line_point(GENAntiSInteractionVertex,GENAntiSMomentumVertex,beamspot);
	double GEN_dz_antiS = AnalyzerAllSteps::dz_line_point(GENAntiSInteractionVertex,GENAntiSMomentumVertex,beamspot);

	//now look at the deltaR between the sum of the antiS daughters and antiS momentum. If the neutron momentum is zero then this should be a delta peak at zero, but if the neutron does not have zero momentum, then the sum of the momenta of the daughters might not be pointing to the antiS
	reco::Candidate::Vector vGENAntiSSumDaug0and1Momentum = vGENAntiSDaug0Momentum + vGENAntiSDaug1Momentum;
	double deltaPhi_sumDaughterMomenta_antiSMomentum = reco::deltaPhi(vGENAntiSSumDaug0and1Momentum.Phi(), GENAntiSMomentumVertex.Phi());
	double deltaEta_sumDaughterMomenta_antiSMomentum = vGENAntiSSumDaug0and1Momentum.Eta() - GENAntiSMomentumVertex.Eta();
	double deltaR_sumDaughterMomenta_antiSMomentum = pow(deltaPhi_sumDaughterMomenta_antiSMomentum*deltaPhi_sumDaughterMomenta_antiSMomentum + deltaEta_sumDaughterMomenta_antiSMomentum*deltaEta_sumDaughterMomenta_antiSMomentum ,0.5); 

	//openingsangle between the displacement and the momentum of the Ks and Lambda: are they pointing back to the beamspot?
	reco::Candidate::Vector GENAntiSInteractionVertex_beamspot(genParticle->daughter(0)->vx()-beamspot.X(),genParticle->daughter(0)->vy()-beamspot.Y(),genParticle->daughter(0)->vz()-beamspot.Z());
	reco::Candidate::Vector GENAntiSDaug0Momentum_vector(genParticle->daughter(0)->px(),genParticle->daughter(0)->py(),genParticle->daughter(0)->pz());
	reco::Candidate::Vector GENAntiSDaug1Momentum_vector(genParticle->daughter(1)->px(),genParticle->daughter(1)->py(),genParticle->daughter(1)->pz());
	double Ks_openings_angle_displacement_momentum = AnalyzerAllSteps::openings_angle(GENAntiSInteractionVertex_beamspot,GENAntiSDaug0Momentum_vector);
	double Lambda_openings_angle_displacement_momentum = AnalyzerAllSteps::openings_angle(GENAntiSInteractionVertex_beamspot,GENAntiSDaug1Momentum_vector);

	//now look at the granddaughters:
	//For the Ks:
	double GEN_Ks_daughter0_px = genParticle->daughter(0)->daughter(0)->px();
	double GEN_Ks_daughter0_py = genParticle->daughter(0)->daughter(0)->py();
	double GEN_Ks_daughter0_pz = genParticle->daughter(0)->daughter(0)->pz();
	TVector3 GEN_Ks_daughter0_momentum(genParticle->daughter(0)->daughter(0)->px(),genParticle->daughter(0)->daughter(0)->py(),genParticle->daughter(0)->daughter(0)->pz());
	double GEN_Ks_daughter0_pt = genParticle->daughter(0)->daughter(0)->pt();
	double GEN_Ks_daughter0_eta = genParticle->daughter(0)->daughter(0)->eta();
	double GEN_Ks_daughter0_phi = genParticle->daughter(0)->daughter(0)->phi();
	double GEN_Ks_daughter0_vx = genParticle->daughter(0)->daughter(0)->vx();
	double GEN_Ks_daughter0_vy = genParticle->daughter(0)->daughter(0)->vy();
	double GEN_Ks_daughter0_vz = genParticle->daughter(0)->daughter(0)->vz();
	TVector3 GEN_Ks_daughter0_creation_vertex(genParticle->daughter(0)->daughter(0)->vx(),genParticle->daughter(0)->daughter(0)->vy(),genParticle->daughter(0)->daughter(0)->vz());
	double GEN_Ks_daughter0_lxy = AnalyzerAllSteps::lxy(beamspot,GEN_Ks_daughter0_creation_vertex);
	double GEN_Ks_daughter0_dxy = AnalyzerAllSteps::dxy_signed_line_point(GEN_Ks_daughter0_creation_vertex,GEN_Ks_daughter0_momentum,beamspot);
	double GEN_Ks_daughter0_dz = AnalyzerAllSteps::dz_line_point(GEN_Ks_daughter0_creation_vertex,GEN_Ks_daughter0_momentum,beamspot);
	reco::Candidate::Vector GEN_Ks_daughter0_creation_vertex_beamspot_vector(genParticle->daughter(0)->daughter(0)->vx()-beamspot.X(),genParticle->daughter(0)->daughter(0)->vy()-beamspot.Y(),genParticle->daughter(0)->daughter(0)->vz()-beamspot.Z());
	reco::Candidate::Vector GEN_Ks_daughter0_momentum_vertex_vector(genParticle->daughter(0)->daughter(0)->px(),genParticle->daughter(0)->daughter(0)->py(),genParticle->daughter(0)->daughter(0)->pz());
	double GEN_Ks_daughter0_openings_angle_displacement_momentum = AnalyzerAllSteps::openings_angle(GEN_Ks_daughter0_creation_vertex_beamspot_vector,GEN_Ks_daughter0_momentum_vertex_vector);

	double GEN_Ks_daughter1_px = genParticle->daughter(0)->daughter(1)->px();
	double GEN_Ks_daughter1_py = genParticle->daughter(0)->daughter(1)->py();
	double GEN_Ks_daughter1_pz = genParticle->daughter(0)->daughter(1)->pz();
	double GEN_Ks_daughter1_pt = genParticle->daughter(0)->daughter(1)->pt();
	double GEN_Ks_daughter1_eta = genParticle->daughter(0)->daughter(1)->eta();
	double GEN_Ks_daughter1_phi = genParticle->daughter(0)->daughter(1)->phi();
	TVector3 GEN_Ks_daughter1_momentum(genParticle->daughter(0)->daughter(1)->px(),genParticle->daughter(0)->daughter(1)->py(),genParticle->daughter(0)->daughter(1)->pz());
	double GEN_Ks_daughter1_vx = genParticle->daughter(0)->daughter(1)->vx();
	double GEN_Ks_daughter1_vy = genParticle->daughter(0)->daughter(1)->vy();
	double GEN_Ks_daughter1_vz = genParticle->daughter(0)->daughter(1)->vz();
	TVector3 GEN_Ks_daughter1_creation_vertex(genParticle->daughter(0)->daughter(1)->vx(),genParticle->daughter(0)->daughter(1)->vy(),genParticle->daughter(0)->daughter(1)->vz());
	double GEN_Ks_daughter1_lxy = AnalyzerAllSteps::lxy(beamspot,GEN_Ks_daughter1_creation_vertex);
	double GEN_Ks_daughter1_dxy = AnalyzerAllSteps::dxy_signed_line_point(GEN_Ks_daughter1_creation_vertex,GEN_Ks_daughter1_momentum,beamspot);
	double GEN_Ks_daughter1_dz = AnalyzerAllSteps::dz_line_point(GEN_Ks_daughter1_creation_vertex,GEN_Ks_daughter1_momentum,beamspot);
	reco::Candidate::Vector GEN_Ks_daughter1_creation_vertex_beamspot_vector(genParticle->daughter(0)->daughter(1)->vx()-beamspot.X(),genParticle->daughter(0)->daughter(1)->vy()-beamspot.Y(),genParticle->daughter(0)->daughter(1)->vz()-beamspot.Z());
	reco::Candidate::Vector GEN_Ks_daughter1_momentum_vertex_vector(genParticle->daughter(0)->daughter(1)->px(),genParticle->daughter(0)->daughter(1)->py(),genParticle->daughter(0)->daughter(1)->pz());
	double GEN_Ks_daughter1_openings_angle_displacement_momentum = AnalyzerAllSteps::openings_angle(GEN_Ks_daughter1_creation_vertex_beamspot_vector,GEN_Ks_daughter1_momentum_vertex_vector);

	//For the AntiLambda: first find the antiproton by finding the one with the highest momentum
	const reco::Candidate  * AntiLambdaAntiProton = genParticle->daughter(1)->daughter(1);	
	const reco::Candidate  * AntiLambdaPion = genParticle->daughter(1)->daughter(0);
	if(genParticle->daughter(1)->daughter(0)->p() > genParticle->daughter(1)->daughter(1)->p()){
		AntiLambdaAntiProton = genParticle->daughter(1)->daughter(0);
		AntiLambdaPion = genParticle->daughter(1)->daughter(1);
	}
	double GEN_AntiLambda_AntiProton_px = AntiLambdaAntiProton->px();
	double GEN_AntiLambda_AntiProton_py = AntiLambdaAntiProton->py();
	double GEN_AntiLambda_AntiProton_pz = AntiLambdaAntiProton->pz();
	double GEN_AntiLambda_AntiProton_pt = AntiLambdaAntiProton->pt();
	double GEN_AntiLambda_AntiProton_eta = AntiLambdaAntiProton->eta();
	double GEN_AntiLambda_AntiProton_phi = AntiLambdaAntiProton->phi();
	TVector3 GEN_AntiLambda_AntiProton_momentum(AntiLambdaAntiProton->px(),AntiLambdaAntiProton->py(),AntiLambdaAntiProton->pz());
	double GEN_AntiLambda_AntiProton_vx = AntiLambdaAntiProton->vx();
	double GEN_AntiLambda_AntiProton_vy = AntiLambdaAntiProton->vy();
	double GEN_AntiLambda_AntiProton_vz = AntiLambdaAntiProton->vz();
	TVector3 GEN_AntiLambda_AntiProton_creation_vertex(AntiLambdaAntiProton->vx(),AntiLambdaAntiProton->vy(),AntiLambdaAntiProton->vz());
	double GEN_AntiLambda_AntiProton_lxy = AnalyzerAllSteps::lxy(beamspot,GEN_AntiLambda_AntiProton_creation_vertex);
	double GEN_AntiLambda_AntiProton_dxy = AnalyzerAllSteps::dxy_signed_line_point(GEN_AntiLambda_AntiProton_creation_vertex,GEN_AntiLambda_AntiProton_momentum,beamspot);
	double GEN_AntiLambda_AntiProton_dz = AnalyzerAllSteps::dz_line_point(GEN_AntiLambda_AntiProton_creation_vertex,GEN_AntiLambda_AntiProton_momentum,beamspot);
	reco::Candidate::Vector GEN_AntiLambda_AntiProton_creation_vertex_beamspot_vector(AntiLambdaAntiProton->vx()-beamspot.X(),AntiLambdaAntiProton->vy()-beamspot.Y(),AntiLambdaAntiProton->vz()-beamspot.Z());
	reco::Candidate::Vector GEN_AntiLambda_AntiProton_momentum_vector(AntiLambdaAntiProton->px(),AntiLambdaAntiProton->py(),AntiLambdaAntiProton->pz());
	double GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum = AnalyzerAllSteps::openings_angle(  GEN_AntiLambda_AntiProton_creation_vertex_beamspot_vector,GEN_AntiLambda_AntiProton_momentum_vector);


	double GEN_AntiLambda_Pion_px = AntiLambdaPion->px();
	double GEN_AntiLambda_Pion_py = AntiLambdaPion->py();
	double GEN_AntiLambda_Pion_pz = AntiLambdaPion->pz();
	double GEN_AntiLambda_Pion_pt = AntiLambdaPion->pt();
	double GEN_AntiLambda_Pion_eta = AntiLambdaPion->eta();
	double GEN_AntiLambda_Pion_phi = AntiLambdaPion->phi();
	TVector3 GEN_AntiLambda_Pion_momentum(AntiLambdaPion->px(),AntiLambdaPion->py(),AntiLambdaPion->pz());
	double GEN_AntiLambda_Pion_vx = AntiLambdaPion->vx();
	double GEN_AntiLambda_Pion_vy = AntiLambdaPion->vy();
	double GEN_AntiLambda_Pion_vz = AntiLambdaPion->vz();
	TVector3 GEN_AntiLambda_Pion_creation_vertex(AntiLambdaPion->vx(),AntiLambdaPion->vy(),AntiLambdaPion->vz());
	double GEN_AntiLambda_Pion_lxy = AnalyzerAllSteps::lxy(beamspot,GEN_AntiLambda_Pion_creation_vertex);
	double GEN_AntiLambda_Pion_dxy = AnalyzerAllSteps::dxy_signed_line_point(GEN_AntiLambda_Pion_creation_vertex,GEN_AntiLambda_Pion_momentum,beamspot);
	double GEN_AntiLambda_Pion_dz = AnalyzerAllSteps::dz_line_point(GEN_AntiLambda_Pion_creation_vertex,GEN_AntiLambda_Pion_momentum,beamspot);
	reco::Candidate::Vector GEN_AntiLambda_Pion_creation_vertex_beamspot_vector(AntiLambdaPion->vx()-beamspot.X(),AntiLambdaPion->vy()-beamspot.Y(),AntiLambdaPion->vz()-beamspot.Z());
	reco::Candidate::Vector GEN_AntiLambda_Pion_momentum_vector(AntiLambdaPion->px(),AntiLambdaPion->py(),AntiLambdaPion->pz());
	double GEN_AntiLambda_Pion_openings_angle_displacement_momentum = AnalyzerAllSteps::openings_angle(  GEN_AntiLambda_Pion_creation_vertex_beamspot_vector,GEN_AntiLambda_Pion_momentum_vector);  


	//loop over the trackingparticles and find the one with the same px, py, pz as the granddaughters of this AntiS. For these trackingparticles it would be interesting to get the numberOfTrackerLayers() and numberOfTrackerHits() to have an idea how many layers the granddaughters cross
	int tp_Ks_daughter0_found = -1;
	int tp_Ks_daughter1_found = -1;
	int tp_AntiLambda_AntiProton_found = -1;
	int tp_AntiLambda_Pion_found = -1;
	if(h_TP.isValid()){
		for(unsigned int i = 0; i < h_TP->size(); ++i){	

			//for the Ks daughters: match GEN particle to tp by looking at the vertices, and only try to find the match if it has not been found yet
			if((float)h_TP->at(i).vx()==(float)GEN_Ks_daughter0_vx && (float)h_TP->at(i).vy()==(float)GEN_Ks_daughter0_vy && (float)h_TP->at(i).vz()==(float)GEN_Ks_daughter0_vz && tp_Ks_daughter0_found < 0) {
				tp_Ks_daughter0_found  = i;
			}//only do the 2nd daughter if the first one has been found	
			else if((float)h_TP->at(i).vx()==(float)GEN_Ks_daughter1_vx && (float)h_TP->at(i).vy()==(float)GEN_Ks_daughter1_vy && (float)h_TP->at(i).vz()==(float)GEN_Ks_daughter1_vz && tp_Ks_daughter0_found > -1 && tp_Ks_daughter1_found < 0) {
				tp_Ks_daughter1_found = i;
			}	

			//for the AntiLambda daughter: do the same as for the Ks
			if((float)h_TP->at(i).vx()==(float)GEN_AntiLambda_AntiProton_vx && (float)h_TP->at(i).vy()==(float)GEN_AntiLambda_AntiProton_vy && (float)h_TP->at(i).vz()==(float)GEN_AntiLambda_AntiProton_vz && tp_AntiLambda_AntiProton_found < 0) {
				tp_AntiLambda_AntiProton_found = i;
			}//only do the 2nd daughter if the first one has been found 
			else if((float)h_TP->at(i).vx()==(float)GEN_AntiLambda_Pion_vx && (float)h_TP->at(i).vy()==(float)GEN_AntiLambda_Pion_vy && (float)h_TP->at(i).vz()==(float)GEN_AntiLambda_Pion_vz && tp_AntiLambda_Pion_found < 0 && tp_AntiLambda_AntiProton_found > -1) {
				tp_AntiLambda_Pion_found = i;
			}	
		}
	}
	else{
		std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!trackingparticle collection is not valid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	}

	if(tp_AntiLambda_AntiProton_found > -1 && tp_AntiLambda_Pion_found > -1){
		if(h_TP->at(tp_AntiLambda_AntiProton_found).p() < h_TP->at(tp_AntiLambda_Pion_found).p()){
			int tmp = tp_AntiLambda_AntiProton_found ;
			tp_AntiLambda_AntiProton_found = tp_AntiLambda_Pion_found;
			tp_AntiLambda_Pion_found = tmp;
		} 
	}

	int Ks_daughter0_numberOfTrackerLayers = 0;
	int Ks_daughter1_numberOfTrackerLayers = 0;
	int AntiLambda_AntiProton_numberOfTrackerLayers = 0;
	int AntiLambda_Pion_numberOfTrackerLayers = 0;

	int Ks_daughter0_numberOfTrackerHits = 0;
	int Ks_daughter1_numberOfTrackerHits = 0;
	int AntiLambda_AntiProton_numberOfTrackerHits = 0;
	int AntiLambda_Pion_numberOfTrackerHits = 0;

	if(tp_Ks_daughter0_found > -1){Ks_daughter0_numberOfTrackerLayers = h_TP->at(tp_Ks_daughter0_found).numberOfTrackerLayers(); Ks_daughter0_numberOfTrackerHits = h_TP->at(tp_Ks_daughter0_found).numberOfTrackerHits();}
	if(tp_Ks_daughter1_found > -1){Ks_daughter1_numberOfTrackerLayers = h_TP->at(tp_Ks_daughter1_found).numberOfTrackerLayers(); Ks_daughter1_numberOfTrackerHits = h_TP->at(tp_Ks_daughter1_found).numberOfTrackerHits();}
	if(tp_AntiLambda_AntiProton_found > -1){AntiLambda_AntiProton_numberOfTrackerLayers = h_TP->at(tp_AntiLambda_AntiProton_found).numberOfTrackerLayers(); AntiLambda_AntiProton_numberOfTrackerHits = h_TP->at(tp_AntiLambda_AntiProton_found).numberOfTrackerHits();}
	if(tp_AntiLambda_Pion_found > -1){AntiLambda_Pion_numberOfTrackerLayers = h_TP->at(tp_AntiLambda_Pion_found).numberOfTrackerLayers(); AntiLambda_Pion_numberOfTrackerHits = h_TP->at(tp_AntiLambda_Pion_found).numberOfTrackerHits();}

	Init(); 

	_S_n_loops.push_back(v_antiS_momenta_and_itt[itDuplicateAntiS][1]);
	_S_charge.push_back(genParticle->charge());
	_S_nGoodPV.push_back(nGoodPV);
	_S_event_weighting_factor.push_back(AnalyzerAllSteps::EventWeightingFactor(genParticle->theta()));
	_S_event_weighting_factor_PU.push_back(weight_PU);

	_S_lxy_interaction_vertex.push_back(GENLxy_interactionVertex);
	_S_lxy_interaction_vertex_beamspot.push_back(GENLxy_interactionVertex_beamspot);
	_S_lxy_interaction_vertex_beampipeCenterData.push_back(GENLxy_interactionVertex_beampipeCenterData);
	_S_lxyz_interaction_vertex.push_back(GENLxyz_interactionVertex);
	_S_mass.push_back(GEN_Smass);
	_S_Mt.push_back(GEN_Smass_trans);

	_n_M.push_back(GEN_n_invM);
	_n_p.push_back(GEN_n_p);

	_S_daughters_deltaphi.push_back(GENDeltaPhiDaughters);
	_S_daughters_deltaeta.push_back(GENDeltaEtaDaughters);
	_S_daughters_openingsangle.push_back(GENOpeningsAngleDaughters);
	_S_Ks_openingsangle.push_back(GENOpeningsAngleAntiSKs);
	_S_Lambda_openingsangle.push_back(GENOpeningsAngleAntiSLambda);
	_S_sumDaughters_openingsangle.push_back(S_sumDaughters_openingsangle);
	_S_sumDaughters_deltaPhi.push_back(S_sumDaughters_deltaPhi);
	_S_sumDaughters_deltaEta.push_back(S_sumDaughters_deltaEta);
	_S_sumDaughters_deltaR.push_back(S_sumDaughters_deltaR);
	_S_daughters_DeltaR.push_back(GENDeltaRDaughters);
	_S_eta.push_back(genParticle->eta());
	_Ks_eta.push_back(genParticle->daughter(0)->eta());
	_Lambda_eta.push_back(genParticle->daughter(1)->eta());

	_S_dxy.push_back(GEN_dxy_antiS);
	_Ks_dxy.push_back(GEN_dxy_daughter0);
	_Lambda_dxy.push_back(GEN_dxy_daughter1);

	_S_dxy_over_lxy.push_back(GEN_dxy_antiS/GENLxy_interactionVertex);
	_Ks_dxy_over_lxy.push_back(GEN_dxy_daughter0/GENLxy_interactionVertex);
	_Lambda_dxy_over_lxy.push_back(GEN_dxy_daughter1/GENLxy_interactionVertex);

	_S_dz.push_back(GEN_dz_antiS);
	_Ks_dz.push_back(GEN_dz_daughter0);
	_Lambda_dz.push_back(GEN_dz_daughter1);

	_deltaR_sumDaughterMomenta_antiSMomentum.push_back(deltaR_sumDaughterMomenta_antiSMomentum);

	_Ks_openings_angle_displacement_momentum.push_back(Ks_openings_angle_displacement_momentum);
	_Lambda_openings_angle_displacement_momentum.push_back(Lambda_openings_angle_displacement_momentum);
	

	_S_pt.push_back(genParticle->pt());
	_Ks_pt.push_back(genParticle->daughter(0)->pt());
	_Lambda_pt.push_back(genParticle->daughter(1)->pt());
	
	_S_pz.push_back(genParticle->pz());
	_Ks_pz.push_back(genParticle->daughter(0)->pz());
	_Lambda_pz.push_back(genParticle->daughter(1)->pz());

	_S_vx_interaction_vertex.push_back(GENAntiSInteractionVertex.X());
	_S_vy_interaction_vertex.push_back(GENAntiSInteractionVertex.Y());
	_S_vz_interaction_vertex.push_back(GENAntiSInteractionVertex.Z());

	_S_vx.push_back(genParticle->vx());	
	_S_vy.push_back(genParticle->vy());	
	_S_vz.push_back(genParticle->vz());	

	_GEN_Ks_daughter0_px.push_back(GEN_Ks_daughter0_px); 
	_GEN_Ks_daughter0_py.push_back(GEN_Ks_daughter0_py); 
	_GEN_Ks_daughter0_pz.push_back(GEN_Ks_daughter0_pz); 
	_GEN_Ks_daughter0_pt.push_back(GEN_Ks_daughter0_pt); 
	_GEN_Ks_daughter0_eta.push_back(GEN_Ks_daughter0_eta); 
	_GEN_Ks_daughter0_phi.push_back(GEN_Ks_daughter0_phi); 
	_GEN_Ks_daughter0_vx.push_back(GEN_Ks_daughter0_vx); 
	_GEN_Ks_daughter0_vy.push_back(GEN_Ks_daughter0_vy); 
	_GEN_Ks_daughter0_vz.push_back(GEN_Ks_daughter0_vz); 
	_GEN_Ks_daughter0_lxy.push_back(GEN_Ks_daughter0_lxy);
	_GEN_Ks_daughter0_dxy.push_back(GEN_Ks_daughter0_dxy);
	_GEN_Ks_daughter0_dz.push_back(GEN_Ks_daughter0_dz);
	_GEN_Ks_daughter0_openings_angle_displacement_momentum.push_back(GEN_Ks_daughter0_openings_angle_displacement_momentum);

	_GEN_Ks_daughter1_px.push_back(GEN_Ks_daughter1_px); 
	_GEN_Ks_daughter1_py.push_back(GEN_Ks_daughter1_py); 
	_GEN_Ks_daughter1_pz.push_back(GEN_Ks_daughter1_pz); 
	_GEN_Ks_daughter1_pt.push_back(GEN_Ks_daughter1_pt); 
	_GEN_Ks_daughter1_eta.push_back(GEN_Ks_daughter1_eta); 
	_GEN_Ks_daughter1_phi.push_back(GEN_Ks_daughter1_phi); 
	_GEN_Ks_daughter1_vx.push_back(GEN_Ks_daughter1_vx); 
	_GEN_Ks_daughter1_vy.push_back(GEN_Ks_daughter1_vy); 
	_GEN_Ks_daughter1_vz.push_back(GEN_Ks_daughter1_vz); 
	_GEN_Ks_daughter1_lxy.push_back(GEN_Ks_daughter1_lxy);
	_GEN_Ks_daughter1_dxy.push_back(GEN_Ks_daughter1_dxy);
	_GEN_Ks_daughter1_dz.push_back(GEN_Ks_daughter1_dz);
	_GEN_Ks_daughter1_openings_angle_displacement_momentum.push_back(GEN_Ks_daughter1_openings_angle_displacement_momentum);

	_GEN_AntiLambda_AntiProton_px.push_back(GEN_AntiLambda_AntiProton_px);
	_GEN_AntiLambda_AntiProton_py.push_back(GEN_AntiLambda_AntiProton_py);
	_GEN_AntiLambda_AntiProton_pz.push_back(GEN_AntiLambda_AntiProton_pz);
	_GEN_AntiLambda_AntiProton_pt.push_back(GEN_AntiLambda_AntiProton_pt);
	_GEN_AntiLambda_AntiProton_eta.push_back(GEN_AntiLambda_AntiProton_eta);
	_GEN_AntiLambda_AntiProton_phi.push_back(GEN_AntiLambda_AntiProton_phi);
	_GEN_AntiLambda_AntiProton_vx.push_back(GEN_AntiLambda_AntiProton_vx);
	_GEN_AntiLambda_AntiProton_vy.push_back(GEN_AntiLambda_AntiProton_vy);
	_GEN_AntiLambda_AntiProton_vz.push_back(GEN_AntiLambda_AntiProton_vz);
	_GEN_AntiLambda_AntiProton_lxy.push_back(GEN_AntiLambda_AntiProton_lxy);
	_GEN_AntiLambda_AntiProton_dxy.push_back(GEN_AntiLambda_AntiProton_dxy);
	_GEN_AntiLambda_AntiProton_dz.push_back(GEN_AntiLambda_AntiProton_dz);
	_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum.push_back(GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum);

	_GEN_AntiLambda_Pion_px.push_back(GEN_AntiLambda_Pion_px);
	_GEN_AntiLambda_Pion_py.push_back(GEN_AntiLambda_Pion_py);
	_GEN_AntiLambda_Pion_pz.push_back(GEN_AntiLambda_Pion_pz);
	_GEN_AntiLambda_Pion_pt.push_back(GEN_AntiLambda_Pion_pt);
	_GEN_AntiLambda_Pion_eta.push_back(GEN_AntiLambda_Pion_eta);
	_GEN_AntiLambda_Pion_phi.push_back(GEN_AntiLambda_Pion_phi);
	_GEN_AntiLambda_Pion_vx.push_back(GEN_AntiLambda_Pion_vx);
	_GEN_AntiLambda_Pion_vy.push_back(GEN_AntiLambda_Pion_vy);
	_GEN_AntiLambda_Pion_vz.push_back(GEN_AntiLambda_Pion_vz);
	_GEN_AntiLambda_Pion_lxy.push_back(GEN_AntiLambda_Pion_lxy);
	_GEN_AntiLambda_Pion_dxy.push_back(GEN_AntiLambda_Pion_dxy);
	_GEN_AntiLambda_Pion_dz.push_back(GEN_AntiLambda_Pion_dz);
	_GEN_AntiLambda_Pion_openings_angle_displacement_momentum.push_back(GEN_AntiLambda_Pion_openings_angle_displacement_momentum);

	_GEN_Ks_daughter0_numberOfTrackerLayers.push_back(Ks_daughter0_numberOfTrackerLayers);
	_GEN_Ks_daughter1_numberOfTrackerLayers.push_back(Ks_daughter1_numberOfTrackerLayers);
	_GEN_AntiLambda_AntiProton_numberOfTrackerLayers.push_back(AntiLambda_AntiProton_numberOfTrackerLayers);
	_GEN_AntiLambda_Pion_numberOfTrackerLayers.push_back(AntiLambda_Pion_numberOfTrackerLayers);

	_GEN_Ks_daughter0_numberOfTrackerHits.push_back(Ks_daughter0_numberOfTrackerHits);
	_GEN_Ks_daughter1_numberOfTrackerHits.push_back(Ks_daughter1_numberOfTrackerHits);
	_GEN_AntiLambda_AntiProton_numberOfTrackerHits.push_back(AntiLambda_AntiProton_numberOfTrackerHits);
	_GEN_AntiLambda_Pion_numberOfTrackerHits.push_back(AntiLambda_Pion_numberOfTrackerHits);

  	_tree->Fill();

	int cutNumberOfTrackerHits = 7;
	if(	Ks_daughter0_numberOfTrackerHits >= cutNumberOfTrackerHits && 
		Ks_daughter1_numberOfTrackerHits >= cutNumberOfTrackerHits && 
		AntiLambda_AntiProton_numberOfTrackerHits >= cutNumberOfTrackerHits && 
		AntiLambda_Pion_numberOfTrackerHits >= cutNumberOfTrackerHits){
		return true;
	}
	else return false;

}


void FlatTreeProducerGENSIM::endJob()
{
}

void
FlatTreeProducerGENSIM::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void
FlatTreeProducerGENSIM::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
FlatTreeProducerGENSIM::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
FlatTreeProducerGENSIM::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatTreeProducerGENSIM::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

FlatTreeProducerGENSIM::~FlatTreeProducerGENSIM()
{
	string particle = "S";
	if(m_lookAtAntiS){
		particle = "anti-S";
	}
	

	std::cout << "The total number GEN " << particle << " that were found is (this includes duplicates): " << nTotalGENS << std::endl; 
	std::cout << "The total number of unique " << particle << ": " << nTotalUniqueGenS  << std::endl; 
	std::cout << "The total number GEN " << particle << " giving 2 correct daughters: " << nTotalGiving2DaughtersGENS << std::endl; 
	std::cout << "The total number GEN " << particle << " giving 2 correct daughters and 4 granddaughters: " << nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS << std::endl; 
	std::cout << "The total number GEN " << particle << " decaying to ALL correct granddaughters is: " << nTotalCorrectGENS << std::endl;
	std::cout << "so the branching ratio you extract for all charged decay is the division of the 2 previous lines: " << (double)nTotalCorrectGENS/(double)nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS << std::endl;

	std::cout << std::endl;

	std::cout << "The efficiency factor with NUM= number of antiS with 4 correct granddaughters and DENOM = number of unique antiS: " << (double)nTotalCorrectGENS/(double)nTotalUniqueGenS << std::endl;

	std::cout << std::endl;

	std::cout << "The total number GEN " << particle << " decaying to ALL correct particles and interacting in the beampipe (lxy < 2.26cm) is: " << nTotalCorrectGENSInteractingInBeampipe << std::endl;
	std::cout << "The total number GEN " << particle << " decaying to ALL correct particles, that were found with pos eta is: " << nTotalGENSPosEta << std::endl;
        std::cout << "The total number GEN " << particle << " decaying to ALL correct particles, that were found with neg eta is: " << nTotalGENSNegEta << std::endl; 

	std::cout << "nTotalRecoconstructableGENS_posEta " << nTotalRecoconstructableGENS_posEta << std::endl;
	std::cout << "nTotalRecoconstructableGENS_negEta " << nTotalRecoconstructableGENS_negEta << std::endl;

	std::cout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;	
	std::cout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;	
	std::cout << "-BUT THESE SHOULD BE REWEIGHTED BASED ON HOW MUCH MATERAIAL THE ANTIS PASSES AND for the PU distribution (both vz and #PU)--------------------" << std::endl;	
	std::cout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;	
	std::cout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;	

	std::cout << "The total number of unique " << particle << ": " << nTotalUniqueGenS_weighted  << std::endl; 
	std::cout << "The total number GEN " << particle << " giving 2 correct daughters: " << nTotalGiving2DaughtersGENS_weighted << std::endl; 
	std::cout << "The total number GEN " << particle << " giving 2 correct daughters and 4 granddaughters: " << nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS_weighted << std::endl; 
	std::cout << "The total number GEN " << particle << " decaying to ALL correct granddaughters is: " << nTotalCorrectGENS_weighted << std::endl;
	std::cout << "so the branching ratio you extract for all charged decay is the division of the 2 previous lines: " << (double)nTotalCorrectGENS_weighted/(double)nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS_weighted << std::endl;

	std::cout << std::endl;

	std::cout << "The efficiency factor with NUM= number of antiS with 4 correct granddaughters and DENOM = number of unique antiS: " << (double)nTotalCorrectGENS_weighted/(double)nTotalUniqueGenS_weighted << std::endl;

	std::cout << std::endl;

	std::cout << "The total number GEN " << particle << " decaying to ALL correct particles and interacting in the beampipe (lxy < 2.26cm) is: " << nTotalCorrectGENSInteractingInBeampipe_weighted << std::endl;
	std::cout << "The total number GEN " << particle << " decaying to ALL correct particles, that were found with pos eta is: " << nTotalGENSPosEta_weighted << std::endl;
        std::cout << "The total number GEN " << particle << " decaying to ALL correct particles, that were found with neg eta is: " << nTotalGENSNegEta_weighted << std::endl; 

	std::cout << "The efficiency factor with NUM = number of antiS with 4 correct granddaughters which are all reconstructable and DENOM = number of unique antiS: " << (double)nTotalCorrectGENS_Reconstructable_weighted/(double)nTotalUniqueGenS_weighted << std::endl;

}


void
FlatTreeProducerGENSIM::Init()
{


    	_S_n_loops.clear();
    	_S_charge.clear();
    	_S_nGoodPV.clear();
    	_S_event_weighting_factor.clear();
    	_S_event_weighting_factor_PU.clear();

    	_S_lxy_interaction_vertex.clear();
    	_S_lxy_interaction_vertex_beamspot.clear();
    	_S_lxy_interaction_vertex_beampipeCenterData.clear();
    	_S_lxyz_interaction_vertex.clear();
        _S_error_lxy_interaction_vertex.clear();
        _S_mass.clear();
        _S_Mt.clear();
        _S_chi2_ndof.clear();

	_n_M.clear();
	_n_p.clear();

	_S_daughters_deltaphi.clear();
	_S_daughters_deltaeta.clear();
	_S_daughters_openingsangle.clear();
	_S_Ks_openingsangle.clear();
	_S_Lambda_openingsangle.clear();
	_S_sumDaughters_openingsangle.clear();
	_S_sumDaughters_deltaPhi.clear();
	_S_sumDaughters_deltaEta.clear();
	_S_sumDaughters_deltaR.clear();
	_S_daughters_DeltaR.clear();
	_S_eta.clear();
	_Ks_eta.clear();
	_Lambda_eta.clear();

	_S_dxy.clear();
	_Ks_dxy.clear();
	_Lambda_dxy.clear();

	_S_dxy_over_lxy.clear();
	_Ks_dxy_over_lxy.clear();
	_Lambda_dxy_over_lxy.clear();

	_S_dz.clear();
	_Ks_dz.clear();
	_Lambda_dz.clear();
	_S_dz_min.clear();
	_Ks_dz_min.clear();
	_Lambda_dz_min.clear();

	_deltaR_sumDaughterMomenta_antiSMomentum.clear();

	_Ks_openings_angle_displacement_momentum.clear();
	_Lambda_openings_angle_displacement_momentum.clear();

	_S_pt.clear();
	_Ks_pt.clear();
	_Lambda_pt.clear();

	_S_pz.clear();
	_Ks_pz.clear();
	_Lambda_pz.clear();

	_S_vx_interaction_vertex.clear();
	_S_vy_interaction_vertex.clear();
	_S_vz_interaction_vertex.clear();

	_S_vx.clear();
	_S_vy.clear();
	_S_vz.clear();


	_GEN_Ks_daughter0_px.clear(); 
	_GEN_Ks_daughter0_py.clear(); 
	_GEN_Ks_daughter0_pz.clear(); 
	_GEN_Ks_daughter0_pt.clear(); 
	_GEN_Ks_daughter0_eta.clear(); 
	_GEN_Ks_daughter0_phi.clear(); 
	_GEN_Ks_daughter0_vx.clear(); 
	_GEN_Ks_daughter0_vy.clear(); 
	_GEN_Ks_daughter0_vz.clear(); 
	_GEN_Ks_daughter0_lxy.clear();
	_GEN_Ks_daughter0_dxy.clear();
	_GEN_Ks_daughter0_dz.clear();
	_GEN_Ks_daughter0_openings_angle_displacement_momentum.clear();

	_GEN_Ks_daughter1_px.clear();
	_GEN_Ks_daughter1_py.clear();
	_GEN_Ks_daughter1_pz.clear();
	_GEN_Ks_daughter1_pt.clear();
	_GEN_Ks_daughter1_eta.clear();
	_GEN_Ks_daughter1_phi.clear();
	_GEN_Ks_daughter1_vx.clear();
	_GEN_Ks_daughter1_vy.clear();
	_GEN_Ks_daughter1_vz.clear();
	_GEN_Ks_daughter1_lxy.clear();
	_GEN_Ks_daughter1_dxy.clear();
	_GEN_Ks_daughter1_dz.clear();
	_GEN_Ks_daughter1_openings_angle_displacement_momentum.clear();

	_GEN_AntiLambda_AntiProton_px.clear();
	_GEN_AntiLambda_AntiProton_py.clear();
	_GEN_AntiLambda_AntiProton_pz.clear();
	_GEN_AntiLambda_AntiProton_pt.clear();
	_GEN_AntiLambda_AntiProton_eta.clear();
	_GEN_AntiLambda_AntiProton_phi.clear();
	_GEN_AntiLambda_AntiProton_vx.clear();
	_GEN_AntiLambda_AntiProton_vy.clear();
	_GEN_AntiLambda_AntiProton_vz.clear();
	_GEN_AntiLambda_AntiProton_lxy.clear();
	_GEN_AntiLambda_AntiProton_dxy.clear();
	_GEN_AntiLambda_AntiProton_dz.clear();
	_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum.clear();

	_GEN_AntiLambda_Pion_px.clear();
	_GEN_AntiLambda_Pion_py.clear();
	_GEN_AntiLambda_Pion_pz.clear();
	_GEN_AntiLambda_Pion_pt.clear();
	_GEN_AntiLambda_Pion_eta.clear();
	_GEN_AntiLambda_Pion_phi.clear();
	_GEN_AntiLambda_Pion_vx.clear();
	_GEN_AntiLambda_Pion_vy.clear();
	_GEN_AntiLambda_Pion_vz.clear();
	_GEN_AntiLambda_Pion_lxy.clear();
	_GEN_AntiLambda_Pion_dxy.clear();
	_GEN_AntiLambda_Pion_dz.clear();
	_GEN_AntiLambda_Pion_openings_angle_displacement_momentum.clear();

        _GEN_Ks_daughter0_numberOfTrackerLayers.clear();
        _GEN_Ks_daughter1_numberOfTrackerLayers.clear();
        _GEN_AntiLambda_AntiProton_numberOfTrackerLayers.clear();
        _GEN_AntiLambda_Pion_numberOfTrackerLayers.clear();

        _GEN_Ks_daughter0_numberOfTrackerHits.clear();
        _GEN_Ks_daughter1_numberOfTrackerHits.clear();
        _GEN_AntiLambda_AntiProton_numberOfTrackerHits.clear();
        _GEN_AntiLambda_Pion_numberOfTrackerHits.clear();

}

DEFINE_FWK_MODULE(FlatTreeProducerGENSIM);
