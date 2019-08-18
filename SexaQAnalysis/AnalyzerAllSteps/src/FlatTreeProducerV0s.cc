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
	_tree_Ks->Branch("_Ks_dxy",&_Ks_dxy);
	_tree_Ks->Branch("_Ks_dz",&_Ks_dz);
	_tree_Ks->Branch("_Ks_dz_min",&_Ks_dz_min);


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
	_tree_Lambda->Branch("_Lambda_dz",&_Lambda_dz);
	_tree_Lambda->Branch("_Lambda_dz_min",&_Lambda_dz_min);


}

void FlatTreeProducerV0s::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


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
  if(h_muons.isValid()){
	for(unsigned int i = 0; i < h_muons->size(); ++i){

		//now check if this muon i is tight and isolated, has a minmal pt and is in the acceptance, if not do not consider this muon
		//bool muon1Tight = muon::isTightMuon(h_muons->at(i), h_offlinePV->at(0));
		//if(!muon1Tight) continue;
		if(h_muons->at(i).pt() < 20) continue;	
		if(fabs(h_muons->at(i).eta()) > 2.4) continue;	

		for(unsigned int j = i + 1; j < h_muons->size(); ++j){
			
			//now check if this muon j is tight and isolated, has a minmal pt and is in the acceptance, if not do not consider this muon
			bool muon2Tight = muon::isTightMuon(h_muons->at(j), h_offlinePV->at(0));
			if(!muon2Tight) continue;

			if(h_muons->at(j).pt() < 20) continue;	
			if(fabs(h_muons->at(j).eta()) > 2.4) continue;	

			LorentzVector p4ZCandidate = h_muons->at(i).p4() + h_muons->at(j).p4();

			//now check if there is no jet (with pt>30GeV) going back to back with the Z 
			bool backToBackHighPtJetFound = false;
			for(unsigned int k = 0; k < h_jets->size(); ++k){
				if(h_jets->at(k).pt() < 30) continue;
				double deltaPhiJetZ = reco::deltaPhi(p4ZCandidate.phi(),h_jets->at(k).phi());
				if(abs(deltaPhiJetZ - TMath::Pi()) < 0.1) backToBackHighPtJetFound = true;

			}

			double invDiMuonMass = p4ZCandidate.mass();

			//ask for the inv mass to be within 15GeV from the known Z mass
			if(abs(invDiMuonMass - 91.1876) < 15/2 && !backToBackHighPtJetFound){ZCandidateMass = invDiMuonMass;  ZCandidatePresent = true; ZCandidatePhi = p4ZCandidate.phi();}
		}


	}
  }

  if(ZCandidatePresent) std::cout << "ZCandMass = " << ZCandidateMass << " ZCandidatePhi " << ZCandidatePhi  << std::endl;

  //now that you have the phi of the Z candidate you know in which direction the hard event goes, so now look at the V0s in the underlying event, which are the V0s which are not in the deltaPhi cone of the hard event, this cone is defined as 60° around the Z and 60° around the Z in the back to back

  if(ZCandidatePresent){

	  if(h_V0Ks.isValid()){
	      for(unsigned int i = 0; i < h_V0Ks->size(); ++i){//loop all RECO Ks
		const reco::VertexCompositeCandidate * Ks = &h_V0Ks->at(i);
		double deltaPhiKsHardCone = reco::deltaPhi( Ks->phi() , ZCandidatePhi ); 
		double deltaPhiKsBackToBackHardCone = reco::deltaPhi( Ks->phi() , -TMath::Pi() + ZCandidatePhi ); 
		if(abs(deltaPhiKsHardCone) <  TMath::Pi()/3 ||  abs(deltaPhiKsBackToBackHardCone) <  TMath::Pi()/3)FillBranchesKs(Ks, beamspot, beamspotVariance, h_offlinePV);
	      }
	  }

	  if(h_V0L.isValid()){
	      for(unsigned int i = 0; i < h_V0L->size(); ++i){//loop all RECO Lambdas
		const reco::VertexCompositeCandidate * L = &h_V0L->at(i);
		double deltaPhiLHardCone = reco::deltaPhi( L->phi() , ZCandidatePhi );
		double deltaPhiLBackToBackHardCone = reco::deltaPhi( L->phi() , -TMath::Pi() + ZCandidatePhi );
		if(abs(deltaPhiLHardCone) <  TMath::Pi()/3 ||  abs(deltaPhiLBackToBackHardCone) <  TMath::Pi()/3)FillBranchesLambda(L, beamspot, beamspotVariance, h_offlinePV);
	      }
	  }

 }


 } //end of analyzer




void FlatTreeProducerV0s::FillBranchesKs(const reco::VertexCompositeCandidate * RECOKs, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV){

        TVector3 KsCreationVertex(RECOKs->vx(),RECOKs->vy(),RECOKs->vz());
        double Lxy = AnalyzerAllSteps::lxy(beamspot,KsCreationVertex);
        TVector3 KsMomentum(RECOKs->px(),RECOKs->py(),RECOKs->pz());
        double dxy = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,beamspot);
        double dz = AnalyzerAllSteps::dz_line_point(KsCreationVertex,KsMomentum,beamspot);
        TVector3 PVmin = AnalyzerAllSteps::dz_line_point_min(KsCreationVertex,KsMomentum,h_offlinePV);
        double dz_min = AnalyzerAllSteps::dz_line_point(KsCreationVertex,KsMomentum,PVmin);





	InitKs();	

	_Ks_mass.push_back(RECOKs->mass());	

	_Ks_pt.push_back(RECOKs->pt());	
	_Ks_pz.push_back(RECOKs->pz());	
	_Ks_Lxy.push_back(Lxy);	
	_Ks_vz.push_back(RECOKs->vz());	

	_Ks_eta.push_back(RECOKs->eta());
	_Ks_phi.push_back(RECOKs->phi());

	_Ks_dxy.push_back(dxy);	
	_Ks_dz.push_back(dz);	
	_Ks_dz_min.push_back(dz_min);	

  	_tree_Ks->Fill();

}

void FlatTreeProducerV0s::FillBranchesLambda(const reco::VertexCompositeCandidate * RECOLambda, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV){


	TVector3 LambdaCreationVertex(RECOLambda->vx(),RECOLambda->vy(),RECOLambda->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,LambdaCreationVertex);
	TVector3 LambdaMomentum(RECOLambda->px(),RECOLambda->py(),RECOLambda->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(LambdaCreationVertex,LambdaMomentum,beamspot);
	double dz = AnalyzerAllSteps::dz_line_point(LambdaCreationVertex,LambdaMomentum,beamspot);
	TVector3 PVmin = AnalyzerAllSteps::dz_line_point_min(LambdaCreationVertex,LambdaMomentum,h_offlinePV);
	double dz_min = AnalyzerAllSteps::dz_line_point(LambdaCreationVertex,LambdaMomentum,PVmin);


	InitLambda();	

	_Lambda_mass.push_back(RECOLambda->mass());	

	_Lambda_pt.push_back(RECOLambda->pt());	
	_Lambda_pz.push_back(RECOLambda->pz());	
	_Lambda_Lxy.push_back(Lxy);	
	_Lambda_vz.push_back(RECOLambda->vz());	

	_Lambda_eta.push_back(RECOLambda->eta());
	_Lambda_phi.push_back(RECOLambda->phi());

	_Lambda_dxy.push_back(dxy);	
	_Lambda_dz.push_back(dz);	
	_Lambda_dz_min.push_back(dz_min);	

  	_tree_Lambda->Fill();

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

	_Ks_dxy.clear();	
	_Ks_dz.clear();
	_Ks_dz_min.clear();	
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
	_Lambda_dz.clear();
	_Lambda_dz_min.clear();	
}


DEFINE_FWK_MODULE(FlatTreeProducerV0s);
