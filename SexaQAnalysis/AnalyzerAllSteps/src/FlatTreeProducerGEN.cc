#include "../interface/FlatTreeProducerGEN.h"
#include <typeinfo>

FlatTreeProducerGEN::FlatTreeProducerGEN(edm::ParameterSet const& pset):
  m_lookAtAntiS(pset.getUntrackedParameter<bool>("lookAtAntiS")),
  m_runningOnData(pset.getUntrackedParameter<bool>("runningOnData")),
  m_bsTag(pset.getParameter<edm::InputTag>("beamspot")),
  m_genParticlesTag_GEN(pset.getParameter<edm::InputTag>("genCollection_GEN")),

  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN))
{

}


void FlatTreeProducerGEN::beginJob() {

    
        // Initialize when class is created
        edm::Service<TFileService> fs ;

	//very basic info on the charged pions in events
	_tree_pi = fs->make <TTree>("FlatTreeGENLevelPi","treePi");
	_tree_pi->Branch("_pi_eta",&_pi_eta);

	//some GEN Sbar kinematics
        _tree = fs->make <TTree>("FlatTreeGENLevel","tree");
	_tree->Branch("_S_charge",&_S_charge);
	_tree->Branch("_S_mass",&_S_mass);
	_tree->Branch("_S_eta",&_S_eta);
	_tree->Branch("_S_pt",&_S_pt);
	_tree->Branch("_S_pz",&_S_pz);
	_tree->Branch("_S_vx",&_S_vx);
	_tree->Branch("_S_vy",&_S_vy);
	_tree->Branch("_S_vz",&_S_vz);



}

void FlatTreeProducerGEN::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

 

  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  //iEvent.getByToken(m_bsToken, h_bs);

  //SIM particles: normal Gen particles or PlusGEANT
  edm::Handle<vector<reco::GenParticle>> h_genParticles;
  iEvent.getByToken(m_genParticlesToken_GEN, h_genParticles);


  //beamspot
  TVector3 beamspot(0,0,0);
  TVector3 beamspotVariance(0,0,0);
  if(h_bs.isValid()){  
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));			
  }

  nEvents++;

  //first find the GEN particles which are proper antiS, and also save some info on charged pions
  if(!m_runningOnData && m_lookAtAntiS){
	  if(h_genParticles.isValid()){
	      int nPionsThisEvent = 0;
	      int nPionsThisEventEtaSmaller4 = 0;
	      InitPi();
	      for(unsigned int i = 0; i < h_genParticles->size(); ++i){//loop all genparticlesPlusGEANT

			const reco::Candidate * genParticle = &h_genParticles->at(i);
			//count the number of charged pions
			if( abs(genParticle->pdgId() ) == 211 ){
				 nPionsThisEvent++;
				 if(abs(genParticle->eta()) < 4) nPionsThisEventEtaSmaller4++;
				 FillBranchesPion(genParticle->eta());
			}

			//and now for antiS
			if(genParticle->pdgId() != AnalyzerAllSteps::pdgIdAntiS) continue;
			nTotalGENS++;	
			if(genParticle->eta()>0)nTotalGENSPosEta++;	
			if(genParticle->eta()<0)nTotalGENSNegEta++;
				
			FillBranchesGENAntiS(genParticle,beamspot, beamspotVariance);


	      }//for(unsigned int i = 0; i < h_genParticles->size(); ++i)
	      _tree_pi->Fill();
	      nPions = nPions + nPionsThisEvent;
	      nPionsEtaSmaller4 = nPionsEtaSmaller4 + nPionsThisEventEtaSmaller4;
	  }//if(h_genParticles.isValid())
  }


 } //end of analyzer


void FlatTreeProducerGEN::FillBranchesPion(double eta){
	_pi_eta.push_back(eta);	
}

void FlatTreeProducerGEN::FillBranchesGENAntiS(const reco::Candidate  * genParticle, TVector3 beamspot, TVector3 beamspotVariance){
	
	Init(); 

	_S_charge.push_back(genParticle->charge());
	_S_mass.push_back(genParticle->mass());
	_S_eta.push_back(genParticle->eta());
	_S_pt.push_back(genParticle->pt());
	_S_pz.push_back(genParticle->pz());
	_S_vx.push_back(genParticle->vx());	
	_S_vy.push_back(genParticle->vy());	
	_S_vz.push_back(genParticle->vz());	

  	_tree->Fill();

}


void FlatTreeProducerGEN::endJob()
{
}

void
FlatTreeProducerGEN::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void
FlatTreeProducerGEN::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
FlatTreeProducerGEN::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
FlatTreeProducerGEN::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatTreeProducerGEN::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

FlatTreeProducerGEN::~FlatTreeProducerGEN()
{
	string particle = "S";
	if(m_lookAtAntiS){
		particle = "anti-S";
	}

	std::cout << "The total number of events: " << nEvents << std::endl;
	std::cout << "The total number of pions: " << nPions << std::endl;
	std::cout << "The total number of pions with |eta| smaller than 4: " << nPionsEtaSmaller4 << std::endl;
	

	std::cout << "The total number GEN " << particle << " that were found is: " << nTotalGENS << std::endl; 
	std::cout << "The total number GEN " << particle << " that were found with pos eta is: " << nTotalGENSPosEta << std::endl; 
	std::cout << "The total number GEN " << particle << " that were found with neg eta is: " << nTotalGENSNegEta << std::endl; 
}


void FlatTreeProducerGEN::InitPi()
{
	_pi_eta.clear();
}

void FlatTreeProducerGEN::Init()
{

    	_S_charge.clear();
    	_S_mass.clear();
	_S_eta.clear();
	_S_pt.clear();
	_S_pz.clear();
	_S_vx.clear();
	_S_vy.clear();
	_S_vz.clear();

}

DEFINE_FWK_MODULE(FlatTreeProducerGEN);
