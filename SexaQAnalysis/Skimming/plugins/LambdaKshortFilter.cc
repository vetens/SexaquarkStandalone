
#include "SexaQAnalysis/Skimming/plugins/LambdaKshortFilter.h"


LambdaKshortFilter::LambdaKshortFilter(edm::ParameterSet const& pset):
  lambdaCollectionTag_(pset.getParameter<edm::InputTag>("lambdaCollection")),
  kshortCollectionTag_(pset.getParameter<edm::InputTag>("kshortCollection")),
  genCollectionTag_   (pset.getParameter<edm::InputTag>("genCollection")),
  isData_       (pset.getParameter<bool>  ("isData")),
  minNrLambda_  (pset.getParameter<unsigned int>("minNrLambda")),
  minNrKshort_  (pset.getParameter<unsigned int>("minNrKshort")),
  minPtLambda_  (pset.getParameter<double>("minPtLambda")),
  minPtKshort_  (pset.getParameter<double>("minPtKshort")),
  maxEtaLambda_ (pset.getParameter<double>("maxEtaLambda")),
  maxEtaKshort_ (pset.getParameter<double>("maxEtaKshort")),
  minMassLambda_(pset.getParameter<double>("minMassLambda")),
  minMassKshort_(pset.getParameter<double>("minMassKshort")),
  maxMassLambda_(pset.getParameter<double>("maxMassLambda")),
  maxMassKshort_(pset.getParameter<double>("maxMassKshort")),
  checkLambdaDaughters_(pset.getParameter<bool>("checkLambdaDaughters")),
  prescaleFalse_(pset.getParameter<unsigned int>("prescaleFalse"))
{
  lambdaCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(lambdaCollectionTag_);
  kshortCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(kshortCollectionTag_);
  genCollectionToken_    = consumes<std::vector<reco::GenParticle> >             (genCollectionTag_);
  nreject_ = 0;
  produces<reco::CandidatePtrVector>("kshort");
  produces<reco::CandidatePtrVector>("lambda");
}


bool LambdaKshortFilter::filter(edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  auto kshorts = std::make_unique<reco::CandidatePtrVector>();
  auto lambdas = std::make_unique<reco::CandidatePtrVector>();

  // select on reco lambdas and kaons
  if (isData_) {
    // read out lambdas
    edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_lambda;
    iEvent.getByToken(lambdaCollectionToken_, h_lambda);
    if(!h_lambda.isValid()) {
      std::cout << "Missing collection during LambdaKshortFilter : " << lambdaCollectionTag_ << " ... skip entry !" << std::endl;
      return false;
    }
    // read out kaons
    edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_kshort;
    iEvent.getByToken(kshortCollectionToken_ , h_kshort);
    if(!h_kshort.isValid()) {
      std::cout << "Missing collection during LambdaKshortFilter : " << kshortCollectionTag_ << " ... skip entry !" << std::endl;
      return false;
    }
    std::cout << "LambdaKshortFilter: Starting a new event" << std::endl;
    // select the lambdas passing kinematic cuts
    for (unsigned int l = 0; l < h_lambda->size(); ++l) {
      if (h_lambda->at(l).pt()       > minPtLambda_   &&
          fabs(h_lambda->at(l).eta()) < maxEtaLambda_ &&
          h_lambda->at(l).mass()     > minMassLambda_ &&
          h_lambda->at(l).mass()     < maxMassLambda_) {
        edm::Ptr<reco::VertexCompositeCandidate> lptr(h_lambda,l);
        lambdas->push_back(std::move(lptr));
      }
    }

    // select the kshorts passing kinematic cuts and non-overlapping with lambdas
    for (unsigned int k = 0; k < h_kshort->size(); ++k) {
      if (h_kshort->at(k).pt()       > minPtKshort_   &&
          fabs(h_kshort->at(k).eta()) < maxEtaKshort_ &&
	  h_kshort->at(k).mass()     > minMassKshort_ &&
	  h_kshort->at(k).mass()     < maxMassKshort_) {
        edm::Ptr<reco::VertexCompositeCandidate> kptr(h_kshort,k);
	// check for overlaps with the lambdas, and keep the lambda in case
        bool overlap = false;
        for (auto lptr : *lambdas) {
          for (unsigned int li = 0; li < lptr->numberOfDaughters() && !overlap; ++li) {
            for (unsigned int ki = 0; ki < kptr->numberOfDaughters() && !overlap; ++ki) {
	      if (lptr->daughter(li)->px() == kptr->daughter(ki)->px() &&
	          lptr->daughter(li)->py() == kptr->daughter(ki)->py() &&
	          lptr->daughter(li)->pz() == kptr->daughter(ki)->pz()) {
                overlap = true;
	      }
	    }
          }
          if (overlap){
		 //std::cout << "LambdaKshortFilter: OVERLAP FOUND" << std::endl;
		 break;
	  }	
        }
        if (!overlap) kshorts->push_back(std::move(kptr));
      //  kshorts->push_back(std::move(kptr));
      }
    }


  // if not data, then select on gen particles
  } else {

    // read out genparticles
    edm::Handle<std::vector<reco::GenParticle> > h_genparts;
    iEvent.getByToken(genCollectionToken_, h_genparts);
    if(!h_genparts.isValid()) {
      std::cout << "Missing collection : " << genCollectionTag_ << " ... skip entry !" << std::endl;
      return false;
    }

    // find lambdas and kshorts
    for (unsigned int p = 0; p < h_genparts->size(); ++p) {
      if (fabs(h_genparts->at(p).pdgId()) == 3122 &&
          fabs(h_genparts->at(p).eta()) < maxEtaLambda_ &&
	  h_genparts->at(p).pt() > minPtLambda_) {
        // check the decay to be to a (anti)proton (and a pion) -> need Geant collected genparticles to do that
        if (checkLambdaDaughters_) { // make sure to collect geant genparticles first
          if (!h_genparts->at(p).daughter(0) || !h_genparts->at(p).daughter(1)) {
            continue;
          }
          if (fabs(h_genparts->at(p).daughter(0)->pdgId()) != 2212 &&
              fabs(h_genparts->at(p).daughter(1)->pdgId()) != 2212) continue;
        }
        edm::Ptr<reco::GenParticle> lptr(h_genparts,p);
        lambdas->push_back(std::move(lptr));
      }
      if ((fabs(h_genparts->at(p).pdgId()) == 310) &&
          fabs(h_genparts->at(p).eta()) < maxEtaKshort_ &&
	  h_genparts->at(p).pt() > minPtKshort_) {
	edm::Ptr<reco::GenParticle> kptr(h_genparts,p);
        kshorts->push_back(std::move(kptr));
      }
    }


  }

  // get the vector sizes before they disappear when putting in the event
  unsigned int nl = lambdas->size(), nk = kshorts->size();

  iEvent.put(std::move(lambdas),"lambda");
  iEvent.put(std::move(kshorts),"kshort");
/*  if(kshorts != NULL){
  for(int k = 0; k<(int)kshorts->size(); k++){
        std::cout << "LambdaKshortFilter: kshort momenta " << k << "px,py,pz: " << (*kshorts)[k]->px() << "," << (*kshorts)[k]->py() << "," << (*kshorts)[k]->pz() << std::endl; 
  }
  }

  if(lambdas != NULL){
  for(int l = 0; l<(int)lambdas->size(); l++){
        std::cout << "LambdaKshortFilter: lambda momenta " << l << "px,py,pz: " << (*lambdas)[l]->px() << "," << (*lambdas)[l]->py() << "," << (*lambdas)[l]->pz() << std::endl; 
  }
  }
*/
  //std::cout << "LambdaKshortFilter: n lambdas put in events " << nl << std::endl;
  //std::cout << "LambdaKshortFilter: n kshorts put in events " << nk << std::endl;

  // throw away events on data without sufficient lambdas or kshorts
  if (nl < minNrLambda_ || nk < minNrKshort_) {
    ++nreject_;
    return (prescaleFalse_ ? !(nreject_ % prescaleFalse_) : false);
  }
  // if we reach here there's a sufficient number of good lambdas and kshorts
 
  //for debugging: check that the lambda and kshort which you save are really different

/*  for(unsigned int l = 0; l < lambdas->size(); ++l){
 	std::cout << "put a lambda" << std::endl;
  }
  for(unsigned int k = 0; k < kshorts->size(); ++k){
	std::cout << "put a kshort" << std::endl;
  }

  for(unsigned int l = 0; l < lambdas->size(); ++l){
  	for(unsigned int k = 0; k < kshorts->size(); ++k){
	 if((*lambdas)[l]->px() == (*kshorts)[k]->px()){
  		std::cout << "------------------------------------------" << std::endl;
      		std::cout << "Lambda number " << l << " momenta: "<< (*lambdas)[l]->px()  << ", " << (*lambdas)[l]->py() << ", " <<  (*lambdas)[l]->pz() << std::endl;
      		std::cout << "Kshort number " << k << " momenta: "<<  (*kshorts)[k]->px() << ", " << (*kshorts)[k]->py() << ", " <<  (*kshorts)[k]->pz() << std::endl;
         }
  	}
  }
*/


   return true;

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(LambdaKshortFilter);
