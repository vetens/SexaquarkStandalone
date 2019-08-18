#ifndef LambdaKshortFilter_h
#define LambdaKshortFilter_h
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <vector>

  
class LambdaKshortFilter : public edm::EDFilter {

  public:

    explicit LambdaKshortFilter(edm::ParameterSet const& cfg);
    virtual ~LambdaKshortFilter() {}
    virtual bool filter(edm::Event & iEvent, edm::EventSetup const & iSetup);

  private:
  
    edm::InputTag lambdaCollectionTag_;
    edm::InputTag kshortCollectionTag_;
    edm::InputTag genCollectionTag_;
    edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > lambdaCollectionToken_;
    edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > kshortCollectionToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle> >              genCollectionToken_;
    bool isData_;
    unsigned int minNrLambda_,   minNrKshort_;
    double       minPtLambda_,   minPtKshort_;
    double       maxEtaLambda_,  maxEtaKshort_;
    double       minMassLambda_, minMassKshort_;
    double       maxMassLambda_, maxMassKshort_;
    bool checkLambdaDaughters_;
    unsigned int prescaleFalse_, nreject_;

};


#endif
