#ifndef MassFilter_h
#define MassFilter_h
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

  
class MassFilter : public edm::EDFilter {

  public:

    explicit MassFilter(edm::ParameterSet const& cfg);
    virtual ~MassFilter() {}
    virtual bool filter(edm::Event & iEvent, edm::EventSetup const & iSetup);

  private:

    edm::InputTag lkPairCollectionTag_;
    edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > lkPairCollectionToken_;
    double minMass_, maxMass_, targetMass_;
    unsigned int prescaleFalse_, nreject_;
    reco::LeafCandidate::LorentzVector n_;

};


#endif

