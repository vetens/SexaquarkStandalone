
#include "SexaQAnalysis/Skimming/plugins/MassFilter.h"


MassFilter::MassFilter(edm::ParameterSet const & pset) :
  lkPairCollectionTag_(pset.getParameter<edm::InputTag>("lambdakshortCollection")),
  minMass_            (pset.getParameter<double>("minMass")),
  maxMass_            (pset.getParameter<double>("maxMass")),
  // mass of a potential fixed target particle, to estimate the "missing" mass,
  // i.e. the mass of the incoming particle. Set to zero to get invariant mass.
  targetMass_         (pset.getParameter<double>("targetMass")),
  prescaleFalse_      (pset.getParameter<unsigned int>("prescaleFalse"))
{
  lkPairCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(lkPairCollectionTag_);
  n_ = reco::LeafCandidate::LorentzVector(0,0,0,targetMass_);
  nreject_ = 0;
  produces<std::vector<reco::VertexCompositePtrCandidate> >("sVertexCompositePtrCandidate");
}


bool MassFilter::filter(edm::Event & iEvent, edm::EventSetup const & iSetup)
{

  edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_lkPair;
  iEvent.getByToken(lkPairCollectionToken_, h_lkPair);
  if(!h_lkPair.isValid()) {
    std::cout << "Missing collection during MassFilter: " << lkPairCollectionTag_ << " ... skip entry !" << std::endl;
    return false;
  }

  auto lkPairs = std::make_unique<std::vector<reco::VertexCompositePtrCandidate> >();

  // find any Lambda - Kshort pair that matches the mass window
  for (auto lk : *h_lkPair) {
    // calculate the missing mass
    reco::LeafCandidate::LorentzVector m = lk.daughter(0)->p4() + lk.daughter(1)->p4() - n_;
    // impose the mass window
    if (m.M() > minMass_ && m.M() < maxMass_) {
      reco::VertexCompositePtrCandidate  lkPass(lk.charge(), lk.p4(), lk.vertex(), lk.vertexCovariance(), lk.vertexChi2(), lk.vertexNdof(), 0, 0, true);
      lkPass.setP4(m);
      lkPass.addDaughter(lk.sourceCandidatePtr(0));
      lkPass.addDaughter(lk.sourceCandidatePtr(1));
      lkPairs->push_back(std::move(lkPass));
    }
  }

  // get the vector size before they disappear when putting in the event
  unsigned int n = lkPairs->size();

  iEvent.put(std::move(lkPairs),"sVertexCompositePtrCandidate");

  // throw away events on data without good lambda-kshort pairs
  if (n == 0) {
    ++nreject_;
    return (prescaleFalse_ ? !(nreject_ % prescaleFalse_) : false);
  }
  // if we reach here there's a good lambda-kshort pair
  return true;

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MassFilter);
