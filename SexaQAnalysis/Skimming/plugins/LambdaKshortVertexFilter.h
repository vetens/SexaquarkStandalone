#ifndef LambdaKshortVertexFilter_h
#define LambdaKshortVertexFilter_h
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "RecoVertex/KalmanVertexFit/interface/SimpleVertexTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include <vector>
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"  
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Vector.h"

class LambdaKshortVertexFilter : public edm::EDFilter {

  public:
    typedef math::XYZPoint Point;
    typedef math::XYZVector Vector;

    explicit LambdaKshortVertexFilter(edm::ParameterSet const& cfg);
    virtual ~LambdaKshortVertexFilter() {}
    virtual bool filter(edm::Event & iEvent, edm::EventSetup const & iSetup);
    ParticleMass charged_pi_mass = 0.13957061;
    ParticleMass KshortMass = 0.497611;
    ParticleMass proton_mass = 0.9382720813;
    ParticleMass LambdaMass =  1.115683;
    //set the sigma on some of the masses 1000 times higher then the real world average values to give the fit some liberty
    float charged_pi_mass_sigma = 0.00000024;
    float KshortMassSigma = 0.013/1000;
    float proton_mass_sigma = 0.001;
    float LambdaMassSigma = 0.006;
    //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
    float chi = 0.;
    float ndf = 0.;

 private:

    edm::InputTag lambdaCollectionTag_;
    edm::InputTag kshortCollectionTag_;
    edm::InputTag genCollectionTag_;
  //  edm::InputTag beamspotCollectionTag_;
    
    edm::EDGetTokenT<reco::CandidatePtrVector> lambdaCollectionToken_;
    edm::EDGetTokenT<reco::CandidatePtrVector> kshortCollectionToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genCollectionToken_;
  //  edm::EDGetTokenT<reco::BeamSpot> beamspotCollectionToken_;

    double maxchi2ndofVertexFit_;

    //functions 
    bool allCollectionValid(edm::Handle<reco::CandidatePtrVector> h_lambda,edm::Handle<reco::CandidatePtrVector> h_kshort);
    bool checkRefCountedKinematicTree(RefCountedKinematicTree Tree); 
    RefCountedKinematicTree KinfitTwoTTracks(reco::TransientTrack ttrack1, reco::TransientTrack ttrack2, ParticleMass trackMass1, float trackMassSigma1, ParticleMass trackMass2, float trackMassSigma2, ParticleMass combinedMass, float combinedMassSigma);
    RefCountedKinematicParticle getTopParticleFromTree(RefCountedKinematicTree Tree);
    RefCountedKinematicVertex returnVertexFromTree(const RefCountedKinematicTree& myTree) const;
    
};


#endif

