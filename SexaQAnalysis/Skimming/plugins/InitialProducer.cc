// -*- C++ -*-
//
// Package:    SexaQAnalysis/Skimming
// Class:      InitialProducer
// 
/**\class InitialProducer InitialProducer.cc SexaQAnalysis/Skimming/plugins/InitialProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  localusers user
//         Created:  Tue, 24 Jul 2018 10:26:06 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
//
// class declaration
//

class InitialProducer : public edm::stream::EDProducer<> {
   public:
      explicit InitialProducer(const edm::ParameterSet&);
      ~InitialProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      //typedefs
      typedef math::XYZTLorentzVector LorentzVector;

   private:
      //************remove 2 lines below
      edm::InputTag trackCollectionTag_;
      edm::InputTag lambdaCollectionTag_;
      edm::InputTag kshortCollectionTag_;
      edm::InputTag offlinePrimaryVerticesCollectionTag_;
      edm::InputTag ak4PFJetsCollectionTag_;
      edm::InputTag muonsCollectionTag_;
      edm::InputTag electronsCollectionTag_;
      edm::InputTag METCollectionTag_;
      edm::EDGetTokenT<std::vector<reco::Track> > tracksCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > lambdaCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > kshortCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > offlinePrimaryVerticesCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::PFJet> > ak4PFJetsCollectionToken_;
      edm::EDGetTokenT<std::vector<edm::FwdPtr<reco::PFCandidate> > > muonsCollectionToken_;
      edm::EDGetTokenT<std::vector<edm::FwdPtr<reco::PFCandidate> > > electronsCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::PFMET>  > METCollectionToken_;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::InputTag src_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
//InitialProducer::InitialProducer(const edm::ParameterSet& iConfig)
InitialProducer::InitialProducer(edm::ParameterSet const& pset):
trackCollectionTag_(pset.getParameter<edm::InputTag>("TrackCollection")),
lambdaCollectionTag_(pset.getParameter<edm::InputTag>("lambdaCollection")),
kshortCollectionTag_(pset.getParameter<edm::InputTag>("kshortCollection")),
offlinePrimaryVerticesCollectionTag_(pset.getParameter<edm::InputTag>("offlinePrimaryVerticesCollection")),
ak4PFJetsCollectionTag_(pset.getParameter<edm::InputTag>("ak4PFJetsCollection")),
muonsCollectionTag_(pset.getParameter<edm::InputTag>("muonsCollection")),
electronsCollectionTag_(pset.getParameter<edm::InputTag>("electronsCollection")),
METCollectionTag_(pset.getParameter<edm::InputTag>("METCollection"))
{
   //register your products
   tracksCollectionToken_ = consumes<std::vector<reco::Track> >(trackCollectionTag_);
   lambdaCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(lambdaCollectionTag_);
   kshortCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(kshortCollectionTag_);
   offlinePrimaryVerticesCollectionToken_ = consumes<std::vector<reco::Vertex> >(offlinePrimaryVerticesCollectionTag_);
   ak4PFJetsCollectionToken_ = consumes<std::vector<reco::PFJet> >(ak4PFJetsCollectionTag_);
   muonsCollectionToken_ = consumes<std::vector<edm::FwdPtr<reco::PFCandidate> > >(muonsCollectionTag_);
   electronsCollectionToken_ = consumes<std::vector<edm::FwdPtr<reco::PFCandidate> > >(electronsCollectionTag_);
   METCollectionToken_ = consumes<std::vector<reco::PFMET> >(METCollectionTag_);
   produces<std::vector<int>>("ntracks");
   produces<std::vector<int>>("nlambdas");
   produces<std::vector<int>>("nkshorts");
   produces<std::vector<int>>("nkshortsAndNlambdas");
   produces<std::vector<int>>("nPVs");
   produces<std::vector<int>>("njets");
   produces<std::vector<reco::Particle::LorentzVector>>("TwoTopJets");
   produces<std::vector<int>>("nmuons");
   produces<std::vector<int>>("nelectrons");
   produces<std::vector<double>>("HT");
   produces<std::vector<double>>("TKHT");
   produces<std::vector<LorentzVector>>("MET");
   produces<std::vector<math::XYZVector>>("TKMET");
  
}


InitialProducer::~InitialProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
InitialProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco; 
   using namespace std;
  //std::cout << "InitialProducer: starting for a new event" << std::endl; 
  //ntracks part
  edm::Handle<std::vector<reco::Track >> h_tracks;
  iEvent.getByToken(tracksCollectionToken_, h_tracks);
  if(!h_tracks.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << trackCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto ntracks = std::make_unique<std::vector<int>>();
  ntracks->push_back((int)h_tracks->size());
  iEvent.put(std::move(ntracks), "ntracks");

  //nlambdas part
  edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_lambdas;
  iEvent.getByToken(lambdaCollectionToken_, h_lambdas);
  if(!h_lambdas.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << lambdaCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nlambdas = std::make_unique<std::vector<int>>();
  nlambdas->push_back((int)h_lambdas->size());
  iEvent.put(std::move(nlambdas), "nlambdas");

  //nkshorts part
  edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_kshorts;
  iEvent.getByToken(kshortCollectionToken_, h_kshorts);
  if(!h_kshorts.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << kshortCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nkshorts = std::make_unique<std::vector<int>>();
  nkshorts->push_back((int)h_kshorts->size());
  iEvent.put(std::move(nkshorts), "nkshorts");

  //nkshorts and nlambdas part: if at least 1 kshort and at least 1 lambda is present then there is the potential to reconstruct an S
  auto nkshortsAndNlambdas = std::make_unique<std::vector<int>>();
  int nkshort_and_nlambdas_present = 0;
/*  if((int)h_kshorts->size() >= 1) std::cout << "InitialProducer: Found " << (int)h_kshorts->size() << "kshorts" << std::endl;
  if((int)h_lambdas->size() >= 1) std::cout << "InitialProduer: Found " << (int)h_lambdas->size() << "lambdas" << std::endl;
  if((int)h_kshorts->size() >= 1 && (int)h_lambdas->size() >= 1) std::cout << "InitialProducer: Found both a kshort and a lambda: " << std::endl;
  for(int k = 0; k<(int)h_kshorts->size(); k++){
	std::cout << "Initialproducer: kshort " << k << "px,py,pz: " << h_kshorts->at(k).px() << "," << h_kshorts->at(k).py() << "," << h_kshorts->at(k).pz() << std::endl; 
  }
  for(int l = 0; l<(int)h_lambdas->size(); l++){
	std::cout << "Initialproducer: lambda " << l << "px,py,pz: " << h_lambdas->at(l).px() << "," << h_lambdas->at(l).py() << "," << h_lambdas->at(l).pz() << std::endl; 
  }*/
  nkshortsAndNlambdas->push_back((int) nkshort_and_nlambdas_present);
  iEvent.put(std::move(nkshortsAndNlambdas), "nkshortsAndNlambdas");

  //nPVs part
  edm::Handle<std::vector<reco::Vertex> > h_PVs;
  iEvent.getByToken(offlinePrimaryVerticesCollectionToken_, h_PVs);
  if(!h_PVs.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << offlinePrimaryVerticesCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nPVs = std::make_unique<std::vector<int>>();
  nPVs->push_back((int)h_PVs->size());
  iEvent.put(std::move(nPVs), "nPVs");

  //njets part and HT part
  edm::Handle<std::vector<reco::PFJet> > h_jets;
  iEvent.getByToken(ak4PFJetsCollectionToken_, h_jets);
  if(!h_jets.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << ak4PFJetsCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto njets = std::make_unique<std::vector<int>>();
  njets->push_back((int)h_jets->size());
  iEvent.put(std::move(njets), "njets");

  //select the 2 jets with highest momentum
  auto TwoTopJets = std::make_unique<std::vector<reco::PFJet>>();
  reco::PFJet dummyPFJet; 
  TwoTopJets->push_back(dummyPFJet); 
  TwoTopJets->push_back(dummyPFJet);
  double sumJetpT = 0;
  for (unsigned int j = 0; j < h_jets->size(); ++j) {

	double thisJetMomentum = h_jets->at(j).p();
	if(thisJetMomentum > TwoTopJets->at(0).p()) TwoTopJets->at(0)  = h_jets->at(j);
	else if(thisJetMomentum > TwoTopJets->at(1).p()) TwoTopJets->at(1)  = h_jets->at(j);
        
	sumJetpT = sumJetpT + h_jets->at(j).pt();
 
  }
 
  auto highestMomentJetsLorentzVectors = std::make_unique<std::vector<reco::Particle::LorentzVector>>();
  highestMomentJetsLorentzVectors->push_back(TwoTopJets->at(0).p4());
  highestMomentJetsLorentzVectors->push_back(TwoTopJets->at(1).p4());
  iEvent.put(std::move(highestMomentJetsLorentzVectors), "TwoTopJets");

  auto HT = std::make_unique<std::vector<double>>();
  HT->push_back(sumJetpT);
  iEvent.put(std::move(HT), "HT");

  //nmuons part
  edm::Handle<std::vector<edm::FwdPtr<reco::PFCandidate>> > h_muons;
  iEvent.getByToken(muonsCollectionToken_, h_muons);
  if(!h_muons.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << muonsCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nmuons = std::make_unique<std::vector<int>>();
  nmuons->push_back((int)h_muons->size());
  iEvent.put(std::move(nmuons), "nmuons");

  //nelectrons part
  edm::Handle<std::vector<edm::FwdPtr<reco::PFCandidate>> > h_electrons;
  iEvent.getByToken(electronsCollectionToken_, h_electrons);
  if(!h_electrons.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << electronsCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nelectrons = std::make_unique<std::vector<int>>();
  nelectrons->push_back((int)h_electrons->size());
  iEvent.put(std::move(nelectrons), "nelectrons");

  //MET part
  edm::Handle<vector<reco::PFMET> > h_MET;
  iEvent.getByToken(METCollectionToken_, h_MET);
  if(!h_MET.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << METCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto MET = std::make_unique<std::vector<LorentzVector>>();
  MET->push_back(h_MET->at(0).p4()); 
  iEvent.put(std::move(MET), "MET");

  //TKMET and TKHT part 
  //TKMET: calculate the missing pT from the tracks: so sum all the tracks and take the minus sign
  //TKHT: scalar sim of all track pts
  math::XYZVector sumTK(0,0,0);
  double sumTrackpT = 0;
  for(unsigned int j = 0; j < h_tracks->size(); ++j) {
	//sumTK = sumTK - h_tracks->at(j).innerMomentum();
	sumTrackpT = sumTrackpT + h_tracks->at(j).pt();
  }	
  auto TKMET = std::make_unique<std::vector<math::XYZVector>>();
  TKMET->push_back(sumTK);
  iEvent.put(std::move(TKMET), "TKMET");

  auto TKHT = std::make_unique<std::vector<double>>();
  TKHT->push_back(sumTrackpT);
  iEvent.put(std::move(TKHT), "TKHT");


}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
InitialProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
InitialProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
InitialProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
InitialProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
InitialProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
InitialProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
InitialProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(InitialProducer);
