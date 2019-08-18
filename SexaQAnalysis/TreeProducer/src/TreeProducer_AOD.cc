#include "SexaQAnalysis/TreeProducer/interface/TreeProducer_AOD.h"


//
// constructors and destructor
//

TreeProducer_AOD::TreeProducer_AOD(edm::ParameterSet const& pset):
trackCollectionTag_(pset.getParameter<edm::InputTag>("TrackCollection")),
lambdaCollectionTag_(pset.getParameter<edm::InputTag>("lambdaCollection")),
kshortCollectionTag_(pset.getParameter<edm::InputTag>("kshortCollection")),
offlinePrimaryVerticesCollectionTag_(pset.getParameter<edm::InputTag>("offlinePrimaryVerticesCollection")),
ak4PFJetsCollectionTag_(pset.getParameter<edm::InputTag>("ak4PFJetsCollection")),
muonsCollectionTag_(pset.getParameter<edm::InputTag>("muonsCollection")),
electronsCollectionTag_(pset.getParameter<edm::InputTag>("electronsCollection")),
METCollectionTag_(pset.getParameter<edm::InputTag>("METCollection")),
tracksCollectionToken_(consumes<std::vector<reco::Track> >(trackCollectionTag_)),
lambdaCollectionToken_(consumes<std::vector<reco::VertexCompositeCandidate> >(lambdaCollectionTag_)),
kshortCollectionToken_(consumes<std::vector<reco::VertexCompositeCandidate> >(kshortCollectionTag_)),
offlinePrimaryVerticesCollectionToken_(consumes<std::vector<reco::Vertex> >(offlinePrimaryVerticesCollectionTag_)),
ak4PFJetsCollectionToken_(consumes<std::vector<reco::PFJet> >(ak4PFJetsCollectionTag_)),
muonsCollectionToken_(consumes<std::vector<edm::FwdPtr<reco::PFCandidate> > >(muonsCollectionTag_)),
electronsCollectionToken_(consumes<std::vector<edm::FwdPtr<reco::PFCandidate> > >(electronsCollectionTag_)),
METCollectionToken_(consumes<std::vector<reco::PFMET> >(METCollectionTag_))
{
}

TreeProducer_AOD::~TreeProducer_AOD()
{
}

//
// member functions
//

// ------------ method called for each event  ------------
void
TreeProducer_AOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Initialize branches
  Init();

  // HANDLES //
  // Get collections
/*  edm::Handle<vector<reco::GenParticle> > H_partons;
  if(!_isData){
    iEvent.getByToken(m_partons, H_partons);
    if(!H_partons.isValid()) {
      if(verbose>0) cout << "Missing collection during TreeProducer_AOD: genParticles ... skip entry !" << endl;
      return;
    }
  }
*/
  //ntracks part
  edm::Handle<std::vector<reco::Track >> h_tracks;
  iEvent.getByToken(tracksCollectionToken_, h_tracks);
  if(!h_tracks.isValid()) {
      std::cout << "Missing collection during TreeProducer_AOD : " << trackCollectionTag_ << " ... skip entry !" << std::endl;
  }
  _nTrack = h_tracks->size();

  //nlambdas part
  edm::Handle<std::vector<reco::VertexCompositeCandidate >> h_lambdas;
  iEvent.getByToken(lambdaCollectionToken_, h_lambdas);
  if(!h_lambdas.isValid()) {
      std::cout << "Missing collection during TreeProducer_AOD : " << lambdaCollectionTag_ << " ... skip entry !" << std::endl;
  }
  _nLambdas = h_lambdas->size();

  //nkshorts part
  edm::Handle<std::vector<reco::VertexCompositeCandidate >> h_kshorts;
  iEvent.getByToken(kshortCollectionToken_, h_kshorts);
  if(!h_kshorts.isValid()) {
      std::cout << "Missing collection during TreeProducer_AOD : " << kshortCollectionTag_ << " ... skip entry !" << std::endl;
  }
  _nKshorts = h_kshorts->size();

  //nPV part
  edm::Handle<std::vector<reco::Vertex >> h_PV;
  iEvent.getByToken(offlinePrimaryVerticesCollectionToken_, h_PV);
  if(!h_PV.isValid()) {
      std::cout << "Missing collection during TreeProducer_AOD : " << offlinePrimaryVerticesCollectionTag_ << " ... skip entry !" << std::endl;
  }
  _nPV = h_PV->size();

  //nmuons part
  edm::Handle<std::vector<edm::FwdPtr<reco::PFCandidate> >> h_muons;
  iEvent.getByToken(muonsCollectionToken_, h_muons);
  if(!h_muons.isValid()) {
      std::cout << "Missing collection during TreeProducer_AOD : " << muonsCollectionTag_ << " ... skip entry !" << std::endl;
  }
  _nmuons = h_muons->size();

  //nelectrons part
  edm::Handle<std::vector<edm::FwdPtr<reco::PFCandidate> >> h_electrons;
  iEvent.getByToken(electronsCollectionToken_, h_electrons);
  if(!h_electrons.isValid()) {
      std::cout << "Missing collection during TreeProducer_AOD : " << electronsCollectionTag_ << " ... skip entry !" << std::endl;
  }
  _nelectrons = h_electrons->size();


  // GLOBAL EVENT INFORMATIONS //

  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  _nEvent = iEvent.id().event();





  _tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
TreeProducer_AOD::beginJob()
{
        // Initialize when class is created
        edm::Service<TFileService> fs ;
        _tree = fs->make <TTree>("SexaQAnalysis","tree");

        // Declare tree's branches
        // Event
        _tree->Branch("nEvent",&_nEvent,"nEvent/I");
        _tree->Branch("nRun",&_nRun,"nRun/I");
        _tree->Branch("nLumi",&_nLumi,"nLumi/I");
        //
        _tree->Branch("nTrack",&_nTrack,"nTrack/I");
        _tree->Branch("nLambdas",&_nLambdas,"nLambdas/I");
        _tree->Branch("nKshorts",&_nKshorts,"nKshorts/I");
        _tree->Branch("nPV",&_nPV,"nPV/I");
        _tree->Branch("nmuons",&_nmuons,"nmuons/I");
        _tree->Branch("nelectrons",&_nelectrons,"nelectrons/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void
TreeProducer_AOD::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
TreeProducer_AOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TreeProducer_AOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TreeProducer_AOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TreeProducer_AOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
TreeProducer_AOD::Init()
{
  _nEvent = _nRun = _nLumi = 0;

  //Tracks
  _nTrack = 0;
  _nLambdas = 0;
  _nKshorts = 0;
  _nPV = 0;
  _nmuons = 0;
  _nelectrons = 0;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer_AOD);
