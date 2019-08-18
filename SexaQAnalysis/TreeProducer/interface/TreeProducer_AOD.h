// Description: EDAnalyzer produce flat trees from AOD for HexaAnalysis

// C++ lib
#include <vector>

// ROOT
#include "TTree.h"
#include "TMatrixD.h"
//#include "TLorentzVector.h"
#include "TPRegexp.h"

// CMSSW standard lib
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// CMSSW specific lib
//#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

// others
using namespace std;
int verbose=1;

//
// class declaration
//

class TreeProducer_AOD : public edm::EDAnalyzer {
 public:
  explicit TreeProducer_AOD(const edm::ParameterSet&);
  ~TreeProducer_AOD();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static bool ptSorter(const reco::Track & i, const reco::Track & j);

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void Init();

  // ----------member data ---------------------------
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
 

//   GlobalPoint vertexPosition;

  // Tree and its branches
  TTree* _tree;

  // Global quantities
  int _nEvent, _nRun, _nLumi;

  int _nTrack, _nLambdas, _nKshorts, _nPV, _nmuons,  _nelectrons;

};

namespace reco {
	template<typename T>
	class RecoPtSorter{
	public:
		bool operator ()(const T & i, const T & j) const {
			return (i->pt() > j->pt());
		}
	};
}

//
// constants, enums and typedefs
//
namespace reco {
  typedef std::vector<Track> TrackCollection;
  typedef edm::Ref<TrackCollection> TrackRef;
}
//
// static data member definitions
//
