#ifndef FlatTreeProducerV0s_h
#define FlatTreeProducerV0s_h
 
#include "AnalyzerAllSteps.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/LorentzVector.h"

using namespace edm;
using namespace std; 
class FlatTreeProducerV0s : public edm::EDAnalyzer
 {
  public:
    explicit FlatTreeProducerV0s(edm::ParameterSet const& cfg);
    virtual ~FlatTreeProducerV0s();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void FillBranchesKs(const reco::VertexCompositeCandidate * Ks, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV);    
    void FillBranchesLambda(const reco::VertexCompositeCandidate * Lambda, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV);    

  private:
    bool m_lookAtAntiS;
    bool m_runningOnData; 

    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    void InitKs();
    void InitLambda();

    edm::Service<TFileService> m_fs;
 
    edm::InputTag m_bsTag;
    edm::InputTag m_offlinePVTag;
    edm::InputTag m_genParticlesTag_GEN;
    edm::InputTag m_genParticlesTag_SIM_GEANT;
    edm::InputTag m_generalTracksTag;
    edm::InputTag m_sCandsTag;
    edm::InputTag m_V0KsTag;
    edm::InputTag m_V0LTag;
    edm::InputTag m_muonsTag;
    edm::InputTag m_jetsTag;


    edm::EDGetTokenT<reco::BeamSpot> m_bsToken;
    edm::EDGetTokenT<vector<reco::Vertex>> m_offlinePVToken;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_GEN; 
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_SIM_GEANT; 
    edm::EDGetTokenT<View<reco::Track>> m_generalTracksToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_sCandsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0KsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0LToken;
    edm::EDGetTokenT<vector<reco::Muon> > m_muonsToken;
    edm::EDGetTokenT<vector<reco::PFJet> > m_jetsToken;

   
    TTree* _tree_Ks;   
    TTree* _tree_Lambda;   

    //definition of variables which should go to tree
    std::vector<float> _Ks_mass,_Ks_pt,_Ks_pz,_Ks_Lxy,_Ks_vz,_Ks_eta,_Ks_phi,_Ks_dxy,_Ks_dz,_Ks_dz_min;
    std::vector<float> _Lambda_mass,_Lambda_pt,_Lambda_pz,_Lambda_Lxy,_Lambda_vz,_Lambda_eta,_Lambda_phi,_Lambda_dxy,_Lambda_dz,_Lambda_dz_min;

     };

#endif

