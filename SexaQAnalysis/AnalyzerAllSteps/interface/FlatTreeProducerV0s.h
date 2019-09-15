#ifndef FlatTreeProducerV0s_h
#define FlatTreeProducerV0s_h
 
#include "AnalyzerAllSteps.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
using namespace edm;
using namespace std; 
class FlatTreeProducerV0s : public edm::EDAnalyzer
 {
  public:
    explicit FlatTreeProducerV0s(edm::ParameterSet const& cfg);
    virtual ~FlatTreeProducerV0s();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    bool IsolationCriterium(reco::Muon muon);
    void FillBranchesV0(const reco::VertexCompositeCandidate * V0, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV, edm::Handle<vector<reco::GenParticle>> h_genParticles, edm::Handle<View<reco::Track>> h_generalTracks, std::string V0Type);    


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
    void InitZ();
    void InitPV();
    void InitBeamspot();

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
    TTree* _tree_Z;   
    TTree* _tree_PV;
    TTree* _tree_beamspot;

    //definition of variables which should go to tree
    //Ks
    std::vector<float> _Ks_mass,_Ks_pt,_Ks_pz,_Ks_Lxy,_Ks_vz,_Ks_eta,_Ks_phi,_Ks_dxy_beamspot,_Ks_dxy_min_PV,_Ks_dxy_PV0,_Ks_dxy_000,_Ks_dz_beamspot,_Ks_dz_min_PV,_Ks_dz_PV0,_Ks_dz_000,_Ks_vz_dz_min_PV;
    std::vector<float> _Ks_deltaRBestMatchingGENParticle,_Ks_trackPair_mindeltaR,_Ks_trackPair_mass,_Ks_Track1Track2_openingsAngle, _Ks_Track1Track2_deltaR,_Ks_Track1_openingsAngle,_Ks_Track2_openingsAngle,_Ks_Track1_deltaR,_Ks_Track2_deltaR;
    std::vector<float> _Ks_daughterTrack1_eta,_Ks_daughterTrack1_phi,_Ks_daughterTrack1_pt,_Ks_daughterTrack1_pz,_Ks_daughterTrack1_dxy,_Ks_daughterTrack1_dz,_Ks_daughterTrack1_charge,_Ks_daughterTrack1_chi2,_Ks_daughterTrack1_ndof,_Ks_daughterTrack1_dxy_beamspot,_Ks_daughterTrack1_dz_beamspot,_Ks_daughterTrack1_dz_min_PV,_Ks_daughterTrack1_dz_PV0,_Ks_daughterTrack1_dz_000;
    std::vector<float> _Ks_daughterTrack2_eta,_Ks_daughterTrack2_phi,_Ks_daughterTrack2_pt,_Ks_daughterTrack2_pz,_Ks_daughterTrack2_charge,_Ks_daughterTrack2_chi2,_Ks_daughterTrack2_ndof,_Ks_daughterTrack2_dxy_beamspot,_Ks_daughterTrack2_dz_beamspot,_Ks_daughterTrack2_dz_min_PV,_Ks_daughterTrack2_dz_PV0,_Ks_daughterTrack2_dz_000;

   //Lambda
    std::vector<float> _Lambda_mass,_Lambda_pt,_Lambda_pz,_Lambda_Lxy,_Lambda_vz,_Lambda_eta,_Lambda_phi,_Lambda_dxy_beamspot,_Lambda_dxy_min_PV,_Lambda_dxy_PV0,_Lambda_dxy_000,_Lambda_dz_beamspot,_Lambda_dz_min_PV,_Lambda_dz_PV0,_Lambda_dz_000,_Lambda_vz_dz_min_PV;
    std::vector<float> _Lambda_deltaRBestMatchingGENParticle,_Lambda_trackPair_mindeltaR,_Lambda_trackPair_mass,_Lambda_Track1Track2_openingsAngle, _Lambda_Track1Track2_deltaR,_Lambda_Track1_openingsAngle,_Lambda_Track2_openingsAngle,_Lambda_Track1_deltaR,_Lambda_Track2_deltaR;
    std::vector<float> _Lambda_daughterTrack1_eta,_Lambda_daughterTrack1_phi,_Lambda_daughterTrack1_pt,_Lambda_daughterTrack1_pz,_Lambda_daughterTrack1_dxy,_Lambda_daughterTrack1_dz,_Lambda_daughterTrack1_charge,_Lambda_daughterTrack1_chi2,_Lambda_daughterTrack1_ndof,_Lambda_daughterTrack1_dxy_beamspot,_Lambda_daughterTrack1_dz_beamspot,_Lambda_daughterTrack1_dz_min_PV,_Lambda_daughterTrack1_dz_PV0,_Lambda_daughterTrack1_dz_000;
    std::vector<float> _Lambda_daughterTrack2_eta,_Lambda_daughterTrack2_phi,_Lambda_daughterTrack2_pt,_Lambda_daughterTrack2_pz,_Lambda_daughterTrack2_charge,_Lambda_daughterTrack2_chi2,_Lambda_daughterTrack2_ndof,_Lambda_daughterTrack2_dxy_beamspot,_Lambda_daughterTrack2_dz_beamspot,_Lambda_daughterTrack2_dz_min_PV,_Lambda_daughterTrack2_dz_PV0,_Lambda_daughterTrack2_dz_000;



    std::vector<float> _Z_mass,_Z_dz_PV_muon1,_Z_dz_PV_muon2,_Z_ptMuMu;
    std::vector<float> _PV_n,_PV0_lxy,_PV0_vz;
    std::vector<float> _beampot_lxy,_beampot_vz;

   




     };

#endif

