#ifndef FlatTreeProducerGENSIM_h
#define FlatTreeProducerGENSIM_h
 
#include "AnalyzerAllSteps.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
using namespace edm;
using namespace std; 
class FlatTreeProducerGENSIM : public edm::EDAnalyzer
 {
  public:
    explicit FlatTreeProducerGENSIM(edm::ParameterSet const& cfg);
    virtual ~FlatTreeProducerGENSIM();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    bool FillBranchesGENAntiS(const reco::Candidate  * genParticle, TVector3 beamspot, TVector3 beamspotVariance, vector<vector<float>> v_antiS_momenta_and_itt, edm::Handle<TrackingParticleCollection>  h_TP, unsigned int nGoodPV);

  private:
    int nTotalGENS=0;
    int nTotalUniqueGenS=0;
    int nTotalGiving2DaughtersGENS=0;
    int nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS=0;
    int nTotalCorrectGENS=0;
    int nTotalCorrectGENSInteractingInBeampipe=0;
    int nTotalGENSPosEta=0;
    int nTotalGENSNegEta=0;
    int nTotalRecoconstructableGENS_negEta=0;
    int nTotalRecoconstructableGENS_posEta=0;
    bool m_lookAtAntiS;
    bool m_runningOnData; 

    double  nTotalUniqueGenS_weighted =0.;
    double  nTotalGiving2DaughtersGENS_weighted = 0.;
    double  nTotalGivingCorrectDaughtersAnd4GrandDaughtersGENS_weighted = 0.;
    double  nTotalCorrectGENS_weighted = 0.;
    double  nTotalCorrectGENSInteractingInBeampipe_weighted = 0.;
    double  nTotalGENSPosEta_weighted = 0.;
    double  nTotalGENSNegEta_weighted = 0.;
    double  nTotalCorrectGENS_Reconstructable_weighted = 0.;

    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    void Init();

    edm::Service<TFileService> m_fs;
 
    edm::InputTag m_bsTag;
    edm::InputTag m_offlinePVTag;
    edm::InputTag m_genParticlesTag_GEN;
    edm::InputTag m_genParticlesTag_SIM_GEANT;
    edm::InputTag m_TPTag;

    edm::EDGetTokenT<reco::BeamSpot> m_bsToken;
    edm::EDGetTokenT<vector<reco::Vertex>> m_offlinePVToken;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_GEN; 
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_SIM_GEANT; 
    edm::EDGetTokenT<vector<TrackingParticle> > m_TPToken;   

    //Two trees: one with info on all the AntiS and one with info on the AntiS which give rise to our signal (so the correct granddaughters)
    TTree* _treeAllAntiS;
    TTree* _tree;   

    //definition of variables which should go to _treeAllAntiS
    std::vector<float> _S_eta_all,_S_event_weighting_factor_all,_S_event_weighting_factor_PU_all,_S_vz_creation_vertex_all,_S_pt_all,_S_pz_all;
    std::vector<int> _S_reconstructable_all, _S_nGoodPV_all;

    //definition of variables which should go to _tree
    std::vector<float> _S_n_loops;
    std::vector<float> _S_charge;
    std::vector<float> _S_nGoodPV;
    std::vector<float> _S_event_weighting_factor,_S_event_weighting_factor_PU;
    std::vector<float> _S_lxy_interaction_vertex, _S_lxy_interaction_vertex_beamspot, _S_lxy_interaction_vertex_beampipeCenterData, _S_lxyz_interaction_vertex,  _S_error_lxy_interaction_vertex,_S_mass,_S_Mt,_S_chi2_ndof;
    std::vector<float> _n_M, _n_p;
    std::vector<float> _S_daughters_deltaphi,_S_daughters_deltaeta,_S_daughters_openingsangle,_S_Ks_openingsangle,_S_Lambda_openingsangle,_S_sumDaughters_openingsangle,_S_sumDaughters_deltaPhi,_S_sumDaughters_deltaEta,_S_sumDaughters_deltaR,_S_daughters_DeltaR,_S_eta,_Ks_eta,_Lambda_eta;
    std::vector<float> _S_dxy,_Ks_dxy,_Lambda_dxy;
    std::vector<float> _S_dxy_over_lxy,_Ks_dxy_over_lxy,_Lambda_dxy_over_lxy;
    std::vector<float> _S_dz,_Ks_dz,_Lambda_dz,_S_dz_min,_Ks_dz_min,_Lambda_dz_min;
    std::vector<float> _deltaR_sumDaughterMomenta_antiSMomentum,_Ks_openings_angle_displacement_momentum,_Lambda_openings_angle_displacement_momentum; 
    std::vector<float> _S_pt,_Ks_pt,_Lambda_pt;
    std::vector<float> _S_pz,_Ks_pz,_Lambda_pz;
    std::vector<float> _S_vx_interaction_vertex,_S_vy_interaction_vertex,_S_vz_interaction_vertex;
    std::vector<float> _S_vx,_S_vy,_S_vz;

    std::vector<float> _GEN_Ks_daughter0_px,_GEN_Ks_daughter0_py,_GEN_Ks_daughter0_pz,_GEN_Ks_daughter0_pt,_GEN_Ks_daughter0_eta,_GEN_Ks_daughter0_phi;
    std::vector<float> _GEN_Ks_daughter0_vx,_GEN_Ks_daughter0_vy,_GEN_Ks_daughter0_vz,_GEN_Ks_daughter0_lxy,_GEN_Ks_daughter0_dxy,_GEN_Ks_daughter0_dz,_GEN_Ks_daughter0_openings_angle_displacement_momentum;
    std::vector<float> _GEN_Ks_daughter1_px,_GEN_Ks_daughter1_py,_GEN_Ks_daughter1_pz,_GEN_Ks_daughter1_pt,_GEN_Ks_daughter1_eta,_GEN_Ks_daughter1_phi;
    std::vector<float> _GEN_Ks_daughter1_vx,_GEN_Ks_daughter1_vy,_GEN_Ks_daughter1_vz,_GEN_Ks_daughter1_lxy,_GEN_Ks_daughter1_dxy,_GEN_Ks_daughter1_dz,_GEN_Ks_daughter1_openings_angle_displacement_momentum;
    std::vector<float> _GEN_AntiLambda_AntiProton_px,_GEN_AntiLambda_AntiProton_py,_GEN_AntiLambda_AntiProton_pz,_GEN_AntiLambda_AntiProton_pt,_GEN_AntiLambda_AntiProton_eta,_GEN_AntiLambda_AntiProton_phi;
    std::vector<float> _GEN_AntiLambda_AntiProton_vx,_GEN_AntiLambda_AntiProton_vy,_GEN_AntiLambda_AntiProton_vz,_GEN_AntiLambda_AntiProton_lxy,_GEN_AntiLambda_AntiProton_dxy,_GEN_AntiLambda_AntiProton_dz,_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum;
    std::vector<float> _GEN_AntiLambda_Pion_px,_GEN_AntiLambda_Pion_py,_GEN_AntiLambda_Pion_pz,_GEN_AntiLambda_Pion_pt,_GEN_AntiLambda_Pion_phi,_GEN_AntiLambda_Pion_eta;
    std::vector<float> _GEN_AntiLambda_Pion_vx,_GEN_AntiLambda_Pion_vy,_GEN_AntiLambda_Pion_vz,_GEN_AntiLambda_Pion_lxy,_GEN_AntiLambda_Pion_dxy,_GEN_AntiLambda_Pion_dz,_GEN_AntiLambda_Pion_openings_angle_displacement_momentum;
    
    std::vector<float> _GEN_Ks_daughter0_numberOfTrackerLayers,_GEN_Ks_daughter1_numberOfTrackerLayers,_GEN_AntiLambda_AntiProton_numberOfTrackerLayers,_GEN_AntiLambda_Pion_numberOfTrackerLayers;
    std::vector<float> _GEN_Ks_daughter0_numberOfTrackerHits,_GEN_Ks_daughter1_numberOfTrackerHits,_GEN_AntiLambda_AntiProton_numberOfTrackerHits,_GEN_AntiLambda_Pion_numberOfTrackerHits;
     };

#endif

