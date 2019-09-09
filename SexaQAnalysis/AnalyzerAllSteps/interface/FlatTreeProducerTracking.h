#ifndef AnalyzerAllSteps_h
#define AnalyzerAllSteps_h
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "AnalyzerAllSteps.h"
using namespace edm;
using namespace std; 
class FlatTreeProducerTracking : public edm::EDAnalyzer
 {
  public:
    explicit FlatTreeProducerTracking(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    virtual void endJob();


    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    virtual ~FlatTreeProducerTracking();


    void FillTreesTracks(const TrackingParticle& tp, TVector3 beamspot, int nPVs, const reco::Track *matchedTrackPointer, bool matchingTrackFound, std::vector<double> tpIsGrandDaughterAntiS);



 //   void FillHistosNonAntiSTracksAll(const TrackingParticle& tp, TVector3 beamspot, int nPVs, int matchedTrackQuality);
    void FillHistosAntiSTracks(const TrackingParticle& tp, TVector3 beamspot, TrackingParticleCollection const & TPColl, edm::Handle<TrackingParticleCollection> h_TP, edm::Handle< reco::TrackToTrackingParticleAssociator> h_trackAssociator, edm::Handle<View<reco::Track>> h_generalTracks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L);
    void FillHistosAntiSKsDaughterTracksRECO(const TrackingParticle& tp, TVector3 beamspot);
    void FillHistosAntiSKsDaughterTracksAll(const TrackingParticle& tp, TVector3 beamspot);
    void FillHistosAntiSAntiLAntiProtonDaughterTracksRECO(const TrackingParticle& tp, TVector3 beamspot); 
    void FillHistosAntiSAntiLAntiProtonDaughterTracksAll(const TrackingParticle& tp, TVector3 beamspot); 
    void FillHistosAntiSAntiLPosPionDaughterTracksRECO(const TrackingParticle& tp, TVector3 beamspot); 
    void FillHistosAntiSAntiLPosPionDaughterTracksAll(const TrackingParticle& tp, TVector3 beamspot); 
    void FillMajorEfficiencyPlot(std::vector<bool>granddaughterTrackMatched,const reco::Track *matchedTrackPointerKsPosPion,const reco::Track *matchedTrackPointerKsNegPion,const reco::Track *matchedTrackPointerAntiLPosPion,const reco::Track *matchedTrackPointerAntiLNegProton, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L);

    void FillHistosGENAntiS(const reco::Candidate * Candidate, TVector3 beamspot);
    void FillHistosGENInteractingAntiS(const reco::Candidate * Candidate, TVector3 beamspot);
    void FillHistosGENKsNonAntiS(const reco::Candidate *Candidate, TVector3 beamspot );
    void FillHistosGENAntiLambdaNonAntiS(const reco::Candidate *Candidate,TVector3 beamspot );
    void FillHistosGENKsAntiS(const reco::Candidate *Candidate, TVector3 beamspot );
    void FillHistosGENAntiLambdaAntiS(const reco::Candidate *Candidate,TVector3 beamspot );
    void RecoEvaluationKsNonAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks,  TVector3 beamspot);
    void RecoEvaluationKsAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks,  TVector3 beamspot);
    void RecoEvaluationAntiLambdaNonAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L,  TVector3 beamspot);
    void RecoEvaluationAntiLambdaAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L,  TVector3 beamspot);
    void RecoEvaluationAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands, TVector3 beamspot, TVector3 beamspotVariance, edm::Handle<vector<reco::Vertex>> h_offlinePV);


  private:
    //---- configurable parameters --------
    bool m_lookAtAntiS;

    void InitTracking();
    
    edm::Service<TFileService> m_fs;
 
    edm::InputTag m_bsTag;
    edm::InputTag m_offlinePVTag;
    edm::InputTag m_genParticlesTag_GEN;
    edm::InputTag m_genParticlesTag_SIM_GEANT;
    edm::InputTag m_generalTracksTag;
    edm::InputTag m_sCandsTag;
    edm::InputTag m_V0KsTag;
    edm::InputTag m_V0LTag;
    edm::InputTag m_trackAssociatorTag;
    edm::InputTag m_TPTag;
//    edm::InputTag m_PileupInfoTag;


    edm::EDGetTokenT<reco::BeamSpot> m_bsToken;
    edm::EDGetTokenT<vector<reco::Vertex>> m_offlinePVToken;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_GEN;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_SIM_GEANT;
    //edm::EDGetTokenT<vector<reco::Track>> m_generalTracksToken;
    edm::EDGetTokenT<View<reco::Track>> m_generalTracksToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_sCandsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0KsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0LToken;
    edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator>  m_trackAssociatorToken;
    edm::EDGetTokenT<vector<TrackingParticle> > m_TPToken;
//    edm::EDGetTokenT<vector<PileupSummaryInfo> > m_PileupInfoToken;
   
    TTree* _tree_tracks;   
    std::vector<float> _tp_pt,_tp_eta,_tp_phi,_tp_pz,_tp_Lxy_beamspot,_tp_vz_beamspot,_tp_dxy_beamspot,_tp_dz_beamspot,_tp_etaOfGrandMotherAntiS;
    std::vector<int> _tp_numberOfTrackerLayers,_tp_charge,_tp_reconstructed,_tp_isAntiSTrack;
    std::vector<float> _matchedTrack_pt,_matchedTrack_eta,_matchedTrack_phi,_matchedTrack_pz,_matchedTrack_chi2,_matchedTrack_ndof,_matchedTrack_dxy_beamspot,_matchedTrack_dz_beamspot;
    std::vector<int> _matchedTrack_trackQuality,_matchedTrack_charge,_matchedTrack_isLooper;



     };

#endif

