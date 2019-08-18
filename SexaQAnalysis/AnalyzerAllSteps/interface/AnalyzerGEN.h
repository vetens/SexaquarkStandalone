#ifndef AnalyzerAllSteps_h
#define AnalyzerAllSteps_h
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "AnalyzerAllSteps.h"
using namespace edm;
using namespace std; 
class AnalyzerGEN : public edm::EDAnalyzer
 {
  public:
    explicit AnalyzerGEN(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~AnalyzerGEN();

    bool isTpGrandDaughterAntiS(TrackingParticleCollection const & TPColl, const TrackingParticle& tp);

    void FillHistosNonAntiSTracksRECO(const TrackingParticle& tp, TVector3 beamspot, int nPVs, int matchedTrackQuality);
    void FillHistosNonAntiSTracksAll(const TrackingParticle& tp, TVector3 beamspot, int nPVs, int matchedTrackQuality);
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
    int m_nEvent, m_nRun, m_nLumi;
    
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
   
 
    int verbose=1;
    
    TString a = ""; //"WjetsMC" "ZeroBias" "MET" "MinBiasMC" "SingleMuon"
    TString b = a + "";
    
    //--------- Histogram Declaration --------------------//
    std::map<TString, TH1F *> histos_th1f;
    std::map<TString, TH2I *> histos_th2i;
    std::map<TString, TH2F *> histos_th2f;
    std::map<TString, TH3F *> histos_th3f;
    std::map<TString, TEfficiency *> histos_teff;
    std::map<TString, TProfile *> histos_TProfile;

     };

#endif

