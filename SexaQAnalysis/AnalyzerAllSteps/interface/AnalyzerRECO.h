#ifndef AnalyzerRECO_h
#define AnalyzerRECO_h
 
#include "AnalyzerAllSteps.h"
using namespace edm;
using namespace std; 
class AnalyzerRECO : public edm::EDAnalyzer
 {
  public:
    explicit AnalyzerRECO(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~AnalyzerRECO();


    void FillHistosPV(reco::Vertex PrimVertex, TVector3 beamspot);
    void FillHistosRECOKs(const reco::VertexCompositeCandidate * RECOKs, TVector3 beamspot, edm::Handle<vector<reco::Vertex>> h_offlinePV);
    void FillHistosRECOLambda(const reco::VertexCompositeCandidate * RECOAntiLambda, TVector3 beamspot,edm::Handle<vector<reco::Vertex>> h_offlinePV);
    void FillHistosRECOTracks(const reco::Track *track, TVector3 beamspot);
    void FillHistosRECOAntiS(const reco::VertexCompositeCandidate * RECOAntiS, TVector3 beamspot,TVector3 beamspotVariance, int eventId, edm::Handle<vector<reco::Vertex>> h_offlinePV);


  private:
    //---- configurable parameters --------
    bool m_lookAtAntiS;
    int m_nEvent, m_nRun, m_nLumi;
    
    edm::Service<TFileService> m_fs;
 
    edm::InputTag m_bsTag;
    edm::InputTag m_offlinePVTag;
    edm::InputTag m_generalTracksTag;
    edm::InputTag m_sCandsTag;
    edm::InputTag m_V0KsTag;
    edm::InputTag m_V0LTag;


    edm::EDGetTokenT<reco::BeamSpot> m_bsToken;
    edm::EDGetTokenT<vector<reco::Vertex>> m_offlinePVToken;
    edm::EDGetTokenT<View<reco::Track>> m_generalTracksToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_sCandsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0KsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0LToken;
   
 
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

