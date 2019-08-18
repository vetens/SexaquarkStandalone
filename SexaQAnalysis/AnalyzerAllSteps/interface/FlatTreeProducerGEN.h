#ifndef FlatTreeProducerGEN_h
#define FlatTreeProducerGEN_h
 
#include "AnalyzerAllSteps.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
using namespace edm;
using namespace std; 
class FlatTreeProducerGEN : public edm::EDAnalyzer
 {
  public:
    explicit FlatTreeProducerGEN(edm::ParameterSet const& cfg);
    virtual ~FlatTreeProducerGEN();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void FillBranchesGENAntiS(const reco::Candidate  * genParticle, TVector3 beamspot, TVector3 beamspotVariance);

  private:
    int nTotalGENS=0;
    int nTotalGENSNegEta=0;
    int nTotalGENSPosEta=0;
    bool m_lookAtAntiS;
    bool m_runningOnData; 

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
    edm::InputTag m_genParticlesTag_GEN;

    edm::EDGetTokenT<reco::BeamSpot> m_bsToken;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_GEN; 
 
    TTree* _tree;   

    //definition of variables which should go to tree
    std::vector<float> _S_charge;
    std::vector<float> _S_dxy;
    std::vector<float> _S_eta;
    std::vector<float> _S_dz,_S_dz_min;
    std::vector<float> _S_pt;
    std::vector<float> _S_pz;
    std::vector<float> _S_vx,_S_vy,_S_vz;


     };

#endif

