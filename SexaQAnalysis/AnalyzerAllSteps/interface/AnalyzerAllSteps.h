#ifndef My_azimuthal_MC_h
#define My_azimuthal_MC_h
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>
#include <math.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TString.h>

#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
using namespace edm;
using namespace std;

#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//matching on hits specific:
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
 
class AnalyzerAllSteps : public edm::EDAnalyzer
 {
  public:
    static constexpr double deltaRCutV0RECOKs = 0.1;
    static constexpr double deltaRCutV0RECOLambda = 0.1;
    static constexpr double deltaRCutRECOAntiS = 0.1;


    //definition of the background cuts
    static constexpr double MinLxyCut = 1.9;
    static constexpr double MaxErrorLxyCut = 0.1;
    static constexpr double MinErrorMassCut = 0.;
    static constexpr double MaxNormChi2Cut = 4.;

    static constexpr double MinAbsDeltaPhiDaughtersCut = 1.;
    static constexpr double MaxAbsDeltaPhiDaughtersCut = 3.;
    static constexpr double MaxOpeningsAngleCut = 1.2;
    static constexpr double MaxAbsDeltaEtaDaughCut = 1.5;
    static constexpr double MinDxyOverLxyCut = 0.;
    static constexpr double MaxDxyOverLxyCut = 0.1;

    static constexpr double DxyKsExclusionRangeMinCut = 0.;
    static constexpr double DxyKsExclusionRangeMaxCut = 0.2;
    static constexpr double DxyAntiLExclusionRangeMinCut = 0.;
    static constexpr double DxyAntiLExclusionRangeMaxCut = 0.2;

    static constexpr double MinLambdaPtCut = 1.2;//zero background: 1.5
    static constexpr double dzAntiSPVminCut = 2.;
    static constexpr double vzAntiSInteractionVertexCut = 4.;//zero background: 5
    static constexpr double antiSEtaCut = 1.3;//zero background: 1.5
	

    //some pdgIds
    const static int pdgIdAntiS = -1020000020;
    const static int pdgIdKs = 310;
    const static int pdgIdAntiLambda = -3122;
    const static int pdgIdAntiProton = -2212;
    const static int pdgIdPosPion = 211;
    const static int pdgIdNegPion = -211;

    double static openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2);
    double static deltaR(double phi1, double eta1, double phi2, double eta2);
    double static lxy(TVector3 v1, TVector3 v2);
    double static lxyz(TVector3 v1, TVector3 v2);
    TVector3 static PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    TVector3 static vec_dxy_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double static dxy_signed_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double static dxyz_signed_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in);
    double static std_dev_lxy(double vx, double vy, double vx_var, double vy_var, double bx_x, double bx_y, double bx_x_var, double bx_y_var);
    double static XYpointingAngle(const reco::Candidate  * particle,TVector3 beamspot);
    double static CosOpeningsAngle(TVector3 vec1, TVector3 vec2);
    double static dz_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in);
    TVector3 static dz_line_point_min(TVector3 Point_line_in, TVector3 Vector_along_line_in, edm::Handle<vector<reco::Vertex>> h_offlinePV);
    double static sgn(double input);
    int static getDaughterParticlesTypes(const reco::Candidate * genParticle);
    int static trackQualityAsInt(const reco::Track *track);
     };

#endif

