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
    //cuts for defining matching between GEN and RECO V0s and AntiS
    static constexpr double deltaRCutV0RECOKs = 0.03;
    static constexpr double deltaLCutV0RECOKs = 2.;
    static constexpr double deltaRCutV0RECOLambda = 0.03;
    static constexpr double deltaLCutV0RECOLambda = 3.;
    static constexpr double deltaRCutRECOAntiS = 0.5;
    static constexpr double deltaLCutInteractionVertexAntiSmin = 2.;

    //definition of the background cuts
    static constexpr double MinLxyCut = 1.9;

    //the location of the center of the beampipe in x and y, from https://arxiv.org/pdf/1807.03289.pdf
    static constexpr double center_beampipe_x = 0.124; //cm
    static constexpr double center_beampipe_y = 0.27; //cm

    //some pdgIds
    const static int pdgIdAntiS = -1020000020;
    const static int pdgIdKs = 310;
    const static int pdgIdAntiLambda = -3122;
    const static int pdgIdAntiProton = -2212;
    const static int pdgIdPosPion = 211;
    const static int pdgIdNegPion = -211;

    //some pdg masses
    static constexpr double pdgMassChargedPion = 0.13957061;

    //some functions to calculate kinematic variables which are used everywhere
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
    std::vector<double> static isTpGrandDaughterAntiS(TrackingParticleCollection const & TPColl, const TrackingParticle& tp);
    double static EventWeightingFactor(double etaAntiS);
    double static PUReweighingFactor(map<double,double>,double MC_PV_vz);

    //definitions of the maps for the PU reweighing. These maps are obtained with /user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/PUReweighing each map is for a certain PU, starting at PU0. Each map has as key the absolute z location and as value the 2D (PU and z location of the valid PVs) reweighing parameter. At the end it includes a vector with all the maps, so when you need the reweighing parameter you need to go to a certain location in the vector, given by the PU, and then find the best mathing z value and get that key's value --> see src/AnalyzerAllSteps.cc for the actual values
    static  map<double,double> mapPU0,mapPU1,mapPU2,mapPU3,mapPU4,mapPU5,mapPU6,mapPU7,mapPU8,mapPU9,mapPU10,mapPU11,mapPU12,mapPU13,mapPU14,mapPU15,mapPU16,mapPU17,mapPU18,mapPU19,mapPU20,mapPU21,mapPU22,mapPU23,mapPU24,mapPU25,mapPU26,mapPU27,mapPU28,mapPU29,mapPU30,mapPU31,mapPU32,mapPU33,mapPU34,mapPU35,mapPU36,mapPU37,mapPU38,mapPU39,mapPU40,mapPU41,mapPU42,mapPU43,mapPU44,mapPU45,mapPU46,mapPU47,mapPU48,mapPU49,mapPU50,mapPU51,mapPU52,mapPU53,mapPU54,mapPU55,mapPU56,mapPU57,mapPU58,mapPU59,mapPU60; 
    static vector<map<double, double>> v_mapPU; 

    };

#endif

