import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * 
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import * 
FlatTreeProducerTracking= cms.EDAnalyzer('FlatTreeProducerTracking',
    beamspot = cms.InputTag("offlineBeamSpot","",""),
    offlinePV = cms.InputTag("offlinePrimaryVertices","",""),
    genCollection_GEN =  cms.InputTag("genParticles","","GEN"),
    genCollection_SIM_GEANT =  cms.InputTag("genParticlesPlusGEANT","","SIM"),
    generalTracksCollection =  cms.InputTag("generalTracks","",""),
    sexaqCandidates = cms.InputTag("lambdaKshortVertexFilter", "sParticles",""),
    V0KsCollection = cms.InputTag("generalV0Candidates","Kshort",""), #can also be SEXAQ
    V0LCollection = cms.InputTag("generalV0Candidates","Lambda",""), #can also be SEXAQ
    trackAssociators = cms.InputTag("quickTrackAssociatorByHits"),
    TrackingParticles = cms.InputTag("mix","MergedTrackTruth"),
#    PileupInfo = cms.InputTag("addPileupInfo","","HLT")

    #################
    #for the V0Fitter
    #################
   useVertex = cms.bool(False),
   # which V0s to reconstruct
   doKShorts = cms.bool(True),
   doLambdas = cms.bool(True),
  
    # which vertex fitting algorithm to use
   # True -> KalmanVertexFitter (recommended)
   # False -> AdaptiveVertexFitter (not recommended)
    vertexFitter = cms.bool(True),

   # use the refitted tracks returned from the KVF for V0Candidate kinematics
   # this is automatically set to False if using the AdaptiveVertexFitter
   useRefTracks = cms.bool(True),

   # -- cuts on initial track collection --
   # Track normalized Chi2 <
   tkChi2Cut = cms.double(10.),
   # Number of valid hits on track >=
   tkNHitsCut = cms.int32(7),
   # Pt of track >
   tkPtCut = cms.double(0.35),
   # Track impact parameter significance >
   tkIPSigXYCut = cms.double(2.),
   tkIPSigZCut = cms.double(-1.),

   # -- cuts on the vertex --
   # Vertex chi2 <
   vtxChi2Cut = cms.double(15.),
   # XY decay distance significance >
   vtxDecaySigXYCut = cms.double(10.),
   # XYZ decay distance significance >
   vtxDecaySigXYZCut = cms.double(-1.),

   # -- miscellaneous cuts --
   # POCA distance between tracks <
   tkDCACut = cms.double(2.),
   # invariant mass of track pair - assuming both tracks are charged pions <
   mPiPiCut = cms.double(0.6),
   # check if either track has a hit radially inside the vertex position minus this number times the sigma of the vertex fit
   # note: Set this to -1 to disable this cut, which MUST be done if you want to run V0Producer on the AOD track collection!
   innerHitPosCut = cms.double(4.),
   # cos(angleXY) between x and p of V0 candidate >
   cosThetaXYCut = cms.double(0.9998),
   # cos(angleXYZ) between x and p of V0 candidate >
   cosThetaXYZCut = cms.double(-2.),

   # -- cuts on the V0 candidate mass --
   # V0 mass window +- pdg value
   kShortMassCut = cms.double(0.07),
   lambdaMassCut = cms.double(0.05)

)
