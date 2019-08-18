import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * 
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import * 
AnalyzerGEN = cms.EDAnalyzer('AnalyzerGEN',
    beamspot = cms.InputTag("offlineBeamSpot"),
    offlinePV = cms.InputTag("offlinePrimaryVertices","","RECO"),
    genCollection_GEN =  cms.InputTag("genParticles","","GEN"),
    genCollection_SIM_GEANT =  cms.InputTag("genParticlesPlusGEANT","","SIM"),
    generalTracksCollection =  cms.InputTag("generalTracks","","RECO"),
    sexaqCandidates = cms.InputTag("lambdaKshortVertexFilter", "sParticles","SEXAQ"),
    V0KsCollection = cms.InputTag("generalV0Candidates","Kshort","RECO"),
    V0LCollection = cms.InputTag("generalV0Candidates","Lambda","RECO"),
    trackAssociators = cms.InputTag("quickTrackAssociatorByHits"),
    TrackingParticles = cms.InputTag("mix","MergedTrackTruth")
#    PileupInfo = cms.InputTag("addPileupInfo","","HLT")
)
