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
    TrackingParticles = cms.InputTag("mix","MergedTrackTruth")
#    PileupInfo = cms.InputTag("addPileupInfo","","HLT")
)
