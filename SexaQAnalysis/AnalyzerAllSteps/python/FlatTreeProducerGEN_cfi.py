import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * 
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import * 
FlatTreeProducerGEN = cms.EDAnalyzer('FlatTreeProducerGEN',
    lookAtAntiS = cms.untracked.bool(False),
    runningOnData = cms.untracked.bool(False),
    beamspot = cms.InputTag("offlineBeamSpot"),
    genCollection_GEN =  cms.InputTag("genParticles","","GEN")
)
