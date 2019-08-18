import FWCore.ParameterSet.Config as cms

lambdaKshortFilter = cms.EDFilter(
    'LambdaKshortFilter',
    lambdaCollection = cms.InputTag("generalV0Candidates","Lambda"),
    kshortCollection = cms.InputTag("generalV0Candidates","Kshort"),
    #genCollection    = cms.InputTag("genCollection"),
    genCollection = cms.InputTag("genParticles","","HLT"),
    isData = cms.bool(True),
    minNrLambda = cms.uint32(1),
    minNrKshort = cms.uint32(1),
    minPtLambda = cms.double(0),
    minPtKshort = cms.double(0),
    maxEtaLambda = cms.double(99999), #no eta cut any more
    maxEtaKshort = cms.double(99999), #no eta cut any more
    minMassLambda = cms.double(1.106), # -3sigma arXiv:1102.4282
    minMassKshort = cms.double(0.473), # +3sigma arXiv:1102.4282
    maxMassLambda = cms.double(1.126), # -3sigma arXiv:1102.4282
    maxMassKshort = cms.double(0.522), # +3sigma arXiv:1102.4282
    checkLambdaDaughters = cms.bool(False),
    prescaleFalse = cms.uint32(0) # 0 means no prescale, reject all
)
