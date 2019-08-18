import FWCore.ParameterSet.Config as cms

massFilter = cms.EDFilter(
    'MassFilter',
    lambdakshortCollection = cms.InputTag("lambdaKshortVertexFilter", "sParticles"),
    minMass = cms.double(-10000),
    maxMass = cms.double(10000),
    targetMass = cms.double(0),  # neutron mass 0.939565
    prescaleFalse = cms.uint32(0) # 0 means no prescale, reject all
)
