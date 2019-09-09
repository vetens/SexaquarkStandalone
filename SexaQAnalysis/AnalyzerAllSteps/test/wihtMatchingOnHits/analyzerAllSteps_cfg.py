import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

runningOnData = False #this flag is used only for choosing to run both GEN and RECO analyzers (if flag is False) or only the RECO analyzer (when flag is True). It needs to be true for data, because the GEN analyzer cannot run on data as it requires some special collections
lookAtAntiS = True  #This flag should be False if you are running on data unless you want to unblind. If you are running on MC it should be True as you want to see the signal.

options = VarParsing ('analysis')
options.parseArguments()
## data or MC options
options.register(
	'isData',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC')

options.register(
	'maxEvts',-1,VarParsing.multiplicity.singleton,VarParsing.varType.int,
	'flag to indicate max events to process')
	
options.isData==True

process = cms.Process("SEXAQDATAANA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

if(options.isData==True):
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_miniAODv2_v1', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(options.inputFiles),
  duplicateCheckMode = cms.untracked.string ("noDuplicateCheck")
)

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')


### validation-specific includes
process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
#process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("Validation.RecoTrack.cuts_cff")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("Validation.Configuration.postValidation_cff")
process.quickTrackAssociatorByHits.SimToRecoDenominator = 'reco'

process.quickTrackAssociatorByHits.useClusterTPAssociation = True
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")


#TrackRECO on hits:
process.validation = cms.Sequence(
    process.tpClusterProducer *
    process.quickTrackAssociatorByHits
)



# Analyzers
#only run this analyzer when on data, cause you do not have the collections which are required for the top one which is for GEN. For GEN also run this one
process.load("SexaQAnalysis.AnalyzerAllSteps.AnalyzerRECO_cfi")
process.AnalyzerRECO.lookAtAntiS = lookAtAntiS
if(runningOnData==True):
	process.analyzerallsteps = cms.Path(process.AnalyzerRECO)

#run this analyer only on GEN:
process.load("SexaQAnalysis.AnalyzerAllSteps.AnalyzerGEN_cfi")
if(runningOnData==False):
	process.analyzerallsteps = cms.Path(process.validation*process.AnalyzerGEN*process.AnalyzerRECO)


process.p = cms.Schedule(
  process.analyzerallsteps
)


# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string(options.outputFile)
)


#Keep edm output file --> used in the analyzer
#process.out = cms.OutputModule("PoolOutputModule",
#  outputCommands = cms.untracked.vstring(
#     'keep *'
#  ),
#   fileName = cms.untracked.string("AOD_test_matchingHits.root"),
# # SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring('p') )
#)

#process.output_step = cms.EndPath(process.out)

