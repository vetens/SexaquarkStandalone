import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

runningOnData = False #should be False as you will only run this on MC  
lookAtAntiS =   True

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
#again, not sure why this would be here since this is MC only but I will update if needed
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#inlist = open("gensimlist_trial3.txt", "r")
inlist = open("EDM_SkimmedTrial3.txt", "r")
process.source = cms.Source("PoolSource",
	#fileNames = cms.untracked.vstring(options.inputFiles),
	fileNames = cms.untracked.vstring(*(inlist.readlines())),
  duplicateCheckMode = cms.untracked.string ("noDuplicateCheck")
)

# testing fix to beamspot issue
#process.load('RecoVertex.BeamSpotProducer.BeamSpot_cfi')
# checking what gen info there is
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printTree = cms.EDAnalyzer("ParticleListDrawer",
#  maxEventsToPrint = cms.untracked.int32(10),
#  printVertex = cms.untracked.bool(False),
#  printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
#  src = cms.InputTag("genParticlesPlusGEANT")
#)

#process.Troubleshoot = cms.Path(process.offlineBeamSpot * process.printTree)
#without beamspot fix (should now be implemented at gensim level)
#process.Troubleshoot = cms.Path(process.printTree)


process.load("SexaQAnalysis.AnalyzerAllSteps.FlatTreeProducerGENSIM_cfi")
process.FlatTreeProducerGENSIM.runningOnData = runningOnData
process.FlatTreeProducerGENSIM.lookAtAntiS = lookAtAntiS
process.flattreeproducer = cms.Path(process.FlatTreeProducerGENSIM)

process.p = cms.Schedule(
# testing fix to beamspot issue & checking what gen info there is
#  process.Troubleshoot,

  process.flattreeproducer
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

#print process.dumpPython()
