# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --filein file:SUS-RunIISummer16DR80Premix-00068_step1.root --fileout file:SUS-RunIISummer16DR80Premix-00068.root --mc --eventcontent AODSIM --runUnscheduled --datatier AODSIM --conditions 80X_mcRun2_asymptotic_2016_TrancheIV_v6 --step RAW2DIGI,RECO,EI --nThreads 4 --era Run2_2016 --python_filename SUS-RunIISummer16DR80Premix-00068_2_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 82
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2016)
options = VarParsing ('analysis')
options.parseArguments()
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("RecoVertex.V0Producer.generalV0Candidates_cfi")
process.generalV0Candidates.innerHitPosCut = -1
process.generalV0Candidates.cosThetaXYCut = -1
process.generalV0Candidates.kShortMassCut = 0.03
process.generalV0Candidates.lambdaMassCut = 0.015


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)

# Input source
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(options.inputFiles),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string ("noDuplicateCheck")
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:82'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.AODSIMEventContent.outputCommands
)

#process.genParticles.src = cms.InputTag("source")

# Additional output definition
process.AODSIMoutput.outputCommands += (
"keep *_genParticlesPlusGEANT_*_*",
"keep *_g4SimHits_*_*",
"keep *_simSiPixelDigis_*_*",
"keep *_simMuonRPCDigis_*_*",
"keep *_simSiStripDigis_*_*",
"keep *_mix_MergedTrackTruth_*",
"keep *_siPixelDigis_*_*",
"keep *_siStripDigis_*_*",
"keep *_siStripDigis_*_*",
"keep *_siPixelClusters_*_*",
"keep *_siStripClusters_*_*",
"keep *_generalTracks_*_*",
"keep *_mixData_MergedTrackTruth_HLT",
"keep *_mix_MergedTrackTruth_*"
)
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.eventinterpretaion_step,process.endjob_step,process.AODSIMoutput_step)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)


