from WMCore.Configuration import Configuration

day = "23082019"
version = "v1"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'FlatTreeProducer_trial10_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'FlatTreeProducerV0s_cfg.py' 
config.JobType.priority = 100

config.section_('Data') 
config.Data.unitsPerJob = 25 
config.Data.totalUnits = 239
config.Data.publication = False
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/MC/FlatTreeV0s' 
config.Data.userInputFiles = open('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/inputFiles_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt').readlines()
config.Data.outputPrimaryDataset = "crab_FlatTree_trial10_"


config.section_('User') 
config.User.voGroup = 'becms'

config.section_('Site') 
#config.Site.whitelist =['T2_BE_IIHE'] 
config.Site.storageSite = 'T2_BE_IIHE'
