from WMCore.Configuration import Configuration

day = "19082019"
version = "v1"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'FlatTreeProducerV0s_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'FlatTreeProducerV0s_cfg.py' 
config.JobType.priority = 100

config.section_('Data') 
config.Data.unitsPerJob = 50 
config.Data.totalUnits = 1121
config.Data.publication = False
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/MC/FlatTreeV0s' 
config.Data.userInputFiles = open('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/inputFiles_ttZJets_13TeV_madgraphMLM.txt').readlines()
#config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
#config.Data.inputDBS = 'global'
config.Data.outputPrimaryDataset = "crab_FlatTreeV0_ttZJets_13TeV_madgraphMLM"
#config.Data.ignoreLocality = True


config.section_('User') 
config.User.voGroup = 'becms'

config.section_('Site') 
#config.Site.whitelist =['T2_BE_IIHE'] 
config.Site.storageSite = 'T2_BE_IIHE'
