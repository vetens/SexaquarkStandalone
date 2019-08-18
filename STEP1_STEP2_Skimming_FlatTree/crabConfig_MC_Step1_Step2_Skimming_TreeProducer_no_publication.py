from WMCore.Configuration import Configuration

day = "05082019"
version = "v1"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'Step1_Step2_Skimming_FlatTree_Sexaq_trial14_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'master_STEP1_STEP2_Skimming_FlatTree_SLIMMED_cfg.py' 
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2500

config.section_('Data') 
config.Data.unitsPerJob = 1 
config.Data.totalUnits = 1000
config.Data.publication = False 
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/crmc_Sexaq/Skimmed' 
config.Data.userInputFiles = open('/user/jdeclerc/CMSSW_8_0_30/src/STEP1_STEP2_Skimming_FlatTree/inputFiles_trial14.txt').readlines() 
config.Data.outputPrimaryDataset = "CRAB_Step1Sexaq_trial14"

config.section_('User') 
config.User.voGroup = 'becms'

config.section_('Site') 
config.Site.whitelist = ['T2_BE_IIHE'] 
config.Site.storageSite = 'T2_BE_IIHE'
