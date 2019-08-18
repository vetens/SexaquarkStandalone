from WMCore.Configuration import Configuration

day = "23072019"
version = "v1"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'Step1Sexaq_trial13_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'step1_noSexaq_bad_efficiency_tryToFix_8.py' 
config.JobType.numCores = 1

config.section_('Data') 
config.Data.unitsPerJob = 1 
config.Data.totalUnits = 10000
config.Data.publication = False 
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/crmc_Sexaq/Step1' 
config.Data.userInputFiles = open('/user/jdeclerc/CMSSW_8_0_30/src/STEP1_Sexaq/crab/inputFiles_Trial13.txt').readlines() 
config.Data.outputPrimaryDataset = "CRAB_Step1Sexaq_trial13"

config.section_('User') 
config.User.voGroup = 'becms'

config.section_('Site') 
config.Site.whitelist = ['T2_BE_IIHE'] 
config.Site.storageSite = 'T2_BE_IIHE'
