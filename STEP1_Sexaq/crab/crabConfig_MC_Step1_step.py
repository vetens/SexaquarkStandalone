from WMCore.Configuration import Configuration

day = "26082019"
version = "v1"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'Step1Sexaq_trial15_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'step1_noSexaq_bad_efficiency_tryToFix_8.py' 
config.JobType.numCores = 1

config.section_('Data') 
config.Data.unitsPerJob = 1 
config.Data.totalUnits = 1000
config.Data.publication = False 
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/crmc_Sexaq/Step1' 
config.Data.userInputFiles = open('/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/STEP1_Sexaq/crab/inputFiles_SIMSexaq_trial15_11082019_v1.txt').readlines() 
config.Data.outputPrimaryDataset = "CRAB_Step1Sexaq_trial15"

config.section_('User') 
config.User.voGroup = 'becms'

config.section_('Site') 
config.Site.whitelist = ['T2_BE_IIHE'] 
config.Site.storageSite = 'T2_BE_IIHE'
