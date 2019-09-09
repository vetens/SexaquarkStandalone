from WMCore.Configuration import Configuration

day = "27082019"
version = "v1"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'SkimmingSexaq_trial15_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'treeproducer_data_cfg.py' 

config.section_('Data') 
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1000
config.Data.publication = False 
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/crmc_Sexaq/Skimmed' 
config.Data.userInputFiles = open('/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/RunManySkimming/crab/inputFiles_Step2_trial15.txt').readlines() 
config.Data.outputPrimaryDataset = "CRAB_SimSexaq_Skimming_trial15_"+day+"_"+version

config.section_('User') 
config.User.voGroup = 'becms'

config.section_('Site') 
config.Site.whitelist = ['T2_BE_IIHE'] 
config.Site.storageSite = 'T2_BE_IIHE'
