from WMCore.Configuration import Configuration

day = "10082019"
version = "v1"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'Step1_Step2_Skimming_FlatTree_trial15_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'master_STEP1_STEP2_Skimming_FlatTree_SLIMMED_cfg.py' 
config.JobType.maxMemoryMB = 2500

config.section_('Data') 
config.Data.unitsPerJob = 1 
config.Data.totalUnits = 10000
config.Data.publication = True 
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/crmc_Sexaq/Skimmed' 
#config.Data.userInputFiles = open('/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30/src/STEP2_Sexaq/crabSlimmedFEVT/inputFiles.txt').readlines()
config.Data.inputDataset = '/CRAB_SimSexaq_trial15/jdeclerc-crab_SIMSexaq_trial15_08082019_v1-33065c3a8a979241c5e811e20eaee23d/USER'
config.Data.inputDBS = 'phys03'
#config.Data.outputPrimaryDataset = "crab_Step2SexaqWithPU2016NeutrinoGun_tryToFix_8_trial11"
config.Data.ignoreLocality = True


config.section_('User') 
#config.User.voGroup = 'becms'

config.section_('Site') 
config.Site.whitelist =['T2_BE_IIHE','T2_AT_Vienna','T2_BE_UCL','T2_CH_CERN_HLT','T2_DE_DESY','T2_DE_RWTH','T2_FR_IPHC','T1_DE_KIT','T1_UK_RAL','T2_HU_Budapest','T2_IT_Bari','T2_IT_Legnaro','T2_IT_Pisa'] 
config.Site.storageSite = 'T2_BE_IIHE'
