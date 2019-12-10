from WMCore.Configuration import Configuration

day = "02102019"
version = "v1"
trial = "21"
mass = "1p8GeV"


config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'FlatTreeProducerBDT_trial'+trial+'_'+mass+'_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = '../FlatTreeProducerBDT_cfg.py' 
config.JobType.priority = 103

config.section_('Data') 
config.Data.unitsPerJob = 20
config.Data.totalUnits = 8910
config.Data.publication = True
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/crmc_Sexaq/FlatTree_BDT' 
#config.Data.userInputFiles = open('/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/crab/inputFiles_Skimming_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_27102019_v8_WITH_SEXAQ.txt').readlines() 
config.Data.inputDataset = '/CRAB_SimSexaq_trial21/lowette-crab_Step1_Step2_Skimming_FlatTree_trial21_1p8GeV_23102019_v1-8925145305413877174dac643a893255/USER'
config.Data.inputDBS = 'phys03'
#config.Data.outputPrimaryDataset = "FlatTreeProducerBDT_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_"+day+'_'+version

config.section_('User') 
config.User.voGroup = 'becms'

config.section_('Site') 
#config.Site.whitelist =['T2_BE_IIHE','T2_AT_Vienna','T2_BE_UCL','T2_CH_CERN_HLT','T2_DE_DESY','T2_DE_RWTH','T2_FR_IPHC','T1_DE_KIT','T1_UK_RAL','T2_HU_Budapest','T2_IT_Bari','T2_IT_Legnaro','T2_IT_Pisa']
config.Site.whitelist =['T2_BE_IIHE']
config.Site.storageSite = 'T2_BE_IIHE'
