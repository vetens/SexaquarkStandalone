from WMCore.Configuration import Configuration

day = "20072019"
version = "v1"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'FlatTreeProducerV0s_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'FlatTreeProducerV0s_cfg.py' 

config.section_('Data') 
config.Data.unitsPerJob = 10 
config.Data.totalUnits = 2513
config.Data.publication = False
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/MC/FlatTreeV0s' 
#config.Data.userInputFiles = open('/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30/src/STEP2_Sexaq/crabSlimmedFEVT/inputFiles.txt').readlines()
config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
config.Data.inputDBS = 'global'
#config.Data.outputPrimaryDataset = "crab_Step2SexaqWithPU2016NeutrinoGun_tryToFix_8_trial11"
#config.Data.ignoreLocality = True


config.section_('User') 
#config.User.voGroup = 'becms'

config.section_('Site') 
#config.Site.whitelist =['T2_AT_Vienna','T2_BE_UCL','T2_CH_CERN_HLT','T2_DE_DESY','T2_DE_RWTH','T2_FR_IPHC'] 
config.Site.storageSite = 'T2_BE_IIHE'
