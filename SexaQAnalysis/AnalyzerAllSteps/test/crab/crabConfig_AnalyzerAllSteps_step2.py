from WMCore.Configuration import Configuration

day = "02052019"
version = "v2"

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'AnalyzerAllStepsSexaqStep2_'+day+'_'+version 

config.section_('JobType') 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'analyzerAllSteps_cfg.py' 

config.section_('Data') 
config.Data.unitsPerJob = 1 
config.Data.totalUnits = 2448
config.Data.publication = False 
config.Data.splitting = 'FileBased' 
config.Data.outLFNDirBase = '/store/user/jdeclerc/crmc_Sexaq/Analyzed_Step2' 
config.Data.userInputFiles = open('/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/AnalyzerAllSteps/test/crab/inputFilesAllStep2_changedXYCosThetaCut.txt').readlines() 
config.Data.outputPrimaryDataset = "CRAB_AnalyzerAllStepsStep2_Background_DisabledCosThetaXYCutV0Fitter"+day+'_'+version

config.section_('User') 
config.User.voGroup = 'becms'

config.section_('Site') 
config.Site.whitelist = ['T2_BE_IIHE'] 
config.Site.storageSite = 'T2_BE_IIHE'
