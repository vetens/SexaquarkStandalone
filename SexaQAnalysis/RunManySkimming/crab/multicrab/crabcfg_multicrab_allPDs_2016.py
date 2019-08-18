from CRABClient.UserUtilities import config
trial = "R"
config = config()

pyCfgParams = ['isData=True']

#config.General.requestName = 'SexaQ_SingleMuon'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../treeproducer_data_cfg.py'


config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 200000
#config.Data.splitting = 'Automatic'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.runRange = ''
config.Data.outLFNDirBase = '/store/user/jdeclerc/data_Sexaq/trial'+trial
config.Data.publication = False

config.Site.storageSite = 'T2_BE_IIHE'
#if only using the below whitelist then a lot of jobs which cannot find the datasets
#config.Site.whitelist = ['T2_BE_IIHE','T2_HU_Budapest','T2_ES_CIEMAT','T2_EE_Estonia','T1_IT_CNAF','T2_US_Caltech','T2_US_MIT','T2_US_Nebraska','T2_US_Florida','T2_US_Wisconsin','T2_US_Vanderbilt','T2_CH_CERN','T2_PK_NCP','T1_US_FNAL','T2_UK_London_IC','T2_UK_SGrid_RALPP','T2_CN_Beijing','T2_DE_RWTH','T2_US_Purdue','T2_AT_Vienna']
#very thight blacklist:
#config.Site.blacklist = ['T1_DE_KIT','T2_GR_Ioannina','T2_IT_Legnaro','T1_UK_RAL','T1_RU_JINR','T1_FR_CCIN2P3','T2_IT_Bari','T2_IN_TIFR','T2_US_UCSD','T2_IT_Rome','T2_BR_UERJ','T2_FR_GRIF_LLR','T1_FR_CCIN2P3']
#loosened blacklist
#config.Site.blacklist = ['T2_GR_Ioannina','T2_IT_Legnaro','T2_IT_Bari','T2_IN_TIFR','T2_IT_Rome','T2_BR_UERJ','T2_FR_GRIF_LLR']

#blacklist IIHE, because you want to keep your quota free at this site to run the MC production
#config.Site.blacklist = ['T2_BE_IIHE']

def submit_single_run(arg1, arg2):

	#config.General.requestName = 'BTagCSV_Run2016B-07Aug17_ver2-v1'
	word1 =  arg1 + "_" + arg2 + "_trial"+trial
	print word1
	config.General.requestName = word1
	config.Data.outputDatasetTag = word1
	#config.Data.inputDataset = '/BTagCSV/Run2016B-07Aug17_ver2-v1/AOD'
	word2 =  "/" + arg1 + "/" + arg2 + "/AOD" 
	print word2
	config.Data.inputDataset =  word2
	submit(config)


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects16_trial'+trial

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
#submit_single_run("ZeroBias","Run2016B-07Aug17-v1")
#submit_single_run("ZeroBias","Run2016C-07Aug17-v1")
#submit_single_run("ZeroBias","Run2016E-07Aug17-v1")
#submit_single_run("ZeroBias","Run2016F-07Aug17-v1")
#submit_single_run("ZeroBias","Run2016H-07Aug17-v1")

###submit_single_run("BTagCSV","Run2016B-07Aug17_ver2-v1") #--> only availalbe on TAPE
#submit_single_run("BTagCSV","Run2016C-07Aug17-v1")
####submit_single_run("BTagCSV","Run2016D-07Aug17-v1") #--> only availalbe on TAPE
#submit_single_run("BTagCSV","Run2016E-07Aug17-v1")
#submit_single_run("BTagCSV","Run2016F-07Aug17-v1")
####submit_single_run("BTagCSV","Run2016G-07Aug17-v1") #--> only availalbe on TAPE
submit_single_run("BTagCSV","Run2016H-07Aug17-v1")
####submit_single_run("BTagMu","Run2016B-07Aug17_ver2-v1") #--> only availalbe on TAPE
#submit_single_run("BTagMu","Run2016C-07Aug17-v1")
#submit_single_run("BTagMu","Run2016D-07Aug17-v1")
#submit_single_run("BTagMu","Run2016E-07Aug17-v1")
#submit_single_run("BTagMu","Run2016F-07Aug17-v1")
#submit_single_run("BTagMu","Run2016G-07Aug17-v1")
#submit_single_run("BTagMu","Run2016H-07Aug17-v1")
#submit_single_run("Charmonium","Run2016B-07Aug17_ver2-v1")
#submit_single_run("Charmonium","Run2016C-07Aug17-v1")
#submit_single_run("Charmonium","Run2016D-07Aug17-v1")
#submit_single_run("Charmonium","Run2016E-07Aug17-v1")
#submit_single_run("Charmonium","Run2016F-07Aug17-v1")
#submit_single_run("Charmonium","Run2016G-07Aug17-v1")
#submit_single_run("Charmonium","Run2016H-07Aug17-v1")
#submit_single_run("DisplacedJet","Run2016B-07Aug17_ver2-v1")
####submit_single_run("DisplacedJet","Run2016C-07Aug17-v1") #--> only availalbe on TAPE
#submit_single_run("DisplacedJet","Run2016D-07Aug17-v1")
#submit_single_run("DisplacedJet","Run2016E-07Aug17-v1")
#submit_single_run("DisplacedJet","Run2016F-07Aug17-v1")
####submit_single_run("DisplacedJet","Run2016G-07Aug17-v1") #--> only availalbe on TAPE
#submit_single_run("DisplacedJet","Run2016H-07Aug17-v1")
#submit_single_run("DoubleEG","Run2016B-07Aug17_ver2-v2")
####submit_single_run("DoubleEG","Run2016C-07Aug17-v1") #--> only availalbe on TAPE
#submit_single_run("DoubleEG","Run2016D-07Aug17-v1")
#submit_single_run("DoubleEG","Run2016E-07Aug17-v1")
#submit_single_run("DoubleEG","Run2016F-07Aug17-v1")
####submit_single_run("DoubleEG","Run2016G-07Aug17-v1") #--> only availalbe on TAPE
####submit_single_run("DoubleEG","Run2016H-07Aug17-v1") #--> only availalbe on TAPE
#submit_single_run("DoubleMuon","Run2016B-07Aug17_ver2-v1")
#submit_single_run("DoubleMuon","Run2016C-07Aug17-v1")
#submit_single_run("DoubleMuon","Run2016D-07Aug17-v1")
#submit_single_run("DoubleMuon","Run2016E-07Aug17-v1")
#submit_single_run("DoubleMuon","Run2016F-07Aug17-v1")
#submit_single_run("DoubleMuon","Run2016G-07Aug17-v1")
#submit_single_run("DoubleMuon","Run2016H-07Aug17-v1")
#submit_single_run("DoubleMuonLowMass","Run2016B-07Aug17_ver2-v1")
#submit_single_run("DoubleMuonLowMass","Run2016C-07Aug17-v1")
#submit_single_run("DoubleMuonLowMass","Run2016D-07Aug17-v2")
#submit_single_run("DoubleMuonLowMass","Run2016E-07Aug17-v1")
#submit_single_run("DoubleMuonLowMass","Run2016F-07Aug17-v1")
####submit_single_run("DoubleMuonLowMass","Run2016G-07Aug17-v1") #--> only availalbe on TAPE
#submit_single_run("DoubleMuonLowMass","Run2016H-07Aug17-v1")
#submit_single_run("HTMHT","Run2016B-07Aug17_ver2-v1")
#submit_single_run("HTMHT","Run2016C-07Aug17-v1")
#submit_single_run("HTMHT","Run2016D-07Aug17-v1")
#submit_single_run("HTMHT","Run2016E-07Aug17-v1")
#submit_single_run("HTMHT","Run2016F-07Aug17-v1")
#submit_single_run("HTMHT","Run2016G-23Sep2016-v2")
#submit_single_run("HTMHT","Run2016H-07Aug17-v1")
#submit_single_run("JetHT","Run2016B-07Aug17_ver2-v1")
#submit_single_run("JetHT","Run2016C-07Aug17-v1")
####submit_single_run("JetHT","Run2016D-07Aug17-v1") #--> only availalbe on TAPE
#submit_single_run("JetHT","Run2016E-07Aug17-v1")
#submit_single_run("JetHT","Run2016F-07Aug17-v1")
#submit_single_run("JetHT","Run2016G-07Aug17-v1")
#submit_single_run("JetHT","Run2016H-07Aug17-v1")
#submit_single_run("MET","Run2016B-07Aug17_ver2-v1")
#submit_single_run("MET","Run2016C-07Aug17-v1")
#submit_single_run("MET","Run2016D-07Aug17-v1")
#submit_single_run("MET","Run2016E-07Aug17-v1")
#submit_single_run("MET","Run2016F-07Aug17-v1")
#submit_single_run("MET","Run2016G-07Aug17-v1")
#submit_single_run("MET","Run2016H-07Aug17-v1")
#submit_single_run("MuOnia","Run2016B-07Aug17_ver2-v1")
#submit_single_run("MuOnia","Run2016C-07Aug17-v1")
#submit_single_run("MuOnia","Run2016D-07Aug17-v1")
#submit_single_run("MuOnia","Run2016E-07Aug17-v1")
#submit_single_run("MuOnia","Run2016F-07Aug17-v1")
#submit_single_run("MuOnia","Run2016G-07Aug17-v2")
#submit_single_run("MuOnia","Run2016H-07Aug17-v1")
#submit_single_run("MuonEG","Run2016B-07Aug17_ver2-v1")
#submit_single_run("MuonEG","Run2016C-07Aug17-v1")
#submit_single_run("MuonEG","Run2016D-07Aug17-v1")
#submit_single_run("MuonEG","Run2016E-07Aug17-v1")
#submit_single_run("MuonEG","Run2016F-07Aug17-v1")
#submit_single_run("MuonEG","Run2016G-07Aug17-v1")
#submit_single_run("MuonEG","Run2016H-07Aug17-v1")
#submit_single_run("SingleElectron","Run2016B-07Aug17_ver2-v2")
####submit_single_run("SingleElectron","Run2016C-07Aug17-v1") #--> only availalbe on TAPE
#submit_single_run("SingleElectron","Run2016D-07Aug17-v1")
#submit_single_run("SingleElectron","Run2016E-07Aug17-v1")
#submit_single_run("SingleElectron","Run2016F-07Aug17-v1")
#submit_single_run("SingleElectron","Run2016G-07Aug17-v1")
###submit_single_run("SingleElectron","Run2016H-07Aug17-v1") #--> only availalbe on TAPE

#submit_single_run("SingleMuon","Run2016B-07Aug17_ver2-v1")
#submit_single_run("SingleMuon","Run2016C-07Aug17-v1")
#submit_single_run("SingleMuon","Run2016D-07Aug17-v1")
#submit_single_run("SingleMuon","Run2016E-07Aug17-v1")
#submit_single_run("SingleMuon","Run2016F-07Aug17-v1")
#submit_single_run("SingleMuon","Run2016G-07Aug17-v1")
#submit_single_run("SingleMuon","Run2016H-07Aug17-v1")

#submit_single_run("SinglePhoton","Run2016B-07Aug17_ver2-v1")
#submit_single_run("SinglePhoton","Run2016C-07Aug17-v1")
#submit_single_run("SinglePhoton","Run2016D-07Aug17-v1")
#submit_single_run("SinglePhoton","Run2016E-07Aug17-v1")
#submit_single_run("SinglePhoton","Run2016F-07Aug17-v1")
#submit_single_run("SinglePhoton","Run2016G-07Aug17-v1")
#submit_single_run("SinglePhoton","Run2016H-07Aug17-v1")
#submit_single_run("Tau","Run2016B-07Aug17_ver2-v1")
#submit_single_run("Tau","Run2016C-07Aug17-v1")
#submit_single_run("Tau","Run2016D-07Aug17-v1")
#submit_single_run("Tau","Run2016E-07Aug17-v1")
#submit_single_run("Tau","Run2016F-07Aug17-v1")
#submit_single_run("Tau","Run2016G-07Aug17-v1")
#submit_single_run("Tau","Run2016H-07Aug17-v1")

