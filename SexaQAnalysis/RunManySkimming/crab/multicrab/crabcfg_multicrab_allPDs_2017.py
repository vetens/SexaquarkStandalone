from CRABClient.UserUtilities import config
trial = "K"
config = config()

pyCfgParams = ['isData=True']

#config.General.requestName = 'SexaQ_SingleMuon'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../treeproducer_data_cfg.py'

#config.Data.inputDataset = '/SingleMuon/Run2016G-23Sep2016-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 500000
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
config.Data.runRange = ''
config.Data.outLFNDirBase = '/store/user/jdeclerc/data_Sexaq/trial'+trial
config.Data.publication = False
#config.Data.outputDatasetTag = 'allPDs_2016_multicrab'
#config.Data.splitting = 'Automatic'

config.Site.storageSite = 'T2_BE_IIHE'

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
    config.General.workArea = 'crab_projects17_trial'+trial

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
#submit_single_run("BTagCSV","Run2017B-17Nov2017-v1")
#submit_single_run("BTagCSV","Run2017C-17Nov2017-v1")
#submit_single_run("BTagCSV","Run2017D-17Nov2017-v1")
#submit_single_run("BTagCSV","Run2017E-17Nov2017-v1")
#submit_single_run("BTagCSV","Run2017F-17Nov2017-v1")
#submit_single_run("BTagMu","Run2017B-17Nov2017-v1")
#submit_single_run("BTagMu","Run2017C-17Nov2017-v1")
#submit_single_run("BTagMu","Run2017D-17Nov2017-v1")
#submit_single_run("BTagMu","Run2017E-17Nov2017-v1")
#submit_single_run("BTagMu","Run2017F-17Nov2017-v1")
#submit_single_run("Charmonium","Run2017B-17Nov2017-v1")
#submit_single_run("Charmonium","Run2017C-17Nov2017-v1")
#submit_single_run("Charmonium","Run2017D-17Nov2017-v1")
#submit_single_run("Charmonium","Run2017E-17Nov2017-v1")
#submit_single_run("Charmonium","Run2017F-17Nov2017-v1")
#submit_single_run("DisplacedJet","Run2017C-17Nov2017-v1")
#submit_single_run("DisplacedJet","Run2017D-17Nov2017-v1")
#submit_single_run("DisplacedJet","Run2017E-17Nov2017-v1")
#submit_single_run("DisplacedJet","Run2017F-17Nov2017-v1")
#submit_single_run("DoubleEG","Run2017B-17Nov2017-v1")
#submit_single_run("DoubleEG","Run2017C-17Nov2017-v1")
#submit_single_run("DoubleEG","Run2017D-17Nov2017-v1")
#submit_single_run("DoubleEG","Run2017E-17Nov2017-v1")
#submit_single_run("DoubleEG","Run2017F-09May2018-v1")
#submit_single_run("DoubleMuon","Run2017B-17Nov2017-v1")
#submit_single_run("DoubleMuon","Run2017C-17Nov2017-v1")
#submit_single_run("DoubleMuon","Run2017D-17Nov2017-v1")
#submit_single_run("DoubleMuon","Run2017E-17Nov2017-v1")
#submit_single_run("DoubleMuon","Run2017F-17Nov2017-v1")
#submit_single_run("DoubleMuon","Run2017G-17Nov2017-v1")
#submit_single_run("DoubleMuon","Run2017H-17Nov2017-v1")
#submit_single_run("DoubleMuonLowMass","Run2017B-17Nov2017-v1")
#submit_single_run("DoubleMuonLowMass","Run2017C-17Nov2017-v1")
#submit_single_run("DoubleMuonLowMass","Run2017D-17Nov2017-v1")
#submit_single_run("DoubleMuonLowMass","Run2017E-17Nov2017-v1")
#submit_single_run("DoubleMuonLowMass","Run2017F-17Nov2017-v1")
#submit_single_run("HTMHT","Run2017B-17Nov2017-v1")
#submit_single_run("HTMHT","Run2017C-17Nov2017-v1")
#submit_single_run("HTMHT","Run2017D-17Nov2017-v1")
#submit_single_run("HTMHT","Run2017E-17Nov2017-v1")
#submit_single_run("HTMHT","Run2017F-17Nov2017-v1")
#submit_single_run("JetHT","Run2017A-12Sep2017-v1")
#submit_single_run("JetHT","Run2017B-17Nov2017-v1")
#submit_single_run("JetHT","Run2017C-17Nov2017-v1")
#submit_single_run("JetHT","Run2017D-17Nov2017-v1")
#submit_single_run("JetHT","Run2017E-17Nov2017-v1")
#submit_single_run("JetHT","Run2017F-17Nov2017-v1")
#submit_single_run("MET","Run2017B-17Nov2017-v1")
#submit_single_run("MET","Run2017C-17Nov2017-v1")
#submit_single_run("MET","Run2017D-17Nov2017-v1")
#submit_single_run("MET","Run2017E-17Nov2017-v1")
#submit_single_run("MET","Run2017F-17Nov2017-v1")
#submit_single_run("MuOnia","Run2017B-17Nov2017-v1")
#submit_single_run("MuOnia","Run2017C-17Nov2017-v1")
#submit_single_run("MuOnia","Run2017D-17Nov2017-v1")
#submit_single_run("MuOnia","Run2017E-17Nov2017-v1")
#submit_single_run("MuOnia","Run2017F-17Nov2017-v1")
#submit_single_run("MuonEG","Run2017B-17Nov2017-v1")
#submit_single_run("MuonEG","Run2017C-17Nov2017-v1")
#submit_single_run("MuonEG","Run2017D-17Nov2017-v1")
#submit_single_run("MuonEG","Run2017E-17Nov2017-v1")
#submit_single_run("MuonEG","Run2017F-17Nov2017-v1")
#submit_single_run("SingleElectron","Run2017B-17Nov2017-v1")
#submit_single_run("SingleElectron","Run2017C-17Nov2017-v1")
#submit_single_run("SingleElectron","Run2017D-17Nov2017-v1")
#submit_single_run("SingleElectron","Run2017E-17Nov2017-v1")
#submit_single_run("SingleElectron","Run2017F-17Nov2017-v1")
submit_single_run("SingleMuon","Run2017B-17Nov2017-v1")
submit_single_run("SingleMuon","Run2017C-17Nov2017-v1")
submit_single_run("SingleMuon","Run2017D-17Nov2017-v1")
submit_single_run("SingleMuon","Run2017E-17Nov2017-v1")
submit_single_run("SingleMuon","Run2017F-17Nov2017-v1")
submit_single_run("SingleMuon","Run2017G-17Nov2017-v1")
submit_single_run("SingleMuon","Run2017H-17Nov2017-v2")
#submit_single_run("SinglePhoton","Run2017B-17Nov2017-v1")
#submit_single_run("SinglePhoton","Run2017C-17Nov2017-v1")
#submit_single_run("SinglePhoton","Run2017D-17Nov2017-v1")
#submit_single_run("SinglePhoton","Run2017E-17Nov2017-v1")
#submit_single_run("SinglePhoton","Run2017F-17Nov2017-v1")
#submit_single_run("Tau","Run2017B-17Nov2017-v1")
#submit_single_run("Tau","Run2017C-17Nov2017-v1")
#submit_single_run("Tau","Run2017D-17Nov2017-v1")
#submit_single_run("Tau","Run2017E-17Nov2017-v1")
#submit_single_run("Tau","Run2017F-17Nov2017-v1")


