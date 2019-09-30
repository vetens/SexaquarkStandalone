import os 

l_inputFiles = ["inputFiles_FlatTreeProducerBDT_SingleMuon_Run2016G.txt"]
l_outputFiles = ["FlatTreeBDT_SingleMuon_Run2016G.root"]

inputListDir = "/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/inputFilesLists"

i = 0
for inputList in os.listdir(inputListDir):
	file = open("shell/"+str(i)+".sh","w") 
#	file.write("#!/bin/bash"+"\n")
#	file.write("cd /afs/cern.ch/user/j/jdeclerc/CMSSW_9_4_7/src/"+"\n")
#	file.write("cmsenv"+"\n")

	file.write("source $VO_CMS_SW_DIR/cmsset_default.sh"+"\n")
	file.write("export SCRAM_ARCH=slc6_amd64_gcc481"+"\n")
	file.write("cd /user/jdeclerc/CMSSW_8_0_30_bis/src ; eval `scram runtime -sh` ; cd - >/dev/null" + "\n")

#	file.write("export X509_USER_PROXY=/tmp/x509up_u20641"+"\n")

        file.write("export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)"+"\n")
        file.write("export X509_USER_PROXY=/tmp/x509up_u20641"+"\n")
        file.write("voms-proxy-init -pwstdin < /storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30/src/STEP1_STEP2_Skimming_FlatTree/qsub/password"+"\n")
	


	file.write("cd /user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT" +"\n")

	#RunG
	file.write("cmsRun /user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/FlatTreeProducerBDT_cfg.py print inputFiles_load="+inputListDir+"/"+inputList + " outputFile=file:/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/Results/FlatTreeBDT_"+inputList[31:-4]+".root")

	#RunH
	#file.write("cmsRun /storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTreeProducerV0s_cfg.py print inputFiles=file:/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/DataMCComparison/DoubleMuon/crab_SkimmingSexaq_DoubleMuonData_18082019_v6_runH/190818_063505/000"+str(int(i/1000))+"/events_skimmed_"+str(i)+".root outputFile=file:/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTreeV0_DoubleMuon_v6_runH_"+str(i)+".root")

	#DY MC
	#file.write("cmsRun /storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTreeProducerV0s_cfg.py print inputFiles=file:/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/DataMCComparison/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_SkimmingSexaq_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_19082019_v1/190819_150350/000"+str(int(i/1000))+"/events_skimmed_"+str(i)+".root outputFile=file:/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_"+str(i)+".root")
	file.close()
	i+=1

