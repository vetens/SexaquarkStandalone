import os 


inputListDir = "/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/inputFilesLists"

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

i = 0
for inputList in os.listdir(inputListDir):
	name = inputList[31:-11]

	#before making the file count the number  of lines in the input txt file
	nLines = file_len(inputListDir+'/'+inputList)
	print inputList,': ', nLines

	file = open("shell/crabConfig_"+name+".py","w") 
#	file.write("#!/bin/bash"+"\n")
	file.write('from WMCore.Configuration import Configuration \n')
	file.write('day=\"03102019\"\n')
	file.write('version=\"v1\"\n')
	file.write('trial = \"'+name+'\"\n')
	file.write('mass = \"\"\n')

	file.write('config = Configuration()\n')
	file.write('config.section_(\'General\')\n')
	file.write('config.General.transferOutputs = True\n')
	file.write('config.General.transferLogs = True\n')
	file.write('config.General.requestName = \'FlatTreeProducerBDT_trial\'+trial+\'_\'+mass+\'_\'+day+\'_\'+version\n')

	file.write('config.section_(\'JobType\')\n')	
	file.write('config.JobType.pluginName = \'Analysis\'\n')	
	file.write('config.JobType.psetName = \'../../FlatTreeProducerBDT_cfg.py\'\n')	
	file.write('config.JobType.priority = 103\n')	

	file.write('config.section_(\'Data\')\n')
	file.write('config.Data.unitsPerJob = 50\n')
	file.write('config.Data.totalUnits = '+str(nLines)+'\n')
	file.write('config.Data.publication = False\n')
	file.write('config.Data.splitting = \'FileBased\'\n')
	file.write('config.Data.outLFNDirBase = \'/store/user/jdeclerc/data_Sexaq/trialR/ALL/ALL_v7\'\n')
	file.write('config.Data.userInputFiles = open(\''+inputListDir+'/'+inputList+'\').readlines()\n')
	file.write('config.Data.outputPrimaryDataset = \"FlatTreeProducerBDT_\"+trial+\"_\"+day+\'_\'+version\n')

	file.write('config.section_(\'User\')\n')		
	file.write('config.User.voGroup = \'becms\'\n')		

	file.write('config.section_(\'Site\')\n')
	file.write('config.Site.whitelist =[\'T2_BE_IIHE\']\n')
	file.write('config.Site.storageSite = \'T2_BE_IIHE\'\n')

	file.close()
	i+=1

