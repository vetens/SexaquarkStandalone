
import os

dir_bash_scripts = "BashScriptsToQSub"
os.mkdir(dir_bash_scripts)

input_files_directory = "/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/GENSIM/CRAB_SimSexaq_7000GENInputFiles/crab_SIMSexaq_16042019_v1/190416_214431/"
output_files_directory = "/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Step1/CRAB_SimSexaq_7000Step1InputFilesQsub/"

idx = 0
for subdir, dirs, files in (os.walk(input_files_directory)):
	for file in files:
		if(file[-4:] == "root"):
			path_and_file_name = "file://" + os.path.join(subdir, file)
			print path_and_file_name
			f = open(dir_bash_scripts+"/"+str(idx)+".sh","w")	
			outputfile = "./Step1"+file

			f.write("#!/bin/bash"+'\n')
			f.write("source $VO_CMS_SW_DIR/cmsset_default.sh"+'\n')
			f.write("export SCRAM_ARCH=slc6_amd64_gcc530"+'\n')
			f.write("cd /user/jdeclerc/CMSSW_8_0_30/src ; eval `scram runtime -sh` ; cd - >/dev/null"+'\n')
			f.write("export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)"+'\n')
			f.write("export X509_USER_PROXY=/tmp/x509up_u20641"+'\n')
			f.write("voms-proxy-init -pwstdin < /user/jdeclerc/CMSSW_8_0_30/src/STEP1_Sexaq/password"+'\n')
			f.write("tmpdir=$TMPDIR"+'\n')
			f.write("mkdir -p $tmpdir"+'\n')
			f.write("cmsRun /user/jdeclerc/CMSSW_8_0_30/src/STEP1_Sexaq_NO_PU/bigsubmission/SUS-RunIISummer16DR80Premix-00068_IIDD_cfg.py inputFiles="+path_and_file_name +" " + "outputFile="+outputfile+'\n')
			f.write("gfal-mkdir srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Step1/`basename $basedir` 2>/dev/null"+'\n')
			f.write("for f in `ls *root` ; do"+'\n')
			f.write('\t'+"gfal-copy "+outputfile+" srm://maite.iihe.ac.be:8443"+output_files_directory+'\n')
			f.write("done"+'\n')
			f.write("rm "+outputfile + '\n')

			f.close()
			
			

			idx +=1
