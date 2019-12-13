source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc481
cd /user/jdeclerc/CMSSW_8_0_30_bis/src ; eval `scram runtime -sh` ; cd - >/dev/null
export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
export X509_USER_PROXY=/tmp/x509up_u20641
voms-proxy-init -pwstdin < /storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30/src/STEP1_STEP2_Skimming_FlatTree/qsub/password
cd /user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT
cmsRun /user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/FlatTreeProducerBDT_cfg.py print inputFiles_load=/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/inputFilesLists/inputFiles_FlatTreeProducerBDT_MuOnia_Run2016G-07Aug17-v2_trialR.txt outputFile=file:/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/Results/FlatTreeBDT_MuOnia_Run2016G-07Aug17-v2_trialR.root