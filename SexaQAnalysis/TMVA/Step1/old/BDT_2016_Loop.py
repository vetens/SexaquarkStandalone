import ROOT
from ROOT import *
import os
import shutil
import glob
import sys

gROOT.SetBatch(ROOT.kTRUE)

############################################
######HERE DEFINE ALL INPUT#################
############################################

SignFile1 = ROOT.TFile.Open(config_dict["config_SignalFile"]) 

BkgFile  = ROOT.TFile.Open(config_dict["config_BkgFile"])

# Get signal and background trees from file
SignalTree1     = SignFile1.Get("FlatTreeProducerBDT/FlatTree")
#I only want the antiS which match a GEN antiS in lxyz of the interaction vertex and the charge should also be negative
gROOT.cd()
selectedSignalTree1 = SignalTree1.CopyTree(config_dict["config_SelectionSignalAntiS"]) 

BkgTree        = BkgFile.Get("FlatTreeProducerBDT/FlatTree")
#I only want the S with positive charge
gROOT.cd()
#for selecting Data S BKG: ---> standard
selectedBkgTree = BkgTree.CopyTree(config_dict["config_SelectionBkgS"])

trainTestSplit = 0.8



###########################################
####NOW RUN THE REAL STUFF#################
###########################################
localoutputdir = sys.argv[2]
#os.mkdir(localoutputdir)

#outputdir = sys.argv[3]
#os.mkdir(outputdir)



#helper function to move all files from one dir to another
def moveAllFilesinDir(srcDir, dstDir):
    # Check if both the are directories
    if os.path.isdir(srcDir) and os.path.isdir(dstDir) :
        # Iterate over all the files in source directory
        for filePath in glob.glob(srcDir + '\*'):
            # Move each file to destination Directory
            shutil.move(filePath, dstDir);
    else:
        print("srcDir & dstDir should be Directories")


def SingleTraining(kin_variables,iteration):
	# define the dataloader
	dataloader = ROOT.TMVA.DataLoader(localoutputdir+'/dataset_BDT_v'+str(iteration))

	# Add variables to dataloader
	for variable in kin_variables:
		dataloader.AddVariable(variable)

	dataloader.AddSignalTree(selectedSignalTree1, 1)
	dataloader.AddBackgroundTree(selectedBkgTree, 1)

	#adding the cuts and so
	dataloader.PrepareTrainingAndTestTree(ROOT.TCut(config_dict["config_pre_BDT_cuts"]),\
		'TrainTestSplit_Signal={}:'.format(trainTestSplit)+\
		'TrainTestSplit_Background={}:'.format(trainTestSplit)+'SplitMode=Random')

	# Setup TMVA
	ROOT.TMVA.Tools.Instance()
	ROOT.TMVA.PyMethodBase.PyInitialize()

	outputFile = ROOT.TFile.Open(localoutputdir+'/plots_BDT_v'+str(iteration)+'.root', 'RECREATE')
	factory = ROOT.TMVA.Factory('TMVAClassification', outputFile,
		'!V:!Silent:Color:Transformations=I:'+\
		'AnalysisType=Classification')


	factory.BookMethod(dataloader,'BDT', 'BDT',
		'H:!V:VarTransform=None:'+\
		'NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=12:UseBaggedBoost=True')

	factory.TrainAllMethods()

	factory.TestAllMethods()

	factory.EvaluateAllMethods()

	# Enable Javascript for ROOT so that we can draw the canvas
	#ROOT.enableJSVis()
	#%jsroot on
	# Print ROC
	#canvas = factory.GetROCCurve(dataloader)
	#canvas.Draw()
	#canvas.SaveAs("BDT_2016_v"+str(iteration)+".root")
	outputFile.Close()
	#moveAllFilesinDir(localoutputdir, outputdir)

#list of variables to cut on
AllKinVariables = [
"_S_vz_interaction_vertex",
#"_S_lxy_interaction_vertex",

#"_S_daughters_deltaphi",
#"_S_daughters_deltaeta",
"_S_daughters_openingsangle",
#"_S_daughters_DeltaR",
"_S_Ks_openingsangle",
#"_S_Lambda_openingsangle",

#"_S_eta",
#"_Ks_eta",
"_Lambda_eta",

"_S_dxy_over_lxy",
#"_Ks_dxy_over_lxy",
"_Lambda_dxy_over_lxy",

#"_S_dz_min",
"_Ks_dz_min",
#"_Lambda_dz_min",

"_Ks_pt",

#"_Lambda_lxy_decay_vertex",
"_S_chi2_ndof",
#"_S_pz"

]

#run SingleTraining a few times by always dropping one of the KinVariables, the variables which will be dropped is the one corresponding to the index 'i_var'
#for i_var in range(0,len(AllKinVariables)) :
#	kinVariablesSubset = []
#	for i_var_subset in range(0,len(AllKinVariables)):
#		if(i_var_subset is not i_var):#the only variable not to add is the one corresponding to the index 'i_var'
#			kinVariablesSubset.append(AllKinVariables[i_var_subset])
#	print "-----> The full Kin Variable list is: "
#	print AllKinVariables
#	print "-----> But will for now only train on: "
#	print kinVariablesSubset
#	SingleTraining(kinVariablesSubset,i_var) 
	

kinVariablesSubset = []
i_var = int(sys.argv[1])
for i_var_subset in range(0,len(AllKinVariables)):
	if(i_var_subset is not i_var):#the only variable not to add is the one corresponding to the index 'i_var'
		kinVariablesSubset.append(AllKinVariables[i_var_subset])
print "-----> The full Kin Variable list is: "
print AllKinVariables
print "-----> But will for now only train on: "
print kinVariablesSubset
SingleTraining(kinVariablesSubset,i_var) 
