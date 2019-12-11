import os
import array
import numpy

import ROOT
from ROOT import *

import sys
sys.path.insert(1, './../..')
import configBDT as config

config_dict = config.config_dict

config = "bkgReference" #"partialUnblinding" (use antiS data, but do not look at stuff with BDT > 0.1), "bkgReference" (use the S as bkg reference) or unblind (use the antiS data fully) or unblindMC (use the antiS data fully, because it is MC so it is fine), 10%Unblind (unblind 10% of the data)
performOverlapCheck = False

#pointer to the results of the training:
dataset = "dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing"
#for the data:
input_directory_data = "/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/ALL/ALL_v7"
#if you want to apply the BDT on the signal MC
#input_directory_data = "/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Final/1p8GeV/FlatTreeBDTMCSignal"
#if you want to apply the BDT on MC which for sure does not contain signal (the DYJets sample also used for systematics study is also used here):
#input_directory_data = "/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Final/1p8GeV/FlatTreeBDTMCBackground"
#for data, SingleMuon Run2016H only
#input_directory_data = "/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/ALL/SingleMuonRunHOnly_v7"
#for data, SingleMuon Run2016H but with S and antiS reconstructed X events, so this one you can unblind
#input_directory_data = "/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialtrialXEventSReco/SingleMuon/SingleMuon_Run2016H-07Aug17-v1_trialtrialXEventSReco/191123_211600/FlatTreeBDT/XEventSReco"

overlapList = []
def checkInOverlapList(eta,S_lxy_interaction_vertex,BDT,fileIn):
	entryFound = False
	for e in overlapList:#loop over the overlap list
		if (e[0] == eta and e[1] == S_lxy_interaction_vertex): #if the eta has already been found count plus 1 and append the file list to this entry
			e[3]+=1
			print "found a ", e[3], " th/nd/rd duplicate!!, with eta: ", eta 
			e.append(fileIn)
			entryFound = True
	
	if(not entryFound): 
		overlapList.append([eta,S_lxy_interaction_vertex,BDT,1,fileIn])
	return entryFound


class TreeCloner(object):
		
	#make a separate directory to store the result
	dirname = "BDTApplied_"+config+"_"+dataset+"_OverlapCheck"+str(performOverlapCheck)
	if not os.path.exists(dirname):
    		os.makedirs(dirname)
	cwd = os.getcwd()

	#get the reader
	getBDTSexaqReader    = TMVA.Reader();

	#just define some variables
        var1 = array.array('f',[0])
        var2 = array.array('f',[0])
        var3 = array.array('f',[0])
        var4 = array.array('f',[0])
        var5 = array.array('f',[0])
        var6 = array.array('f',[0])
        var7 = array.array('f',[0])
        var8 = array.array('f',[0])
        var9 = array.array('f',[0])
        var10 = array.array('f',[0])
        var11 = array.array('f',[0])
        var12 = array.array('f',[0])
        var13 = array.array('f',[0])
        var14 = array.array('f',[0])
        var15 = array.array('f',[0])
        var16 = array.array('f',[0])
        var17 = array.array('f',[0])
        var18 = array.array('f',[0])
        var19 = array.array('f',[0])
        var20 = array.array('f',[0])
        var21 = array.array('f',[0])
	#add these variables to the reader 
        getBDTSexaqReader.AddVariable("_S_vz_interaction_vertex",           (var1))   
        getBDTSexaqReader.AddVariable("_S_lxy_interaction_vertex_beampipeCenter",           (var2))   
        getBDTSexaqReader.AddVariable("_S_daughters_deltaphi",      (var3))   
        getBDTSexaqReader.AddVariable("_S_daughters_deltaeta",      (var4))   
        getBDTSexaqReader.AddVariable("_S_daughters_openingsangle",      (var5))   
        getBDTSexaqReader.AddVariable("_S_daughters_DeltaR",      (var6))   
        getBDTSexaqReader.AddVariable("_S_Ks_openingsangle",      (var7))   
        getBDTSexaqReader.AddVariable("_S_Lambda_openingsangle",      (var8))   
        getBDTSexaqReader.AddVariable("_S_eta",     (var9))   
        getBDTSexaqReader.AddVariable("_Ks_eta",     (var10))   
        #getBDTSexaqReader.AddVariable("_Lambda_eta",     (var11))   
        getBDTSexaqReader.AddVariable("_S_dxy_over_lxy",          (var12))   
        getBDTSexaqReader.AddVariable("_Ks_dxy_over_lxy", (var13))   
        getBDTSexaqReader.AddVariable("_Lambda_dxy_over_lxy",         (var14))   
        getBDTSexaqReader.AddVariable("_S_dz_min",  (var15))   
        getBDTSexaqReader.AddVariable("_Ks_dz_min",  (var16))   
        getBDTSexaqReader.AddVariable("_Lambda_dz_min",          (var17))   
        getBDTSexaqReader.AddVariable("_Ks_pt",      (var18))   
        getBDTSexaqReader.AddVariable("_Lambda_lxy_decay_vertex",      (var19))   
        getBDTSexaqReader.AddVariable("_S_chi2_ndof",      (var20))   
        #getBDTSexaqReader.AddVariable("_S_pz",      (var21))   



	#book the weights from the training to the Reader
        getBDTSexaqReader.BookMVA("BDT","/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step1/"+dataset+"/weights/TMVAClassification_BDT.weights.xml")
	#for each file in the directory now evaluate the weights
	nEntriesinTreeSelectedTotal = 0
	iFile = 0
	for fileIn in sorted(os.listdir(input_directory_data)):
		#if(iFile not in [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]):
		#	iFile+=1
		#	continue
		if(not fileIn.endswith(".root")):continue
		print "the input file: ",fileIn

		fileIn = input_directory_data+"/"+fileIn
		fileH  = TFile.Open(fileIn)
		inTree = fileH.Get('FlatTree')
		#inTree = fileH.Get('FlatTreeProducerBDT/FlatTree')
		gROOT.cd()
		##The below line is very important. If I change Alt$(_S_charge,0) == 1 to Alt$(_S_charge,0) == -1 I am looking at signal and unblinding OMG

		#to select from the input trees the S bkg and apply the pre BDT cuts
		inTreeSelected = inTree.CopyTree(config_dict["config_SelectionBkgS"] + ' && ' + config_dict["config_pre_BDT_cuts"]  )

		nEntriesinTreeSelected = inTreeSelected.GetEntries()
                print "number of entries in the inTreeSelected: ", nEntriesinTreeSelected
		nEntriesinTreeSelectedTotal = nEntriesinTreeSelectedTotal + nEntriesinTreeSelected
		iFile+=1
	print "Total number of entries in the inTreeSelected :", nEntriesinTreeSelectedTotal

nOverlaps = 0
for e in overlapList:
	if(e[3] > 1):
		nOverlaps = nOverlaps + e[3]-1 #the overlap is the total times this events occurs (e[2]), with a -1 because one you should keep as the original one
		print e
print 'number of unique events: ', len(overlapList)
print 'number of overlaps: ', nOverlaps
