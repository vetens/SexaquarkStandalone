#from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, gROOT
import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

import sys
sys.path.insert(1, '../../../TMVA')
import configBDT as config

config_dict = config.config_dict

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = ""
tdrstyle.setTDRStyle()

colours = [4,2,35,1,38,41]

maxEvents = 1e99


#you can use this script to compared different collections of S and antiS to eachother. The most common one is comparing Data-S-BKG to MC-AntiS-Signal. The other option is comparing MC-S-BKG to MC-AntiS-BKG to Data-S-BKG


# Open file
SignFile= ROOT.TFile.Open(config_dict["config_SignalFile"]) 
BkgFile_Data_S  = ROOT.TFile.Open(config_dict["config_BkgFileData"]) 

#for when you have the BDT parameter in your tree:
#select the trees
# Get signal and background trees from file
SignalTree = SignFile.Get("FlatTreeProducerBDT/FlatTree")

BkgTree_Data_S = BkgFile_Data_S.Get("FlatTreeProducerBDT/FlatTree")


#now apply 1 by 1 the cuts and see how many events remain:
i_cut = 0
l_cut_applied = ""
for cut in ["pre_BDT_noCut","fiducial_region_cuts","pre_BDT_cut4","pre_BDT_cut1","pre_BDT_cut3"]:
	continue
	#first append the cut to l_cut_applied
	if(i_cut > 0): l_cut_applied = l_cut_applied + " && " + config_dict[cut]
	else: l_cut_applied = config_dict[cut]
	print "----> for cut = ", l_cut_applied
	print ""
	#I only want the antiS which match a GEN antiS in lxyz of the interaction vertex and the charge should also be negative
	gROOT.cd()
	selectedSignalTree = SignalTree.CopyTree(config_dict["config_SelectionSignalAntiS"]+ ' && ' + l_cut_applied)

	#for the MC signal the efficiency calculation is not so easy as for the background samples below as I have to reweigh the events
	nom = 0.
	for entry in  selectedSignalTree:
		nom += selectedSignalTree._S_event_weighting_factorALL[0]
	print "The remaining events for signal: ", nom


	gROOT.cd()
	#for selecting Data S BKG: ---> standard
	selectedBkgTree_Data_S_BKG = BkgTree_Data_S.CopyTree(config_dict["config_SelectionBkgS"] + ' && ' + l_cut_applied)
	print "The remaining events for background: ", selectedBkgTree_Data_S_BKG.GetEntries()
	print ""
	print ""
	i_cut += 1


#And now just evalueate the efficiency of the BDT parameter
SignFile = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_unblind_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing/DiscrApplied_combined_FlatTreeBDT_trial17AND21_1p8GeV_02112019_v1.root")
BkgFile_Data_S  = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")

SignalTree = SignFile.Get("FlatTree")
BkgTree_Data_S = BkgFile_Data_S.Get("FlatTree")

nom = 0
for entry in  SignalTree:
	nom += SignalTree._S_event_weighting_factorALL[0]
print '(just as a check: these are number of events in the signaltree and bkgtree in the tree with the BDT Discriminator attached (should be the same as the previous lines): ', nom, ' and ', BkgTree_Data_S.GetEntries()

print 'now applying cut on BDT classifier. The cut is at: ', config_dict["BDT_classifier_cut"]
nom = 0
for entry in  SignalTree:
                if(SignalTree.SexaqBDT>config_dict["BDT_classifier_cut"]): nom += SignalTree._S_event_weighting_factorALL[0]
print "The remaining events for signal: ", nom

BDT_cut = "SexaqBDT >"+ str(config_dict["BDT_classifier_cut"])
selectedBkgTree_Data_S = BkgTree_Data_S.CopyTree(BDT_cut)
print "The remaining events for background: ", selectedBkgTree_Data_S.GetEntries()



