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

colours = [4,2,35,30,43,1,6]
markerStyle = [20,21,22,23,33,34,35]

maxEvents = 1e99


#you can use this script to compared different collections of S and antiS to eachother. The most common one is comparing Data-S-BKG to MC-AntiS-Signal. The other option is comparing MC-S-BKG to MC-AntiS-BKG to Data-S-BKG
configuration = "all" # "MC-AntiS-Signal" "Data-S-BKG_to_MC-AntiS-Signal" or "MC-S-BKG_to_MC-AntiS-BKG_to_Data-S-BKG" or "all" 
#cut on the BDT parameter to select a minimal BDT, to go and look in the tail of the BDT distribution (towards the signal), put to -999 if you want all the events
min_BDT_classifier = -999. #-999 or -0.15


#################################################################
#With FiducialRegionCut applied
#################################################################
#MC_S_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_Cut/MC_S_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__31102019_v1.root")
#MC_S_Bkg_Tree = MC_S_Bkg_File.Get("FlatTree")
#
##MC AntiS BKG from DYJets sample
#MC_AntiS_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_Cut/MC_AntiS_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__31102019_v1.root")
#MC_AntiS_Bkg_Tree = MC_AntiS_Bkg_File.Get("FlatTree")
#
##Data S Bkg from SingleMuon Run2016H
#Data_S_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_Cut/Data_S_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
#Data_S_Bkg_Tree = Data_S_Bkg_File.Get("FlatTree")
#
##Data antiS Bkg from SingleMuon Run2016H, this is BKG because it is for a small BDT parameter
#Data_AntiS_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_Cut/Data_AntiS_Bkg/BDTApplied_partialUnblinding_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
#Data_AntiS_Bkg_Tree = Data_AntiS_Bkg_File.Get("FlatTree")
#
##The signal MC
#MC_AntiS_Sgn_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_Cut/MC_AntiS_Signal/BDTApplied_unblindMC_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_combined_FlatTreeBDT_trial17AND21_1p8GeV_02112019_v1.root")
#MC_AntiS_Sgn_Tree = MC_AntiS_Sgn_File.Get("FlatTree")
#
##Data antiS from SingleElectron Run2016H with antiS reconstructed X events
#Data_AntiS_XEvent_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_Cut/Data_AntiS_BKG_XEvent/BDTApplied_unblind_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_test_SingleMuonRun2016H_XEventAntiS.root")
#Data_AntiS_XEvent_Tree = Data_AntiS_XEvent_File.Get("FlatTree")

#################################################################
#With FiducialRegionCut and DeltaPhi cut applied
#################################################################
#MC_S_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Cut/MC_S_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__31102019_v1.root")
#MC_S_Bkg_Tree = MC_S_Bkg_File.Get("FlatTree")
#
##MC AntiS BKG from DYJets sample
#MC_AntiS_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Cut/MC_AntiS_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__31102019_v1.root")
#MC_AntiS_Bkg_Tree = MC_AntiS_Bkg_File.Get("FlatTree")
#
##Data S Bkg from SingleMuon Run2016H
#Data_S_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Cut/Data_S_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
#Data_S_Bkg_Tree = Data_S_Bkg_File.Get("FlatTree")
#
##Data antiS Bkg from SingleMuon Run2016H, this is BKG because it is for a small BDT parameter
#Data_AntiS_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Cut/Data_AntiS_Bkg/BDTApplied_partialUnblinding_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
#Data_AntiS_Bkg_Tree = Data_AntiS_Bkg_File.Get("FlatTree")
#
##The signal MC
#MC_AntiS_Sgn_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Cut/MC_AntiS_Signal/BDTApplied_unblindMC_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_combined_FlatTreeBDT_trial17AND21_1p8GeV_02112019_v1.root")
#MC_AntiS_Sgn_Tree = MC_AntiS_Sgn_File.Get("FlatTree")
#
##Data antiS from SingleElectron Run2016H with antiS reconstructed X events
#Data_AntiS_XEvent_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Cut/Data_AntiS_BKG_XEvent/BDTApplied_unblind_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_test_SingleMuonRun2016H_XEventAntiS.root")
#Data_AntiS_XEvent_Tree = Data_AntiS_XEvent_File.Get("FlatTree")



##################################################################
##With FiducialRegionCut and DeltaPhi cut and lxy cut applied
##################################################################
#MC_S_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Lxy_Cut/MC_S_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__31102019_v1.root")
#MC_S_Bkg_Tree = MC_S_Bkg_File.Get("FlatTree")
#
##MC AntiS BKG from DYJets sample
#MC_AntiS_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Lxy_Cut/MC_AntiS_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__31102019_v1.root")
#MC_AntiS_Bkg_Tree = MC_AntiS_Bkg_File.Get("FlatTree")
#
##Data S Bkg from SingleMuon Run2016H
#Data_S_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Lxy_Cut/Data_S_BKG/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
#Data_S_Bkg_Tree = Data_S_Bkg_File.Get("FlatTree")
#
##Data antiS Bkg from SingleMuon Run2016H, this is BKG because it is for a small BDT parameter
#Data_AntiS_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Lxy_Cut/Data_AntiS_Bkg/BDTApplied_partialUnblinding_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
#Data_AntiS_Bkg_Tree = Data_AntiS_Bkg_File.Get("FlatTree")
#
##The signal MC
#MC_AntiS_Sgn_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Lxy_Cut/MC_AntiS_Signal/BDTApplied_unblindMC_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_combined_FlatTreeBDT_trial17AND21_1p8GeV_02112019_v1.root")
#MC_AntiS_Sgn_Tree = MC_AntiS_Sgn_File.Get("FlatTree")
#
##Data antiS from SingleElectron Run2016H with antiS reconstructed X events
#Data_AntiS_XEvent_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied_CutFlow/FiducialRegion_DeltaPhi_Lxy_Cut/Data_AntiS_BKG_XEvent/BDTApplied_unblind_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapFalse/DiscrApplied_test_SingleMuonRun2016H_XEventAntiS.root")
#Data_AntiS_XEvent_Tree = Data_AntiS_XEvent_File.Get("FlatTree")


##################################################################
##With all pre BDT cuts applied
##################################################################
#Open files and trees:
#MC S BKG from DYJets sample
MC_S_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied/BkgReference/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapCheckFalse/ziscrApplied_FlatTreeBDT_MC-S-BKG-DYJets_blablablablablabl.root")
MC_S_Bkg_Tree = MC_S_Bkg_File.Get("FlatTree")

#MC AntiS BKG from DYJets sample
MC_AntiS_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied/BkgReference/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapCheckFalse/ziscrApplied_FlatTreeBDT_MC-AntiS-BKG-DYJets_blablablablablabl.root")
MC_AntiS_Bkg_Tree = MC_AntiS_Bkg_File.Get("FlatTree")

#Data S Bkg from SingleMuon Run2016H
Data_S_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied/BkgReference/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapCheckFalse/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
Data_S_Bkg_Tree = Data_S_Bkg_File.Get("FlatTree")

#Data antiS Bkg from SingleMuon Run2016H, this is BKG because it is for a small BDT parameter
Data_AntiS_Bkg_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied/partialUnblinding/BDTApplied_partialUnblinding_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapCheckFalse/DiscrApplied_FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
Data_AntiS_Bkg_Tree = Data_AntiS_Bkg_File.Get("FlatTree")

#The signal MC
MC_AntiS_Sgn_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied/Unblinded/SignalMC/DiscrApplied_combined_FlatTreeBDT_trial17AND21_1p8GeV_02112019_v1.root")
MC_AntiS_Sgn_Tree = MC_AntiS_Sgn_File.Get("FlatTree")

#Data antiS from SingleElectron Run2016H with antiS reconstructed X events
Data_AntiS_XEvent_File = ROOT.TFile.Open("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDTApplied/Unblinded/Data_AntiS_BKG_XEvent/DiscrApplied_test_SingleMuonRun2016H_XEventAntiS.root")
Data_AntiS_XEvent_Tree = Data_AntiS_XEvent_File.Get("FlatTree")

l_y_axis_ranges = [
0.08,
14.,
1.1,
1,
2.5,
2,
2,
4.5,
0.6,
0.6,
0.5,
15,
10,
12,
5,
10,
1,
1.2,
0.20,
10.,
0.05,
250.,
0.95,
50,
50
]

Legend = ["MC-S-BKG"   ,"MC-#bar{S}-BKG" ,"Data-S-BKG"   ,"Data-#bar{S}-BKG (BDT < 0.1)  ","Data-#bar{S}-X event BKG" ,"MC-#bar{S}-Signal"]
l_tree = [MC_S_Bkg_Tree,MC_AntiS_Bkg_Tree,Data_S_Bkg_Tree,Data_AntiS_Bkg_Tree             , Data_AntiS_XEvent_Tree    , MC_AntiS_Sgn_Tree]

plots_output_dir = "plots_BackgroundVsSignal/"+configuration+"/"

TH1_ll = [] #list of list of 1D histos 
TH2_ll = [] #list of list of 2D histos

iTree = 0
for inTree in l_tree:

	gROOT.cd()
	nEntries1 = inTree.GetEntries()
#	tree = inTree.CopyTree(config_dict["config_pre_BDT_cuts"])
	tree = inTree.CopyTree(config_dict["pre_BDT_noCut"])
	nEntries2 = tree.GetEntries()
	print "---------------------------------------------------"
	print "running for ", Legend[iTree]
	print "ratio of the number of events surviving the config_pre_BDT_cuts, just as a check. This should be = 1: ", float(nEntries2) ,"/",float(nEntries1), " = ", float(nEntries2)/float(nEntries1)

	h_S_vz_interaction_vertex= TH1F('h_S_vz_interaction_vertex','; absolute v_{z} iv ^{(}#bar{S} ^{)} (cm); 1/N_{ev} Events/3cm',20,-40,40)
	h_S_lxy_interaction_vertex = TH1F('h_S_lxy_interaction_vertex','; l_{0,bpc} iv ^{(}#bar{S} ^{)} (cm); 1/N_{ev} Events/0.2mm',31,1.9,2.52)

	h_S_daughters_deltaphi = TH1F('h_S_daughters_deltaphi','; #Delta#phi( ^{(} #bar{#Lambda} ^{)} ^{0} , K_{S}^{0}) (rad); 1/N_{ev} Events/0.2rad',32,-3.2,3.2)
	h_S_daughters_deltaeta = TH1F('h_S_daughters_deltaeta','; #Delta#eta( ^{(} #bar{#Lambda} ^{)} ^{0} , K_{S}^{0}) ; 1/N_{ev} Events/0.2rad',31,-3.1,3.1)
	h_S_daughters_openingsangle = TH1F('h_S_daughters_openingsangle','; openings angle( ^{(} #bar{#Lambda} ^{)} ^{0} , K_{S}^{0} ) (rad); 1/N_{ev} Events/0.1rad',35,0,3.5)
	h_S_daughters_DeltaR = TH1F('h_S_daughters_DeltaR','; #DeltaR( ^{(} #bar{#Lambda} ^{)} ^{0} , K_{S}^{0} ); 1/N_{ev} Events/0.1#DeltaR',60,0.5,6.5)
	h_S_Ks_openingsangle = TH1F('h_S_Ks_openingsangle','; openings angle( ^{(} #bar{S} ^{)} , K_{S}^{0}) (rad); 1/N_{ev} Events/0.1rad',20,0,2)
	h_S_Lambda_openingsangle = TH1F('h_S_Lambda_openingsangle','; openings angle( ^{(} #bar{#Lambda} ^{)} ^{0} ,  ^{(} #bar{S} ^{)} ) (rad); 1/N_{ev} Events/0.1rad',20,0,2)

	h_S_eta = TH1F('h_S_eta','; #eta( ^{(} #bar{S} ^{)} ); 1/N_{ev} Events/0.5#eta',20,-5,5)
	h_Ks_eta = TH1F('h_Ks_eta','; #eta(K_{S}^{0}) ; 1/N_{ev} Events/0.5#eta',20,-5,5)
	h_Lambda_eta = TH1F('h_Lambda_eta','; #eta( ^{(} #bar{#Lambda} ^{)} ^{0}); 1/N_{ev} Events/0.5#eta',20,-5,5)

	h_S_dxy_over_lxy = TH1F('h_S_dxy_over_lxy','; d_{0,bs}/l_{0,bs} ( ^{(} #bar{S} ^{)} ); 1/N_{ev} Events/0.1',20,-1,1)
	h_Ks_dxy_over_lxy = TH1F('h_Ks_dxy_over_lxy','; d_{0,bs}/l_{0,bs} (K_{S}^{0}); 1/N_{ev} Events/0.1',20,-1,1)
	h_Lambda_dxy_over_lxy = TH1F('h_Lambda_dxy_over_lxy','; d_{0,bs}/l_{0,bs} ( ^{(} #bar{#Lambda} ^{)} ^{0}) ; 1/N_{ev} Events/0.1',20,-1,1)

	h_S_dz_min = TH1F('h_S_dz_min','; min d_{z,bs}  ^{(} #bar{S} ^{)}  (cm); 1/N_{ev} Events/cm',20,-10,10)
	h_Ks_dz_min = TH1F('h_Ks_dz_min','; min d_{z,bs} K_{S}^{0} (cm); 1/N_{ev} Events/cm',60,-30,30)
	h_Lambda_dz_min = TH1F('h_Lambda_dz_min','; min d_{z,bs}  ^{(} #bar{#Lambda} ^{)} ^{0}  (cm); 1/N_{ev} Events/cm',60,-30,30)

	h_Ks_pt = TH1F('h_Ks_pt','; p_{t} K_{S}^{0} (GeV/c); 1/N_{ev} Events/0.4GeV/c',20,0,8)

	h_Lambda_lxy_decay_vertex = TH1F('h_Lambda_lxy_decay_vertex','; l_{0} ^{(} #bar{#Lambda} ^{)} ^{0} decay vertex (cm); 1/N_{ev} Events/cm',20,1.9,21.9)
	h_S_chi2_ndof = TH1F('h_S_chi2_ndof','; #chi^{2}/ndof ^{(}#bar{S} ^{)} annihilation vertex; 1/N_{ev} Events',44,0,11)

	h_S_pz = TH1F('h_S_pz','; p_{z}  ( ^{(} #bar{S} ^{)} ) (GeV/c); 1/N_{ev} Events/5GeV/c',16,-40,40)

	h_S_error_lxy_interaction_vertex = TH1F('h_S_error_lxy_interaction_vertex','; #sigma(l_{0,bpc} iv ^{(}#bar{S} ^{)} ) (cm); 1/N_{ev} Events/0.004mm',10,0,0.04)
	h_S_mass = TH1F('h_S_mass','; m_{ ^{(} #bar{S} ^{)} ,obs} (GeV/c^{2}); 1/N_{ev} Events/0.25GeV/c^{2}',40,-5,5)

	h_S_BDT = TH1F('h_S_BDT','; BDT classifier; 1/N_{ev} Events/0.05 BDT class.',40,-1,1)
	tprof_reweighing_factor = TProfile('tprof_reweighing_factor',';#eta ^{(}#bar{S} ^{)};reweighing factor',20,-5,5,0,50)

	for i in range(0,nEntries2):
		if(i==0): print 'charge of the S/antiS: ', tree._S_charge[0]
		if(i==maxEvents):
			break
		if(i%1e4 == 0):
			print "reached entry: ", i
		tree.GetEntry(i)
		
		#for when you have the BDT parameter in your tree:
		if(tree.SexaqBDT <= min_BDT_classifier):continue

		#need to reweigh the MC signal events, because the ones with high eta are more important, because they will pass more material
		reweighing_factor = config.calc_reweighing_factor(tree._S_eta[0],'MC-#bar{S}-Signal' in Legend[iTree])
		reweighing_factor = 1
		if('MC-#bar{S}-Signal' in Legend[iTree]): #if MC signal reweigh both for the pathlength through the beampipe and for the z location of the PV and PU
			reweighing_factor = tree._S_event_weighting_factor[0]*tree._S_event_weighting_factorPU[0]
		elif('MC-S-BKG' in Legend[iTree] or 'MC-#bar{S}-BKG' in Legend[iTree]):#if MC background only reweigh for the z location of the PV and PU
			reweighing_factor = tree._S_event_weighting_factorPU[0]
			#reweighing_factor = 1.


		h_S_vz_interaction_vertex.Fill(tree._S_vz_interaction_vertex[0],reweighing_factor)
		h_S_lxy_interaction_vertex.Fill(tree._S_lxy_interaction_vertex_beampipeCenter[0],reweighing_factor)

		h_S_daughters_deltaphi.Fill(tree._S_daughters_deltaphi[0],reweighing_factor)
		h_S_daughters_deltaeta.Fill(tree._S_daughters_deltaeta[0],reweighing_factor)
		h_S_daughters_openingsangle.Fill(tree._S_daughters_openingsangle[0],reweighing_factor)
		h_S_daughters_DeltaR.Fill(tree._S_daughters_DeltaR[0],reweighing_factor)
		h_S_Ks_openingsangle.Fill(tree._S_Ks_openingsangle[0],reweighing_factor)
		h_S_Lambda_openingsangle.Fill(tree._S_Lambda_openingsangle[0],reweighing_factor)

		h_S_eta.Fill(tree._S_eta[0],reweighing_factor)
		h_Ks_eta.Fill(tree._Ks_eta[0],reweighing_factor)
		h_Lambda_eta.Fill(tree._Lambda_eta[0],reweighing_factor)

		h_S_dxy_over_lxy.Fill(tree._S_dxy_over_lxy[0],reweighing_factor)
		h_Ks_dxy_over_lxy.Fill(tree._Ks_dxy_over_lxy[0],reweighing_factor)
		h_Lambda_dxy_over_lxy.Fill(tree._Lambda_dxy_over_lxy[0],reweighing_factor)

		h_S_dz_min.Fill(tree._S_dz_min[0],reweighing_factor)
		h_Ks_dz_min.Fill(tree._Ks_dz_min[0],reweighing_factor)
		h_Lambda_dz_min.Fill(tree._Lambda_dz_min[0],reweighing_factor)

		h_Ks_pt.Fill(tree._Ks_pt[0],reweighing_factor)
		
		h_Lambda_lxy_decay_vertex.Fill(tree._Lambda_lxy_decay_vertex[0],reweighing_factor)
		h_S_chi2_ndof.Fill(tree._S_chi2_ndof[0],reweighing_factor)

		h_S_pz.Fill(tree._S_pz[0],reweighing_factor)

		h_S_error_lxy_interaction_vertex.Fill(tree._S_error_lxy_interaction_vertex_beampipeCenter[0],reweighing_factor)  
		h_S_mass.Fill(tree._S_mass[0],reweighing_factor)
		h_S_BDT.Fill(tree.SexaqBDT,reweighing_factor)

		tprof_reweighing_factor.Fill(tree._S_eta[0],reweighing_factor)	

	TH1_l = [h_S_vz_interaction_vertex,h_S_lxy_interaction_vertex,h_S_daughters_deltaphi,h_S_daughters_deltaeta,h_S_daughters_openingsangle,h_S_daughters_DeltaR,h_S_Ks_openingsangle,h_S_Lambda_openingsangle,h_S_eta,h_Ks_eta,h_Lambda_eta,h_S_dxy_over_lxy,h_Ks_dxy_over_lxy,h_Lambda_dxy_over_lxy,h_S_dz_min,h_Ks_dz_min,h_Lambda_dz_min,h_Ks_pt,h_Lambda_lxy_decay_vertex,h_S_chi2_ndof,h_S_pz,h_S_error_lxy_interaction_vertex,h_S_mass,h_S_BDT,tprof_reweighing_factor]
	for h in TH1_l:
		h.SetDirectory(0) 
	TH1_ll.append(TH1_l)

	TH2_l = []
	for h in TH2_l:
		h.SetDirectory(0) 
	TH2_ll.append(TH2_l)

	iTree+=1

fOut = TFile(plots_output_dir+'macro_FlatTree_BDT_trial17.root','RECREATE')

nHistos = len(TH1_ll[0])
nSamples = len(TH1_ll)

c_all  = TCanvas("c_all","c_all",1800,1400)
c_all_log  = TCanvas("c_all_log","c_all_log",1600,1200)
c_all.Divide(6,5)
for i in range(0,nHistos):#each list contains a list of histograms. Each list represents background or signal. The histos need to be overlaid one list to the other
	h = TH1_ll[0][i]
	c_name = "c_"+h.GetName()
	c = TCanvas(c_name,"")
	legend = TLegend(0.65,0.77,0.99,0.99)
	for j in range(0,nSamples):
		h = TH1_ll[j][i]
		if(h.GetName() != "tprof_reweighing_factor"):
			if(h.GetSumw2N() == 0):
				h.Sumw2(kTRUE)
			if(h.Integral() != 0):
				h.Scale(1./h.Integral(), "width");
		h.GetYaxis().SetRangeUser(0.,l_y_axis_ranges[i])
		if("dxy_over_lxy" in h.GetName()):
			h.GetYaxis().SetRangeUser(1e-2,l_y_axis_ranges[i])
		if("dz_min" in h.GetName()):
			h.GetYaxis().SetRangeUser(1e-5,l_y_axis_ranges[i])
		if("chi2_ndof" in h.GetName()):
			h.GetYaxis().SetRangeUser(1e-3,l_y_axis_ranges[i])
		if("BDT" in h.GetName()):
			h.GetYaxis().SetRangeUser(1e-3,l_y_axis_ranges[i])
		if j == 0:
			h.Draw("L")
			c.Update()
		else:
			h.Draw("PCE1same")
			c.Update()
		if("dxy_over_lxy" in h.GetName() or "dz_min" in h.GetName() or "chi2_ndof" in h.GetName() or "BDT" in h.GetName()):
			c.SetLogy()
		h.SetLineColor(colours[j])
		h.SetLineWidth(2)
		h.SetMarkerStyle(markerStyle[j])
		h.SetMarkerColor(colours[j])
		h.SetMarkerSize(0.8)
		h.SetStats(0)
		legend.AddEntry(h, Legend[j],"ep")
	legend.Draw()
	CMS_lumi.CMS_lumi(c, 0, 11)
	c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
	c.Write()

	c_all.cd(i+1)
	c.SetGridx(0)
	c.SetGridy(0)
	c.DrawClonePad()
	

c_all.Write()
c_all.SaveAs(plots_output_dir+"c_all.pdf")
c_all_log.Write()
c_all_log.SaveAs(plots_output_dir+"c_all_log.pdf")


i = 0
for l in TH2_ll:
	for h in l:
		c_name = "c_"+h.GetName()+"_"+Legend[i]
		c = TCanvas(c_name,"");
		c.SetRightMargin(0.2) #make room for the tile of the z scale
		if(h.GetSumw2N() == 0):
			h.Sumw2(kTRUE)
		h.Scale(1./h.Integral(), "width");
		h.Draw("colz")
		h.SetStats(0)
		CMS_lumi.CMS_lumi(c, 0, 11)
		c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
		c.Write()
	i += 1

#c_antiS_eta_pt = TCanvas("c_antiS_eta_pt","")
#h_antiS_eta_pt.DrawNormalized()
#c_antiS_eta_pt.Write()
#
#c_antiS_eta_pz = TCanvas("c_antiS_eta_pz","")
#h_antiS_eta_pz.DrawNormalized()
#c_antiS_eta_pz.Write()

fOut.Write()
fOut.Close()
