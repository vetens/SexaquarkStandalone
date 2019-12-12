import ROOT
from ROOT import *

from os import environ
#give the output classifier a meaningful name
version = "dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing"

import sys
sys.path.insert(1, './..')
import configBDT as config

config_dict = config.config_dict

# Open files
SignFile1 = ROOT.TFile.Open(config_dict["config_SignalFile"])

#BkgFile  = ROOT.TFile.Open("/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/Results/FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
BkgFile  = ROOT.TFile.Open(config_dict["config_BkgFileData"])

# Get signal and background trees from file

SignalTree1     = SignFile1.Get("FlatTreeProducerBDT/FlatTree")
#apply the pre-BDT cuts and also select from the ntuples the Sbar (based on charge) and should be signal (based on matching with GEN particles)
gROOT.cd()
selectedSignalTree1 = SignalTree1.CopyTree( config_dict["config_SelectionSignalAntiS"] + ' && ' + config_dict["config_pre_BDT_cuts"] )


BkgTree        = BkgFile.Get("FlatTreeProducerBDT/FlatTree")
gROOT.cd()
#apply pre-BDT cuts and select data BKG S (look at the charge)
selectedBkgTree = BkgTree.CopyTree(config_dict["config_SelectionBkgS"] + ' && ' + config_dict["config_pre_BDT_cuts"] )

trainTestSplit = 0.8

# Add variables to dataloader
dataloader = ROOT.TMVA.DataLoader('dataset_BDT_2016'+version) 
#dataloader.AddVariable("_S_error_lxy_interaction_vertex_beampipeCenter") #selected
dataloader.AddVariable("_S_vz_interaction_vertex") 
dataloader.AddVariable("_S_lxy_interaction_vertex_beampipeCenter") 

dataloader.AddVariable("_S_daughters_deltaphi")
dataloader.AddVariable("_S_daughters_deltaeta") 
dataloader.AddVariable("_S_daughters_openingsangle")
dataloader.AddVariable("_S_daughters_DeltaR") 
dataloader.AddVariable("_S_Ks_openingsangle") 
dataloader.AddVariable("_S_Lambda_openingsangle") 

dataloader.AddVariable("_S_eta") 
dataloader.AddVariable("_Ks_eta") 
#dataloader.AddVariable("_Lambda_eta")

dataloader.AddVariable("_S_dxy_over_lxy") 
dataloader.AddVariable("_Ks_dxy_over_lxy") 
dataloader.AddVariable("_Lambda_dxy_over_lxy") 

#dataloader.AddVariable("_S_dxy_dzPVmin")
#dataloader.AddVariable("_Ks_dxy_dzPVmin")
#dataloader.AddVariable("_Lambda_dxy_dzPVmin")
#dataloader.AddVariable("_S_dxy")
#dataloader.AddVariable("_Ks_dxy")
#dataloader.AddVariable("_Lambda_dxy")

dataloader.AddVariable("_S_dz_min")
dataloader.AddVariable("_Ks_dz_min") 
dataloader.AddVariable("_Lambda_dz_min") 

#dataloader.AddVariable("_S_pt") 
dataloader.AddVariable("_Ks_pt")
#dataloader.AddVariable("_Lambda_pt") 

dataloader.AddVariable("_Lambda_lxy_decay_vertex")
dataloader.AddVariable("_S_chi2_ndof")
#dataloader.AddVariable("_S_pz")


# Add trees to dataloader
dataloader.AddSignalTree(selectedSignalTree1, 1)
dataloader.AddBackgroundTree(selectedBkgTree, 1)

#do the event by event reweighing for signal (more weight goes to the events with larger distance to travel through the beampipe and reweighing also done for different PV distribution)
dataloader.SetSignalWeightExpression("_S_event_weighting_factorALL")


dataloader.PrepareTrainingAndTestTree(ROOT.TCut(config_dict["config_pre_BDT_cuts"]),\
	'TrainTestSplit_Signal={}:'.format(trainTestSplit)+\
	'TrainTestSplit_Background={}:'.format(trainTestSplit)+'SplitMode=Random')


# Setup TMVA
ROOT.TMVA.Tools.Instance()
ROOT.TMVA.PyMethodBase.PyInitialize()

outputFile = ROOT.TFile.Open('BDTOutput_2016_'+version+'.root', 'RECREATE')
factory = ROOT.TMVA.Factory('TMVAClassification', outputFile,
        '!V:!Silent:Color:Transformations=I:'+\
        'AnalysisType=Classification')


# BDT method
factory.BookMethod(dataloader,'BDT', 'BDT',
                'H:!V:VarTransform=None:'+\
                'NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=12:UseBaggedBoost=True')
factory.TrainAllMethods()

factory.TestAllMethods()

factory.EvaluateAllMethods()

canvas = factory.GetROCCurve(dataloader)
canvas.Draw()
canvas.SaveAs("BDT_2016_"+version+".root")
