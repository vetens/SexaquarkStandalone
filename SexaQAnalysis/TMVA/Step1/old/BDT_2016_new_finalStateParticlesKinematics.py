import ROOT
from ROOT import *

# Select Theano as backend for Keras
from os import environ
version = "vSelected14Parameters_CutLxy_CutErrorLxy_CutDeltaPhi_CutDxyOverLxy_CutTightMaxLxyInteractionVertex_newParameters_LambdaLxyDecayVertex_SChi2Ndof_Spz"

# Open file
#SignFile1 = ROOT.TFile.Open("/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/test_FlatTreeBDT_trial15.root")
#SignFile1 = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/lowette/crmc_Sexaq/Skimmed/CRAB_SimSexaq_trial17/crab_Step1_Step2_Skimming_FlatTree_trial17_18092019_v1/190918_051631/combined_FlatTreeBDT_Skimmed_trial17_21.root")
SignFile1 = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/lowette/crmc_Sexaq/Skimmed/CRAB_SimSexaq_trial17/crab_Step1_Step2_Skimming_FlatTree_trial17_18092019_v1/190918_051631/combined_FlatTreeBDT_Skimmed_trial17_21_v2.root")

#BkgFile  = ROOT.TFile.Open("/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/Results/FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")
BkgFile  = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/ALL_v2/FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root")

# Get signal and background trees from file
SignalTree1     = SignFile1.Get("FlatTreeProducerBDT/FlatTree")
#I only want the antiS which match a GEN antiS in lxyz of the interaction vertex and the charge should also be negative
gROOT.cd()
selectedSignalTree1 = SignalTree1.CopyTree('Alt$(_S_charge,0) == -1 && Alt$(_S_deltaLInteractionVertexAntiSmin,0) < 0.5 && Alt$(_S_lxy_interaction_vertex,0) < 2.5')


BkgTree        = BkgFile.Get("FlatTreeProducerBDT/FlatTree")
#I only want the S with positive charge
gROOT.cd()
#for selecting Data S BKG: ---> standard
selectedBkgTree = BkgTree.CopyTree('Alt$(_S_charge,0) == 1 && Alt$(_S_lxy_interaction_vertex,0) < 2.5')

trainTestSplit = 0.8

MasterCut = ROOT.TCut("Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) &&  Alt$(_S_dxy_over_lxy,0) >= 0 &&\
Alt$(_Lambda_vz_decay_vertex,0) >= -105 && Alt$(_Lambda_vz_decay_vertex,0) <= 105 && Alt$(_Lambda_lxy_decay_vertex,0) <= 49.5 &&\
Alt$(_Ks_vz_decay_vertex,0) >= -105 && Alt$(_Ks_vz_decay_vertex,0) <= 105 && Alt$(_Ks_lxy_decay_vertex,0) <= 49.5 &&\
Alt$(_RECO_Lambda_daughter0_pt,0) >= 0.35 && Alt$(_RECO_Lambda_daughter0_pz,0) <= 22 && Alt$(_RECO_Lambda_daughter0_dxy_beamspot,0) >= 0 && Alt$(_RECO_Lambda_daughter0_dxy_beamspot,0) <= 8.5 && Alt$(_RECO_Lambda_daughter0_dz_beamspot,0) >= -27 && Alt$(_RECO_Lambda_daughter0_dz_beamspot,0) <= 27 &&\
Alt$(_RECO_Lambda_daughter1_pt,0) >= 0.35 && Alt$(_RECO_Lambda_daughter1_pz,0) <= 22 && Alt$(_RECO_Lambda_daughter1_dxy_beamspot,0) >= 0 && Alt$(_RECO_Lambda_daughter1_dxy_beamspot,0) <= 8.5 && Alt$(_RECO_Lambda_daughter1_dz_beamspot,0) >= -27 && Alt$(_RECO_Lambda_daughter1_dz_beamspot,0) <= 27 &&\
Alt$(_RECO_Ks_daughter0_pt,0) >= 0.35 && Alt$(_RECO_Ks_daughter0_pz,0) <= 22 && Alt$(_RECO_Ks_daughter0_dxy_beamspot,0) >= 0 && Alt$(_RECO_Ks_daughter0_dxy_beamspot,0) <= 8.5 && Alt$(_RECO_Ks_daughter0_dz_beamspot,0) >= -27 && Alt$(_RECO_Ks_daughter0_dz_beamspot,0) <= 27 &&\
Alt$(_RECO_Ks_daughter1_pt,0) >= 0.35 && Alt$(_RECO_Ks_daughter1_pz,0) <= 22 && Alt$(_RECO_Ks_daughter1_dxy_beamspot,0) >= 0 && Alt$(_RECO_Ks_daughter1_dxy_beamspot,0) <= 8.5 && Alt$(_RECO_Ks_daughter1_dz_beamspot,0) >= -27 && Alt$(_RECO_Ks_daughter1_dz_beamspot,0) <= 27")



# Add variables to dataloader
dataloader = ROOT.TMVA.DataLoader('dataset_BDT_2016') 
#decay vertex Lambda
dataloader.AddVariable("_Lambda_vz_decay_vertex") 
dataloader.AddVariable("_Lambda_lxy_decay_vertex") 
#decay vertex Ks
dataloader.AddVariable("_Ks_vz_decay_vertex") 
dataloader.AddVariable("_Ks_lxy_decay_vertex") 
#track1 Lambda
dataloader.AddVariable("_RECO_Lambda_daughter0_pt") 
dataloader.AddVariable("_RECO_Lambda_daughter0_pz") 
dataloader.AddVariable("_RECO_Lambda_daughter0_dxy_beamspot") 
dataloader.AddVariable("_RECO_Lambda_daughter0_dz_beamspot") 
#track2 Lambda
dataloader.AddVariable("_RECO_Lambda_daughter1_pt") 
dataloader.AddVariable("_RECO_Lambda_daughter1_pz") 
dataloader.AddVariable("_RECO_Lambda_daughter1_dxy_beamspot") 
dataloader.AddVariable("_RECO_Lambda_daughter1_dz_beamspot") 
#track1 Ks
dataloader.AddVariable("_RECO_Ks_daughter0_pt") 
dataloader.AddVariable("_RECO_Ks_daughter0_pz") 
dataloader.AddVariable("_RECO_Ks_daughter0_dxy_beamspot") 
dataloader.AddVariable("_RECO_Ks_daughter0_dz_beamspot") 
#track2 Ks
dataloader.AddVariable("_RECO_Ks_daughter1_pt") 
dataloader.AddVariable("_RECO_Ks_daughter1_pz") 
dataloader.AddVariable("_RECO_Ks_daughter1_dxy_beamspot") 
dataloader.AddVariable("_RECO_Ks_daughter1_dz_beamspot") 

# Add trees to dataloader
dataloader.AddSignalTree(selectedSignalTree1, 1)
dataloader.AddBackgroundTree(selectedBkgTree, 1)

dataloader.PrepareTrainingAndTestTree(MasterCut,\
	'TrainTestSplit_Signal={}:'.format(trainTestSplit)+\
	'TrainTestSplit_Background={}:'.format(trainTestSplit)+'SplitMode=Random')


# Setup TMVA
ROOT.TMVA.Tools.Instance()
ROOT.TMVA.PyMethodBase.PyInitialize()

outputFile = ROOT.TFile.Open('BDTOutput_2016_'+version+'.root', 'RECREATE')
factory = ROOT.TMVA.Factory('TMVAClassification', outputFile,
        '!V:!Silent:Color:Transformations=I:'+\
        'AnalysisType=Classification')

#'!V:!Silent:Color:Transformations=I;D;P;G,D:'+\



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
