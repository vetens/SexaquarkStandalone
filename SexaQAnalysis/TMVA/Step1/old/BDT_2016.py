import ROOT
from ROOT import *

# Select Theano as backend for Keras
from os import environ
version = "v_18_variables_AdaBoost_04092019_MCBkg_noInvMassCutInTreeProduction_NoInvMassCutInBDTScript"

# Open file
#SignFile = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Analyzed_Skimmed/CRAB_AnalyzerAllSkimmed_WithPU2016NeutrinoGun_tryToFix_8_10072019_v1/crab_AnalyzerAllStepsSkimmedSexaqWithPU2016NeutrinoGun_tryToFix_8_10072019_v1/190710_104816/FlatTree_WithPU2016NeutrinoGun_tryToFix_8_10072019_v1.root")
SignFile1 = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaqWithPU2016NeutrinoGun_tryToFix_8_10072019_v1/FlatTree_Skimmed_trial10.root")
SignFile2 = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaqWithPU2016NeutrinoGun_tryToFix_8_trial11_17072019_v1/crab_SkimmingSexaqWithPU2016NeutrinoGun_tryToFix_8_trial11_17072019_v1/190717_083219/FlatTree/FlatTree_Skimmed_trial11.root")
SignFile3 = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaq_Skimming_trial13_25072019_v1/crab_SkimmingSexaq_trial13_25072019_v1/190725_051522/FlatTree_Skimmed_trial13.root")
SignFile4 = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_Step1Sexaq_trial14/crab_Step1_Step2_Skimming_FlatTree_Sexaq_trial14_04082019_v1/190804_115510/FlatTree_Skimmed_trial14.root")
SignFile5 = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaq_trial16/crab_Step1_Step2_Skimming_FlatTree_trial16_26082019_v1/190826_190802/combined_FlatTree_Skimmed_trial16.root")

#BkgFile  = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SingleElectron/FlatTree_SingleElectron2016_Background.root")
#BkgFile  = ROOT.TFile.Open("/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SingleMuon/FlatTree_SingleMuonRun2016H_Background.root")
#Try Bkg sample without applying the inv_mass > 0 cut in the filling of the tree:
BkgFile  = ROOT.TFile.Open("/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducer/FlatTree_SingleMuonRun2016H_Background_NoMinMassCut.root")
#Try as Bkg some background extracted from MC:
#BkgFile  = ROOT.TFile.Open("/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducer/test_trial11_S_candidates_inMC.root")


# Get signal and background trees from file
SignalTree1     = SignFile1.Get("FlatTreeProducer/FlatTree")
SignalTree2     = SignFile2.Get("FlatTreeProducer/FlatTree")
SignalTree3     = SignFile3.Get("FlatTreeProducer/FlatTree")
SignalTree4     = SignFile4.Get("FlatTreeProducer/FlatTree")
SignalTree5     = SignFile5.Get("FlatTreeProducer/FlatTree")



BkgTree        = BkgFile.Get("FlatTreeProducer/FlatTree")

# Add variables to dataloader
dataloader = ROOT.TMVA.DataLoader('dataset_BDT_2016') 
#dataloader.AddVariable("_S_error_lxy_interaction_vertex") #selected
dataloader.AddVariable("_Ks_vz_decay_vertex") #selected  --> might still be interesting
dataloader.AddVariable("_S_lxy_interaction_vertex") #selected


dataloader.AddVariable("_S_daughters_deltaphi")
dataloader.AddVariable("_S_daughters_deltaeta") # selected
dataloader.AddVariable("_S_daughters_openingsangle")
dataloader.AddVariable("_S_daughters_DeltaR") 
dataloader.AddVariable("_S_Ks_openingsangle")
dataloader.AddVariable("_S_Lambda_openingsangle")

dataloader.AddVariable("_S_eta") 
dataloader.AddVariable("_Ks_eta") # selected
dataloader.AddVariable("_Lambda_eta")

dataloader.AddVariable("_S_dxy_over_lxy") #selected
dataloader.AddVariable("_Ks_dxy_over_lxy") #selected
dataloader.AddVariable("_Lambda_dxy_over_lxy") #selected


#don't use following dxy variables as dxy_over_lxy seems the one which is most discriminating
#dataloader.AddVariable("_S_dxy_dzPVmin")
#dataloader.AddVariable("_Ks_dxy_dzPVmin")
#dataloader.AddVariable("_Lambda_dxy_dzPVmin")
#dataloader.AddVariable("_S_dxy")
#dataloader.AddVariable("_Ks_dxy")
#dataloader.AddVariable("_Lambda_dxy")

dataloader.AddVariable("_S_dz_min")
dataloader.AddVariable("_Ks_dz_min") # selected
dataloader.AddVariable("_Lambda_dz_min") #selected

#dataloader.AddVariable("_S_pt") 
dataloader.AddVariable("_Ks_pt")# --> might still be interesting 
#dataloader.AddVariable("_Lambda_pt") 



# Add trees to dataloader
dataloader.AddSignalTree(SignalTree1, 1)
dataloader.AddSignalTree(SignalTree2, 1)
dataloader.AddSignalTree(SignalTree3, 1)
dataloader.AddSignalTree(SignalTree4, 1)
dataloader.AddSignalTree(SignalTree5, 1)
dataloader.AddBackgroundTree(BkgTree, 1)
trainTestSplit = 0.8

MasterCut = ROOT.TCut("\
Alt$(_S_error_lxy_interaction_vertex,0) < 0.1 &&\
(Alt$(_S_lxy_interaction_vertex,0) > 1.9 && Alt$(_S_lxy_interaction_vertex,0) < 12) &&\
Alt$(_S_chi2_ndof,0) < 4. && \
(Alt$(_S_daughters_deltaphi,0) < -0.6 || Alt$(_S_daughters_deltaphi,0) > 0.6) && \
Alt$(_S_daughters_openingsangle,0) < 1.6 && \
(Alt$(_S_daughters_deltaeta,0) > -2 && Alt$(_S_daughters_deltaeta,0) < 2) && \
Alt$(_S_Ks_openingsangle,0) < 1.4 && \
Alt$(_S_Lambda_openingsangle,0) < 1 && Alt$(_S_daughters_DeltaR,0) < 3.5 && \
(Alt$(_S_dxy_over_lxy,0) > 0 && Alt$(_S_dxy_over_lxy,0) < 0.25) && \
(Alt$(_S_dz_min,0) > -5 && Alt$(_S_dz_min,0) < 5)"
)
#Alt$(_S_mass,0) > 0. && Alt$(_S_chi2_ndof,0) < 4. && \
#Alt$(_S_daughters_deltaeta,0) < 2.5 && Alt$(_S_daughters_deltaeta,0) > -2.5 && \
#MasterCut = ROOT.TCut("Alt$(_S_error_lxy_interaction_vertex,0) < 0.1 && Alt$(_S_lxy_interaction_vertex,0) > 1.9 && Alt$(_S_mass,0) > 0. && Alt$(_S_chi2_ndof,0) < 4." )

dataloader.PrepareTrainingAndTestTree(MasterCut,\
	'TrainTestSplit_Signal={}:'.format(trainTestSplit)+\
	'TrainTestSplit_Background={}:'.format(trainTestSplit)+'SplitMode=Random')


# Setup TMVA
ROOT.TMVA.Tools.Instance()
ROOT.TMVA.PyMethodBase.PyInitialize()

outputFile = ROOT.TFile.Open('BDTOutput_2016_'+version+'.root', 'RECREATE')
factory = ROOT.TMVA.Factory('TMVAClassification', outputFile,
        '!V:!Silent:Color:Transformations=I;D;P;G,D:'+\
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
