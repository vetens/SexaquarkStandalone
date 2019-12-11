import numpy as np
from ROOT import *
import pandas as pd

import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle 

sys.path.insert(1, '/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA')
import configBDT as config
config_dict = config.config_dict

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

colours = [1,2,4,35,38,41]

maxEvents = 1e99

verbose = False

plots_output_dir = "plots_syst_evaluation/"

#inFiles = [TFile("/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerTracking/test_FlatTreeSkimming_Step1_Step2_Skimming_FlatTree_trial17_1p8GeV_17102019_v1_191017_220444_numberOfTrackerHits.root",'read')]
inFiles = [TFile("file:/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/FlatTree_Skimmed/CRAB_SimSexaq_trial21/crab_FlatTreeProducerTracking_trial21_02112019_v1_1p8GeV/191102_062811/0000/combined/combined_FlatTreeTracking_trial21_02112019_v1_1p8GeV.root",'read')]

fOut = TFile(plots_output_dir+'macro_syst_evaluation_antiS_RECO_eff2.root','RECREATE')

def printProgress(i):
	if(i%10000 == 0):
		print 'reached track: ', i, ': ', float(i)/float(min(maxEvents,tree.GetEntries()))*100, '%'

#define the fiducial region
FidReg_minPt = config_dict["config_fidRegion_FiducialRegionptMin"] 
FidReg_maxPt = config_dict["config_fidRegion_FiducialRegionptMax"]

FidReg_minPz = config_dict["config_fidRegion_FiducialRegionpzMin"]
FidReg_maxPz = config_dict["config_fidRegion_FiducialRegionpzMax"]

FidReg_minvz = config_dict["config_fidRegion_FiducialRegionvzMin"]
FidReg_maxvz = config_dict["config_fidRegion_FiducialRegionvzMax"]

FidReg_maxlxy = config_dict["config_fidRegion_FiducialRegionlxyMax"]

FidReg_mindxy = config_dict["config_fidRegion_FiducialRegiondxyMin"]
FidReg_maxdxy = config_dict["config_fidRegion_FiducialRegiondxyMax"]

FidReg_mindz = config_dict["config_fidRegion_FiducialRegiondzMin"]
FidReg_maxdz = config_dict["config_fidRegion_FiducialRegiondzMax"]


#the files with the correction factors
corr_factors_pt = pd.read_csv("../MCToData/Data_MC_plots_8_final/h_RECO_Ks_pt_tracks1_and_2.dat") 
corr_factors_pz = pd.read_csv("../MCToData/Data_MC_plots_8_final/h_RECO_Ks_pz_tracks1_and_2.dat") 
corr_factors_dxy = pd.read_csv("../MCToData/Data_MC_plots_8_final/h_RECO_Ks_Track1_and_2_dxy_beamspot.dat") 
corr_factors_dz = pd.read_csv("../MCToData/Data_MC_plots_8_final/h_RECO_Ks_Track1_and_2_dz_min_PV.dat") 
corr_factors_lxy = pd.read_csv("../MCToData/Data_MC_plots_8_final/h_RECO_Ks_lxy.dat") 
corr_factors_vz = pd.read_csv("../MCToData/Data_MC_plots_8_final/h_RECO_Ks_vz.dat") 

#histograms:
h_corr_factor_pt = TH1F("corr_factor_pt",";correction factor pt;",100,0.2,1.8)
h_corr_factor_pz = TH1F("corr_factor_pz",";correction factor pz;",100,0.2,1.8)
h_corr_factor_dxy = TH1F("corr_factor_dxy",";correction factor dxy;",100,0.2,1.8)
h_corr_factor_dz = TH1F("corr_factor_dz",";correction factor dz;",100,0.2,1.8)
h_corr_factor_lxy = TH1F("corr_factor_lxy",";correction factor lxy;",100,0.2,1.8)
h_corr_factor_vz = TH1F("corr_factor_vz",";correction factor vz;",100,0.2,1.8)
h_corr_factor_this_antiS = TH1F("h_corr_factor_this_antiS",";correction factor (C_{k});#Events;",16,0.2,1.8)

#for each kinematic variable plot the correction factor of that kinematic variable versus the correction parameters of the others.
nbins_corr_par = 60
min_bin_corr_par = 0.7
max_bin_corr_par = 1.3
#pt against the others
h_corr_pt_pt = TH2F("h_corr_pt_pt",";pt correction factor; pt correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pt_pz = TH2F("h_corr_pt_pz",";pt correction factor; pz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pt_dxy = TH2F("h_corr_pt_dxy",";pt correction factor; dxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pt_dz = TH2F("h_corr_pt_dz",";pt correction factor; dz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pt_lxy = TH2F("h_corr_pt_lxy",";pt correction factor; lxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pt_vz = TH2F("h_corr_pt_vz",";pt correction factor; vz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)

#pz against the others
h_corr_pz_pt = TH2F("h_corr_pz_pt",";pz correction factor;pt correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pz_pz = TH2F("h_corr_pz_pz",";pz correction factor;pz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pz_dxy = TH2F("h_corr_pz_dxy",";pz correction factor;dxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pz_dz = TH2F("h_corr_pz_dz",";pz correction factor;dz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pz_lxy = TH2F("h_corr_pz_lxy",";pz correction factor;lxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_pz_vz = TH2F("h_corr_pz_vz",";pz correction factor;vz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)

#dxy against the others
h_corr_dxy_pt = TH2F("h_corr_dxy_pt",";dxy correction factor;pt correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dxy_pz = TH2F("h_corr_dxy_pz",";dxy correction factor;pz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dxy_dxy = TH2F("h_corr_dxy_dxy",";dxy correction factor;dxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dxy_dz = TH2F("h_corr_dxy_dz",";dxy correction factor;dz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dxy_lxy = TH2F("h_corr_dxy_lxy",";dxy correction factor;lxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dxy_vz = TH2F("h_corr_dxy_vz",";dxy correction factor;vz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)

#dz against the others
h_corr_dz_pt = TH2F("h_corr_dz_pt",";dz correction factor;pt correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dz_pz = TH2F("h_corr_dz_pz",";dz correction factor;pz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dz_dxy = TH2F("h_corr_dz_dxy",";dz correction factor;dxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dz_dz = TH2F("h_corr_dz_dz",";dz correction factor;dz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dz_lxy = TH2F("h_corr_dz_lxy",";dz correction factor;lxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_dz_vz = TH2F("h_corr_dz_vz",";dz correction factor;vz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)

#lxy against the others
h_corr_lxy_pt = TH2F("h_corr_lxy_pt",";lxy correction factor;pt correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_lxy_pz = TH2F("h_corr_lxy_pz",";lxy correction factor;pz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_lxy_dxy = TH2F("h_corr_lxy_dxy",";lxy correction factor;dxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_lxy_dz = TH2F("h_corr_lxy_dz",";lxy correction factor;dz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_lxy_lxy = TH2F("h_corr_lxy_lxy",";lxy correction factor;lxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_lxy_vz = TH2F("h_corr_lxy_vz",";lxy correction factor;vz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)

#dxy against the others
h_corr_vz_pt = TH2F("h_corr_vz_pt",";vz correction factor;pt correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_vz_pz = TH2F("h_corr_vz_pz",";vz correction factor;pz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_vz_dxy = TH2F("h_corr_vz_dxy",";vz correction factor;dxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_vz_dz = TH2F("h_corr_vz_dz",";vz correction factor;dz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_vz_lxy = TH2F("h_corr_vz_lxy",";vz correction factor;lxy correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)
h_corr_vz_vz = TH2F("h_corr_vz_vz",";vz correction factor;vz correction factor;",nbins_corr_par,min_bin_corr_par,max_bin_corr_par,nbins_corr_par,min_bin_corr_par,max_bin_corr_par)

ll_corr_eff = [
[h_corr_pt_pt,h_corr_pt_pz,h_corr_pt_dxy,h_corr_pt_dz,h_corr_pt_lxy,h_corr_pt_vz],
[h_corr_pz_pt,h_corr_pz_pz,h_corr_pz_dxy,h_corr_pz_dz,h_corr_pz_lxy,h_corr_pz_vz],
[h_corr_dxy_pt,h_corr_dxy_pz,h_corr_dxy_dxy,h_corr_dxy_dz,h_corr_dxy_lxy,h_corr_dxy_vz],
[h_corr_dz_pt,h_corr_dz_pz,h_corr_dz_dxy,h_corr_dz_dz,h_corr_dz_lxy,h_corr_dz_vz],
[h_corr_lxy_pt,h_corr_lxy_pz,h_corr_lxy_dxy,h_corr_lxy_dz,h_corr_lxy_lxy,h_corr_lxy_vz],
[h_corr_vz_pt,h_corr_vz_pz,h_corr_vz_dxy,h_corr_vz_dz,h_corr_vz_lxy,h_corr_vz_vz]
]

#for each final state particle check which are the cuts for the fiducial region where it gets killed
h_FiducialRegionCuts_1 = TH1I("h_FiducialRegionCuts_1",";Failing fiducial region cuts;",13,-0.5,12.5)
h_FiducialRegionCuts_2 = TH1I("h_FiducialRegionCuts_2",";Failing fiducial region cuts;",13,-0.5,12.5)
h_FiducialRegionCuts_3 = TH1I("h_FiducialRegionCuts_3",";Failing fiducial region cuts;",13,-0.5,12.5)
h_FiducialRegionCuts_4 = TH1I("h_FiducialRegionCuts_4",";Failing fiducial region cuts;",13,-0.5,12.5)
l_h_FiducialRegionCuts = [h_FiducialRegionCuts_1,h_FiducialRegionCuts_2,h_FiducialRegionCuts_3,h_FiducialRegionCuts_4]


#a list with counters for the reconstructed particles, so there are 7 entries for each of the 7 particles
nAntiS_reconstructable = 0.
nAntiSRecoAlsoOutsideFiducialRegion = 0.
nAntiSRecoInsideFiducialRegion = 0.
#and count separately the reconstruction efficiency of the Ks, antiLambda and antiS if their respective daughters got reconstructed.
nKsRECOIfBothDaughtersReco = 0.
nKsTOTALIfBothDaughtersReco = 0.
nAntiLambdaRECOIfBothDaughtersReco = 0.
nAntiLambdaTOTALIfBothDaughtersReco = 0.
nAntiSRECOIfBothDaughtersReco = 0.
nAntiSTOTALIfBothDaughtersReco = 0.

for iFile, fIn in enumerate(inFiles,start = 1):
	print "Starting with inputFile: ", str(iFile) ,"/",str(len(inFiles)), ':', fIn.GetName()
	tree = fIn.Get('FlatTreeProducerTracking/FlatTreeTpsAntiS') 
	for i in range(0,tree.GetEntries()):
		if(i>maxEvents):
			break
		tree.GetEntry(i)
		printProgress(i)

		weightFactor = tree._tpsAntiS_event_weighting_factor[0]*tree._tpsAntiS_event_weighting_factorPU[0]

		boolNGrandDaughtersWithTrackerHitsLargerThan6 = False
#		if(tree._tpsAntiS_numberOfTrackerHits[3] >= 7 and tree._tpsAntiS_numberOfTrackerHits[4] >= 7 and tree._tpsAntiS_numberOfTrackerHits[5] >= 7 and tree._tpsAntiS_numberOfTrackerHits[6] >= 7): 
#			boolNGrandDaughtersWithTrackerHitsLargerThan6 = True


		NGrandDaughtersWithTrackerHitsLargerThan6 = 0
		for j in range(0,len(tree._tpsAntiS_type)):
			#now look at _tpAntiS_type from 3 to 6, this are the granddaughters:
			if(tree._tpsAntiS_type[j] >= 3 and tree._tpsAntiS_type[j] <= 6 ):
				if(tree._tpsAntiS_numberOfTrackerHits[j] >= 7):
					NGrandDaughtersWithTrackerHitsLargerThan6 += 1

		boolNGrandDaughtersWithTrackerHitsLargerThan6 =  NGrandDaughtersWithTrackerHitsLargerThan6==4

		antiSReconstructed = False
		KsReconstructed = False
		antiLambdaReconstructed = False
		if(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[1] < config_dict["GENRECO_matcher_Ks_deltaL"] and tree._tpsAntiS_bestDeltaRWithRECO[1] < config_dict["GENRECO_matcher_Ks_deltaR"]):
                        KsReconstructed = True
                if(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[2] < config_dict["GENRECO_matcher_AntiL_deltaL"] and tree._tpsAntiS_bestDeltaRWithRECO[2] < config_dict["GENRECO_matcher_AntiL_deltaR"]):
                        antiLambdaReconstructed = True
                if(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[0] < config_dict["GENRECO_matcher_AntiS_deltaL"] and tree._tpsAntiS_bestDeltaRWithRECO[0] < config_dict["GENRECO_matcher_AntiS_deltaR"] and tree._tpsAntiS_reconstructed[3] == 1 and tree._tpsAntiS_reconstructed[4] == 1 and tree._tpsAntiS_reconstructed[5] == 1 and tree._tpsAntiS_reconstructed[6] == 1 and KsReconstructed and antiLambdaReconstructed):
                        antiSReconstructed = True

		#only evaluate the below for reconstructable antiS
		if(not boolNGrandDaughtersWithTrackerHitsLargerThan6): continue

		nAntiS_reconstructable+=weightFactor


		if(antiSReconstructed): nAntiSRecoAlsoOutsideFiducialRegion+=weightFactor	


		#define a fiducial region based on the kinematics of the final state particles. Don't look at antiS which fall outside this fiducial region because outside this region systematic uncertainties become too large
		insideSystUncFiducialRegion = True
		for i_par in [3,4,5,6]:
			pt = tree._tpsAntiS_pt[i_par]
			pz = tree._tpsAntiS_pz[i_par]	
			dxy = tree._tpsAntiS_dxy_beamspot[i_par]	
			dz = tree._tpsAntiS_dz_beamspot[i_par]	
			lxy = tree._tpsAntiS_Lxy_beamspot[i_par]	
			vz = tree._tpsAntiS_vz_beamspot[i_par]	
			#if(pt > 10 or pz > 22 or lxy > 60 or abs(vz) > 150  or dxy > 12 or abs(dz) > 35 ):
			if(pt < FidReg_minPt or pt > FidReg_maxPt or pz < FidReg_minPz or pz > FidReg_maxPz or lxy > FidReg_maxlxy or vz < FidReg_minvz or vz > FidReg_maxvz or dxy < FidReg_mindxy  or dxy > FidReg_maxdxy or dz < FidReg_mindz or dz > FidReg_maxdz ):
				insideSystUncFiducialRegion = False
			if(antiSReconstructed):
				if(pt < FidReg_minPt):
					l_h_FiducialRegionCuts[i_par-3].Fill(0)
				if(pt > FidReg_maxPt):
					l_h_FiducialRegionCuts[i_par-3].Fill(1)
				if(pz < FidReg_minPz):
					l_h_FiducialRegionCuts[i_par-3].Fill(3)
				if(pz > FidReg_maxPz):
					l_h_FiducialRegionCuts[i_par-3].Fill(4)
				if(lxy > FidReg_maxlxy):
					l_h_FiducialRegionCuts[i_par-3].Fill(5)
				if(vz < FidReg_minvz):
					l_h_FiducialRegionCuts[i_par-3].Fill(6)
				if(vz > FidReg_maxvz):
					l_h_FiducialRegionCuts[i_par-3].Fill(7)
				if(dxy < FidReg_mindxy):
					l_h_FiducialRegionCuts[i_par-3].Fill(8)
				if(dxy > FidReg_maxdxy):
					l_h_FiducialRegionCuts[i_par-3].Fill(9)
				if(dz < FidReg_mindz):
					l_h_FiducialRegionCuts[i_par-3].Fill(10)
				if(dz > FidReg_maxdz):
					l_h_FiducialRegionCuts[i_par-3].Fill(11)

		#only proceed with antiS events which have final state particles in the fiducial region
		if(not insideSystUncFiducialRegion): continue

		if(antiSReconstructed): nAntiSRecoInsideFiducialRegion += weightFactor		

		#count the reconstructed particles with the requirement that there daughters got reconstructed
		if(tree._tpsAntiS_reconstructed[3] == 1 and tree._tpsAntiS_reconstructed[4] == 1):
			if(tree._tpsAntiS_reconstructed[1] == 1):
				nKsRECOIfBothDaughtersReco += weightFactor
			nKsTOTALIfBothDaughtersReco += weightFactor
		if(tree._tpsAntiS_reconstructed[5] == 1 and tree._tpsAntiS_reconstructed[6] == 1):
			if(tree._tpsAntiS_reconstructed[2] == 1):
				nAntiLambdaRECOIfBothDaughtersReco += weightFactor
			nAntiLambdaTOTALIfBothDaughtersReco += weightFactor
		if(tree._tpsAntiS_reconstructed[1] == 1 and tree._tpsAntiS_reconstructed[2] == 1):
			if(antiSReconstructed):
				nAntiSRECOIfBothDaughtersReco += weightFactor
			nAntiSTOTALIfBothDaughtersReco += weightFactor

		#now for events where the antiS got reconstructed get the correction factor for each of the kinematic variables for each of the granddaughters:
		if(antiSReconstructed):
			#for each granddaughter:
			if(verbose):print "#################################################################"
			corr_factor_this_antiS = 1
			error_corr_factor_this_antiS = 0
			for i_par in [3,4,5,6]:
				best_matching_variables = []
				corrections_factors = []
				corrections_factors_error = []
				

				pt = tree._tpsAntiS_pt[i_par]
				best_matching_variables.append(corr_factors_pt.iloc[(corr_factors_pt['kinVariable']-pt).abs().argsort()[:1]]["kinVariable"].values[0])
				corrections_factors.append(corr_factors_pt.iloc[(corr_factors_pt['kinVariable']-pt).abs().argsort()[:1]]["DataOverMC"].values[0])
				corrections_factors_error.append(corr_factors_pt.iloc[(corr_factors_pt['kinVariable']-pt).abs().argsort()[:1]]["DataOverMCError"].values[0])

				pz = tree._tpsAntiS_pz[i_par]
				best_matching_variables.append(corr_factors_pz.iloc[(corr_factors_pz['kinVariable']-pz).abs().argsort()[:1]]["kinVariable"].values[0])
				corrections_factors.append(corr_factors_pz.iloc[(corr_factors_pz['kinVariable']-pz).abs().argsort()[:1]]["DataOverMC"].values[0])
				corrections_factors_error.append(corr_factors_pz.iloc[(corr_factors_pz['kinVariable']-pz).abs().argsort()[:1]]["DataOverMCError"].values[0])

				dxy = tree._tpsAntiS_dxy_beamspot[i_par]	
				best_matching_variables.append(corr_factors_dxy.iloc[(corr_factors_dxy['kinVariable']-dxy).abs().argsort()[:1]]["kinVariable"].values[0])
				corrections_factors.append(corr_factors_dxy.iloc[(corr_factors_dxy['kinVariable']-dxy).abs().argsort()[:1]]["DataOverMC"].values[0])
				corrections_factors_error.append(corr_factors_dxy.iloc[(corr_factors_dxy['kinVariable']-dxy).abs().argsort()[:1]]["DataOverMCError"].values[0])

				dz = tree._tpsAntiS_dz_beamspot[i_par]	
				best_matching_variables.append(corr_factors_dz.iloc[(corr_factors_dz['kinVariable']-dz).abs().argsort()[:1]]["kinVariable"].values[0])
				corrections_factors.append(corr_factors_dz.iloc[(corr_factors_dz['kinVariable']-dz).abs().argsort()[:1]]["DataOverMC"].values[0])
				corrections_factors_error.append(corr_factors_dz.iloc[(corr_factors_dz['kinVariable']-dz).abs().argsort()[:1]]["DataOverMCError"].values[0])

				lxy = tree._tpsAntiS_Lxy_beamspot[i_par]	
				best_matching_variables.append(corr_factors_lxy.iloc[(corr_factors_lxy['kinVariable']-lxy).abs().argsort()[:1]]["kinVariable"].values[0])
				corrections_factors.append(corr_factors_lxy.iloc[(corr_factors_lxy['kinVariable']-lxy).abs().argsort()[:1]]["DataOverMC"].values[0])
				corrections_factors_error.append(corr_factors_lxy.iloc[(corr_factors_lxy['kinVariable']-lxy).abs().argsort()[:1]]["DataOverMCError"].values[0])

				vz = tree._tpsAntiS_vz_beamspot[i_par]	
				best_matching_variables.append(corr_factors_vz.iloc[(corr_factors_vz['kinVariable']-vz).abs().argsort()[:1]]["kinVariable"].values[0])
				corrections_factors.append(corr_factors_vz.iloc[(corr_factors_vz['kinVariable']-vz).abs().argsort()[:1]]["DataOverMC"].values[0])
				corrections_factors_error.append(corr_factors_vz.iloc[(corr_factors_vz['kinVariable']-vz).abs().argsort()[:1]]["DataOverMCError"].values[0])
		
				#fill the correlation plots of the correction parameters:
				for i in range(0,len(ll_corr_eff)):
					for j in range(0,len(ll_corr_eff[0])):
						ll_corr_eff[i][j].Fill(corrections_factors[i],corrections_factors[j],weightFactor)
				

				corr_factor_this_particle = 1	
				error_corr_factor_this_particle = 0
				for i_cor in range(0,len(corrections_factors)):
					corr_factor_this_particle = corr_factor_this_particle*corrections_factors[i_cor]
					error_corr_factor_this_particle = error_corr_factor_this_particle +  np.power(corrections_factors_error[i_cor],2)/np.power(corrections_factors[i_cor],2)
				corr_factor_this_antiS = corr_factor_this_antiS * corr_factor_this_particle
				error_corr_factor_this_antiS = error_corr_factor_this_antiS + error_corr_factor_this_particle 
				if(verbose):
					print "particle has pt: ", pt, " so the best found value is: ", best_matching_variables[0]," and the correction factor is: ", corrections_factors[0], " and the error is: ", corrections_factors_error[0]
					print "particle has pz: ", pz, " so the best found value is: ", best_matching_variables[1]," and the correction factor is: ", corrections_factors[1], " and the error is: ", corrections_factors_error[1]
					print "particle has dxy: ", dxy, " so the best found value is: ", best_matching_variables[2]," and the correction factor is: ", corrections_factors[2], " and the error is: ", corrections_factors_error[2]
					print "particle has dz: ", dz, " so the best found value is: ", best_matching_variables[3]," and the correction factor is: ", corrections_factors[3], " and the error is: ", corrections_factors_error[3]
					print "particle has lxy: ", lxy, " so the best found value is: ", best_matching_variables[4]," and the correction factor is: ", corrections_factors[4], " and the error is: ", corrections_factors_error[4]
					print "particle has vz: ", vz, " so the best found value is: ", best_matching_variables[5]," and the correction factor is: ", corrections_factors[5], " and the error is: ", corrections_factors_error[5]
					print "The overall correction factor for this particle: ", corr_factor_this_particle
					print "-----------------------------------------------"
			corr_factor_this_antiS = corr_factor_this_antiS
			error_corr_factor_this_antiS = corr_factor_this_antiS*np.sqrt(error_corr_factor_this_antiS)
			print "The overall correction factor for this antiS: ", corr_factor_this_antiS, " with an error of: ", error_corr_factor_this_antiS
			h_corr_factor_this_antiS.Fill(corr_factor_this_antiS,weightFactor)
			h_corr_factor_pt.Fill(corrections_factors[0])
			h_corr_factor_pz.Fill(corrections_factors[1])
			h_corr_factor_dxy.Fill(corrections_factors[2])
			h_corr_factor_dz.Fill(corrections_factors[3])
			h_corr_factor_lxy.Fill(corrections_factors[4])
			h_corr_factor_vz.Fill(corrections_factors[5])

			
		#now select tracks which have the antiS properly reconstructed and check for these events how the kinematics look like of the final state particles.
		#if(antiSReconstructed):
h_corr_factor_this_antiS.Write()

c_name = "c_"+h_corr_factor_this_antiS.GetName()
c = TCanvas(c_name,"")
CMS_lumi.CMS_lumi(c, 0, 11)
h_corr_factor_this_antiS.SetLineColor(colours[0])
h_corr_factor_this_antiS.SetLineWidth(2)
h_corr_factor_this_antiS.SetMarkerStyle(22+0)
h_corr_factor_this_antiS.SetMarkerColor(colours[0])
h_corr_factor_this_antiS.Fit("gaus")
h_corr_factor_this_antiS.Draw("")
c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
c.Write()


h_corr_factor_pt.Write()
h_corr_factor_pz.Write()
h_corr_factor_dxy.Write()
h_corr_factor_dz.Write()
h_corr_factor_lxy.Write()
h_corr_factor_vz.Write()

l_dir = ["pt","pz","dxy","dz","lxy","vz"]
#list of list for the correlation factors between the correction parameters
ll_correlation_correction_factors = []
i = 0
for l in ll_corr_eff:
	corr_parameters_dir = fOut.mkdir("corr_parameters_"+l_dir[i])
	corr_parameters_dir.cd()
	c = TCanvas("corr_parameters_"+l_dir[i],"corr_parameters_"+l_dir[i],800,600)
	gStyle.SetOptStat(0)
	c.Divide(2,3)
	j = 1
	l_correlation_correction_factors = []
	for h in l:
		c.cd(j)
		h.Draw("colz")
		h.GetXaxis().SetTitleSize(0.05)
		h.GetYaxis().SetTitleSize(0.05)
		print "Correlation factors for ", h.GetName(), " is: ", h.GetCorrelationFactor()
		h.Write()	
		j+=1
		l_correlation_correction_factors.append(h.GetCorrelationFactor())
	ll_correlation_correction_factors.append(l_correlation_correction_factors)
	c.Write()
	c.SaveAs(plots_output_dir+"corr_parameters_"+l_dir[i]+".pdf")
	i+=1
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print ll_correlation_correction_factors
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

fOut.cd()
text_labels = ["p_{t} correction factor","p_{z} correction factor", "d_{0} correction factor", "d_{z} correction factor", "l_{0} correction factor", "v_{z} correction factor"]
h2_correlation_factors_corr_parameters = TH2D("h2_correlation_factors_corr_parameters",";;;correlation factor",6,-0.5,5.5,6,-0.5,5.5)

c_name = "c_"+h2_correlation_factors_corr_parameters.GetName()
c = TCanvas(c_name,"")

for i in range(0,len(ll_correlation_correction_factors)):
	for j in range(0,len(ll_correlation_correction_factors[0])):
		h2_correlation_factors_corr_parameters.SetBinContent(i+1,j+1,ll_correlation_correction_factors[i][j])

xax = h2_correlation_factors_corr_parameters.GetXaxis()
yax = h2_correlation_factors_corr_parameters.GetYaxis()
for i in range(0,len(text_labels)):
	xax.SetBinLabel(i+1,text_labels[i])
	yax.SetBinLabel(i+1,text_labels[i])

xax.LabelsOption("v")
#CMS_lumi.CMS_lumi(c, 0, 11)
c.SetBottomMargin(3)
c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
h2_correlation_factors_corr_parameters.Draw("colztext")
c.Write()
h2_correlation_factors_corr_parameters.Write()		

FiducialCuts_dir = fOut.mkdir("FiducialCuts")
FiducialCuts_dir.cd()
for h in l_h_FiducialRegionCuts:
	h.Write()



fOut.Close()

#print some conclusions:

#summary of reconstruction efficiencies:
particles = ["antiS","Ks","AntiLambda","Ks-piplus","Ks-pineg","AntiLambda-pion","AntiLambda-antiproton"]
print "###################################################################################################"
print "########################################RECO EFF###################################################"
print "###################################################################################################"

print "AntiS reconstruction efficiency (weighted for the path length through the beampipe and PV) for antiS with final state particles which are reconstructable: ", nAntiSRecoAlsoOutsideFiducialRegion,"/",nAntiS_reconstructable, " = ", nAntiSRecoAlsoOutsideFiducialRegion/nAntiS_reconstructable, "       --> so 1 out of: " , nAntiS_reconstructable/nAntiSRecoAlsoOutsideFiducialRegion

print "AntiS reconstruction efficiency (weighted for the path length through the beampipe and PV) for antiS with final state particles which are reconstructable and have to be in the fiducial region: ", nAntiSRecoInsideFiducialRegion,"/",nAntiS_reconstructable, " = ", nAntiSRecoInsideFiducialRegion/nAntiS_reconstructable, "       --> so 1 out of: " , nAntiS_reconstructable/nAntiSRecoInsideFiducialRegion


print "------------------------------------------"
print "Ks if both daughters got reconstructed: ", nKsRECOIfBothDaughtersReco,"/",nKsTOTALIfBothDaughtersReco," = ", nKsRECOIfBothDaughtersReco/nKsTOTALIfBothDaughtersReco, "              --> so 1 out of :", nKsTOTALIfBothDaughtersReco/nKsRECOIfBothDaughtersReco
print "AntiLambda if both daughters got reconstructed: ", nAntiLambdaRECOIfBothDaughtersReco,"/",nAntiLambdaTOTALIfBothDaughtersReco," = ", nAntiLambdaRECOIfBothDaughtersReco/nAntiLambdaTOTALIfBothDaughtersReco, "              --> so 1 out of :", nAntiLambdaTOTALIfBothDaughtersReco/nAntiLambdaRECOIfBothDaughtersReco
print "AntiS if both daughters got reconstructed: ", nAntiSRECOIfBothDaughtersReco,"/",nAntiSTOTALIfBothDaughtersReco," = ", nAntiSRECOIfBothDaughtersReco/nAntiSTOTALIfBothDaughtersReco, "              --> so 1 out of :", nAntiSTOTALIfBothDaughtersReco/nAntiSRECOIfBothDaughtersReco

