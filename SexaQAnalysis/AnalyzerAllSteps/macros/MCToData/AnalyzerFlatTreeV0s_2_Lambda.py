from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, gROOT
import numpy as np
import copy
from array import array

fIns = [
#the standard V0s:
#'/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_standardV0s/combined_FlatTreeV0_DoubleMuonData_RunC_D_E_G_H.root',
#'/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_standardV0s/combined_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root',

#the adapted V0s:
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/combined_FlatTreeV0_DoubleMuonData_Run_C_H.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunC_v1/combined_FlatTreeV0_DoubleMuonData_Run_C.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunD_v2/combined_FlatTreeV0_DoubleMuonData_Run_D.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunE_v3/combined_FlatTreeV0_DoubleMuonData_Run_E.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunH_v6/combined_FlatTreeV0_DoubleMuonData_Run_H.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/Results/FlatTrees_DoubleMuon_adaptedV0s_extendedFlatTree6/RunG_v5/combined_FlatTreeV0_DoubleMuonData_Run_G.root',
'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/Results/FlatTrees_DoubleMuon_adaptedV0s_extendedFlatTree6/combined_FlatTreeV0_DoubleMuonData_Run_G_H.root',
'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/Results/FlatTrees_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_adaptedV0s_extendedFlatTree6/combined_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root'
]

def eff_error(data,MC):
	return data/MC*np.sqrt(1/MC+1/data)	

#checks which one is the largest in absolute terms and then returns the signed value
def maxSigned(x1, x2):
	if(abs(x1) > abs(x2)):
		return x1
	else:
		return x2

#fist find which of the trees has the smallest number of entries. Data or MC?
fIn1 = TFile(fIns[0],'read')
fIn2 = TFile(fIns[1],'read')
treeLambda1 = fIn1.Get('FlatTreeProducerV0s/FlatTreeLambda')
treeLambda2 = fIn2.Get('FlatTreeProducerV0s/FlatTreeLambda')
min_entries = min( treeLambda1.GetEntries(), treeLambda2.GetEntries() )
#min_entries = 100000 
print "the nentries in Data: " , treeLambda1.GetEntries()
print "the nentries in MC: " , treeLambda2.GetEntries()
print "the min_enties: ", min_entries

#get the means of the PV_z distribution for data and MC. Use these values to offset the distributions later
#treePV1 = fIn1.Get('FlatTreeProducerV0s/FlatTreePV')
#treePV1.Draw("_PV0_vz>>hist_PVz_Data")
#hist_PVz_Data = gROOT.FindObject('hist_PVz_Data')
#mean_hist_PVz_Data = hist_PVz_Data.GetMean()
#print "mean PV_z for Data ", mean_hist_PVz_Data
#
#treePV2 = fIn2.Get('FlatTreeProducerV0s/FlatTreePV')
#treePV2.Draw("_PV0_vz>>hist_PVz_MC")
#hist_PVz_MC = gROOT.FindObject('hist_PVz_MC')
#mean_hist_PVz_MC = hist_PVz_MC.GetMean()
#print "mean PV_z for MC ", mean_hist_PVz_MC
#
##get the means of the Lambda_vz distributions for data and MC. Use these values to offset the distributions later. 
#treeLambda1.Draw("_Lambda_vz>>hist_Lambda_vz_Data")
#hist_Lambda_vz_Data = gROOT.FindObject('hist_Lambda_vz_Data')
#mean_hist_Lambda_vz_Data = hist_Lambda_vz_Data.GetMean()
#print "mean Lambda_vz for Data ", mean_hist_Lambda_vz_Data 
#
#treeLambda2.Draw("_Lambda_vz>>hist_Lambda_vz_MC")
#hist_Lambda_vz_MC = gROOT.FindObject('hist_Lambda_vz_MC')
#mean_hist_Lambda_vz_MC = hist_Lambda_vz_MC.GetMean()
#print "mean Lambda_vz for MC ", mean_hist_Lambda_vz_MC 
#
##get the means of the _Lambda_dz_PV distributions for data and MC. Use these values to offset the distributions later. 
#treeLambda1.Draw("_Lambda_dz>>hist_Lambda_dz_PV_data")
#hist_Lambda_dz_PV_data= gROOT.FindObject('hist_Lambda_dz_PV_data')
#mean_hist_Lambda_dz_PV_Data = hist_Lambda_dz_PV_data.GetMean()
#print "mean Lambda_dz_PV for Data ", mean_hist_Lambda_dz_PV_Data
#
#treeLambda2.Draw("_Lambda_vz>>hist_Lambda_dz_PV_MC")
#hist_Lambda_dz_PV_MC = gROOT.FindObject('hist_Lambda_dz_PV_MC')
#mean_hist_Lambda_dz_PV_MC = hist_Lambda_dz_PV_MC.GetMean()
#print "mean Lambda_dz_PV for MC ", mean_hist_Lambda_dz_PV_MC 


#first loop over the trees and look at the PV0 distribution in z both for DATA and MC. This is important as these distributions are different for data and MC, so you have to correct for this when calculating e.g. the vz
#iFile = 0
#h_RECO_PV0_vz_Data_shifted = TH1F('h_RECO_PV0_vz_Data_shifted'," ;PV vz (cm) - mean PV vz; #entries",200,-100,100)
#h_RECO_PV0_vz_MC_shifted = TH1F('h_RECO_PV0_vz_MC_shifted'," ;PV vz - mean PV vz (cm); #entries",200,-100,100)
#for fIn in fIns:
#	fIn = TFile(fIn,'read')
#	treePV  = fIn.Get('FlatTreeProducerV0s/FlatTreePV')
#	for i in range(0,min_entries):
#		treePV.GetEntry(i)
#		if(iFile == 0):
#			h_RECO_PV0_vz_Data_shifted.Fill(treePV._PV0_vz[0]-mean_hist_PVz_Data)	
#		else:	
#			h_RECO_PV0_vz_MC_shifted.Fill(treePV._PV0_vz[0]-mean_hist_PVz_MC)
#	iFile = iFile+1	
#
#mean_vz_data = h_RECO_PV0_vz_Data_shifted.GetMean()
#RMS_vz_data  = h_RECO_PV0_vz_Data_shifted.GetRMS()
#mean_vz_MC   = h_RECO_PV0_vz_MC_shifted.GetMean()
#RMS_vz_MC    = h_RECO_PV0_vz_MC_shifted.GetRMS()
#
##normalize the histogram of PV in data to later use it for reweighting 
#h_RECO_PV0_vz_Data_shifted.Scale(1./h_RECO_PV0_vz_Data_shifted.Integral(), "width")
#h_RECO_PV0_vz_MC_shifted.Scale(1./h_RECO_PV0_vz_MC_shifted.Integral(), "width")


#print "Mean and RMS of PV Data: ", mean_vz_data, " ", RMS_vz_data
#print "Mean and RMS of PV MC: ", mean_vz_MC, " ", RMS_vz_MC

iFile = 0
h_RECO_PV0_vz_Data= TH1F('h_RECO_PV0_vz_Data'," ;PV vz (cm); #entries",20000,-100,100)
h_RECO_PV0_vz_MC= TH1F('h_RECO_PV0_vz_MC'," ;PV vz (cm); #entries",20000,-100,100)
for fIn in fIns:
	fIn = TFile(fIn,'read')
	treePV  = fIn.Get('FlatTreeProducerV0s/FlatTreePV')
	for i in range(0,min_entries):
		treePV.GetEntry(i)
		if(iFile == 0):
			h_RECO_PV0_vz_Data.Fill(treePV._PV0_vz[0])	
		else:	
			h_RECO_PV0_vz_MC.Fill(treePV._PV0_vz[0])
	iFile = iFile+1	


#normalize the histogram of PV in data to later use it for reweighting 
#h_RECO_PV0_vz_Data.Scale(1./h_RECO_PV0_vz_Data.Integral(), "width")
#h_RECO_PV0_vz_MC.Scale(1./h_RECO_PV0_vz_MC.Integral(), "width")



fOut = TFile('Results/output_DataToMC_RunG_H_with_dxy_dz_min_PV_cut_reweighing_on_Lambda_vz_Dz_min_PV.root','RECREATE')
h_dummy =  TH1F('h_dummy'," x; y",20,0,10)
histos_Data = [h_dummy]*9
histos_MC = [h_dummy]*9


nLambdahortPerEventData = []
nLambdahortPerEventMC = []

iFile = 0
for fIn in fIns:

	print "creating plots for file: ", fIn

	fIn = TFile(fIn,'read')

	treeLambda = fIn.Get('FlatTreeProducerV0s/FlatTreeLambda') 
	treeZ = fIn.Get('FlatTreeProducerV0s/FlatTreeZ')
	treePV = fIn.Get('FlatTreeProducerV0s/FlatTreePV')

	nEntriesLambda = treeLambda.GetEntries()
	nEntriesZ = treeZ.GetEntries()
	nEntriesPV = treePV.GetEntries()

	print 'Number of entries in the treeLambda: ', nEntriesLambda
	print 'Number of entries in the treeZ: ', nEntriesZ
	print 'Number of entries in the treePV: ', nEntriesPV

	bins1 = np.arange(-300,-200,20)
	bins2 = np.arange(-200,-100,10)
	bins3 = np.arange(-100,-80,5)
	bins4 = np.arange(-80,-70,2)
	bins5 = np.arange(-70,-64,1.5)
	bins6 = np.arange(-64,64,1)
	bins7 = np.arange(64,70,1.5)
	bins8 = np.arange(70,80,2)
	bins9 = np.arange(80,100,5)
	bins10 = np.arange(100,200,10)
	bins11 = np.arange(200,300,20)
	args = (bins1,bins2,bins3,bins4,bins5,bins6,bins7,bins8,bins9,bins10,bins11)
	bins = np.concatenate(args)
	print bins
	h_RECO_Lambda_vz = TH1F('h_RECO_Lambda_vz'," ;Lambda v_{z} (cm); #entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_vz_NOTWeighted = TH1F('h_RECO_Lambda_vz_NOTWeighted'," ;Lambda v_{z} (cm); #entries",len(bins)-1,array('d',bins))
	bins1 = np.arange(0,50,1)
	bins2 = np.arange(50,56,1.5)
	bins3 = np.arange(56,65,1.8)
	bins4 = np.arange(65,70,2.5)
	bins5 = np.arange(70,90,5)
	args = (bins1,bins2,bins3,bins4,bins5)
	bins = np.concatenate(args)
	h_RECO_Lambda_lxy = TH1F('h_RECO_Lambda_lxy'," ;Lambda l_{0} (cm);#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_lxy_NOTWeighted = TH1F('h_RECO_Lambda_lxy_NOTWeighted'," ;Lambda l_{0} (cm);#entries",len(bins)-1,array('d',bins))
	bins1 = np.arange(0,8,0.2)
	bins2 = np.arange(8,10.1,0.3)
	args = (bins1,bins2)
	bins = np.concatenate(args)
	h_RECO_Lambda_pt = TH1F('h_RECO_Lambda_pt',"; Lambda p_{t} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_pt_NOTWeighted = TH1F('h_RECO_Lambda_pt_NOTWeighted',"; Lambda p_{t} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_pt_tracks1 = TH1F('h_RECO_Lambda_pt_tracks1',"; Lambda daughter 1 p_{t} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_pt_tracks2 = TH1F('h_RECO_Lambda_pt_tracks2',"; Lambda daughter 2 p_{t} (GeV);#entries",len(bins)-1,array('d',bins))
	bins1 = np.arange(0,15,1)
	bins2 = np.arange(15,21,1.5)
	bins3 = np.arange(21,25,2)
	bins4 = np.arange(25,40,5)
	args = (bins1,bins2,bins3,bins4)
	bins = np.concatenate(args)
	h_RECO_Lambda_pz= TH1F('h_RECO_Lambda_pz',"; Lambda p_{z} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_pz_NOTWeighted= TH1F('h_RECO_Lambda_pz_NOTWeighted',"; Lambda p_{z} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_pz_tracks1 = TH1F('h_RECO_Lambda_pz_tracks1',"; Lambda daughter 1 p_{z} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_pz_tracks2 = TH1F('h_RECO_Lambda_pz_tracks2',"; Lambda daughter 2 p_{z} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_dxy_PV= TH1F('h_RECO_Lambda_dxy_PV',"; Lambda d_{0} (cm); #entries",100,-10,10)
	h_RECO_Lambda_dxy_PV_NOTWeighted= TH1F('h_RECO_Lambda_dxy_PV_NOTWeighted',"; Lambda d_{0} (cm); #entries",100,-10,10)
	bins = [ -30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30 ]
	bins1 = np.arange(-30,-20,2)
	bins2 = np.arange(-20,-14,1.5)
	bins3 = np.arange(-14,14,1)
	bins4 = np.arange(14,20,1.5)
	bins5 = np.arange(20,30,2)
	args = (bins1,bins2,bins3,bins4,bins5)
	bins = np.concatenate(args)
	h_RECO_Lambda_dz= TH1F('h_RECO_Lambda_dz',";  Lambda d_{z}(0,0,0) (cm); #entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_dz_NOTWeighted = TH1F('h_RECO_Lambda_dz_NOTWeighted',";  Lambda d_{z}(0,0,0) (cm); #entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_eta= TH1F('h_RECO_Lambda_eta',"; Lambda #eta (cm);#entries",100,-4,4)
	h_RECO_Lambda_eta_NOTWeighted= TH1F('h_RECO_Lambda_eta_NOTWeighted',"; Lambda #eta (cm);#entries",100,-4,4)
	h_RECO_Lambda_phi= TH1F('h_RECO_Lambda_phi',"; Lambda #phi (cm);#entries",100,-4,4)
	h_RECO_Lambda_phi_NOTWeighted= TH1F('h_RECO_Lambda_phi_NOTWeighted',"; Lambda #phi (cm);#entries",100,-4,4)
	h_RECO_Lambda_Track1Track2_openingsAngle= TH1F('h_RECO_Lambda_Track1Track2_openingsAngle',"; openingsangle between track1 and track2 of the Lambda (rad) ;#entries",100,35,3.5)
	h_RECO_Lambda_Track1Track2_openingsAngle_NOTWeighted = TH1F('h_RECO_Lambda_Track1Track2_openingsAngle_NOTWeighted',"; openingsangle between track1 and track2 of the Lambda (rad) ;#entries",100,35,3.5)
	h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p2 = TH1F('h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p2',"; openingsangle between track1 and track2 of the Lambda (rad) ;#entries",100,35,3.5)
	h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p5= TH1F('h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p5',"; openingsangle between track1 and track2 of the Lambda (rad) ;#entries",100,35,3.5)
	bins = [ 0, .5, 1, 1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,8,9,10,12,14,16 ]
	h_RECO_Lambda_Track1Track2_max_dxy_beamspot= TH1F('h_RECO_Lambda_Track1Track2_max_dxy_beamspot',"; maxSigned[d_{0}(beamspot) track1, d_{0}(beamspot) track2]   Lambda (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_Track1Track2_max_dxy_beamspot_ptCut1p2= TH1F('h_RECO_Lambda_Track1Track2_max_dxy_beamspot_ptCut1p2',"; maxSigned[d_{0}(beamspot) track1, d_{0}(beamspot) track2]   Lambda (cm) ;#entries",len(bins)-1,array('d',bins))
	bins = [ -40,-35,-30,-28,-26,-24,-22,-20,-18,-16,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,30,35,40 ]
	h_RECO_Lambda_Track1Track2_max_dz_min_PV = TH1F('h_RECO_Lambda_Track1Track2_max_dz_min_PV',"; maxSigned[d_{z}(best PV Lambda) track1, d_{z}(best PV Lambda) track2]   Lambda (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Lambda_Track1Track2_max_dz_beamspot = TH1F('h_RECO_Lambda_Track1Track2_max_dz_beamspot',"; maxSigned[d_{z}(beamspot) track1, d_{z}(beamspot) track2]   Lambda (cm) ;#entries",len(bins)-1,array('d',bins))

	histos= [h_RECO_Lambda_vz, h_RECO_Lambda_vz_NOTWeighted, h_RECO_Lambda_lxy, h_RECO_Lambda_lxy_NOTWeighted, h_RECO_Lambda_pt, h_RECO_Lambda_pt_NOTWeighted, h_RECO_Lambda_pt_tracks1, h_RECO_Lambda_pt_tracks2, h_RECO_Lambda_pz, h_RECO_Lambda_pz_NOTWeighted,h_RECO_Lambda_pz_tracks1, h_RECO_Lambda_pz_tracks2, h_RECO_Lambda_dxy_PV, h_RECO_Lambda_dxy_PV_NOTWeighted, h_RECO_Lambda_dz, h_RECO_Lambda_dz_NOTWeighted, h_RECO_Lambda_eta, h_RECO_Lambda_eta_NOTWeighted, h_RECO_Lambda_phi, h_RECO_Lambda_phi_NOTWeighted, h_RECO_Lambda_Track1Track2_openingsAngle, h_RECO_Lambda_Track1Track2_openingsAngle_NOTWeighted,h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p2,h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p5,h_RECO_Lambda_Track1Track2_max_dxy_beamspot,h_RECO_Lambda_Track1Track2_max_dxy_beamspot_ptCut1p2,h_RECO_Lambda_Track1Track2_max_dz_min_PV,h_RECO_Lambda_Track1Track2_max_dz_beamspot] 
	#check the reweighting
	h_RECO_PV0_vz_MC_reweighted_to_data = TH1F('h_RECO_PV0_vz_MC_reweighted_to_data'," ;PV v_{z} (cm); #entries",200,-100,100)
	h_RECO_PV0_vz_MC_NOTreweighted_to_data = TH1F('h_RECO_PV0_vz_MC_NOTreweighted_to_data',"; PV v_{z} (cm); #entries",200,-100,100)


	for i in range(0,min_entries):
		if(i%100000 == 0):
			print 'reached event ', i
		
		treeZ.GetEntry(i)
		treeLambda.GetEntry(i)
		treePV.GetEntry(i)
		
		nLambdahortThisEvent = 0
	

		#only select a certain hard scale
#		if(treeZ._Z_ptMuMu[0] < 5):
#		       continue
		
		#now loop over all the Lambda in this event:
		for j in range(0, len(treeLambda._Lambda_mass)):
			
			#require only Lambda coming from the first PV
			#if(treeLambda._Lambda_mass[j] > 0.48 and treeLambda._Lambda_mass[j] < 0.52 and  abs(treeLambda._Lambda_eta[j]) < 2 and treeLambda._Lambda_dxy_beamspot[j] < 0.1 and treeLambda._Lambda_dxy_beamspot[j] > 0. and abs(treeLambda._Lambda_dz_PV0[j]) < 0.2) :
			#require only Lambda which point to a 'best PV'
			if(treeLambda._Lambda_mass[j] > 1.115-0.01 and treeLambda._Lambda_mass[j] < 1.115+0.01 and  abs(treeLambda._Lambda_eta[j]) < 2 and treeLambda._Lambda_dxy_beamspot[j] < 0.1 and treeLambda._Lambda_dxy_beamspot[j] > 0. and abs(treeLambda._Lambda_dz_min_PV[j]) < 0.2) :
			#require only Lambda which point away in dz
			#if(treeLambda._Lambda_mass[j] > 0.48 and treeLambda._Lambda_mass[j] < 0.52 and  abs(treeLambda._Lambda_eta[j]) < 2 and treeLambda._Lambda_dxy_beamspot[j] < 0.1 and treeLambda._Lambda_dxy_beamspot[j] > 0. and abs(treeLambda._Lambda_dz_min_PV[j]) > 0.2) :
			#require 'PU Lambda' by asking for a Lambda with large _Lambda_dz_PV0, but small _Lambda_dz_min_PV
			#if(treeLambda._Lambda_mass[j] > 0.48 and treeLambda._Lambda_mass[j] < 0.52 and  abs(treeLambda._Lambda_eta[j]) < 2 and treeLambda._Lambda_dxy_beamspot[j] < 0.1 and treeLambda._Lambda_dxy_beamspot[j] > 0. and abs(treeLambda._Lambda_dz_PV0[j]) > 0.2 and abs(treeLambda._Lambda_dz_min_PV[j]) < 0.2) :

			#if(treeLambda._Lambda_mass[j] > 0.48 and treeLambda._Lambda_mass[j] < 0.52 and  abs(treeLambda._Lambda_eta[j]) < 2 and treeLambda._Lambda_dxy_beamspot[j] < 0.1 and treeLambda._Lambda_dxy_beamspot[j] > 0.) :
				nLambdahortThisEvent+=1
				x_axis_h_RECO_PV0_vz_Data = h_RECO_PV0_vz_Data.GetXaxis()
				#bin_prob_h_RECO_PV0_vz_Data = x_axis_h_RECO_PV0_vz_Data.FindBin(treePV._PV0_vz[0]) #the bin which contains the prob of finding treePV._PV0_vz[0] in the data
				bin_prob_h_RECO_PV0_vz_Data = x_axis_h_RECO_PV0_vz_Data.FindBin(treeLambda._Lambda_vz_dz_min_PV[j]) 
				prob_h_RECO_PV0_vz_Data = h_RECO_PV0_vz_Data.GetBinContent(bin_prob_h_RECO_PV0_vz_Data) #the actual probability

				x_axis_h_RECO_PV0_vz_MC= h_RECO_PV0_vz_MC.GetXaxis()
				#bin_prob_h_RECO_PV0_vz_MC= x_axis_h_RECO_PV0_vz_MC.FindBin(treePV._PV0_vz[0]) #the bin which contains the prob of finding treePV._PV0_vz[0] in the data
				bin_prob_h_RECO_PV0_vz_MC= x_axis_h_RECO_PV0_vz_MC.FindBin(treeLambda._Lambda_vz_dz_min_PV[j]) 
				prob_h_RECO_PV0_vz_MC= h_RECO_PV0_vz_MC.GetBinContent(bin_prob_h_RECO_PV0_vz_MC) #the actual probability

				PV_Z_reweighting_factor = 1
				if(prob_h_RECO_PV0_vz_MC > 0 and iFile == 1): #only reweigh for MC
					PV_Z_reweighting_factor = prob_h_RECO_PV0_vz_Data/prob_h_RECO_PV0_vz_MC

				#print "------------------------------------------------------------------------------------"
				#print "PV_z: ", treePV._PV0_vz[0]
				#print "PV_z shifted Data, prob: ", treePV._PV0_vz[0]-mean_hist_PVz_Data, " ", prob_h_RECO_PV0_vz_Data_shifted  
				#print "PV_z shifted MC, prob: ", treePV._PV0_vz[0]-mean_hist_PVz_MC, " ", prob_h_RECO_PV0_vz_MC_shifted  
				#print "reweighting factor: ",PV_Z_reweighting_factor
				if(iFile == 0):
					histos[0].Fill(treeLambda._Lambda_vz[j], PV_Z_reweighting_factor)
				if(iFile == 1):
					histos[0].Fill(treeLambda._Lambda_vz[j], PV_Z_reweighting_factor)
				histos[1].Fill(treeLambda._Lambda_vz[j])

				histos[2].Fill(treeLambda._Lambda_Lxy[j], PV_Z_reweighting_factor)
				histos[3].Fill(treeLambda._Lambda_Lxy[j])

				histos[4].Fill(treeLambda._Lambda_pt[j], PV_Z_reweighting_factor)
				histos[5].Fill(treeLambda._Lambda_pt[j])
				histos[6].Fill(treeLambda._Lambda_daughterTrack1_pt[j]) #daughters of the Lambda are identical objects so can put them in same plot
				histos[7].Fill(treeLambda._Lambda_daughterTrack2_pt[j])

				histos[8].Fill(treeLambda._Lambda_pz[j], PV_Z_reweighting_factor)
				histos[9].Fill(treeLambda._Lambda_pz[j])
				histos[10].Fill(treeLambda._Lambda_daughterTrack1_pz[j])
				histos[11].Fill(treeLambda._Lambda_daughterTrack2_pz[j])

				histos[12].Fill(treeLambda._Lambda_dxy_beamspot[j], PV_Z_reweighting_factor)
				histos[13].Fill(treeLambda._Lambda_dxy_beamspot[j])

				histos[14].Fill(treeLambda._Lambda_dz_000[j], PV_Z_reweighting_factor)
				histos[15].Fill(treeLambda._Lambda_dz_000[j])

				histos[16].Fill(treeLambda._Lambda_eta[j], PV_Z_reweighting_factor)
				histos[17].Fill(treeLambda._Lambda_eta[j])
		
				histos[18].Fill(treeLambda._Lambda_phi[j], PV_Z_reweighting_factor)
				histos[19].Fill(treeLambda._Lambda_phi[j])

				histos[20].Fill(treeLambda._Lambda_Track1Track2_openingsAngle[j], PV_Z_reweighting_factor)
				histos[21].Fill(treeLambda._Lambda_Track1Track2_openingsAngle[j])

				if(treeLambda._Lambda_pt[j] > 1.2):
					histos[22].Fill(treeLambda._Lambda_Track1Track2_openingsAngle[j], PV_Z_reweighting_factor)
				if(treeLambda._Lambda_pt[j] > 1.5):
					histos[23].Fill(treeLambda._Lambda_Track1Track2_openingsAngle[j], PV_Z_reweighting_factor)

				histos[24].Fill(maxSigned(treeLambda._Lambda_daughterTrack1_dxy_beamspot[j],treeLambda._Lambda_daughterTrack2_dxy_beamspot[j]), PV_Z_reweighting_factor)
				if(treeLambda._Lambda_pt[j] > 1.2):
					histos[25].Fill(maxSigned(treeLambda._Lambda_daughterTrack1_dxy_beamspot[j],treeLambda._Lambda_daughterTrack2_dxy_beamspot[j]), PV_Z_reweighting_factor)
				histos[26].Fill(maxSigned(treeLambda._Lambda_daughterTrack1_dz_min_PV[j],treeLambda._Lambda_daughterTrack2_dz_min_PV[j]), PV_Z_reweighting_factor)
				histos[27].Fill(maxSigned(treeLambda._Lambda_daughterTrack1_dz_beamspot[j],treeLambda._Lambda_daughterTrack2_dz_beamspot[j]), PV_Z_reweighting_factor)

				if(iFile == 1): #for MC
					h_RECO_PV0_vz_MC_reweighted_to_data.Fill(treePV._PV0_vz[0], PV_Z_reweighting_factor)
					h_RECO_PV0_vz_MC_NOTreweighted_to_data.Fill(treePV._PV0_vz[0], 1)


				#for V0 pointing at the PV0
				#if(abs(treeLambda._Lambda_dxy[j]) < 0.015  and abs(treeLambda._Lambda_dz_PV[j]) < 0.01 ):

		if(iFile == 0):#for Data
			nLambdahortPerEventData.append(nLambdahortThisEvent)
		if(iFile == 1): #for MC
			nLambdahortPerEventMC.append(nLambdahortThisEvent)
			
			
	for h in histos:
		h.SetDirectory(0)
	if(iFile == 0):
		histos_Data = histos
	else:
		histos_MC = histos

	iFile += 1


fOut.mkdir('PV_Histos')
fOut.cd('PV_Histos')
h_RECO_PV0_vz_Data.Rebin(100)
h_RECO_PV0_vz_Data.Scale(1./h_RECO_PV0_vz_Data.Integral(), "width")
h_RECO_PV0_vz_Data.Write()

h_RECO_PV0_vz_MC.Rebin(100)
h_RECO_PV0_vz_MC.Scale(1./h_RECO_PV0_vz_MC.Integral(), "width")
h_RECO_PV0_vz_MC.Write()

h_RECO_PV0_vz_MC_reweighted_to_data.Scale(1./h_RECO_PV0_vz_MC_reweighted_to_data.Integral(), "width")
h_RECO_PV0_vz_MC_reweighted_to_data.Write()
h_RECO_PV0_vz_MC_NOTreweighted_to_data.Scale(1./h_RECO_PV0_vz_MC_NOTreweighted_to_data.Integral(), "width")
h_RECO_PV0_vz_MC_NOTreweighted_to_data.Write()

fOut.mkdir('Data_Histos')
fOut.cd('Data_Histos')
for h in histos_Data:
#	h.Scale(1./h.Integral(),"Width")
	h.Write()

fOut.mkdir('MC_Histos')
fOut.cd('MC_Histos')
for h in histos_MC:
#	h.Scale(1./h.Integral(),"Width")
	h.Write()

fOut.mkdir('Data_TO_MC_Histos')
fOut.cd('Data_TO_MC_Histos')
i = 0
for i in range(0,len(histos_Data)):
	#histos_Data[i].Divide(histos_MC[i])
	#histos_Data[i].SetName(histos_Data[i].GetName()+'DataToMC')
	#histos_Data[i].Write()

	histo_Divide = histos_Data[i].Clone()
	for j in range(1,histo_Divide.GetNbinsX()):
		data = histos_Data[i].GetBinContent(j) 
		MC = histos_MC[i].GetBinContent(j)
		ratio = 0
		error = 0
		if(MC != 0 and data != 0):
			ratio = data/MC
			error = eff_error(data,MC)
		histo_Divide.SetBinContent(j, ratio)
		histo_Divide.SetBinError(j, error)

	histo_Divide.SetName(histo_Divide.GetName()+'_DataToMC')
	histo_Divide.GetYaxis().SetTitle('Data/MC')
	histo_Divide.Write()

print "number of Lambdahort going into the histos for Data: ", sum(nLambdahortPerEventData)
print "number of Lambdahort going into the histos for MC: ", sum(nLambdahortPerEventMC)

print "mean number of Lambdahort per event in data (x1e3): ", sum(nLambdahortPerEventData)*1e3/len(nLambdahortPerEventData)
print "mean number of Lambdahort per event in MC (x1e3): ", sum(nLambdahortPerEventMC)*1e3/len(nLambdahortPerEventMC)
print "ratio of the nLambdahortData/nLambdahortMC: ", float(sum(nLambdahortPerEventData))/float(sum(nLambdahortPerEventMC))


fOut.Close()
