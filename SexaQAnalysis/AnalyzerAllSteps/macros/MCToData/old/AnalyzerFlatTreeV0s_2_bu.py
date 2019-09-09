from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, gROOT
import numpy as np
import copy

fIns = [
#the standard V0s:
#'/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_standardV0s/combined_FlatTreeV0_DoubleMuonData_RunC_D_E_G_H.root',
#'/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_standardV0s/combined_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root',

#the adapted V0s:
'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/combined_FlatTreeV0_DoubleMuonData_Run_C_H.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunC_v1/combined_FlatTreeV0_DoubleMuonData_Run_C.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunD_v2/combined_FlatTreeV0_DoubleMuonData_Run_D.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunE_v3/combined_FlatTreeV0_DoubleMuonData_Run_E.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunH_v6/combined_FlatTreeV0_DoubleMuonData_Run_H.root',
#'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/RunG_v5/combined_FlatTreeV0_DoubleMuonData_Run_G.root',
'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_adaptedV0s/combined_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root'
]

def eff_error(data,MC):
	return data/MC*np.sqrt(1/MC+1/data)	


#fist find which of the trees has the smallest number of entries. Data or MC?
fIn1 = TFile(fIns[0],'read')
fIn2 = TFile(fIns[1],'read')
treeKs1 = fIn1.Get('FlatTreeProducerV0s/FlatTreeKs')
treeKs2 = fIn2.Get('FlatTreeProducerV0s/FlatTreeKs')
min_entries = min( treeKs1.GetEntries(), treeKs2.GetEntries() )
#min_entries = 100000 
print "the nentries in Data: " , treeKs1.GetEntries()
print "the nentries in MC: " , treeKs2.GetEntries()
print "the min_enties: ", min_entries

#get the means of the PV_z distribution for data and MC. Use these values to offset the distributions later
treePV1 = fIn1.Get('FlatTreeProducerV0s/FlatTreePV')
treePV1.Draw("_PV0_vz>>hist_PVz_Data")
hist_PVz_Data = gROOT.FindObject('hist_PVz_Data')
mean_hist_PVz_Data = hist_PVz_Data.GetMean()
print "mean PV_z for Data ", mean_hist_PVz_Data

treePV2 = fIn2.Get('FlatTreeProducerV0s/FlatTreePV')
treePV2.Draw("_PV0_vz>>hist_PVz_MC")
hist_PVz_MC = gROOT.FindObject('hist_PVz_MC')
mean_hist_PVz_MC = hist_PVz_MC.GetMean()
print "mean PV_z for MC ", mean_hist_PVz_MC

#get the means of the Ks_vz distributions for data and MC. Use these values to offset the distributions later. 
treeKs1.Draw("_Ks_vz>>hist_Ks_vz_Data")
hist_Ks_vz_Data = gROOT.FindObject('hist_Ks_vz_Data')
mean_hist_Ks_vz_Data = hist_Ks_vz_Data.GetMean()
print "mean Ks_vz for Data ", mean_hist_Ks_vz_Data 

treeKs2.Draw("_Ks_vz>>hist_Ks_vz_MC")
hist_Ks_vz_MC = gROOT.FindObject('hist_Ks_vz_MC')
mean_hist_Ks_vz_MC = hist_Ks_vz_MC.GetMean()
print "mean Ks_vz for MC ", mean_hist_Ks_vz_MC 

#get the means of the _Ks_dz_PV distributions for data and MC. Use these values to offset the distributions later. 
treeKs1.Draw("_Ks_dz>>hist_Ks_dz_PV_data")
hist_Ks_dz_PV_data= gROOT.FindObject('hist_Ks_dz_PV_data')
mean_hist_Ks_dz_PV_Data = hist_Ks_dz_PV_data.GetMean()
print "mean Ks_dz_PV for Data ", mean_hist_Ks_dz_PV_Data

treeKs2.Draw("_Ks_vz>>hist_Ks_dz_PV_MC")
hist_Ks_dz_PV_MC = gROOT.FindObject('hist_Ks_dz_PV_MC')
mean_hist_Ks_dz_PV_MC = hist_Ks_dz_PV_MC.GetMean()
print "mean Ks_dz_PV for MC ", mean_hist_Ks_dz_PV_MC 


#first loop over the trees and look at the PV0 distribution in z both for DATA and MC. This is important as these distributions are different for data and MC, so you have to correct for this when calculating e.g. the vz
iFile = 0
h_RECO_PV0_vz_Data_shifted = TH1F('h_RECO_PV0_vz_Data_shifted'," ;PV vz (cm) - mean PV vz; #entries",200,-100,100)
h_RECO_PV0_vz_MC_shifted = TH1F('h_RECO_PV0_vz_MC_shifted'," ;PV vz - mean PV vz (cm); #entries",200,-100,100)
for fIn in fIns:
	fIn = TFile(fIn,'read')
	treePV  = fIn.Get('FlatTreeProducerV0s/FlatTreePV')
	for i in range(0,min_entries):
		treePV.GetEntry(i)
		if(iFile == 0):
			h_RECO_PV0_vz_Data_shifted.Fill(treePV._PV0_vz[0]-mean_hist_PVz_Data)	
		else:	
			h_RECO_PV0_vz_MC_shifted.Fill(treePV._PV0_vz[0]-mean_hist_PVz_MC)
	iFile = iFile+1	

mean_vz_data = h_RECO_PV0_vz_Data_shifted.GetMean()
RMS_vz_data  = h_RECO_PV0_vz_Data_shifted.GetRMS()
mean_vz_MC   = h_RECO_PV0_vz_MC_shifted.GetMean()
RMS_vz_MC    = h_RECO_PV0_vz_MC_shifted.GetRMS()

#normalize the histogram of PV in data to later use it for reweighting 
h_RECO_PV0_vz_Data_shifted.Scale(1./h_RECO_PV0_vz_Data_shifted.Integral(), "width")
h_RECO_PV0_vz_MC_shifted.Scale(1./h_RECO_PV0_vz_MC_shifted.Integral(), "width")


print "Mean and RMS of PV Data: ", mean_vz_data, " ", RMS_vz_data
print "Mean and RMS of PV MC: ", mean_vz_MC, " ", RMS_vz_MC

fOut = TFile('output_DataToMC_RunC_H.root','RECREATE')
h_dummy =  TH1F('h_dummy'," x; y",20,0,10)
histos_Data = [h_dummy]*9
histos_MC = [h_dummy]*9


nKshortPerEventData = []
nKshortPerEventMC = []

iFile = 0
for fIn in fIns:

	print "creating plots for file: ", fIn

	fIn = TFile(fIn,'read')

	treeKs = fIn.Get('FlatTreeProducerV0s/FlatTreeKs') 
	treeLambda = fIn.Get('FlatTreeProducerV0s/FlatTreeLambda') 
	treeZ = fIn.Get('FlatTreeProducerV0s/FlatTreeZ')
	treePV = fIn.Get('FlatTreeProducerV0s/FlatTreePV')

	nEntriesKs = treeKs.GetEntries()
	nEntriesLambda = treeLambda.GetEntries()
	nEntriesZ = treeZ.GetEntries()
	nEntriesPV = treePV.GetEntries()

	print 'Number of entries in the treeKs: ', nEntriesKs
	print 'Number of entries in the treeLambda: ', nEntriesLambda
	print 'Number of entries in the treeZ: ', nEntriesZ
	print 'Number of entries in the treePV: ', nEntriesPV

	h_RECO_Ks_vz = TH1F('h_RECO_Ks_vz'," ;Ks vz (cm); #entries",600,-300,300)
	h_RECO_Ks_vz_NOTWeighted = TH1F('h_RECO_Ks_vz_NOTWeighted'," ;Ks vz (cm); #entries",600,-300,300)
	h_RECO_Ks_lxy = TH1F('h_RECO_Ks_lxy'," ;Ks lxy (cm);#entries",120,0,120)
	h_RECO_Ks_pt = TH1F('h_RECO_Ks_pt',"; Ks pt (GeV);#entries",100,0,10)
	h_RECO_Ks_pz= TH1F('h_RECO_Ks_pz',"; Ks pz (GeV);#entries",100,0,40)
	h_RECO_Ks_dxy_PV= TH1F('h_RECO_Ks_dxy_PV',"; Ks dxy (cm); #entries",100,-10,10)
	h_RECO_Ks_dz_PV= TH1F('h_RECO_Ks_dz_PV',";  Ks dz_PV (cm); #entries",100,-40,40)
	h_RECO_Ks_dz_PV_NOTWeighted = TH1F('h_RECO_Ks_dz_PV_NOTWeighted',";  Ks dz_PV (cm); #entries",100,-40,40)
	h_RECO_Ks_eta= TH1F('h_RECO_Ks_eta',"; Ks #eta (cm);#entries",100,-4,4)
	h_RECO_Ks_phi= TH1F('h_RECO_Ks_phi',"; Ks #phi (cm);#entries",100,-4,4)
	histos= [h_RECO_Ks_vz, h_RECO_Ks_vz_NOTWeighted, h_RECO_Ks_lxy, h_RECO_Ks_pt, h_RECO_Ks_pz, h_RECO_Ks_dxy_PV, h_RECO_Ks_dz_PV, h_RECO_Ks_dz_PV_NOTWeighted, h_RECO_Ks_eta, h_RECO_Ks_phi] 
	#check the reweighting
	h_RECO_PV0_vz_MC_shifted_reweighted_to_data = TH1F('h_RECO_PV0_vz_MC_shifted_reweighted_to_data'," PV vz (cm); #entries",200,-100,100)
	h_RECO_PV0_vz_MC_shifted_NOTreweighted_to_data = TH1F('h_RECO_PV0_vz_MC_shifted_NOTreweighted_to_data'," PV vz (cm); #entries",200,-100,100)


	for i in range(0,min_entries):
		if(i%100000 == 0):
			print 'reached event ', i
		
		treeZ.GetEntry(i)
		treeKs.GetEntry(i)
		treePV.GetEntry(i)
		
		nKshortThisEvent = 0
	

		#only select a certain hard scale
#		if(treeZ._Z_ptMuMu[0] < 5):
#		       continue
		
		#now loop over DATA the Ks in this event:
		for j in range(0, len(treeKs._Ks_mass)):
			
			#standard cuts on the V0s
			if(treeKs._Ks_mass[j] > 0.48 and treeKs._Ks_mass[j] < 0.52 and  abs(treeKs._Ks_eta[j]) < 2 ) :
				nKshortThisEvent+=1
				x_axis_h_RECO_PV0_vz_Data_shifted = h_RECO_PV0_vz_Data_shifted.GetXaxis()
				bin_prob_h_RECO_PV0_vz_Data_shifted = x_axis_h_RECO_PV0_vz_Data_shifted.FindBin(treePV._PV0_vz[0]-mean_hist_PVz_Data) #the bin which contains the prob of finding treePV._PV0_vz[0] in the data
				prob_h_RECO_PV0_vz_Data_shifted = h_RECO_PV0_vz_Data_shifted.GetBinContent(bin_prob_h_RECO_PV0_vz_Data_shifted) #the actual probability

				x_axis_h_RECO_PV0_vz_MC_shifted = h_RECO_PV0_vz_MC_shifted.GetXaxis()
				bin_prob_h_RECO_PV0_vz_MC_shifted = x_axis_h_RECO_PV0_vz_MC_shifted.FindBin(treePV._PV0_vz[0]-mean_hist_PVz_MC) #the bin which contains the prob of finding treePV._PV0_vz[0] in the data
				prob_h_RECO_PV0_vz_MC_shifted = h_RECO_PV0_vz_MC_shifted.GetBinContent(bin_prob_h_RECO_PV0_vz_MC_shifted) #the actual probability

				PV_Z_reweighting_factor = prob_h_RECO_PV0_vz_Data_shifted/prob_h_RECO_PV0_vz_MC_shifted

				#print "------------------------------------------------------------------------------------"
				#print "PV_z: ", treePV._PV0_vz[0]
				#print "PV_z shifted Data, prob: ", treePV._PV0_vz[0]-mean_hist_PVz_Data, " ", prob_h_RECO_PV0_vz_Data_shifted  
				#print "PV_z shifted MC, prob: ", treePV._PV0_vz[0]-mean_hist_PVz_MC, " ", prob_h_RECO_PV0_vz_MC_shifted  
				#print "reweighting factor: ",PV_Z_reweighting_factor

				if(iFile == 0): #for data do not reweigh
					histos[0].Fill(treeKs._Ks_vz[j]-mean_hist_Ks_vz_Data)
					histos[1].Fill(treeKs._Ks_vz[j]-mean_hist_Ks_vz_Data)
					histos[6].Fill(treeKs._Ks_dz_PV[j]-mean_hist_Ks_dz_PV_Data)
				if(iFile == 1):#for MC do reweigh
					histos[0].Fill(treeKs._Ks_vz[j]-mean_hist_Ks_vz_MC, PV_Z_reweighting_factor)
					histos[1].Fill(treeKs._Ks_vz[j]-mean_hist_Ks_vz_MC)
					histos[6].Fill(treeKs._Ks_dz_PV[j]-mean_hist_Ks_dz_PV_MC  , PV_Z_reweighting_factor)
				histos[2].Fill(treeKs._Ks_Lxy[j])
				histos[3].Fill(treeKs._Ks_pt[j])
				histos[4].Fill(treeKs._Ks_pz[j])
				histos[5].Fill(treeKs._Ks_dxy[j])
				histos[7].Fill(treeKs._Ks_dz_PV[j])
				histos[8].Fill(treeKs._Ks_eta[j])
				histos[9].Fill(treeKs._Ks_phi[j])

				if(iFile == 1): #for MC
					h_RECO_PV0_vz_MC_shifted_reweighted_to_data.Fill(treePV._PV0_vz[0]-mean_hist_PVz_Data, PV_Z_reweighting_factor)
					h_RECO_PV0_vz_MC_shifted_NOTreweighted_to_data.Fill(treePV._PV0_vz[0]-mean_hist_PVz_MC, 1)

				#for V0 pointing at the PV0
				#if(abs(treeKs._Ks_dxy[j]) < 0.015  and abs(treeKs._Ks_dz_PV[j]) < 0.01 ):

		if(iFile == 0):#for Data
			nKshortPerEventData.append(nKshortThisEvent)
		if(iFile == 1): #for MC
			nKshortPerEventMC.append(nKshortThisEvent)
			
			
	for h in histos:
		h.SetDirectory(0)
	if(iFile == 0):
		histos_Data = histos
	else:
		histos_MC = histos

	iFile += 1


fOut.mkdir('PV_Histos')
fOut.cd('PV_Histos')
h_RECO_PV0_vz_Data_shifted.Write()
h_RECO_PV0_vz_MC_shifted.Write()
h_RECO_PV0_vz_MC_shifted_reweighted_to_data.Scale(1./h_RECO_PV0_vz_MC_shifted_reweighted_to_data.Integral(), "width")
h_RECO_PV0_vz_MC_shifted_reweighted_to_data.Write()
h_RECO_PV0_vz_MC_shifted_NOTreweighted_to_data.Scale(1./h_RECO_PV0_vz_MC_shifted_NOTreweighted_to_data.Integral(), "width")
h_RECO_PV0_vz_MC_shifted_NOTreweighted_to_data.Write()

fOut.mkdir('Data_Histos')
fOut.cd('Data_Histos')
for h in histos_Data:
	h.Write()

fOut.mkdir('MC_Histos')
fOut.cd('MC_Histos')
for h in histos_MC:
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

print "number of Kshort going into the histos for Data: ", sum(nKshortPerEventData)
print "number of Kshort going into the histos for MC: ", sum(nKshortPerEventMC)

print "mean number of Kshort per event in data (x1e3): ", sum(nKshortPerEventData)*1e3/len(nKshortPerEventData)
print "mean number of Kshort per event in MC (x1e3): ", sum(nKshortPerEventMC)*1e3/len(nKshortPerEventMC)
print "ratio of the nKshortData/nKshortMC: ", float(sum(nKshortPerEventData))/float(sum(nKshortPerEventMC))


fOut.Close()
