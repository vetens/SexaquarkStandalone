from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas

fIns = [
#the standard V0s:
#'/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_standardV0s/combined_FlatTreeV0_DoubleMuonData_RunC_D_E_G_H.root',
#'/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_standardV0s/combined_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root',

#the adapted V0s:
'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_adaptedV0s/combined_FlatTreeV0_DoubleMuonData_Run_C_H.root',
'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_adaptedV0s/combined_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root'
]


fOuts  = ['analyzedFlatTreeV0s_combined_FlatTreeV0_DoubleMuonData_Run_C_H.root',
	  'analyzedFlatTreeV0s_test_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root']

#fist find which of the trees has the smallest number of entries. Data or MC?
fIn1 = TFile(fIns[0],'read')
fIn2 = TFile(fIns[1],'read')
treeKs1 = fIn1.Get('FlatTreeProducerV0s/FlatTreeKs')
treeKs2 = fIn2.Get('FlatTreeProducerV0s/FlatTreeKs')
min_entries = min( treeKs1.GetEntries(), treeKs2.GetEntries() )
print "the nentries in Data: " , treeKs1.GetEntries()
print "the nentries in MC: " , treeKs2.GetEntries()
print "the min_enties: ", min_entries

#first loop over the trees and look at the PV0 distribution in z both for DATA and MC. This is important as these distributions are different for data and MC, so you have to correct for this when calculating e.g. the vz
iFile = 0
h_RECO_PV0_vz_Data = TH1F('h_RECO_PV0_vz_Data'," PV vz (cm); #entries",200,-100,100)
h_RECO_PV0_vz_MC = TH1F('h_RECO_PV0_vz_MC'," PV vz (cm); #entries",200,-100,100)
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

mean_vz_data = h_RECO_PV0_vz_Data.GetMean()
RMS_vz_data  = h_RECO_PV0_vz_Data.GetRMS()
mean_vz_MC   = h_RECO_PV0_vz_MC.GetMean()
RMS_vz_MC    = h_RECO_PV0_vz_MC.GetRMS()

#normalize the histogram of PV in data to later use it for normalization
h_RECO_PV0_vz_Data.Scale(1./h_RECO_PV0_vz_Data.Integral(), "width")
h_RECO_PV0_vz_MC.Scale(1./h_RECO_PV0_vz_MC.Integral(), "width")
fOut = TFile("beampsot_dist.root",'RECREATE')
h_RECO_PV0_vz_Data.Write()
h_RECO_PV0_vz_MC.Write()
fOut.Close()


print "Mean and RMS of PV Data: ", mean_vz_data, " ", RMS_vz_data
print "Mean and RMS of PV MC: ", mean_vz_MC, " ", RMS_vz_MC



iFile = 0
for fIn in fIns:

	fIn = TFile(fIn,'read')

	treeKs = fIn.Get('FlatTreeProducerV0s/FlatTreeKs') 
	treeLambda = fIn.Get('FlatTreeProducerV0s/FlatTreeLambda') 
	treeZ = fIn.Get('FlatTreeProducerV0s/FlatTreeZ')
	treePV = fIn.Get('FlatTreeProducerV0s/FlatTreePV')

	nEntriesKs = treeKs.GetEntries()
	nEntriesLambda = treeLambda.GetEntries()
	nEntriesZ = treeZ.GetEntries()


	print 'Number of entries in the treeKs: ', nEntriesKs
	print 'Number of entries in the treeLambda: ', nEntriesLambda
	print 'Number of entries in the treeZ: ', nEntriesZ

	#check the reweighting
	h_RECO_PV0_vz_MC_reweighted_to_data = TH1F('h_RECO_PV0_vz_MC_reweighted_to_data'," PV vz (cm); #entries",200,-100,100)

	#for all V0s
	h2_RECO_Ks_lxy_vz_all = TH2F('h2_RECO_Ks_lxy_vz_all',"; Ks vz (cm); Ks l_{0} (cm)",300,-300,300,120,0,120)
	h2_RECO_Ks_pt_pz_all= TH2F('h2_RECO_Ks_pt_pz_all',"; Ks pt (GeV); Ks |pz| (GeV)",100,0,10,100,0,40)
	h2_RECO_Ks_dxy_dz_PV_all= TH2F('h2_RECO_Ks_dxy_dz_PV_all',"; Ks dxy (cm); Ks dz_PV (cm)",100,-10,10,100,-40,40)
	h2_RECO_Ks_eta_phi_all= TH2F('h2_RECO_Ks_eta_phi_all',"; Ks #eta (cm); Ks #phi (cm)",100,-4,4,100,-4,4)

	#for V0s pointing at the PV
	h2_RECO_Ks_lxy_vz_pointing = TH2F('h2_RECO_Ks_lxy_vz_pointing',"; Ks vz (cm); Ks l_{0} (cm)",300,-300,300,120,0,120)
	h2_RECO_Ks_pt_pz_pointing = TH2F('h2_RECO_Ks_pt_pz_pointing',"; Ks pt (GeV); Ks |pz| (GeV)",100,0,10,100,0,40)
	h2_RECO_Ks_dxy_dz_PV_pointing = TH2F('h2_RECO_Ks_dxy_dz_PV_pointing',"; Ks dxy (cm); Ks dz_PV (cm)",400,-0.1,0.1,1600,-4,4)
	h2_RECO_Ks_eta_phi_pointing = TH2F('h2_RECO_Ks_eta_phi_pointing',"; Ks #eta (cm); Ks #phi (cm)",100,-4,4,100,-4,4)

	for i in range(0,min_entries):
		if(i%100000 == 0):
			print 'reached event ', i
		
		treeZ.GetEntry(i)
		treeKs.GetEntry(i)
		treePV.GetEntry(i)
		#only select a certain hard scale
		#if(treeZ._Z_ptMuMu[0] < 5):
		#	continue
		#now loop over all the Ks in this event:
		for j in range(0, len(treeKs._Ks_mass)):
			#standard cuts on the V0s
			if(treeKs._Ks_mass[j] > 0.48 and treeKs._Ks_mass[j] < 0.52 and  abs(treeKs._Ks_eta[j]) < 2):

				x_axis_h_RECO_PV0_vz_Data = h_RECO_PV0_vz_Data.GetXaxis()
				bin_prob_h_RECO_PV0_vz_Data = x_axis_h_RECO_PV0_vz_Data.FindBin(treePV._PV0_vz[0]) #the bin which contains the prob of finding treePV._PV0_vz[0] in the data
				prob_h_RECO_PV0_vz_Data = h_RECO_PV0_vz_Data.GetBinContent(bin_prob_h_RECO_PV0_vz_Data) #the actual probability
				prob_h_RECO_PV0_vz_MC = h_RECO_PV0_vz_MC.GetBinContent(bin_prob_h_RECO_PV0_vz_Data) #the actual probability

				if(iFile == 0):#data
                                	h2_RECO_Ks_lxy_vz_all.Fill(treeKs._Ks_vz[j]-treePV._PV0_vz[0],treeKs._Ks_Lxy[j])
                                else:#MC
                                	h2_RECO_Ks_lxy_vz_all.Fill(treeKs._Ks_vz[j]-treePV._PV0_vz[0],treeKs._Ks_Lxy[j],prob_h_RECO_PV0_vz_Data/prob_h_RECO_PV0_vz_MC)
					h_RECO_PV0_vz_MC_reweighted_to_data.Fill(treePV._PV0_vz[0], prob_h_RECO_PV0_vz_Data/prob_h_RECO_PV0_vz_MC)

				h2_RECO_Ks_pt_pz_all.Fill(treeKs._Ks_pt[j],treeKs._Ks_pz[j])
				h2_RECO_Ks_dxy_dz_PV_all.Fill(treeKs._Ks_dxy[j],treeKs._Ks_dz_PV[j])
				h2_RECO_Ks_eta_phi_all.Fill(treeKs._Ks_eta[j],treeKs._Ks_phi[j])

				#for V0 pointing at the PV0
				if(abs(treeKs._Ks_dxy[j]) < 0.015  and abs(treeKs._Ks_dz_PV[j]) < 0.01 ):
					if(iFile == 0):#data
						h2_RECO_Ks_lxy_vz_pointing.Fill(treeKs._Ks_vz[j]-treePV._PV0_vz[0],treeKs._Ks_Lxy[j])
					else:#MC
						h2_RECO_Ks_lxy_vz_pointing.Fill(treeKs._Ks_vz[j]-treePV._PV0_vz[0],treeKs._Ks_Lxy[j],prob_h_RECO_PV0_vz_Data/prob_h_RECO_PV0_vz_MC)
					h2_RECO_Ks_pt_pz_pointing.Fill(treeKs._Ks_pt[j],treeKs._Ks_pz[j])
					h2_RECO_Ks_dxy_dz_PV_pointing.Fill(treeKs._Ks_dxy[j],treeKs._Ks_dz_PV[j])
					h2_RECO_Ks_eta_phi_pointing.Fill(treeKs._Ks_eta[j],treeKs._Ks_phi[j])
		


	fOut = TFile(fOuts[iFile],'RECREATE')
	if(iFile == 1):
		h_RECO_PV0_vz_MC_reweighted_to_data.Scale(1./h_RECO_PV0_vz_MC_reweighted_to_data.Integral(), "width")
		h_RECO_PV0_vz_MC_reweighted_to_data.Write()

	h2_RECO_Ks_lxy_vz_all.Write()
	h2_RECO_Ks_pt_pz_all.Write()
	h2_RECO_Ks_dxy_dz_PV_all.Write()
	h2_RECO_Ks_eta_phi_all.Write()
	
	h2_RECO_Ks_lxy_vz_pointing.Write()
	h2_RECO_Ks_pt_pz_pointing.Write()
	h2_RECO_Ks_dxy_dz_PV_pointing.Write()
	h2_RECO_Ks_eta_phi_pointing.Write()

	iFile += 1

#h2_antiS_vz_lxy.Write()
#
#c_antiS_mass = TCanvas("c_antiS_mass","");
#h_antiS_mass.DrawNormalized()
#h_antiS_mass.Fit("gaus")
#h_antiS_mass.DrawNormalized()
#c_antiS_mass.Write()


	fOut.Write()
	fOut.Close()
