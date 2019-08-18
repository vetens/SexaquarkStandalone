from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas

fIns = [
'/user/jdeclerc/CMSSW_8_0_30/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/test_SingleMuon_FlatTreeV0s.root'
#'/user/jdeclerc/CMSSW_8_0_30/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/test_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root'
]

OutputFileName = 'analyzedFlatTreeV0s_t_SingleMuon_FlatTreeV0s.root'
#OutputFileName = 'analyzedFlatTreeV0s_test_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root'

h2_RECO_Ks_lxy_vz_kinregA = TH2F('h2_RECO_Ks_lxy_vz_kinregA',"; Ks |vz| (cm); Ks l_{0} (cm)",300,0,300,120,0,120)
h2_RECO_Ks_pt_pz_kinregA= TH2F('h2_RECO_Ks_pt_pz_kinregA',"; Ks pt (GeV); Ks |pz| (GeV)",100,0,10,100,0,40)
h2_RECO_Ks_dxy_dz_min_kinregA= TH2F('h2_RECO_Ks_dxy_dz_min_kinregA',"; Ks dxy (cm); Ks dz_min (cm)",100,-2,2,100,-40,40)
h2_RECO_Ks_eta_phi_kinregA= TH2F('h2_RECO_Ks_eta_phi_kinregA',"; Ks #eta (cm); Ks #phi (cm)",100,-4,4,100,-4,4)


for fIn in fIns:

	fIn = TFile(fIn,'read')

	treeKs = fIn.Get('FlatTreeProducerV0s/FlatTreeKs') 
	treeLambda = fIn.Get('FlatTreeProducerV0s/FlatTreeLambda') 

	nEntriesKs = treeKs.GetEntries()
	nEntriesLambda = treeLambda.GetEntries()

	print 'Number of entries in the treeKs: ', nEntriesKs
	print 'Number of entries in the treeLambda: ', nEntriesLambda

	for i in range(0,nEntriesKs):
		treeKs.GetEntry(i)
		#only mass cuts:
		#if(treeKs._Ks_mass[0] > 0.48 and treeKs._Ks_mass[0] < 0.52):
		#trying to get lxy, vz flat:
		#if(treeKs._Ks_mass[0] > 0.48 and treeKs._Ks_mass[0] < 0.52 and 1 < treeKs._Ks_pt[0] and treeKs._Ks_pt[0] < 4 and 1 < abs(treeKs._Ks_pz[0]) and abs(treeKs._Ks_pz[0]) < 10  and abs(treeKs._Ks_dxy[0]) < 0.05 and 0 < abs(treeKs._Ks_dz_min[0]) and abs(treeKs._Ks_dz_min[0]) < 0.1 and 0 < treeKs._Ks_eta[0] and treeKs._Ks_eta[0] < 2 and -1.5 < treeKs._Ks_phi[0] and treeKs._Ks_phi[0] < 0):
		#trying to get dxy flat:
		if(treeKs._Ks_mass[0] > 0.48 and treeKs._Ks_mass[0] < 0.52 and 2 < treeKs._Ks_pt[0] and treeKs._Ks_pt[0] < 3 and 0 < abs(treeKs._Ks_pz[0]) and abs(treeKs._Ks_pz[0]) < 10  and 10 < treeKs._Ks_Lxy[0] and treeKs._Ks_Lxy[0] < 20 and 30 < abs(treeKs._Ks_vz[0]) and abs(treeKs._Ks_vz[0]) < 80 and abs(treeKs._Ks_eta[0]) < 2 ):
			h2_RECO_Ks_lxy_vz_kinregA.Fill(treeKs._Ks_vz[0],treeKs._Ks_Lxy[0])
			h2_RECO_Ks_pt_pz_kinregA.Fill(treeKs._Ks_pt[0],treeKs._Ks_pz[0])
			h2_RECO_Ks_dxy_dz_min_kinregA.Fill(treeKs._Ks_dxy[0],treeKs._Ks_dz_min[0])
			h2_RECO_Ks_eta_phi_kinregA.Fill(treeKs._Ks_eta[0],treeKs._Ks_phi[0])
		


fOut = TFile(OutputFileName,'RECREATE')

h2_RECO_Ks_lxy_vz_kinregA.Write()
h2_RECO_Ks_pt_pz_kinregA.Write()
h2_RECO_Ks_dxy_dz_min_kinregA.Write()
h2_RECO_Ks_eta_phi_kinregA.Write()

#h2_antiS_vz_lxy.Write()
#
#c_antiS_mass = TCanvas("c_antiS_mass","");
#h_antiS_mass.DrawNormalized()
#h_antiS_mass.Fit("gaus")
#h_antiS_mass.DrawNormalized()
#c_antiS_mass.Write()


fOut.Write()
fOut.Close()
