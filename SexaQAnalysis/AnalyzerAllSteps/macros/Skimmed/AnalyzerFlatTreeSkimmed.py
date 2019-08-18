from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas

#fIn = TFile('/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/AnalyzerAllSteps/test/wihtMatchingOnHits/test_TrackMatchingOnHits.root', 'read')
fIns = [
#'/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaqWithPU2016NeutrinoGun_tryToFix_8_10072019_v1/FlatTree_WithPU2016NeutrinoGun_tryToFix_8_10072019_v1.root',
#'/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaqWithPU2016NeutrinoGun_tryToFix_8_trial11_17072019_v1/crab_SkimmingSexaqWithPU2016NeutrinoGun_tryToFix_8_trial11_17072019_v1/190717_083219/FlatTree/FlatTree_Skimmed_trial11.root',
'/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaq_Skimming_trial13_25072019_v1/crab_SkimmingSexaq_trial13_25072019_v1/190725_051522/FlatTree_Skimmed_trial13.root'
]

OutputFileName = 'macroFlatTree.root'

h_antiS_mass = TH1F('h_antiS_mass','; #bar{S} mass (GeV); #Entries',200,0,20)
h_antiS_pt = TH1F('h_antiS_pt','; #bar{S} p_{t} (GeV); #Entries',200,0,20)
h_antiS_pz = TH1F('h_antiS_pz','; #bar{S} |p_{z}| (GeV); #Entries',100,0,100)
h2_antiS_pt_pz = TH2F('h2_antiS_pt_pz','; #bar{S} p_{t} (GeV);#bar{S} |p_{z}| (GeV); #Entries',100,0,10,100,0,100)
h_antiS_eta = TH1F('h_antiS_eta','; #bar{S} #eta ; #Entries',160,-8,8)
h_antiS_lxy = TH1F('h_antiS_lxy','; #bar{S} interaction vertex l_{0} (cm) ; #Entries',100,0,10)
h_antiS_vz = TH1F('h_antiS_vz','; #bar{S} interaction vertex |v_{z}| (cm); #Entries',50,0,200)

h_antiS_eta_pt = TH2F('h_antiS_eta_pt',';#bar{S} #eta; #bar{S} p_{t} (GeV); #Entries',160,-8,8,200,0,20)
h_antiS_eta_pz = TH2F('h_antiS_eta_pz',';#bar{S} #eta; #bar{S} |p_{z}| (GeV); #Entries',160,-8,8,200,0,100)

h_antiS_vx_vy = TH2F('h_antiS_vx_vy',';#bar{S} interaction vertex v_{x} (cm); #bar{S} interaction vertex v_{y} (cm); #Entries',240,-120,120,240,-120,120)
h_antiS_vz_vx = TH2F('h_antiS_vz_vx',';#bar{S} interaction vertex v_{z} (cm); #bar{S} interaction vertex v_{x} (cm); #Entries',240,-120,120,240,-120,120)
h2_antiS_vz_lxy = TH2F('h2_antiS_vz_lxy',';#bar{S} interaction vertex v_{z} (cm); #bar{S} interaction vertex l_{0} (cm); #Entries',240,-120,120,120,0,120)
h2_Ks_vz_lxy = TH2F('h2_Ks_vz_lxy',';Ks decay vertex v_{z} (cm); Ks decay vertex l_{0} (cm); #Entries',240,-120,120,120,0,120)
h2_AntiLambda_vz_lxy = TH2F('h2_AntiLambda_vz_lxy',';#bar{#Lambda}^{0} decay vertex v_{z} (cm); #bar{#Lambda}^{0} decay vertex l_{0} (cm); #Entries',240,-120,120,120,0,120)


h2_Ks_pt_pz = TH2F('h2_Ks_pt_pz','; K^{0}_{S} p_{t} (GeV); K^{0}_{S} |p_{z}| (GeV); #Entries',100,0,10,100,0,100)
h2_AntiLambda_pt_pz = TH2F('h2_AntiLambda_pt_pz','; #bar{#Lambda}^{0} p_{t} (GeV); #bar{#Lambda}^{0} |p_{z}| (GeV); #Entries',100,0,10,100,0,100)
h2_Ks_dxy_dz = TH2F('h2_Ks_dxy_dz','; K^{0}_{S} d_{xy} (cm); K^{0}_{S} d_{z} min (cm); #Entries',100,-20,20,100,-40,40)
h2_AntiLambda_dxy_dz = TH2F('h2_AntiLambda_dxy_dz','; #bar{#Lambda}^{0} d_{xy} (cm); #bar{#Lambda}^{0} d_{z} min (cm); #Entries',100,-20,20,100,-40,40)

for fIn in fIns:
	fIn = TFile(fIn,'read')
	tree = fIn.Get('FlatTreeProducer/FlatTree') 
	nEntries = tree.GetEntries()
	print 'Number of entries in the tree: ', nEntries
	for i in range(0,nEntries):
		tree.GetEntry(i)
		h_antiS_mass.Fill(tree._S_mass[0])
		h_antiS_pt.Fill(tree._S_pt[0])
		h_antiS_pz.Fill(abs(tree._S_pz[0]))
		h2_antiS_pt_pz.Fill(tree._S_pt[0],abs(tree._S_pz[0]))
		h_antiS_eta.Fill(tree._S_eta[0])
		h_antiS_lxy.Fill(tree._S_lxy_interaction_vertex[0])	
		h_antiS_vz.Fill(abs(tree._S_vz[0]))	

		h_antiS_eta_pt.Fill(tree._S_eta[0],tree._S_pt[0])
		h_antiS_eta_pz.Fill(tree._S_eta[0],abs(tree._S_pz[0]))
		h_antiS_vx_vy.Fill(tree._S_vx[0],tree._S_vy[0])
		h_antiS_vz_vx.Fill(tree._S_vz[0],tree._S_vx[0])
		h2_antiS_vz_lxy.Fill(tree._S_vz[0],tree._S_lxy_interaction_vertex[0])
		h2_Ks_vz_lxy.Fill(tree._Ks_vz_decay_vertex[0],tree._Ks_lxy_decay_vertex[0])
		h2_AntiLambda_vz_lxy.Fill(tree._Lambda_vz_decay_vertex[0],tree._Lambda_lxy_decay_vertex[0])

		h2_Ks_pt_pz.Fill(tree._Ks_pt[0] , tree._Ks_pz[0] )
		h2_AntiLambda_pt_pz.Fill(tree._Lambda_pt[0] , tree._Lambda_pz[0] )
		h2_Ks_dxy_dz.Fill(tree._Ks_dxy[0],tree._Ks_dz_min[0])
		h2_AntiLambda_dxy_dz.Fill(tree._Lambda_dxy[0],tree._Lambda_dz_min[0])

fOut = TFile(OutputFileName,'RECREATE')

h2_antiS_vz_lxy.Write()
h2_Ks_vz_lxy.Write()
h2_AntiLambda_vz_lxy.Write()
h2_Ks_pt_pz.Write()
h2_AntiLambda_pt_pz.Write()
h2_Ks_dxy_dz.Write()
h2_AntiLambda_dxy_dz.Write()

c_antiS_mass = TCanvas("c_antiS_mass","");
h_antiS_mass.DrawNormalized()
h_antiS_mass.Fit("gaus")
h_antiS_mass.DrawNormalized()
c_antiS_mass.Write()

c_antiS_pt = TCanvas("c_antiS_pt","");
h_antiS_pt.DrawNormalized()
c_antiS_pt.Write()

c_antiS_pz = TCanvas("c_antiS_pz","");
h_antiS_pz.DrawNormalized()
c_antiS_pz.Write()

c_antiS_eta = TCanvas("c_antiS_eta","");
h_antiS_eta.DrawNormalized()
c_antiS_eta.Write()

c_antiS_lxy = TCanvas("c_antiS_lxy","");
h_antiS_lxy.DrawNormalized()
c_antiS_lxy.Write()

c_antiS_vz = TCanvas("c_antiS_vz","");
h_antiS_vz.DrawNormalized()
c_antiS_vz.Write()

c_antiS_eta_pt = TCanvas("c_antiS_eta_pt","")
h_antiS_eta_pt.DrawNormalized()
c_antiS_eta_pt.Write()

c_antiS_eta_pz = TCanvas("c_antiS_eta_pz","")
h_antiS_eta_pz.DrawNormalized()
c_antiS_eta_pz.Write()

fOut.Write()
fOut.Close()
