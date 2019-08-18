from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, TLegend, TProfile2D


fIn = TFile('/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Step1/CRAB_SimSexaqWithPU2016NeutrinoGun_tryToFix_8_trial10/crab_Step1SexaqWithPU2016NeutrinoGun_tryToFix_8_14072019_v4/190714_185701/FlatTreeGENSIM_WithPU2016NeutrinoGun_tryToFix_8.root', 'read')
tree = fIn.Get('FlatTreeProducerGENSIM/FlatTreeGENLevel') 
treeAllAntiS = fIn.Get('FlatTreeProducerGENSIM/FlatTreeGENLevelAllAntiS') 

fOut = TFile('macroFlatTreeGENSIM_WithPU2016NeutrinoGun_tryToFix_8.root','RECREATE')

eta_distribution_antiS_dir = fOut.mkdir("eta_distribution_antiS")
eta_distribution_antiS_dir.cd()

h_eta_all_AntiS = TH1F('h_eta_all_AntiS','; #eta #bar{S}; #Entries',160,-8,8)
h_eta_all_AntiS_cut_eta_smaller4 = TH1F('h_eta_all_AntiS_cut_eta_smaller4','; #eta #bar{S}; #Entries',160,-8,8)


#first do this small loop which runs over the tree containing all the antiS (i.e. also the ones which do not go the correct granddaughters)
for i in range(0,treeAllAntiS.GetEntries()):
	treeAllAntiS.GetEntry(i)
	h_eta_all_AntiS.Fill(treeAllAntiS._S_eta_all[0])
	if(abs(treeAllAntiS._S_eta_all[0])<4):
		h_eta_all_AntiS_cut_eta_smaller4.Fill(treeAllAntiS._S_eta_all[0])
c_eta_all_AntiS = TCanvas("c_eta_all_AntiS","");
h_eta_all_AntiS.DrawNormalized()
h_eta_all_AntiS.Write()
c_eta_all_AntiS_cut_eta_smaller4 = TCanvas("c_eta_all_AntiS_cut_eta_smaller4","");
h_eta_all_AntiS_cut_eta_smaller4.DrawNormalized()
h_eta_all_AntiS_cut_eta_smaller4.Write()


l_h_lxy_eta = []
l_h_vz_eta = []
l_h_lxy_n_loops = []
l_h_eta_n_loops = []
l_h_vz_n_loops = []
eta_unit = 0.2
for i in range(0,400):
	lower_eta = i*0.2
	higher_eta = (i+1)*0.2
	h_lxy_eta = TH1F('h_lxy_eta'+str(lower_eta)+'_to_'+str(higher_eta),'h_lxy_eta; l_{0} #bar{S} interaction vertex (cm); #Entries',4*140,0,140)
	l_h_lxy_eta.append(h_lxy_eta)
	h_vz_eta = TH1F('h_vz_eta'+str(lower_eta)+'_to_'+str(higher_eta),'h_vz_eta; vz #bar{S} interaction vertex (cm); #Entries',600,-300,300)
	l_h_vz_eta.append(h_vz_eta)
	h_lxy_n_loops = TH1F('h_lxy_n_loops'+str(i),'h_lxy_n_loops; l_{0} #bar{S} interaction vertex (cm); #Entries',4*140,0,140)
	h_eta_n_loops = TH1F('h_eta_n_loops'+str(i),'h_eta_n_loops; #eta #bar{S} interaction vertex (cm); #Entries',90,-4.5,4.5)
	h_vz_n_loops = TH1F('h_vz_n_loops'+str(i),'h_vz_n_loops; vz #bar{S} interaction vertex (cm); #Entries',600,-300,300)
	l_h_lxy_n_loops.append(h_lxy_n_loops)
	l_h_eta_n_loops.append(h_eta_n_loops)
	l_h_vz_n_loops.append(h_vz_n_loops)


	

nEntries = tree.GetEntries()
print 'Number of entries in the tree: ', nEntries
for i in range(0,nEntries):
	tree.GetEntry(i)
	for j in range(0,len(tree._S_lxy_interaction_vertex)):
		eta_range = int(abs(tree._S_eta[j])/eta_unit)
		l_h_lxy_eta[eta_range].Fill(tree._S_lxy_interaction_vertex[j])
		l_h_vz_eta[eta_range].Fill(tree._S_vz_interaction_vertex[j])
		if(int(tree._S_n_loops[j])<400):
			if(abs(tree._S_eta[j])<1.1):
				l_h_lxy_n_loops[int(tree._S_n_loops[j]/100)].Fill(tree._S_lxy_interaction_vertex[j])
			else:
				l_h_lxy_n_loops[10].Fill(tree._S_lxy_interaction_vertex[j])
			l_h_eta_n_loops[int(tree._S_n_loops[j])].Fill(tree._S_eta[j])
			l_h_vz_n_loops[int(tree._S_n_loops[j])].Fill(tree._S_vz_interaction_vertex[j])
		else:
			l_h_lxy_n_loops[10].Fill(tree._S_lxy_interaction_vertex[j])


lxy_eta_dir = fOut.mkdir("lxy_eta")
lxy_eta_dir.cd()
for h in l_h_lxy_eta:
	h.Write()

vz_eta_dir = fOut.mkdir("vz_eta")
vz_eta_dir.cd()
for h in l_h_vz_eta:
	h.Write()

lxy_n_loops_dir = fOut.mkdir("lxy_n_loops")
lxy_n_loops_dir.cd()
for h in l_h_lxy_n_loops:
	h.Write()

eta_n_loops_dir = fOut.mkdir("eta_n_loops")
eta_n_loops_dir.cd()
for h in l_h_eta_n_loops:
	h.Write()

vz_n_loops_dir = fOut.mkdir("vz_n_loops")
vz_n_loops_dir.cd()
for h in l_h_vz_n_loops:
	h.Write()

overall_dir = fOut.mkdir("overall")
overall_dir.cd()
#plot the n_loops versus eta:
h2_n_loops_vs_eta = TH2F("h2_n_loops_vs_eta","; #eta(#bar{S}); #Loops; #Entries",80,-4,4,600,0,600)
#plot the number of interactions in each loop
h_n_interactions_vs_n_loops = TH1F("h_n_interactions_vs_n_loops","; #Loops; #bar{S} interacting after #Loops",1600,0,1600)
h_n_interactions_vs_n_loops_barrel = TH1F("h_n_interactions_vs_n_loops_barrel","; #Loops; #bar{S} interacting after #Loops",1600,0,1600)
h_n_interactions_vs_n_loops_endcap = TH1F("h_n_interactions_vs_n_loops_endcap","; #Loops; #bar{S} interacting after #Loops",1600,0,1600)
h_n_interactions_vs_n_loops_endcap_more_forward_than_endcap = TH1F("h_n_interactions_vs_n_loops_endcap_more_forward_than_endcap","; #Loops; #bar{S} interacting after #Loops",1600,0,1600)
#2D plots of the interaction vertex location
h2_interaction_vertex_vx_vy = TH2F("h2_interaction_vertex_vx_vy","; v_{x} #bar{S} interaction vertex (cm); v_{y} #bar{S} interaction vertex (cm);#Entries",2500,-125,125,2500,-125,125)
h2_interaction_vertex_vz_vx = TH2F("h2_interaction_vertex_vz_vx","; v_{z} #bar{S} interaction vertex (cm); v_{x} #bar{S} interaction vertex (cm);#Entries",2500,-125,125,2500,-125,125)
h2_interaction_vertex_vz_vy = TH2F("h2_interaction_vertex_vz_vy","; v_{z} #bar{S} interaction vertex (cm); v_{y} #bar{S} interaction vertex (cm);#Entries",2500,-125,125,2500,-125,125)
h2_interaction_vertex_vz_lxy = TH2F("h2_interaction_vertex_vz_lxy","; v_{z} #bar{S} interaction vertex (cm); l_{0} #bar{S} interaction vertex (cm);#Entries",2500,-125,125,1250,0,125)

for i in range(0,nEntries):
	tree.GetEntry(i)
	for j in range(0,len(tree._S_lxy_interaction_vertex)):
		h2_n_loops_vs_eta.Fill(tree._S_eta[j],tree._S_n_loops[j])
		h2_interaction_vertex_vx_vy.Fill(tree._S_vx_interaction_vertex[j],tree._S_vy_interaction_vertex[j])
		h2_interaction_vertex_vz_vx.Fill(tree._S_vz_interaction_vertex[j],tree._S_vx_interaction_vertex[j])
		h2_interaction_vertex_vz_vy.Fill(tree._S_vz_interaction_vertex[j],tree._S_vy_interaction_vertex[j])
		h2_interaction_vertex_vz_lxy.Fill(tree._S_vz_interaction_vertex[j],tree._S_lxy_interaction_vertex[j])
		h_n_interactions_vs_n_loops.Fill(tree._S_n_loops[j])
		if(abs(tree._S_eta[j])<1.1):
			h_n_interactions_vs_n_loops_barrel.Fill(tree._S_n_loops[j])
		elif(abs(tree._S_eta[j])>1.1 and abs(tree._S_eta[j])<2.5):
			h_n_interactions_vs_n_loops_endcap.Fill(tree._S_n_loops[j])
		else:
			h_n_interactions_vs_n_loops_endcap_more_forward_than_endcap.Fill(tree._S_n_loops[j])



momenta_daughters_and_grandaughters_dir = fOut.mkdir("momenta_daughters_and_grandaughters")
momenta_daughters_and_grandaughters_dir.cd()
h_pt_Ks = TH1F("h_pt_Ks",";p_{t} (GeV); #Entries",80,0,8)
h_pt_Ks_daug0 = TH1F("h_pt_Ks_daug0",";p_{t} (GeV); #Entries",80,0,8)
h_pt_Ks_daug1 = TH1F("h_pt_Ks_daug1",";p_{t} (GeV); #Entries",80,0,8)
h_pz_Ks = TH1F("h_pz_Ks",";p_{z} (GeV); #Entries",800,0,80)
h_pz_Ks_daug0 = TH1F("h_pz_Ks_daug0",";p_{z} (GeV); #Entries",800,0,80)
h_pz_Ks_daug1 = TH1F("h_pz_Ks_daug1",";p_{z} (GeV); #Entries",800,0,80)
h_dxy_Ks_daug0 = TH1F("h_dxy_Ks_daug0",";d_{0} (cm); #Entries",200,-10,10)
h_dxy_Ks_daug1 = TH1F("h_dxy_Ks_daug1",";d_{0} (cm); #Entries",200,-10,10)
h_dz_Ks_daug0 = TH1F("h_dz_Ks_daug0",";d_{z} (cm); #Entries",400,-100,100)
h_dz_Ks_daug1 = TH1F("h_dz_Ks_daug1",";d_{z} (cm); #Entries",400,-100,100)

h_pt_AntiLambda = TH1F("h_pt_AntiLambda",";p_{t} #bar{#Lambda}^{0} (GeV); #Entries",80,0,8)
h_pt_AntiLambda_AntiProton = TH1F("h_pt_AntiLambda_AntiProton",";p_{t} #bar{#Lambda}^{0}-#bar(p)(GeV); #Entries",80,0,8)
h_pt_AntiLambda_Pion = TH1F("h_pt_AntiLambda_Pion",";p_{t} #bar{#Lambda}^{0}(GeV)-#pi^{+}; #Entries",80,0,8)
h_pz_AntiLambda = TH1F("h_pz_AntiLambda",";p_{z} #bar{#Lambda}^{0}(GeV); #Entries",800,0,80)
h_pz_AntiLambda_AntiProton = TH1F("h_pz_AntiLambda_AntiProton",";p_{z} #bar{#Lambda}^{0}-#bar(p)(GeV); #Entries",800,0,80)
h_pz_AntiLambda_Pion = TH1F("h_pz_AntiLambda_Pion",";p_{z} #bar{#Lambda}^{0}(GeV)-#pi^{+}; #Entries",800,0,80)

PCA_granddaughters_dir = fOut.mkdir("PCA_granddaughters")
PCA_granddaughters_dir.cd()

h_dxy_AntiLambda_AntiProton = TH1F("h_dxy_AntiLambda_AntiProton",";d_{0} #bar{#Lambda}^{0}-#bar{p} (cm); #Entries",200,-10,10)
h_dxy_AntiLambda_Pion = TH1F("h_dxy_AntiLambda_Pion",";d_{0} #bar{#Lambda}^{0}(GeV)-#pi^{+} (cm); #Entries",200,-10,10)
h_dz_AntiLambda_AntiProton = TH1F("h_dz_AntiLambda_AntiProton",";d_{z} #bar{#Lambda}^{0}-#bar(p) (cm); #Entries",400,-100,100)
h_dz_AntiLambda_Pion = TH1F("h_dz_AntiLambda_Pion",";d_{z} #bar{#Lambda}^{0}-#pi^{+} (cm); #Entries",400,-100,100)

decay_vertex_dir = fOut.mkdir("decay_vertex_V0s")
decay_vertex_dir.cd()
h_lxy_creation_vertex_Ks_daughters = TH1F("h_lxy_creation_vertex_Ks_daughters",";l_{0} (cm) ",200,0,200)
h_vz_creation_vertex_Ks_daughters = TH1F("h_vz_creation_vertex_Ks_daughters",";v_{z} (cm)",200,-400,400)
h_lxy_creation_vertex_AntiLambda_daughters = TH1F("h_lxy_creation_vertex_AntiLambda_daughters",";l_{0} (cm)",200,0,200)
h_vz_creation_vertex_AntiLambda_daughters = TH1F("h_vz_creation_vertex_AntiLambda_daughters",";v_{z} (cm)",200,-400,400)

h2_vx_vy_creation_vertex_Ks_daughters = TH2F("h2_vx_vy_creation_vertex_Ks_daughters",";v_{x} decay vertex K_{S}^{0} (cm); v_{y} decay vertex K_{S}^{0} (cm);#Entries",240,-120,120,240,-120,120)
h2_vz_vx_creation_vertex_Ks_daughters = TH2F("h2_vz_vx_creation_vertex_Ks_daughters",";v_{z} decay vertex K_{S}^{0} (cm);v_{x} decay vertex K_{S}^{0} (cm);#Entries",600,-300,300,240,-120,120)
h2_vz_lxy_creation_vertex_Ks_daughters = TH2F("h2_vz_lxy_creation_vertex_Ks_daughters",";v_{z} decay vertex K_{S}^{0} (cm);l_{0} decay vertex K_{S}^{0} (cm);#Entries",600,-300,300,120,0,120)
h2_vx_vy_creation_vertex_AntiLambda_daughters = TH2F("h2_vx_vy_creation_vertex_AntiLambda_daughters",";v_{x} decay vertex #bar{#Lambda}^{0} (cm);v_{y} decayvertex #bar{#Lambda}^{0} (cm);#Entries",240,-120,120,240,-120,120)
h2_vz_vx_creation_vertex_AntiLambda_daughters = TH2F("h2_vz_vx_creation_vertex_AntiLambda_daughters",";v_{z} decay vertex #bar{#Lambda}^{0} (cm);v_{x} decay vertex #bar{#Lambda}^{0} (cm);#Entries",600,-300,300,240,-120,120)
h2_vz_lxy_creation_vertex_AntiLambda_daughters = TH2F("h2_vz_lxy_creation_vertex_AntiLambda_daughters",";v_{z} decay vertex #bar{#Lambda}^{0} (cm);l_{0} decay vertex #bar{#Lambda}^{0} (cm);#Entries",600,-300,300,120,0,120)

prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt",";v_{x} decay vertex K_{S}^{0} (cm); v_{y} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",240,-120,120,240,-120,120)
prof2_vz_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vz_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt",";v_{z} decay vertex K_{S}^{0} (cm); v_{x} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",600,-300,300,240,-120,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",600,-300,300,120,0,120)
prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt",";v_{x} decay vertex #bar{#Lambda}^{0} (cm); v_{y} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",240,-120,120,240,-120,120)
prof2_vz_vx_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vz_vx_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); v_{x} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",600,-300,300,240,-120,120)
prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); l_{0} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",600,-300,300,120,0,120)
prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt",";v_{x} decay vertex #bar{#Lambda}^{0} (cm); v_{y} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",240,-120,120,240,-120,120)
prof2_vz_vx_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vz_vx_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); v_{x} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",600,-300,300,240,-120,120)
prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); l_{0} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",600,-300,300,120,0,120)

prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt = TProfile2D("prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt",";v_{x} decay vertex K_{S}^{0} (cm); v_{y} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",240,-120,120,240,-120,120)
prof2_vz_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt",";v_{z} decay vertex K_{S}^{0} (cm); v_{x} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",600,-300,300,240,-120,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",600,-300,300,120,0,120)
prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt = TProfile2D("prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt",";v_{x} decay vertex #bar{#Lambda}^{0} (cm); v_{y} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",240,-120,120,240,-120,120)
prof2_vz_vx_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_vx_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); v_{x} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",600,-300,300,240,-120,120)
prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); l_{0} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",600,-300,300,120,0,120)
prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt = TProfile2D("prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt",";v_{x} decay vertex #bar{#Lambda}^{0} (cm); v_{y} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",240,-120,120,240,-120,120)
prof2_vz_vx_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_vx_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); v_{x} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",600,-300,300,240,-120,120)
prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); l_{0} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",600,-300,300,120,0,120)

h_antiS_lxy = TH1F('h_antiS_lxy','; #bar{S} interaction vertex l_{0} (cm) ; #Entries',120,0,120)
h_antiS_vz = TH1F('h_antiS_vz','; #bar{S} interaction vertex |v_{z}| (cm); #Entries',50,0,200)

for i in range(0,nEntries):
	tree.GetEntry(i)
	h_antiS_lxy.Fill(tree._S_lxy_interaction_vertex[0])	
	h_antiS_vz.Fill(abs(tree._S_vz_interaction_vertex[0]))	
	h_pt_Ks.Fill(tree._Ks_pt[0])
	h_pt_Ks_daug0.Fill(tree._GEN_Ks_daughter0_pt[0])
	h_pt_Ks_daug1.Fill(tree._GEN_Ks_daughter1_pt[0])
	h_pz_Ks.Fill(abs(tree._Ks_pz[0]))
	h_pz_Ks_daug0.Fill(abs(tree._GEN_Ks_daughter0_pz[0]))
	h_pz_Ks_daug1.Fill(abs(tree._GEN_Ks_daughter1_pz[0]))
	h_dxy_Ks_daug0.Fill(tree._GEN_Ks_daughter0_dxy[0])
	h_dxy_Ks_daug1.Fill(tree._GEN_Ks_daughter1_dxy[0])
	h_dz_Ks_daug0.Fill(tree._GEN_Ks_daughter0_dz[0])
	h_dz_Ks_daug1.Fill(tree._GEN_Ks_daughter1_dz[0])


	h_pt_AntiLambda.Fill(tree._Lambda_pt[0])
	h_pt_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_pt[0])
	h_pt_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_pt[0])
	h_pz_AntiLambda.Fill(abs(tree._Lambda_pz[0]))
	h_pz_AntiLambda_AntiProton.Fill(abs(tree._GEN_AntiLambda_AntiProton_pz[0]))
	h_pz_AntiLambda_Pion.Fill(abs(tree._GEN_AntiLambda_Pion_pz[0]))
	h_dxy_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_dxy[0])
	h_dxy_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_dxy[0])
	h_dz_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_dz[0])
	h_dz_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_dz[0])
	
	#need only one of the granddaughter of each V0, as their creation vertex is the same
	h_lxy_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_lxy[0])
	h_vz_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_vz[0])
	h_lxy_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_lxy[0])
	h_vz_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vz[0])

	h2_vx_vy_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_vx[0],tree._GEN_Ks_daughter0_vy[0])
	h2_vz_vx_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_vx[0])
	h2_vz_lxy_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0])
	h2_vx_vy_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_vy[0])
	h2_vz_vx_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_vx[0])
	h2_vz_lxy_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_lxy[0])

	#fill the tprofiles
	#the Ks daughters are both pions, so can put them in same tprofiles
	if(tree._GEN_Ks_daughter0_pt[0]>0.1 and tree._GEN_Ks_daughter0_pt[0]<0.9):
		prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter0_vx[0],tree._GEN_Ks_daughter0_vy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
		prof2_vz_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_vx[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
	if(tree._GEN_Ks_daughter1_pt[0]>0.1 and tree._GEN_Ks_daughter1_pt[0]<0.9):
		prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter1_vx[0],tree._GEN_Ks_daughter1_vy[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])
		prof2_vz_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter1_vz[0],tree._GEN_Ks_daughter1_vx[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter1_vz[0],tree._GEN_Ks_daughter1_lxy[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])

	if(tree._GEN_AntiLambda_AntiProton_pt[0]>0.1 and tree._GEN_AntiLambda_AntiProton_pt[0]<0.9):
		prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_vy[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
		prof2_vz_vx_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_lxy[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
	if(tree._GEN_AntiLambda_Pion_pt[0]>0.1 and tree._GEN_AntiLambda_Pion_pt[0]<0.9):
		prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_Pion_vx[0],tree._GEN_AntiLambda_Pion_vy[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])
		prof2_vz_vx_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_Pion_vz[0],tree._GEN_AntiLambda_Pion_vx[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_Pion_vz[0],tree._GEN_AntiLambda_Pion_lxy[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])

	if(tree._GEN_Ks_daughter0_pt[0]>0.9):
		prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter0_vx[0],tree._GEN_Ks_daughter0_vy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
		prof2_vz_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_vx[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
	if(tree._GEN_Ks_daughter1_pt[0]>0.9):
		prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter1_vx[0],tree._GEN_Ks_daughter1_vy[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])
		prof2_vz_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter1_vz[0],tree._GEN_Ks_daughter1_vx[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter1_vz[0],tree._GEN_Ks_daughter1_lxy[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])
	
	if(tree._GEN_AntiLambda_AntiProton_pt[0]>0.9):
		prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_vy[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
		prof2_vz_vx_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_lxy[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
	if(tree._GEN_AntiLambda_Pion_pt[0]>0.9):
		prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_Pion_vx[0],tree._GEN_AntiLambda_Pion_vy[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])
		prof2_vz_vx_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_Pion_vz[0],tree._GEN_AntiLambda_Pion_vx[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_Pion_vz[0],tree._GEN_AntiLambda_Pion_lxy[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])

momenta_daughters_and_grandaughters_dir.cd()

c_antiS_lxy = TCanvas("c_antiS_lxy","");
h_antiS_lxy.DrawNormalized()
c_antiS_lxy.Write()

c_antiS_vz = TCanvas("c_antiS_vz","");
h_antiS_vz.DrawNormalized()
c_antiS_vz.Write()



c_ptKs = TCanvas("c_ptKs","");
h_pt_Ks_daug0.SetLineColor(38)
h_pt_Ks_daug0.SetFillColorAlpha(38,0.5)
h_pt_Ks_daug0.SetLineWidth(3)
h_pt_Ks.SetLineColor(46)
h_pt_Ks.SetFillColorAlpha(46,0.5)
h_pt_Ks.SetLineWidth(3)
h_pt_Ks.DrawNormalized();
h_pt_Ks_daug0.DrawNormalized("same");
legend_ptKs = TLegend(0.1,0.7,0.48,0.9);
legend_ptKs.AddEntry(h_pt_Ks_daug0,"K_{S}^{0} daughters","l");
legend_ptKs.AddEntry(h_pt_Ks,"K_{S}^{0}","l");
legend_ptKs.Draw();
c_ptKs.Write()

c_pzKs = TCanvas("c_pzKs","");
h_pz_Ks_daug0.SetLineColor(38)
h_pz_Ks_daug0.SetLineWidth(3)
h_pz_Ks_daug0.SetFillColorAlpha(38,0.5)
h_pz_Ks.SetLineColor(46)
h_pz_Ks.SetFillColorAlpha(46,0.5)
h_pz_Ks.SetLineWidth(3)
h_pz_Ks.DrawNormalized();
h_pz_Ks_daug0.DrawNormalized("same");
legend_pzKs = TLegend(0.1,0.7,0.48,0.9);
legend_pzKs.AddEntry(h_pz_Ks_daug0,"K_{S}^{0} daughters","l");
legend_pzKs.AddEntry(h_pz_Ks,"K_{S}^{0}","l");
legend_pzKs.Draw();
c_pzKs.Write()

c_pzAntiLambda = TCanvas("c_pzAntiLambda","");
h_pz_AntiLambda_AntiProton.SetLineColor(38)
h_pz_AntiLambda_AntiProton.SetLineWidth(3)
h_pz_AntiLambda_AntiProton.SetFillColorAlpha(38,0.5)
h_pz_AntiLambda_Pion.SetLineColor(46)
h_pz_AntiLambda_Pion.SetLineWidth(3)
h_pz_AntiLambda_Pion.SetFillColorAlpha(46,0.5)
h_pz_AntiLambda.SetLineColor(41)
h_pz_AntiLambda.SetLineWidth(3)
h_pz_AntiLambda.SetFillColorAlpha(41,0.5)
h_pz_AntiLambda.DrawNormalized();
h_pz_AntiLambda_AntiProton.DrawNormalized("same");
h_pz_AntiLambda_Pion.DrawNormalized("same");
legend_pzAntiLambda = TLegend(0.1,0.7,0.48,0.9);
legend_pzAntiLambda.AddEntry(h_pz_AntiLambda,"#bar{#Lambda}^{0}","l");
legend_pzAntiLambda.AddEntry(h_pz_AntiLambda_AntiProton,"#bar{#Lambda}^{0}-#bar{p}","l");
legend_pzAntiLambda.AddEntry(h_pz_AntiLambda_Pion,"#bar{#Lambda}^{0}-#pi^{+}","l");
legend_pzAntiLambda.Draw();
c_pzAntiLambda.Write()

c_ptAntiLambda = TCanvas("c_ptAntiLambda","");
h_pt_AntiLambda_AntiProton.SetLineColor(38)
h_pt_AntiLambda_AntiProton.SetLineWidth(3)
h_pt_AntiLambda_AntiProton.SetFillColorAlpha(38,0.5)
h_pt_AntiLambda_Pion.SetLineColor(46)
h_pt_AntiLambda_Pion.SetLineWidth(3)
h_pt_AntiLambda_Pion.SetFillColorAlpha(46,0.5)
h_pt_AntiLambda.SetLineColor(41)
h_pt_AntiLambda.SetLineWidth(3)
h_pt_AntiLambda.SetFillColorAlpha(41,0.5)
h_pt_AntiLambda.DrawNormalized();
h_pt_AntiLambda_AntiProton.DrawNormalized("same");
h_pt_AntiLambda_Pion.DrawNormalized("same");
legend_ptAntiLambda = TLegend(0.1,0.7,0.48,0.9);
legend_ptAntiLambda.AddEntry(h_pt_AntiLambda,"#bar{#Lambda}^{0}","l");
legend_ptAntiLambda.AddEntry(h_pt_AntiLambda_AntiProton,"#bar{#Lambda}^{0}-#bar{p}","l");
legend_ptAntiLambda.AddEntry(h_pt_AntiLambda_Pion,"#bar{#Lambda}^{0}-#pi^{+}","l");
legend_ptAntiLambda.Draw();
c_ptAntiLambda.Write()

PCA_granddaughters_dir.cd()
c_dxyGranddaugthers = TCanvas("c_dxyGranddaugthers","");
h_dxy_AntiLambda_AntiProton.SetLineColor(38)
h_dxy_AntiLambda_AntiProton.SetLineWidth(3)
h_dxy_AntiLambda_AntiProton.SetFillColorAlpha(38,0.5)
h_dxy_AntiLambda_Pion.SetLineColor(46)
h_dxy_AntiLambda_Pion.SetLineWidth(3)
h_dxy_AntiLambda_Pion.SetFillColorAlpha(46,0.5)
h_dxy_Ks_daug0.SetLineColor(41)
h_dxy_Ks_daug0.SetLineWidth(3)
h_dxy_Ks_daug0.SetFillColorAlpha(41,0.5)
h_dxy_AntiLambda_AntiProton.DrawNormalized()
h_dxy_AntiLambda_Pion.DrawNormalized("same")
h_dxy_Ks_daug0.DrawNormalized("same")
legend_dxyGranddaughters = TLegend(0.1,0.7,0.48,0.9);
legend_dxyGranddaughters.AddEntry(h_dxy_AntiLambda_AntiProton,"#bar{#Lambda}^{0}-#bar{p}","l")
legend_dxyGranddaughters.AddEntry(h_dxy_AntiLambda_Pion,"#bar{#Lambda}^{0}-#pi^{+}","l")
legend_dxyGranddaughters.AddEntry(h_dxy_Ks_daug0,"K_{S}^{0} daughters","l")
legend_dxyGranddaughters.Draw()
c_dxyGranddaugthers.Write()

c_dzGranddaugthers = TCanvas("c_dzGranddaugthers","");
h_dz_AntiLambda_AntiProton.SetLineColor(38)
h_dz_AntiLambda_AntiProton.SetLineWidth(3)
h_dz_AntiLambda_AntiProton.SetFillColorAlpha(38,0.5)
h_dz_AntiLambda_Pion.SetLineColor(46)
h_dz_AntiLambda_Pion.SetLineWidth(3)
h_dz_AntiLambda_Pion.SetFillColorAlpha(46,0.5)
h_dz_Ks_daug0.SetLineColor(41)
h_dz_Ks_daug0.SetLineWidth(3)
h_dz_Ks_daug0.SetFillColorAlpha(41,0.5)
h_dz_AntiLambda_AntiProton.DrawNormalized()
h_dz_AntiLambda_Pion.DrawNormalized("same")
h_dz_Ks_daug0.DrawNormalized("same")
legend_dzGranddaughters = TLegend(0.1,0.7,0.48,0.9);
legend_dzGranddaughters.AddEntry(h_dz_AntiLambda_AntiProton,"#bar{#Lambda}^{0}-#bar{p}","l")
legend_dzGranddaughters.AddEntry(h_dz_AntiLambda_Pion,"#bar{#Lambda}^{0}-#pi^{+}","l")
legend_dzGranddaughters.AddEntry(h_dz_Ks_daug0,"K_{S}^{0} daughters","l")
legend_dzGranddaughters.Draw()
c_dzGranddaugthers.Write()


decay_vertex_dir.cd()

c_lxy_decayVertexV0s = TCanvas("c_lxy_decayVertexV0s","");
h_lxy_creation_vertex_Ks_daughters.SetLineColor(38)
h_lxy_creation_vertex_Ks_daughters.SetFillColorAlpha(38,0.5)
h_lxy_creation_vertex_Ks_daughters.SetLineWidth(3)
h_lxy_creation_vertex_AntiLambda_daughters.SetLineColor(46)
h_lxy_creation_vertex_AntiLambda_daughters.SetFillColorAlpha(46,0.5)
h_lxy_creation_vertex_AntiLambda_daughters.SetLineWidth(3)
h_lxy_creation_vertex_AntiLambda_daughters.DrawNormalized();
h_lxy_creation_vertex_Ks_daughters.DrawNormalized("same");
legend_lxy_decayVertexV0s = TLegend(0.1,0.7,0.48,0.9);
legend_lxy_decayVertexV0s.AddEntry(h_lxy_creation_vertex_Ks_daughters,"K_{S}^{0} decay vertex","l");
legend_lxy_decayVertexV0s.AddEntry(h_lxy_creation_vertex_AntiLambda_daughters,"#bar{#Lambda}^{0} decay vertex","l");
legend_lxy_decayVertexV0s.Draw();
c_lxy_decayVertexV0s.Write()

c_vz_decayVertexV0s = TCanvas("c_vz_decayVertexV0s","");
h_vz_creation_vertex_Ks_daughters.SetLineColor(38)
h_vz_creation_vertex_Ks_daughters.SetFillColorAlpha(38,0.5)
h_vz_creation_vertex_Ks_daughters.SetLineWidth(3)
h_vz_creation_vertex_AntiLambda_daughters.SetLineColor(46)
h_vz_creation_vertex_AntiLambda_daughters.SetFillColorAlpha(46,0.5)
h_vz_creation_vertex_AntiLambda_daughters.SetLineWidth(3)
h_vz_creation_vertex_AntiLambda_daughters.DrawNormalized();
h_vz_creation_vertex_Ks_daughters.DrawNormalized("same");
legend_vz_decayVertexV0s = TLegend(0.1,0.7,0.48,0.9);
legend_vz_decayVertexV0s.AddEntry(h_vz_creation_vertex_Ks_daughters,"K_{S}^{0} decay vertex","l");
legend_vz_decayVertexV0s.AddEntry(h_vz_creation_vertex_AntiLambda_daughters,"#bar{#Lambda}^{0} decay vertex","l");
legend_vz_decayVertexV0s.Draw();
c_vz_decayVertexV0s.Write()


fOut.Write()
fOut.Close()
