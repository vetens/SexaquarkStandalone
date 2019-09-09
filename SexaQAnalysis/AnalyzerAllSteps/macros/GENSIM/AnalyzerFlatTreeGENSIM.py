import numpy as np
from ROOT import *

import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

colours = [1,2,4,30,38,41]

maxNEntries = 1e5

plots_output_dir = "plots_GENSIM/"

fIn = TFile('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerGENSIM/test_FlatTree_GENSIM_trial17.root', 'read')
tree = fIn.Get('FlatTreeProducerGENSIM/FlatTreeGENLevel') 
treeAllAntiS = fIn.Get('FlatTreeProducerGENSIM/FlatTreeGENLevelAllAntiS') 

fOut = TFile('macro_test_FlatTree_GENSIM_trial17.root','RECREATE')

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

#now make some plots which validate the looping mechanism for the antiS

#first create some list containing some plots
l_h_lxy_eta = []
l_h_vz_eta = []
l_h_lxy_n_loops = []
l_h_eta_n_loops = []
l_h_vz_n_loops = []
eta_unit = 0.2
for i in range(0,400):
	lower_eta = i*eta_unit
	higher_eta = (i+1)*eta_unit
	h_lxy_eta = TH1F('h_lxy_eta'+str(lower_eta)+'_to_'+str(higher_eta),'h_lxy_eta; l_{0} #bar{S} interaction vertex (cm); #Entries',4*140,0,140)
	l_h_lxy_eta.append(h_lxy_eta)
	h_vz_eta = TH1F('h_vz_eta'+str(lower_eta)+'_to_'+str(higher_eta),'h_vz_eta; vz #bar{S} interaction vertex (cm); #Entries',600,-300,300)
	l_h_vz_eta.append(h_vz_eta)
	h_lxy_n_loops = TH1F('h_lxy_n_loops'+str(i),'; l_{0} #bar{S} interaction vertex (cm); #Events/cm',140,0,140)
	h_eta_n_loops = TH1F('h_eta_n_loops'+str(i),'; #eta #bar{S} interaction vertex (cm); #Entries',90,-4.5,4.5)
	h_vz_n_loops = TH1F('h_vz_n_loops'+str(i),'; vz #bar{S} interaction vertex (cm); #Entries',600,-300,300)
	l_h_lxy_n_loops.append(h_lxy_n_loops)
	l_h_eta_n_loops.append(h_eta_n_loops)
	l_h_vz_n_loops.append(h_vz_n_loops)


	
#fill plots with for different eta ranges the lxy and vz of the interaction vertex
#fill plots for the first 100 loops, 2nd 100 loops , 3rd 100 loops, ... the lxy interaction vertex for eta within the barrel. Each plot should look alike. 
nEntries = tree.GetEntries()
print 'Number of entries in the tree: ', nEntries
for i in range(0,nEntries):
	if i > maxNEntries:
		break
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

#overlap the first five of l_h_lxy_n_loops
c_lxy_n_loops = TCanvas("lxy_n_loops","");
legend_c_lxy_n_loops = TLegend(0.6,0.85,0.99,0.99)
legend_c_lxy_n_loops_text = ["#bar{S} doing [0,100[ loops","#bar{S} doing [100,200[ loops","#bar{S} doing [200,300[ loops","#bar{S} doing [300,400[ loops"]
for i in range(0,4):
	h = l_h_lxy_n_loops[i]
	#if(h.GetSumw2N() == 0):
	#	h.Sumw2(kTRUE)
	#h.Scale(1./h.Integral(), "width");
	if i == 0:
		h.Draw("")
	else:
		h.Draw("same")
	h.SetLineColor(colours[i])
	h.SetFillColorAlpha(colours[i],0.9)
	#h.SetMarkerStyle(22+i)
	h.SetMarkerColor(colours[i])
	legend_c_lxy_n_loops.AddEntry(h,legend_c_lxy_n_loops_text[i],"f")
	h.SetStats(0)
legend_c_lxy_n_loops.Draw()
CMS_lumi.CMS_lumi(c_lxy_n_loops, 0, 11)
c_lxy_n_loops.SaveAs(plots_output_dir+c_lxy_n_loops.GetName()+".pdf")
c_lxy_n_loops.Write()
	

eta_n_loops_dir = fOut.mkdir("eta_n_loops")
eta_n_loops_dir.cd()
for h in l_h_eta_n_loops:
	h.Write()

vz_n_loops_dir = fOut.mkdir("vz_n_loops")
vz_n_loops_dir.cd()
for h in l_h_vz_n_loops:
	h.Write()



#some more plots for the looping mechanism and for the interaction vertex location of the antiS


overall_dir = fOut.mkdir("overall")
overall_dir.cd()
#plot the n_loops versus eta:
h2_n_loops_vs_eta = TH2F("h2_n_loops_vs_eta","; #eta(#bar{S}); #Loops; #Entries",20,-4,4,200,0,600)
#plot the number of antiS interacting after certain amount of loops for all, barrel and endcap
h_n_interactions_vs_n_loops = TH1F("h_n_interactions_vs_n_loops","; #Loops; #bar{S} interacting after #Loops",400,0,1600)
h_n_interactions_vs_n_loops_barrel = TH1F("h_n_interactions_vs_n_loops_barrel","; #Loops; #bar{S} interacting after #Loops",400,0,1600)
h_n_interactions_vs_n_loops_endcap = TH1F("h_n_interactions_vs_n_loops_endcap","; #Loops; #bar{S} interacting after #Loops",400,0,1600)
h_n_interactions_vs_n_loops_endcap_more_forward_than_endcap = TH1F("h_n_interactions_vs_n_loops_endcap_more_forward_than_endcap","; #Loops; #bar{S} interacting after #Loops",400,0,1600)
#2D plots of the interaction vertex location
h2_interaction_vertex_vx_vy = TH2F("h2_interaction_vertex_vx_vy","; v_{x} #bar{S} interaction vertex (cm); v_{y} #bar{S} interaction vertex (cm);Events/mm^{2}",2500,-125,125,2500,-125,125)
h2_interaction_vertex_vx_vy_zoom = TH2F("h2_interaction_vertex_vx_vy_zoom","; v_{x} #bar{S} interaction vertex (cm); v_{y} #bar{S} interaction vertex (cm);Events/4.10^{-4}mm^{2}",2500,-25,25,2500,-25,25)
h2_interaction_vertex_vz_lxy = TH2F("h2_interaction_vertex_vz_lxy","; v_{z} #bar{S} interaction vertex (cm); l_{0} #bar{S} interaction vertex (cm);Events/mm^{2}",2500,-125,125,1250,0,125)
h2_interaction_vertex_vz_lxy_zoom = TH2F("h2_interaction_vertex_vz_lxy_zoom","; v_{z} #bar{S} interaction vertex (cm); l_{0} #bar{S} interaction vertex (cm);Events/mm^{2}",2500,-125,125,200,0,20)

for i in range(0,nEntries):
	if i > maxNEntries:
		break
	tree.GetEntry(i)
	for j in range(0,len(tree._S_lxy_interaction_vertex)):
		h2_n_loops_vs_eta.Fill(tree._S_eta[j],tree._S_n_loops[j])
		h2_interaction_vertex_vx_vy.Fill(tree._S_vx_interaction_vertex[j],tree._S_vy_interaction_vertex[j])
		h2_interaction_vertex_vx_vy_zoom.Fill(tree._S_vx_interaction_vertex[j],tree._S_vy_interaction_vertex[j])
		h2_interaction_vertex_vz_lxy.Fill(tree._S_vz_interaction_vertex[j],tree._S_lxy_interaction_vertex[j])
		h2_interaction_vertex_vz_lxy_zoom.Fill(tree._S_vz_interaction_vertex[j],tree._S_lxy_interaction_vertex[j])
		h_n_interactions_vs_n_loops.Fill(tree._S_n_loops[j])
		if(abs(tree._S_eta[j])<1.1):
			h_n_interactions_vs_n_loops_barrel.Fill(tree._S_n_loops[j])
		elif(abs(tree._S_eta[j])>1.1 and abs(tree._S_eta[j])<2.5):
			h_n_interactions_vs_n_loops_endcap.Fill(tree._S_n_loops[j])
		else:
			h_n_interactions_vs_n_loops_endcap_more_forward_than_endcap.Fill(tree._S_n_loops[j])

l_TH1F = [h_n_interactions_vs_n_loops,h_n_interactions_vs_n_loops_barrel,h_n_interactions_vs_n_loops_endcap,h_n_interactions_vs_n_loops_endcap_more_forward_than_endcap]
for h in l_TH1F:
        h.SetDirectory(0)

c_n_interactions_vs_n_loops = TCanvas("c_n_interactions_vs_n_loops","");
legend_c_n_interactions_vs_n_loops = TLegend(0.7,0.85,0.99,0.99)
legend_text_c_n_interactions_vs_n_loops = ["All #bar{S}","#bar{S} in the barrel","#bar{S} in the endcap","#bar{S} |#eta| > 2.5"]
i_l_TH1F = 0
for h in l_TH1F:
	if(i_l_TH1F==0):
		h.Draw("L")
	else:	
		h.Draw("sameL")
	h.SetLineColor(colours[i_l_TH1F])
	h.SetMarkerStyle(22+i)
	h.SetMarkerColor(colours[i_l_TH1F])
	h.SetStats(0)
	legend_c_n_interactions_vs_n_loops.AddEntry(h,legend_text_c_n_interactions_vs_n_loops[i_l_TH1F],"lep")
	i_l_TH1F+=1

legend_c_n_interactions_vs_n_loops.Draw()
c_n_interactions_vs_n_loops.SetLogy()
CMS_lumi.CMS_lumi(c_n_interactions_vs_n_loops, 0, 11)
c_n_interactions_vs_n_loops.SaveAs(plots_output_dir+c_n_interactions_vs_n_loops.GetName()+".pdf")
c_n_interactions_vs_n_loops.Write()


l_TH2F = [h2_interaction_vertex_vz_lxy,h2_interaction_vertex_vz_lxy_zoom,h2_interaction_vertex_vx_vy,h2_interaction_vertex_vx_vy_zoom]
for h in l_TH2F:
	h.SetDirectory(0)
for h in l_TH2F:
	c= TCanvas(h.GetName(),"");
	c.SetRightMargin(0.2) #make room for the tile of the z scale
	if(h.GetSumw2N() == 0):
		h.Sumw2(kTRUE)
	h.Scale(1./h.Integral(), "width");
	h.Draw("colz")
	h.SetStats(0)
	CMS_lumi.CMS_lumi(c, 0, 11)
	c.SetLogz()
	c.SaveAs(plots_output_dir+h.GetName()+".pdf")
	c.Write()

antiS_properties_dir = fOut.mkdir("antiS_properties")

#properties of the antiS, so these are the ones which interact and go the correct granddaughters:
antiS_properties_dir.cd()
h_antiS_lxy = TH1F('h_antiS_lxy','; #bar{S} interaction vertex l_{0} (cm) ; Events/mm',1200,0,120)
h_antiS_vz = TH1F('h_antiS_vz','; #bar{S} interaction vertex |v_{z}| (cm); Events/cm',120,0,120)
h_neutron_momentum = TH1F('h_neutron_momentum','; p neutron (GeV); Events/0.01GeV',85,0,0.85)
h2_antiS_inv_mass_p = TH2F("h2_antiS_inv_mass_p","; mass #bar{S} (GeV); p #bar{S} (GeV);#Entries",100,-10,10,60,0,60)
h2_antiS_inv_mass_p_Ks_plus_Lambda = TH2F("h2_antiS_inv_mass_p_Ks_plus_Lambda","; m_{#bar{S},obs} (GeV); |#vec{p}_{K_{s}} + #vec{p}_{#bar{#Lambda}}| (GeV);Events/GeV^{2}",110,-5,6,60,0,40)
h_antiS_sumDaughters_openingsangle = TH1F('h_antiS_sumDaughters_openingsangle','; openings angle(#vec{p}_{K_{s}}+#vec{p}_{#bar{#Lambda}},#vec{p}_{#bar{S}}) (rad); Events/mrad',200,0,0.2)

#properties of the antiS daughters and granddaughters
momenta_daughters_and_grandaughters_dir = fOut.mkdir("momenta_daughters_and_grandaughters")
momenta_daughters_and_grandaughters_dir.cd()
h_pt_Ks = TH1F("h_pt_Ks",";p_{t} (GeV); #Entries",80,0,8)
h_pt_Ks_daug0 = TH1F("h_pt_Ks_daug0",";p_{t} (GeV); Events/0.1GeV",80,0,8)
h_pt_Ks_daug1 = TH1F("h_pt_Ks_daug1",";p_{t} (GeV); Events/0.1GeV",80,0,8)
h_pz_Ks = TH1F("h_pz_Ks",";p_{z} (GeV); Events/1GeV",80,0,80)
h_pz_Ks_daug0 = TH1F("h_pz_Ks_daug0",";p_{z} (GeV); Events/1GeV",80,0,80)
h_pz_Ks_daug1 = TH1F("h_pz_Ks_daug1",";p_{z} (GeV); Events/1GeV",80,0,80)
h2_pt_pz_Ks_daug0 = TH2F("h2_pt_pz_Ks_daug0",";p_{t} (GeV);|p_{z}| (GeV); Events/(1GeV*0.1GeV)",80,0,8,80,0,80)
h2_pt_pz_Ks_daug1 = TH2F("h2_pt_pz_Ks_daug1",";p_{t} (GeV);|p_{z}| (GeV); Events/(1GeV*0.1GeV)",80,0,8,80,0,80)

h_pt_AntiLambda = TH1F("h_pt_AntiLambda",";p_{t} (GeV); Events/0.1GeV",80,0,8)
h_pt_AntiLambda_AntiProton = TH1F("h_pt_AntiLambda_AntiProton",";p_{t} (GeV); Events/0.1GeV",80,0,8)
h_pt_AntiLambda_Pion = TH1F("h_pt_AntiLambda_Pion",";p_{t} (GeV); Events/0.1GeV",80,0,8)
h_pz_AntiLambda = TH1F("h_pz_AntiLambda",";p_{z} (GeV); Events/1GeV",80,0,80)
h_pz_AntiLambda_AntiProton = TH1F("h_pz_AntiLambda_AntiProton",";p_{z} (GeV); Events/1GeV",80,0,80)
h_pz_AntiLambda_Pion = TH1F("h_pz_AntiLambda_Pion",";p_{z} (GeV); Events/1GeV",80,0,80)
h2_pt_pz_AntiLambda_AntiProton = TH2F("h2_pt_pz_AntiLambda_AntiProton",";p_{t} (GeV);|p_{z}| (GeV); Events/(1GeV*0.1GeV)",80,0,8,80,0,80)
h2_pt_pz_AntiLambda_Pion = TH2F("h2_pt_pz_AntiLambda_Pion",";p_{t} (GeV);|p_{z}| (GeV); Events/(1GeV*0.1GeV)",80,0,8,80,0,80)

#PCAs of the granddughters of the antiS
PCA_granddaughters_dir = fOut.mkdir("PCA_granddaughters")
PCA_granddaughters_dir.cd()

h_dxy_Ks_daug0 = TH1F("h_dxy_Ks_daug0",";d_{0} (cm); Events/mm",200,-10,10)
h_dxy_Ks_daug1 = TH1F("h_dxy_Ks_daug1",";d_{0} (cm); Events/mm",200,-10,10)
h_dz_Ks_daug0 = TH1F("h_dz_Ks_daug0",";d_{z} (cm); Events/cm",200,-100,100)
h_dz_Ks_daug1 = TH1F("h_dz_Ks_daug1",";d_{z} (cm); Events/cm",200,-100,100)

h_dxy_AntiLambda_AntiProton = TH1F("h_dxy_AntiLambda_AntiProton",";d_{0} (cm); Events/mm",200,-10,10)
h_dxy_AntiLambda_Pion = TH1F("h_dxy_AntiLambda_Pion",";d_{0} (cm); Events/mm",200,-10,10)
h_dz_AntiLambda_AntiProton = TH1F("h_dz_AntiLambda_AntiProton",";d_{z} (cm); Events/cm",200,-100,100)
h_dz_AntiLambda_Pion = TH1F("h_dz_AntiLambda_Pion",";d_{z} (cm); Events/cm",200,-100,100)


#info on the decay vertices of the V0s
decay_vertex_dir = fOut.mkdir("decay_vertex_V0s")
decay_vertex_dir.cd()
h_lxy_creation_vertex_Ks_daughters = TH1F("h_lxy_creation_vertex_Ks_daughters",";l_{0} creation vertex (cm);Events/cm ",120,0,120)
h_vz_creation_vertex_Ks_daughters = TH1F("h_vz_creation_vertex_Ks_daughters",";v_{z} creation vertex (cm);Events/10cm",80,-400,400)
h_lxy_creation_vertex_AntiLambda_daughters = TH1F("h_lxy_creation_vertex_AntiLambda_daughters;Events/cm",";l_{0} creation vertex (cm)",120,0,120)
h_vz_creation_vertex_AntiLambda_daughters = TH1F("h_vz_creation_vertex_AntiLambda_daughters;Events/10cm",";v_{z} creation vertex (cm)",80,-400,400)

h2_vx_vy_creation_vertex_Ks_daughters = TH2F("h2_vx_vy_creation_vertex_Ks_daughters",";v_{x} decay vertex K_{S}^{0} (cm); v_{y} decay vertex K_{S}^{0} (cm);Events/cm^{2}",240,-120,120,240,-120,120)
h2_vz_lxy_creation_vertex_Ks_daughters = TH2F("h2_vz_lxy_creation_vertex_Ks_daughters",";v_{z} decay vertex K_{S}^{0} (cm);l_{0} decay vertex K_{S}^{0} (cm);Events/cm^{2}",600,-300,300,120,0,120)
h2_vx_vy_creation_vertex_AntiLambda_daughters = TH2F("h2_vx_vy_creation_vertex_AntiLambda_daughters",";v_{x} decay vertex #bar{#Lambda}^{0} (cm);v_{y} decayvertex #bar{#Lambda}^{0} (cm);Events/cm^{2}",240,-120,120,240,-120,120)
h2_vz_lxy_creation_vertex_AntiLambda_daughters = TH2F("h2_vz_lxy_creation_vertex_AntiLambda_daughters",";v_{z} decay vertex #bar{#Lambda}^{0} (cm);l_{0} decay vertex #bar{#Lambda}^{0} (cm);Events/cm^{2}",600,-300,300,120,0,120)

#count the number of tracker layers for a track produced at a certain location
#NOTE: IMPORTANT: THE numberOfTrackerLayers INFORMATION WILL ONLY BE AVAILABLE IN THE TREE IF THE TREE WAS CREATED BY RUNNING ON A SAMPLE WHICH CONTAINS TRACKINGPARTICLES, THESE TRACKINGPARTICLES ARE ONLY AVAILABLE AFTER THE RECONSTRUCTION, SO ACTUALLY IT IS NOT REALLY CORRECT THAT THESE PLOTS ARE HERE IN A SCRIPT WHICH IS ONLY SUPPOSED TO ANALYZE THE GENSIM STEP
prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt",";v_{x} decay vertex K_{S}^{0} (cm); v_{y} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",240,-120,120,240,-120,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",600,-300,300,120,0,120)
prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt",";v_{x} decay vertex #bar{#Lambda}^{0} (cm); v_{y} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",240,-120,120,240,-120,120)
prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); l_{0} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",600,-300,300,120,0,120)
prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt",";v_{x} decay vertex #bar{#Lambda}^{0} (cm); v_{y} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",240,-120,120,240,-120,120)
prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt = TProfile2D("prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); l_{0} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",600,-300,300,120,0,120)

prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt = TProfile2D("prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt",";v_{x} decay vertex K_{S}^{0} (cm); v_{y} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",240,-120,120,240,-120,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",600,-300,300,120,0,120)
prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt = TProfile2D("prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt",";v_{x} decay vertex #bar{#Lambda}^{0} (cm); v_{y} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",240,-120,120,240,-120,120)
prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); l_{0} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #bar{p} daughter",600,-300,300,120,0,120)
prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt = TProfile2D("prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt",";v_{x} decay vertex #bar{#Lambda}^{0} (cm); v_{y} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",240,-120,120,240,-120,120)
prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt = TProfile2D("prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt",";v_{z} decay vertex #bar{#Lambda}^{0} (cm); l_{0} decay vertex #bar{#Lambda}^{0} (cm);mean #tracker layers hit by the track from #pi^{+} daughter",600,-300,300,120,0,120)

#study a bit more in detail for one of the track collections (in the end it are all just tracks here, so if you constrain the kinematics it does not matter much any more rather these are tracks from Ks or antiLambda)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_lowPz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_lowPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_middlePz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_middlePz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_highPz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_highPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_lowPz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_lowPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_middlePz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_middlePz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_highPz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_highPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_lowPz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_lowPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_middlePz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_middlePz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_highPz = TProfile2D("prof2_lxy_vx_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_highPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)


for i in range(0,nEntries):

	if i > maxNEntries:
		break

	tree.GetEntry(i)


	h_antiS_lxy.Fill(tree._S_lxy_interaction_vertex[0])	
	h_antiS_vz.Fill(abs(tree._S_vz_interaction_vertex[0]))	
	h_neutron_momentum.Fill(abs(tree._n_p[0]))
	h_antiS_sumDaughters_openingsangle.Fill(tree._S_sumDaughters_openingsangle[0])	
	h2_antiS_inv_mass_p.Fill(tree._S_mass[0], np.sqrt( tree._S_pz[0]*tree._S_pz[0] + tree._S_pt[0]*tree._S_pt[0] ) )
	h2_antiS_inv_mass_p_Ks_plus_Lambda.Fill(tree._S_mass[0], np.sqrt( pow(tree._Ks_pz[0] + tree._Lambda_pz[0] , 2) + pow(tree._Ks_pt[0] + tree._Lambda_pt[0] , 2) ) )
	
	h_pt_Ks.Fill(tree._Ks_pt[0])
	h_pt_Ks_daug0.Fill(tree._GEN_Ks_daughter0_pt[0])
	h_pt_Ks_daug1.Fill(tree._GEN_Ks_daughter1_pt[0])
	h_pz_Ks.Fill(abs(tree._Ks_pz[0]))
	h_pz_Ks_daug0.Fill(abs(tree._GEN_Ks_daughter0_pz[0]))
	h_pz_Ks_daug1.Fill(abs(tree._GEN_Ks_daughter1_pz[0]))
	h2_pt_pz_Ks_daug0.Fill(tree._GEN_Ks_daughter0_pt[0], abs(tree._GEN_Ks_daughter0_pz[0]))
	h2_pt_pz_Ks_daug1.Fill(tree._GEN_Ks_daughter1_pt[0], abs(tree._GEN_Ks_daughter1_pz[0]))

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
	h2_pt_pz_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_pt[0],abs(tree._GEN_AntiLambda_AntiProton_pz[0]))
	h2_pt_pz_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_pt[0],abs(tree._GEN_AntiLambda_Pion_pz[0]))

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
	h2_vz_lxy_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0])
	h2_vx_vy_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_vy[0])
	h2_vz_lxy_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_lxy[0])

	#fill the tprofiles
	#the Ks daughters are both pions, so can put them in same tprofiles
	#low pt
	if(tree._GEN_Ks_daughter0_pt[0]>0.1 and tree._GEN_Ks_daughter0_pt[0]<0.9):
		prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter0_vx[0],tree._GEN_Ks_daughter0_vy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
	if(tree._GEN_Ks_daughter1_pt[0]>0.1 and tree._GEN_Ks_daughter1_pt[0]<0.9):
		prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter1_vx[0],tree._GEN_Ks_daughter1_vy[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt.Fill(tree._GEN_Ks_daughter1_vz[0],tree._GEN_Ks_daughter1_lxy[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])

	if(tree._GEN_AntiLambda_AntiProton_pt[0]>0.1 and tree._GEN_AntiLambda_AntiProton_pt[0]<0.9):
		prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_vy[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_lxy[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
	if(tree._GEN_AntiLambda_Pion_pt[0]>0.1 and tree._GEN_AntiLambda_Pion_pt[0]<0.9):
		prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_Pion_vx[0],tree._GEN_AntiLambda_Pion_vy[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_lowPt.Fill(tree._GEN_AntiLambda_Pion_vz[0],tree._GEN_AntiLambda_Pion_lxy[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])

	#high pt
	if(tree._GEN_Ks_daughter0_pt[0]>0.9):
		prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter0_vx[0],tree._GEN_Ks_daughter0_vy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])
	if(tree._GEN_Ks_daughter1_pt[0]>0.9):
		prof2_vx_vy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter1_vx[0],tree._GEN_Ks_daughter1_vy[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt.Fill(tree._GEN_Ks_daughter1_vz[0],tree._GEN_Ks_daughter1_lxy[0],tree._GEN_Ks_daughter1_numberOfTrackerLayers[0])
	
	if(tree._GEN_AntiLambda_AntiProton_pt[0]>0.9):
		prof2_vx_vy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_vy[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_AntiLambda_AntiProton_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_lxy[0],tree._GEN_AntiLambda_AntiProton_numberOfTrackerLayers[0])
	if(tree._GEN_AntiLambda_Pion_pt[0]>0.9):
		prof2_vx_vy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_Pion_vx[0],tree._GEN_AntiLambda_Pion_vy[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])
		prof2_vz_lxy_creation_vertex_AntiLambda_Pion_numberOfTrackerLayers_highPt.Fill(tree._GEN_AntiLambda_Pion_vz[0],tree._GEN_AntiLambda_Pion_lxy[0],tree._GEN_AntiLambda_Pion_numberOfTrackerLayers[0])


	#low pt Ks
	if(tree._GEN_Ks_daughter0_pt[0]>0.35 and tree._GEN_Ks_daughter0_pt[0]<0.5):
		if(abs(tree._GEN_Ks_daughter0_pz[0]) < 1):
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_lowPz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	
		elif(abs(tree._GEN_Ks_daughter0_pz[0]) > 1 and abs(tree._GEN_Ks_daughter0_pz[0]) < 3):
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_middlePz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	
		else:
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_lowPt_highPz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	
	elif(tree._GEN_Ks_daughter0_pt[0]>0.5 and tree._GEN_Ks_daughter0_pt[0]<1.):
		if(abs(tree._GEN_Ks_daughter0_pz[0]) < 1):
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_lowPz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	
		elif(abs(tree._GEN_Ks_daughter0_pz[0]) > 1 and abs(tree._GEN_Ks_daughter0_pz[0]) < 3):
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_middlePz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	
		else:
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_middlePt_highPz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	
	else:
		if(abs(tree._GEN_Ks_daughter0_pz[0]) < 1):
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_lowPz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	
		elif(abs(tree._GEN_Ks_daughter0_pz[0]) > 1 and abs(tree._GEN_Ks_daughter0_pz[0]) < 3):
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_middlePz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	
		else:
			prof2_vz_lxy_creation_vertex_Ks_daughters_numberOfTrackerLayers_highPt_highPz.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],tree._GEN_Ks_daughter0_numberOfTrackerLayers[0])	

	


momenta_daughters_and_grandaughters_dir.cd()

c_antiS_lxy = TCanvas("c_antiS_lxy","");
h_antiS_lxy.DrawNormalized()
c_antiS_lxy.Write()

c_antiS_vz = TCanvas("c_antiS_vz","");
h_antiS_vz.DrawNormalized()
c_antiS_vz.Write()

l_TH1F = [h_antiS_lxy,h_antiS_vz]
for h in l_TH1F:
        h.SetDirectory(0)
i_l_TH1F = 0
for h in l_TH1F:
	c= TCanvas(h.GetName(),"");
	if(h.GetSumw2N() == 0):
		h.Sumw2(kTRUE)
	#h.Scale(1./h.Integral(), "width");
	h.Draw("CL")
	h.SetStats(0)
	CMS_lumi.CMS_lumi(c, 0, 11)
	if i_l_TH1F == 0:
		c.SetLogy()
	c.SaveAs(plots_output_dir+h.GetName()+".pdf")
	c.Write()
	i_l_TH1F+=1

h2_antiS_inv_mass_p.Write()
h2_antiS_inv_mass_p_Ks_plus_Lambda.Write()


c_h_neutron_momentum = TCanvas("c_"+h_neutron_momentum.GetName(),"");
h_neutron_momentum.Draw("PCE1")
if(h_neutron_momentum.GetSumw2N() == 0):
	h_neutron_momentum.Sumw2(kTRUE)
h_neutron_momentum.Scale(1./h_neutron_momentum.Integral(), "width");
h_neutron_momentum.SetLineColor(colours[0])
h_neutron_momentum.SetMarkerStyle(22)
h_neutron_momentum.SetMarkerColor(colours[0])
h_neutron_momentum.SetStats(0)
CMS_lumi.CMS_lumi(c_h_neutron_momentum, 0, 11)
c_h_neutron_momentum.SetLogy()
c_h_neutron_momentum.SaveAs(plots_output_dir+c_h_neutron_momentum.GetName()+".pdf")
c_h_neutron_momentum.Write()


h_antiS_sumDaughters_openingsangle
c_h_antiS_sumDaughters_openingsangle = TCanvas("c_"+h_antiS_sumDaughters_openingsangle.GetName(),"");
h_antiS_sumDaughters_openingsangle.Draw("PCE1")
if(h_antiS_sumDaughters_openingsangle.GetSumw2N() == 0):
	h_antiS_sumDaughters_openingsangle.Sumw2(kTRUE)
h_antiS_sumDaughters_openingsangle.Scale(1./h_antiS_sumDaughters_openingsangle.Integral(), "width");
h_antiS_sumDaughters_openingsangle.SetLineColor(colours[0])
h_antiS_sumDaughters_openingsangle.SetMarkerStyle(22)
h_antiS_sumDaughters_openingsangle.SetMarkerColor(colours[0])
h_antiS_sumDaughters_openingsangle.SetStats(0)
CMS_lumi.CMS_lumi(c_h_antiS_sumDaughters_openingsangle, 0, 11)
c_h_antiS_sumDaughters_openingsangle.SetLogy()
c_h_antiS_sumDaughters_openingsangle.SaveAs(plots_output_dir+c_h_antiS_sumDaughters_openingsangle.GetName()+".pdf")
c_h_antiS_sumDaughters_openingsangle.Write()


c_h2_antiS_inv_mass_p_Ks_plus_Lambda= TCanvas(h2_antiS_inv_mass_p_Ks_plus_Lambda.GetName(),"");
c_h2_antiS_inv_mass_p_Ks_plus_Lambda.SetRightMargin(0.2) #make room for the tile of the z scale
if(h2_antiS_inv_mass_p_Ks_plus_Lambda.GetSumw2N() == 0):
	h2_antiS_inv_mass_p_Ks_plus_Lambda.Sumw2(kTRUE)
h2_antiS_inv_mass_p_Ks_plus_Lambda.Scale(1./h2_antiS_inv_mass_p_Ks_plus_Lambda.Integral(), "width");
h2_antiS_inv_mass_p_Ks_plus_Lambda.Draw("SURF1")
h2_antiS_inv_mass_p_Ks_plus_Lambda.SetStats(0)
h2_antiS_inv_mass_p_Ks_plus_Lambda.GetXaxis().SetTitleOffset(1.5)
h2_antiS_inv_mass_p_Ks_plus_Lambda.GetYaxis().SetTitleOffset(2)
c_h2_antiS_inv_mass_p_Ks_plus_Lambda.SetLogz()
CMS_lumi.CMS_lumi(c, 0, 11)
c_h2_antiS_inv_mass_p_Ks_plus_Lambda.SaveAs(plots_output_dir+h2_antiS_inv_mass_p_Ks_plus_Lambda.GetName()+".pdf")
c_h2_antiS_inv_mass_p_Ks_plus_Lambda.Write()


ll_TH1F = [
[h_pt_Ks_daug0,h_pt_Ks],
[h_pz_Ks_daug0,h_pz_Ks],
[h_pz_AntiLambda_Pion,h_pz_AntiLambda_AntiProton,h_pz_AntiLambda],
[h_pt_AntiLambda_Pion,h_pt_AntiLambda_AntiProton,h_pt_AntiLambda]
]

ll_legend_text =  [
["K_{S}^{0} daughters","K_{S}^{0}"],
["K_{S}^{0} daughters","K_{S}^{0}"],
["#bar{#Lambda}-#pi^{+}","#bar{#Lambda}-#bar{p}","#bar{#Lambda}"],
["#bar{#Lambda}-#pi^{+}","#bar{#Lambda}-#bar{p}","#bar{#Lambda}"]
]
for l in ll_TH1F:
	for h in l:
		h.SetDirectory(0)

i_ll_TH1F = 0
for l in ll_TH1F:
	legend = TLegend(0.6,0.85,0.99,0.99)
	c = TCanvas("c_"+l[0].GetName(),"");
	i_l_TH1F = 0
	for h in l:
		if(i_l_TH1F==0):
			h.Draw("PCE1")
		else:	
			h.Draw("PCE1same")
		if(h.GetSumw2N() == 0):
        	       h.Sumw2(kTRUE)
        	h.Scale(1./h.Integral(), "width");
		h.SetLineColor(colours[i_l_TH1F])
		h.SetMarkerStyle(22+i_l_TH1F)
		h.SetMarkerColor(colours[i_l_TH1F])
		h.SetStats(0)
		legend.AddEntry(h,ll_legend_text[i_ll_TH1F][i_l_TH1F],"lep")
		i_l_TH1F+=1

	legend.Draw()
	CMS_lumi.CMS_lumi(c, 0, 11)
	c.SaveAs(plots_output_dir+c.GetName()+".pdf")
	c.Write()
	i_ll_TH1F +=1

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

ll_TH1F = [
[h_lxy_creation_vertex_Ks_daughters,h_lxy_creation_vertex_AntiLambda_daughters],
[h_vz_creation_vertex_Ks_daughters,h_vz_creation_vertex_AntiLambda_daughters],
[h_dxy_AntiLambda_AntiProton,h_dxy_AntiLambda_Pion,h_dxy_Ks_daug0],
[h_dz_AntiLambda_AntiProton,h_dz_AntiLambda_Pion,h_dz_Ks_daug0]
]

ll_legend_text =  [
["K_{S}^{0} daughters","#bar{#Lambda} daughters"],
["K_{S}^{0} daughters","#bar{#Lambda} daughters"],
["d_{0} #bar{#Lambda}-#bar{p}","d_{0} #bar{#Lambda}-#pi^+","d_{0} K_{S} daughters"],
["d_{z} #bar{#Lambda}-#bar{p}","d_{z} #bar{#Lambda}-#pi^+","d_{z} K_{S} daughters"],
]
for l in ll_TH1F:
	for h in l:
		h.SetDirectory(0)

i_ll_TH1F = 0
for l in ll_TH1F:
	legend = TLegend(0.6,0.85,0.99,0.99)
	c = TCanvas("c_"+l[0].GetName(),"");
	i_l_TH1F = 0
	for h in l:
		if(i_l_TH1F==0):
			h.Draw("PCE1")
		else:	
			h.Draw("PCE1same")
		if(h.GetSumw2N() == 0):
        	       h.Sumw2(kTRUE)
        	h.Scale(1./h.Integral(), "width");
		h.SetLineColor(colours[i_l_TH1F])
		h.SetMarkerStyle(22+i_l_TH1F)
		h.SetMarkerColor(colours[i_l_TH1F])
		h.SetStats(0)
		legend.AddEntry(h,ll_legend_text[i_ll_TH1F][i_l_TH1F],"lep")
		i_l_TH1F+=1

	legend.Draw()
	CMS_lumi.CMS_lumi(c, 0, 11)
	c.SaveAs(plots_output_dir+c.GetName()+".pdf")
	c.Write()
	i_ll_TH1F +=1


l_TH2F = [h2_vx_vy_creation_vertex_Ks_daughters,h2_vz_lxy_creation_vertex_Ks_daughters,h2_vx_vy_creation_vertex_AntiLambda_daughters,h2_vz_lxy_creation_vertex_AntiLambda_daughters]
for h in l_TH2F:
        h.SetDirectory(0)
for h in l_TH2F:
        c= TCanvas(h.GetName(),"");
        c.SetRightMargin(0.2) #make room for the tile of the z scale
        if(h.GetSumw2N() == 0):
                h.Sumw2(kTRUE)
        h.Scale(1./h.Integral(), "width");
        h.Draw("colz")
        h.SetStats(0)
	c.SetLogz()
        CMS_lumi.CMS_lumi(c, 0, 11)
        c.SaveAs(plots_output_dir+h.GetName()+".pdf")
        c.Write()


fOut.Write()
fOut.Close()
