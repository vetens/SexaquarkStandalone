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

colours = [1,2,4,35,38,41]

maxNEntries = 1e99

plots_output_dir = "plots_GENSIM/"

inputFile = '/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerGENSIM/test_FlatTreeGENSIM_Skimmed_trial25_1p8GeV_27102019_v1.root'
fIn = TFile(inputFile,'read')
tree = fIn.Get('FlatTreeProducerGENSIM/FlatTreeGENLevel') 
treeAllAntiS = fIn.Get('FlatTreeProducerGENSIM/FlatTreeGENLevelAllAntiS') 

fOut = TFile('macro_'+inputFile.rsplit('/', 1)[-1],'RECREATE')

all_antiS_dir = fOut.mkdir("all_antiS")
all_antiS_dir.cd()

h_eta_all_AntiS_non_weighted = TH1F('h_eta_all_AntiS_non_weighted','; #eta #bar{S}; #Entries/0.1#eta',160,-8,8)
h_eta_all_AntiS = TH1F('h_eta_all_AntiS','; #eta #bar{S}; #Entries/0.1#eta',160,-8,8)
h_vz_interaction_all_AntiS_non_weighted = TH1F('h_vz_interaction_all_AntiS_non_weighted','; absolute v_{z} interaction vertex #bar{S}; #Entries/5cm',100,-250,250)
h_vz_interaction_all_AntiS = TH1F('h_vz_interaction_all_AntiS','; absolute v_{z} interaction vertex #bar{S}; #Entries/5cm',100,-250,250)
h_vz_interaction_all_AntiS_zoom = TH1F('h_vz_interaction_all_AntiS_zoom','; absolute v_{z} interaction vertex #bar{S}; #Entries/0.5cm',100,-25,25)
tprof_eta_weighting_factor = TProfile('tprof_eta_weighting_factor','; #eta #bar{S}; Event weighting factor/0.1#eta',90,-4.5,4.5)
tprof_vz_weighting_factor = TProfile('tprof_vz_weighting_factor','; absolute v_{Z} #bar{S} interaction vertex; Event weighting factor/5cm',100,-250,250)
tprof_vz_weighting_factor_zoom = TProfile('tprof_vz_weighting_factor_zoom','; absolute v_{Z} #bar{S} interaction vertex; Event weighting factor/2mm',100,-10,10)
h_pt_all_AntiS = TH1F('h_pt_all_AntiS','; p_{t} #bar{S}; #Entries/0.1GeV',100,0,10)
h_pz_all_AntiS = TH1F('h_pz_all_AntiS','; |p_{z}| #bar{S}; #Entries/1GeV',80,0,80)
h_vz_creation_vertex_all_AntiS = TH1F('h_vz_creation_vertex_all_AntiS','; v_{z} creation vertex #bar{S}; #Entries/cm',60,-30,30)
h_nPV_all_AntiS = TH1F('h_nPV_all_AntiS','; #PV; #Entries',60,0,60)


#investigate the reconstructability: have to use histograms for nominator and denominator, because teff does not support weights and using tprofile the bins with y = 0 are not displayed
h_nom_vz_antiS_reconstructable = TH1F('h_nom_vz_antiS_reconstructable','; absolute v_{z} #bar{S} interaction vertex  (cm); Reconstructability/5cm',60,-150,150)
h_nom_eta_antiS_reconstructable = TH1F('h_nom_eta_antiS_reconstructable','; #eta #bar{S} interaction vertex  (cm); Reconstructability/0.1#eta',160,-8,8)
h_nom_pt_antiS_reconstructable = TH1F('h_nom_pt_antiS_reconstructable','; p_{t} #bar{S} (cm); Reconstructability/0.2GeV',50,0,10)
h_nom_pz_antiS_reconstructable = TH1F('h_nom_pz_antiS_reconstructable','; p_{z} #bar{S} (cm); Reconstructability/1GeV',80,0,80)

h_denom_vz_antiS_reconstructable = TH1F('h_denom_vz_antiS_reconstructable','; absolute v_{z} #bar{S} interaction vertex  (cm); Reconstructability/5cm',60,-150,150)
h_denom_eta_antiS_reconstructable = TH1F('h_denom_eta_antiS_reconstructable','; #eta #bar{S} interaction vertex  (cm); Reconstructability/0.1#eta',160,-8,8)
h_denom_pt_antiS_reconstructable = TH1F('h_denom_pt_antiS_reconstructable','; p_{t} #bar{S} (cm); Reconstructability/0.2GeV',50,0,10)
h_denom_pz_antiS_reconstructable = TH1F('h_denom_pz_antiS_reconstructable','; p_{z} #bar{S} (cm); Reconstructability/1GeV',80,0,80)


#first do this small loop which runs over the tree containing all the antiS (i.e. also the ones which do not go the correct granddaughters)
nAntiSReconstructable = 0.
nAntiSTotal = 0.
for i in range(0,treeAllAntiS.GetEntries()):
	treeAllAntiS.GetEntry(i)

	#the AntiS does not necessarily have two daughters, so I cannot get the vz through that for all antiS, so have to calculate the vz from eta assuming the beampipe is infinitely thin and at a radius of 2.21cm 
	vz_interaction_antiS = 2.21/np.tan( 2*np.arctan( np.exp(-treeAllAntiS._S_eta_all[0]) ) ) 
	weight_factor = treeAllAntiS._S_event_weighting_factor_all[0]*treeAllAntiS._S_event_weighting_factor_PU_all[0]

	h_eta_all_AntiS_non_weighted.Fill(treeAllAntiS._S_eta_all[0],1)
	h_eta_all_AntiS.Fill(treeAllAntiS._S_eta_all[0],weight_factor)
	tprof_eta_weighting_factor.Fill(treeAllAntiS._S_eta_all[0],weight_factor)
	tprof_vz_weighting_factor.Fill(vz_interaction_antiS,weight_factor)
	tprof_vz_weighting_factor_zoom.Fill(vz_interaction_antiS,weight_factor)

	h_vz_interaction_all_AntiS_non_weighted.Fill(vz_interaction_antiS,1)
	h_vz_interaction_all_AntiS.Fill(vz_interaction_antiS,weight_factor)
	h_vz_interaction_all_AntiS_zoom.Fill(vz_interaction_antiS,weight_factor)

	h_pt_all_AntiS.Fill(treeAllAntiS._S_pt_all[0],weight_factor)
	h_pz_all_AntiS.Fill(abs(treeAllAntiS._S_pz_all[0]),weight_factor)

	h_vz_creation_vertex_all_AntiS.Fill(treeAllAntiS._S_vz_creation_vertex_all[0],weight_factor)
	h_nPV_all_AntiS.Fill(treeAllAntiS._S_nGoodPV_all[0],weight_factor)

	if(treeAllAntiS._S_reconstructable_all[0] == 1):
		nAntiSReconstructable += weight_factor
		h_nom_vz_antiS_reconstructable.Fill(vz_interaction_antiS,weight_factor)
		h_nom_eta_antiS_reconstructable.Fill(treeAllAntiS._S_eta_all[0],weight_factor)
		h_nom_pt_antiS_reconstructable.Fill(treeAllAntiS._S_pt_all[0],weight_factor)
		h_nom_pz_antiS_reconstructable.Fill(treeAllAntiS._S_pz_all[0],weight_factor)

	h_denom_vz_antiS_reconstructable.Fill(vz_interaction_antiS,weight_factor)
	h_denom_eta_antiS_reconstructable.Fill(treeAllAntiS._S_eta_all[0],weight_factor)
	h_denom_pt_antiS_reconstructable.Fill(treeAllAntiS._S_pt_all[0],weight_factor)
	h_denom_pz_antiS_reconstructable.Fill(treeAllAntiS._S_pz_all[0],weight_factor)
	nAntiSTotal += weight_factor

print 'the overall antiS reconstructability: ', nAntiSReconstructable, '/', nAntiSTotal, '=', float(nAntiSReconstructable)/float(nAntiSTotal)

l_nom = [h_nom_vz_antiS_reconstructable,h_nom_eta_antiS_reconstructable,h_nom_pt_antiS_reconstructable,h_nom_pz_antiS_reconstructable]
l_denom = [h_denom_vz_antiS_reconstructable,h_denom_eta_antiS_reconstructable,h_denom_pt_antiS_reconstructable,h_denom_pz_antiS_reconstructable]
for i in range(0,len(l_nom)):
	teff = TEfficiency(l_nom[i],l_denom[i])
	teff.SetName('teff_'+l_nom[i].GetName())
	c = TCanvas("c_"+teff.GetName(),"");
	teff.Draw("")
	CMS_lumi.CMS_lumi(c, 0, 11)
	c.SaveAs(plots_output_dir+c.GetName()+".pdf")
	c.Write()
	teff.Write()

	
l_tprof  = [h_eta_all_AntiS,h_vz_interaction_all_AntiS,h_vz_interaction_all_AntiS_zoom,h_pt_all_AntiS,h_pz_all_AntiS,tprof_eta_weighting_factor,tprof_vz_weighting_factor,tprof_vz_weighting_factor_zoom]
for h in l_tprof:
        c = TCanvas("c_"+h.GetName(),"");
        h.Draw("PCE1")
        h.SetLineColor(colours[0])
        h.SetMarkerStyle(22)
        h.SetMarkerColor(colours[0])
        h.SetStats(0)
        CMS_lumi.CMS_lumi(c, 0, 11)
        #c_h_antiS_sumDaughters_openingsangle.SetLogy()
        c.SaveAs(plots_output_dir+c.GetName()+".pdf")
        c.Write()	

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
	h_lxy_eta = TH1F('h_lxy_eta'+str(lower_eta)+'_to_'+str(higher_eta),'h_lxy_eta; l_{0} #bar{S} interaction vertex (cm); #Entries/0.1mm',300,0,3)
	l_h_lxy_eta.append(h_lxy_eta)
	h_vz_eta = TH1F('h_vz_eta'+str(lower_eta)+'_to_'+str(higher_eta),'h_vz_eta; vz #bar{S} interaction vertex (cm); #Entries',600,-300,300)
	l_h_vz_eta.append(h_vz_eta)
	h_lxy_n_loops = TH1F('h_lxy_n_loops'+str(i),'; l_{0} #bar{S} interaction vertex (cm); #Events/0.1mm',10,2.16,2.26)
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
	eta_range = int(abs(tree._S_eta[0])/eta_unit)
	l_h_lxy_eta[eta_range].Fill(tree._S_lxy_interaction_vertex[0])
	l_h_vz_eta[eta_range].Fill(tree._S_vz_interaction_vertex[0])
	if(int(tree._S_n_loops[0])<400):
		if(abs(tree._S_eta[0])<1.1):
			l_h_lxy_n_loops[int(tree._S_n_loops[0]/100)].Fill(tree._S_lxy_interaction_vertex[0])
		else:
			l_h_lxy_n_loops[10].Fill(tree._S_lxy_interaction_vertex[0])
		l_h_eta_n_loops[int(tree._S_n_loops[0])].Fill(tree._S_eta[0])
		l_h_vz_n_loops[int(tree._S_n_loops[0])].Fill(tree._S_vz_interaction_vertex[0])
	else:
		l_h_lxy_n_loops[10].Fill(tree._S_lxy_interaction_vertex[0])


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
#2D plots of the interaction vertex location--> do not weigh these yet, it is just to show that the interaction probability is flat. 
h_interaction_vertex_lxy_unweighted = TH1F("h_interaction_vertex_lxy_unweighted","; l_{xy} #bar{S} interaction vertex (cm); Events/0.1mm",10,2.16,2.26)
h_interaction_vertex_lxy_absolute_unweighted = TH1F("h_interaction_vertex_lxy_absolute_unweighted","; absolute l_{xy} #bar{S} interaction vertex (cm); Events/0.1mm",10,2.16,2.26)
h2_interaction_vertex_vx_vy_unweighted = TH2F("h2_interaction_vertex_vx_vy_unweighted","; v_{x} #bar{S} interaction vertex (cm); v_{y} #bar{S} interaction vertex (cm);Events/(0.1mm x 0.1mm)",2000,-10,10,2000,-10,10)
h2_interaction_vertex_vx_vy_zoom_unweighted = TH2F("h2_interaction_vertex_vx_vy_zoom_unweighted","; absolute v_{x} #bar{S} interaction vertex (cm); absolute v_{y} #bar{S} interaction vertex (cm);Events/(0.1mm x 0.1mm)",600,-3,3,600,-3,3)
h2_interaction_vertex_vz_lxy_unweighted = TH2F("h2_interaction_vertex_vz_lxy_unweighted","; v_{z} #bar{S} interaction vertex (cm); l_{0} #bar{S} interaction vertex (cm);Events/(10cm x 0.1mm)",34,-170,170,1000,0,10)
h2_interaction_vertex_vz_lxy_zoom_unweighted = TH2F("h2_interaction_vertex_vz_lxy_zoom_unweighted","; v_{z} #bar{S} interaction vertex (cm); l_{0} #bar{S} interaction vertex (cm);Events/(10cm x 0.1mm)",34,-170,170,10,2.16,2.26)

for i in range(0,tree.GetEntries()):
	if i > maxNEntries:
		break
	tree.GetEntry(i)
	h2_n_loops_vs_eta.Fill(tree._S_eta[0],tree._S_n_loops[0])
	h_interaction_vertex_lxy_unweighted.Fill(tree._S_lxy_interaction_vertex[0])
	h_interaction_vertex_lxy_absolute_unweighted.Fill( np.sqrt( np.power(tree._S_vx_interaction_vertex[0],2) + np.power(tree._S_vy_interaction_vertex[0],2) ) )
	h2_interaction_vertex_vx_vy_unweighted.Fill(tree._S_vx_interaction_vertex[0],tree._S_vy_interaction_vertex[0])
	h2_interaction_vertex_vx_vy_zoom_unweighted.Fill(tree._S_vx_interaction_vertex[0],tree._S_vy_interaction_vertex[0])
	h2_interaction_vertex_vz_lxy_unweighted.Fill(tree._S_vz_interaction_vertex[0],tree._S_lxy_interaction_vertex[0])
	h2_interaction_vertex_vz_lxy_zoom_unweighted.Fill(tree._S_vz_interaction_vertex[0],tree._S_lxy_interaction_vertex[0])
	h_n_interactions_vs_n_loops.Fill(tree._S_n_loops[0])
	if(abs(tree._S_eta[0])<1.1):
		h_n_interactions_vs_n_loops_barrel.Fill(tree._S_n_loops[0])
	elif(abs(tree._S_eta[0])>1.1 and abs(tree._S_eta[0])<2.5):
		h_n_interactions_vs_n_loops_endcap.Fill(tree._S_n_loops[0])
	else:
		h_n_interactions_vs_n_loops_endcap_more_forward_than_endcap.Fill(tree._S_n_loops[0])

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

c_interaction_vertex_lxy_unweighted = TCanvas("c_interaction_vertex_lxy_absolute_unweighted","");
h_interaction_vertex_lxy_absolute_unweighted.SetLineColor(colours[0])
h_interaction_vertex_lxy_absolute_unweighted.SetMarkerStyle(22)
h_interaction_vertex_lxy_absolute_unweighted.SetMarkerColor(colours[0])
h_interaction_vertex_lxy_absolute_unweighted.SetStats(0)
h_interaction_vertex_lxy_absolute_unweighted.Draw("PCE")
c_interaction_vertex_lxy_unweighted.Write()
c_interaction_vertex_lxy_unweighted.SaveAs(plots_output_dir+c_interaction_vertex_lxy_unweighted.GetName()+'.pdf')

l_TH2F = [h2_interaction_vertex_vz_lxy_unweighted,h2_interaction_vertex_vz_lxy_zoom_unweighted,h2_interaction_vertex_vx_vy_unweighted,h2_interaction_vertex_vx_vy_zoom_unweighted]
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
	#c.SetLogz()
	c.SaveAs(plots_output_dir+h.GetName()+".pdf")
	c.Write()

#enough looping validation, now do the actual parameters of the S and antiS. So these should all be weighted.
antiS_properties_dir = fOut.mkdir("antiS_properties")

antiS_properties_dir.cd()
#interaction vertex 2D plots, this time weighted
h2_interaction_vertex_vx_vy= TH2F("h2_interaction_vertex_vx_vy","; v_{x} #bar{S} interaction vertex (cm); v_{y} #bar{S} interaction vertex (cm);Events/(0.1mm x 0.1mm)",2000,-10,10,2000,-10,10)
h2_interaction_vertex_vx_vy_zoom= TH2F("h2_interaction_vertex_vx_vy_zoom","; v_{x} #bar{S} interaction vertex (cm); v_{y} #bar{S} interaction vertex (cm);Events/(0.1mm x 0.1mm)",600,-3,3,600,-3,3)
h2_interaction_vertex_vz_lxy= TH2F("h2_interaction_vertex_vz_lxy","; v_{z} #bar{S} interaction vertex (cm); l_{0} #bar{S} interaction vertex (cm);Events/(10cm x 0.1mm)",34,-170,170,1000,0,10)
h2_interaction_vertex_vz_lxy_zoom= TH2F("h2_interaction_vertex_vz_lxy_zoom","; v_{z} #bar{S} interaction vertex (cm); l_{0} #bar{S} interaction vertex (cm);Events/(10cm x 0.1mm)",34,-170,170,10,2.16,2.26)

#properties of the antiS, so these are the ones which interact and go the correct granddaughters:
h_antiS_eta = TH1F('h_antiS_eta','; #bar{S} #eta ; Events/mm',160,-8,8)
h_antiS_lxy_creation_vertex = TH1F('h_antiS_lxy_creation_vertex','; #bar{S} cv l_{0} (cm) ; Events/0.1mm',2000,-1,1)
h_antiS_vz_creation_vertex = TH1F('h_antiS_vz_creation_vertex','; #bar{S} cv absolute v_{z} (cm); Events/cm',60,-30,30)
h_antiS_lxy_interaction_vertex = TH1F('h_antiS_lxy_interaction_vertex','; l_{0,bs} interaction vertex #bar{S} (cm) ; Events/0.1mm',10,2.16,2.26)
h_antiS_vz_interaction_vertex = TH1F('h_antiS_vz_interaction_vertex','; absolute v_{z}  interaction vertex #bar{S} (cm); Events/5cm',68,-170,170)
h_neutron_momentum = TH1F('h_neutron_momentum','; p neutron (GeV); Events/0.01GeV',85,0,0.85)
h2_antiS_inv_mass_p = TH2F("h2_antiS_inv_mass_p","; mass #bar{S} (GeV); p #bar{S} (GeV);#Entries",100,-10,10,60,0,60)
h2_antiS_inv_mass_p_Ks_plus_Lambda = TH2F("h2_antiS_inv_mass_p_Ks_plus_Lambda","; m_{#bar{S},obs} (GeV); |#vec{p}_{K_{s}} + #vec{p}_{#bar{#Lambda}^{0}}| (GeV);Events/(0.1GeVx1GeV)",110,-5,6,40,0,40)
h_antiS_sumDaughters_openingsangle = TH1F('h_antiS_sumDaughters_openingsangle','; openings angle(#vec{p}_{K_{s}}+#vec{p}_{#bar{#Lambda}^{0}},#vec{p}_{#bar{S}}) (rad); Events/mrad',200,0,0.2)
h_antiS_sumDaughters_deltaPhi = TH1F('h_antiS_sumDaughters_deltaPhi','; #Delta#phi(#vec{p}_{K_{s}}+#vec{p}_{#bar{#Lambda}^{0}},#vec{p}_{#bar{S}}) (rad); Events/10 mrad',30,0,0.3)
h_antiS_sumDaughters_deltaEta = TH1F('h_antiS_sumDaughters_deltaEta','; #Delta#eta(#vec{p}_{K_{s}}+#vec{p}_{#bar{#Lambda}^{0}},#vec{p}_{#bar{S}}); Events/0.01#eta',30,0,0.3)
h_antiS_sumDaughters_deltaR = TH1F('h_antiS_sumDaughters_deltaR','; #DeltaR(#vec{p}_{K_{s}}+#vec{p}_{#bar{#Lambda}^{0}},#vec{p}_{#bar{S}}); Events/0.01#DeltaR',30,0,0.3)

#properties of the antiS daughters and granddaughters
momenta_daughters_and_grandaughters_dir = fOut.mkdir("momenta_daughters_and_grandaughters")
momenta_daughters_and_grandaughters_dir.cd()
h_pt_Ks = TH1F("h_pt_Ks",";p_{t} (GeV); #Entries",60,0,6)
h_p_Ks = TH1F("h_p_Ks",";p (GeV); #Entries",60,0,6)
h_pt_Ks_daug0 = TH1F("h_pt_Ks_daug0",";p_{t} (GeV); Events/0.1GeV",60,0,6)
h_pt_Ks_daug1 = TH1F("h_pt_Ks_daug1",";p_{t} (GeV); Events/0.1GeV",60,0,6)
h_p_Ks_daug0 = TH1F("h_p_Ks_daug0",";p (GeV); Events/0.1GeV",60,0,6)
h_p_Ks_daug1 = TH1F("h_p_Ks_daug1",";p (GeV); Events/0.1GeV",60,0,6)
h_pz_Ks = TH1F("h_pz_Ks",";p_{z} (GeV); Events/1GeV",30,0,30)
h_pz_Ks_daug0 = TH1F("h_pz_Ks_daug0",";p_{z} (GeV); Events/1GeV",30,0,30)
h_pz_Ks_daug1 = TH1F("h_pz_Ks_daug1",";p_{z} (GeV); Events/1GeV",30,0,30)
h2_pt_pz_Ks_daug0 = TH2F("h2_pt_pz_Ks_daug0",";p_{t} (GeV);|p_{z}| (GeV); Events/(1GeV*0.1GeV)",60,0,6,30,0,30)
h2_pt_pz_Ks_daug1 = TH2F("h2_pt_pz_Ks_daug1",";p_{t} (GeV);|p_{z}| (GeV); Events/(1GeV*0.1GeV)",60,0,6,30,0,30)
h2_eta_pt_Ks_daug0 = TH2F("h2_eta_pt_Ks_daug0",";#eta K_{s} daughter;p_{t} K_{s} daughter (GeV); Events/(0.1#eta*0.1GeV)",100,-5,5,60,0,6)
h2_eta_pt_Ks_daug1 = TH2F("h2_eta_pt_Ks_daug1",";#eta K_{s} daughter;p_{t} K_{s} daughter (GeV);  Events/(0.1#eta*0.1GeV)",100,-5,5,60,0,6)
h2_eta_pz_Ks_daug0 = TH2F("h2_eta_pz_Ks_daug0",";#eta K_{s} daughter;|p_{z}| K_{s} daughter (GeV); Events/(0.1#eta*1GeV)",100,-5,5,60,0,60)
h2_eta_pz_Ks_daug1 = TH2F("h2_eta_pz_Ks_daug1",";#eta K_{s} daughter;|p_{z}| K_{s} daughter (GeV);  Events/(0.1#eta*1GeV)",100,-5,5,60,0,60)
h2_eta_p_Ks_daug0 = TH2F("h2_eta_p_Ks_daug0",";#eta K_{s} daughter;p K_{s} daughter (GeV); Events/(0.1#eta*0.1GeV)",100,-5,5,600,0,60)
h2_eta_p_Ks_daug1 = TH2F("h2_eta_p_Ks_daug1",";#eta K_{s} daughter;p K_{s} daughter (GeV);  Events/(0.1#eta*0.1GeV)",100,-5,5,600,0,60)


h_pt_AntiLambda = TH1F("h_pt_AntiLambda",";p_{t} (GeV); Events/0.1GeV",60,0,6)
h_p_AntiLambda = TH1F("h_p_AntiLambda",";p (GeV); Events/0.001GeV",6000,0,6)
h_pt_AntiLambda_AntiProton = TH1F("h_pt_AntiLambda_AntiProton",";p_{t} (GeV); Events/0.1GeV",60,0,6)
h_pt_AntiLambda_Pion = TH1F("h_pt_AntiLambda_Pion",";p_{t} (GeV); Events/0.1GeV",60,0,6)
h_p_AntiLambda_AntiProton = TH1F("h_p_AntiLambda_AntiProton",";} (GeV); Events/0.01GeV",6000,0,6)
h_p_AntiLambda_Pion = TH1F("h_p_AntiLambda_Pion",";p (GeV); Events/0.01GeV",6000,0,6)
h_pz_AntiLambda = TH1F("h_pz_AntiLambda",";p_{z} (GeV); Events/1GeV",30,0,30)
h_pz_AntiLambda_AntiProton = TH1F("h_pz_AntiLambda_AntiProton",";p_{z} (GeV); Events/1GeV",30,0,30)
h_pz_AntiLambda_Pion = TH1F("h_pz_AntiLambda_Pion",";p_{z} (GeV); Events/1GeV",30,0,30)
h2_pt_pz_AntiLambda_AntiProton = TH2F("h2_pt_pz_AntiLambda_AntiProton",";p_{t} (GeV);|p_{z}| (GeV); Events/(1GeV*0.1GeV)",60,0,6,30,0,30)
h2_pt_pz_AntiLambda_Pion = TH2F("h2_pt_pz_AntiLambda_Pion",";p_{t} (GeV);|p_{z}| (GeV); Events/(1GeV*0.1GeV)",60,0,6,30,0,30)
h2_eta_pt_AntiLambda_AntiProton = TH2F("h2_eta_pt_AntiLambda_AntiProton",";#eta #bar{#Lambda}^{0}-#bar{p};p_{t} #bar{#Lambda}^{0}-#bar{p} (GeV); Events/(0.1#eta*0.1GeV)",100,-5,5,60,0,6)
h2_eta_pt_AntiLambda_Pion = TH2F("h2_eta_pt_AntiLambda_Pion",";#eta #bar{#Lambda}^{0}-#pi^{+};p_{t} #bar{#Lambda}^{0}-#pi^{+} (GeV); Events/(0.1#eta*0.1GeV)",80,-4,4,15,0,1.5)
h2_eta_pz_AntiLambda_AntiProton = TH2F("h2_eta_pz_AntiLambda_AntiProton",";#eta #bar{#Lambda}^{0}-#bar{p};p_{z}| #bar{#Lambda}^{0}-#bar{p} (GeV); Events/(0.1#eta*1GeV)",100,-5,5,60,0,60)
h2_eta_pz_AntiLambda_Pion = TH2F("h2_eta_pz_AntiLambda_Pion",";#eta #bar{#Lambda}^{0}-#pi^{+};|p_{z}| #bar{#Lambda}^{0}-#pi^{+} (GeV); Events/(0.1#eta*0.1GeV)",80,-4,4,60,0,6)
h2_eta_p_AntiLambda_AntiProton = TH2F("h2_eta_p_AntiLambda_AntiProton",";#eta #bar{#Lambda}^{0}-#bar{p};p #bar{#Lambda}^{0}-#bar{p} (GeV); Events/(0.1#eta*0.1GeV)",100,-5,5,600,0,60)
h2_eta_p_AntiLambda_Pion = TH2F("h2_eta_p_AntiLambda_Pion",";#eta #bar{#Lambda}^{0}-#pi^{+};p #bar{#Lambda}^{0}-#pi^{+} (GeV); Events/(0.1#eta*0.1GeV)",80,-4,4,60,0,6)

#the 3D angle(displacement, momentum) of the Ks, Lambda and 4 granddaugthers. This learns what is the the fraction that is actually pointing backwards
h_Ks_openings_angle_displacement_momentum = TH1F("h_Ks_openings_angle_displacement_momentum",";openings angle (#vec{p},#vec{l}_{bs,cv}) (rad);Events/0.1rad ",70,-3.5,3.5)
h_Lambda_openings_angle_displacement_momentum = TH1F("h_Lambda_openings_angle_displacement_momentum",";openings angle (#vec{p},#vec{l}_{bs,cv}) (rad);Events/0.1rad ",70,-3.5,3.5)
h_GEN_Ks_daughters_openings_angle_displacement_momentum = TH1F("h_GEN_Ks_daughters_openings_angle_displacement_momentum",";openings angle (#vec{p},#vec{l}_{bs,cv}) (rad);Events/0.1rad ",70,-3.5,3.5)
h_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum = TH1F("h_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum",";openings angle (#vec{p},#vec{l}_{bs,cv}) (rad);Events/0.1rad ",70,-3.5,3.5)
h_GEN_AntiLambda_Pion_openings_angle_displacement_momentum = TH1F("h_GEN_AntiLambda_Pion_openings_angle_displacement_momentum",";openings angle (#vec{p},#vec{l}_{bs,cv}) (rad);Events/0.1rad ",70,-3.5,3.5)

#PCAs of the granddughters of the antiS
PCA_granddaughters_dir = fOut.mkdir("PCA_granddaughters")
PCA_granddaughters_dir.cd()

h_dxy_Ks_daug0 = TH1F("h_dxy_Ks_daug0",";d_{0} (cm); Events/mm",200,-10,10)
h_dxy_Ks_daug1 = TH1F("h_dxy_Ks_daug1",";d_{0} (cm); Events/mm",200,-10,10)
h_dz_Ks_daug0 = TH1F("h_dz_Ks_daug0",";d_{z} (cm); Events/5cm",60,-150,150)
h_dz_Ks_daug1 = TH1F("h_dz_Ks_daug1",";d_{z} (cm); Events/5cm",60,-150,150)

h_dxy_AntiLambda_AntiProton = TH1F("h_dxy_AntiLambda_AntiProton",";d_{0} (cm); Events/mm",200,-10,10)
h_dxy_AntiLambda_Pion = TH1F("h_dxy_AntiLambda_Pion",";d_{0} (cm); Events/mm",200,-10,10)
h_dz_AntiLambda_AntiProton = TH1F("h_dz_AntiLambda_AntiProton",";d_{z} (cm); Events/5cm",60,-150,150)
h_dz_AntiLambda_Pion = TH1F("h_dz_AntiLambda_Pion",";d_{z} (cm); Events/5cm",60,-150,150)


#info on the decay vertices of the V0s
decay_vertex_dir = fOut.mkdir("decay_vertex_V0s")
decay_vertex_dir.cd()
h_lxy_creation_vertex_Ks_daughters = TH1F("h_lxy_creation_vertex_Ks_daughters",";l_{0,bs} cv (cm);Events/cm ",80,0,80)
h_vz_creation_vertex_Ks_daughters = TH1F("h_vz_creation_vertex_Ks_daughters",";absolute v_{z} cv (cm);Events/10cm",50,-250,250)
h_lxy_creation_vertex_AntiLambda_daughters = TH1F("h_lxy_creation_vertex_AntiLambda_daughters",";l_{0,bs} cv (cm);Events/cm",80,0,80)
h_vz_creation_vertex_AntiLambda_daughters = TH1F("h_vz_creation_vertex_AntiLambda_daughters",";absoolute v_{z} cv (cm);Events/10cm",50,-250,250)

h2_vx_vy_creation_vertex_Ks_daughters = TH2F("h2_vx_vy_creation_vertex_Ks_daughters",";absolute v_{x} decay vertex K_{S}^{0} (cm);absolute v_{y} decay vertex K_{S}^{0} (cm);Events/cm^{2}",160,-80,80,160,-80,80)
h2_vz_lxy_creation_vertex_Ks_daughters = TH2F("h2_vz_lxy_creation_vertex_Ks_daughters",";absolute v_{z} decay vertex K_{S}^{0} (cm);l_{0,bs} decay vertex K_{S}^{0} (cm);Events/cm^{2}",500,-250,250,80,0,80)
h2_vx_vy_creation_vertex_AntiLambda_daughters = TH2F("h2_vx_vy_creation_vertex_AntiLambda_daughters",";absolute v_{x} decay vertex #bar{#Lambda}^{0} (cm); absolute v_{y} decay vertex #bar{#Lambda}^{0} (cm);Events/cm^{2}",160,-80,80,160,-80,80)
h2_vz_lxy_creation_vertex_AntiLambda_daughters = TH2F("h2_vz_lxy_creation_vertex_AntiLambda_daughters",";absolute v_{z} decay vertex #bar{#Lambda}^{0} (cm);l_{0,bs} decay vertex #bar{#Lambda}^{0} (cm);Events/cm^{2}",500,-250,250,80,0,80)


for i in range(0,tree.GetEntries()):

	if i > maxNEntries:
		break

	tree.GetEntry(i)

	#only make these plots for reconstructable antiS which means that all their granddaughters should have 7 hits
	if(tree._GEN_Ks_daughter0_numberOfTrackerHits[0] < 7 or tree._GEN_Ks_daughter1_numberOfTrackerHits[0] < 7 or tree._GEN_AntiLambda_AntiProton_numberOfTrackerHits[0] < 7 or tree._GEN_AntiLambda_Pion_numberOfTrackerHits[0] < 7 ):
		continue
	weight_factor = tree._S_event_weighting_factor[0]*tree._S_event_weighting_factor_PU[0]

	h2_interaction_vertex_vx_vy.Fill(tree._S_vx_interaction_vertex[0],tree._S_vy_interaction_vertex[0],weight_factor)
	h2_interaction_vertex_vx_vy_zoom.Fill(tree._S_vx_interaction_vertex[0],tree._S_vy_interaction_vertex[0],weight_factor)
	h2_interaction_vertex_vz_lxy.Fill(tree._S_vz_interaction_vertex[0],tree._S_lxy_interaction_vertex[0],weight_factor)
	h2_interaction_vertex_vz_lxy_zoom.Fill(tree._S_vz_interaction_vertex[0],tree._S_lxy_interaction_vertex[0],weight_factor)

	h_antiS_eta.Fill( tree._S_eta[0] )	
	h_antiS_lxy_creation_vertex.Fill( np.sqrt( np.power(tree._S_vx[0],2) + np.power(tree._S_vy[0],2) ),weight_factor )	
	h_antiS_vz_creation_vertex.Fill(tree._S_vz[0],weight_factor)	
	h_antiS_lxy_interaction_vertex.Fill(tree._S_lxy_interaction_vertex[0],weight_factor)
	h_antiS_vz_interaction_vertex.Fill(tree._S_vz_interaction_vertex[0],weight_factor)
	h_neutron_momentum.Fill(abs(tree._n_p[0]),weight_factor)
	h_antiS_sumDaughters_openingsangle.Fill(tree._S_sumDaughters_openingsangle[0],weight_factor)	
	h_antiS_sumDaughters_deltaPhi.Fill(tree._S_sumDaughters_deltaPhi[0],weight_factor)	
	h_antiS_sumDaughters_deltaEta.Fill(tree._S_sumDaughters_deltaEta[0],weight_factor)	
	h_antiS_sumDaughters_deltaR.Fill(tree._S_sumDaughters_deltaR[0],weight_factor)	
	h2_antiS_inv_mass_p.Fill(tree._S_mass[0], np.sqrt( tree._S_pz[0]*tree._S_pz[0] + tree._S_pt[0]*tree._S_pt[0] ) ,weight_factor)
	h2_antiS_inv_mass_p_Ks_plus_Lambda.Fill(tree._S_mass[0], np.sqrt( pow(tree._Ks_pz[0] + tree._Lambda_pz[0] , 2) + pow(tree._Ks_pt[0] + tree._Lambda_pt[0] , 2) ) ,weight_factor)
	
	h_pt_Ks.Fill(tree._Ks_pt[0],weight_factor)
	h_p_Ks.Fill( np.sqrt( tree._Ks_pt[0]*tree._Ks_pt[0] + tree._Ks_pz[0]*tree._Ks_pz[0] ) ,weight_factor)
	h_pt_Ks_daug0.Fill(tree._GEN_Ks_daughter0_pt[0],weight_factor)
	h_pt_Ks_daug1.Fill(tree._GEN_Ks_daughter1_pt[0],weight_factor)
	h_p_Ks_daug0.Fill(np.sqrt( np.power(tree._GEN_Ks_daughter0_pt[0],2) + np.power(tree._GEN_Ks_daughter0_pz[0],2) ) ,weight_factor)
	h_p_Ks_daug1.Fill(np.sqrt( np.power(tree._GEN_Ks_daughter1_pt[0],2) + np.power(tree._GEN_Ks_daughter1_pz[0],2) ) ,weight_factor)
	h_pz_Ks.Fill(abs(tree._Ks_pz[0]),weight_factor)
	h_pz_Ks_daug0.Fill(abs(tree._GEN_Ks_daughter0_pz[0]),weight_factor)
	h_pz_Ks_daug1.Fill(abs(tree._GEN_Ks_daughter1_pz[0]),weight_factor)
	h2_pt_pz_Ks_daug0.Fill(tree._GEN_Ks_daughter0_pt[0], abs(tree._GEN_Ks_daughter0_pz[0]),weight_factor)
	h2_pt_pz_Ks_daug1.Fill(tree._GEN_Ks_daughter1_pt[0], abs(tree._GEN_Ks_daughter1_pz[0]),weight_factor)
	h2_eta_pt_Ks_daug0.Fill(tree._GEN_Ks_daughter0_eta[0],tree._GEN_Ks_daughter0_pt[0],weight_factor)
	h2_eta_pt_Ks_daug1.Fill(tree._GEN_Ks_daughter1_eta[0],tree._GEN_Ks_daughter1_pt[0],weight_factor)
	h2_eta_pz_Ks_daug0.Fill(tree._GEN_Ks_daughter0_eta[0],abs(tree._GEN_Ks_daughter0_pz[0]),weight_factor)
	h2_eta_pz_Ks_daug1.Fill(tree._GEN_Ks_daughter1_eta[0],abs(tree._GEN_Ks_daughter1_pz[0]),weight_factor)
	h2_eta_p_Ks_daug0.Fill(tree._GEN_Ks_daughter0_eta[0],np.sqrt(np.power(tree._GEN_Ks_daughter0_pt[0],2)+np.power(tree._GEN_Ks_daughter0_pz[0],2)),weight_factor)
	h2_eta_p_Ks_daug1.Fill(tree._GEN_Ks_daughter1_eta[0],np.sqrt(np.power(tree._GEN_Ks_daughter1_pt[0],2)+np.power(tree._GEN_Ks_daughter1_pz[0],2)),weight_factor)

	h_dxy_Ks_daug0.Fill(tree._GEN_Ks_daughter0_dxy[0],weight_factor)
	h_dxy_Ks_daug1.Fill(tree._GEN_Ks_daughter1_dxy[0],weight_factor)
	h_dz_Ks_daug0.Fill(tree._GEN_Ks_daughter0_dz[0],weight_factor)
	h_dz_Ks_daug1.Fill(tree._GEN_Ks_daughter1_dz[0],weight_factor)


	h_pt_AntiLambda.Fill(tree._Lambda_pt[0],weight_factor)
	h_p_AntiLambda.Fill( np.sqrt( tree._Lambda_pt[0]*tree._Lambda_pt[0] + tree._Lambda_pz[0]*tree._Lambda_pz[0] ) ,weight_factor)
	h_pt_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_pt[0],weight_factor)
	h_pt_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_pt[0],weight_factor)
	h_p_AntiLambda_AntiProton.Fill( np.sqrt( np.power(tree._GEN_AntiLambda_AntiProton_pt[0],2) + np.power(tree._GEN_AntiLambda_AntiProton_pz[0],2) ) ,weight_factor)
	h_p_AntiLambda_Pion.Fill( np.sqrt( np.power(tree._GEN_AntiLambda_Pion_pt[0],2) +  np.power(tree._GEN_AntiLambda_Pion_pz[0],2) ) ,weight_factor)
	h_pz_AntiLambda.Fill(abs(tree._Lambda_pz[0]),weight_factor)
	h_pz_AntiLambda_AntiProton.Fill(abs(tree._GEN_AntiLambda_AntiProton_pz[0]),weight_factor)
	h_pz_AntiLambda_Pion.Fill(abs(tree._GEN_AntiLambda_Pion_pz[0]),weight_factor)
	h2_pt_pz_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_pt[0],abs(tree._GEN_AntiLambda_AntiProton_pz[0]),weight_factor)
	h2_pt_pz_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_pt[0],abs(tree._GEN_AntiLambda_Pion_pz[0]),weight_factor)
	h2_eta_pt_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_eta[0],tree._GEN_AntiLambda_AntiProton_pt[0],weight_factor)
	h2_eta_pt_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_eta[0],tree._GEN_AntiLambda_Pion_pt[0],weight_factor)
	h2_eta_pz_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_eta[0],abs(tree._GEN_AntiLambda_AntiProton_pz[0]),weight_factor)
	h2_eta_pz_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_eta[0],abs(tree._GEN_AntiLambda_Pion_pz[0]),weight_factor)
	h2_eta_p_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_eta[0],np.sqrt(np.power(tree._GEN_AntiLambda_AntiProton_pz[0],2)+np.power(tree._GEN_AntiLambda_AntiProton_pt[0],2)),weight_factor)
	h2_eta_p_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_eta[0],np.sqrt(np.power(tree._GEN_AntiLambda_Pion_pz[0],2)+np.power(tree._GEN_AntiLambda_Pion_pt[0],2)),weight_factor)


	h_Ks_openings_angle_displacement_momentum.Fill(tree._Ks_openings_angle_displacement_momentum[0],weight_factor)
	h_Lambda_openings_angle_displacement_momentum.Fill(tree._Lambda_openings_angle_displacement_momentum[0],weight_factor)
	h_GEN_Ks_daughters_openings_angle_displacement_momentum.Fill(tree._GEN_Ks_daughter0_openings_angle_displacement_momentum[0],weight_factor)
	h_GEN_Ks_daughters_openings_angle_displacement_momentum.Fill(tree._GEN_Ks_daughter1_openings_angle_displacement_momentum[0],weight_factor)
	h_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum.Fill(tree._GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum[0],weight_factor)
	h_GEN_AntiLambda_Pion_openings_angle_displacement_momentum.Fill(tree._GEN_AntiLambda_Pion_openings_angle_displacement_momentum[0],weight_factor)

	h_dxy_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_dxy[0],weight_factor)
	h_dxy_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_dxy[0],weight_factor)
	h_dz_AntiLambda_AntiProton.Fill(tree._GEN_AntiLambda_AntiProton_dz[0],weight_factor)
	h_dz_AntiLambda_Pion.Fill(tree._GEN_AntiLambda_Pion_dz[0],weight_factor)
	
	#need only one of the granddaughter of each V0, as their creation vertex is the same
	h_lxy_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_lxy[0],weight_factor)
	h_vz_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_vz[0],weight_factor)
	h_lxy_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_lxy[0],weight_factor)
	h_vz_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],weight_factor)

	h2_vx_vy_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_vx[0],tree._GEN_Ks_daughter0_vy[0],weight_factor)
	h2_vz_lxy_creation_vertex_Ks_daughters.Fill(tree._GEN_Ks_daughter0_vz[0],tree._GEN_Ks_daughter0_lxy[0],weight_factor)
	h2_vx_vy_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vx[0],tree._GEN_AntiLambda_AntiProton_vy[0],weight_factor)
	h2_vz_lxy_creation_vertex_AntiLambda_daughters.Fill(tree._GEN_AntiLambda_AntiProton_vz[0],tree._GEN_AntiLambda_AntiProton_lxy[0],weight_factor)

	
antiS_properties_dir.cd()

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
	#c.SetLogz()
	c.SaveAs(plots_output_dir+h.GetName()+".pdf")
	c.Write()

h_antiS_eta.Write()

h_antiS_lxy_creation_vertex.Write()
h_antiS_vz_creation_vertex.Write()

c_antiS_lxy_interaction_vertex = TCanvas("c_antiS_lxy_interaction_vertex","");
h_antiS_lxy_interaction_vertex.DrawNormalized()
c_antiS_lxy_interaction_vertex.Write()

c_antiS_vz_interaction_vertex = TCanvas("c_antiS_vz_interaction_vertex","");
h_antiS_vz_interaction_vertex.DrawNormalized()
c_antiS_vz_interaction_vertex.Write()

l_TH1F = [h_antiS_vz_creation_vertex,h_antiS_lxy_interaction_vertex,h_antiS_vz_interaction_vertex]
for h in l_TH1F:
        h.SetDirectory(0)
i_l_TH1F = 0
for h in l_TH1F:
	c= TCanvas(h.GetName(),"");
	if(h.GetSumw2N() == 0):
		h.Sumw2(kTRUE)
	h.Scale(1./h.Integral(), "width");
	h.Draw("CL")
	h.SetStats(0)
	CMS_lumi.CMS_lumi(c, 0, 11)
	#if i_l_TH1F == 1:
	#	c.SetLogy()
	c.SaveAs(plots_output_dir+h.GetName()+".pdf")
	c.Write()
	i_l_TH1F+=1

fermiMomentum_dir = fOut.mkdir("fermiMomentum")
fermiMomentum_dir.cd()

h2_antiS_inv_mass_p.Write()
h2_antiS_inv_mass_p_Ks_plus_Lambda.Write()

l_TH2F = [h2_eta_pt_AntiLambda_Pion,h2_eta_pz_AntiLambda_Pion,h2_eta_p_AntiLambda_Pion]
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


l_TH1F = [h_antiS_sumDaughters_openingsangle, h_antiS_sumDaughters_deltaPhi, h_antiS_sumDaughters_deltaEta, h_antiS_sumDaughters_deltaR]
for h in l_TH1F:
	c = TCanvas("c_"+h.GetName(),"");
	h.Draw("PCE1")
	if(h.GetSumw2N() == 0):
		h.Sumw2(kTRUE)
	h.Scale(1./h.Integral(), "width");
	h.SetLineColor(colours[0])
	h.SetMarkerStyle(22)
	h.SetMarkerColor(colours[0])
	h.SetStats(0)
	CMS_lumi.CMS_lumi(c, 0, 11)
	#c_h_antiS_sumDaughters_openingsangle.SetLogy()
	c.SaveAs(plots_output_dir+c.GetName()+".pdf")
	c.Write()


c_h2_antiS_inv_mass_p_Ks_plus_Lambda= TCanvas(h2_antiS_inv_mass_p_Ks_plus_Lambda.GetName(),"");
c_h2_antiS_inv_mass_p_Ks_plus_Lambda.SetRightMargin(0.2) #make room for the tile of the z scale
if(h2_antiS_inv_mass_p_Ks_plus_Lambda.GetSumw2N() == 0):
	h2_antiS_inv_mass_p_Ks_plus_Lambda.Sumw2(kTRUE)
h2_antiS_inv_mass_p_Ks_plus_Lambda.Scale(1./h2_antiS_inv_mass_p_Ks_plus_Lambda.Integral(), "width");
#h2_antiS_inv_mass_p_Ks_plus_Lambda.Draw("SURF1") #for 3D plot
h2_antiS_inv_mass_p_Ks_plus_Lambda.Draw("colz")
h2_antiS_inv_mass_p_Ks_plus_Lambda.SetStats(0)
h2_antiS_inv_mass_p_Ks_plus_Lambda.GetXaxis().SetTitleOffset(1.5)
h2_antiS_inv_mass_p_Ks_plus_Lambda.GetYaxis().SetTitleOffset(2)
c_h2_antiS_inv_mass_p_Ks_plus_Lambda.SetLogz()
CMS_lumi.CMS_lumi(c_h2_antiS_inv_mass_p_Ks_plus_Lambda, 0, 11)
c_h2_antiS_inv_mass_p_Ks_plus_Lambda.SaveAs(plots_output_dir+h2_antiS_inv_mass_p_Ks_plus_Lambda.GetName()+".pdf")
c_h2_antiS_inv_mass_p_Ks_plus_Lambda.Write()

momenta_daughters_and_grandaughters_dir.cd()

ll_TH1F = [
[h_pt_Ks_daug0,h_pt_Ks],
[h_pz_Ks_daug0,h_pz_Ks],
[h_pz_AntiLambda_Pion,h_pz_AntiLambda_AntiProton,h_pz_AntiLambda],
[h_pt_AntiLambda_Pion,h_pt_AntiLambda_AntiProton,h_pt_AntiLambda],
[h_GEN_Ks_daughters_openings_angle_displacement_momentum,h_Ks_openings_angle_displacement_momentum],
[h_GEN_AntiLambda_Pion_openings_angle_displacement_momentum,h_GEN_AntiLambda_AntiProton_openings_angle_displacement_momentum,h_Lambda_openings_angle_displacement_momentum]
]


ll_legend_text =  [
["K_{S}^{0} daughters","K_{S}^{0}"],
["K_{S}^{0} daughters","K_{S}^{0}"],
["#bar{#Lambda}^{0}-#pi^{+}","#bar{#Lambda}^{0}-#bar{p}","#bar{#Lambda}^{0}"],
["#bar{#Lambda}^{0}-#pi^{+}","#bar{#Lambda}^{0}-#bar{p}","#bar{#Lambda}^{0}"],
["K_{S}^{0} daughters","K_{S}^{0}"],
["#bar{#Lambda}^{0}-#pi^{+}","#bar{#Lambda}^{0}-#bar{p}","#bar{#Lambda}^{0}"]
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
legend_dxyGranddaughters = TLegend(0.1,0.7,0.55,0.9);
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
["K_{S}^{0} daughters","#bar{#Lambda}^{0} daughters"],
["K_{S}^{0} daughters","#bar{#Lambda}^{0} daughters"],
["#bar{#Lambda}^{0}-#bar{p}","#bar{#Lambda}^{0}-#pi^{+}","K_{S} daughters"],
["#bar{#Lambda}^{0}-#bar{p}","#bar{#Lambda}^{0}-#pi^{+}","K_{S} daughters"],
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
