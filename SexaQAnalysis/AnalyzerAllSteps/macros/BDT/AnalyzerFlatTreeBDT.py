#from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, gROOT
import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

import sys
sys.path.insert(1, '/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA')
import configBDT as config

config_dict = config.config_dict

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

colours = [1,2,4,35,38,41]

maxEvents = 5e4



#you can use this script to compared different collections of S and antiS to eachother. The most common one is comparing Data-S-BKG to MC-AntiS-Signal. The other option is comparing MC-S-BKG to MC-AntiS-BKG to Data-S-BKG
configuration = "all" # "Data-S-BKG_to_MC-AntiS-Signal" or "MC-S-BKG_to_MC-AntiS-BKG_to_Data-S-BKG" or "all" 

# Open file
SignFile1 = ROOT.TFile.Open(config_dict["config_SignalFile"]) 

BkgFile  = ROOT.TFile.Open(config_dict["config_BkgFile"]) 

# Get signal and background trees from file
SignalTree1     = SignFile1.Get("FlatTreeProducerBDT/FlatTree")
#I only want the antiS which match a GEN antiS in lxyz of the interaction vertex and the charge should also be negative
gROOT.cd()
selectedSignalTree1 = SignalTree1.CopyTree(config_dict["config_SelectionSignalAntiS"]+ ' && ' + config_dict["config_pre_BDT_cuts"])



BkgTree        = BkgFile.Get("FlatTreeProducerBDT/FlatTree")
gROOT.cd()
#for selecting MC S BKG: 
selectedBkgTree_MC_S_BKG = SignalTree1.CopyTree(config_dict["config_SelectionBkgS"] + ' && ' + config_dict["config_pre_BDT_cuts"])
#for selecting MC antiS BKG: 
selectedBkgTree_MC_AntiS_BKG = SignalTree1.CopyTree(config_dict["config_SelectionBkgAntiS"] +' &&' + config_dict["config_pre_BDT_cuts"])
#for selecting Data S BKG: ---> standard
selectedBkgTree_Data_S_BKG = BkgTree.CopyTree(config_dict["config_SelectionBkgS"] + ' && ' + config_dict["config_pre_BDT_cuts"])

l_y_axis_ranges = [
0.08,
5.,
1.1,
1,
2.,
1.4,
3,
3,
1.,
0.9,
0.8,
15,
4,
12,
3,
10,
1,
1.2,
0.15,
1.4,
0.05,
20.,
0.7
]

Legend_Data_S_BKG_to_MC_AntiS_Signal = ["MC-#bar{S}-Signal", "Data-S-BKG"]
Legend_MC_S_BKG_to_MC_AntiS_BKG_to_Data_S_BKG = ["MC-S-BKG","MC-#bar{S}-BKG","Data-S-BKG"]
Legend_all = ["MC-#bar{S}-Signal","MC-S-BKG","MC-#bar{S}-BKG","Data-S-BKG"]

plots_output_dir = "plots_BackgroundVsSignal/"+configuration+"/"

Legend = Legend_Data_S_BKG_to_MC_AntiS_Signal
l_tree = [selectedSignalTree1,selectedBkgTree_Data_S_BKG]

if(configuration == "MC-S-BKG_to_MC-AntiS-BKG_to_Data-S-BKG"):
	Legend = Legend_MC_S_BKG_to_MC_AntiS_BKG_to_Data_S_BKG
	l_tree = [selectedBkgTree_MC_S_BKG, selectedBkgTree_MC_AntiS_BKG, selectedBkgTree_Data_S_BKG]

if(configuration == "all"):
	Legend = Legend_all
	l_tree = [selectedSignalTree1,selectedBkgTree_MC_S_BKG, selectedBkgTree_MC_AntiS_BKG, selectedBkgTree_Data_S_BKG]


TH1_ll = [] #list of list of 1D histos 
TH2_ll = [] #list of list of 2D histos

iFile = 0
for tree in l_tree:

	h_S_vz_interaction_vertex= TH1F('h_S_vz_interaction_vertex','; v_{z} ^{(}#bar{S} ^{)}) interaction vertex (cm); 1/N_{ev} Events/5cm',20,-50,50)
	h_S_lxy_interaction_vertex = TH1F('h_S_lxy_interaction_vertex','; l_{0} interaction vertex ^{(}#bar{S} ^{)} (cm); 1/N_{ev} Events/0.5mm',16,1.8,2.6)

	h_S_daughters_deltaphi = TH1F('h_S_daughters_deltaphi','; #Delta#phi(K_{S}, ^{(} #bar{#Lambda} ^{)} ) (rad); 1/N_{ev} Events/0.1rad',70,-3.5,3.5)
	h_S_daughters_deltaeta = TH1F('h_S_daughters_deltaeta','; #Delta#eta(K_{S}, ^{(} #bar{#Lambda} ^{)} ) ; 1/N_{ev} Events/0.1rad',60,-3,3)
	h_S_daughters_openingsangle = TH1F('h_S_daughters_openingsangle','; openings angle(K_{S},^{(} #bar{#Lambda} ^{)} ) (rad); 1/N_{ev} Events/0.1rad',35,0,3.5)
	h_S_daughters_DeltaR = TH1F('h_S_daughters_DeltaR','; #DeltaR(K_{S}, ^{(} #bar{#Lambda} ^{)} ); 1/N_{ev} Events/0.1#DeltaR',60,0.5,6.5)
	h_S_Ks_openingsangle = TH1F('h_S_Ks_openingsangle','; openings angle( ^{(} #bar{S} ^{)} ,K_{S}) (rad); 1/N_{ev} Events/0.1rad',20,0,2)
	h_S_Lambda_openingsangle = TH1F('h_S_Lambda_openingsangle','; openings angle( ^{(} #bar{S} ^{)},^{(} #bar{#Lambda} ^{)} ) (rad); 1/N_{ev} Events/',20,0,2)

	h_S_eta = TH1F('h_S_eta','; #eta( ^{(} #bar{S} ^{)} ); 1/N_{ev} Events/0.1#eta',100,-5,5)
	h_Ks_eta = TH1F('h_Ks_eta','; #eta(K_{S}) ; 1/N_{ev} Events/0.1#eta',100,-5,5)
	h_Lambda_eta = TH1F('h_Lambda_eta','; #eta( ^{(} #bar{#Lambda} ^{)} ); 1/N_{ev} Events/0.1#eta',100,-5,5)

	h_S_dxy_over_lxy = TH1F('h_S_dxy_over_lxy','; d_{0}/l_{0} ( ^{(} #bar{S} ^{)} ); 1/N_{ev} Events/0.1',20,-1,1)
	h_Ks_dxy_over_lxy = TH1F('h_Ks_dxy_over_lxy','; d_{0}/l_{0} (K_{S}); 1/N_{ev} Events/0.1',20,-1,1)
	h_Lambda_dxy_over_lxy = TH1F('h_Lambda_dxy_over_lxy','; d_{0}/l_{0} ( ^{(} #bar{#Lambda} ^{)}) ; 1/N_{ev} Events/0.1',20,-1,1)

	h_S_dz_min = TH1F('h_S_dz_min','; min d_{z}  ^{(} #bar{S} ^{)}  (cm); 1/N_{ev} Events/cm',20,-10,10)
	h_Ks_dz_min = TH1F('h_Ks_dz_min','; min d_{z} K_{S} (cm); 1/N_{ev} Events/cm',60,-30,30)
	h_Lambda_dz_min = TH1F('h_Lambda_dz_min','; min d_{z}  ^{(} #bar{#Lambda} ^{)}  (cm); 1/N_{ev} Events/cm',60,-30,30)

	h_Ks_pt = TH1F('h_Ks_pt','; p_{t} K_{S} (GeV); 1/N_{ev} Events/0.1GeV',80,0,8)

	h_Lambda_lxy_decay_vertex = TH1F('h_Lambda_lxy_decay_vertex','; l_{0} ^{(} #bar{#Lambda} ^{)} decay vertex (cm); 1/N_{ev} Events/cm',20,1.9,21.9)
	h_S_chi2_ndof = TH1F('h_S_chi2_ndof','; #chi^{2}/ndof ^{(}#bar{S} ^{)} annihilation vertex; 1/N_{ev} Events',110,0,11)

	h_S_pz = TH1F('h_S_pz','; p_{z}  ( ^{(} #bar{S} ^{)} ) (GeV); 1/N_{ev} Events/5GeV',16,-40,40)

	h_S_error_lxy_interaction_vertex = TH1F('h_S_error_lxy_interaction_vertex','; #sigma(l_{0} interaction vertex ^{(}#bar{S} ^{)}) (cm); 1/N_{ev} Events/0.1mm',15,0,0.15)
	h_S_mass = TH1F('h_S_mass','; m_{ ^{(} #bar{S} ^{)} ,obs} (GeV); 1/N_{ev} Events/0.1GeV',100,-5,5)

	nEntries = tree.GetEntries()
	print 'Number of entries in the tree: ', nEntries
	for i in range(0,nEntries):
		if(i==maxEvents):
			break
		if(i%1e4 == 0):
			print "reached entry: ", i
		tree.GetEntry(i)
#		if(tree._S_mass[0] < 0):
#			continue
		h_S_vz_interaction_vertex.Fill(tree._S_vz_interaction_vertex[0])
		h_S_lxy_interaction_vertex.Fill(tree._S_lxy_interaction_vertex[0])

		h_S_daughters_deltaphi.Fill(tree._S_daughters_deltaphi[0])
		h_S_daughters_deltaeta.Fill(tree._S_daughters_deltaeta[0])
		h_S_daughters_openingsangle.Fill(tree._S_daughters_openingsangle[0])
		h_S_daughters_DeltaR.Fill(tree._S_daughters_DeltaR[0])
		h_S_Ks_openingsangle.Fill(tree._S_Ks_openingsangle[0])
		h_S_Lambda_openingsangle.Fill(tree._S_Lambda_openingsangle[0])

		h_S_eta.Fill(tree._S_eta[0])
		h_Ks_eta.Fill(tree._Ks_eta[0])
		h_Lambda_eta.Fill(tree._Lambda_eta[0])

		h_S_dxy_over_lxy.Fill(tree._S_dxy_over_lxy[0])
		h_Ks_dxy_over_lxy.Fill(tree._Ks_dxy_over_lxy[0])
		h_Lambda_dxy_over_lxy.Fill(tree._Lambda_dxy_over_lxy[0])

		h_S_dz_min.Fill(tree._S_dz_min[0])
		h_Ks_dz_min.Fill(tree._Ks_dz_min[0])
		h_Lambda_dz_min.Fill(tree._Lambda_dz_min[0])

		h_Ks_pt.Fill(tree._Ks_pt[0])
		
		h_Lambda_lxy_decay_vertex.Fill(tree._Lambda_lxy_decay_vertex[0])
		h_S_chi2_ndof.Fill(tree._S_chi2_ndof[0])

		h_S_pz.Fill(tree._S_pz[0])

		h_S_error_lxy_interaction_vertex.Fill(tree._S_error_lxy_interaction_vertex[0])  
		h_S_mass.Fill(tree._S_mass[0])	

	TH1_l = [h_S_vz_interaction_vertex,h_S_lxy_interaction_vertex,h_S_daughters_deltaphi,h_S_daughters_deltaeta,h_S_daughters_openingsangle,h_S_daughters_DeltaR,h_S_Ks_openingsangle,h_S_Lambda_openingsangle,h_S_eta,h_Ks_eta,h_Lambda_eta,h_S_dxy_over_lxy,h_Ks_dxy_over_lxy,h_Lambda_dxy_over_lxy,h_S_dz_min,h_Ks_dz_min,h_Lambda_dz_min,h_Ks_pt,h_Lambda_lxy_decay_vertex,h_S_chi2_ndof,h_S_pz,h_S_error_lxy_interaction_vertex,h_S_mass]
	for h in TH1_l:
		h.SetDirectory(0) 
	TH1_ll.append(TH1_l)

	TH2_l = []
	for h in TH2_l:
		h.SetDirectory(0) 
	TH2_ll.append(TH2_l)

	iFile+=1

fOut = TFile(plots_output_dir+'macro_FlatTree_BDT_trial17.root','RECREATE')

nHistos = len(TH1_ll[0])
nSamples = len(TH1_ll)

c_all  = TCanvas("c_all","c_all",1800,1400)
c_all_log  = TCanvas("c_all_log","c_all_log",1600,1200)
c_all.Divide(6,5)
for i in range(0,nHistos):#each list contains a list of histograms. Each list represents background or signal. The histos need to be overlaid one list to the other
	h = TH1_ll[0][i]
	c_name = "c_"+h.GetName()
	c = TCanvas(c_name,"")
	legend = TLegend(0.7,0.85,0.99,0.99)
	for j in range(0,nSamples):
		h = TH1_ll[j][i]
		if(h.GetSumw2N() == 0):
			h.Sumw2(kTRUE)
		if(h.Integral() != 0):
			h.Scale(1./h.Integral(), "width");
		h.GetYaxis().SetRangeUser(0.,l_y_axis_ranges[i])
		if("dxy_over_lxy" in h.GetName()):
			h.GetYaxis().SetRangeUser(1e-2,l_y_axis_ranges[i])
		if("dz_min" in h.GetName()):
			h.GetYaxis().SetRangeUser(1e-5,l_y_axis_ranges[i])
		if("chi2_ndof" in h.GetName()):
			h.GetYaxis().SetRangeUser(1e-3,l_y_axis_ranges[i])
		if j == 0:
			h.Draw("L")
			c.Update()
		else:
			h.Draw("PCE1same")
			c.Update()
		if("dxy_over_lxy" in h.GetName() or "dz_min" in h.GetName() or "chi2_ndof" in h.GetName()):
			c.SetLogy()
		h.SetLineColor(colours[j])
		h.SetLineWidth(1)
		h.SetMarkerStyle(22+j)
		h.SetMarkerColor(colours[j])
		h.SetMarkerSize(0.5)
		h.SetStats(0)
		legend.AddEntry(h, Legend[j],"lep")
	legend.Draw()
	CMS_lumi.CMS_lumi(c, 0, 11)
	c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
	c.Write()

	c_all.cd(i+1)
	c.SetGridx(0)
	c.SetGridy(0)
	c.DrawClonePad()
	

c_all.Write()
c_all.SaveAs(plots_output_dir+"c_all.pdf")
c_all_log.Write()
c_all_log.SaveAs(plots_output_dir+"c_all_log.pdf")


i = 0
for l in TH2_ll:
	for h in l:
		c_name = "c_"+h.GetName()+"_"+Legend[i]
		c = TCanvas(c_name,"");
		c.SetRightMargin(0.2) #make room for the tile of the z scale
		if(h.GetSumw2N() == 0):
			h.Sumw2(kTRUE)
		h.Scale(1./h.Integral(), "width");
		h.Draw("colz")
		h.SetStats(0)
		CMS_lumi.CMS_lumi(c, 0, 11)
		c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
		c.Write()
	i += 1

#c_antiS_eta_pt = TCanvas("c_antiS_eta_pt","")
#h_antiS_eta_pt.DrawNormalized()
#c_antiS_eta_pt.Write()
#
#c_antiS_eta_pz = TCanvas("c_antiS_eta_pz","")
#h_antiS_eta_pz.DrawNormalized()
#c_antiS_eta_pz.Write()

fOut.Write()
fOut.Close()
