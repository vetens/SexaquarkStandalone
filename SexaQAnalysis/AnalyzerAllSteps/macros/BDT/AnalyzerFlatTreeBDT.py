#from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, gROOT
from ROOT import *
import numpy as np
import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

colours = [1,2,4,35,38,41]

#you can use this script to compared different collections of S and antiS to eachother. The most common one is comparing Data-S-BKG to MC-AntiS-Signal. The other option is comparing MC-S-BKG to MC-AntiS-BKG to Data-S-BKG
configuration = "compare_NEG_to_POS_InvMass" # "Data-S-BKG_to_MC-AntiS-Signal" or "MC-S-BKG_to_MC-AntiS-BKG_to_Data-S-BKG" or "all" or "compare_NEG_to_POS_InvMass" 


# Open file
#SignFile1 = ROOT.TFile.Open("/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/test_FlatTreeBDT_trial15.root")
SignFile1 = TFile.Open("/pnfs/iihe/cms/store/user/lowette/crmc_Sexaq/Skimmed/CRAB_SimSexaq_trial17/crab_Step1_Step2_Skimming_FlatTree_trial17_18092019_v1/190918_051631/combined_FlatTreeBDT_Skimmed_trial17_21.root")

BkgFile  = TFile.Open("/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerBDT/test_FlatTreeBDT_SingleMuon_Run2016H.root")

# Get signal and background trees from file
SignalTree1     = SignFile1.Get("FlatTreeProducerBDT/FlatTree")
#I only want the antiS which match a GEN antiS in lxyz of the interaction vertex and the charge should also be negative
gROOT.cd()
selectedSignalTree1 = SignalTree1.CopyTree('Alt$(_S_charge,0) == -1 && Alt$(_S_deltaLInteractionVertexAntiSmin,0) < 0.5 && Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) && Alt$(_S_dxy_over_lxy,0) >= 0')
selectedSignalTree1_posMass = SignalTree1.CopyTree('Alt$(_S_charge,0) == -1 && Alt$(_S_deltaLInteractionVertexAntiSmin,0) < 0.5 && Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) && Alt$(_S_dxy_over_lxy,0) >= 0 &&   Alt$(_S_mass,0) > 0' )
selectedSignalTree1_negMass = SignalTree1.CopyTree('Alt$(_S_charge,0) == -1 && Alt$(_S_deltaLInteractionVertexAntiSmin,0) < 0.5 && Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) && Alt$(_S_dxy_over_lxy,0) >= 0 && Alt$(_S_mass,0) < 0')


BkgTree        = BkgFile.Get("FlatTreeProducerBDT/FlatTree")
gROOT.cd()
#for selecting MC S BKG: 
selectedBkgTree_MC_S_BKG = SignalTree1.CopyTree('Alt$(_S_charge,0) == 1 && Alt$(_S_deltaLInteractionVertexAntiSmin,0) > 0 && Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) && Alt$(_S_dxy_over_lxy,0) >= 0')
#for selecting MC antiS BKG: 
selectedBkgTree_MC_AntiS_BKG = SignalTree1.CopyTree('Alt$(_S_charge,0) == -1 && Alt$(_S_deltaLInteractionVertexAntiSmin,0) > 1 && Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) && Alt$(_S_dxy_over_lxy,0) >= 0')
#for selecting Data S BKG: ---> standard
selectedBkgTree_Data_S_BKG = BkgTree.CopyTree('Alt$(_S_charge,0) == 1 && Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) && Alt$(_S_dxy_over_lxy,0) >= 0')
selectedBkgTree_Data_S_BKG_posMass = BkgTree.CopyTree('Alt$(_S_charge,0) == 1 && Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) && Alt$(_S_dxy_over_lxy,0) >= 0 && Alt$(_S_mass,0) > 0')
selectedBkgTree_Data_S_BKG_negMass = BkgTree.CopyTree('Alt$(_S_charge,0) == 1 && Alt$(_S_error_lxy_interaction_vertex,0) < 0.13 && (Alt$(_S_daughters_deltaphi,0) < -1 || Alt$(_S_daughters_deltaphi,0) > 1) && Alt$(_S_dxy_over_lxy,0) >= 0 && Alt$(_S_mass,0) < 0')

l_y_axis_ranges = [0.03,1,20,1.2,1.2,3,1.5,3.5,4.5,0.7,0.7,0.7,12,3,10,1,0.7,1,1.4,0.6]

Legend_Data_S_BKG_to_MC_AntiS_Signal = ["MC-#bar{S}-Signal", "Data-S-BKG"]
Legend_MC_S_BKG_to_MC_AntiS_BKG_to_Data_S_BKG = ["MC-S-BKG","MC-#bar{S}-BKG","Data-S-BKG"]
Legend_all = ["MC-#bar{S}-Signal","MC-S-BKG","MC-#bar{S}-BKG","Data-S-BKG"]
Legend_compare_NEG_to_POS_InvMass = ["MC-#bar{S}-Signal M > 0","MC-#bar{S}-Signal M < 0","Data-S-BKG M > 0","Data-S-BKG M < 0"]

plots_output_dir = "plots_BackgroundVsSignal/"+configuration+"/"

Legend = Legend_Data_S_BKG_to_MC_AntiS_Signal
l_tree = [selectedSignalTree1,selectedBkgTree_Data_S_BKG]

if(configuration == "MC-S-BKG_to_MC-AntiS-BKG_to_Data-S-BKG"):
	Legend = Legend_MC_S_BKG_to_MC_AntiS_BKG_to_Data_S_BKG
	l_tree = [selectedBkgTree_MC_S_BKG, selectedBkgTree_MC_AntiS_BKG, selectedBkgTree_Data_S_BKG]

if(configuration == "all"):
	Legend = Legend_all
	l_tree = [selectedSignalTree1,selectedBkgTree_MC_S_BKG, selectedBkgTree_MC_AntiS_BKG, selectedBkgTree_Data_S_BKG]

if(configuration == "compare_NEG_to_POS_InvMass"):
	Legend = Legend_compare_NEG_to_POS_InvMass
	l_tree = [selectedSignalTree1_posMass,selectedSignalTree1_negMass,selectedBkgTree_Data_S_BKG_posMass,selectedBkgTree_Data_S_BKG_negMass]

TH1_ll = [] #list of list of 1D histos 
TH2_ll = [] #list of list of 2D histos

iFile = 0
for tree in l_tree:

	h_Ks_vz_decay_vertex = TH1F('h_Ks_vz_decay_vertex','; v_{z} K_{s} decay vertex (cm); 1/N_{ev} Events/5GeV',84,-120,120)
	h_S_lxy_interaction_vertex = TH1F('h_S_lxy_interaction_vertex','; l_{0} interaction vertex ^{(}#bar{S} ^{)} (cm); 1/N_{ev} Events/cm',81,1.9,10)
	h_S_error_lxy_interaction_vertex = TH1F('h_S_error_lxy_interaction_vertex','; #sigma(l_{0} interaction vertex ^{(}#bar{S} ^{)}) (cm); 1/N_{ev} Events/0.1mm',300,0,3)

	h_S_daughters_deltaphi = TH1F('h_S_daughters_deltaphi','; #Delta#phi(K_{S}, ^{(} #bar{#Lambda} ^{)} ) (rad); 1/N_{ev} Events/0.1rad',70,-3.5,3.5)
	h_S_daughters_deltaeta = TH1F('h_S_daughters_deltaeta','; #Delta#eta(K_{S}, ^{(} #bar{#Lambda} ^{)} ) ; 1/N_{ev} Events/0.1rad',60,-3,3)
	h_S_daughters_openingsangle = TH1F('h_S_daughters_openingsangle','; openings angle(K_{S},^{(} #bar{#Lambda} ^{)} ) (rad); 1/N_{ev} Events/0.1rad',35,0,3.5)
	h_S_daughters_DeltaR = TH1F('h_S_daughters_DeltaR','; #DeltaR(K_{S}, ^{(} #bar{#Lambda} ^{)} ); 1/N_{ev} Events/0.1#DeltaR',40,0,4)
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
	h_S_mass = TH1F('h_S_mass','; m_{ ^{(} #bar{S} ^{)} ,obs} (GeV); 1/N_{ev} Events/0.1GeV',100,-5,5)
 


	nEntries = tree.GetEntries()
	print 'Number of entries in the tree: ', nEntries
	for i in range(0,nEntries):
		if(i==5e4):
			break
		if(i%1e4 == 0):
			print "reached entry: ", i
		tree.GetEntry(i)
#		if(tree._S_mass[0] < 0):
#			continue
		h_Ks_vz_decay_vertex.Fill(tree._Ks_vz_decay_vertex[0])
		h_S_lxy_interaction_vertex.Fill(tree._S_lxy_interaction_vertex[0])
		h_S_error_lxy_interaction_vertex.Fill(tree._S_error_lxy_interaction_vertex[0])  

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
		h_S_mass.Fill(tree._S_mass[0])	

	TH1_l = [h_Ks_vz_decay_vertex,h_S_lxy_interaction_vertex,h_S_error_lxy_interaction_vertex,h_S_daughters_deltaphi,h_S_daughters_deltaeta,h_S_daughters_openingsangle,h_S_daughters_DeltaR,h_S_Ks_openingsangle,h_S_Lambda_openingsangle,h_S_eta,h_Ks_eta,h_Lambda_eta,h_S_dxy_over_lxy,h_Ks_dxy_over_lxy,h_Lambda_dxy_over_lxy,h_S_dz_min,h_Ks_dz_min,h_Lambda_dz_min,h_Ks_pt,h_S_mass]
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

c_all  = TCanvas("c_all","c_all",1600,1200)
c_all_log  = TCanvas("c_all_log","c_all_log",1600,1200)
c_all.Divide(5,4)
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
		if j == 0:
			h.Draw("L")
			c.Update()
		else:
			h.Draw("PCE1same")
			c.Update()
		if("dxy_over_lxy" in h.GetName() or "dz_min" in h.GetName()):
			c.SetLogy()
		h.SetLineColor(colours[j])
		h.SetMarkerStyle(22+j)
		h.SetMarkerColor(colours[j])
		h.SetMarkerSize(0.6)
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
