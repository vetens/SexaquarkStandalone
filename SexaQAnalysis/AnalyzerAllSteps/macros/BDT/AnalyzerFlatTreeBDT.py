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

colours = [2,4,35,1,38,41]

maxEvents = 5e4


#you can use this script to compared different collections of S and antiS to eachother. The most common one is comparing Data-S-BKG to MC-AntiS-Signal. The other option is comparing MC-S-BKG to MC-AntiS-BKG to Data-S-BKG
configuration = "all" # "MC-AntiS-Signal" "Data-S-BKG_to_MC-AntiS-Signal" or "MC-S-BKG_to_MC-AntiS-BKG_to_Data-S-BKG" or "all" 

# Open file
SignFile1 = ROOT.TFile.Open(config_dict["config_SignalFile"]) 

BkgFile  = ROOT.TFile.Open(config_dict["config_BkgFile"]) 

# Get signal and background trees from file
SignalTree1     = SignFile1.Get("FlatTreeProducerBDT/FlatTree")
#I only want the antiS which match a GEN antiS in lxyz of the interaction vertex and the charge should also be negative
gROOT.cd()
selectedSignalTree1_preBDT_cuts = SignalTree1.CopyTree(config_dict["config_SelectionSignalAntiS"])
selectedSignalTree1 = SignalTree1.CopyTree(config_dict["config_SelectionSignalAntiS"]+ ' && ' + config_dict["config_pre_BDT_cuts"])
#for the MC signal the efficiency calculation is not so easy as for the MC samples below as I have to reweigh the events
den = 0.
nom = 0.
for entry in selectedSignalTree1_preBDT_cuts:
	den += config.calc_reweighing_factor(selectedSignalTree1_preBDT_cuts._S_eta[0],True)	
for entry in  selectedSignalTree1:
	nom += config.calc_reweighing_factor(selectedSignalTree1._S_eta[0],True)
print "The efficiency for signal for the pre-BDT cuts: ", nom/den


BkgTree        = BkgFile.Get("FlatTreeProducerBDT/FlatTree")
gROOT.cd()
#for selecting MC S BKG: 
nEntries_MC_S_BKG_preBDTCuts = SignalTree1.GetEntries(config_dict["config_SelectionBkgS"])
selectedBkgTree_MC_S_BKG = SignalTree1.CopyTree(config_dict["config_SelectionBkgS"] + ' && ' + config_dict["config_pre_BDT_cuts"])
print "The efficiency for MC_S_BKG for the pre-BDT cuts: ", float(selectedBkgTree_MC_S_BKG.GetEntries())/float(nEntries_MC_S_BKG_preBDTCuts)
#for selecting MC antiS BKG: 
nEntries_MC_AntiS_BKG_preBDTCuts = SignalTree1.GetEntries(config_dict["config_SelectionBkgAntiS"])
selectedBkgTree_MC_AntiS_BKG = SignalTree1.CopyTree(config_dict["config_SelectionBkgAntiS"] +' &&' + config_dict["config_pre_BDT_cuts"])
print "The efficiency for MC_AntiS_BKG for the pre-BDT cuts: ", float(selectedBkgTree_MC_AntiS_BKG.GetEntries())/float(nEntries_MC_AntiS_BKG_preBDTCuts)
#for selecting Data S BKG: ---> standard
nEntries_Data_S_BKG_preBDTCuts = BkgTree.GetEntries(config_dict["config_SelectionBkgS"])
selectedBkgTree_Data_S_BKG = BkgTree.CopyTree(config_dict["config_SelectionBkgS"] + ' && ' + config_dict["config_pre_BDT_cuts"])
print "The efficiency for Data_S_BKG for the pre-BDT cuts: ", float(selectedBkgTree_Data_S_BKG.GetEntries())/float(nEntries_Data_S_BKG_preBDTCuts)

l_y_axis_ranges = [
0.08,
14.,
1.3,
1,
3.,
2,
3,
4.5,
1,
1,
1,
15,
40,
12,
5,
10,
1,
1.2,
0.20,
2.,
0.05,
250.,
0.7,
50
]

Legend_MC_AntiS_Signal = ["MC-S-BKG","MC-#bar{S}-Signal"]
Legend_Data_S_BKG_to_MC_AntiS_Signal = ["MC-#bar{S}-Signal", "Data-S-BKG"]
Legend_MC_S_BKG_to_MC_AntiS_BKG_to_Data_S_BKG = ["MC-S-BKG","MC-#bar{S}-BKG","Data-S-BKG"]
Legend_all = ["MC-S-BKG","MC-#bar{S}-BKG","Data-S-BKG","MC-#bar{S}-Signal"]

plots_output_dir = "plots_BackgroundVsSignal/"+configuration+"/"

Legend = Legend_Data_S_BKG_to_MC_AntiS_Signal
l_tree = [selectedSignalTree1,selectedBkgTree_Data_S_BKG]

if(configuration == "MC-AntiS-Signal"):
	Legend = Legend_MC_AntiS_Signal
	l_tree = [selectedSignalTree1]

if(configuration == "MC-S-BKG_to_MC-AntiS-BKG_to_Data-S-BKG"):
	Legend = Legend_MC_S_BKG_to_MC_AntiS_BKG_to_Data_S_BKG
	l_tree = [selectedBkgTree_MC_S_BKG, selectedBkgTree_MC_AntiS_BKG, selectedBkgTree_Data_S_BKG]

if(configuration == "all"):
	Legend = Legend_all
	l_tree = [selectedBkgTree_MC_S_BKG, selectedBkgTree_MC_AntiS_BKG, selectedBkgTree_Data_S_BKG,selectedSignalTree1]

#first loop over all the entries to see how much each of the pre-BDT cuts is excluding
numerators = [[0,0,0,0]*len(l_tree)] #for each tree there are 4 cuts to evaluate
denominators = [[0,0,0,0]*len(l_tree)] #for each tree there are 4 cuts to evaluate

TH1_ll = [] #list of list of 1D histos 
TH2_ll = [] #list of list of 2D histos

iTree = 0
for tree in l_tree:

	h_S_vz_interaction_vertex= TH1F('h_S_vz_interaction_vertex','; v_{z} ^{(}#bar{S} ^{)}) interaction vertex w.r.t beampipe center (cm); 1/N_{ev} Events/3cm',20,-40,40)
	h_S_lxy_interaction_vertex = TH1F('h_S_lxy_interaction_vertex','; l_{0} interaction vertex ^{(}#bar{S} ^{)} (cm); 1/N_{ev} Events/0.2mm',16,2.05,2.37)

	h_S_daughters_deltaphi = TH1F('h_S_daughters_deltaphi','; #Delta#phi(K_{S}, ^{(} #bar{#Lambda} ^{)} ) (rad); 1/N_{ev} Events/0.1rad',70,-3.5,3.5)
	h_S_daughters_deltaeta = TH1F('h_S_daughters_deltaeta','; #Delta#eta(K_{S}, ^{(} #bar{#Lambda} ^{)} ) ; 1/N_{ev} Events/0.1rad',60,-3,3)
	h_S_daughters_openingsangle = TH1F('h_S_daughters_openingsangle','; openings angle(K_{S},^{(} #bar{#Lambda} ^{)} ) (rad); 1/N_{ev} Events/0.1rad',35,0,3.5)
	h_S_daughters_DeltaR = TH1F('h_S_daughters_DeltaR','; #DeltaR(K_{S}, ^{(} #bar{#Lambda} ^{)} ); 1/N_{ev} Events/0.1#DeltaR',60,0.5,6.5)
	h_S_Ks_openingsangle = TH1F('h_S_Ks_openingsangle','; openings angle( ^{(} #bar{S} ^{)} ,K_{S}) (rad); 1/N_{ev} Events/0.1rad',20,0,2)
	h_S_Lambda_openingsangle = TH1F('h_S_Lambda_openingsangle','; openings angle( ^{(} #bar{S} ^{)},^{(} #bar{#Lambda} ^{)} ) (rad); 1/N_{ev} Events/0.1rad',20,0,2)

	h_S_eta = TH1F('h_S_eta','; #eta( ^{(} #bar{S} ^{)} ); 1/N_{ev} Events/0.5#eta',20,-5,5)
	h_Ks_eta = TH1F('h_Ks_eta','; #eta(K_{S}) ; 1/N_{ev} Events/0.5#eta',20,-5,5)
	h_Lambda_eta = TH1F('h_Lambda_eta','; #eta( ^{(} #bar{#Lambda} ^{)} ); 1/N_{ev} Events/0.5#eta',20,-5,5)

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

	h_S_error_lxy_interaction_vertex = TH1F('h_S_error_lxy_interaction_vertex','; #sigma(l_{0} interaction vertex ^{(}#bar{S} ^{)} w.r.t beampipe center) (cm); 1/N_{ev} Events/0.004mm',10,0,0.04)
	h_S_mass = TH1F('h_S_mass','; m_{ ^{(} #bar{S} ^{)} ,obs} (GeV); 1/N_{ev} Events/0.1GeV',100,-5,5)
	tprof_reweighing_factor = TProfile('tprof_reweighing_factor',';#eta ^{(}#bar{S} ^{)};reweighing factor',20,-5,5,0,50)

	nEntries = tree.GetEntries()
	print 'Number of entries in the tree: ', nEntries
	for i in range(0,nEntries):
		if(i==maxEvents):
			break
		if(i%1e4 == 0):
			print "reached entry: ", i
		tree.GetEntry(i)

		#need to reweigh the MC signal events, because the ones with high eta are more important, because they will pass more material
		reweighing_factor = config.calc_reweighing_factor(tree._S_eta[0],'MC-#bar{S}-Signal' in Legend[iTree])

		h_S_vz_interaction_vertex.Fill(tree._S_vz_interaction_vertex[0],reweighing_factor)
		h_S_lxy_interaction_vertex.Fill(tree._S_lxy_interaction_vertex_beampipeCenter[0],reweighing_factor)

		h_S_daughters_deltaphi.Fill(tree._S_daughters_deltaphi[0],reweighing_factor)
		h_S_daughters_deltaeta.Fill(tree._S_daughters_deltaeta[0],reweighing_factor)
		h_S_daughters_openingsangle.Fill(tree._S_daughters_openingsangle[0],reweighing_factor)
		h_S_daughters_DeltaR.Fill(tree._S_daughters_DeltaR[0],reweighing_factor)
		h_S_Ks_openingsangle.Fill(tree._S_Ks_openingsangle[0],reweighing_factor)
		h_S_Lambda_openingsangle.Fill(tree._S_Lambda_openingsangle[0],reweighing_factor)

		h_S_eta.Fill(tree._S_eta[0],reweighing_factor)
		h_Ks_eta.Fill(tree._Ks_eta[0],reweighing_factor)
		h_Lambda_eta.Fill(tree._Lambda_eta[0],reweighing_factor)

		h_S_dxy_over_lxy.Fill(tree._S_dxy_over_lxy[0],reweighing_factor)
		h_Ks_dxy_over_lxy.Fill(tree._Ks_dxy_over_lxy[0],reweighing_factor)
		h_Lambda_dxy_over_lxy.Fill(tree._Lambda_dxy_over_lxy[0],reweighing_factor)

		h_S_dz_min.Fill(tree._S_dz_min[0],reweighing_factor)
		h_Ks_dz_min.Fill(tree._Ks_dz_min[0],reweighing_factor)
		h_Lambda_dz_min.Fill(tree._Lambda_dz_min[0],reweighing_factor)

		h_Ks_pt.Fill(tree._Ks_pt[0],reweighing_factor)
		
		h_Lambda_lxy_decay_vertex.Fill(tree._Lambda_lxy_decay_vertex[0],reweighing_factor)
		h_S_chi2_ndof.Fill(tree._S_chi2_ndof[0],reweighing_factor)

		h_S_pz.Fill(tree._S_pz[0],reweighing_factor)

		h_S_error_lxy_interaction_vertex.Fill(tree._S_error_lxy_interaction_vertex_beampipeCenter[0],reweighing_factor)  
		h_S_mass.Fill(tree._S_mass[0],reweighing_factor)

		tprof_reweighing_factor.Fill(tree._S_eta[0],reweighing_factor)	

	TH1_l = [h_S_vz_interaction_vertex,h_S_lxy_interaction_vertex,h_S_daughters_deltaphi,h_S_daughters_deltaeta,h_S_daughters_openingsangle,h_S_daughters_DeltaR,h_S_Ks_openingsangle,h_S_Lambda_openingsangle,h_S_eta,h_Ks_eta,h_Lambda_eta,h_S_dxy_over_lxy,h_Ks_dxy_over_lxy,h_Lambda_dxy_over_lxy,h_S_dz_min,h_Ks_dz_min,h_Lambda_dz_min,h_Ks_pt,h_Lambda_lxy_decay_vertex,h_S_chi2_ndof,h_S_pz,h_S_error_lxy_interaction_vertex,h_S_mass,tprof_reweighing_factor]
	for h in TH1_l:
		h.SetDirectory(0) 
	TH1_ll.append(TH1_l)

	TH2_l = []
	for h in TH2_l:
		h.SetDirectory(0) 
	TH2_ll.append(TH2_l)

	iTree+=1

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
		if(h.GetName() != "tprof_reweighing_factor"):
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
