#macro to plot the kinemaics of the generated Sbar of different masses

#from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, gROOT
from ROOT import *
import numpy as np
import sys
#sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
sys.path.append('/afs/cern.ch/work/w/wvetens/Sexaquarks/CMSSW_10_2_26/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

colours = [1,2,4,35,38,41]

#fIn = [
#TFile('/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/crmc/Sexaquark_13TeV_trial20/FlatTree/FlatTree_GEN_trial20.root', 'read'),
#TFile('/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/crmc/Sexaquark_13TeV_trial19/FlatTree/FlatTree_GEN_trial19.root', 'read'),
#TFile('/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/crmc/Sexaquark_13TeV_trial17/FlatTree/FlatTree_GEN_trial17.root', 'read'),
#TFile('/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/crmc/Sexaquark_13TeV_trial18/FlatTree/FlatTree_GEN_trial18.root', 'read'),
#]

#mass = ["1.2","1.5","1.8","2.1"]
fIn = [
TFile('/afs/cern.ch/work/w/wvetens/Sexaquarks/CMSSW_10_2_26/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerGEN/output.root', 'read')
]
mass = ["1.8"]

plots_output_dir = "plots_GEN/"

maxEvents = 1e6

TH1_ll = [] #list of list of 1D histos 
TH2_ll = [] #list of list of 2D histos

iFile = 0
for f in fIn:
	tree = f.Get('FlatTreeProducerGEN/FlatTreeGENLevel') 

	#h_antiS_pt = TH1F('h_antiS_pt','; #bar{S} p_{T} (GeV/c); 1/N_{ev} Events/0.1GeV/c',100,0,10)
	#h_antiS_pz = TH1F('h_antiS_pz','; #bar{S} |p_{z}| (GeV/c); 1/N_{ev} Events/1GeV/c',80,0,80)
	#h_antiS_eta = TH1F('h_antiS_eta','; #bar{S} #eta ; 1/N_{ev} Events/0.1#eta',160,-8,8)
	#h_antiS_vz = TH1F('h_antiS_vz','; #bar{S} absolute creation vertex z (cm) ; 1/N_{ev} Events/cm',40,-20,20)
	#h_antiS_vz_interaction_vertex = TH1F('h_antiS_vz_interaction_vertex','; #bar{S} beampipe crossing location absolute z (cm) ; 1/N_{ev} Events/5cm',100,-250,250)
	#h_antiS_lxy = TH1F('h_antiS_lxy','; #bar{S} absolute creation vertex l_{0} (cm) ; 1/N_{ev} Events/mm',20,-1,1)

	#h_antiS_eta_pt = TH2F('h_antiS_eta_pt',';#bar{S} #eta; #bar{S} p_{T} (GeV/c); 1/N_{ev} Events/0.1#eta/0.1GeV/c',160,-8,8,100,0,10)
	#h_antiS_eta_pz = TH2F('h_antiS_eta_pz',';#bar{S} #eta; #bar{S} |p_{z}| (GeV/c); 1/N_{ev} Events/0.1#eta/1GeV/c',160,-8,8,100,0,100)
	h_antiS_pt = TH1F('h_antiS_pt','; #bar{S} p_{T} (GeV/c); 1/N_{ev} Events/0.1GeV/c',100,0,10)
	h_antiS_pz = TH1F('h_antiS_pz','; #bar{S} |p_{z}| (GeV/c); 1/N_{ev} Events/1GeV/c',80,0,80)
	h_antiS_eta = TH1F('h_antiS_eta','; #bar{S} #eta ; 1/N_{ev} Events/0.1#eta',160,-8,8)
	h_antiS_vz = TH1F('h_antiS_vz','; #bar{S} absolute creation vertex z (cm) ; 1/N_{ev} Events/cm',40,-20,20)
	h_antiS_vz_interaction_vertex = TH1F('h_antiS_vz_interaction_vertex','; #bar{S} beampipe crossing location absolute z (cm) ; 1/N_{ev} Events/5cm',100,-250,250)
	h_antiS_lxy = TH1F('h_antiS_lxy','; #bar{S} absolute creation vertex l_{0} (cm) ; 1/N_{ev} Events/mm',20,-1,1)

	h_antiS_eta_pt = TH2F('h_antiS_eta_pt',';#bar{S} #eta; #bar{S} p_{T} (GeV/c); 1/N_{ev} Events/0.1#eta/0.1GeV/c',160,-8,8,100,0,10)
	h_antiS_eta_pz = TH2F('h_antiS_eta_pz',';#bar{S} #eta; #bar{S} |p_{z}| (GeV/c); 1/N_{ev} Events/0.1#eta/1GeV/c',160,-8,8,100,0,100)

	nEntries = tree.GetEntries()
	print 'Number of entries in the tree: ', nEntries
	for i in range(0,nEntries):
		if(i==maxEvents):
			break
		if(i%1e4 == 0):
			print "reached entry: ", i
		tree.GetEntry(i)
		h_antiS_pt.Fill(tree._S_pt[0])
		h_antiS_pz.Fill(tree._S_pz[0])
		h_antiS_eta.Fill(tree._S_eta[0])
		h_antiS_vz.Fill(tree._S_vz[0])
		vz_interaction_antiS = 2.21/np.tan( 2*np.arctan( np.exp(-tree._S_eta[0]) ) )
		h_antiS_vz_interaction_vertex.Fill(vz_interaction_antiS)
		h_antiS_lxy.Fill( np.sqrt( np.power(tree._S_vx[0],2) + np.power(tree._S_vy[0],2) ))

		h_antiS_eta_pt.Fill(tree._S_eta[0],tree._S_pt[0])
		h_antiS_eta_pz.Fill(tree._S_eta[0],abs(tree._S_pz[0]))

	TH1_l = [h_antiS_pt,h_antiS_pz,h_antiS_eta,h_antiS_vz,h_antiS_vz_interaction_vertex,h_antiS_lxy]
	for h in TH1_l:
		h.SetDirectory(0) 
	TH1_ll.append(TH1_l)

	TH2_l = [h_antiS_eta_pt, h_antiS_eta_pz]
	for h in TH2_l:
		h.SetDirectory(0) 
	TH2_ll.append(TH2_l)

	iFile+=1

fOut = TFile('macro_FlatTree_GEN_trial17.root','RECREATE')

nHistos = len(TH1_ll[0])
nMasses = len(TH1_ll)

for i in range(0,nHistos):#each list contains a list of histograms. Each list represents a specific mass. The histos need to be overlaid one list to the other
	h = TH1_ll[0][i]
	c_name = "c_"+h.GetName()
	c = TCanvas(c_name,"")
	legend = TLegend(0.7,0.85,0.99,0.99)
	for j in range(0,nMasses):
		h = TH1_ll[j][i]
		if(h.GetSumw2N() == 0):
			h.Sumw2(kTRUE)
		h.Scale(1./h.Integral(), "width");
		if j == 0:
			h.Draw("PCE1")
		else:
			h.Draw("same")
		h.SetLineColor(colours[j])
		h.SetMarkerStyle(22+j)
		h.SetMarkerColor(colours[j])
		h.SetStats(0)
		legend.AddEntry(h,"#bar{S} mass = "+ mass[j] +" GeV/c^{2}  ","lep")
	legend.Draw()
	CMS_lumi.CMS_lumi(c, 0, 11)
        gPad.Update()
	c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
	c.Write()


iMass = 0
for l in TH2_ll:
	for h in l:
		c_name = "c_"+h.GetName()+"_"+mass[iMass]
		c = TCanvas(c_name,"");
		c.SetRightMargin(0.2) #make room for the tile of the z scale
		if(h.GetSumw2N() == 0):
			h.Sumw2(kTRUE)
		h.Scale(1./h.Integral(), "width");
		h.Draw("colz")
		h.SetStats(0)
		CMS_lumi.CMS_lumi(c, 0, 11)
                gPad.Update()
		c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
		c.Write()
	iMass += 1


fOut.Write()
fOut.Close()
