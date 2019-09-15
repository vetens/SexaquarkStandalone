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

inputFileLocation = '/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/FlatTree_Skimmed/CRAB_SimSexaq_trial17/crab_FlatTreeProducerTracking_trial17_14092019_v2/190914_134721/'
nInputFiles = 11
maxEvents = 5e99

plots_output_dir = "plots_Tracking_AntiS_specific/"

inFiles =  []

for i in range(1,nInputFiles+1):
	inFiles.append(TFile(inputFileLocation+'combined_FlatTree_Tracking_Skimmed_trial17_part'+str(i)+'.root','read'))	

fOut = TFile(plots_output_dir+'macro_combined_FlatTree_Tracking_Skimmed_trial17.root','RECREATE')

#histos to plot the number of granddaughters with >= 7 tracker hits versus a certain parameter of the antiS
tprof_numberGranddaughters_7hits_eta_antiS = TProfile('tprof_numberGranddaughters_7hits_eta_antiS',";#eta #bar{S};#final state particles with >= 7 numberOfTrackerLayers",20,-5,5) 
tprof_numberGranddaughters_7hits_vz_interaction_antiS = TProfile('tprof_numberGranddaughters_7hits_vz_interaction_antiS',";vz(beamspot) interaction vertex #bar{S} (cm);#final state particles with >= 7 numberOfTrackerLayers",50,-100,100) 
tprof_numberGranddaughters_7hits_lxy_interaction_antiS = TProfile('tprof_numberGranddaughters_7hits_lxy_interaction_antiS',";l_{0}(beamspot) interaction vertex #bar{S} (cm);#final state particles with >= 7 numberOfTrackerLayers",30,0,120) 
teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS',";#eta #bar{S};event fracion #geq 7 numberOfTrackerLayers each final state particle",100,-5,5) 
teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS',";#phi #bar{S};event fracion #geq 7 numberOfTrackerLayers each final state particle",100,-5,5) 
teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS',";vz(beamspot) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerLayers each final state particle",200,-100,100) 
teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS',";l_{0}(beamspot) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerLayers each final state particle",30,0,120) 
tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS = TProfile2D('tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS',";v_{z}(beamspot) interaction vertex #bar{S} (cm);l_{0}(beamspot) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerLayers each final state particle",240,-120,120,1200,0,120) 
teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS',";p_{t} #bar{S} (GeV);event fracion #geq 7 numberOfTrackerLayers each final state particle",100,0,10) 
teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS',";p_{z} #bar{S} (GeV);event fracion #geq 7 numberOfTrackerLayers each final state particle",100,0,100) 
t2_allTracksMoreThan6Hits_efficiency = TH2I("t2_allTracksMoreThan6Hits_efficiency",";all tracks >= 7 hits?;#bar{S} reconstructed;",2,-0.5,1.5,2,-0.5,1.5)

#reconstruction efficiencies of the antiS, based on the lxyz between the RECO and GEN interaction vertex
h_AntiS_deltaLInteractionVertexAntiSmin = TH1F("h_AntiS_deltaLInteractionVertexAntiSmin",";min #DeltaL_{xyz}(RECO #bar{S},simulated #bar{S}) interaction vertex;Events/mm",1000,0,100) 
teff_AntiS_RECO_eff_eta_antiS= TEfficiency('teff_AntiS_RECO_eff_eta_antiS',";#eta simulated #bar{S};#bar{S} RECO eff",100,-5,5) 
teff_AntiS_RECO_eff_vz_antiS= TEfficiency('teff_AntiS_RECO_eff_vz_antiS',";vz(beamspot) interaction vertex simulated #bar{S} (cm); #bar{S} RECO eff",200,-100,100) 
teff_AntiS_RECO_eff_lxy_antiS= TEfficiency('teff_AntiS_RECO_eff_lxy_antiS',";l_{0}(beamspot) interaction vertex simulated #bar{S} (cm); #bar{S} RECO eff",40,0,40) 
teff_AntiS_RECO_eff_pt_antiS= TEfficiency('teff_AntiS_RECO_eff_pt_antiS',";p_{t} simulated #bar{S} (GeV); #bar{S} RECO eff",100,0,10) 
teff_AntiS_RECO_eff_pz_antiS= TEfficiency('teff_AntiS_RECO_eff_pz_antiS',";p_{z} simulated #bar{S} (GeV); #bar{S} RECO eff",100,0,100) 

#reconstruction accuracy of the antiS
h1_AntiS_RECO_Acc_eta = TH1F("h1_AntiS_RECO_Acc_eta","; #eta(simulated #bar{S}) - #eta(reconstructed #bar{S});Events",80,-0.4,0.4)
h1_AntiS_RECO_Acc_phi = TH1F("h1_AntiS_RECO_Acc_phi","; #phi(simulated #bar{S}) - #phi(reconstructed #bar{S});Events",80,-0.4,0.4)
h1_AntiS_RECO_Acc_vz = TH1F("h1_AntiS_RECO_Acc_vz","; #v_{z}(interaction vertex, beamspot)(simulated #bar{S}) - #v_{z}(interaction vertex, beamspot)(reconstructed #bar{S});Events",100,-5,5)
h1_AntiS_RECO_Acc_lxy = TH1F("h1_AntiS_RECO_Acc_lxy","; #l_{0}(interaction vertex, beamspot)(simulated #bar{S}) - #l_{0}(interaction vertex, beamspot)(reconstructed #bar{S});Events",100,-5,5)
h1_AntiS_RECO_Acc_pt = TH1F("h1_AntiS_RECO_Acc_pt","; #p_{t}(simulated #bar{S}) - #p_{t}(reconstructed #bar{S});Events",200,-2,2)
h1_AntiS_RECO_Acc_pz = TH1F("h1_AntiS_RECO_Acc_pz","; #p_{z}(simulated #bar{S}) - #p_{z}(reconstructed #bar{S});Events",200,-2,2)

#invariant mass of the antiS
h_AntiS_massMinusNeutron = TH1F("h_AntiS_massMinusNeutron",";m_{#bar{S},obs} (GeV);Events/0.01GeV",2000,-10,10) 
h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda = TH2F("h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda","; m_{#bar{S},obs} (GeV); |#vec{p}_{K_{s}} + #vec{p}_{#bar{#Lambda}}| (GeV);Events/GeV^{2}",110,-5,6,60,0,40)

#list of list for histograms containging kinematics of the granddaughters of the AntiS
h_AntiSGrandDaughter_lxy_1 = TH1F("h_AntiSGrandDaughter_lxy_1",";l_{0}(creation vertex, beamspot) (cm);Events/cm",60,0,60) 
h_AntiSGrandDaughter_vz_1 = TH1F("h_AntiSGrandDaughter_vz_1",";v_{z}(creation vertex, beamspot) (cm);Events/cm",40,-200,200) 
h_AntiSGrandDaughter_dxy_1 = TH1F("h_AntiSGrandDaughter_dxy_1",";d_{0}(beamspot) (cm);Events/mm",40,-20,20) 
h_AntiSGrandDaughter_dz_1 = TH1F("h_AntiSGrandDaughter_dz_1",";d_{z}(beamspot) (cm);Events/cm",30,-60,60) 
h_AntiSGrandDaughter_pt_1 = TH1F("h_AntiSGrandDaughter_pt_1",";p_{t} (GeV);Events/0.1GeV",50,0,5) 
h_AntiSGrandDaughter_pz_1 = TH1F("h_AntiSGrandDaughter_pz_1",";p_{z} (GeV);Events/GeV",20,0,20) 

h_AntiSGrandDaughter_lxy_2 = TH1F("h_AntiSGrandDaughter_lxy_2",";l_{0}(creation vertex, beamspot) (cm);Events/cm",60,0,60) 
h_AntiSGrandDaughter_vz_2 = TH1F("h_AntiSGrandDaughter_vz_2",";v_{z}(creation vertex, beamspot) (cm);Events/cm",40,-200,200) 
h_AntiSGrandDaughter_dxy_2 = TH1F("h_AntiSGrandDaughter_dxy_2",";d_{0}(beamspot) (cm);Events/mm",40,-20,20)
h_AntiSGrandDaughter_dz_2 = TH1F("h_AntiSGrandDaughter_dz_2",";d_{z}(beamspot) (cm);Events/cm",30,-60,60)  
h_AntiSGrandDaughter_pt_2 = TH1F("h_AntiSGrandDaughter_pt_2",";p_{t} (GeV);Events/0.1GeV",50,0,5)
h_AntiSGrandDaughter_pz_2 = TH1F("h_AntiSGrandDaughter_pz_2",";p_{z} (GeV);Events/GeV",20,0,20) 

h_AntiSGrandDaughter_lxy_3 = TH1F("h_AntiSGrandDaughter_lxy_3",";l_{0}(creation vertex, beamspot) (cm);Events/cm",60,0,60)
h_AntiSGrandDaughter_vz_3 = TH1F("h_AntiSGrandDaughter_vz_3",";v_{z}(creation vertex, beamspot) (cm);Events/cm",40,-200,200)  
h_AntiSGrandDaughter_dxy_3 = TH1F("h_AntiSGrandDaughter_dxy_3",";d_{0}(beamspot) (cm);Events/mm",40,-20,20) 
h_AntiSGrandDaughter_dz_3 = TH1F("h_AntiSGrandDaughter_dz_3",";d_{z}(beamspot) (cm);Events/cm",30,-60,60)    
h_AntiSGrandDaughter_pt_3 = TH1F("h_AntiSGrandDaughter_pt_3",";p_{t} (GeV);Events/0.1GeV",50,0,5)
h_AntiSGrandDaughter_pz_3 = TH1F("h_AntiSGrandDaughter_pz_3",";p_{z} (GeV);Events/GeV",20,0,20)   

h_AntiSGrandDaughter_lxy_4 = TH1F("h_AntiSGrandDaughter_lxy_4",";l_{0}(creation vertex, beamspot) (cm);Events/cm",60,0,60)
h_AntiSGrandDaughter_vz_4 = TH1F("h_AntiSGrandDaughter_vz_4",";v_{z}(creation vertex, beamspot) (cm);Events/cm",40,-200,200) 
h_AntiSGrandDaughter_dxy_4 = TH1F("h_AntiSGrandDaughter_dxy_4",";d_{0}(beamspot) (cm);Events/mm",40,-20,20)
h_AntiSGrandDaughter_dz_4 = TH1F("h_AntiSGrandDaughter_dz_4",";d_{z}(beamspot) (cm);Events/cm",30,-60,60)   
h_AntiSGrandDaughter_pt_4 = TH1F("h_AntiSGrandDaughter_pt_4",";p_{t} (GeV);Events/0.1GeV",50,0,5)
h_AntiSGrandDaughter_pz_4 = TH1F("h_AntiSGrandDaughter_pz_4",";p_{z} (GeV);Events/GeV",20,0,20) 

ll_kinematics_granddaughters = [
[h_AntiSGrandDaughter_lxy_1,h_AntiSGrandDaughter_vz_1,h_AntiSGrandDaughter_dxy_1,h_AntiSGrandDaughter_dz_1,h_AntiSGrandDaughter_pt_1,h_AntiSGrandDaughter_pz_1],
[h_AntiSGrandDaughter_lxy_2,h_AntiSGrandDaughter_vz_2,h_AntiSGrandDaughter_dxy_2,h_AntiSGrandDaughter_dz_2,h_AntiSGrandDaughter_pt_2,h_AntiSGrandDaughter_pz_2],
[h_AntiSGrandDaughter_lxy_3,h_AntiSGrandDaughter_vz_3,h_AntiSGrandDaughter_dxy_3,h_AntiSGrandDaughter_dz_3,h_AntiSGrandDaughter_pt_3,h_AntiSGrandDaughter_pz_3],
[h_AntiSGrandDaughter_lxy_4,h_AntiSGrandDaughter_vz_4,h_AntiSGrandDaughter_dxy_4,h_AntiSGrandDaughter_dz_4,h_AntiSGrandDaughter_pt_4,h_AntiSGrandDaughter_pz_4],
]

for iFile, fIn in enumerate(inFiles,start = 1):
	print "Starting with inputFile: ", str(iFile) ,"/",str(len(inFiles)), ':', fIn.GetName()
	tree = fIn.Get('FlatTreeProducerTracking/FlatTreeTpsAntiS') 
	for i in range(0,tree.GetEntries()):
		tree.GetEntry(i)

		#my requirement for having reconstructed antiS
		antiSReconstructed = tree._tpsAntiS_deltaLInteractionVertexAntiSmin[0] < 0.5

		if(i>maxEvents):
			break

		if(i%10000 == 0):
			print 'reached track: ', i, ': ', float(i)/float(min(maxEvents,tree.GetEntries()))*100, '%'

		NGrandDaughtersWithTrackerLayerHitsLargerThan6 = 0
		for j in range(0,len(tree._tpsAntiS_type)):
			#now look at _tpAntiS_type from 3 to 6, that are the granddaughters:
			if(tree._tpsAntiS_type[j] >= 3 and tree._tpsAntiS_type[j] <= 6 ):
				#print "particle type ", tree._tpsAntiS_type[j], " has numberOfTrackerLayerHits = " , tree._tpsAntiS_numberOfTrackerLayers[j]
				if(tree._tpsAntiS_numberOfTrackerLayers[j] >= 7):
					NGrandDaughtersWithTrackerLayerHitsLargerThan6 += 1

		tprof_numberGranddaughters_7hits_eta_antiS.Fill(tree._tpsAntiS_eta[0],NGrandDaughtersWithTrackerLayerHitsLargerThan6)
		tprof_numberGranddaughters_7hits_vz_interaction_antiS.Fill(tree._tpsAntiS_vz_beamspot[1],NGrandDaughtersWithTrackerLayerHitsLargerThan6) #interaction vertex is the creation vertex of the daughter
		tprof_numberGranddaughters_7hits_lxy_interaction_antiS.Fill(tree._tpsAntiS_Lxy_beamspot[1],NGrandDaughtersWithTrackerLayerHitsLargerThan6) #interaction vertex is the creation vertex of the daughter

		boolNGrandDaughtersWithTrackerLayerHitsLargerThan6 =  NGrandDaughtersWithTrackerLayerHitsLargerThan6==4
		teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_eta[0])
		teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_phi[0])
		teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_vz_beamspot[1])
		teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_Lxy_beamspot[1])
		tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.Fill(tree._tpsAntiS_vz_beamspot[1],tree._tpsAntiS_Lxy_beamspot[1],boolNGrandDaughtersWithTrackerLayerHitsLargerThan6)
		teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_pt[0])
		teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_pz[0])

		#now instead of looking at the NGrandDaughtersWithTrackerLayerHitsLargerThan6 as a proxy for tracking efficiency now look at the real tracking efficiency for antiS
		teff_AntiS_RECO_eff_eta_antiS.Fill(antiSReconstructed,tree._tpsAntiS_eta[0])
		teff_AntiS_RECO_eff_vz_antiS.Fill(antiSReconstructed,tree._tpsAntiS_vz_beamspot[1])
		teff_AntiS_RECO_eff_lxy_antiS.Fill(antiSReconstructed,tree._tpsAntiS_Lxy_beamspot[1])
		teff_AntiS_RECO_eff_pt_antiS.Fill(antiSReconstructed,tree._tpsAntiS_pt[0])
		teff_AntiS_RECO_eff_pz_antiS.Fill(antiSReconstructed,tree._tpsAntiS_pz[0])
		#just a check of the tree: if a antiS got reconstructed then also all the daughters should be reconstructed. 
		if(antiSReconstructed):
			print tree._tpsAntiS_reconstructed[0], " ", tree._tpsAntiS_reconstructed[1], " ",tree._tpsAntiS_reconstructed[2], " ",tree._tpsAntiS_bestRECO_pt[3], " ",tree._tpsAntiS_bestRECO_pt[4], " ",tree._tpsAntiS_bestRECO_pt[5], " ",tree._tpsAntiS_bestRECO_pt[6]
		#accuracies:
		if(antiSReconstructed):
			h1_AntiS_RECO_Acc_eta.Fill(tree._tpsAntiS_eta[0]-tree._tpsAntiS_bestRECO_eta[0])
			h1_AntiS_RECO_Acc_phi.Fill(tree._tpsAntiS_phi[0]-tree._tpsAntiS_bestRECO_phi[0])
			h1_AntiS_RECO_Acc_vz.Fill(tree._tpsAntiS_vz_beamspot[1]-tree._tpsAntiS_bestRECO_vz_beamspot[0])
			h1_AntiS_RECO_Acc_lxy.Fill(tree._tpsAntiS_Lxy_beamspot[1]-tree._tpsAntiS_bestRECO_Lxy_beamspot[0])
			h1_AntiS_RECO_Acc_pt.Fill(tree._tpsAntiS_pt[0]-tree._tpsAntiS_bestRECO_pt[0])
			h1_AntiS_RECO_Acc_pz.Fill(tree._tpsAntiS_pz[0]-tree._tpsAntiS_bestRECO_pz[0])

		#check now how well boolNGrandDaughtersWithTrackerLayerHitsLargerThan6 is a proxy for the efficiency
		t2_allTracksMoreThan6Hits_efficiency.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,antiSReconstructed)
		h_AntiS_deltaLInteractionVertexAntiSmin.Fill(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[0])
			
		#now select tracks which have the antiS properly reconstructed and check for these events how the kinematics look like.
		if(antiSReconstructed):
			h_AntiS_massMinusNeutron.Fill(tree._tpsAntiS_bestRECO_massMinusNeutron[0])
			momentum_Ks_plus_Lambda = np.sqrt(  np.power(tree._tpsAntiS_pt[1]+tree._tpsAntiS_pt[2],2) + np.power(tree._tpsAntiS_pz[1]+tree._tpsAntiS_pz[2],2) )
			h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda.Fill(tree._tpsAntiS_bestRECO_massMinusNeutron[0],momentum_Ks_plus_Lambda)
			for j in range(0,len(tree._tpsAntiS_type)):
				#now look at _tpAntiS_type from 3 to 6, that are the granddaughters:
				if(tree._tpsAntiS_type[j] >= 3 and tree._tpsAntiS_type[j] <= 6 ):
					ll_kinematics_granddaughters[j-3][0].Fill(tree._tpsAntiS_Lxy_beamspot[j])	
					ll_kinematics_granddaughters[j-3][1].Fill(tree._tpsAntiS_vz_beamspot[j])	
					ll_kinematics_granddaughters[j-3][2].Fill(tree._tpsAntiS_dxy_beamspot[j])	
					ll_kinematics_granddaughters[j-3][3].Fill(tree._tpsAntiS_dz_beamspot[j])	
					ll_kinematics_granddaughters[j-3][4].Fill(tree._tpsAntiS_pt[j])	
					ll_kinematics_granddaughters[j-3][5].Fill(tree._tpsAntiS_pz[j])	


nHistos = len(ll_kinematics_granddaughters[0])
n_granddaughters = len(ll_kinematics_granddaughters)
leg_granddaugthers = ["K_{S}-#pi^{+}","K_{S}-#pi^{-}","#bar{#Lambda}-#pi^{+}","#bar{#Lambda}-#bar{p}"]
for i in range(0,nHistos):#each list contains a list of histograms. he histos need to be overlaid one list to the other,
	h = ll_kinematics_granddaughters[0][i]
	c_name = "c_"+h.GetName()
	c = TCanvas(c_name,"")
	legend = TLegend(0.8,0.85,0.99,0.99)
	for j in [2,3,0,1]: #start with the soft pion
		h = ll_kinematics_granddaughters[j][i]
		#if(h.GetSumw2N() == 0):
		#	h.Sumw2(kTRUE)
		#h.Scale(1./h.Integral(), "width");
		if j == 2:
			h.Draw("")
		else:
			h.Draw("same")
		h.SetLineColor(colours[j])
		h.SetFillColorAlpha(colours[j],0.15)
		h.SetLineWidth(2)
		h.SetMarkerStyle(22+j)
		h.SetMarkerColor(colours[j])
		h.SetStats(0)
		legend.AddEntry(h,leg_granddaugthers[j],"l")
	legend.Draw()
	CMS_lumi.CMS_lumi(c, 0, 11)
	c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
	c.Write()

RECO_AntiS_dir = fOut.mkdir("RECO_AntiS")
RECO_AntiS_dir.cd()
h_AntiS_massMinusNeutron.Write()
h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda.Write()


numberTrackerHits_dir = fOut.mkdir("numberTrackerHits")
numberTrackerHits_dir.cd()
tprof_numberGranddaughters_7hits_eta_antiS.Write()
tprof_numberGranddaughters_7hits_vz_interaction_antiS.Write()
tprof_numberGranddaughters_7hits_lxy_interaction_antiS.Write()
teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS.Write()
teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS.Write()
teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS.Write()

l_teff = [teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS,teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS]
for teff in l_teff:
	teff.SetDirectory(0)
for i in range(0,len(l_teff)):
	h = l_teff[i]
	c_name = "c_"+h.GetName()
	c = TCanvas(c_name,"")
	CMS_lumi.CMS_lumi(c, 0, 11)
	h.Draw()
	#c.Update()
	#h.GetPaintedGraph().GetYaxis().SetTitleSize(0.1)
	#c.Update()
	#c.GetHistogram().GetYaxis().SetTitleSize(0.1)
	#c.Update()
	c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
	c.Write()

teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS.Write()

tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.Write()
c_tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS= TCanvas(tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.GetName(),"");
c_tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.SetRightMargin(0.2) #make room for the tile of the z scale
#tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.GetZaxis().SetTitleSize(1)
tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.Draw("colz")
CMS_lumi.CMS_lumi(c_tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS, 0, 11)
c_tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.SetLogz()
c_tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.SaveAs(plots_output_dir+tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.GetName()+".pdf")
c_tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.Write()

teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS.Write()
teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS.Write()
t2_allTracksMoreThan6Hits_efficiency.Write()


antiS_RECO_eff_dir = fOut.mkdir("antiS_RECO_eff")
antiS_RECO_eff_dir.cd()
h_AntiS_deltaLInteractionVertexAntiSmin.Write()
teff_AntiS_RECO_eff_eta_antiS.Write()
teff_AntiS_RECO_eff_eta_antiS.GetPassedHistogram().Write()
teff_AntiS_RECO_eff_eta_antiS.GetTotalHistogram().Write()
teff_AntiS_RECO_eff_vz_antiS.Write()
teff_AntiS_RECO_eff_vz_antiS.GetPassedHistogram().Write()
teff_AntiS_RECO_eff_vz_antiS.GetTotalHistogram().Write()
teff_AntiS_RECO_eff_lxy_antiS.Write()
teff_AntiS_RECO_eff_lxy_antiS.GetPassedHistogram().Write()
teff_AntiS_RECO_eff_lxy_antiS.GetTotalHistogram().Write()
teff_AntiS_RECO_eff_pt_antiS.Write()
teff_AntiS_RECO_eff_pt_antiS.GetPassedHistogram().Write()
teff_AntiS_RECO_eff_pt_antiS.GetTotalHistogram().Write()
teff_AntiS_RECO_eff_pz_antiS.Write()
teff_AntiS_RECO_eff_pz_antiS.GetPassedHistogram().Write()
teff_AntiS_RECO_eff_pz_antiS.GetTotalHistogram().Write()


antiS_RECO_accuracy_dir = fOut.mkdir("antiS_RECO_accuracy")
antiS_RECO_accuracy_dir.cd()
h1_AntiS_RECO_Acc_eta.Write()
h1_AntiS_RECO_Acc_phi.Write()
h1_AntiS_RECO_Acc_vz.Write()
h1_AntiS_RECO_Acc_lxy.Write()
h1_AntiS_RECO_Acc_pt.Write()
h1_AntiS_RECO_Acc_pz.Write()

granddaughterkinematics_dir = fOut.mkdir("granddaughterkinematics")
granddaughterkinematics_dir.cd()
for l in ll_kinematics_granddaughters:
	for h in l:
		h.Write()

fOut.Close()

