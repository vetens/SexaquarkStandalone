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

inputFileLocation = '/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/FlatTree_Skimmed/CRAB_SimSexaq_trial17/crab_FlatTreeProducerTracking_trial17_14092019_v2/190914_134721/'
nInputFiles = 11
maxEvents = 1e99

plots_output_dir = "plots_Tracking_AntiS_specific/"

inFiles =  []

for i in range(1,nInputFiles+1):
	inFiles.append(TFile(inputFileLocation+'combined_FlatTree_Tracking_Skimmed_trial17_part'+str(i)+'.root','read'))	
inFiles = [TFile("/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerTracking/crab/hadd/combined_test_Be_4p3eta_FlatTreeProducerTracking_trial17.root",'read')]

fOut = TFile(plots_output_dir+'macro_combined_FlatTree_Tracking_Skimmed_trial17_acceptance.root','RECREATE')

def printProgress(i):
	if(i%10000 == 0):
		print 'reached track: ', i, ': ', float(i)/float(min(maxEvents,tree.GetEntries()))*100, '%'

#histos to plot the number of granddaughters with >= 7 tracker hits versus a certain parameter of the antiS
tprof_numberGranddaughters_7hits_eta_antiS = TProfile('tprof_numberGranddaughters_7hits_eta_antiS',";#eta #bar{S};#final state particles with >= 7 numberOfTrackerLayers",20,-5,5) 
tprof_numberGranddaughters_7hits_vz_interaction_antiS = TProfile('tprof_numberGranddaughters_7hits_vz_interaction_antiS',";vz(beamspot) interaction vertex #bar{S} (cm);#final state particles with >= 7 numberOfTrackerLayers",50,-100,100) 
tprof_numberGranddaughters_7hits_lxy_interaction_antiS = TProfile('tprof_numberGranddaughters_7hits_lxy_interaction_antiS',";l_{0}(beamspot) interaction vertex #bar{S} (cm);#final state particles with >= 7 numberOfTrackerLayers",30,0,120) 
teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS',";#eta #bar{S};event fracion #geq 7 numberOfTrackerLayers each final state particle",100,-5,5) 
teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS',";#phi #bar{S};event fracion #geq 7 numberOfTrackerLayers each final state particle",100,-5,5) 
teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS',";vz(beamspot) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerLayers each final state particle",200,-100,100) 
teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS',";l_{0}(beamspot) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerLayers each final state particle",30,0,120) 
tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS = TProfile2D('tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS',";v_{z}(beamspot) interaction vertex #bar{S} (cm);l_{0}(beamspot) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerLayers each final state particle",250,-125,125,1250,0,125) 
teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS',";p_{t} #bar{S} (GeV);event fracion #geq 7 numberOfTrackerLayers each final state particle",100,0,10) 
teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS= TEfficiency('teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS',";p_{z} #bar{S} (GeV);event fracion #geq 7 numberOfTrackerLayers each final state particle",100,0,100) 
#check the correlation between antiS being reconstructed and all tracks having larger than a certain amount of hits
h2_allTracksMoreThan4Hits_efficiency = TH2I("h2_allTracksMoreThan4Hits_efficiency",";all tracks >= 5 hits?;#bar{S} reconstructed;",2,-0.5,1.5,2,-0.5,1.5)
h2_allTracksMoreThan5Hits_efficiency = TH2I("h2_allTracksMoreThan5Hits_efficiency",";all tracks >= 6 hits?;#bar{S} reconstructed;",2,-0.5,1.5,2,-0.5,1.5)
h2_allTracksMoreThan6Hits_efficiency = TH2I("h2_allTracksMoreThan6Hits_efficiency",";all tracks >= 7 hits?;#bar{S} reconstructed;",2,-0.5,1.5,2,-0.5,1.5)


#reconstruction efficiencies of the antiS, based on the lxyz between the RECO and GEN interaction vertex
h_AntiS_deltaLInteractionVertexAntiSmin = TH1F("h_AntiS_deltaLInteractionVertexAntiSmin",";min #DeltaL_{xyz}(RECO #bar{S},simulated #bar{S}) interaction vertex;Events/mm",1000,0,100) 
teff_AntiS_RECO_eff_eta_antiS= TEfficiency('teff_AntiS_RECO_eff_eta_antiS',";#eta simulated #bar{S};",100,-5,5) 
teff_AntiS_RECO_eff_eta_antiS_antiS= TEfficiency('teff_AntiS_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};",100,-5,5) 
teff_AntiS_RECO_eff_vz_antiS= TEfficiency('teff_AntiS_RECO_eff_vz_antiS',";v_{z}(beamspot) interaction vertex simulated #bar{S} (cm);",40,-100,100) 
teff_AntiS_RECO_eff_lxy_antiS= TEfficiency('teff_AntiS_RECO_eff_lxy_antiS',";l_{0}(beamspot) interaction vertex simulated #bar{S} (cm);",60,0,60) 
teff_AntiS_RECO_eff_pt_antiS= TEfficiency('teff_AntiS_RECO_eff_pt_antiS',";p_{t} simulated #bar{S} (GeV); ",100,0,10) 
teff_AntiS_RECO_eff_pz_antiS= TEfficiency('teff_AntiS_RECO_eff_pz_antiS',";p_{z} simulated #bar{S} (GeV); ",80,0,80) 
teff_AntiS_RECO_eff_dxy_antiS= TEfficiency('teff_AntiS_RECO_eff_dxy_antiS',";d_{0}(beamspot) simulated #bar{S} (cm); ",40,-10,10) 
teff_AntiS_RECO_eff_dz_antiS= TEfficiency('teff_AntiS_RECO_eff_dz_antiS',";d_{z}(beamspot) simulated #bar{S} (cm); ",40,-100,100) 
teff_AntiS_RECO_eff_numberOfTrackerLayerHits_antiS= TEfficiency('teff_AntiS_RECO_eff_numberOfTrackerLayerHits_antiS',";numberOfTrackerLayerHits simulated #bar{S} ; ",40,0-0.5,40-0.5) 
#reco eff for the Ks
teff_Ks_RECO_eff_eta_antiS= TEfficiency('teff_Ks_RECO_eff_eta_antiS',";#eta simulated K_{S};",100,-5,5) 
teff_Ks_RECO_eff_eta_antiS_antiS= TEfficiency('teff_Ks_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};",100,-5,5) 
teff_Ks_RECO_eff_vz_antiS= TEfficiency('teff_Ks_RECO_eff_vz_antiS',";v_{z}(beamspot) decay vertex simulated K_{S} (cm);",40,-100,100) 
teff_Ks_RECO_eff_lxy_antiS= TEfficiency('teff_Ks_RECO_eff_lxy_antiS',";l_{0}(beamspot) decay vertex simulated K_{S} (cm);",60,0,60) 
teff_Ks_RECO_eff_pt_antiS= TEfficiency('teff_Ks_RECO_eff_pt_antiS',";p_{t} simulated K_{S} (GeV); ",100,0,10) 
teff_Ks_RECO_eff_pz_antiS= TEfficiency('teff_Ks_RECO_eff_pz_antiS',";p_{z} simulated K_{S} (GeV); ",80,0,80) 
teff_Ks_RECO_eff_dxy_antiS= TEfficiency('teff_Ks_RECO_eff_dxy_antiS',";d_{0}(beamspot) simulated K_{S} (cm); ",40,-10,10) 
teff_Ks_RECO_eff_dz_antiS= TEfficiency('teff_Ks_RECO_eff_dz_antiS',";d_{z}(beamspot) simulated K_{S} (cm); ",40,-100,100) 
teff_Ks_RECO_eff_numberOfTrackerLayerHits_antiS= TEfficiency('teff_Ks_RECO_eff_numberOfTrackerLayerHits_antiS',";numberOfTrackerLayerHits simulated K_{S} ; ",40,0-0.5,40-0.5) 
#reco eff for the AntiLambda
teff_AntiLambda_RECO_eff_eta_antiS= TEfficiency('teff_AntiLambda_RECO_eff_eta_antiS',";#eta simulated #bar{#Lambda};",100,-5,5) 
teff_AntiLambda_RECO_eff_eta_antiS_antiS= TEfficiency('teff_AntiLambda_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};",100,-5,5) 
teff_AntiLambda_RECO_eff_vz_antiS= TEfficiency('teff_AntiLambda_RECO_eff_vz_antiS',";v_{z}(beamspot) decay vertex simulated #bar{#Lambda} (cm);",40,-100,100) 
teff_AntiLambda_RECO_eff_lxy_antiS= TEfficiency('teff_AntiLambda_RECO_eff_lxy_antiS',";l_{0}(beamspot) decay vertex simulated #bar{#Lambda} (cm);",60,0,60) 
teff_AntiLambda_RECO_eff_pt_antiS= TEfficiency('teff_AntiLambda_RECO_eff_pt_antiS',";p_{t} simulated #bar{#Lambda} (GeV); ",100,0,10) 
teff_AntiLambda_RECO_eff_pz_antiS= TEfficiency('teff_AntiLambda_RECO_eff_pz_antiS',";p_{z} simulated #bar{#Lambda} (GeV); ",80,0,80) 
teff_AntiLambda_RECO_eff_dxy_antiS= TEfficiency('teff_AntiLambda_RECO_eff_dxy_antiS',";d_{0}(beamspot) simulated #bar{#Lambda} (cm); ",40,-10,10) 
teff_AntiLambda_RECO_eff_dz_antiS= TEfficiency('teff_AntiLambda_RECO_eff_dz_antiS',";d_{z}(beamspot) simulated #bar{#Lambda} (cm); ",40,-100,100) 
teff_AntiLambda_RECO_eff_numberOfTrackerLayerHits_antiS= TEfficiency('teff_AntiLambda_RECO_eff_numberOfTrackerLayerHits_antiS',";numberOfTrackerLayerHits simulated #bar{#Lambda} (cm); ",40,0-0.5,40-0.5) 
#reco eff for the Ks daughter 0 and 1
teff_Ksdaugthers_RECO_eff_eta_antiS= TEfficiency('teff_Ksdaugthers_RECO_eff_eta_antiS',";#eta simulated K_{S} daughters;",100,-5,5) 
teff_Ksdaughters_RECO_eff_eta_antiS_antiS= TEfficiency('teff_Ksdaughters_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};",100,-5,5) 
teff_Ksdaugthers_RECO_eff_vz_antiS= TEfficiency('teff_Ksdaugthers_RECO_eff_vz_antiS',";v_{z}(beamspot) creation vertex simulated K_{S} daughters (cm);",40,-100,100) 
teff_Ksdaugthers_RECO_eff_lxy_antiS= TEfficiency('teff_Ksdaugthers_RECO_eff_lxy_antiS',";l_{0}(beamspot) creation vertex simulated K_{S} daughters (cm);",60,0,60) 
teff_Ksdaugthers_RECO_eff_pt_antiS= TEfficiency('teff_Ksdaugthers_RECO_eff_pt_antiS',";p_{t} simulated K_{S} daughters (GeV); ",100,0,10) 
teff_Ksdaugthers_RECO_eff_pz_antiS= TEfficiency('teff_Ksdaugthers_RECO_eff_pz_antiS',";p_{z} simulated K_{S} daughters (GeV); ",80,0,80) 
teff_Ksdaugthers_RECO_eff_dxy_antiS= TEfficiency('teff_Ksdaugthers_RECO_eff_dxy_antiS',";d_{0}(beamspot) simulated K_{S} daughters (cm); ",40,-10,10) 
teff_Ksdaugthers_RECO_eff_dz_antiS= TEfficiency('teff_Ksdaugthers_RECO_eff_dz_antiS',";d_{z}(beamspot) simulated K_{S} daughters (cm); ",40,-100,100) 
teff_Ksdaugthers_RECO_eff_numberOfTrackerLayerHits_antiS= TEfficiency('teff_Ksdaugthers_RECO_eff_numberOfTrackerLayerHits_antiS',";numberOfTrackerLayerHits simulated K_{S} daughters; ",40,0-0.5,40-0.5) 
#reco eff for the soft pion from the AntiLambda
teff_AntiLambdaPion_RECO_eff_eta_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_eta_antiS',";#eta simulated #bar{#Lambda}-#pi^{+};",100,-5,5) 
teff_AntiLambdaPion_RECO_eff_eta_antiS_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};",100,-5,5) 
teff_AntiLambdaPion_RECO_eff_vz_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_vz_antiS',";v_{z}(beamspot) creation vertex simulated #bar{#Lambda}-#pi^{+} (cm);",40,-100,100) 
teff_AntiLambdaPion_RECO_eff_lxy_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_lxy_antiS',";l_{0}(beamspot) creation vertex simulated #bar{#Lambda}-#pi^{+} (cm);",60,0,60) 
teff_AntiLambdaPion_RECO_eff_pt_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_pt_antiS',";p_{t} simulated #bar{#Lambda}-#pi^{+} (GeV); ",100,0,10) 
teff_AntiLambdaPion_RECO_eff_pz_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_pz_antiS',";p_{z} simulated #bar{#Lambda}-#pi^{+} (GeV); ",80,0,80) 
teff_AntiLambdaPion_RECO_eff_dxy_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_dxy_antiS',";d_{0}(beamspot) simulated #bar{#Lambda}-#pi^{+} (cm); ",40,-10,10) 
teff_AntiLambdaPion_RECO_eff_dz_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_dz_antiS',";d_{z}(beamspot) simulated #bar{#Lambda}-#pi^{+} (cm); ",40,-100,100) 
teff_AntiLambdaPion_RECO_eff_numberOfTrackerLayerHits_antiS= TEfficiency('teff_AntiLambdaPion_RECO_eff_numberOfTrackerLayerHits_antiS',";numberOfTrackerLayerHits simulated #bar{#Lambda}-#pi^{+}; ",40,0-0.5,40-0.5) 
#reco eff for the antiproton from the AntiLambda
teff_AntiLambdaAntiProton_RECO_eff_eta_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_eta_antiS',";#eta simulated #bar{#Lambda}-#pi^{+};",100,-5,5) 
teff_AntiLambdaAntiProton_RECO_eff_eta_antiS_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};",100,-5,5) 
teff_AntiLambdaAntiProton_RECO_eff_vz_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_vz_antiS',";v_{z}(beamspot) creation vertex simulated #bar{#Lambda}-#pi^{+} (cm);",40,-100,100) 
teff_AntiLambdaAntiProton_RECO_eff_lxy_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_lxy_antiS',";l_{0}(beamspot) creation vertex simulated #bar{#Lambda}-#pi^{+} (cm);",60,0,60) 
teff_AntiLambdaAntiProton_RECO_eff_pt_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_pt_antiS',";p_{t} simulated #bar{#Lambda}-#pi^{+} (GeV); ",100,0,10) 
teff_AntiLambdaAntiProton_RECO_eff_pz_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_pz_antiS',";p_{z} simulated #bar{#Lambda}-#pi^{+} (GeV); ",80,0,80) 
teff_AntiLambdaAntiProton_RECO_eff_dxy_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_dxy_antiS',";d_{0}(beamspot) simulated #bar{#Lambda}-#pi^{+} (cm); ",40,-10,10) 
teff_AntiLambdaAntiProton_RECO_eff_dz_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_dz_antiS',";d_{z}(beamspot) simulated #bar{#Lambda}-#pi^{+} (cm); ",40,-100,100) 
teff_AntiLambdaAntiProton_RECO_eff_numberOfTrackerLayerHits_antiS= TEfficiency('teff_AntiLambdaAntiProton_RECO_eff_numberOfTrackerLayerHits_antiS',";numberOfTrackerLayerHits simulated #bar{#Lambda}-#pi^{+}; ",40,0-0.5,40-0.5) 

ll_efficiencies = [
[teff_AntiS_RECO_eff_eta_antiS,teff_AntiS_RECO_eff_eta_antiS_antiS,teff_AntiS_RECO_eff_vz_antiS,teff_AntiS_RECO_eff_lxy_antiS,teff_AntiS_RECO_eff_pt_antiS,teff_AntiS_RECO_eff_pz_antiS,teff_AntiS_RECO_eff_dxy_antiS,teff_AntiS_RECO_eff_dz_antiS,teff_AntiS_RECO_eff_numberOfTrackerLayerHits_antiS],
[teff_Ks_RECO_eff_eta_antiS,teff_Ks_RECO_eff_eta_antiS_antiS,teff_Ks_RECO_eff_vz_antiS,teff_Ks_RECO_eff_lxy_antiS,teff_Ks_RECO_eff_pt_antiS,teff_Ks_RECO_eff_pz_antiS,teff_Ks_RECO_eff_dxy_antiS,teff_Ks_RECO_eff_dz_antiS,teff_Ks_RECO_eff_numberOfTrackerLayerHits_antiS],
[teff_AntiLambda_RECO_eff_eta_antiS,teff_AntiLambda_RECO_eff_eta_antiS_antiS,teff_AntiLambda_RECO_eff_vz_antiS,teff_AntiLambda_RECO_eff_lxy_antiS,teff_AntiLambda_RECO_eff_pt_antiS,teff_AntiLambda_RECO_eff_pz_antiS,teff_AntiLambda_RECO_eff_dxy_antiS,teff_AntiLambda_RECO_eff_dz_antiS,teff_AntiLambda_RECO_eff_numberOfTrackerLayerHits_antiS],
[teff_Ksdaugthers_RECO_eff_eta_antiS,teff_Ksdaughters_RECO_eff_eta_antiS_antiS,teff_Ksdaugthers_RECO_eff_vz_antiS,teff_Ksdaugthers_RECO_eff_lxy_antiS,teff_Ksdaugthers_RECO_eff_pt_antiS,teff_Ksdaugthers_RECO_eff_pz_antiS,teff_Ksdaugthers_RECO_eff_dxy_antiS,teff_Ksdaugthers_RECO_eff_dz_antiS,teff_Ksdaugthers_RECO_eff_numberOfTrackerLayerHits_antiS],
[teff_AntiLambdaPion_RECO_eff_eta_antiS,teff_AntiLambdaPion_RECO_eff_eta_antiS_antiS,teff_AntiLambdaPion_RECO_eff_vz_antiS,teff_AntiLambdaPion_RECO_eff_lxy_antiS,teff_AntiLambdaPion_RECO_eff_pt_antiS,teff_AntiLambdaPion_RECO_eff_pz_antiS,teff_AntiLambdaPion_RECO_eff_dxy_antiS,teff_AntiLambdaPion_RECO_eff_dz_antiS,teff_AntiLambdaPion_RECO_eff_numberOfTrackerLayerHits_antiS],
[teff_AntiLambdaAntiProton_RECO_eff_eta_antiS,teff_AntiLambdaAntiProton_RECO_eff_eta_antiS_antiS,teff_AntiLambdaAntiProton_RECO_eff_vz_antiS,teff_AntiLambdaAntiProton_RECO_eff_lxy_antiS,teff_AntiLambdaAntiProton_RECO_eff_pt_antiS,teff_AntiLambdaAntiProton_RECO_eff_pz_antiS,teff_AntiLambdaAntiProton_RECO_eff_dxy_antiS,teff_AntiLambdaAntiProton_RECO_eff_dz_antiS,teff_AntiLambdaAntiProton_RECO_eff_numberOfTrackerLayerHits_antiS]
]

for l in ll_efficiencies:
	for h in l:
		h.Clone()
		h.SetDirectory(0)

#create the same list of histograms as in ll_efficiencies but fill these only for antiS which have all granddaughters within tracker acceptance, where this is based on having >=7 hits on all tracks.
ll_efficiencies_acceptance = []
for l in range(0,len(ll_efficiencies)):
	ll_efficiencies_acceptance.append([])
	for h in range(0,len(ll_efficiencies[0])):
		histo = ll_efficiencies[l][h]
		histoClone = histo.Clone()
		histoClone.SetName(histo.GetName()+"_acceptance")
		histoClone.SetDirectory(0)
		ll_efficiencies_acceptance[l].append(histoClone)

for l in ll_efficiencies_acceptance:
	for h in l:
		h.Clone()
		h.SetDirectory(0)

#reconstruction accuracy of the antiS
h1_AntiS_RECO_Acc_eta = TH1F("h1_AntiS_RECO_Acc_eta","; #eta_{sim #bar{S}} - #eta_{reco #bar{S}};Events/0.01#eta",80,-0.4,0.4)
h1_AntiS_RECO_Acc_phi = TH1F("h1_AntiS_RECO_Acc_phi","; #phi_{sim #bar{S}} - #phi_{reco #bar{S}} (rad);Events/0.01rad",80,-0.4,0.4)
h1_AntiS_RECO_Acc_vz = TH1F("h1_AntiS_RECO_Acc_vz",";  v_{z, sim #bar{S}}(int vertex, beamspot) - v_{z, reco #bar{S}}(int vertex, beamspot) (cm);Events/mm",10,-1,1)
h1_AntiS_RECO_Acc_lxy = TH1F("h1_AntiS_RECO_Acc_lxy","; l_{0, sim #bar{S}}(int vertex, beamspot) - l_{0, reco #bar{S}}(int vertex, beamspot) (cm);Events/mm",10,-1,1)
h1_AntiS_RECO_Acc_pt = TH1F("h1_AntiS_RECO_Acc_pt","; p_{t, sim #bar{S}} - p_{t, reco #bar{S}} (GeV);Events/0.01GeV",400,-2,2)
h1_AntiS_RECO_Acc_pz = TH1F("h1_AntiS_RECO_Acc_pz","; p_{z, sim #bar{S}} - p_{z, reco #bar{S}} (GeV);Events/0.1GeV",40,-2,2)

#invariant mass of the antiS
h_AntiS_massMinusNeutron = TH1F("h_AntiS_massMinusNeutron",";m_{#bar{S},obs} RECO (GeV);Events/0.1GeV",180,-6,12) 
h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda = TH2F("h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda","; m_{#bar{S},obs} (GeV); |#vec{p}_{K_{s}, RECO} + #vec{p}_{#bar{#Lambda}, RECO}| (GeV);Events/GeV^{2}",110,-5,6,60,0,40)

#list of list for histograms containging kinematics of the granddaughters of the AntiS, for AntiS which got reconstructed. These plots are necessary because they tell where you need to be sure that your tracking is properly described by MC. 
h_AntiSGrandDaughter_lxy_1 = TH1F("h_AntiSGrandDaughter_lxy_1",";Simulated track l_{0}(creation vertex, beamspot) (cm);Events/cm",60,0,60) 
h_AntiSGrandDaughter_vz_1 = TH1F("h_AntiSGrandDaughter_vz_1",";Simulated track v_{z}(creation vertex, beamspot) (cm);Events/10cm",40,-200,200) 
h_AntiSGrandDaughter_dxy_1 = TH1F("h_AntiSGrandDaughter_dxy_1",";Simulated track d_{0}(beamspot) (cm);Events/cm",40,-20,20) 
h_AntiSGrandDaughter_dz_1 = TH1F("h_AntiSGrandDaughter_dz_1",";Simulated track d_{z}(beamspot) (cm);Events/4cm",30,-60,60) 
h_AntiSGrandDaughter_pt_1 = TH1F("h_AntiSGrandDaughter_pt_1",";Simulated track p_{t} (GeV);Events/0.1GeV",50,0,5) 
h_AntiSGrandDaughter_pz_1 = TH1F("h_AntiSGrandDaughter_pz_1",";Simulated track p_{z} (GeV);Events/GeV",20,0,20) 

h_AntiSGrandDaughter_lxy_3 = TH1F("h_AntiSGrandDaughter_lxy_3",";Simulated track l_{0}(creation vertex, beamspot) (cm);Events/cm",60,0,60)
h_AntiSGrandDaughter_vz_3 = TH1F("h_AntiSGrandDaughter_vz_3",";Simulated track v_{z}(creation vertex, beamspot) (cm);Events/10cm",40,-200,200)  
h_AntiSGrandDaughter_dxy_3 = TH1F("h_AntiSGrandDaughter_dxy_3",";Simulated track d_{0}(beamspot) (cm);Events/cm",40,-20,20) 
h_AntiSGrandDaughter_dz_3 = TH1F("h_AntiSGrandDaughter_dz_3",";Simulated track d_{z}(beamspot) (cm);Events/4cm",30,-60,60)    
h_AntiSGrandDaughter_pt_3 = TH1F("h_AntiSGrandDaughter_pt_3",";Simulated track p_{t} (GeV);Events/0.1GeV",50,0,5)
h_AntiSGrandDaughter_pz_3 = TH1F("h_AntiSGrandDaughter_pz_3",";Simulated track p_{z} (GeV);Events/GeV",20,0,20)   

h_AntiSGrandDaughter_lxy_4 = TH1F("h_AntiSGrandDaughter_lxy_4",";Simulated track l_{0}(creation vertex, beamspot) (cm);Events/cm",60,0,60)
h_AntiSGrandDaughter_vz_4 = TH1F("h_AntiSGrandDaughter_vz_4",";Simulated track v_{z}(creation vertex, beamspot) (cm);Events/10cm",40,-200,200) 
h_AntiSGrandDaughter_dxy_4 = TH1F("h_AntiSGrandDaughter_dxy_4",";Simulated track d_{0}(beamspot) (cm);Events/cm",40,-20,20)
h_AntiSGrandDaughter_dz_4 = TH1F("h_AntiSGrandDaughter_dz_4",";Simulated track d_{z}(beamspot) (cm);Events/4cm",30,-60,60)   
h_AntiSGrandDaughter_pt_4 = TH1F("h_AntiSGrandDaughter_pt_4",";Simulated track p_{t} (GeV);Events/0.1GeV",50,0,5)
h_AntiSGrandDaughter_pz_4 = TH1F("h_AntiSGrandDaughter_pz_4",";Simulated track p_{z} (GeV);Events/GeV",20,0,20) 

ll_kinematics_granddaughters_of_RECO_AntiS = [
[h_AntiSGrandDaughter_lxy_1,h_AntiSGrandDaughter_vz_1,h_AntiSGrandDaughter_dxy_1,h_AntiSGrandDaughter_dz_1,h_AntiSGrandDaughter_pt_1,h_AntiSGrandDaughter_pz_1],
[h_AntiSGrandDaughter_lxy_3,h_AntiSGrandDaughter_vz_3,h_AntiSGrandDaughter_dxy_3,h_AntiSGrandDaughter_dz_3,h_AntiSGrandDaughter_pt_3,h_AntiSGrandDaughter_pz_3],
[h_AntiSGrandDaughter_lxy_4,h_AntiSGrandDaughter_vz_4,h_AntiSGrandDaughter_dxy_4,h_AntiSGrandDaughter_dz_4,h_AntiSGrandDaughter_pt_4,h_AntiSGrandDaughter_pz_4],
]

#a list with counters for the reconstructed particles, so there are 7 entries for each of the 7 particles
nAntiS = 0.
nTotalReconstructed = [0.,0.,0.,0.,0.,0.,0.]
#and count separately the reconstruction efficiency of the Ks, antiLambda and antiS if their respective daughters got reconstructed.
nKsRECOIfBothDaughtersReco = 0.
nKsTOTALIfBothDaughtersReco = 0.
nAntiLambdaRECOIfBothDaughtersReco = 0.
nAntiLambdaTOTALIfBothDaughtersReco = 0.
nAntiSRECOIfBothDaughtersReco = 0.
nAntiSTOTALIfBothDaughtersReco = 0.

nAntiSWithAllGranddaughtersMoreThan6Hits = 0.
for iFile, fIn in enumerate(inFiles,start = 1):
	print "Starting with inputFile: ", str(iFile) ,"/",str(len(inFiles)), ':', fIn.GetName()
	tree = fIn.Get('FlatTreeProducerTracking/FlatTreeTpsAntiS') 
	for i in range(0,tree.GetEntries()):
		tree.GetEntry(i)

		nAntiS+=1

		#my requirement for having reconstructed antiS: based on the 3D distance between the GEN and RECO interaction vertex of the antiS.
		antiSReconstructed = tree._tpsAntiS_deltaLInteractionVertexAntiSmin[0] < 0.5

		if(i>maxEvents):
			break
		
		printProgress(i)

                #just a check of the tree: if a antiS got reconstructed then also all the daughters should be reconstructed.  
                if(antiSReconstructed): 
                        print tree._tpsAntiS_reconstructed[0], " ", tree._tpsAntiS_reconstructed[1], " ",tree._tpsAntiS_reconstructed[2], " ",tree._tpsAntiS_reconstructed[3], " ",tree._tpsAntiS_reconstructed[4], " ",tree._tpsAntiS_reconstructed[5], " ",tree._tpsAntiS_reconstructed[6]

		#count the reconstructed particles with the requirement that there daughters got reconstructed
		if(tree._tpsAntiS_reconstructed[3] == 1 and tree._tpsAntiS_reconstructed[4] == 1):
			if(tree._tpsAntiS_reconstructed[1] == 1):
				nKsRECOIfBothDaughtersReco += 1
			nKsTOTALIfBothDaughtersReco +=1
		if(tree._tpsAntiS_reconstructed[5] == 1 and tree._tpsAntiS_reconstructed[6] == 1):
			if(tree._tpsAntiS_reconstructed[2] == 1):
				nAntiLambdaRECOIfBothDaughtersReco += 1
			nAntiLambdaTOTALIfBothDaughtersReco += 1
		if(tree._tpsAntiS_reconstructed[1] == 1 and tree._tpsAntiS_reconstructed[2] == 1):
			if(antiSReconstructed):
				nAntiSRECOIfBothDaughtersReco += 1
			nAntiSTOTALIfBothDaughtersReco += 1

		#count the number of granddaughters of this antiS which have >= a certain amount of hits
		NGrandDaughtersWithTrackerLayerHitsLargerThan4 = 0
		NGrandDaughtersWithTrackerLayerHitsLargerThan5 = 0
		NGrandDaughtersWithTrackerLayerHitsLargerThan6 = 0
		for j in range(0,len(tree._tpsAntiS_type)):
			#now look at _tpAntiS_type from 3 to 6, this are the granddaughters:
			if(tree._tpsAntiS_type[j] >= 3 and tree._tpsAntiS_type[j] <= 6 ):
				#print "particle type ", tree._tpsAntiS_type[j], " has numberOfTrackerLayerHits = " , tree._tpsAntiS_numberOfTrackerLayers[j]
				if(tree._tpsAntiS_numberOfTrackerLayers[j] >= 5):
					NGrandDaughtersWithTrackerLayerHitsLargerThan4 += 1
				if(tree._tpsAntiS_numberOfTrackerLayers[j] >= 6):
					NGrandDaughtersWithTrackerLayerHitsLargerThan5 += 1
				if(tree._tpsAntiS_numberOfTrackerLayers[j] >= 7):
					NGrandDaughtersWithTrackerLayerHitsLargerThan6 += 1

		tprof_numberGranddaughters_7hits_eta_antiS.Fill(tree._tpsAntiS_eta[0],NGrandDaughtersWithTrackerLayerHitsLargerThan6)
		tprof_numberGranddaughters_7hits_vz_interaction_antiS.Fill(tree._tpsAntiS_vz_beamspot[1],NGrandDaughtersWithTrackerLayerHitsLargerThan6) #interaction vtx is the creation vtx of the daughter
		tprof_numberGranddaughters_7hits_lxy_interaction_antiS.Fill(tree._tpsAntiS_Lxy_beamspot[1],NGrandDaughtersWithTrackerLayerHitsLargerThan6) #interaction vtx is the creation vtx of the daughter

		boolNGrandDaughtersWithTrackerLayerHitsLargerThan4 =  NGrandDaughtersWithTrackerLayerHitsLargerThan4==4
		boolNGrandDaughtersWithTrackerLayerHitsLargerThan5 =  NGrandDaughtersWithTrackerLayerHitsLargerThan5==4
		boolNGrandDaughtersWithTrackerLayerHitsLargerThan6 =  NGrandDaughtersWithTrackerLayerHitsLargerThan6==4

		#plots having the fraction of all events with all granddaughters with >= 7 hits in function of a kinematic variable of the antiS
		teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_eta[0])
		teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_phi[0])
		teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_vz_beamspot[1])
		teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_Lxy_beamspot[1])
		tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.Fill(tree._tpsAntiS_vz_beamspot[1],tree._tpsAntiS_Lxy_beamspot[1],boolNGrandDaughtersWithTrackerLayerHitsLargerThan6)
		teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_pt[0])
		teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,tree._tpsAntiS_pz[0])
		#counter for all the antiS which have all granddaughters producing >= 7 hits
		if(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6):
			nAntiSWithAllGranddaughtersMoreThan6Hits+=1

		#now instead of looking at the NGrandDaughtersWithTrackerLayerHitsLargerThan6 as a proxy for tracking efficiency now look at the real tracking efficiency for antiS and its daughters
		for i in range(0,7):#for all of the 7 particles
			index = i
			if (i == 3 or i == 4):#take the two daughters of the Ks together
				index = 3
			elif(i == 5):#because I took the previous one together
				index = 4
			elif(i == 6):
				index = 5

			particleReconstructed = int(tree._tpsAntiS_reconstructed[i])
			if(index == 0): #for the antiS use the definition on top to say if it was reconstructed, not the one from the tree
				particleReconstructed = antiSReconstructed

			#make a global count of what got reconstructed
			if(particleReconstructed): #here I still count the Ks daughters separately
				nTotalReconstructed[i] += 1
			
			ll_efficiencies[index][0].Fill(particleReconstructed,tree._tpsAntiS_eta[i])
			ll_efficiencies[index][1].Fill(particleReconstructed,tree._tpsAntiS_eta[0]) #always in function of eta of the original antiS
			if(i == 0): #if antiS use interaction vertex as point of reference
				ll_efficiencies[index][2].Fill(particleReconstructed,tree._tpsAntiS_vz_beamspot[1])
				ll_efficiencies[index][3].Fill(particleReconstructed,tree._tpsAntiS_Lxy_beamspot[1])
			if(i == 1): #if Ks use the  decay vertex as point of reference
				ll_efficiencies[index][2].Fill(particleReconstructed,tree._tpsAntiS_vz_beamspot[3])
				ll_efficiencies[index][3].Fill(particleReconstructed,tree._tpsAntiS_Lxy_beamspot[3])
			if(i == 2): #if AntiL use the  decay vertex as point of reference
				ll_efficiencies[index][2].Fill(particleReconstructed,tree._tpsAntiS_vz_beamspot[5])
				ll_efficiencies[index][3].Fill(particleReconstructed,tree._tpsAntiS_Lxy_beamspot[5])
			else: #for the tracks use there creation vertex as point of reference
				ll_efficiencies[index][2].Fill(particleReconstructed,tree._tpsAntiS_vz_beamspot[i])
				ll_efficiencies[index][3].Fill(particleReconstructed,tree._tpsAntiS_Lxy_beamspot[i])
			ll_efficiencies[index][4].Fill(particleReconstructed,tree._tpsAntiS_pt[i])
			ll_efficiencies[index][5].Fill(particleReconstructed,tree._tpsAntiS_pz[i])
			ll_efficiencies[index][6].Fill(particleReconstructed,tree._tpsAntiS_dxy_beamspot[i])
			ll_efficiencies[index][7].Fill(particleReconstructed,tree._tpsAntiS_dz_beamspot[i])
			ll_efficiencies[index][8].Fill(particleReconstructed,tree._tpsAntiS_numberOfTrackerLayers[i])
			
			#it is almost a requirement for the antiS to be reconstructed that boolNGrandDaughtersWithTrackerLayerHitsLargerThan6 so now look for the events which have all final state particles falling in the acceptance what are the remaining inneficiencies
			if(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6):
				ll_efficiencies_acceptance[index][0].Fill(particleReconstructed,tree._tpsAntiS_eta[i])
				ll_efficiencies_acceptance[index][1].Fill(particleReconstructed,tree._tpsAntiS_eta[0]) #always in function of eta of the original antiS
				if(i == 0): #if antiS use interaction vertex as point of reference
					ll_efficiencies_acceptance[index][2].Fill(particleReconstructed,tree._tpsAntiS_vz_beamspot[1])
					ll_efficiencies_acceptance[index][3].Fill(particleReconstructed,tree._tpsAntiS_Lxy_beamspot[1])
				if(i == 1): #if Ks use the  decay vertex as point of reference
					ll_efficiencies_acceptance[index][2].Fill(particleReconstructed,tree._tpsAntiS_vz_beamspot[3])
					ll_efficiencies_acceptance[index][3].Fill(particleReconstructed,tree._tpsAntiS_Lxy_beamspot[3])
				if(i == 2): #if AntiL use the  decay vertex as point of reference
					ll_efficiencies_acceptance[index][2].Fill(particleReconstructed,tree._tpsAntiS_vz_beamspot[5])
					ll_efficiencies_acceptance[index][3].Fill(particleReconstructed,tree._tpsAntiS_Lxy_beamspot[5])
				else: #for the tracks use there creation vertex as point of reference
					ll_efficiencies_acceptance[index][2].Fill(particleReconstructed,tree._tpsAntiS_vz_beamspot[i])
					ll_efficiencies_acceptance[index][3].Fill(particleReconstructed,tree._tpsAntiS_Lxy_beamspot[i])
				ll_efficiencies_acceptance[index][4].Fill(particleReconstructed,tree._tpsAntiS_pt[i])
				ll_efficiencies_acceptance[index][5].Fill(particleReconstructed,tree._tpsAntiS_pz[i])
				ll_efficiencies_acceptance[index][6].Fill(particleReconstructed,tree._tpsAntiS_dxy_beamspot[i])
				ll_efficiencies_acceptance[index][7].Fill(particleReconstructed,tree._tpsAntiS_dz_beamspot[i])
				ll_efficiencies_acceptance[index][8].Fill(particleReconstructed,tree._tpsAntiS_numberOfTrackerLayers[i])
					
		#accuracies:
		if(antiSReconstructed):
			h1_AntiS_RECO_Acc_eta.Fill(tree._tpsAntiS_eta[0]-tree._tpsAntiS_bestRECO_eta[0])
			h1_AntiS_RECO_Acc_phi.Fill(tree._tpsAntiS_phi[0]-tree._tpsAntiS_bestRECO_phi[0])
			h1_AntiS_RECO_Acc_vz.Fill(tree._tpsAntiS_vz_beamspot[1]-tree._tpsAntiS_bestRECO_vz_beamspot[0])
			h1_AntiS_RECO_Acc_lxy.Fill(tree._tpsAntiS_Lxy_beamspot[1]-tree._tpsAntiS_bestRECO_Lxy_beamspot[0])
			h1_AntiS_RECO_Acc_pt.Fill(tree._tpsAntiS_pt[0]-tree._tpsAntiS_bestRECO_pt[0])
			h1_AntiS_RECO_Acc_pz.Fill(tree._tpsAntiS_pz[0]-tree._tpsAntiS_bestRECO_pz[0])

		#check now how well boolNGrandDaughtersWithTrackerLayerHitsLargerThan6 is a proxy for the efficiency
		h2_allTracksMoreThan4Hits_efficiency.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan4,antiSReconstructed)
		h2_allTracksMoreThan5Hits_efficiency.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan5,antiSReconstructed)
		h2_allTracksMoreThan6Hits_efficiency.Fill(boolNGrandDaughtersWithTrackerLayerHitsLargerThan6,antiSReconstructed)
		h_AntiS_deltaLInteractionVertexAntiSmin.Fill(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[0])
			
		#now select tracks which have the antiS properly reconstructed and check for these events how the kinematics look like of the final state particles.
		if(antiSReconstructed):
			h_AntiS_massMinusNeutron.Fill(tree._tpsAntiS_bestRECO_massMinusNeutron[0])
			momentum_Ks_plus_Lambda = np.sqrt(  np.power(tree._tpsAntiS_bestRECO_pt[1]+tree._tpsAntiS_bestRECO_pt[2],2) + np.power(tree._tpsAntiS_bestRECO_pz[1]+tree._tpsAntiS_bestRECO_pz[2],2) )
			h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda.Fill(tree._tpsAntiS_bestRECO_massMinusNeutron[0],momentum_Ks_plus_Lambda)
			for j in range(0,len(tree._tpsAntiS_type)):
				#now look at _tpAntiS_type from 3 to 6, that are the granddaughters:
				if(tree._tpsAntiS_type[j] >= 3 and tree._tpsAntiS_type[j] <= 6 ):
					index = 0
					if(tree._tpsAntiS_type[j]==3 or tree._tpsAntiS_type[j]==4): #both daughters of the Ks can go in the same plot
						index = 0
					elif(tree._tpsAntiS_type[j]==5):
						index = 1
					elif(tree._tpsAntiS_type[j]==6):
						index = 2
					ll_kinematics_granddaughters_of_RECO_AntiS[index][0].Fill(tree._tpsAntiS_Lxy_beamspot[j])	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][1].Fill(tree._tpsAntiS_vz_beamspot[j])	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][2].Fill(tree._tpsAntiS_dxy_beamspot[j])	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][3].Fill(tree._tpsAntiS_dz_beamspot[j])	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][4].Fill(tree._tpsAntiS_pt[j])	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][5].Fill(tree._tpsAntiS_pz[j])	


#plot kinematics of final state particles for which the antiS got reconstructed
granddaughterkinematics_dir_of_RECO_AntiS = fOut.mkdir("granddaughterkinematics_of_RECO_AntiS") 
granddaughterkinematics_dir_of_RECO_AntiS.cd() 
for l in ll_kinematics_granddaughters_of_RECO_AntiS: 
        for h in l: 
                h.Write()

nHistos = len(ll_kinematics_granddaughters_of_RECO_AntiS[0])
n_granddaughters = len(ll_kinematics_granddaughters_of_RECO_AntiS)
leg_granddaugthers = ["K_{S} daughters","#bar{#Lambda}-#pi^{+}","#bar{#Lambda}-#bar{p}"]
for i in range(0,nHistos):#each list contains a list of histograms. the histos need to be overlaid one list to the other,
	c_name = "c_"+ll_kinematics_granddaughters_of_RECO_AntiS[0][i].GetName()
	c = TCanvas(c_name,"")
	legend = TLegend(0.8,0.85,0.99,0.99)
	for j in [1,2,0]: #first the soft pion of the antiLambda, then the antiproton and then the Ks daughters (this is also the order you are using at GEN level) 
		h = ll_kinematics_granddaughters_of_RECO_AntiS[j][i]
		if(h.GetSumw2N() == 0):
			h.Sumw2(kTRUE)
		h.Scale(1./h.Integral(), "width");
		if j == 1:
			h.Draw("PCE1")
		else:
			h.Draw("PCE1same")
		h.SetLineColor(colours[j])
		h.SetFillColorAlpha(colours[j],0.15)
		h.SetLineWidth(2)
		h.SetMarkerStyle(22+j)
		h.SetMarkerColor(colours[j])
		h.SetStats(0)
		legend.AddEntry(h,leg_granddaugthers[j],"lep")
	legend.Draw()
	CMS_lumi.CMS_lumi(c, 0, 11)
	c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
	c.Write()

#some kinematics of the AntiS
RECO_AntiS_dir = fOut.mkdir("RECO_AntiS")
RECO_AntiS_dir.cd()
h_AntiS_massMinusNeutron.Write()
c_name = "c_"+h_AntiS_massMinusNeutron.GetName()
c = TCanvas(c_name,"")
if(h_AntiS_massMinusNeutron.GetSumw2N() == 0):
	h.Sumw2(kTRUE)
h_AntiS_massMinusNeutron.Scale(1./h.Integral(), "width");
h_AntiS_massMinusNeutron.Draw("PCE1")
h_AntiS_massMinusNeutron.SetLineColor(colours[0])
h_AntiS_massMinusNeutron.SetLineWidth(2)
h_AntiS_massMinusNeutron.SetMarkerStyle(22+0)
h_AntiS_massMinusNeutron.SetMarkerColor(colours[0])
h_AntiS_massMinusNeutron.SetStats(0)
CMS_lumi.CMS_lumi(c, 0, 11)
c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
c.Write()
h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda.Write()

#for the number of tracker hits study
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
h2_allTracksMoreThan4Hits_efficiency.Write()
h2_allTracksMoreThan5Hits_efficiency.Write()
h2_allTracksMoreThan6Hits_efficiency.Write()

#The reconstruction efficiencies of all the particles
#list of directory names:
l_dir_names = ["RECO_eff_antiS", "RECO_eff_Ks", "RECO_eff_antiLambda", "RECO_eff_Ks_daughters", "RECO_eff_antiLambda_pion", "RECO_eff_antiLambda_antiProton"]
for il, l in enumerate(ll_efficiencies,start = 0):
	directory = fOut.mkdir(l_dir_names[il])
	directory.cd()
	for h in l:
		h.Write()
		h.GetPassedHistogram().Write()
		h.GetTotalHistogram().Write()
		hc_PassedHistogram = h.GetPassedHistogram().Clone()
	
		c_name = "c_"+hc_PassedHistogram.GetName()
		c = TCanvas(c_name,"")
		if(hc_PassedHistogram.GetSumw2N() == 0):
			hc_PassedHistogram.Sumw2(kTRUE)
		hc_PassedHistogram.Scale(1./hc_PassedHistogram.Integral(), "width");

		if("pt" in hc_PassedHistogram.GetName()):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/0.1GeV")	
		elif("pz" in hc_PassedHistogram.GetName()):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/1GeV")	
		elif("eta" in hc_PassedHistogram.GetName()):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/0.1#eta")	
		elif("lxy" in hc_PassedHistogram.GetName()):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/cm")	
		elif("vz" in hc_PassedHistogram.GetName()):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/5cm")	
		hc_PassedHistogram.SetLineColor(colours[0])
		hc_PassedHistogram.SetLineWidth(2)
		hc_PassedHistogram.SetMarkerStyle(22+0)
		hc_PassedHistogram.SetMarkerColor(colours[0])
		hc_PassedHistogram.Draw("PCE1")
		CMS_lumi.CMS_lumi(c, 0, 11)
		c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")

#reconstruction efficiencies of all the particles which have all final state particles within the accepatance
l_dir_names = ["RECO_eff_antiS_acceptance", "RECO_eff_Ks_acceptance", "RECO_eff_antiLambda_acceptance", "RECO_eff_Ks_daughters_acceptance", "RECO_eff_antiLambda_pion_acceptance", "RECO_eff_antiLambda_antiProton_acceptance"]
#l_particle = ["#bar{S}","Ks","Lambda","KsDaughters","AntiLambdaPion","AntiLambdaAntiProton"]
for il, l in enumerate(ll_efficiencies_acceptance,start = 0):
	directory = fOut.mkdir(l_dir_names[il])
	directory.cd()
	for h in l:
		h.Write()
		hc = h.Clone()
		hc.GetPassedHistogram().Write()
		hc.GetTotalHistogram().Write()

		h_TotalHistogram = hc.GetTotalHistogram()
		hc_TotalHistogram = h_TotalHistogram.Clone()
		#now make for each particle and for each variable a single plot which contains on the same plot the normalised distribution of a certain kin variable and the reconstruction eff wrt that kin variable
		c_name = "c_"+hc.GetName()+"_distribution_and_efficiency"
		c = TCanvas(c_name,"")
		legend = TLegend(0.7,0.9,0.99,0.99)
		if(hc_TotalHistogram.GetSumw2N() == 0):
			hc_TotalHistogram.Sumw2(kTRUE)
		hc_TotalHistogram.Scale(1./hc_TotalHistogram.Integral(), "width");
		
		hc.SetLineColor(colours[1])
		hc.SetLineWidth(2)
		hc.SetMarkerStyle(22+1)
		hc.SetMarkerColor(colours[1])
		legend.AddEntry(hc,"Reconstruction efficiency","lep")
		hc.Draw("")

		#scale the distribution a bit to make it visible on the scale of the efficiency
		hc_TotalHistogram.SetLineColor(colours[0])
		hc_TotalHistogram.SetLineWidth(2)
		hc_TotalHistogram.SetMarkerStyle(22+0)
		hc_TotalHistogram.SetMarkerColor(colours[0])
		hc_TotalHistogram.SetStats(0)
		#hc_TotalHistogram.DrawNormalized("samePCE1")
		#scale the distribution a bit to make it visible on the scale of the efficiency
		if("pz" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale()
			legend.AddEntry(hc_TotalHistogram,"Distribution","lep")
		elif("eta" in hc_TotalHistogram.GetName() or "dxy" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(1)
                        legend.AddEntry(hc_TotalHistogram,"Distribution","lep")
		elif("vz" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(30)
                        legend.AddEntry(hc_TotalHistogram,"Distribution x30","lep")
		elif("dz" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(20)
                        legend.AddEntry(hc_TotalHistogram,"Distribution x20","lep")
		elif("lxy" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(8)
                        legend.AddEntry(hc_TotalHistogram,"Distribution x8","lep")
		elif("pt" in hc_TotalHistogram.GetName() ):
			hc_TotalHistogram.Scale(0.5)
                        legend.AddEntry(hc_TotalHistogram,"Distribution x0.25","lep")
		else:	
			hc_TotalHistogram.Scale(50)
			legend.AddEntry(hc_TotalHistogram,"Distribution x50","lep")
		hc_TotalHistogram.Draw("samePCE1")

		legend.Draw()
		CMS_lumi.CMS_lumi(c, 0, 11)
		c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
		c.Write()
	
#the accuracy plots
antiS_RECO_accuracy_dir = fOut.mkdir("antiS_RECO_accuracy")
antiS_RECO_accuracy_dir.cd()

l_h = [h1_AntiS_RECO_Acc_eta,h1_AntiS_RECO_Acc_phi,h1_AntiS_RECO_Acc_vz,h1_AntiS_RECO_Acc_lxy,h1_AntiS_RECO_Acc_pt,h1_AntiS_RECO_Acc_pz]
for h in l_h:
	h.SetDirectory(0)
for i in range(0,len(l_h)):
	h = l_h[i]
	c_name = "c_"+h.GetName()
	c = TCanvas(c_name,"")
	CMS_lumi.CMS_lumi(c, 0, 11)
	if(h.GetSumw2N() == 0):
		h.Sumw2(kTRUE)
	h.Scale(1./h.Integral(), "width");
	
	h.SetLineColor(colours[0])
	h.SetLineWidth(2)
	h.SetMarkerStyle(22+0)
	h.SetMarkerColor(colours[0])
	h.Fit("gaus")
	h.Draw("")
	c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
	c.Write()


fOut.Close()




#print some conclusions:

print "Fraction of antiS with all granddaughters generating >= 7 tracker hits=: ", nAntiSWithAllGranddaughtersMoreThan6Hits, "/",nAntiS, " = ", nAntiSWithAllGranddaughtersMoreThan6Hits/nAntiS 

#summary of reconstruction efficiencies:
particles = ["antiS","Ks","AntiLambda","Ks-piplus","Ks-pineg","AntiLambda-pion","AntiLambda-antiproton"]
print "###################################################################################################"
print "########################################RECO EFF###################################################"
print "###################################################################################################"
for i_par, nominator in enumerate(nTotalReconstructed,start = 0):
	if i_par == 1 or i_par == 3:
		print "------------------------------------------"
	print particles[i_par],": ", nominator,"/",nAntiS, " = ", nominator/nAntiS, "       --> so 1 out of: " , nAntiS/nominator

print "------------------------------------------"
print "Ks if both daughters got reconstructed: ", nKsRECOIfBothDaughtersReco,"/",nKsTOTALIfBothDaughtersReco," = ", nKsRECOIfBothDaughtersReco/nKsTOTALIfBothDaughtersReco, "              --> so 1 out of :", nKsTOTALIfBothDaughtersReco/nKsRECOIfBothDaughtersReco
print "AntiLambda if both daughters got reconstructed: ", nAntiLambdaRECOIfBothDaughtersReco,"/",nAntiLambdaTOTALIfBothDaughtersReco," = ", nAntiLambdaRECOIfBothDaughtersReco/nAntiLambdaTOTALIfBothDaughtersReco, "              --> so 1 out of :", nAntiLambdaTOTALIfBothDaughtersReco/nAntiLambdaRECOIfBothDaughtersReco
print "AntiS if both daughters got reconstructed: ", nAntiSRECOIfBothDaughtersReco,"/",nAntiSTOTALIfBothDaughtersReco," = ", nAntiSRECOIfBothDaughtersReco/nAntiSTOTALIfBothDaughtersReco, "              --> so 1 out of :", nAntiSTOTALIfBothDaughtersReco/nAntiSRECOIfBothDaughtersReco

