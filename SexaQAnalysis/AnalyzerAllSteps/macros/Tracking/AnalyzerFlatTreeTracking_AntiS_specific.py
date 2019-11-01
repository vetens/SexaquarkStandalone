import numpy as np
from ROOT import *
import random
import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle 
import collections

sys.path.insert(1, '/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA')
import configBDT as config
config_dict = config.config_dict

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

colours = [1,2,4,35,38,41]

maxEvents1 = 2e99
maxEvents2 = 2e99


plots_output_dir = "plots_Tracking_AntiS_specific/"

inFiles = [TFile("file:/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerTracking/test_FlatTreeTracking_Step1_Step2_Skimming_FlatTree_trial17_1p8GeV_17102019_v1.root",'read')]


fOut = TFile(plots_output_dir+'macro_combined_FlatTree_Tracking_Skimmed_trial17.root','RECREATE')

def printProgress(i):
	if(i%10000 == 0):
		print 'reached track: ', i, ': ', float(i)/float(min(maxEvents2,tree.GetEntries()))*100, '%'

def FillHistosEfficiency(tree,ll_efficiencies,i_particle,index,weightFactor,i):

	tree.GetEntry(i)

	ll_efficiencies[index][0].Fill(tree._tpsAntiS_eta[i_particle],weightFactor)
	ll_efficiencies[index][1].Fill(tree._tpsAntiS_eta[0],weightFactor) #always in function of eta of the original antiS
	if(i_particle == 0): #if antiS use interaction vertex as point of reference
		ll_efficiencies[index][2].Fill(tree._tpsAntiS_vz[1],weightFactor)
		ll_efficiencies[index][3].Fill(tree._tpsAntiS_Lxy_beampipeCenter[1],weightFactor)
	elif(i_particle == 1): #if Ks use the  decay vertex as point of reference
		ll_efficiencies[index][2].Fill(tree._tpsAntiS_vz[3],weightFactor)
		ll_efficiencies[index][3].Fill(tree._tpsAntiS_Lxy_beamspot[3],weightFactor)
	elif(i_particle == 2): #if AntiL use the  decay vertex as point of reference
		ll_efficiencies[index][2].Fill(tree._tpsAntiS_vz[5],weightFactor)
		ll_efficiencies[index][3].Fill(tree._tpsAntiS_Lxy_beamspot[5],weightFactor)
	elif(i_particle == 3 or i_particle == 4 or i_particle == 5 or i_particle == 6): #for the tracks use there creation vertex as point of reference
		ll_efficiencies[index][2].Fill(tree._tpsAntiS_vz[i_particle],weightFactor)
		ll_efficiencies[index][3].Fill(tree._tpsAntiS_Lxy_beamspot[i_particle],weightFactor)
	ll_efficiencies[index][4].Fill(tree._tpsAntiS_pt[i_particle],weightFactor)
	ll_efficiencies[index][5].Fill(tree._tpsAntiS_pz[i_particle],weightFactor)
	ll_efficiencies[index][6].Fill(np.sqrt(tree._tpsAntiS_pz[i_particle]*tree._tpsAntiS_pz[i_particle]+tree._tpsAntiS_pt[i_particle]*tree._tpsAntiS_pt[i_particle]),weightFactor)
	ll_efficiencies[index][7].Fill(tree._tpsAntiS_dxyTrack_beamspot[i_particle],weightFactor)
	ll_efficiencies[index][8].Fill(tree._tpsAntiS_dzTrack_beamspot[i_particle],weightFactor)
	ll_efficiencies[index][9].Fill(tree._tpsAntiS_numberOfTrackerHits[i_particle],weightFactor)
	ll_efficiencies[index][10].Fill(tree._tpsAntiS_phi[i_particle],weightFactor)


#first make a few plots on the PV distribution to check the reweighing
PV_dir = fOut.mkdir("PV")
PV_dir.cd()
h_PVz_non_weighed = TH1F("h_PVz_non_weighed",";PV absolute z;Events/cm",600,-30,30)
h_PVz_PUweighed = TH1F("h_PVz_PUweighed",";PV absolute z;Events/cm",600,-30,30)
h_weight_parameter_PVz = TH1F("h_weight_parameter_PVz",";weight parameter PVz;",200,0,2)

h_nPV_MC_non_weighed = TH1I('h_nPV_MC_non_weighed','; # valid PV; Events',60,-0.5,59.5)
h_nPV_MC_PUweighed = TH1I('h_nPV_MC_PUweighed','; # valid PV; Events',60,-0.5,59.5)
h_weight_parameter_nPV = TH1F("h_weight_parameter_nPV",";weight parameter PVz;",200,0,2)

h2_nPV_vzPV_MC = TH2F('h2_nPV_vzPV_MC','; #PV; absolute v_{z} PV (cm);  Events',60,-0.5,59.5,600,-30,30)

for iFile, fIn in enumerate(inFiles,start = 1):
        print "Starting with inputFile: ", str(iFile) ,"/",str(len(inFiles)), ':', fIn.GetName()
	treePV = fIn.Get('FlatTreeProducerTracking/FlatTreePV')
	for i in range(0,treePV.GetEntries()):
		if(i>maxEvents1):
			break

		treePV.GetEntry(i)

		#draw a random vertex (you should be able to reweigh using a random vertex as a random vertex is representative for the event and this is an event by event reweighing) from the distribution and do nPV reweighing on this vertex
		randomVertexIndex = random.randint(0,len(treePV._goodPV_weightPU)-1)
		sum_weight = 0.
		for j in range(0,len(treePV._goodPV_weightPU)):
			weightPU = treePV._goodPV_weightPU[j]
			sum_weight = weightPU + sum_weight
			#print 'weightPU: ',weightPU

			h_PVz_non_weighed.Fill(treePV._goodPVzPOG[j],1.)
			h_PVz_PUweighed.Fill(treePV._goodPVzPOG[j],weightPU)
			h2_nPV_vzPV_MC.Fill(len(treePV._goodPV_weightPU),treePV._goodPVzPOG[j])	
			h_weight_parameter_PVz.Fill(weightPU)
			if j==randomVertexIndex: random_weight = weightPU

			#dic_PVz_weight.update({treePV._goodPVzPOG[j]:weightPU})
		

		#random weight
		h_nPV_MC_non_weighed.Fill(len(treePV._goodPV_weightPU),1.)
		h_nPV_MC_PUweighed.Fill(len(treePV._goodPV_weightPU),random_weight*len(treePV._goodPV_weightPU)/18.479)#/17.88*1.5) #the last two factors are just to get the overall distribution in the same ballpark as the unweighed one for MC. These are constant weights so they just represent an overall scaling to the histogram and in the end only relative weights are important
		h_weight_parameter_nPV.Fill(random_weight)


h_PVz_non_weighed.Write()
h_PVz_PUweighed.Write()
h_weight_parameter_PVz.Write()

h_nPV_MC_non_weighed.Write()
h_nPV_MC_PUweighed.Write()
h_weight_parameter_nPV.Write()

h2_nPV_vzPV_MC.Write()

#histos to plot the number of granddaughters with >= 7 tracker hits versus a certain parameter of the antiS
tprof_numberGranddaughters_7hits_eta_antiS = TProfile('tprof_numberGranddaughters_7hits_eta_antiS',";#eta #bar{S};#final state particles with >= 7 numberOfTrackerHits",20,-5,5) 
tprof_numberGranddaughters_7hits_vz_interaction_antiS = TProfile('tprof_numberGranddaughters_7hits_vz_interaction_antiS',";vz(bs) interaction vertex #bar{S} (cm);#final state particles with >= 7 numberOfTrackerHits",50,-100,100) 
tprof_numberGranddaughters_7hits_lxy_interaction_antiS = TProfile('tprof_numberGranddaughters_7hits_lxy_interaction_antiS',";l_{0}(bs) interaction vertex #bar{S} (cm);#final state particles with >= 7 numberOfTrackerHits",30,0,120) 
teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS= TProfile('teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS',";#eta #bar{S};event fracion #geq 7 numberOfTrackerHits each final state particle",100,-5,5) 
teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS= TProfile('teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS',";#phi #bar{S};event fracion #geq 7 numberOfTrackerHits each final state particle",100,-5,5) 
teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS= TProfile('teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS',";vz(bs) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerHits each final state particle",200,-100,100) 
teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS= TProfile('teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS',";l_{0}(bs) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerHits each final state particle",30,0,120) 
tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS = TProfile2D('tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS',";v_{z}(bs) interaction vertex #bar{S} (cm);l_{0}(bs) interaction vertex #bar{S} (cm);event fracion #geq 7 numberOfTrackerHits each final state particle",250,-125,125,1250,0,125) 
teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS= TProfile('teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS',";p_{t} #bar{S} (GeV);event fracion #geq 7 numberOfTrackerHits each final state particle",100,0,10) 
teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS= TProfile('teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS',";p_{z} #bar{S} (GeV);event fracion #geq 7 numberOfTrackerHits each final state particle",100,0,100) 
#check the correlation between antiS being reconstructed and all tracks having larger than a certain amount of hits
h2_allTracksMoreThan4Hits_efficiency = TH2F("h2_allTracksMoreThan4Hits_efficiency",";all tracks >= 5 hits?;#bar{S} reconstructed;",2,-0.5,1.5,2,-0.5,1.5)
h2_allTracksMoreThan5Hits_efficiency = TH2F("h2_allTracksMoreThan5Hits_efficiency",";all tracks >= 6 hits?;#bar{S} reconstructed;",2,-0.5,1.5,2,-0.5,1.5)
h2_allTracksMoreThan6Hits_efficiency = TH2F("h2_allTracksMoreThan6Hits_efficiency",";all tracks >= 7 hits?;#bar{S} reconstructed;",2,-0.5,1.5,2,-0.5,1.5)

#plot the variables used for matching the antiS and V0 GEN to RECO
h_GENRECO_matcher1_antiS = TH1F("h_GENRECO_matcher1_antiS",";#DeltaR #bar{S} (simulation, RECO);Events/0.1#DeltaR",20,0,2)
h_GENRECO_matcher1_Ks = TH1F("h_GENRECO_matcher1_Ks",";#DeltaL_{xyz,dv} K_{S}^{0} (simulation, RECO);Events/mm",100,0,10)
h_GENRECO_matcher1_AntiLambda = TH1F("h_GENRECO_matcher1_AntiLambda",";#DeltaL_{xyz,dv} #bar{#Lambda}^{0} (simulation, RECO);Events/mm",100,0,10)

h_GENRECO_matcher2_antiS = TH1F("h_GENRECO_matcher2_antiS",";#DeltaL_{xyz,iv}#bar{S} (simulation, RECO);Events/mm",100,0,10)
h_GENRECO_matcher2_Ks = TH1F("h_GENRECO_matcher2_Ks",";#Delta R K_{S}^{0} (simulation, RECO);Events/0.001#DeltaR",200,0,0.2)
h_GENRECO_matcher2_AntiLambda = TH1F("h_GENRECO_matcher2_AntiLambda",";#Delta R #bar{#Lambda}^{0} (simulation, RECO);Events/0.001#DeltaR",200,0,0.2)


#reconstruction efficiencies of the antiS, based on the lxyz between the RECO and GEN interaction vertex
h_AntiS_deltaLInteractionVertexAntiSmin = TH1F("h_AntiS_deltaLInteractionVertexAntiSmin",";min #DeltaL_{xyz}(RECO #bar{S},simulated #bar{S}) interaction vertex;Events/mm",1000,0,100) 


h_teff_nomAntiS_RECO_eff_eta_antiS= TH1F('h_teff_nomAntiS_RECO_eff_eta_antiS',";#eta simulated #bar{S};Efficiency",100,-5,5) 
h_teff_nomAntiS_RECO_eff_eta_antiS_antiS= TH1F('h_teff_nomAntiS_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};Efficiency",100,-5,5) 
h_teff_nomAntiS_RECO_eff_vz_antiS= TH1F('h_teff_nomAntiS_RECO_eff_vz_antiS',";absolute v_{z} interaction vertex simulated #bar{S} (cm);Efficiency",100,-200,200) 
h_teff_nomAntiS_RECO_eff_lxy_antiS= TH1F('h_teff_nomAntiS_RECO_eff_lxy_antiS',";l_{0}(bpc) interaction vertex simulated #bar{S} (cm);Efficiency",50,2,2.5) 
h_teff_nomAntiS_RECO_eff_pt_antiS= TH1F('h_teff_nomAntiS_RECO_eff_pt_antiS',";p_{t} simulated #bar{S} (GeV);Efficiency ",100,0,10) 
h_teff_nomAntiS_RECO_eff_pz_antiS= TH1F('h_teff_nomAntiS_RECO_eff_pz_antiS',";p_{z} simulated #bar{S} (GeV);Efficiency ",80,0,80) 
h_teff_nomAntiS_RECO_eff_p_antiS= TH1F('h_teff_nomAntiS_RECO_eff_p_antiS',";p simulated #bar{S} (GeV);Efficiency ",400,0,40) 
h_teff_nomAntiS_RECO_eff_dxy_antiS= TH1F('h_teff_nomAntiS_RECO_eff_dxy_antiS',";d_{0}(bs) simulated #bar{S} (cm);Efficiency ",40,-20,20) 
h_teff_nomAntiS_RECO_eff_dz_antiS= TH1F('h_teff_nomAntiS_RECO_eff_dz_antiS',";d_{z}(bs) simulated #bar{S} (cm);Efficiency ",40,-100,100) 
h_teff_nomAntiS_RECO_eff_numberOfTrackerHits_antiS= TH1F('h_teff_nomAntiS_RECO_eff_numberOfTrackerHits_antiS',";numberOfTrackerHits simulated #bar{S} ;Efficiency ",40,0-0.5,40-0.5) 
h_teff_nomAntiS_RECO_eff_phi_antiS= TH1F('h_teff_nomAntiS_RECO_eff_phi_antiS',";#phi simulated #bar{S};Efficiency",100,-5,5) 
#reco eff for the Ks
h_teff_nomKs_RECO_eff_eta_antiS= TH1F('h_teff_nomKs_RECO_eff_eta_antiS',";#eta simulated K_{S}^{0};",100,-5,5) 
h_teff_nomKs_RECO_eff_eta_antiS_antiS= TH1F('h_teff_nomKs_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};",100,-5,5) 
h_teff_nomKs_RECO_eff_vz_antiS= TH1F('h_teff_nomKs_RECO_eff_vz_antiS',";absolute v_{z} decay vertex simulated K_{S}^{0} (cm);",100,-200,200)
h_teff_nomKs_RECO_eff_lxy_antiS= TH1F('h_teff_nomKs_RECO_eff_lxy_antiS',";l_{0}(bs) decay vertex simulated K_{S}^{0} (cm);",80,0,80) 
h_teff_nomKs_RECO_eff_pt_antiS= TH1F('h_teff_nomKs_RECO_eff_pt_antiS',";p_{t} simulated K_{S}^{0} (GeV); ",100,0,10) 
h_teff_nomKs_RECO_eff_pz_antiS= TH1F('h_teff_nomKs_RECO_eff_pz_antiS',";p_{z} simulated K_{S}^{0} (GeV); ",80,0,80) 
h_teff_nomKs_RECO_eff_p_antiS= TH1F('h_teff_nomKs_RECO_eff_p_antiS',";p simulated K_{S}^{0} (GeV); ",400,0,40) 
h_teff_nomKs_RECO_eff_dxy_antiS= TH1F('h_teff_nomKs_RECO_eff_dxy_antiS',";d_{0}(bs) simulated K_{S}^{0} (cm); ",40,-20,20) 
h_teff_nomKs_RECO_eff_dz_antiS= TH1F('h_teff_nomKs_RECO_eff_dz_antiS',";d_{z}(bs) simulated K_{S}^{0} (cm); ",40,-100,100) 
h_teff_nomKs_RECO_eff_numberOfTrackerHits_antiS= TH1F('h_teff_nomKs_RECO_eff_numberOfTrackerHits_antiS',";numberOfTrackerHits simulated K_{S}^{0} ; ",40,0-0.5,40-0.5) 
h_teff_nomKs_RECO_eff_phi_antiS= TH1F('h_teff_nomKs_RECO_eff_phi_antiS',";#phi simulated K_{S}^{0};",100,-5,5) 
#reco eff for the AntiLambda
h_teff_nomAntiLambda_RECO_eff_eta_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_eta_antiS',";#eta simulated #bar{#Lambda}^{0};Efficiency",100,-5,5) 
h_teff_nomAntiLambda_RECO_eff_eta_antiS_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};Efficiency",100,-5,5) 
h_teff_nomAntiLambda_RECO_eff_vz_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_vz_antiS',";absolute v_{z} decay vertex simulated #bar{#Lambda}^{0} (cm);Efficiency",100,-200,200) 
h_teff_nomAntiLambda_RECO_eff_lxy_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_lxy_antiS',";l_{0}(bs) decay vertex simulated #bar{#Lambda}^{0} (cm);Efficiency",80,0,80) 
h_teff_nomAntiLambda_RECO_eff_pt_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_pt_antiS',";p_{t} simulated #bar{#Lambda}^{0} (GeV);Efficiency ",100,0,10) 
h_teff_nomAntiLambda_RECO_eff_pz_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_pz_antiS',";p_{z} simulated #bar{#Lambda}^{0} (GeV);Efficiency ",80,0,80) 
h_teff_nomAntiLambda_RECO_eff_p_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_p_antiS',";p simulated #bar{#Lambda}^{0} (GeV);Efficiency ",400,0,40) 
h_teff_nomAntiLambda_RECO_eff_dxy_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_dxy_antiS',";d_{0}(bs) simulated #bar{#Lambda}^{0} (cm);Efficiency ",40,-20,20) 
h_teff_nomAntiLambda_RECO_eff_dz_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_dz_antiS',";d_{z}(bs) simulated #bar{#Lambda}^{0} (cm);Efficiency ",40,-100,100) 
h_teff_nomAntiLambda_RECO_eff_numberOfTrackerHits_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_numberOfTrackerHits_antiS',";numberOfTrackerHits simulated #bar{#Lambda}^{0} (cm);Efficiency ",40,0-0.5,40-0.5) 
h_teff_nomAntiLambda_RECO_eff_phi_antiS= TH1F('h_teff_nomAntiLambda_RECO_eff_phi_antiS',";#phi simulated #bar{#Lambda}^{0};Efficiency",100,-5,5) 
#reco eff for the Ks daughter 0 and 1
h_teff_nomKsdaugthers_RECO_eff_eta_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_eta_antiS',";#eta simulated K_{S}^{0} daughters;Efficiency",100,-5,5) 
h_teff_nomKsdaughters_RECO_eff_eta_antiS_antiS= TH1F('h_teff_nomKsdaughters_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};Efficiency",100,-5,5) 
h_teff_nomKsdaugthers_RECO_eff_vz_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_vz_antiS',";absolute v_{z} cv simulated K_{S}^{0} daughters (cm);Efficiency",100,-200,200)
h_teff_nomKsdaugthers_RECO_eff_lxy_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_lxy_antiS',";l_{0}(bs) cv simulated K_{S}^{0} daughters (cm);Efficiency",80,0,80) 
h_teff_nomKsdaugthers_RECO_eff_pt_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_pt_antiS',";p_{t} simulated K_{S}^{0} daughters (GeV);Efficiency ",100,0,10) 
h_teff_nomKsdaugthers_RECO_eff_pz_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_pz_antiS',";p_{z} simulated K_{S}^{0} daughters (GeV);Efficiency ",80,0,80) 
h_teff_nomKsdaugthers_RECO_eff_p_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_p_antiS',";p simulated K_{S}^{0} daughters (GeV);Efficiency ",400,0,40) 
h_teff_nomKsdaugthers_RECO_eff_dxy_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_dxy_antiS',";d_{0}(bs) simulated K_{S}^{0} daughters (cm);Efficiency ",40,-20,20) 
h_teff_nomKsdaugthers_RECO_eff_dz_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_dz_antiS',";d_{z}(bs) simulated K_{S}^{0} daughters (cm);Efficiency ",40,-100,100) 
h_teff_nomKsdaugthers_RECO_eff_numberOfTrackerHits_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_numberOfTrackerHits_antiS',";numberOfTrackerHits simulated K_{S}^{0} daughters;Efficiency ",40,0-0.5,40-0.5) 
h_teff_nomKsdaugthers_RECO_eff_phi_antiS= TH1F('h_teff_nomKsdaugthers_RECO_eff_phi_antiS',";#phi simulated K_{S}^{0} daughters;Efficiency",100,-5,5) 
#reco eff for the soft pion from the AntiLambda
h_teff_nomAntiLambdaPion_RECO_eff_eta_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_eta_antiS',";#eta simulated #bar{#Lambda}^{0}-#pi^{+};Efficiency",100,-5,5) 
h_teff_nomAntiLambdaPion_RECO_eff_eta_antiS_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};Efficiency",100,-5,5) 
h_teff_nomAntiLambdaPion_RECO_eff_vz_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_vz_antiS',";absolute v_{z} cv simulated #bar{#Lambda}^{0}-#pi^{+} (cm);Efficiency",100,-200,200) 
h_teff_nomAntiLambdaPion_RECO_eff_lxy_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_lxy_antiS',";l_{0}(bs) cv simulated #bar{#Lambda}^{0}-#pi^{+} (cm);Efficiency",80,0,80) 
h_teff_nomAntiLambdaPion_RECO_eff_pt_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_pt_antiS',";p_{t} simulated #bar{#Lambda}^{0}-#pi^{+} (GeV);Efficiency ",100,0,10) 
h_teff_nomAntiLambdaPion_RECO_eff_pz_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_pz_antiS',";p_{z} simulated #bar{#Lambda}^{0}-#pi^{+} (GeV);Efficiency ",80,0,80) 
h_teff_nomAntiLambdaPion_RECO_eff_p_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_p_antiS',";p simulated #bar{#Lambda}^{0}-#pi^{+} (GeV);Efficiency ",400,0,40) 
h_teff_nomAntiLambdaPion_RECO_eff_dxy_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_dxy_antiS',";d_{0}(bs) simulated #bar{#Lambda}^{0}-#pi^{+} (cm);Efficiency ",40,-20,20) 
h_teff_nomAntiLambdaPion_RECO_eff_dz_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_dz_antiS',";d_{z}(bs) simulated #bar{#Lambda}^{0}-#pi^{+} (cm);Efficiency ",40,-100,100) 
h_teff_nomAntiLambdaPion_RECO_eff_numberOfTrackerHits_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_numberOfTrackerHits_antiS',";numberOfTrackerHits simulated #bar{#Lambda}^{0}-#pi^{+};Efficiency ",40,0-0.5,40-0.5) 
h_teff_nomAntiLambdaPion_RECO_eff_phi_antiS= TH1F('h_teff_nomAntiLambdaPion_RECO_eff_phi_antiS',";#phi simulated #bar{#Lambda}^{0}-#pi^{+};Efficiency",100,-5,5) 
#reco eff for the antiproton from the AntiLambda
h_teff_nomAntiLambdaAntiProton_RECO_eff_eta_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_eta_antiS',";#eta simulated #bar{#Lambda}^{0}-#bar{p};Efficiency",100,-5,5) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_eta_antiS_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_eta_antiS_antiS',";#eta simulated #bar{S};Efficiency",100,-5,5) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_vz_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_vz_antiS',";absolute v_{z} cv simulated #bar{#Lambda}^{0}-#bar{p} (cm);Efficiency",100,-200,200) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_lxy_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_lxy_antiS',";l_{0}(bs) cv simulated #bar{#Lambda}^{0}-#bar{p} (cm);Efficiency",80,0,80) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_pt_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_pt_antiS',";p_{t} simulated #bar{#Lambda}^{0}-#bar{p} (GeV);Efficiency ",100,0,10) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_pz_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_pz_antiS',";p_{z} simulated #bar{#Lambda}^{0}-#bar{p} (GeV);Efficiency ",80,0,80) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_p_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_p_antiS',";p simulated #bar{#Lambda}^{0}-#bar{p} (GeV);Efficiency ",400,0,40) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_dxy_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_dxy_antiS',";d_{0}(bs) simulated #bar{#Lambda}^{0}-#bar{p} (cm);Efficiency ",40,-20,20) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_dz_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_dz_antiS',";d_{z}(bs) simulated #bar{#Lambda}^{0}-#bar{p} (cm);Efficiency ",40,-100,100) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_numberOfTrackerHits_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_numberOfTrackerHits_antiS',";numberOfTrackerHits simulated #bar{#Lambda}^{0}-#bar{p};Efficiency ",40,0-0.5,40-0.5) 
h_teff_nomAntiLambdaAntiProton_RECO_eff_phi_antiS= TH1F('h_teff_nomAntiLambdaAntiProton_RECO_eff_phi_antiS',";#phi simulated #bar{#Lambda}^{0}-#bar{p};Efficiency",100,-5,5) 

ll_efficiencies_nom = [
[h_teff_nomAntiS_RECO_eff_eta_antiS,h_teff_nomAntiS_RECO_eff_eta_antiS_antiS,h_teff_nomAntiS_RECO_eff_vz_antiS,h_teff_nomAntiS_RECO_eff_lxy_antiS,h_teff_nomAntiS_RECO_eff_pt_antiS,h_teff_nomAntiS_RECO_eff_pz_antiS,h_teff_nomAntiS_RECO_eff_p_antiS,h_teff_nomAntiS_RECO_eff_dxy_antiS,h_teff_nomAntiS_RECO_eff_dz_antiS,h_teff_nomAntiS_RECO_eff_numberOfTrackerHits_antiS,h_teff_nomAntiS_RECO_eff_phi_antiS],
[h_teff_nomKs_RECO_eff_eta_antiS,h_teff_nomKs_RECO_eff_eta_antiS_antiS,h_teff_nomKs_RECO_eff_vz_antiS,h_teff_nomKs_RECO_eff_lxy_antiS,h_teff_nomKs_RECO_eff_pt_antiS,h_teff_nomKs_RECO_eff_pz_antiS,h_teff_nomKs_RECO_eff_p_antiS,h_teff_nomKs_RECO_eff_dxy_antiS,h_teff_nomKs_RECO_eff_dz_antiS,h_teff_nomKs_RECO_eff_numberOfTrackerHits_antiS,h_teff_nomKs_RECO_eff_phi_antiS],
[h_teff_nomAntiLambda_RECO_eff_eta_antiS,h_teff_nomAntiLambda_RECO_eff_eta_antiS_antiS,h_teff_nomAntiLambda_RECO_eff_vz_antiS,h_teff_nomAntiLambda_RECO_eff_lxy_antiS,h_teff_nomAntiLambda_RECO_eff_pt_antiS,h_teff_nomAntiLambda_RECO_eff_pz_antiS,h_teff_nomAntiLambda_RECO_eff_p_antiS,h_teff_nomAntiLambda_RECO_eff_dxy_antiS,h_teff_nomAntiLambda_RECO_eff_dz_antiS,h_teff_nomAntiLambda_RECO_eff_numberOfTrackerHits_antiS,h_teff_nomAntiLambda_RECO_eff_phi_antiS],
[h_teff_nomKsdaugthers_RECO_eff_eta_antiS,h_teff_nomKsdaughters_RECO_eff_eta_antiS_antiS,h_teff_nomKsdaugthers_RECO_eff_vz_antiS,h_teff_nomKsdaugthers_RECO_eff_lxy_antiS,h_teff_nomKsdaugthers_RECO_eff_pt_antiS,h_teff_nomKsdaugthers_RECO_eff_pz_antiS,h_teff_nomKsdaugthers_RECO_eff_p_antiS,h_teff_nomKsdaugthers_RECO_eff_dxy_antiS,h_teff_nomKsdaugthers_RECO_eff_dz_antiS,h_teff_nomKsdaugthers_RECO_eff_numberOfTrackerHits_antiS,h_teff_nomKsdaugthers_RECO_eff_phi_antiS],
[h_teff_nomAntiLambdaPion_RECO_eff_eta_antiS,h_teff_nomAntiLambdaPion_RECO_eff_eta_antiS_antiS,h_teff_nomAntiLambdaPion_RECO_eff_vz_antiS,h_teff_nomAntiLambdaPion_RECO_eff_lxy_antiS,h_teff_nomAntiLambdaPion_RECO_eff_pt_antiS,h_teff_nomAntiLambdaPion_RECO_eff_pz_antiS,h_teff_nomAntiLambdaPion_RECO_eff_p_antiS,h_teff_nomAntiLambdaPion_RECO_eff_dxy_antiS,h_teff_nomAntiLambdaPion_RECO_eff_dz_antiS,h_teff_nomAntiLambdaPion_RECO_eff_numberOfTrackerHits_antiS,h_teff_nomAntiLambdaPion_RECO_eff_phi_antiS],
[h_teff_nomAntiLambdaAntiProton_RECO_eff_eta_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_eta_antiS_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_vz_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_lxy_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_pt_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_pz_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_p_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_dxy_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_dz_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_numberOfTrackerHits_antiS,h_teff_nomAntiLambdaAntiProton_RECO_eff_phi_antiS]
]

for l in ll_efficiencies_nom:
	for h in l:
		h.Clone()
		h.SetDirectory(0)

#create the same list of histograms as in ll_efficiencies_nom, but now for the denominator, have to do like this, because I have to use weights and teff does not allow weights
ll_efficiencies_denom = []
for l in range(0,len(ll_efficiencies_nom)):
	ll_efficiencies_denom.append([])
	for h in range(0,len(ll_efficiencies_nom[0])):
		histo = ll_efficiencies_nom[l][h]
		histoClone = histo.Clone()
		histoClone.SetName((histo.GetName()).replace('nom','denom'))
		histoClone.SetDirectory(0)
		ll_efficiencies_denom[l].append(histoClone)

for l in ll_efficiencies_denom:
	for h in l:
		h.Clone()
		h.SetDirectory(0)


#create the same list of histograms as in ll_efficiencies but fill these only for antiS which have all granddaughters within tracker acceptance, where this is based on having >=7 hits on all tracks.
ll_efficiencies_acceptance_nom = []
for l in range(0,len(ll_efficiencies_nom)):
	ll_efficiencies_acceptance_nom.append([])
	for h in range(0,len(ll_efficiencies_nom[0])):
		histo = ll_efficiencies_nom[l][h]
		histoClone = histo.Clone()
		histoClone.SetName(histo.GetName()+"_acceptance")
		histoClone.SetDirectory(0)
		ll_efficiencies_acceptance_nom[l].append(histoClone)

for l in ll_efficiencies_acceptance_nom:
	for h in l:
		h.Clone()
		h.SetDirectory(0)

ll_efficiencies_acceptance_denom = []
for l in range(0,len(ll_efficiencies_nom)):
	ll_efficiencies_acceptance_denom.append([])
	for h in range(0,len(ll_efficiencies_nom[0])):
		histo = ll_efficiencies_nom[l][h]
		histoClone = histo.Clone()
		histoClone.SetName((histo.GetName()+"_acceptance").replace('nom','denom'))
		histoClone.SetDirectory(0)
		ll_efficiencies_acceptance_denom[l].append(histoClone)

for l in ll_efficiencies_acceptance_denom:
	for h in l:
		h.Clone()
		h.SetDirectory(0)

#reconstruction accuracy of the antiS
h1_AntiS_RECO_Acc_eta = TH1F("h1_AntiS_RECO_Acc_eta","; #eta_{sim #bar{S}} - #eta_{reco #bar{S}};Events/0.01#eta",80,-0.4,0.4)
h1_AntiS_RECO_Acc_phi = TH1F("h1_AntiS_RECO_Acc_phi","; #phi_{sim #bar{S}} - #phi_{reco #bar{S}} (rad);Events/0.01rad",80,-0.4,0.4)
h1_AntiS_RECO_Acc_vz = TH1F("h1_AntiS_RECO_Acc_vz",";  v_{z, sim #bar{S}}(iv, bs) - v_{z, reco #bar{S}}(iv, bs) (cm);Events/0.1mm",100,-1,1)
h1_AntiS_RECO_Acc_lxy = TH1F("h1_AntiS_RECO_Acc_lxy","; l_{0, sim #bar{S}}(iv, bs) - l_{0, reco #bar{S}}(iv, bs) (cm);Events/0.1mm",100,-1,1)
h1_AntiS_RECO_Acc_pt = TH1F("h1_AntiS_RECO_Acc_pt","; p_{t, sim #bar{S}} - p_{t, reco #bar{S}} (GeV);Events/0.02GeV",100,-1,1)
h1_AntiS_RECO_Acc_pz = TH1F("h1_AntiS_RECO_Acc_pz","; p_{z, sim #bar{S}} - p_{z, reco #bar{S}} (GeV);Events/0.02GeV",150,-1.5,1.5)

#invariant mass of the antiS
h_AntiS_massMinusNeutron = TH1F("h_AntiS_massMinusNeutron",";m_{#bar{S},obs} RECO (GeV);Events/0.1GeV",180,-6,12) 
h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda = TH2F("h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda","; m_{#bar{S},obs} (GeV); |#vec{p}_{K_{s}^{0}, RECO} + #vec{p}_{#bar{#Lambda}^{0}, RECO}| (GeV);Events/GeV^{2}",110,-5,6,60,0,40)

#list of list for histograms containging kinematics of the granddaughters of the AntiS, for AntiS which got reconstructed. These plots are necessary because they tell where you need to be sure that your tracking is properly described by MC. 
h_AntiSGrandDaughter_lxy_1 = TH1F("h_AntiSGrandDaughter_lxy_1",";Simulated track l_{0}(cv, bs) (cm);Events/cm",60,0,60) 
h_AntiSGrandDaughter_vz_1 = TH1F("h_AntiSGrandDaughter_vz_1",";Simulated track v_{z}(cv, bs) (cm);Events/10cm",40,-200,200) 
h_AntiSGrandDaughter_dxy_1 = TH1F("h_AntiSGrandDaughter_dxy_1",";Simulated track d_{0}(bs) (cm);Events/cm",40,-20,20) 
h_AntiSGrandDaughter_dz_1 = TH1F("h_AntiSGrandDaughter_dz_1",";Simulated track d_{z}(bs) (cm);Events/4cm",30,-60,60) 
h_AntiSGrandDaughter_pt_1 = TH1F("h_AntiSGrandDaughter_pt_1",";Simulated track p_{t} (GeV);Events/0.1GeV",50,0,5) 
h_AntiSGrandDaughter_pz_1 = TH1F("h_AntiSGrandDaughter_pz_1",";Simulated track p_{z} (GeV);Events/GeV",20,0,20) 

h_AntiSGrandDaughter_lxy_3 = TH1F("h_AntiSGrandDaughter_lxy_3",";Simulated track l_{0}(cv, bs) (cm);Events/cm",60,0,60)
h_AntiSGrandDaughter_vz_3 = TH1F("h_AntiSGrandDaughter_vz_3",";Simulated track v_{z}(cv, bs) (cm);Events/10cm",40,-200,200)  
h_AntiSGrandDaughter_dxy_3 = TH1F("h_AntiSGrandDaughter_dxy_3",";Simulated track d_{0}(bs) (cm);Events/cm",40,-20,20) 
h_AntiSGrandDaughter_dz_3 = TH1F("h_AntiSGrandDaughter_dz_3",";Simulated track d_{z}(bs) (cm);Events/4cm",30,-60,60)    
h_AntiSGrandDaughter_pt_3 = TH1F("h_AntiSGrandDaughter_pt_3",";Simulated track p_{t} (GeV);Events/0.1GeV",50,0,5)
h_AntiSGrandDaughter_pz_3 = TH1F("h_AntiSGrandDaughter_pz_3",";Simulated track p_{z} (GeV);Events/GeV",20,0,20)   

h_AntiSGrandDaughter_lxy_4 = TH1F("h_AntiSGrandDaughter_lxy_4",";Simulated track l_{0}(cv, bs) (cm);Events/cm",60,0,60)
h_AntiSGrandDaughter_vz_4 = TH1F("h_AntiSGrandDaughter_vz_4",";Simulated track v_{z}(cv, bs) (cm);Events/10cm",40,-200,200) 
h_AntiSGrandDaughter_dxy_4 = TH1F("h_AntiSGrandDaughter_dxy_4",";Simulated track d_{0}(bs) (cm);Events/cm",40,-20,20)
h_AntiSGrandDaughter_dz_4 = TH1F("h_AntiSGrandDaughter_dz_4",";Simulated track d_{z}(bs) (cm);Events/4cm",30,-60,60)   
h_AntiSGrandDaughter_pt_4 = TH1F("h_AntiSGrandDaughter_pt_4",";Simulated track p_{t} (GeV);Events/0.1GeV",50,0,5)
h_AntiSGrandDaughter_pz_4 = TH1F("h_AntiSGrandDaughter_pz_4",";Simulated track p_{z} (GeV);Events/GeV",20,0,20) 

ll_kinematics_granddaughters_of_RECO_AntiS = [
[h_AntiSGrandDaughter_lxy_1,h_AntiSGrandDaughter_vz_1,h_AntiSGrandDaughter_dxy_1,h_AntiSGrandDaughter_dz_1,h_AntiSGrandDaughter_pt_1,h_AntiSGrandDaughter_pz_1],
[h_AntiSGrandDaughter_lxy_3,h_AntiSGrandDaughter_vz_3,h_AntiSGrandDaughter_dxy_3,h_AntiSGrandDaughter_dz_3,h_AntiSGrandDaughter_pt_3,h_AntiSGrandDaughter_pz_3],
[h_AntiSGrandDaughter_lxy_4,h_AntiSGrandDaughter_vz_4,h_AntiSGrandDaughter_dxy_4,h_AntiSGrandDaughter_dz_4,h_AntiSGrandDaughter_pt_4,h_AntiSGrandDaughter_pz_4],
]

h_V0FitterCuts_Ks = TH1I("h_V0FitterCuts_Ks","; does V0 for which both tracks surv V0Fitter track cuts get RECO (0 = reconstructed);",62,-1.5,60.5)
h_V0FitterCuts_AntiLambda = TH1I("h_V0FitterCuts_AntiLambda","; does V0 for which both tracks surv V0Fitter track cuts get RECO (0 = reconstructed);",62,-1.5,60.5)
h_tracks_V0FitterCuts_posPion_Ks = TH1I("h_tracks_V0FitterCuts_posPion_Ks",";does RECO track survive V0Fitter track cuts (0 = reconstructed);",62,-1.5,60.5)
h_tracks_V0FitterCuts_negPion_Ks = TH1I("h_tracks_V0FitterCuts_negPion_Ks",";does RECO track survive V0Fitter track cuts (0 = reconstructed);",62,-1.5,60.5)
h_tracks_V0FitterCuts_posPion_AntiLambda = TH1I("h_tracks_V0FitterCuts_posPion_AntiLambda",";does RECO track survive V0Fitter track cuts (0 = reconstructed);",62,-1.5,60.5)
h_tracks_V0FitterCuts_antiProton_AntiLambda = TH1I("h_tracks_V0FitterCuts_antiProton_AntiLambda",";does RECO track survive V0Fitter track cuts (0 = reconstructed);",62,-1.5,60.5)
l_h_V0FitterCuts = [h_V0FitterCuts_Ks,h_V0FitterCuts_AntiLambda,h_tracks_V0FitterCuts_posPion_Ks,h_tracks_V0FitterCuts_negPion_Ks,h_tracks_V0FitterCuts_posPion_AntiLambda,h_tracks_V0FitterCuts_antiProton_AntiLambda] 


#a list with counters for the reconstructed particles, so there are 7 entries for each of the 7 particles
nAntiS = 0.
nAntiS_reconstructable = 0. 
nTotalReconstructed_reconstructable = [0.,0.,0.,0.,0.,0.,0.]
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

		if(i>maxEvents2):
			break

		tree.GetEntry(i)

		weightFactor = tree._tpsAntiS_event_weighting_factor[0]*tree._tpsAntiS_event_weighting_factorPU[0]

		nAntiS+=weightFactor

		#my requirement for having reconstructed antiS: based on the 3D distance between the GEN and RECO interaction vertex of the antiS.
		antiSReconstructed = False
		KsReconstructed = False
		antiLambdaReconstructed = False
		if(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[1] < config_dict["GENRECO_matcher_Ks_deltaL"] and tree._tpsAntiS_bestDeltaRWithRECO[1] < config_dict["GENRECO_matcher_Ks_deltaR"]):
			KsReconstructed = True
		if(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[2] < config_dict["GENRECO_matcher_AntiL_deltaL"] and tree._tpsAntiS_bestDeltaRWithRECO[2] < config_dict["GENRECO_matcher_AntiL_deltaR"]):
			antiLambdaReconstructed = True
		if(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[0] < config_dict["GENRECO_matcher_AntiS_deltaL"] and tree._tpsAntiS_bestDeltaRWithRECO[0] < config_dict["GENRECO_matcher_AntiS_deltaR"] and tree._tpsAntiS_reconstructed[3] == 1 and tree._tpsAntiS_reconstructed[4] == 1 and tree._tpsAntiS_reconstructed[5] == 1 and tree._tpsAntiS_reconstructed[6] == 1 and KsReconstructed and antiLambdaReconstructed):
			antiSReconstructed = True	
	
		printProgress(i)

                #just a check of the tree: if a antiS got reconstructed then also all the daughters should be reconstructed.  
                if(antiSReconstructed): 
                        print antiSReconstructed," ", KsReconstructed, " ", antiLambdaReconstructed, " ",tree._tpsAntiS_reconstructed[3], " ",tree._tpsAntiS_reconstructed[4], " ",tree._tpsAntiS_reconstructed[5], " ",tree._tpsAntiS_reconstructed[6]


		#count the number of granddaughters of this antiS which have >= a certain amount of hits
		NGrandDaughtersWithTrackerHitsLargerThan4 = 0
		NGrandDaughtersWithTrackerHitsLargerThan5 = 0
		NGrandDaughtersWithTrackerHitsLargerThan6 = 0
		for j in range(0,len(tree._tpsAntiS_type)):
			#now look at _tpAntiS_type from 3 to 6, this are the granddaughters:
			if(tree._tpsAntiS_type[j] >= 3 and tree._tpsAntiS_type[j] <= 6 ):
				#print "particle type ", tree._tpsAntiS_type[j], " has numberOfTrackerHits = " , tree._tpsAntiS_numberOfTrackerHits[j]
				if(tree._tpsAntiS_numberOfTrackerHits[j] >= 5):
					NGrandDaughtersWithTrackerHitsLargerThan4 += 1
				if(tree._tpsAntiS_numberOfTrackerHits[j] >= 6):
					NGrandDaughtersWithTrackerHitsLargerThan5 += 1
				if(tree._tpsAntiS_numberOfTrackerHits[j] >= 7):
					NGrandDaughtersWithTrackerHitsLargerThan6 += 1

		tprof_numberGranddaughters_7hits_eta_antiS.Fill(tree._tpsAntiS_eta[0],NGrandDaughtersWithTrackerHitsLargerThan6,weightFactor)
		tprof_numberGranddaughters_7hits_vz_interaction_antiS.Fill(tree._tpsAntiS_vz_beamspot[1],NGrandDaughtersWithTrackerHitsLargerThan6,weightFactor) #interaction vtx is the creation vtx of the daughter
		tprof_numberGranddaughters_7hits_lxy_interaction_antiS.Fill(tree._tpsAntiS_Lxy_beamspot[1],NGrandDaughtersWithTrackerHitsLargerThan6,weightFactor) #interaction vtx is the creation vtx of the daughter

		boolNGrandDaughtersWithTrackerHitsLargerThan4 =  NGrandDaughtersWithTrackerHitsLargerThan4==4
		boolNGrandDaughtersWithTrackerHitsLargerThan5 =  NGrandDaughtersWithTrackerHitsLargerThan5==4
		boolNGrandDaughtersWithTrackerHitsLargerThan6 =  NGrandDaughtersWithTrackerHitsLargerThan6==4

		#plots having the fraction of all events with all granddaughters with >= 7 hits in function of a kinematic variable of the antiS
		teff_fractionAllEvents_numberGranddaughters_7hits_eta_antiS.Fill(tree._tpsAntiS_eta[0],boolNGrandDaughtersWithTrackerHitsLargerThan6,weightFactor)
		teff_fractionAllEvents_numberGranddaughters_7hits_phi_antiS.Fill(tree._tpsAntiS_phi[0],boolNGrandDaughtersWithTrackerHitsLargerThan6,weightFactor)
		teff_fractionAllEvents_numberGranddaughters_7hits_vz_antiS.Fill(tree._tpsAntiS_vz_beamspot[1],boolNGrandDaughtersWithTrackerHitsLargerThan6,weightFactor)
		teff_fractionAllEvents_numberGranddaughters_7hits_lxy_antiS.Fill(tree._tpsAntiS_Lxy_beamspot[1],boolNGrandDaughtersWithTrackerHitsLargerThan6,weightFactor)
		tprof2_fractionAllEvents_numberGranddaughters_7hits_vz_lxy_antiS.Fill(tree._tpsAntiS_vz_beamspot[1],tree._tpsAntiS_Lxy_beamspot[1],boolNGrandDaughtersWithTrackerHitsLargerThan6,weightFactor)
		teff_fractionAllEvents_numberGranddaughters_7hits_pt_antiS.Fill(tree._tpsAntiS_pt[0],boolNGrandDaughtersWithTrackerHitsLargerThan6,weightFactor)
		teff_fractionAllEvents_numberGranddaughters_7hits_pz_antiS.Fill(tree._tpsAntiS_pz[0],boolNGrandDaughtersWithTrackerHitsLargerThan6,weightFactor)
		#counter for all the antiS which have all granddaughters producing >= 7 hits
		if(boolNGrandDaughtersWithTrackerHitsLargerThan6):
			nAntiSWithAllGranddaughtersMoreThan6Hits+=weightFactor

		#check now how well boolNGrandDaughtersWithTrackerHitsLargerThan6 is a proxy for the efficiency
		h2_allTracksMoreThan4Hits_efficiency.Fill(boolNGrandDaughtersWithTrackerHitsLargerThan4,antiSReconstructed,weightFactor)
		h2_allTracksMoreThan5Hits_efficiency.Fill(boolNGrandDaughtersWithTrackerHitsLargerThan5,antiSReconstructed,weightFactor)
		h2_allTracksMoreThan6Hits_efficiency.Fill(boolNGrandDaughtersWithTrackerHitsLargerThan6,antiSReconstructed,weightFactor)
		h_AntiS_deltaLInteractionVertexAntiSmin.Fill(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[0])

		
		if(tree._tpsAntiS_reconstructed[3] == 1 and tree._tpsAntiS_reconstructed[4] == 1 and tree._tpsAntiS_reconstructed[5] == 1 and tree._tpsAntiS_reconstructed[6] == 1):
			h_GENRECO_matcher1_antiS.Fill(tree._tpsAntiS_bestDeltaRWithRECO[0],weightFactor)
		h_GENRECO_matcher1_Ks.Fill(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[1],weightFactor)
		h_GENRECO_matcher1_AntiLambda.Fill(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[2],weightFactor)

		#then apply the rough cuts and check the cut on the 2nd parameter
		if(tree._tpsAntiS_bestDeltaRWithRECO[0] < config_dict["GENRECO_matcher_AntiS_deltaR"] and tree._tpsAntiS_reconstructed[3] == 1 and tree._tpsAntiS_reconstructed[4] == 1 and tree._tpsAntiS_reconstructed[5] == 1 and tree._tpsAntiS_reconstructed[6] == 1):
			h_GENRECO_matcher2_antiS.Fill(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[0],weightFactor)
		if(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[1] < config_dict["GENRECO_matcher_Ks_deltaL"]):
			h_GENRECO_matcher2_Ks.Fill(tree._tpsAntiS_bestDeltaRWithRECO[1],weightFactor)
		if(tree._tpsAntiS_deltaLInteractionVertexAntiSmin[2] < config_dict["GENRECO_matcher_AntiL_deltaL"]):
			h_GENRECO_matcher2_AntiLambda.Fill(tree._tpsAntiS_bestDeltaRWithRECO[2],weightFactor)


		#investigate the cuts in the V0Fitter:
		for j in [1,2,3,4,5,6]:#for the V0s itself and for the tracks investigate the _tpsAntiS_returnCodeV0Fitter leaf which tells something about where the track or the V0 got cut (or survived the) V0Fitter
			#the code stored in _tpsAntiS_bestRECO_returnCodeV0Fitter for the tracks is a bit special. It is actually a binary number (with a length of 5)
			l_h_V0FitterCuts[j-1].Fill(tree._tpsAntiS_returnCodeV0Fitter[j])
		#quick check, that if the V0 it flagged as reconstructed that then also both tracks are reconstructed:
#		if(tree._tpsAntiS_reconstructed[1]):
#			print tree._tpsAntiS_reconstructed[3], "; ", tree._tpsAntiS_reconstructed[4]




		#only do the below for reconstructable events:
		if(not boolNGrandDaughtersWithTrackerHitsLargerThan6): continue

		nAntiS_reconstructable+=weightFactor

		#count the reconstructed particles with the requirement that there daughters got reconstructed
		if(tree._tpsAntiS_reconstructed[3] == 1 and tree._tpsAntiS_reconstructed[4] == 1):
			if(KsReconstructed == 1):
				nKsRECOIfBothDaughtersReco += weightFactor
			nKsTOTALIfBothDaughtersReco +=weightFactor
		if(tree._tpsAntiS_reconstructed[5] == 1 and tree._tpsAntiS_reconstructed[6] == 1):
			if(antiLambdaReconstructed == 1):
				nAntiLambdaRECOIfBothDaughtersReco += weightFactor
			nAntiLambdaTOTALIfBothDaughtersReco += weightFactor
		if(tree._tpsAntiS_reconstructed[1] == 1 and tree._tpsAntiS_reconstructed[2] == 1):
			if(antiSReconstructed):
				nAntiSRECOIfBothDaughtersReco += weightFactor
			nAntiSTOTALIfBothDaughtersReco += weightFactor

		#now instead of looking at the NGrandDaughtersWithTrackerHitsLargerThan6 as a proxy for tracking efficiency now look at the real tracking efficiency for antiS and its daughters
		for i_particle in range(0,7):#for all of the 7 particles
			index = i_particle
			if (i_particle == 3 or i_particle == 4):#take the two daughters of the Ks together
				index = 3
			elif(i_particle == 5):#because I took the previous one together
				index = 4
			elif(i_particle == 6):
				index = 5

			particleReconstructed = int(tree._tpsAntiS_reconstructed[i_particle])
			if(index == 0): #for the antiS use the definition on top to say if it was reconstructed, not the one from the tree
				particleReconstructed = antiSReconstructed
			if(index == 1): 
				particleReconstructed = KsReconstructed 
			if(index == 2): 
				particleReconstructed = antiLambdaReconstructed 

			#make a global count of what got reconstructed
			if(particleReconstructed): 
				nTotalReconstructed_reconstructable[i_particle] += weightFactor

			#fill reconstruction efficiency plots 	
			if(particleReconstructed):
				FillHistosEfficiency(tree,ll_efficiencies_nom,i_particle,index,weightFactor,i)
			FillHistosEfficiency(tree,ll_efficiencies_denom,i_particle,index,weightFactor,i)
			#fill reconstruction efficiency plots for events within acceptance
			if(particleReconstructed and boolNGrandDaughtersWithTrackerHitsLargerThan6):
				FillHistosEfficiency(tree,ll_efficiencies_acceptance_nom,i_particle,index,weightFactor,i)
			if(boolNGrandDaughtersWithTrackerHitsLargerThan6):
				FillHistosEfficiency(tree,ll_efficiencies_acceptance_denom,i_particle,index,weightFactor,i)
		
		#accuracies:
		if(antiSReconstructed):
			h1_AntiS_RECO_Acc_eta.Fill(tree._tpsAntiS_eta[0]-tree._tpsAntiS_bestRECO_eta[0],weightFactor)
			h1_AntiS_RECO_Acc_phi.Fill(tree._tpsAntiS_phi[0]-tree._tpsAntiS_bestRECO_phi[0],weightFactor)
			h1_AntiS_RECO_Acc_vz.Fill(tree._tpsAntiS_vz_beamspot[1]-tree._tpsAntiS_bestRECO_vz_beamspot[0],weightFactor)
			h1_AntiS_RECO_Acc_lxy.Fill(tree._tpsAntiS_Lxy_beamspot[1]-tree._tpsAntiS_bestRECO_Lxy_beamspot[0],weightFactor)
			h1_AntiS_RECO_Acc_pt.Fill(tree._tpsAntiS_pt[0]-tree._tpsAntiS_bestRECO_pt[0],weightFactor)
			h1_AntiS_RECO_Acc_pz.Fill(tree._tpsAntiS_pz[0]-tree._tpsAntiS_bestRECO_pz[0],weightFactor)

		#now select tracks which have the antiS properly reconstructed and check for these events how the kinematics look like of the final state particles.
		if(antiSReconstructed):
			h_AntiS_massMinusNeutron.Fill(tree._tpsAntiS_bestRECO_massMinusNeutron[0],weightFactor)
			momentum_Ks_plus_Lambda = np.sqrt(  np.power(tree._tpsAntiS_bestRECO_pt[1]+tree._tpsAntiS_bestRECO_pt[2],2) + np.power(tree._tpsAntiS_bestRECO_pz[1]+tree._tpsAntiS_bestRECO_pz[2],2) )
			h2_AntiS_inv_massMinusNeutron_p_Ks_plus_Lambda.Fill(tree._tpsAntiS_bestRECO_massMinusNeutron[0],momentum_Ks_plus_Lambda,weightFactor)
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
					ll_kinematics_granddaughters_of_RECO_AntiS[index][0].Fill(tree._tpsAntiS_Lxy_beamspot[j],weightFactor)	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][1].Fill(tree._tpsAntiS_vz_beamspot[j],weightFactor)	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][2].Fill(tree._tpsAntiS_dxy_beamspot[j],weightFactor)	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][3].Fill(tree._tpsAntiS_dz_beamspot[j],weightFactor)	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][4].Fill(tree._tpsAntiS_pt[j],weightFactor)	
					ll_kinematics_granddaughters_of_RECO_AntiS[index][5].Fill(tree._tpsAntiS_pz[j],weightFactor)	




#plot kinematics of final state particles for which the antiS got reconstructed
granddaughterkinematics_dir_of_RECO_AntiS = fOut.mkdir("granddaughterkinematics_of_RECO_AntiS") 
granddaughterkinematics_dir_of_RECO_AntiS.cd() 
for l in ll_kinematics_granddaughters_of_RECO_AntiS: 
        for h in l: 
                h.Write()

nHistos = len(ll_kinematics_granddaughters_of_RECO_AntiS[0])
n_granddaughters = len(ll_kinematics_granddaughters_of_RECO_AntiS)
leg_granddaugthers = ["K_{S}^{0} daughters","#bar{#Lambda}^{0}-#pi^{+}","#bar{#Lambda}^{0}-#bar{p}"]
for i in range(0,nHistos):#each list contains a list of histograms. the histos need to be overlaid one list to the other,
	c_name = "c_"+ll_kinematics_granddaughters_of_RECO_AntiS[0][i].GetName()
	c = TCanvas(c_name,"")
	legend = TLegend(0.8,0.85,0.99,0.99)
	for j in [1,2,0]: #first the soft pion of the antiLambda, then the antiproton and then the Ks daughters (this is also the order you are using at GEN level) 
		h = ll_kinematics_granddaughters_of_RECO_AntiS[j][i]
		if(i==2 or i==3):
			h.SetMaximum(h.GetMaximum()*2)
		if(i == 1):
			h.SetMaximum(h.GetMaximum()*3)
		if(h.GetSumw2N() == 0):
			h.Sumw2(kTRUE)
		#h.Scale(1./h.Integral(), "width");
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
#h_AntiS_massMinusNeutron.Scale(1./h.Integral(), "width");
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
#plot showing that the acceptance criteria is a good one
h2_allTracksMoreThan6Hits_efficiency.Write()
c= TCanvas(h2_allTracksMoreThan6Hits_efficiency.GetName(),"");
c.SetRightMargin(0.2) #make room for the tile of the z scale
if(h2_allTracksMoreThan6Hits_efficiency.GetSumw2N() == 0):
        h2_allTracksMoreThan6Hits_efficiency.Sumw2(kTRUE)
#h.Scale(1./h.Integral(), "width");
h2_allTracksMoreThan6Hits_efficiency.Draw("text")
h2_allTracksMoreThan6Hits_efficiency.SetStats(0)
CMS_lumi.CMS_lumi(c, 0, 11)
c.SaveAs(plots_output_dir+h2_allTracksMoreThan6Hits_efficiency.GetName()+".pdf")
c.Write()


#the matching criteria between GEN and RECO:
GENRECOMatchingCriteria_dir = fOut.mkdir("GENRECOMatchingCriteria")
GENRECOMatchingCriteria_dir.cd()
l_h = [h_GENRECO_matcher1_antiS,h_GENRECO_matcher1_Ks,h_GENRECO_matcher1_AntiLambda,h_GENRECO_matcher2_antiS,h_GENRECO_matcher2_Ks,h_GENRECO_matcher2_AntiLambda]
for h in l_h:
        h.SetDirectory(0)
for i in range(0,len(l_h)):
        h = l_h[i]
        c_name = "c_"+h.GetName()
        c = TCanvas(c_name,"")
    #    CMS_lumi.CMS_lumi(c, 0, 11)
        if(h.GetSumw2N() == 0):
                h.Sumw2(kTRUE)
        #h.Scale(1./h.Integral(), "width");
        h.SetLineColor(colours[0])
        h.SetLineWidth(2)
        h.SetMarkerStyle(22+0)
        h.SetMarkerColor(colours[0])
        h.SetStats(0)
        h.Draw("")
	c.SetLogy()
	CMS_lumi.CMS_lumi(c, 0, 11)
        c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
        c.Write()                        



#The reconstruction efficiencies of all the particles
ll_efficiencies = []
for i in range(0,len(ll_efficiencies_nom)):
	ll_efficiencies.append([])
	for j in range(0,len(ll_efficiencies_nom[i])):
		ll_efficiencies[i].append(TEfficiency(ll_efficiencies_nom[i][j],ll_efficiencies_denom[i][j]))
		

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
		#hc_PassedHistogram.Scale(1./hc_PassedHistogram.Integral(), "width");

		if("pt" in hc_PassedHistogram.GetName() or "_p_" in hc_PassedHistogram.GetName()):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/0.1GeV")	
		elif("pz" in hc_PassedHistogram.GetName() ):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/1GeV")	
		elif("eta" in hc_PassedHistogram.GetName()):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/0.1#eta")	
		elif("lxy" in hc_PassedHistogram.GetName()):
			hc_PassedHistogram.GetYaxis().SetTitle("Events/0.1mm")
			c.SetLogy()	
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

ll_efficiencies_acceptance = []
for i in range(0,len(ll_efficiencies_acceptance_nom)):#loop over different particles
	ll_efficiencies_acceptance.append([])
	for j in range(0,len(ll_efficiencies_acceptance_nom[i])):#loop over different kin variables
		ll_efficiencies_acceptance[i].append(TEfficiency(ll_efficiencies_acceptance_nom[i][j],ll_efficiencies_acceptance_denom[i][j]))

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

		#draw the efficiency plot	
		hc.SetLineColor(colours[1])
		hc.SetLineWidth(2)
		hc.SetMarkerStyle(22+1)
		hc.SetMarkerColor(colours[1])
		legend.AddEntry(hc,"Reconstruction efficiency","lep")
		#hc.SetMinimum(0)
		hc.Draw("")

		#now draw the distribution of that variable on top of it
		#scale the distribution a bit to make it visible on the scale of the efficiency
		hc_TotalHistogram.SetLineColor(colours[0])
		hc_TotalHistogram.SetLineWidth(2)
		hc_TotalHistogram.SetMarkerStyle(22+0)
		hc_TotalHistogram.SetMarkerColor(colours[0])
		hc_TotalHistogram.SetStats(0)
		#scale the distribution a bit to make it visible on the scale of the efficiency
		#start with the most specific ones
		if("AntiS_RECO_eff_lxy_antiS" in hc_TotalHistogram.GetName() ):
			hc_TotalHistogram.Scale(1/4700.)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x1/4700","lep")
		elif("AntiS_RECO_eff_vz_antiS" in hc_TotalHistogram.GetName() ):
			hc_TotalHistogram.Scale(0.5)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x0.5","lep")
		elif("AntiS_RECO_eff_pt_antiS" in hc_TotalHistogram.GetName() ):
			hc_TotalHistogram.Scale(1/40.)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x1/40","lep")
		elif("AntiS_RECO_eff_eta_antiS" in hc_TotalHistogram.GetName() ):
			hc_TotalHistogram.Scale(1/60.)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x1/60","lep")

		elif("Ks_RECO_eff_lxy_antiS" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(1)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x1","lep")
		elif("Ks_RECO_eff_eta_antiS" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(0.25)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x0.25","lep")

		elif("AntiLambda_RECO_eff_lxy_antiS" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(0.5)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x0.5","lep")
		elif("AntiLambda_RECO_eff_vz_antiS" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(10)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x10","lep")
		elif("AntiLambda_RECO_eff_eta_antiS" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(0.1)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x0.1","lep")
		#proceed with the more generic cases
		elif("pz" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale()
			legend.AddEntry(hc_TotalHistogram,"Normalised distribution","lep")
		elif("eta" in hc_TotalHistogram.GetName()): 
			hc_TotalHistogram.Scale(1.)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution","lep")
		elif("dxy" in hc_TotalHistogram.GetName()):
                        hc_TotalHistogram.Scale(1.)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution","lep")
		elif("vz" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(25)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x25","lep")
		elif("dz" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(20)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x20","lep")
		elif("lxy" in hc_TotalHistogram.GetName()):
			hc_TotalHistogram.Scale(2)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x2","lep")
		elif("pt" in hc_TotalHistogram.GetName() or "_p_" in hc_PassedHistogram.GetName()):
			hc_TotalHistogram.Scale(0.4)
                        legend.AddEntry(hc_TotalHistogram,"Normalised distribution x0.4","lep")

		hc_TotalHistogram.SetMinimum(0)
		hc_TotalHistogram.Draw("Y+samePCE1")

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
	#h.Scale(1./h.Integral(), "width");
	
	h.SetLineColor(colours[0])
	h.SetLineWidth(2)
	h.SetMarkerStyle(22+0)
	h.SetMarkerColor(colours[0])
	h.Fit("gaus")
	h.Draw("")
	c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
	c.Write()

#the cuts in the V0Fitter
V0FitterCuts_dir = fOut.mkdir("V0FitterCuts")
V0FitterCuts_dir.cd()
for h in l_h_V0FitterCuts:
	h.Write()

fOut.Close()




#print some conclusions:

print "Fraction of antiS with all granddaughters generating >= 7 tracker hits=: ", nAntiSWithAllGranddaughtersMoreThan6Hits, "/",nAntiS, " = ", nAntiSWithAllGranddaughtersMoreThan6Hits/nAntiS 

#summary of reconstruction efficiencies:
particles = ["antiS","Ks","AntiLambda","Ks-piplus","Ks-pineg","AntiLambda-pion","AntiLambda-antiproton"]
print "###################################################################################################"
print "########################################RECO EFF FOR RECONSTRUCTABLE ANTIS#########################"
print "###################################################################################################"

print "Ks if both daughters got reconstructed: ", nKsRECOIfBothDaughtersReco,"/",nKsTOTALIfBothDaughtersReco," = ", nKsRECOIfBothDaughtersReco/nKsTOTALIfBothDaughtersReco, "              --> so 1 out of :", nKsTOTALIfBothDaughtersReco/nKsRECOIfBothDaughtersReco
print "AntiLambda if both daughters got reconstructed: ", nAntiLambdaRECOIfBothDaughtersReco,"/",nAntiLambdaTOTALIfBothDaughtersReco," = ", nAntiLambdaRECOIfBothDaughtersReco/nAntiLambdaTOTALIfBothDaughtersReco, "              --> so 1 out of :", nAntiLambdaTOTALIfBothDaughtersReco/nAntiLambdaRECOIfBothDaughtersReco
print "AntiS if both daughters got reconstructed: ", nAntiSRECOIfBothDaughtersReco,"/",nAntiSTOTALIfBothDaughtersReco," = ", nAntiSRECOIfBothDaughtersReco/nAntiSTOTALIfBothDaughtersReco, "              --> so 1 out of :", nAntiSTOTALIfBothDaughtersReco/nAntiSRECOIfBothDaughtersReco
print '\n\n'
for i_par, nominator in enumerate(nTotalReconstructed_reconstructable,start = 0):
	if i_par == 1 or i_par == 3:
		print "------------------------------------------"
	print particles[i_par],": ", nominator,"/",nAntiS_reconstructable, " = ", nominator/nAntiS_reconstructable, "       --> so 1 out of: " , nAntiS_reconstructable/nominator
