#With this script you can make some tracking efficiency plots for generic (not necessarily coming from Sbar events) tracks 

from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, TLegend, TProfile2D, TProfile
import numpy as np

fIn = TFile('/afs/cern.ch/work/w/wvetens/Sexaquarks/CMSSW_10_2_26/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerTracking/TrackingFlatTree.root', 'read')
tree = fIn.Get('FlatTreeProducerTracking/FlatTreeTracks') 

fOut = TFile('analyzed_test_FlatTreeTracking_Step2_trial11.root','RECREATE')



eff_pt  = TEfficiency("eff_pt","; Simulated track p_{T} (GeV) ;efficiency",100,0,20)
eff_pt_cut_eta_lxy  = TEfficiency("eff_pt_cut_eta_lxy","; Simulated track p_{T} (GeV);efficiency",100,0,20)
eff_eta  = TEfficiency("eff_eta","; Simulated track #eta ;efficiency",60,-3,3)
eff_phi  = TEfficiency("eff_phi","; Simulated track #phi ;efficiency",60,-3,3)
eff_pz  = TEfficiency("eff_pz","; Simulated track p_{z} (GeV);efficiency",100,-50,50)
eff_Lxy_beamspot  = TEfficiency("eff_Lxy_beamspot","; Simulated track l_{0}(beamspot) (cm);efficiency",60,0,60)
eff_Lxy_beamspot_cut_eta_pt  = TEfficiency("eff_Lxy_beamspot_cut_eta_pt","; Simulated track l_{0}(beamspot) (cm);efficiency",60,0,60)
eff_vz_beamspot  = TEfficiency("eff_vz_beamspot","; Simulated track v_{z}(beamspot) (cm);efficiency",50,-100,100)
eff_dxy_beamspot  = TEfficiency("eff_dxy_beamspot","; Simulated track d_{0}(beamspot) (cm);efficiency",50,-10,10)
eff_dz_beamspot  = TEfficiency("eff_dz_beamspot","; Simulated track d_{z}(beamspot) (cm);efficiency",50,-100,100)
eff_numberOfTrackerLayers = TEfficiency("eff_numberOfTrackerLayers","; Simulated track number of tracker layers hit ;efficiency",50,0,50)
eff_numberOfTrackerLayers_cut_eta_pt = TEfficiency("eff_numberOfTrackerLayers_cut_eta_pt","; Simulated track number of tracker layers hit ;efficiency",50,0,50)
eff_etaOfGrandMotherAntiS = TEfficiency("eff_etaOfGrandMotherAntiS","; Simulated grandmother #bar{S} #eta;efficiency",100,-5,5)

l_eff_histos = [eff_pt,eff_pt_cut_eta_lxy,eff_eta,eff_phi,eff_pz,eff_Lxy_beamspot,eff_Lxy_beamspot_cut_eta_pt,eff_vz_beamspot,eff_dxy_beamspot,eff_dz_beamspot,eff_numberOfTrackerLayers,eff_numberOfTrackerLayers_cut_eta_pt,eff_etaOfGrandMotherAntiS]

#list of the hisos in the above list, there is one list per 'collection' of particles. The collections are: all tracks not from an antiS, all tracks from antiS-Ks, tracks from antiS-AntiLambda-pion, tracks from antiS-AntiLambda-antiProton
ll_eff_histos = []
for l in range(0,4):
	ll_eff_histos.append([])
	for h in l_eff_histos:
                ll_eff_histos[l].append(h.Clone())

#try to say something about the fact if granddaughters from the antiS loose tracking eff because they fall outside of the tracker acceptance.
tprof_etaOfGrandMotherAntiS_numberOfTrackerLayers = TProfile("tprof_etaOfGrandMotherAntiS_numberOfTrackerLayers","; Simulated grandmother #bar{S} #eta; Simulated track mean number of tracker layers hit",100,-5,5,0,20)
tprof_etaOfGrandMotherAntiS_eff = TProfile("tprof_etaOfGrandMotherAntiS_eff","; Simulated grandmother #bar{S} #eta; Efficiency",100,-5,5,0,1)
tprof2_etaOfGrandMotherAntiS_lxyz_numberOfTrackerLayers = TProfile2D("tprof2_etaOfGrandMotherAntiS_lxyz_numberOfTrackerLayers","; Simulated grandmother #bar{S} #eta;  Simulated track l_{xyz}(beamspot) (cm);Simulated track mean number of tracker layers hit",20,-5,5,26,0,120,0,20)
tprof2_etaOfGrandMotherAntiS_lxyz_eff = TProfile2D("tprof2_etaOfGrandMotherAntiS_lxyz_eff","; Simulated grandmother #bar{S} #eta;  Simulated track l_{xyz}(beamspot) (cm); Efficiency",20,-5,5,26,0,130,0,1)

prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_lowPt_lowPz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_lowPt_lowPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_lowPt_middlePz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_lowPt_middlePz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_lowPt_highPz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_lowPt_highPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_middlePt_lowPz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_middlePt_lowPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_middlePt_middlePz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_middlePt_middlePz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_middlePt_highPz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_middlePt_highPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_highPt_lowPz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_highPt_lowPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_highPt_middlePz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_highPt_middlePz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)
prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_highPt_highPz = TProfile2D("prof2_lxy_vx_creation_vertex_daughters_numberOfTrackerLayers_highPt_highPz",";v_{z} decay vertex K_{S}^{0} (cm); l_{0} decay vertex K_{S}^{0} (cm);mean #tracker layers hit by the track from K_{S}^{0} daughter",60,-300,300,120,0,120)

#also count the fraction of tracks for each type of track and how many have >= 7 hits (>= 7 is the requirement in the V0Fitter)
nTrackerLayerHits =[[0.,0.],[0.,0.],[0.,0.],[0.,0.]]

l_tprof =[tprof_etaOfGrandMotherAntiS_numberOfTrackerLayers,tprof_etaOfGrandMotherAntiS_eff,tprof2_etaOfGrandMotherAntiS_lxyz_numberOfTrackerLayers,tprof2_etaOfGrandMotherAntiS_lxyz_eff,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_lowPt_lowPz,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_lowPt_middlePz,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_lowPt_highPz,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_middlePt_lowPz,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_middlePt_middlePz,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_middlePt_highPz,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_highPt_lowPz,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_highPt_middlePz,prof2_vz_lxy_creation_vertex_daughters_numberOfTrackerLayers_highPt_highPz]
ll_tprof = []
for l in range(0,4):
	ll_tprof.append([])
	for h in l_tprof:
                ll_tprof[l].append(h.Clone())

#first do this small loop which runs over the tree containing all the antiS (i.e. also the ones which do not go the correct granddaughters)
for i in range(0,tree.GetEntries()):
	tree.GetEntry(i)

	if(i>1e5):
		break

	if(i%10000 == 0):
		print 'reached track: ', i

	#not interested in tp which are neutral
	if (tree._tp_charge[0] == 0):
		continue

	isReconstructed = tree._tp_reconstructed[0]
	isAntiSTrack = tree._tp_isAntiSTrack[0] #0 is track not from an antiS, 1 is pion from Ks from antiS, 2 is pion from antiL from antiS, 3 is antiproton from antiLambda from antiS
	lxyz = np.sqrt(tree._tp_Lxy_beamspot[0]**2 + tree._tp_vz_beamspot[0]**2)

	#the efficiency plots
	ll_eff_histos[isAntiSTrack][0].Fill(isReconstructed,tree._tp_pt[0])
	if(tree._tp_Lxy_beamspot[0] < 3.5 and abs(tree._tp_eta[0]) < 2.5):	
		ll_eff_histos[isAntiSTrack][1].Fill(isReconstructed,tree._tp_pt[0])
	ll_eff_histos[isAntiSTrack][2].Fill(isReconstructed,tree._tp_eta[0])
	ll_eff_histos[isAntiSTrack][3].Fill(isReconstructed,tree._tp_phi[0])
	ll_eff_histos[isAntiSTrack][4].Fill(isReconstructed,tree._tp_pz[0])
	ll_eff_histos[isAntiSTrack][5].Fill(isReconstructed,tree._tp_Lxy_beamspot[0])
	if(tree._tp_pt[0] > 0.9 and abs(tree._tp_eta[0]) < 2.5):
		ll_eff_histos[isAntiSTrack][6].Fill(isReconstructed,tree._tp_Lxy_beamspot[0])
	ll_eff_histos[isAntiSTrack][7].Fill(isReconstructed,tree._tp_vz_beamspot[0])
	ll_eff_histos[isAntiSTrack][8].Fill(isReconstructed,tree._tp_dxy_beamspot[0])
	ll_eff_histos[isAntiSTrack][9].Fill(isReconstructed,tree._tp_dz_beamspot[0])
	ll_eff_histos[isAntiSTrack][10].Fill(isReconstructed,tree._tp_numberOfTrackerLayers[0])
	if(tree._tp_pt[0] > 0.9 and abs(tree._tp_eta[0]) < 2.5):
		ll_eff_histos[isAntiSTrack][11].Fill(isReconstructed,tree._tp_numberOfTrackerLayers[0])
	ll_eff_histos[isAntiSTrack][12].Fill(isReconstructed,tree._tp_etaOfGrandMotherAntiS[0])	

	#the tprofiles
	ll_tprof[isAntiSTrack][0].Fill(tree._tp_etaOfGrandMotherAntiS[0],tree._tp_numberOfTrackerLayers[0])
	ll_tprof[isAntiSTrack][1].Fill(tree._tp_etaOfGrandMotherAntiS[0],isReconstructed)
	ll_tprof[isAntiSTrack][2].Fill(tree._tp_etaOfGrandMotherAntiS[0],lxyz,tree._tp_numberOfTrackerLayers[0])
	ll_tprof[isAntiSTrack][3].Fill(tree._tp_etaOfGrandMotherAntiS[0],lxyz,isReconstructed)

        if(tree._tp_pt[0]>0.35 and tree._tp_pt[0]<0.5):
                if(abs(tree._tp_pz[0]) < 1):
                        ll_tprof[isAntiSTrack][4].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])
                elif(abs(tree._tp_pz[0]) > 1 and abs(tree._tp_pz[0]) < 3):
                        ll_tprof[isAntiSTrack][5].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])
                else:
                        ll_tprof[isAntiSTrack][6].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])
        elif(tree._tp_pt[0]>0.5 and tree._tp_pt[0]<1.):
                if(abs(tree._tp_pz[0]) < 1):
                        ll_tprof[isAntiSTrack][7].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])
                elif(abs(tree._tp_pz[0]) > 1 and abs(tree._tp_pz[0]) < 3):
                        ll_tprof[isAntiSTrack][8].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])
                else:
                        ll_tprof[isAntiSTrack][9].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])
        else:
                if(abs(tree._tp_pz[0]) < 1):
                        ll_tprof[isAntiSTrack][10].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])
                elif(abs(tree._tp_pz[0]) > 1 and abs(tree._tp_pz[0]) < 3):
                        ll_tprof[isAntiSTrack][11].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])
                else:
                        ll_tprof[isAntiSTrack][12].Fill(tree._tp_vz_beamspot[0],tree._tp_Lxy_beamspot[0],tree._tp_numberOfTrackerLayers[0])

	nTrackerLayerHits[isAntiSTrack][0]+=1
	if(tree._tp_numberOfTrackerLayers[0]>=7):
		nTrackerLayerHits[isAntiSTrack][1]+=1


print "Track type \t Number of Tracks found \t Number of tracks with >= 7 TrackerLayerHits \t fraction"
print "Non-antiS tracks \t", nTrackerLayerHits[0][0], "\t", nTrackerLayerHits[0][1], "\t", nTrackerLayerHits[0][1]/nTrackerLayerHits[0][0]
print "Ks daughters \t", nTrackerLayerHits[1][0], "\t", nTrackerLayerHits[1][1], "\t", nTrackerLayerHits[1][1]/nTrackerLayerHits[1][0]
print "AntiLambda pos pion daughters \t", nTrackerLayerHits[2][0], "\t", nTrackerLayerHits[2][1], "\t", nTrackerLayerHits[2][1]/nTrackerLayerHits[2][0]
print "AntiLambda antiproton daughters \t", nTrackerLayerHits[3][0], "\t", nTrackerLayerHits[3][1], "\t", nTrackerLayerHits[3][1]/nTrackerLayerHits[3][0]

	
dir_name = ["nonAntiSTracks","AntiSKsPionTracks","AntiSAntiLambdaPionTracks","AntiSAntiLambdaAntiProtonTracks"]

i = 0
for l in ll_eff_histos:
	directory = fOut.mkdir('teff_'+dir_name[i])
	directory.cd()	
	for h in l:
		h.Write()
		PassedHisto = h.GetCopyPassedHisto()
		TotalHisto = h.GetCopyTotalHisto()
		PassedHisto.Write()
		TotalHisto.Write()
		
	i+=1
i = 0
for l in ll_tprof:
	directory = fOut.mkdir('tprof_'+dir_name[i])
	directory.cd()	
	for h in l:
		h.Write()
		
	i+=1



fOut.Close()

