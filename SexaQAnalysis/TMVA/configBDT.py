import numpy as np

#BDT classifier cut
BDT_classifier_cut = 0.35

#matching criteria for the GEN and RECO for the Ks, AntiLambda and AntiS
GENRECO_matcher_AntiS_deltaL = 2.
GENRECO_matcher_AntiS_deltaR = 0.5
GENRECO_matcher_Ks_deltaL = 2.
GENRECO_matcher_Ks_deltaR = 0.03
GENRECO_matcher_AntiL_deltaL = 3.
GENRECO_matcher_AntiL_deltaR = 0.03

#cuts to be applied pre BDT:
#some cuts which are quite obvious:
pre_BDT_noCut = "Alt$(_S_lxy_interaction_vertex_beampipeCenter,0) < 99999999." #essentially no cut
pre_BDT_cut1 =  "Alt$(_S_lxy_interaction_vertex_beampipeCenter,0) > 2.02 && Alt$(_S_lxy_interaction_vertex_beampipeCenter,0) < 2.4"
#pre_BDT_cut2 = "Alt$(_S_error_lxy_interaction_vertex_beampipeCenter,0) < 0.015"
pre_BDT_cut3 = "Alt$(_S_dxy_over_lxy,0) >= 0 && Alt$(_S_dxy_over_lxy,0) <= 0.5"
pre_BDT_cut4 = "(Alt$(_S_daughters_deltaphi,0) < -0.5 || Alt$(_S_daughters_deltaphi,0) > 0.5)"

#part of these cuts are the fiducial region in which we understand the systematics:
FiducialRegionptMin = 0.33
FiducialRegionptMax = 999999999.

FiducialRegionpzMin = -22.
FiducialRegionpzMax = 22.

FiducialRegiondxyMin = 0.
FiducialRegiondxyMax = 9.5

FiducialRegiondzMin = -27.
FiducialRegiondzMax = 27.

FiducialRegionlxyMax = 44.5

FiducialRegionvzMin = -125.
FiducialRegionvzMax = 125.



#apply the above limits on all 4 final state particles
fiducial_region_cuts = "Alt$(_Lambda_vz_decay_vertex,0) >= "+str(FiducialRegionvzMin)\
+" && Alt$(_Lambda_vz_decay_vertex,0) <= "+str(FiducialRegionvzMax)\
+" && Alt$(_Lambda_lxy_decay_vertex,0) <= "+str(FiducialRegionlxyMax)\
+" && Alt$(_Ks_vz_decay_vertex,0) >= "+str(FiducialRegionvzMin)\
+" && Alt$(_Ks_vz_decay_vertex,0) <= "+str(FiducialRegionvzMax)\
+" && Alt$(_Ks_lxy_decay_vertex,0) <= "+str(FiducialRegionlxyMax)\
+" && Alt$(_RECO_Lambda_daughter0_pt,0) >= "+str(FiducialRegionptMin)\
+" && Alt$(_RECO_Lambda_daughter0_pz,0) >= "+str(FiducialRegionpzMin)\
+" && Alt$(_RECO_Lambda_daughter0_pz,0) <= "+str(FiducialRegionpzMax)\
+" && Alt$(_RECO_Lambda_daughter0_dxy_beamspot,0) >= "+str(FiducialRegiondxyMin)\
+" && Alt$(_RECO_Lambda_daughter0_dxy_beamspot,0) <= "+str(FiducialRegiondxyMax)\
+" && Alt$(_RECO_Lambda_daughter0_dz_beamspot,0) >= "+str(FiducialRegiondzMin)\
+" && Alt$(_RECO_Lambda_daughter0_dz_beamspot,0) <= "+str(FiducialRegiondzMax)\
+" && Alt$(_RECO_Lambda_daughter1_pt,0) >= "+str(FiducialRegionptMin)\
+" && Alt$(_RECO_Lambda_daughter1_pz,0) >= "+str(FiducialRegionpzMin)\
+" && Alt$(_RECO_Lambda_daughter1_pz,0) <= "+str(FiducialRegionpzMax)\
+" && Alt$(_RECO_Lambda_daughter1_dxy_beamspot,0) >= "+str(FiducialRegiondxyMin)\
+" && Alt$(_RECO_Lambda_daughter1_dxy_beamspot,0) <= "+str(FiducialRegiondxyMax)\
+" && Alt$(_RECO_Lambda_daughter1_dz_beamspot,0) >= "+str(FiducialRegiondzMin)\
+" && Alt$(_RECO_Lambda_daughter1_dz_beamspot,0) <= "+str(FiducialRegiondzMax)\
+" && Alt$(_RECO_Ks_daughter0_pt,0) >= "+str(FiducialRegionptMin)\
+" && Alt$(_RECO_Ks_daughter0_pz,0) >= "+str(FiducialRegionpzMin)\
+" && Alt$(_RECO_Ks_daughter0_pz,0) <= "+str(FiducialRegionpzMax)\
+" && Alt$(_RECO_Ks_daughter0_dxy_beamspot,0) >= "+str(FiducialRegiondxyMin)\
+" && Alt$(_RECO_Ks_daughter0_dxy_beamspot,0) <= "+str(FiducialRegiondxyMax)\
+" && Alt$(_RECO_Ks_daughter0_dz_beamspot,0) >= "+str(FiducialRegiondzMin)\
+" && Alt$(_RECO_Ks_daughter0_dz_beamspot,0) <= "+str(FiducialRegiondzMax)\
+" && Alt$(_RECO_Ks_daughter1_pt,0) >= "+str(FiducialRegionptMin)\
+" && Alt$(_RECO_Ks_daughter1_pz,0) >= "+str(FiducialRegionpzMin)\
+" && Alt$(_RECO_Ks_daughter1_pz,0) <= "+str(FiducialRegionpzMax)\
+" && Alt$(_RECO_Ks_daughter1_dxy_beamspot,0) >= "+str(FiducialRegiondxyMin)\
+" && Alt$(_RECO_Ks_daughter1_dxy_beamspot,0) <= "+str(FiducialRegiondxyMax)\
+" && Alt$(_RECO_Ks_daughter1_dz_beamspot,0) >= "+str(FiducialRegiondzMin)\
+" && Alt$(_RECO_Ks_daughter1_dz_beamspot,0) <= "+str(FiducialRegiondzMax)


#dictionary to load in several macros
config_dict = {
"GENRECO_matcher_AntiS_deltaL":GENRECO_matcher_AntiS_deltaL,
"GENRECO_matcher_AntiS_deltaR":GENRECO_matcher_AntiS_deltaR,
"GENRECO_matcher_Ks_deltaL":GENRECO_matcher_Ks_deltaL,
"GENRECO_matcher_Ks_deltaR":GENRECO_matcher_Ks_deltaR,
"GENRECO_matcher_AntiL_deltaL":GENRECO_matcher_AntiL_deltaL,
"GENRECO_matcher_AntiL_deltaR":GENRECO_matcher_AntiL_deltaR,
"fiducial_region_cuts":fiducial_region_cuts,
"pre_BDT_noCut":pre_BDT_noCut,
"pre_BDT_cut1":pre_BDT_cut1,
#"pre_BDT_cut2":pre_BDT_cut2,
"pre_BDT_cut3":pre_BDT_cut3,
"pre_BDT_cut4":pre_BDT_cut4,
"config_pre_BDT_cuts":fiducial_region_cuts + " && "+pre_BDT_cut4 + " && "+pre_BDT_cut1 + " && " + pre_BDT_cut3,
"BDT_classifier_cut":BDT_classifier_cut ,
"config_SignalFile":"/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Final/1p8GeV/FlatTreeBDTMCSignal/combined_FlatTreeBDT_trial17AND21_1p8GeV_02112019_v1.root",
"config_BkgFileData":"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/ALL/ALL_v6/FlatTreeBDT_SingleMuon_Run2016H-07Aug17-v1_trialR.root",
"config_BkgFileMC":"/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Final/1p8GeV/FlatTreeBDTMCBackground/FlatTreeBDT_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__31102019_v1.root",
"config_SelectionSignalAntiS":"Alt$(_S_charge,0) == -1 && Alt$(_S_deltaLInteractionVertexAntiSmin,0) < "+str(GENRECO_matcher_AntiS_deltaL)+" && Alt$(_S_deltaRAntiSmin,0) <"+str(GENRECO_matcher_AntiS_deltaR),
"config_SelectionBkgS":"Alt$(_S_charge,0) == 1",
"config_SelectionBkgAntiS":"Alt$(_S_charge,0) == -1 && Alt$(_S_deltaLInteractionVertexAntiSmin,0) > 10",
"config_SelectionAntiS":"Alt$(_S_charge,0) == -1",
"config_fidRegion_FiducialRegionptMin":FiducialRegionptMin,
"config_fidRegion_FiducialRegionptMax":FiducialRegionptMax,
"config_fidRegion_FiducialRegionpzMin":FiducialRegionpzMin,
"config_fidRegion_FiducialRegionpzMax":FiducialRegionpzMax,
"config_fidRegion_FiducialRegiondxyMin":FiducialRegiondxyMin,
"config_fidRegion_FiducialRegiondxyMax":FiducialRegiondxyMax,
"config_fidRegion_FiducialRegiondzMin":FiducialRegiondzMin,
"config_fidRegion_FiducialRegiondzMax":FiducialRegiondzMax,
"config_fidRegion_FiducialRegionlxyMax":FiducialRegionlxyMax,
"config_fidRegion_FiducialRegionvzMin":FiducialRegionvzMin,
"config_fidRegion_FiducialRegionvzMax":FiducialRegionvzMax
}


#to calculate the reweighing factor:
def calc_reweighing_factor(eta_S,isMCSignal):
        reweighing_factor = 1
        if(isMCSignal):
                theta = 2*np.arctan(np.exp(-eta_S))
                reweighing_factor = 1/np.sin(theta) #the reweiging factor has to scale with the path length of the particle through the beampipe which depends on theta
        return reweighing_factor
