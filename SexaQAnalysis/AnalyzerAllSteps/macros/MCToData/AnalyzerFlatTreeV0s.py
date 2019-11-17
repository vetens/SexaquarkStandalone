from ROOT import *
import numpy as np
import copy
from array import array

import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = ""
tdrstyle.setTDRStyle()

colours = [1,2,4,35,38,41]


fIns = [
#the standard V0s:
#'/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DoubleMuon_standardV0s/combined_FlatTreeV0_DoubleMuonData_RunC_D_E_G_H.root',
#'/storage_mnt/storage/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/FlatTrees_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_standardV0s/combined_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root',

#the adapted V0s:
'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/Results/FlatTrees_DoubleMuon_adaptedV0s_extendedFlatTree8/combined_FlatTreeV0_DoubleMuonData_Run_G_H.root',
'/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/test/FlatTreeProducerV0s/qsub/Results/FlatTrees_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_adaptedV0s_extendedFlatTree8/combined/combined_FlatTreeV0_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root'
]

maxEvents = 1e99 #good one for speed versus stats is 1e5


def eff_error(data,MC):
	return data/MC*np.sqrt(1/MC+1/data)	

#checks which one is the largest in absolute terms and then returns the signed value
def maxSigned(x1, x2):
	if(abs(x1) > abs(x2)):
		return x1
	else:
		return x2

massKs = 0.493677
c = 299792458

def lifetime(pt,pz,lxy,vz,m):
	lxyz = np.sqrt(lxy*lxy+vz*vz)
	p = np.sqrt(pt*pt+pz*pz)
	pm2 = np.power(p/m,2)
	v = np.sqrt(pm2/(1 + pm2))
	t = lxyz/v
	return t

def beta(pt,pz,lxy,vz,m):
        p = np.sqrt(pt*pt+pz*pz)
        pm2 = np.power(p/m,2)
        beta = np.sqrt(pm2/(1 + pm2))
	return beta

#fist find which of the trees has the smallest number of entries. Data or MC?
fIn1 = TFile(fIns[0],'read')
fIn2 = TFile(fIns[1],'read')
treeKs1 = fIn1.Get('FlatTreeProducerV0s/FlatTreeKs')
treeKs2 = fIn2.Get('FlatTreeProducerV0s/FlatTreeKs')
min_entries = min( treeKs1.GetEntries(), treeKs2.GetEntries() )
print "the nentries in Data: " , treeKs1.GetEntries()
print "the nentries in MC: " , treeKs2.GetEntries()
print "the min_enties: ", min_entries

iFile = 0
h_RECO_PV0_vz_Data= TH1F('h_RECO_PV0_vz_Data'," ;PV vz (cm); #entries",20000,-100,100)
h_RECO_PV0_vz_MC= TH1F('h_RECO_PV0_vz_MC'," ;PV vz (cm); #entries",20000,-100,100)
for fIn in fIns:
	fIn = TFile(fIn,'read')
	treePV  = fIn.Get('FlatTreeProducerV0s/FlatTreePV')
	for i in range(0,min_entries):
		if(i > maxEvents):
                        continue
		treePV.GetEntry(i)
		if(iFile == 0):
			h_RECO_PV0_vz_Data.Fill(treePV._PV0_vz[0])	
		else:	
			h_RECO_PV0_vz_MC.Fill(treePV._PV0_vz[0])
	iFile = iFile+1	

iFile = 0
h_RECO_beamspot_vz_Data= TH1F('h_RECO_beamspot_vz_Data'," ;PV vz (cm); #entries",2000,-2,2)
h_RECO_beamspot_vz_MC= TH1F('h_RECO_beamspot_vz_MC'," ;PV vz (cm); #entries",2000,-2,2)
h_RECO_beamspot_vz_Data.SetDirectory(0)
h_RECO_beamspot_vz_MC.SetDirectory(0)
for fIn in fIns:
	fIn = TFile(fIn,'read')
	treeBeamspot  = fIn.Get('FlatTreeProducerV0s/FlatTreeBeamspot')
	for i in range(0,min_entries):
		if(i > maxEvents):
                        continue
		treeBeamspot.GetEntry(i)
		if(iFile == 0):
			h_RECO_beamspot_vz_Data.Fill(treeBeamspot._beampot_vz[0])	
		else:	
			h_RECO_beamspot_vz_MC.Fill(treeBeamspot._beampot_vz[0])
	iFile = iFile+1	

iFile = 0
h_RECO_Ks_vz_PV_min_Data= TH1F('h_RECO_Ks_vz_PV_min_Data'," ;PV vz min (cm); #entries",4000,-20,20)
h_RECO_Ks_vz_PV_min_MC= TH1F('h_RECO_Ks_vz_PV_min_MC'," ;PV vz min (cm); #entries",4000,-20,20)
h_RECO_track_multiplicity_Data = TH1F('h_RECO_track_multiplicity_Data',';track multiplicity;#entries',250,0,2000)
h_RECO_track_multiplicity_MC = TH1F('h_RECO_track_multiplicity_MC',';track multiplicity;#entries',250,0,2000)
h2_RECO_Ks_vz_PV_min_track_multiplicity_Data = TH2F('h2_RECO_Ks_vz_PV_min_track_multiplicity_Data',';PV vz min (cm);track multiplicity;#entries',40,-20,20,20,0,2000)
h2_RECO_Ks_vz_PV_min_track_multiplicity_MC = TH2F('h2_RECO_Ks_vz_PV_min_track_multiplicity_MC',';PV vz min (cm);track multiplicity;#entries',40,-20,20,20,0,2000)

h_RECO_Ks_vz_PV_min_Data.SetDirectory(0)
h_RECO_Ks_vz_PV_min_MC.SetDirectory(0)
h_RECO_track_multiplicity_Data.SetDirectory(0)
h_RECO_track_multiplicity_MC.SetDirectory(0)
h2_RECO_Ks_vz_PV_min_track_multiplicity_Data.SetDirectory(0)
h2_RECO_Ks_vz_PV_min_track_multiplicity_MC.SetDirectory(0)
for fIn in fIns:
	fIn = TFile(fIn,'read')
	treeKs = fIn.Get('FlatTreeProducerV0s/FlatTreeKs')
	treeGeneral = fIn.Get('FlatTreeProducerV0s/FlatTreeGeneral')
	for i in range(0,min_entries):
		if(i > maxEvents):
                        continue

		treeKs.GetEntry(i)
		treeGeneral.GetEntry(i)

		if(treeGeneral._general_triggerFired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ[0] == 0 and treeGeneral._general_triggerFired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ[0] == 0):
			continue

		for j in range(0, len(treeKs._Ks_mass)):
			#require only Ks which point to a 'best PV'
                        if(treeKs._Ks_mass[j] > 0.48 and treeKs._Ks_mass[j] < 0.52 and  abs(treeKs._Ks_eta[j]) < 2 and treeKs._Ks_dxy_beamspot[j] < 0.1 and treeKs._Ks_dxy_beamspot[j] > 0. and abs(treeKs._Ks_dz_min_PV[j]) < 0.2) :
				if(iFile == 0):
					h_RECO_Ks_vz_PV_min_Data.Fill(treeKs._Ks_vz_dz_min_PV[j])
		#			h_RECO_track_multiplicity_Data.Fill(treeGeneral._general_eventTrackMultiplicity[0])#some events with 2 or more good Ks will contribute twice or more to this histo, but that is ok. Such events are twice as important
					#h2_RECO_Ks_vz_PV_min_track_multiplicity_Data.Fill(treeKs._Ks_vz_dz_min_PV[j],treeGeneral._general_eventTrackMultiplicity[0])	
				else:	
					h_RECO_Ks_vz_PV_min_MC.Fill(treeKs._Ks_vz_dz_min_PV[j])
		#			h_RECO_track_multiplicity_MC.Fill(treeGeneral._general_eventTrackMultiplicity[0])
		#			h2_RECO_Ks_vz_PV_min_track_multiplicity_MC.Fill(treeKs._Ks_vz_dz_min_PV[j],treeGeneral._general_eventTrackMultiplicity[0])
	iFile = iFile+1	

#normalize the histogram of PV in data to later use it for reweighting 
#h_RECO_PV0_vz_Data.Scale(1./h_RECO_PV0_vz_Data.Integral(), "width")
#h_RECO_PV0_vz_MC.Scale(1./h_RECO_PV0_vz_MC.Integral(), "width")

plots_output_dir = 'Results/'
fOut = TFile(plots_output_dir+'output_DataToMC_RunG_H_with_dxy_dz_min_PV_reweighing_on_Ks_vz_Dz_min_PV.root','RECREATE')
h_dummy =  TH1F('h_dummy'," x; y",20,0,10)
histos_Data = [h_dummy]*9
histos_MC = [h_dummy]*9

histos_lifetimes_Data = []
histos_lifetimes_MC = []

nKshortPerEventData = []
nKshortPerEventMC = []

h2_KsLifetime_vz_lxy_Data = TH2F("h2_KsLifetime_vz_lxy_Data",";v_{z} K_{S} decay (cm);l_{0} K_{S} decay (cm);#entries",50,-100,100,20,0,100)
h2_KsLifetime_vz_lxy_MC = TH2F("h2_KsLifetime_vz_lxy_MC",";v_{z} K_{S} decay (cm);l_{0} K_{S} decay (cm);#entries",50,-100,100,20,0,100)
lxyz_Ks_max = 30
lxyz_Ks_bins = 300
h_Ks_reconstruction_MC_passed_1 = TH1F('h_Ks_reconstruction_MC_passed_1',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_2 = TH1F('h_Ks_reconstruction_MC_passed_2',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_3 = TH1F('h_Ks_reconstruction_MC_passed_3',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_4 = TH1F('h_Ks_reconstruction_MC_passed_4',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_5 = TH1F('h_Ks_reconstruction_MC_passed_5',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_6 = TH1F('h_Ks_reconstruction_MC_passed_6',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_7 = TH1F('h_Ks_reconstruction_MC_passed_7',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_8 = TH1F('h_Ks_reconstruction_MC_passed_8',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_9 = TH1F('h_Ks_reconstruction_MC_passed_9',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_10 = TH1F('h_Ks_reconstruction_MC_passed_10',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_11 = TH1F('h_Ks_reconstruction_MC_passed_11',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_12 = TH1F('h_Ks_reconstruction_MC_passed_12',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_13 = TH1F('h_Ks_reconstruction_MC_passed_13',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_14 = TH1F('h_Ks_reconstruction_MC_passed_14',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_15 = TH1F('h_Ks_reconstruction_MC_passed_15',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_16 = TH1F('h_Ks_reconstruction_MC_passed_16',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_17 = TH1F('h_Ks_reconstruction_MC_passed_17',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_18 = TH1F('h_Ks_reconstruction_MC_passed_18',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_19 = TH1F('h_Ks_reconstruction_MC_passed_19',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_passed_20 = TH1F('h_Ks_reconstruction_MC_passed_20',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
histos_Ks_reconstruction_MC_passed = [h_Ks_reconstruction_MC_passed_1,h_Ks_reconstruction_MC_passed_2,h_Ks_reconstruction_MC_passed_3,h_Ks_reconstruction_MC_passed_4,h_Ks_reconstruction_MC_passed_5,h_Ks_reconstruction_MC_passed_6,h_Ks_reconstruction_MC_passed_7,h_Ks_reconstruction_MC_passed_8,h_Ks_reconstruction_MC_passed_9,h_Ks_reconstruction_MC_passed_10,h_Ks_reconstruction_MC_passed_11,h_Ks_reconstruction_MC_passed_12,h_Ks_reconstruction_MC_passed_13,h_Ks_reconstruction_MC_passed_14,h_Ks_reconstruction_MC_passed_15,h_Ks_reconstruction_MC_passed_16,h_Ks_reconstruction_MC_passed_17,h_Ks_reconstruction_MC_passed_18,h_Ks_reconstruction_MC_passed_19,h_Ks_reconstruction_MC_passed_20]
h_Ks_reconstruction_MC_total_1 = TH1F('h_Ks_reconstruction_MC_total_1',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_2 = TH1F('h_Ks_reconstruction_MC_total_2',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_3 = TH1F('h_Ks_reconstruction_MC_total_3',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_4 = TH1F('h_Ks_reconstruction_MC_total_4',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_5 = TH1F('h_Ks_reconstruction_MC_total_5',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_6 = TH1F('h_Ks_reconstruction_MC_total_6',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_7 = TH1F('h_Ks_reconstruction_MC_total_7',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_8 = TH1F('h_Ks_reconstruction_MC_total_8',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_9 = TH1F('h_Ks_reconstruction_MC_total_9',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_10 = TH1F('h_Ks_reconstruction_MC_total_10',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_11 = TH1F('h_Ks_reconstruction_MC_total_11',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_12 = TH1F('h_Ks_reconstruction_MC_total_12',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_13 = TH1F('h_Ks_reconstruction_MC_total_13',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_14 = TH1F('h_Ks_reconstruction_MC_total_14',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_15 = TH1F('h_Ks_reconstruction_MC_total_15',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_16 = TH1F('h_Ks_reconstruction_MC_total_16',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_17 = TH1F('h_Ks_reconstruction_MC_total_17',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_18 = TH1F('h_Ks_reconstruction_MC_total_18',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_19 = TH1F('h_Ks_reconstruction_MC_total_19',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
h_Ks_reconstruction_MC_total_20 = TH1F('h_Ks_reconstruction_MC_total_20',";c.t K_{s} cm;#entries",lxyz_Ks_bins,0,lxyz_Ks_max)
histos_Ks_reconstruction_MC_total = [h_Ks_reconstruction_MC_total_1,h_Ks_reconstruction_MC_total_2,h_Ks_reconstruction_MC_total_3,h_Ks_reconstruction_MC_total_4,h_Ks_reconstruction_MC_total_5,h_Ks_reconstruction_MC_total_6,h_Ks_reconstruction_MC_total_7,h_Ks_reconstruction_MC_total_8,h_Ks_reconstruction_MC_total_9,h_Ks_reconstruction_MC_total_10,h_Ks_reconstruction_MC_total_11,h_Ks_reconstruction_MC_total_12,h_Ks_reconstruction_MC_total_13,h_Ks_reconstruction_MC_total_14,h_Ks_reconstruction_MC_total_15,h_Ks_reconstruction_MC_total_16,h_Ks_reconstruction_MC_total_17,h_Ks_reconstruction_MC_total_18,h_Ks_reconstruction_MC_total_19,h_Ks_reconstruction_MC_total_20]

#f_poisson = TF1("f_poisson","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, lxyz_Ks_max)
f_poisson = TF1("f_poisson","[0]*(TMath::Exp(-x/[1]))", 0, lxyz_Ks_max)

h2_RECO_Ks_vz_lxy_data = TH2F("h2_RECO_Ks_vz_lxy_data",";v_{z} K_{S} decay (cm);l_{0} K_{S} decay (cm);#entries",50,-100,100,20,0,100)
h2_RECO_Ks_vz_lxy_mc   = TH2F("h2_RECO_Ks_vz_lxy_mc",";v_{z} K_{S} decay (cm);l_{0} K_{S} decay (cm);#entries",50,-100,100,20,0,100)

h_Ks_mass_data = TH1F("h_Ks_mass_data",";mass K_{s} (GeV/c^{2});Entries/1MeV/c^{2}",35,0.48,0.515)
h_Ks_mass_mc = TH1F("h_Ks_mass_mc",";mass K_{s} (GeV/c^{2});Entries/1MeV/c^{2}",35,0.48,0.515)

iFile = 0
for fIn in fIns:

	print "creating plots for file: ", fIn

	fIn = TFile(fIn,'read')

	treeKs = fIn.Get('FlatTreeProducerV0s/FlatTreeKs') 
	treeZ = fIn.Get('FlatTreeProducerV0s/FlatTreeZ')
	treePV = fIn.Get('FlatTreeProducerV0s/FlatTreePV')
	treeBeamspot = fIn.Get('FlatTreeProducerV0s/FlatTreeBeamspot')
	treeGeneral = fIn.Get('FlatTreeProducerV0s/FlatTreeGeneral')

	nEntriesKs = treeKs.GetEntries()
	nEntriesZ = treeZ.GetEntries()
	nEntriesPV = treePV.GetEntries()
	nEntriesBeamspot = treeBeamspot.GetEntries()
	nEntriesGeneral = treeGeneral.GetEntries()

	print 'Number of entries in the treeKs: ', nEntriesKs
	print 'Number of entries in the treeZ: ', nEntriesZ
	print 'Number of entries in the treePV: ', nEntriesPV
	print 'Number of entries in the treeBeamspot: ', nEntriesBeamspot
	print 'Number of entries in the treeGeneral: ', nEntriesGeneral

	bins1 = np.arange(-150,-60,10)
	bins2 = np.arange(-60,60,1)
	bins3 = np.arange(60,150,10)
	args = (bins1,bins2,bins3)
	bins = np.concatenate(args)
	h_RECO_Ks_vz = TH1F('h_RECO_Ks_vz'," ;K_{s} daughter 1 and 2 absolute track v_{z} (cm); #entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_vz_beamspot = TH1F('h_RECO_Ks_vz_beamspot'," ;K_{s} daughter 1 and 2 track v_{z,bs} (cm); #entries",len(bins)-1,array('d',bins))

	bins1 = np.arange(0,50,1)
	bins2 = np.arange(50,56,1.5)
	bins3 = np.arange(56,65,1.8)
	bins4 = np.arange(65,70,2.5)
	bins5 = np.arange(70,90,5)
	args = (bins1,bins2,bins3,bins4,bins5)
	bins = np.concatenate(args)
	h_RECO_Ks_lxy = TH1F('h_RECO_Ks_lxy'," ;K_{s} daughter 1 and 2 track l_{0,bs} (cm);#entries",len(bins)-1,array('d',bins))

	bins0 = np.arange(0,0.8,0.02)
	bins1 = np.arange(0.8,2,0.1)
	bins2 = np.arange(2,8,0.2)
	bins3 = np.arange(8,10.1,0.3)
	args = (bins0,bins1,bins2,bins3)
	bins = np.concatenate(args)
	h_RECO_Ks_pt = TH1F('h_RECO_Ks_pt',"; Ks p_{t} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_pt_tracks1 = TH1F('h_RECO_Ks_pt_tracks1',"; Ks daughter 1 p_{t} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_pt_tracks2 = TH1F('h_RECO_Ks_pt_tracks2',"; Ks daughter 2 p_{t} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_pt_tracks1_and_2 = TH1F('h_RECO_Ks_pt_tracks1_and_2',"; K_{s} daughter 1 and 2 track p_{t} (GeV);#entries",len(bins)-1,array('d',bins))

	bins1 = np.arange(-40,-25,5)
	bins2 = np.arange(-25,-21,2)
	bins3 = np.arange(-21,-15,1.5)
	bins4 = np.arange(-15,15,1)
	bins5 = np.arange(15,21,1.5)
	bins6 = np.arange(21,25,2)
	bins7 = np.arange(25,40,5)
	args = (bins1,bins2,bins3,bins4,bins5,bins6,bins7)
	bins = np.concatenate(args)
	h_RECO_Ks_pz= TH1F('h_RECO_Ks_pz',"; Ks p_{z} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_pz_tracks1 = TH1F('h_RECO_Ks_pz_tracks1',"; Ks daughter 1 p_{z} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_pz_tracks2 = TH1F('h_RECO_Ks_pz_tracks2',"; Ks daughter 2 p_{z} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_pz_tracks1_and_2 = TH1F('h_RECO_Ks_pz_tracks1_and_2',"; K_{s} daughter 1 and 2 track p_{z} (GeV);#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_dxy_PV= TH1F('h_RECO_Ks_dxy_PV',"; Ks d_{0} (cm); #entries",100,-10,10)

	bins = [ -30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30 ]
	bins1 = np.arange(-30,-20,2)
	bins2 = np.arange(-20,-14,1.5)
	bins3 = np.arange(-14,14,1)
	bins4 = np.arange(14,20,1.5)
	bins5 = np.arange(20,30,2)
	args = (bins1,bins2,bins3,bins4,bins5)
	bins = np.concatenate(args)
	h_RECO_Ks_dz= TH1F('h_RECO_Ks_dz',";  Ks d_{z}(0,0,0) (cm); #entries",len(bins)-1,array('d',bins))

	bins = [ 0, .5, 1, 1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,8,9,10,12,14,16 ]
	h_RECO_Ks_Track1_dxy_beamspot= TH1F('h_RECO_Ks_Track1_dxy_beamspot',";d_{0}(beamspot) track1 (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_Track2_dxy_beamspot= TH1F('h_RECO_Ks_Track2_dxy_beamspot',";d_{0}(beamspot) track2 (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_Track1_and_2_dxy_beamspot= TH1F('h_RECO_Ks_Track1_and_2_dxy_beamspot',";K_{s} daughter 1 and 2 track d_{0,bs} (cm) ;#entries",len(bins)-1,array('d',bins))

	bins = [ -40,-35,-30,-28,-26,-24,-22,-20,-18,-16,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,30,35,40 ]
	h_RECO_Ks_Track1Track2_max_dz_min_PV = TH1F('h_RECO_Ks_Track1Track2_max_dz_min_PV',"; maxSigned[d_{z}(best PV Ks) track1, d_{z}(best PV Ks) track2]   Ks (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_Track1Track2_max_dz_beamspot = TH1F('h_RECO_Ks_Track1Track2_max_dz_beamspot',"; maxSigned[d_{z}(beamspot) track1, d_{z}(beamspot) track2]   Ks (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_Track1_dz_min_PV = TH1F('h_RECO_Ks_Track1_dz_min_PV',"; d_{z}(best PV Ks) track1 (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_Track2_dz_min_PV = TH1F('h_RECO_Ks_Track2_dz_min_PV',"; d_{z}(best PV Ks) track2 (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_Track1_and_2_dz_min_PV = TH1F('h_RECO_Ks_Track1_and_2_dz_min_PV',";K_{s} daughter 1 and 2 track d_{z}(best PV K_{s}) (cm) ;#entries",len(bins)-1,array('d',bins))
	h_RECO_Ks_Track1_and_2_dz_beamspot = TH1F('h_RECO_Ks_Track1_and_2_dz_beamspot',";K_{s} daughter 1 and 2 track d_{z,bs} track 1 and 2 (cm);#entries",len(bins)-1,array('d',bins))


	histos= [
	h_RECO_Ks_vz, 
	h_RECO_Ks_vz_beamspot, 
	h_RECO_Ks_lxy, 
	h_RECO_Ks_pt, 
	h_RECO_Ks_pt_tracks1, 
	h_RECO_Ks_pt_tracks2,
	h_RECO_Ks_pt_tracks1_and_2, 
	h_RECO_Ks_pz, 
	h_RECO_Ks_pz_tracks1, 
	h_RECO_Ks_pz_tracks2, 
	h_RECO_Ks_pz_tracks1_and_2,
	h_RECO_Ks_dxy_PV, 
	h_RECO_Ks_dz, 
	h_RECO_Ks_Track1_dxy_beamspot,
	h_RECO_Ks_Track2_dxy_beamspot,
	h_RECO_Ks_Track1_and_2_dxy_beamspot,
	h_RECO_Ks_Track1Track2_max_dz_min_PV,
	h_RECO_Ks_Track1Track2_max_dz_beamspot,
	h_RECO_Ks_Track1_dz_min_PV,
	h_RECO_Ks_Track2_dz_min_PV,
	h_RECO_Ks_Track1_and_2_dz_min_PV,
	h_RECO_Ks_Track1_and_2_dz_beamspot,
	] 
	#check the reweighting
	h_RECO_PV0_vz_MC_reweighted_to_data = TH1F('h_RECO_PV0_vz_MC_reweighted_to_data'," ;PV v_{z} (cm); #entries",200,-100,100)
	h_RECO_PV0_vz_MC_NOTreweighted_to_data = TH1F('h_RECO_PV0_vz_MC_NOTreweighted_to_data',"; PV v_{z} (cm); #entries",200,-100,100)
	#check the reweighing factor versus the _Ks_vz_dz_min_PV
	tprof_Ks_vz_dz_min_PV_reweighing_factor = TProfile('tprof_Ks_vz_dz_min_PV_reweighing_factor',"; _Ks_vz_dz_min_PV (cm); reweighing factor",100,-50,50,0,2)
	tprof_Ks_vz_reweighing_factor = TProfile('tprof_Ks_vz_reweighing_factor',"; _Ks_vz (cm); reweighing factor",300,-150,150,0,2)
	#check the reweighing of the beamspot
	h_RECO_beamspot_vz_MC_reweighted_to_data = TH1F('h_RECO_beamspot_vz_MC_reweighted_to_data'," ;beamspot v_{z} (cm); #entries",200,-100,100)
	h_RECO_beamspot_vz_MC_NOTreweighted_to_data = TH1F('h_RECO_beamspot_vz_MC_NOTreweighted_to_data',"; beamspot v_{z} (cm); #entries",200,-100,100)

	#check the reweighing of the beamspot
	h_RECO_Ks_vz_PV_min_MC_reweighted_to_data = TH1F('h_RECO_Ks_vz_PV_min_MC_reweighted_to_data'," ; v_{z} min PV (cm); #entries",40,-20,20)
	h_RECO_Ks_vz_PV_min_MC_NOTreweighted_to_data = TH1F('h_RECO_Ks_vz_PV_min_MC_NOTreweighted_to_data',"; v_{z} min PV (cm); #entries",40,-20,20)

	#some extra plots to check lifetimes of Ks 
	h_lxyz_Ks_1 = TH1F('h_lxyz_Ks_eta_1'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_2 = TH1F('h_lxyz_Ks_eta_2'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_3 = TH1F('h_lxyz_Ks_eta_3'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_4 = TH1F('h_lxyz_Ks_eta_4'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_5 = TH1F('h_lxyz_Ks_eta_5'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_6 = TH1F('h_lxyz_Ks_eta_6'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_7 = TH1F('h_lxyz_Ks_eta_7'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_8 = TH1F('h_lxyz_Ks_eta_8'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_9 = TH1F('h_lxyz_Ks_eta_9'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_10 = TH1F('h_lxyz_Ks_eta_10'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_11 = TH1F('h_lxyz_Ks_eta_11'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_12 = TH1F('h_lxyz_Ks_eta_12'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_13 = TH1F('h_lxyz_Ks_eta_13'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_14 = TH1F('h_lxyz_Ks_eta_14'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_15 = TH1F('h_lxyz_Ks_eta_15'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_16 = TH1F('h_lxyz_Ks_eta_16'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_17 = TH1F('h_lxyz_Ks_eta_17'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_18 = TH1F('h_lxyz_Ks_eta_18'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_19 = TH1F('h_lxyz_Ks_eta_19'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	h_lxyz_Ks_20 = TH1F('h_lxyz_Ks_eta_20'," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)


	histos_lifetime = [
	h_lxyz_Ks_1,
	h_lxyz_Ks_2,
	h_lxyz_Ks_3,
	h_lxyz_Ks_4,
	h_lxyz_Ks_5,
	h_lxyz_Ks_6,
	h_lxyz_Ks_7,
	h_lxyz_Ks_8,
	h_lxyz_Ks_9,
	h_lxyz_Ks_10,
	h_lxyz_Ks_11,
	h_lxyz_Ks_12,
	h_lxyz_Ks_13,
	h_lxyz_Ks_14,
	h_lxyz_Ks_15,
	h_lxyz_Ks_16,
	h_lxyz_Ks_17,
	h_lxyz_Ks_18,
	h_lxyz_Ks_19,
	h_lxyz_Ks_20
	]



	for i in range(0,min_entries):
		if(i > maxEvents):
			continue
		if(i%1e5 == 0):
			print 'reached event ', i
	
		treeZ.GetEntry(i)
		treeKs.GetEntry(i)
		treePV.GetEntry(i)
		treeBeamspot.GetEntry(i)
		treeGeneral.GetEntry(i)
		
		nKshortThisEvent = 0
	

		#only select a certain hard scale
#		if(treeZ._Z_ptMuMu[0] < 5):
#		       continue

		#trigger requirement
		if(treeGeneral._general_triggerFired_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ[0] == 0 and treeGeneral._general_triggerFired_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ[0] == 0):
			continue

		
		#now loop over all the Ks in this event:
		for j in range(0, len(treeKs._Ks_mass)):
			
			#require only Ks coming from the first PV
			#if(treeKs._Ks_mass[j] > 0.48 and treeKs._Ks_mass[j] < 0.52 and  abs(treeKs._Ks_eta[j]) < 2 and treeKs._Ks_dxy_beamspot[j] < 0.1 and treeKs._Ks_dxy_beamspot[j] > 0. and abs(treeKs._Ks_dz_PV0[j]) < 0.2) :
			#require only Ks which point to a 'best PV'
			if(treeKs._Ks_mass[j] > 0.48 and treeKs._Ks_mass[j] < 0.52 and  abs(treeKs._Ks_eta[j]) < 2 and treeKs._Ks_dxy_beamspot[j] < 0.1 and treeKs._Ks_dxy_beamspot[j] > 0. and abs(treeKs._Ks_dz_min_PV[j]) < 0.2) :
			#require only Ks which point away in dz from any PV
			#if(treeKs._Ks_mass[j] > 0.48 and treeKs._Ks_mass[j] < 0.52 and  abs(treeKs._Ks_eta[j]) < 2 and treeKs._Ks_dxy_beamspot[j] < 0.1 and treeKs._Ks_dxy_beamspot[j] > 0. and abs(treeKs._Ks_dz_min_PV[j]) > 0.2) :
			#require 'PU Ks' by asking for a Ks with large _Ks_dz_PV0, but small _Ks_dz_min_PV
			#if(treeKs._Ks_mass[j] > 0.48 and treeKs._Ks_mass[j] < 0.52 and  abs(treeKs._Ks_eta[j]) < 2 and treeKs._Ks_dxy_beamspot[j] < 0.1 and treeKs._Ks_dxy_beamspot[j] > 0. and abs(treeKs._Ks_dz_PV0[j]) > 0.2 and abs(treeKs._Ks_dz_min_PV[j]) < 0.2) :

				#if( abs(treeKs._Ks_eta[j]) > 1): continue


				#1st reweigihing factor
				x_axis_h_RECO_Ks_vz_PV_min_Data = h_RECO_Ks_vz_PV_min_Data.GetXaxis()
				#y_axis_h_RECO_Ks_vz_PV_min_Data = h2_RECO_Ks_vz_PV_min_track_multiplicity_Data.GetYaxis()
				binx_prob_h_RECO_PV0_vz_Data = x_axis_h_RECO_Ks_vz_PV_min_Data.FindBin(treeKs._Ks_vz_dz_min_PV[j])  #treeKs._Ks_vz_dz_min_PV[j]
				#biny_prob_h_RECO_PV0_vz_Data = y_axis_h_RECO_Ks_vz_PV_min_Data.FindBin(treeGeneral._general_eventTrackMultiplicity[0])  #treeKs._Ks_vz_dz_min_PV[j]
				prob_h_RECO_PV0_vz_Data = h_RECO_Ks_vz_PV_min_Data.GetBinContent(binx_prob_h_RECO_PV0_vz_Data) #the actual probability

				x_axis_h_RECO_Ks_vz_PV_min_MC = h_RECO_Ks_vz_PV_min_MC.GetXaxis()
				#y_axis_h_RECO_Ks_vz_PV_min_MC = h2_RECO_Ks_vz_PV_min_track_multiplicity_MC.GetYaxis()
				binx_prob_h_RECO_PV0_vz_MC= x_axis_h_RECO_Ks_vz_PV_min_MC.FindBin(treeKs._Ks_vz_dz_min_PV[j]) #treeKs._Ks_vz_dz_min_PV[j]
				#biny_prob_h_RECO_PV0_vz_MC= y_axis_h_RECO_Ks_vz_PV_min_MC.FindBin(treeGeneral._general_eventTrackMultiplicity[0]) #treeKs._Ks_vz_dz_min_PV[j]
				prob_h_RECO_PV0_vz_MC= h_RECO_Ks_vz_PV_min_MC.GetBinContent(binx_prob_h_RECO_PV0_vz_MC) #the actual probability

				PV_Z_reweighting_factor = 1
				#if(prob_h_RECO_PV0_vz_MC > 0 and iFile == 1): #only reweigh for MC
				#	PV_Z_reweighting_factor = prob_h_RECO_PV0_vz_Data/prob_h_RECO_PV0_vz_MC

				if(iFile == 1):
					if(prob_h_RECO_PV0_vz_MC > 0): PV_Z_reweighting_factor = prob_h_RECO_PV0_vz_Data/prob_h_RECO_PV0_vz_MC
					else: PV_Z_reweighting_factor = 0

				nKshortThisEvent+=PV_Z_reweighting_factor

				#2nd reweighing factor
				x_axis_h_RECO_beamspot_vz_Data = h_RECO_beamspot_vz_Data.GetXaxis()
				bin_prob_h_RECO_beamspot_vz_Data = x_axis_h_RECO_beamspot_vz_Data.FindBin(treeBeamspot._beampot_vz[0]) 
				prob_h_RECO_beamspot_vz_Data = h_RECO_beamspot_vz_Data.GetBinContent(bin_prob_h_RECO_beamspot_vz_Data) #the actual probability

				x_axis_h_RECO_beamspot_vz_MC= h_RECO_beamspot_vz_MC.GetXaxis()
				bin_prob_h_RECO_beamspot_vz_MC= x_axis_h_RECO_beamspot_vz_MC.FindBin(treeBeamspot._beampot_vz[0]) 
				prob_h_RECO_beamspot_vz_MC= h_RECO_beamspot_vz_MC.GetBinContent(bin_prob_h_RECO_beamspot_vz_MC) #the actual probability

				beamspot_Z_reweighting_factor = 1
				if(prob_h_RECO_beamspot_vz_MC > 0 and iFile == 1): #only reweigh for MC
					beamspot_Z_reweighting_factor = prob_h_RECO_beamspot_vz_Data/prob_h_RECO_beamspot_vz_MC

				#if(iFile == 1):
				#	print 'PV_Z_reweighting_factor: ' , str(PV_Z_reweighting_factor)
				#	print 'beamspot_Z_reweighting_factor: ' , str(beamspot_Z_reweighting_factor)

				if(iFile == 1): #for MC
					h_RECO_PV0_vz_MC_reweighted_to_data.Fill(treeKs._Ks_vz_dz_min_PV[j], PV_Z_reweighting_factor)
					h_RECO_PV0_vz_MC_NOTreweighted_to_data.Fill(treeKs._Ks_vz_dz_min_PV[j], 1)

					tprof_Ks_vz_dz_min_PV_reweighing_factor.Fill(treeKs._Ks_vz_dz_min_PV[j],PV_Z_reweighting_factor)
					tprof_Ks_vz_reweighing_factor.Fill(treeKs._Ks_vz[j],PV_Z_reweighting_factor)

					h_RECO_beamspot_vz_MC_reweighted_to_data.Fill(treeBeamspot._beampot_vz[0], beamspot_Z_reweighting_factor)
					h_RECO_beamspot_vz_MC_NOTreweighted_to_data.Fill(treeBeamspot._beampot_vz[0], 1)

					h_RECO_Ks_vz_PV_min_MC_reweighted_to_data.Fill(treeKs._Ks_vz_dz_min_PV[j], PV_Z_reweighting_factor)
					h_RECO_Ks_vz_PV_min_MC_NOTreweighted_to_data.Fill(treeKs._Ks_vz_dz_min_PV[j], 1)


				histos[0].Fill(treeKs._Ks_vz[j], PV_Z_reweighting_factor)
				histos[1].Fill(treeBeamspot._beampot_vz[0]-treeKs._Ks_vz[j], PV_Z_reweighting_factor)

				histos[2].Fill(treeKs._Ks_Lxy[j], PV_Z_reweighting_factor)

				histos[3].Fill(treeKs._Ks_pt[j], PV_Z_reweighting_factor)
				histos[4].Fill(treeKs._Ks_daughterTrack1_pt[j],PV_Z_reweighting_factor) 
				histos[5].Fill(treeKs._Ks_daughterTrack2_pt[j],PV_Z_reweighting_factor)
				histos[6].Fill(treeKs._Ks_daughterTrack1_pt[j],PV_Z_reweighting_factor)
				histos[6].Fill(treeKs._Ks_daughterTrack2_pt[j],PV_Z_reweighting_factor)

				histos[7].Fill(treeKs._Ks_pz[j], PV_Z_reweighting_factor)
				histos[8].Fill(treeKs._Ks_daughterTrack1_pz[j], PV_Z_reweighting_factor)
				histos[9].Fill(treeKs._Ks_daughterTrack2_pz[j], PV_Z_reweighting_factor)
				histos[10].Fill(treeKs._Ks_daughterTrack1_pz[j], PV_Z_reweighting_factor)
				histos[10].Fill(treeKs._Ks_daughterTrack2_pz[j], PV_Z_reweighting_factor)

				histos[11].Fill(treeKs._Ks_dxy_beamspot[j], PV_Z_reweighting_factor)

				histos[12].Fill(treeKs._Ks_dz_000[j], PV_Z_reweighting_factor)

				histos[13].Fill(treeKs._Ks_daughterTrack1_dxy_beamspot[j], PV_Z_reweighting_factor)
				histos[14].Fill(treeKs._Ks_daughterTrack2_dxy_beamspot[j], PV_Z_reweighting_factor)
				histos[15].Fill(treeKs._Ks_daughterTrack1_dxy_beamspot[j], PV_Z_reweighting_factor)
				histos[15].Fill(treeKs._Ks_daughterTrack2_dxy_beamspot[j], PV_Z_reweighting_factor)

				histos[16].Fill(maxSigned(treeKs._Ks_daughterTrack1_dz_min_PV[j],treeKs._Ks_daughterTrack2_dz_min_PV[j]), PV_Z_reweighting_factor)
				histos[17].Fill(maxSigned(treeKs._Ks_daughterTrack1_dz_beamspot[j],treeKs._Ks_daughterTrack2_dz_beamspot[j]), PV_Z_reweighting_factor)
				histos[18].Fill(treeKs._Ks_daughterTrack1_dz_min_PV[j], PV_Z_reweighting_factor)
				histos[19].Fill(treeKs._Ks_daughterTrack2_dz_min_PV[j], PV_Z_reweighting_factor)
				histos[20].Fill(treeKs._Ks_daughterTrack1_dz_min_PV[j], PV_Z_reweighting_factor)
				histos[20].Fill(treeKs._Ks_daughterTrack2_dz_min_PV[j], PV_Z_reweighting_factor)
				histos[21].Fill(treeKs._Ks_daughterTrack1_dz_beamspot[j], PV_Z_reweighting_factor)
				histos[21].Fill(treeKs._Ks_daughterTrack2_dz_beamspot[j], PV_Z_reweighting_factor)
			
				#the 2D plot of lxy verus vz	
				if(iFile == 0):
					h2_RECO_Ks_vz_lxy_data.Fill(treeKs._Ks_vz[j],treeKs._Ks_Lxy[j], PV_Z_reweighting_factor)
					h_Ks_mass_data.Fill(treeKs._Ks_mass[j],PV_Z_reweighting_factor)
				if(iFile == 1):
					h2_RECO_Ks_vz_lxy_mc.Fill(treeKs._Ks_vz[j],treeKs._Ks_Lxy[j], PV_Z_reweighting_factor)
					h_Ks_mass_mc.Fill(treeKs._Ks_mass[j],PV_Z_reweighting_factor)

				#calculate the lifetimes: first calculate lxyz, also in m, in the LAB frame (Lxyz_prime and Lxyz_prime_m) calculate gamma and beta and v, then lifetime in the LAB frame is given by t_prime_Ks and through time dilation you can get it in the COM frame (t_Ks).
				Lxyz_prime = np.sqrt( treeKs._Ks_Lxy[j]*treeKs._Ks_Lxy[j] + (treeKs._Ks_vz[j]-treePV._PV0_vz[0])*(treeKs._Ks_vz[j]-treePV._PV0_vz[0]) )
				Lxyz_prime_m = Lxyz_prime/100
				Ks_p = np.sqrt( treeKs._Ks_pt[j]*treeKs._Ks_pt[j] + treeKs._Ks_pz[j]*treeKs._Ks_pz[j]  )
				gamma_Ks = np.sqrt(1+np.power(Ks_p/massKs,2))
				beta_Ks = np.sqrt(1-1/(gamma_Ks*gamma_Ks))
				v_Ks = beta_Ks*c
				t_prime_Ks = Lxyz_prime_m/v_Ks
				t_Ks = t_prime_Ks/gamma_Ks

				eta_bin = int(abs(treeKs._Ks_eta[j]*10)) - 1
				#eta_bin = 0
				histos_lifetime[ eta_bin ].Fill(t_Ks*c*100)
				#print 'lifetime(cm): ', t_Ks*c*100
				#!!!!!!!!!!!!!The IDEA was to then divide the above lifetime plots by the efficiency versus lxyz from MC, but this is not possible with this MC sample as I do not have the decay position of Ks as Ks are decayed by Geant and Geant decay info is not stored. 
				if(iFile == 1):
					#have to use the reweighing factor only in the nominator of these histograms otherwise will reweigh twice, because they are used for efficiency calculation in the end
					if(treeKs._Ks_deltaRBestMatchingGENParticle[j] < 0.01): histos_Ks_reconstruction_MC_passed[eta_bin].Fill(Lxyz_prime_m,PV_Z_reweighting_factor)
					histos_Ks_reconstruction_MC_total[eta_bin].Fill(Lxyz_prime_m,1)

				if(iFile == 0):
					h2_KsLifetime_vz_lxy_Data.Fill(treeKs._Ks_vz[j],treeKs._Ks_Lxy[j])
				if(iFile == 1):
					h2_KsLifetime_vz_lxy_MC.Fill(treeKs._Ks_vz[j],treeKs._Ks_Lxy[j])

		if(iFile == 0):#for Data
			nKshortPerEventData.append(nKshortThisEvent)
		if(iFile == 1): #for MC
			nKshortPerEventMC.append(nKshortThisEvent)
			
			
	for h in histos:
		h.SetDirectory(0)
	for h in histos_lifetime:
		h.SetDirectory(0)
	if(iFile == 0):
		histos_Data = histos
		histos_lifetimes_Data = histos_lifetime
	else:
		histos_MC = histos
		histos_lifetimes_MC = histos_lifetime

	iFile += 1


fOut.mkdir('PV_Histos')
fOut.cd('PV_Histos')
h_RECO_PV0_vz_Data.Rebin(10)
h_RECO_PV0_vz_Data.Scale(1./h_RECO_PV0_vz_Data.Integral(), "width")
h_RECO_PV0_vz_Data.Write()

h_RECO_PV0_vz_MC.Rebin(10)
h_RECO_PV0_vz_MC.Scale(1./h_RECO_PV0_vz_MC.Integral(), "width")
h_RECO_PV0_vz_MC.Write()

h_RECO_PV0_vz_MC_reweighted_to_data.Scale(1./h_RECO_PV0_vz_MC_reweighted_to_data.Integral(), "width")
h_RECO_PV0_vz_MC_reweighted_to_data.Write()
h_RECO_PV0_vz_MC_NOTreweighted_to_data.Scale(1./h_RECO_PV0_vz_MC_NOTreweighted_to_data.Integral(), "width")
h_RECO_PV0_vz_MC_NOTreweighted_to_data.Write()

tprof_Ks_vz_dz_min_PV_reweighing_factor.Write()
tprof_Ks_vz_reweighing_factor.Write()

h_RECO_beamspot_vz_Data.Rebin(100)
h_RECO_beamspot_vz_Data.Scale(1./h_RECO_beamspot_vz_Data.Integral(), "width")
h_RECO_beamspot_vz_Data.Write()

h_RECO_beamspot_vz_MC.Rebin(100)
h_RECO_beamspot_vz_MC.Scale(1./h_RECO_beamspot_vz_MC.Integral(), "width")
h_RECO_beamspot_vz_MC.Write()


h_RECO_Ks_vz_PV_min_Data.Rebin(100)
h_RECO_Ks_vz_PV_min_Data.Scale(1./h_RECO_Ks_vz_PV_min_Data.Integral(), "width")
h_RECO_Ks_vz_PV_min_Data.Write()

h_RECO_Ks_vz_PV_min_MC.Rebin(100)
h_RECO_Ks_vz_PV_min_MC.Scale(1./h_RECO_Ks_vz_PV_min_MC.Integral(), "width")
h_RECO_Ks_vz_PV_min_MC.Write()

h_RECO_Ks_vz_PV_min_MC_reweighted_to_data.Scale(1./h_RECO_Ks_vz_PV_min_MC_reweighted_to_data.Integral(),"width")
h_RECO_Ks_vz_PV_min_MC_reweighted_to_data.Write()
h_RECO_Ks_vz_PV_min_MC_NOTreweighted_to_data.Scale(1./h_RECO_Ks_vz_PV_min_MC_NOTreweighted_to_data.Integral(),"width")
h_RECO_Ks_vz_PV_min_MC_NOTreweighted_to_data.Write()

#h_RECO_track_multiplicity_Data.Scale(1./h_RECO_track_multiplicity_Data.Integral(), "width")
h_RECO_track_multiplicity_Data.Write()

#h_RECO_track_multiplicity_MC.Scale(1./h_RECO_track_multiplicity_MC.Integral(), "width")
h_RECO_track_multiplicity_MC.Write()

#h2_RECO_Ks_vz_PV_min_track_multiplicity_Data.Scale(1./h2_RECO_Ks_vz_PV_min_track_multiplicity_Data.Integral(),"width")
h2_RECO_Ks_vz_PV_min_track_multiplicity_Data.Write()

#h2_RECO_Ks_vz_PV_min_track_multiplicity_MC.Scale(1./h2_RECO_Ks_vz_PV_min_track_multiplicity_MC.Integral(),"width")
h2_RECO_Ks_vz_PV_min_track_multiplicity_MC.Write()

#save the validation plots for the reweighting procedure to file
validation_reweighing_dir = fOut.mkdir("validation_reweighing")
validation_reweighing_dir.cd()

l_TH1F = [h_RECO_Ks_vz_PV_min_Data,h_RECO_Ks_vz_PV_min_MC,h_RECO_Ks_vz_PV_min_MC_reweighted_to_data]

l_legend_text = ["Data","MC","MC reweighted to data"]

#for h in l_TH1F:
#	h.SetDirectory(0)

i_l_TH1F = 0
legend = TLegend(0.6,0.85,0.99,0.99)
c_validation_reweighing = TCanvas("c_validation_reweighing","");
for h in l_TH1F:
	if(i_l_TH1F==0):
		h.Draw("")
	else:	
		h.Draw("same")
	h.SetLineColor(colours[i_l_TH1F])
	h.SetMarkerStyle(22+i_l_TH1F)
	h.SetMarkerColor(colours[i_l_TH1F])
	h.SetStats(0)
	legend.AddEntry(h,l_legend_text[i_l_TH1F],"lep")
	i_l_TH1F+=1

legend.Draw()
CMS_lumi.CMS_lumi(c_validation_reweighing, 0, 11)
c_validation_reweighing.SaveAs(plots_output_dir+c_validation_reweighing.GetName()+".pdf")
c_validation_reweighing.Write()


#h_RECO_beamspot_vz_MC_reweighted_to_data.Scale(1./h_RECO_beamspot_vz_MC_reweighted_to_data.Integral(), "width")
#h_RECO_beamspot_vz_MC_reweighted_to_data.Write()
#h_RECO_beamspot_vz_MC_NOTreweighted_to_data.Scale(1./h_RECO_beamspot_vz_MC_NOTreweighted_to_data.Integral(), "width")
#h_RECO_beamspot_vz_MC_NOTreweighted_to_data.Write()


fOut.mkdir('Data_Histos')
fOut.cd('Data_Histos')
for h in histos_Data:
#	h.Scale(1./h.Integral(),"Width")
	h.Write()

fOut.mkdir('MC_Histos')
fOut.cd('MC_Histos')
for h in histos_MC:
#	h.Scale(1./h.Integral(),"Width")
	h.Write()

fOut.mkdir('Data_TO_MC_Histos')
fOut.cd('Data_TO_MC_Histos')
i = 0
for i in range(0,len(histos_Data)):
	#histos_Data[i].Divide(histos_MC[i])
	#histos_Data[i].SetName(histos_Data[i].GetName()+'DataToMC')
	#histos_Data[i].Write()

	histo_Divide = histos_Data[i].Clone()
	for j in range(1,histo_Divide.GetNbinsX()):
		data = histos_Data[i].GetBinContent(j) 
		MC = histos_MC[i].GetBinContent(j)
		ratio = 0
		error = 0
		if(MC != 0 and data != 0):
			ratio = data/MC
			error = eff_error(data,MC)
		histo_Divide.SetBinContent(j, ratio)
		histo_Divide.SetBinError(j, error)

	histo_Divide.SetName(histo_Divide.GetName()+'_DataToMC')
	histo_Divide.GetYaxis().SetTitle('Data/MC')
	histo_Divide.Write()

fOut.mkdir('Ks_lifetimes_Data')
fOut.cd('Ks_lifetimes_Data')
for h in histos_lifetimes_Data:
	f_poisson.SetParameters(100, 1)
        h.Fit("f_poisson","R")
	h.Write()
h2_KsLifetime_vz_lxy_Data.Write()

fOut.mkdir('Ks_lifetimes_MC')
fOut.cd('Ks_lifetimes_MC')
for h in histos_lifetimes_MC:
	f_poisson.SetParameters(100, 1)
        h.Fit("f_poisson","R")
	h.Write()
h2_KsLifetime_vz_lxy_MC.Write()

fOut.cd()
#for the 2D histogram of DATA/MC for vz-lxy
h2_KsLifetime_vz_lxy_Data_over_MC = h2_KsLifetime_vz_lxy_Data.Clone()
h2_KsLifetime_vz_lxy_Data_over_MC.SetName("h2_KsLifetime_vz_lxy_Data_over_MC.SetName")
h2_KsLifetime_vz_lxy_Data_over_MC.Divide(h2_KsLifetime_vz_lxy_MC)
h2_KsLifetime_vz_lxy_Data_over_MC.Write()


for i in range(0,len(histos_Ks_reconstruction_MC_passed)):
	h_eff_Ks_reconstruction_MC = histos_Ks_reconstruction_MC_passed[i].Clone()
	h_eff_Ks_reconstruction_MC.Divide(histos_Ks_reconstruction_MC_total[i])
	h_eff_Ks_reconstruction_MC.SetName("h_eff_Ks_reconstruction_MC_"+str(i))
	h_eff_Ks_reconstruction_MC.Write()

	h_corrected_Ks_lifetime_Data = TH1F('h_corrected_Ks_lifetime_Data'+str(i)," ;Ks ct (cm); #entries",lxyz_Ks_bins,0,lxyz_Ks_max)
	for i_bin in range(1,h_corrected_Ks_lifetime_Data.GetNbinsX()):
		uncorrected_data = histos_lifetimes_Data[i].GetBinContent(i_bin) 
		efficiency_MC = h_eff_Ks_reconstruction_MC.GetBinContent(i_bin)
		corrected_data = uncorrected_data
		if efficiency_MC > 0: corrected_data = uncorrected_data/efficiency_MC
		else: print i_bin, "efficiency_MC = 0...."
		h_corrected_Ks_lifetime_Data.SetBinContent(i_bin, corrected_data)
	#f1 = TF1("f1","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, lxyz_Ks_max)
	f_poisson.SetParameters(100, 1)
	h_corrected_Ks_lifetime_Data.Fit("f_poisson","R")
	h_corrected_Ks_lifetime_Data.Write()

h2_RECO_Ks_vz_lxy_data.Write()
h2_RECO_Ks_vz_lxy_mc.Write()


h2_RECO_Ks_vz_xy_ratio = h2_RECO_Ks_vz_lxy_data.Clone()
h2_RECO_Ks_vz_xy_ratio.SetName("h2_RECO_Ks_vz_xy_ratio") 
h2_RECO_Ks_vz_xy_ratio.Divide(h2_RECO_Ks_vz_lxy_mc)
h2_RECO_Ks_vz_xy_ratio.Write()

#masses of the used Ks in data and MC

c_name = "c_Ks_masses"
c = TCanvas(c_name,"")
legend = TLegend(0.8,0.85,0.99,0.99)
l_h_masses = [h_Ks_mass_data,h_Ks_mass_mc]
l_legend = ["data","MC"] 
for j in range(0,len(l_h_masses)):
	h = l_h_masses[j]
	if j == 0:
		h.Draw("PCE1")
	else:
		h.Draw("same")
	h.SetLineColor(colours[j])
	h.SetMarkerStyle(22+j)
	h.SetMarkerColor(colours[j])
	h.SetStats(0)
	legend.AddEntry(h,l_legend[j],"lep")
legend.Draw()
CMS_lumi.CMS_lumi(c, 0, 11)
c.SaveAs(plots_output_dir+c_name+".pdf")
c.Write()

print "number of Kshort going into the histos for Data: ", sum(nKshortPerEventData)
print "number of Kshort going into the histos for MC: ", sum(nKshortPerEventMC)

print "mean number of Kshort per event in data (x1e3): ", sum(nKshortPerEventData)*1e3/len(nKshortPerEventData)
print "mean number of Kshort per event in MC (x1e3): ", sum(nKshortPerEventMC)*1e3/len(nKshortPerEventMC)
print "ratio of the nKshortData/nKshortMC: ", float(sum(nKshortPerEventData))/float(sum(nKshortPerEventMC))


fOut.Close()
