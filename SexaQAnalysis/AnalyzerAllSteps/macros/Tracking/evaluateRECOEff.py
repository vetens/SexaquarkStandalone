import numpy as np
from ROOT import *
import random
import os

inDir = "/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/FlatTree_Skimmed/Skimmed_trial17_1p8GeV_14102019_v1_191014_132642/crab_FlatTreeProducerTracking_trial17_11122019_v1_1p8GeV/191211_081147/0000/"


#fOut = TFile(plots_output_dir+'macro_combined_FlatTree_Tracking_Skimmed_trial17.root','RECREATE')
#fOut = TFile(plots_output_dir+'test.root','RECREATE')

nom = 0.
denom = 0.
for fI in os.listdir(inDir):
	if fI.endswith(".root"):
		print fI
		fIn = TFile("file:"+inDir+fI,'read') 
		tree = fIn.Get('FlatTreeProducerTracking/FlatTreeCounter') 
		for i in range(0,tree.GetEntries()):
			tree.GetEntry(i)
			nom += tree._nRECOAntiS[0]
			denom += tree._nGENAntiS[0]
			print tree._nRECOAntiS[0]/tree._nGENAntiS[0]*100.,"%"
		fIn.Close()

print nom,"/",denom, "=", nom/denom*100. ,"%"
