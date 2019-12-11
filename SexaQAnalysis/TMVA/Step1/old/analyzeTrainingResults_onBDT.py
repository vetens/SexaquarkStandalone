from ROOT import TFile, TH1F, TH2F, TEfficiency,TCanvas,TLegend,TH1D,TGraph,TMultiGraph
from array import array
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
#fIn = TFile('/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/AnalyzerAllSteps/test/wihtMatchingOnHits/test_TrackMatchingOnHits.root', 'read')
lFilesIn = []

nVariables = 9

for i in range (0,nVariables):
	lFilesIn.append("LoopCut"+str(nVariables)+"Var2016/plots_BDT_v"+str(i)+".root")

print lFilesIn


for i in range(0,len(lFilesIn)):
#for i in range(0,5):
        fIn = TFile(lFilesIn[i], 'read')
        BDT_curve_BKG = "LoopCut"+str(nVariables)+"Var2016/dataset_BDT_v"+str(i)+"/Method_BDT/BDT/MVA_BDT_B_high"
	#h_BDT_BKG = TH1D("h_BDT","h_BDT",1000,-1,1) 
	#fIn.GetObject(BDT_curve,h_BDT_BKG)

	h_BDT_BKG = fIn.Get(BDT_curve_BKG)

	binmax = (h_BDT_BKG.GetXaxis().FindBin(0.4))
	binwidth = h_BDT_BKG.GetBinWidth(binmax)
	nbins = h_BDT_BKG.GetNbinsX()
	last_bin_with_background = binmax
	i_bin = binmax
	for j in range(1,nbins):
	
		if(h_BDT_BKG.GetBinContent(i_bin) > 0):
			last_bin_with_background = i_bin
			break
		i_bin = binmax-j

	if(last_bin_with_background == binmax) : print "--->>>> ERROR::the binmax value you specified is too small"
	#print "the last bin with background at a BDT variable of ", h_BDT_BKG.GetXaxis().GetBinCenter(last_bin_with_background) ," with content: ", h_BDT_BKG.GetBinContent(last_bin_with_background)
			
        BDT_curve_SGN = "LoopCut"+str(nVariables)+"Var2016/dataset_BDT_v"+str(i)+"/Method_BDT/BDT/MVA_BDT_S_high"
	h_BDT_SGN = fIn.Get(BDT_curve_SGN)
	bminSGN = h_BDT_SGN.GetXaxis().FindBin( h_BDT_BKG.GetXaxis().GetBinCenter(last_bin_with_background)  )
	print i, ")The signal you can keep behind the last bin with BKG is (this should be maximised) : ", h_BDT_SGN.Integral(bminSGN,9999999999999999) / h_BDT_SGN.Integral(-9999999999999999,9999999999999999), '=' , h_BDT_SGN.Integral(bminSGN,9999999999999999) ,"/",h_BDT_SGN.Integral(-9999999999999999,9999999999999999)
