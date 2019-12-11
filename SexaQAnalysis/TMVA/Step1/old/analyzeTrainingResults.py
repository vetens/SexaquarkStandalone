from ROOT import TFile, TH1F, TH2F, TEfficiency,TCanvas,TLegend,TH1D,TGraph,TMultiGraph
from array import array
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
#fIn = TFile('/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/AnalyzerAllSteps/test/wihtMatchingOnHits/test_TrackMatchingOnHits.root', 'read')
lFilesIn = []

for i in range (0,14):
	lFilesIn.append("LoopCut14Var2016/plots_BDT_v"+str(i)+".root")

print lFilesIn

#c_ROC_overlay = TCanvas("c_ROC_overlay","c_ROC_overlay");
#legend_ROC_overlay = TLegend(0.1,0.1,0.48,0.9);


n_trans_values = 3
trans_value_1 = 1.
trans_value_2 = 0.9995
trans_value_3 = 0.999
trans_values_x = [[0 for i in range(n_trans_values)] for j in range(len(lFilesIn))]
trans_values_y = [[0 for i in range(n_trans_values)] for j in range(len(lFilesIn))]
trans_values_error_y = [[0 for i in range(n_trans_values)] for j in range(len(lFilesIn))]

for i in range(0,len(lFilesIn)):
#for i in range(0,5):
        fIn = TFile(lFilesIn[i], 'read')
        ROC_curve = "LoopCut14Var2016/dataset_BDT_v"+str(i)+"/Method_BDT/BDT/MVA_BDT_rejBvsS"
	print ROC_curve
	h_ROC = TH1D("h_ROC","h_ROC",100,0,1) 
	fIn.GetObject(ROC_curve,h_ROC)

#	h_ROC.SetLineColor(i+1)
#	legend_ROC_overlay.AddEntry(h_ROC,lFilesIn[i],"l")
#	if i == 0:
#		c_ROC_overlay.Draw("")
#	else:
#		c_ROC_overlay.Draw("same")

	#get a value from each ROC curve and plot this on a single graph
	prev_bin_content = 1
	trans1Found = False
	trans2Found = False
	trans3Found = False
	for i_bins in range(1,h_ROC.GetXaxis().GetNbins()):
		bin_content = h_ROC.GetBinContent(i_bins)
		if(bin_content < trans_value_1 and not trans1Found):
			trans1Found = True
			trans_values_y[i][0] = i_bins
			trans_values_error_y[i][0] = h_ROC.GetBinError(i_bins)
			trans_values_x[i][0] = trans_value_1
			
		if(bin_content < trans_value_2 and not trans2Found):
			trans2Found = True
			trans_values_y[i][1] = i_bins
			trans_values_error_y[i][1] = h_ROC.GetBinError(i_bins)
			trans_values_x[i][1] = trans_value_2
		if(bin_content < trans_value_3 and not trans3Found):
			trans3Found = True
			trans_values_y[i][2] = i_bins
			trans_values_error_y[i][2] = h_ROC.GetBinError(i_bins)
			trans_values_x[i][2] = trans_value_3
#fOut = TFile('AnalyzeTraining.root','RECREATE')
#legend_ROC_overlay.Draw()
#c_ROC_overlay.Write()
#fOut.Close()

print "The requested bkg rejection efficiencies are respectively: ",trans_value_1,",",trans_value_2,",",trans_value_3
print "At these bkg rejection efficienes you have a signal efficiency of:"

color=iter(cm.rainbow(np.linspace(0,1,len(trans_values_y))))
for i in range (0,len(trans_values_y)):
	y = trans_values_y[i]
	y_error = trans_values_error_y[i]
	x = trans_values_x[i]
	print i,"\t", y, "\t", y_error
	ax = plt.gca()
   	#plt.plot(x,y,color = next(color),marker='.',linestyle='-', label = lFilesIn[i])
		
#ax.set_xlim([trans_value_1-0.002,trans_value_1+0.0001])
#ax.legend()
#plt.xlabel('Bkg rejection eff')
#plt.ylabel('Signal eff')
#plt.legend(loc=2, prop={'size': 6})
#plt.savefig('mg.pdf')	

