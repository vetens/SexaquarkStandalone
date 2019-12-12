#simle script to make a bit nicer looking plots for the BDT efficiency and the BDT overtraining check

from ROOT import *
import sys

sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 0
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

colours = [1,2,4,35,38,41]

plots_output_dir = "plots_BDTEvaluationPlots/" 

fIn = TFile('./BDTOutput_2016_dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing.root', 'read')

fIn.cd()

dir_name = 'dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing/Method_BDT/BDT/'

fOut = TFile(plots_output_dir+'plots_BDTEvaluationPlots.root','RECREATE')

l_overlap = [dir_name+"MVA_BDT_S",dir_name+"MVA_BDT_Train_S",dir_name+"MVA_BDT_B",dir_name+"MVA_BDT_Train_B"]
l_legend  = ["Signal (test sample)","Signal (trainging sample)","Background (test sample)","Background (training sample)"]

c = TCanvas("TrainingAndTesting", "", 800, 600)
i = 0
legend = TLegend(0.68,0.15,0.99,0.35)
for h in l_overlap:
        h1 = fIn.Get(h)
	h1.SetTitle("")
	h1.GetXaxis().SetTitle("BDT classifier/BDT class.")
	h1.GetYaxis().SetTitle("1/N_{ev} Events")
	print h1.GetName()
	if(i==0):
        	h1.Draw()
	else:
		h1.Draw("same")
	h1.SetLineColor(colours[i])
	h1.SetMarkerStyle(22+i)
	h1.SetMarkerColor(colours[i])
	h1.SetStats(0)
	legend.AddEntry(h1,l_legend[i],"lep")
	i+=1
c.SetLogy()
legend.Draw()
CMS_lumi.CMS_lumi(c, 0, 11)
c.SaveAs(plots_output_dir+c.GetName()+".pdf")

c.Write()


l_overlap = [dir_name+"MVA_BDT_effS",dir_name+"MVA_BDT_effB"]
l_legend  = ["Signal efficiency","Background efficiency"]

c = TCanvas("Efficiency", "", 800, 600)
i = 0
legend = TLegend(0.2,0.3,0.5,0.5)
for h in l_overlap:
        h1 = fIn.Get(h)
	h1.SetTitle("")
	h1.GetXaxis().SetTitle("BDT classifier")
	h1.GetYaxis().SetTitle("Efficiency")
	print h1.GetName()
	if(i==0):
        	h1.Draw()
	else:
		h1.Draw("same")
	h1.SetLineColor(colours[i])
	h1.SetMarkerStyle(22+i)
	h1.SetMarkerColor(colours[i])
	h1.SetStats(0)
	legend.AddEntry(h1,l_legend[i],"l")
	i+=1
legend.Draw()
c.SetLogy()
CMS_lumi.CMS_lumi(c, 0, 11)
c.SaveAs(plots_output_dir+c.GetName()+".pdf")

c.Write()

fOut.Write()
fOut.Close()
