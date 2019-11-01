#from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas, gROOT
from ROOT import *
import numpy as np
import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

sys.path.insert(1, '/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA')
import configBDT as config
config_dict = config.config_dict

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 0
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

colours = [1,2,4,35,38,41]

MaxEvents = 1e7


plots_output_dir = "Results/"

#have to reweigh on the z location of the PV
n_PVZ = 600
min_PVZ = -30
max_PVZ = 30
h_reweighingFactor_PVz = TH1F('h_reweighingFactor_PVz','; absolute v_{z} PV (cm); Events/mm',n_PVZ,min_PVZ,max_PVZ)
#have to reweigh on the number PVs in each event as well
n_PVn = 60
min_PVn = -0.5
max_PVn = 59.5
h_reweighingFactor_nPV = TH1F('h_reweighingFactor_nPV','; #PV; Events',n_PVn,min_PVn,max_PVn)

h2_reweighingFactor_nPV_PVz = TH2F('h2_reweighingFactor_nPV_PVz','; #PV; absolute v_{z} PV (cm);  Events',n_PVn,min_PVn,max_PVn,n_PVZ,min_PVZ,max_PVZ)

#Get the 2D histograms containing data and MC nPV versus PV_vz
fData = TFile.Open('file:/storage_mnt/storage/user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/src/SexaQAnalysis/PrimaryDataInvestigation/combined_PVDistributionsData_allPrimaryDatasets.root')
h2_nPV_vzPV_Data = fData.Get('h2_nPV_vzPV_Data') 
h2_nPV_Data = fData.Get('h_nPV_Data')
fMC = TFile('file:/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/Tracking/plots_Tracking_AntiS_specific_prePVWeighing_trial17/macro_combined_FlatTree_Tracking_Skimmed_trial17.root')
h2_nPV_vzPV_MC = fMC.Get('PV/h2_nPV_vzPV_MC')

#you first need to scale the data to the number of events in MC
print "mean #PV in data: ", h2_nPV_vzPV_Data.GetMean(1)
print "mean #PV in mc  : ", h2_nPV_vzPV_MC.GetMean(1)
NEventsData = h2_nPV_vzPV_Data.GetEntries()/h2_nPV_vzPV_Data.GetMean(1)
NEventsMC   = h2_nPV_vzPV_MC.GetEntries()/h2_nPV_vzPV_MC.GetMean(1)
h2_nPV_vzPV_Data.Scale(NEventsMC/NEventsData)

print 'Events in the data: ', NEventsData 
print 'Events in the MC:   ', NEventsMC

h_nPV_Data = h2_nPV_vzPV_Data.ProjectionX()
h_nPV_MC = h2_nPV_vzPV_MC.ProjectionX()

h_vzPV_Data = h2_nPV_vzPV_Data.ProjectionY()
h_vzPV_MC = h2_nPV_vzPV_MC.ProjectionY()

if(h2_nPV_vzPV_Data.GetNbinsX() != h2_nPV_vzPV_MC.GetNbinsX()):
	print 'Data and MC plot have different number of bins in X'
if(h2_nPV_vzPV_Data.GetNbinsY() != h2_nPV_vzPV_MC.GetNbinsY()):
	print 'Data and MC plot have different number of bins in Y'



f = open(plots_output_dir+'PUReweighing.txt', "w")

#fill the plots with the reweighing parameter
for i in range(1,h_vzPV_Data.GetNbinsX()+1):
	data_PVvz = h_vzPV_Data.GetBinContent(i)
	mc_PVvz = h_vzPV_MC.GetBinContent(i)
	dataToMC_PVvz = 0.
	if(mc_PVvz > 0.):
		dataToMC_PVvz = data_PVvz/mc_PVvz
	h_reweighingFactor_PVz.SetBinContent(i,dataToMC_PVvz)
	

for i in range(1,h_nPV_Data.GetNbinsX()+1):
	data_nPV = h_nPV_Data.GetBinContent(i)
	mc_nPV = h_nPV_MC.GetBinContent(i)
	dataToMC_nPV = 0.
	if(mc_nPV > 0.):
		dataToMC_nPV = data_nPV/mc_nPV
	h_reweighingFactor_nPV.SetBinContent(i,dataToMC_nPV)	

#f.write(",")
#for j in range(1,h2_nPV_vzPV_Data.GetNbinsY()+1):
#	f.write(str(h2_reweighingFactor_nPV_PVz.GetYaxis().GetBinCenter(j))+",")
#f.write("\n")

#loop over the PU
for i in range(1,h2_nPV_vzPV_Data.GetNbinsX()+1):
	f.write('map<double,double>AnalyzerAllSteps::mapPU'+str(int(h2_reweighingFactor_nPV_PVz.GetXaxis().GetBinCenter(i)))+' = ' + '{' )
	for j in range(1,h2_nPV_vzPV_Data.GetNbinsY()+1):
		data_nPV_vzPV = h2_nPV_vzPV_Data.GetBinContent(i,j)
		mc_nPV_vzPV = h2_nPV_vzPV_MC.GetBinContent(i,j)
		dataToMC_nPV_vzPV = 0.
		if(mc_nPV_vzPV>0.):
			dataToMC_nPV_vzPV = data_nPV_vzPV/mc_nPV_vzPV
		h2_reweighingFactor_nPV_PVz.SetBinContent(i,j,dataToMC_nPV_vzPV)	
		f.write('{'+str(h2_reweighingFactor_nPV_PVz.GetYaxis().GetBinCenter(j)) + ','+str(dataToMC_nPV_vzPV)+"}")
		if(j != h2_nPV_vzPV_Data.GetNbinsY()):
			f.write(',')
	f.write("};\n")
f.write('vector<map<double,double>>AnalyzerAllSteps::v_mapPU{')
for i in range(1,h2_nPV_vzPV_Data.GetNbinsX()+1):
	f.write('mapPU'+str(int(h2_reweighingFactor_nPV_PVz.GetXaxis().GetBinCenter(i))))
	if i != h2_nPV_vzPV_Data.GetNbinsX():
		f.write(',')
f.write('};')
f.close()

#now do the reweighing: loop over the events again and reeigh the mc to data. Do this for the vz distribution of the PV as a test that your reweighing works.


h2_nPV_vzPV_MC_reweighed_2D = h2_nPV_vzPV_MC.Clone()
h2_nPV_vzPV_MC_reweighed_2D.SetName('h2_nPV_vzPV_MC_reweighed_2D')
h2_nPV_vzPV_MC_reweighed_2D.Multiply(h2_reweighingFactor_nPV_PVz)


h_nPV_MC_reweighed_2D  = h2_nPV_vzPV_MC_reweighed_2D.ProjectionX()
h_vzPV_MC_reweighed_2D = h2_nPV_vzPV_MC_reweighed_2D.ProjectionY()


fOut = TFile(plots_output_dir+'PUReweighing.root','RECREATE')

h_reweighingFactor_PVz.Write()

h_vzPV_Data.Write()
h_vzPV_MC.Write()

h_nPV_Data.Write()
h_nPV_MC.Write()

h2_nPV_vzPV_Data.Write()
h2_nPV_vzPV_MC.Write()

h2_reweighingFactor_nPV_PVz.Write()
h_nPV_MC_reweighed_2D.Write()
h_vzPV_MC_reweighed_2D.Write()
h2_nPV_vzPV_MC_reweighed_2D.Write()

TH1_l = [h_vzPV_Data,h_vzPV_MC,h_vzPV_MC_reweighed_2D]
Legend_l = ["Data","MC","MC reweighted"]
Legend_l_type = ["l","l","lep"]
c_name = "c_PVz_2D_PUReweighing"
c = TCanvas(c_name,"")
legend = TLegend(0.8,0.85,0.99,0.99)

for j in [2,1,0]:
	h = TH1_l[j]
	#if(h.GetSumw2N() == 0):
	#	h.Sumw2(kTRUE)
	if j == 0:
		h.Draw("Csame")
	elif j == 1:
		h.Draw("Csame")
	elif j == 2:
		h.Draw("PCE")
	h.SetLineColor(colours[j])
	h.SetLineWidth(2)
	h.SetMarkerStyle(22+j)
	h.SetMarkerColor(colours[j])
	h.SetStats(0)
	legend.AddEntry(h,Legend_l[j],Legend_l_type[j])

legend.Draw()
CMS_lumi.CMS_lumi(c, 0, 11)
c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
c.Write()


c_name = "c_nPV_2D_PUReweighing"
c = TCanvas(c_name,"")
legend = TLegend(0.8,0.85,0.99,0.99)
TH1_l = [h_nPV_Data,h_nPV_MC,h_nPV_MC_reweighed_2D]

for j in [2,1,0]:
	h = TH1_l[j]
	#if(h.GetSumw2N() == 0):
	#	h.Sumw2(kTRUE)
	if j == 0:
		h.Draw("Csame")
	elif j == 1:
		h.Draw("Csame")
	elif j == 2:
		h.Draw("PCE")
	h.SetLineColor(colours[j])
	h.SetLineWidth(2)
	h.SetMarkerStyle(22+j)
	h.SetMarkerColor(colours[j])
	h.SetStats(0)
	legend.AddEntry(h,Legend_l[j],Legend_l_type[j])


legend.Draw()
CMS_lumi.CMS_lumi(c, 0, 11)
c.SaveAs(plots_output_dir+c_name.replace(".", "p")+".pdf")
c.Write()

fOut.Write()
fOut.Close()
