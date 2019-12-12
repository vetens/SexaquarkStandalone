#script to combine the results from running the DiscrApplication.py script over multiple datasets 

from ROOT import * 
import os

import sys
sys.path.append('/user/jdeclerc/CMSSW_8_0_30_bis/src/SexaQAnalysis/AnalyzerAllSteps/macros/tdrStyle')
import  CMS_lumi, tdrstyle

gROOT.SetBatch(kTRUE)
gStyle.SetLegendTextSize(0.08)

CMS_lumi.writeExtraText = 0
CMS_lumi.extraText = "Simulation"
tdrstyle.setTDRStyle()

#define here a proper name for the OutputDir
OutputDir = "Results_bkgReference_overlapCheckApplied/" #Results_bkgReference or Results_bkgReference_overlapCheckApplied or Results_partialUnblinging or Results_partialUnblinging_overlapCheckApplied 

colours = [1,2,4,7,8]

#define here the directory which has the results from DiscrApplication.py which you want to combine
fIn_dir = "BDTApplied/BkgReference/BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapCheckTrue/"
#fIn_dir = "BDTApplied/10%Unblind/BDTApplied_10%Unblind_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapCheckTrue/" 
fOut = TFile(OutputDir+'OverLapBDTDistributions.root','RECREATE')


gROOT.SetBatch(kTRUE)

nFilesToRunOver = 0
for fIn in sorted(os.listdir(fIn_dir)):
	if(not fIn.endswith(".root")): continue
	nFilesToRunOver += 1

zoom = False #to make the BDT classifier plots zoomed in on the tail of the distribution
min_BDT_var = -0.7
max_BDT_var = 0.5
bins_BDT_var = 120
if(zoom):
	min_BDT_var = 0.
	max_BDT_var = 0.4
	bins_BDT_var = 40

#one histo to combine all the data in the BDT classifier distributions from all Data you are running over
h_BDT_ALL_Data = TH1F('h_BDT_ALL_Data','; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
#one histo for each of the input files (31 data files + some MC BKG files)
h_BDT_0 = TH1F('h_BDT_0','h_BDT_0; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_1 = TH1F('h_BDT_1','h_BDT_1; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_2 = TH1F('h_BDT_2','h_BDT_2; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_3 = TH1F('h_BDT_3','h_BDT_3; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_4 = TH1F('h_BDT_4','h_BDT_4; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_5 = TH1F('h_BDT_5','h_BDT_5; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_6 = TH1F('h_BDT_6','h_BDT_6; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_7 = TH1F('h_BDT_7','h_BDT_7; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_8 = TH1F('h_BDT_8','h_BDT_8; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_9 = TH1F('h_BDT_9','h_BDT_9; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_10 = TH1F('h_BDT_10','h_BDT_10; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_11 = TH1F('h_BDT_11','h_BDT_11; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_12 = TH1F('h_BDT_12','h_BDT_12; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_13 = TH1F('h_BDT_13','h_BDT_13; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_14 = TH1F('h_BDT_14','h_BDT_14; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_15 = TH1F('h_BDT_15','h_BDT_15; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_16 = TH1F('h_BDT_16','h_BDT_16; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_17 = TH1F('h_BDT_17','h_BDT_17; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_18 = TH1F('h_BDT_18','h_BDT_18; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_19 = TH1F('h_BDT_19','h_BDT_19; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_20 = TH1F('h_BDT_20','h_BDT_20; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_21 = TH1F('h_BDT_21','h_BDT_21; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_22 = TH1F('h_BDT_22','h_BDT_22; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_23 = TH1F('h_BDT_23','h_BDT_23; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_24 = TH1F('h_BDT_24','h_BDT_24; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_25 = TH1F('h_BDT_25','h_BDT_25; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_26 = TH1F('h_BDT_26','h_BDT_26; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_27 = TH1F('h_BDT_27','h_BDT_27; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_28 = TH1F('h_BDT_28','h_BDT_28; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_29 = TH1F('h_BDT_29','h_BDT_29; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_30 = TH1F('h_BDT_30','h_BDT_30; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_31 = TH1F('h_BDT_31','h_BDT_31; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_32 = TH1F('h_BDT_32','h_BDT_32; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)

l_h_BDT = [
h_BDT_0, 
h_BDT_1, 
h_BDT_2, 
h_BDT_3, 
h_BDT_4, 
h_BDT_5, 
h_BDT_6, 
h_BDT_7, 
h_BDT_8, 
h_BDT_9, 
h_BDT_10,
h_BDT_11,
h_BDT_12,
h_BDT_13,
h_BDT_14,
h_BDT_15,
h_BDT_16,
h_BDT_17,
h_BDT_18,
h_BDT_19,
h_BDT_20,
h_BDT_21,
h_BDT_22,
h_BDT_23,
h_BDT_24,
h_BDT_25,
h_BDT_26,
h_BDT_27,
h_BDT_28,
h_BDT_29,
h_BDT_30,
h_BDT_31,
h_BDT_32
]

l_h_BDT_names = []
l_nevents_total = [0.]*nFilesToRunOver

fOut.cd()

#run over all the input files and do some stuff with the BDT classifier value
i_File = 0
for fIn in sorted(os.listdir(fIn_dir)):
	if(not fIn.endswith(".root")): continue
	boolData = "Run" in fIn
	print i_File, ",", fIn, ' . Is this a data file? ', boolData
	f_fIn  = TFile.Open(fIn_dir + fIn)
	inTree = f_fIn.Get('FlatTree')

 	l_nevents_total[i_File] = inTree.GetEntries()
	l_h_BDT[i_File].SetName(l_h_BDT[i_File].GetName()+'_'+fIn[25:-23]) 

	#plot the BDT classifier versus the mass
	h2_BDT_mass = TH2F('h2_BDT_mass'+fIn,'h2_BDT_mass; S mass (GeV); BDT classifier;',400,-20,20,200,-1,1)
	i = 0
	for i in range(inTree.GetEntries()):
		inTree.GetEntry(i)
		if(i%1e4 == 0):
			print "reached event: ", i, "/",inTree.GetEntries()
		h2_BDT_mass.Fill(inTree._S_mass[0],inTree.SexaqBDT)
		l_h_BDT[i_File].Fill(inTree.SexaqBDT)
		if(boolData): #the last two files are MC
			#if(inTree.SexaqBDT>0.1):
				h_BDT_ALL_Data.Fill(inTree.SexaqBDT)

	fOut.cd()
	h2_BDT_mass.Write()
	l_h_BDT_names.append(fIn)
	i_File+=1


fOut.cd()

h_BDT_ALL_Data.Write()


c_BDT_Normalized = TCanvas("c_BDT_Normalized","c_BDT_Normalized",800,800)
c_BDT_Normalized.Divide(6,6)
print 'len l_h_BDT', len(l_h_BDT)
print 'l_h_BDT_names', len(l_h_BDT_names)


for i in range(0,len(l_h_BDT_names)):
	p = c_BDT_Normalized.cd(i+1)
	p.cd()
	p.SetLogy()
	l_h_BDT[i].SetTitle(l_h_BDT_names[i])
	l_h_BDT[i].SetLineColor(1)
	l_h_BDT[i].DrawNormalized()
	l = TLine(0.2,0,0.2,1)
        l.SetLineColor(2)
        l.Draw()
	p.Update()
	c_BDT_Normalized.Update()


#c_BDT_Normalized.SaveAs("CanvasDivideBDTDistributionsBkg.pdf")
c_BDT_Normalized.Write()
print "total number of events in the datasets: ", l_nevents_total 


#now overlap all BDT curves on one canvas
c_BDT_all1 = TCanvas("c_BDT_all1","c_BDT_all1",800,800)
c_BDT_all2 = TCanvas("c_BDT_all2","c_BDT_all2",800,800)
c_BDT_all3 = TCanvas("c_BDT_all3","c_BDT_all3",800,800)
c_BDT_all4 = TCanvas("c_BDT_all4","c_BDT_all4",800,800)
c_BDT_all5 = TCanvas("c_BDT_all5","c_BDT_all5",800,800)
c_BDT_all6 = TCanvas("c_BDT_all6","c_BDT_all6",800,800)
c_BDT_all7 = TCanvas("c_BDT_all7","c_BDT_all7",800,800)
l_c_BDT_all = [c_BDT_all1,c_BDT_all2,c_BDT_all3,c_BDT_all4,c_BDT_all5,c_BDT_all6,c_BDT_all7]
legend1 = TLegend(0.6,0.6,1,1)
legend2 = TLegend(0.6,0.6,1,1)
legend3 = TLegend(0.6,0.6,1,1)
legend4 = TLegend(0.6,0.6,1,1)
legend5 = TLegend(0.6,0.6,1,1)
legend6 = TLegend(0.6,0.6,1,1)
legend7 = TLegend(0.6,0.6,1,1)
l_legend = [legend1,legend2,legend3,legend4,legend5,legend6,legend7]
for i in range(0,len(l_h_BDT_names)):
	c_i = int(i/5)
	l_c_BDT_all[c_i].cd()
	l_h_BDT[i].SetLineColor(colours[i-5*c_i])
	l_h_BDT[i].SetTitle("")
	l_h_BDT[i].GetXaxis().SetLabelSize(0.052)
	l_h_BDT[i].GetYaxis().SetLabelSize(0.052)
	l_h_BDT[i].GetXaxis().SetTitleSize(0.052)
	l_h_BDT[i].GetYaxis().SetTitleSize(0.052)
	l_h_BDT[i].SetMaximum(5000)
	if(zoom):
		l_h_BDT[i].SetMaximum(500)
	if(i==0 or i%5 == 0):
		l_h_BDT[i].Draw()
	else:
		l_h_BDT[i].Draw("same")
	l_legend[c_i].AddEntry(l_h_BDT[i],l_h_BDT_names[i][25:-23]+" ","l")

for i in range(0,len(l_c_BDT_all)):	
	l_c_BDT_all[i].cd()
	l_c_BDT_all[i].SetLogy()
	CMS_lumi.CMS_lumi(l_c_BDT_all[i], 0, 11)
	l_legend[i].Draw()	
	l_c_BDT_all[i].Write()

c_BDT_all_combined = TCanvas("c_BDT_all_combined","c_BDT_all_combined",800,800)
c_BDT_all_combined.Divide(2,4)
for i in range(0,len(l_c_BDT_all)):
	c_BDT_all_combined.cd(i+1)
	l_c_BDT_all[i].DrawClonePad()
c_BDT_all_combined.Write()
c_BDT_all_combined.SaveAs(OutputDir+c_BDT_all_combined.GetName()+"_zoom"+str(zoom)+".pdf")



fOut.Write()
fOut.Close()
