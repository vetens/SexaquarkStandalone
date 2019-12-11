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

OutputDir = "Results_10%Unblinging_overlapCheckApplied/" #Results_bkgReference or Results_bkgReference_overlapCheckApplied or Results_partialUnblinging or Results_partialUnblinging_overlapCheckApplied 

colours = [1,2,4,7,8]

#fIn_dir = "BDTApplied_bkgReference_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_overlapCheckApplied/" 
fIn_dir = "BDTApplied/10%Unblind/BDTApplied_10%Unblind_dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing_OverlapCheckTrue/" 
fOut = TFile(OutputDir+'OverLapBDTDistributions.root','RECREATE')


nEventsPDs= [
100014471.,
64209667,
27342838,
28329985,
61398862,
70124218,
18006002,
19433403,
76538905,
82569062,
44182144, #DoubleMuon Run 2016H should come behind this one
44835919,
52461380,
30133219,
33060229,
118152949,
120597407,
25726953,
37033146,
34955119,
37575978,
33573015,
28812785,
15098617,
12093812,
14945745,
17137991,
31884899,
33674739,
78733260,
75409819,
-1,
-1
]


gROOT.SetBatch(kTRUE)

nFilesToRunOver = 0
for fIn in sorted(os.listdir(fIn_dir)):
	if(not fIn.endswith(".root")): continue
	nFilesToRunOver += 1

zoom = False
min_BDT_var = -0.7
max_BDT_var = 0.5
bins_BDT_var = 120
if(zoom):
	min_BDT_var = 0.
	max_BDT_var = 0.4
	bins_BDT_var = 40
h_BDT_ALL_Data = TH1F('h_BDT_ALL_Data','; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)

h_BDT_0 = TH1F('h_BDT_0','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_1 = TH1F('h_BDT_1','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_2 = TH1F('h_BDT_2','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_3 = TH1F('h_BDT_3','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_4 = TH1F('h_BDT_4','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_5 = TH1F('h_BDT_5','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_6 = TH1F('h_BDT_6','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_7 = TH1F('h_BDT_7','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_8 = TH1F('h_BDT_8','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_9 = TH1F('h_BDT_9','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_10 = TH1F('h_BDT_10','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_11 = TH1F('h_BDT_11','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_12 = TH1F('h_BDT_12','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_13 = TH1F('h_BDT_13','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_14 = TH1F('h_BDT_14','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_15 = TH1F('h_BDT_15','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_16 = TH1F('h_BDT_16','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_17 = TH1F('h_BDT_17','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_18 = TH1F('h_BDT_18','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_19 = TH1F('h_BDT_19','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_20 = TH1F('h_BDT_20','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_21 = TH1F('h_BDT_21','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_22 = TH1F('h_BDT_22','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_23 = TH1F('h_BDT_23','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_24 = TH1F('h_BDT_24','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_25 = TH1F('h_BDT_25','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_26 = TH1F('h_BDT_26','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_27 = TH1F('h_BDT_27','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_28 = TH1F('h_BDT_28','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_29 = TH1F('h_BDT_29','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_30 = TH1F('h_BDT_30','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_31 = TH1F('h_BDT_31','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)
h_BDT_32 = TH1F('h_BDT_32','h_BDT; BDT classifier;',bins_BDT_var,min_BDT_var,max_BDT_var)

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

i_File = 0
for fIn in sorted(os.listdir(fIn_dir)):
	if(not fIn.endswith(".root")): continue
#	if(i_File > 2): continue
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
			if(inTree.SexaqBDT>0.1):
				h_BDT_ALL_Data.Fill(inTree.SexaqBDT)

	fOut.cd()
	h2_BDT_mass.Write()
	l_h_BDT_names.append(fIn)
	i_File+=1



#h_BDT_ALL_Data_normalized = h_BDT_ALL_Data.Clone()
#h_BDT_ALL_Data_normalized.SetName(h_BDT_ALL_Data_normalized.GetName()+'_normalized')
#h_BDT_ALL_Data_normalized.Scale(1./h_BDT_ALL_Data_normalized.Integral(),"width")

fOut.cd()

h_BDT_ALL_Data.Write()

#h_BDT_ALL_Data_normalized.Write()

c_BDT_Normalized = TCanvas("c_BDT_Normalized","c_BDT_Normalized",800,800)
c_BDT_Normalized.Divide(6,6)
print 'len l_h_BDT', len(l_h_BDT)
print 'l_h_BDT_names', len(l_h_BDT_names)

lCuts = [0.1,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45]
NEventsAboveCut = [[0.]*nFilesToRunOver for x in xrange(len(lCuts))] 

for i in range(0,len(l_h_BDT_names)):
#	c_BDT_Normalized.cd(i+1)
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

	for j in range(0,len(lCuts)):
		bmin = l_h_BDT[i].GetXaxis().FindBin( lCuts[j] )
		NEventsAboveCut[j][i] = l_h_BDT[i].Integral(bmin,9999)

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


h_total_events_PDs = TH1F("h_total_events_PDs", ";;Total N Events PDs", nFilesToRunOver, -0.5, nFilesToRunOver-0.5)
h_total_events_PDs_scaled = TH1F("h_total_events_PDs_scaled", ";;Total N Events PDs / 1000", nFilesToRunOver, -0.5, nFilesToRunOver-0.5)
h_total_events = TH1F("h_total_events", ";;S surviving pre-BDT cuts", nFilesToRunOver, -0.5, nFilesToRunOver-0.5)
for i in range(0,len(l_h_BDT_names)):
	h_total_events_PDs.SetBinContent(i+1,nEventsPDs[i])
	h_total_events_PDs_scaled.SetBinContent(i+1,nEventsPDs[i]/1000)
	h_total_events.SetBinContent(i+1,l_nevents_total[i])

h_totalEvents_TotalEventsPDsScaled = TH1F("h_totalEvents_TotalEventsPDsScaled", ";;(S surviving pre-BDT cut)/(events)", nFilesToRunOver, -0.5, nFilesToRunOver-0.5)
for i in range(0,len(l_h_BDT_names)):
	h_totalEvents_TotalEventsPDsScaled.SetBinContent(i+1,h_total_events[i+1]/h_total_events_PDs[i+1])

xax0 = h_total_events_PDs.GetXaxis()
xax1 = h_total_events_PDs_scaled.GetXaxis()
xax2 = h_total_events.GetXaxis()
xax3 = h_totalEvents_TotalEventsPDsScaled.GetXaxis()
for i in range(0,len(l_h_BDT_names)):
	xax0.SetBinLabel(i+1,l_h_BDT_names[i][25:-23])
	xax1.SetBinLabel(i+1,l_h_BDT_names[i][25:-23])
	xax2.SetBinLabel(i+1,l_h_BDT_names[i][25:-23])
	xax3.SetBinLabel(i+1,l_h_BDT_names[i][25:-23])
xax0.LabelsOption("v")
xax1.LabelsOption("v")
xax2.LabelsOption("v")
xax3.LabelsOption("v")


for j in range(0,len(lCuts)):
	print lCuts[j], " gives following events surviving: ",NEventsAboveCut[j]
	h_events_surviving = TH1F("h_events_surviving_cutBDT_"+str(lCuts[j]), ";;(S surviving pre-BDT and BDT cut)/(events)", nFilesToRunOver, -0.5, nFilesToRunOver-0.5)
	for i in range(0,len(l_h_BDT_names)):
		h_events_surviving.SetBinContent(i+1,NEventsAboveCut[j][i]/nEventsPDs[i])
	xax = h_events_surviving.GetXaxis()
	for i in range(0,len(l_h_BDT_names)):
		xax.SetBinLabel(i+1,l_h_BDT_names[i][25:-23])
	xax.LabelsOption("v")
	h_events_surviving.Write()

fOut.Write()
fOut.Close()
