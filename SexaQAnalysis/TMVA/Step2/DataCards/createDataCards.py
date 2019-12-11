from ROOT import *
import numpy as np

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / float(multiplier)

#extract some info from the histograms. The info which you need to extract is from the BDT distribution on MC how much events lie behind a certain cut. And from the BDT distribution of the S in data also the number of events which pass a certain cut
def extractInfo(BDT_cut):
	f_fIn  = TFile.Open('../Results_OverLapBDTDistributions/Results_bkgReference_overlapCheckApplied/OverLapBDTDistributions.root')
        h_BDT_S_Data = f_fIn.Get('h_BDT_ALL_Data')
        h_BDT_S_MC   = f_fIn.Get('h_BDT_32_MC-S-BKG-DYJets')
        h_BDT_AntiS_MC = f_fIn.Get('h_BDT_31_MC-AntiS-BKG-DYJets')
	l_h = [h_BDT_S_Data,h_BDT_S_MC,h_BDT_AntiS_MC]
	l_counts = []
	for h in l_h:
		binx = h.GetXaxis().FindBin(BDT_cut)
		l_counts.append(h.Integral(binx,99999999))
	print 'for BDT cut: ',BDT_cut,' l_counts: ', l_counts
	f_fIn.Close()
	return l_counts

def extractSignalEff(BDT_cut):
	f_fIn  = TFile.Open('../Results_OverLapBDTDistributions/Results_bkgReference/OverLapBDTDistributions.root')	
	
	fIn = TFile('../../Step1/BDTOutput_2016_dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing.root', 'read')
	fIn.cd()
	dir_name = 'dataset_BDT_2016dataset_BDT_2016vSelected19Parameters_CutFiducialRegion_CutDeltaPhi_CutLxy_CutDxyOverLxy_SignalWeighing/Method_BDT/BDT/'
	h_eff = fIn.Get(dir_name+"MVA_BDT_effS")
	binx = h_eff.GetXaxis().FindBin(BDT_cut)
	eff = h_eff.GetBinContent(binx)
	print 'for BDT cut: ', BDT_cut, ' efficiency: ', eff
	return eff
	

#write a datacard
def writeDataCard(BDTCut, signal_events, countsBkg, syst_S_to_AntiS,syst_BDT_sign_eff):
	ext_factor = 0.903
	bkg_rate = countsBkg*ext_factor
        file = open("DataCardsSexaq/datacard_BDTCUT"+str(BDTCut).replace('.','p')+"_signalRate_"+str(signal_events)+".card","w")
	file.write("#"+str(BDTCut)+"\n")
	file.write("#"+str(signal_events)+"\n")
        file.write("imax 1  number of channels"+"\n")
        file.write("jmax 1  number of backgrounds"+"\n")
        file.write("kmax 6  number of nuisance parameters (sources of systematical uncertainties)"+"\n")
	file.write("------------"+"\n")
	file.write("bin bin1"+"\n")
	file.write("observation "+str(bkg_rate)+"\n")
	file.write("------------"+"\n")
	file.write("bin             bin1            bin1"+"\n")
	file.write("process         Sexaq_1p8       BKG"+"\n")
	file.write("process         0               1"+"\n")
	file.write("rate           "+str(signal_events)+ "    "  +str(bkg_rate)+"\n")
	file.write("------------"+"\n")
	file.write("lumi            lnN    1.05    -"+"\n")
	file.write("tracking        lnN    1.24    -"+"\n")
	file.write("syst_BDT_sign_eff lnN    "+str(syst_BDT_sign_eff)+"    -"+"\n")
	file.write("S_to_AntiS      lnN    -  "+str(syst_S_to_AntiS)+"\n")
	file.write("BKG_norm_sys         lnN "+"   -   "+str(1.08)+"\n")
	file.write("BKG_norm_stat        gmN "+str(int(countsBkg))+"   -   "+str(ext_factor)+"\n")

	file.close()

#loop over the configurations you want to change: different BDT cuts and different signal rates
BDT_cut_min = 0
BDT_cut_max = 0.56 #0.56
BDT_cut_step = 0.05
l_BDT_cut = []
i_e=0.
for e in range(0,int(BDT_cut_max/BDT_cut_step)):
	l_BDT_cut.append(BDT_cut_min+i_e*BDT_cut_step)
	i_e+=1	
l_signal_events = [100]

for BDT_cut in l_BDT_cut:
	l_counts = extractInfo(BDT_cut)
	signal_eff = extractSignalEff(BDT_cut)
	#syst_S_to_AntiS = 1.+l_counts[2]/l_counts[1]*np.sqrt(1./l_counts[2]+1./l_counts[1])
	syst_S_to_AntiS = 1.53
	syst_BDT_sign_eff = 0.5*(1-signal_eff)+1.
	for signal_events in l_signal_events:
		writeDataCard(BDT_cut,signal_events*signal_eff,l_counts[0],syst_S_to_AntiS,syst_BDT_sign_eff)

#write the l_BDT_cut and l_signal_events to a file so that these ranges can be used when running combine later to make a 2D plot
conf_file = open("DataCardsSexaq/bin_ranges.dat","w")
a = BDT_cut_min-BDT_cut_step/2.
conf_file.write(str(a)+"\n")
for i in range(0,len(l_BDT_cut)-1):
	a = (l_BDT_cut[i]+l_BDT_cut[i+1])/2.
	conf_file.write(str(a)+"\n")
a = BDT_cut_max+BDT_cut_step/2.
conf_file.write(str(a)+"\n")

conf_file.write("###\n")

conf_file.write(str(0.1)+"\n")
for e in l_signal_events:
	conf_file.write(str(2*e)+"\n")


conf_file.close()
