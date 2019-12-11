import os
import glob
import subprocess
from ROOT import *
from array import array
import numpy as np


def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / float(multiplier)

#first read the BDT_cuts and the signal rates which you have used
#conf_file = open("DataCardsSexaq/bin_ranges.dat")
l_BDT_cut = []
l_signal_events = []
bool_signal_events = False
with open("DataCardsSexaq/bin_ranges.dat", 'r') as handle:
    for line in handle:
	if(line.startswith('###')): 
		bool_signal_events = True
		continue

	if(not bool_signal_events):l_BDT_cut.append(float(line))
	else: l_signal_events.append(float(line))


for i in range(0,len(l_BDT_cut)):
	l_BDT_cut[i] = truncate(l_BDT_cut[i],3)

print 'truncated'
print l_BDT_cut
print l_signal_events

#now make a th2F with the above lists
a_l_BDT_cut = array('d',l_BDT_cut)
a_l_signal_events = array('d',l_signal_events)

print ""
print a_l_BDT_cut
print a_l_signal_events

h2_BDTcut_SignalEvents = TH2F('h2_BDTcut_SignalEvents'," ; BDT classifier cut; Signal events; Expected limit on signal strength (r)",len(l_BDT_cut)-1,array('d',l_BDT_cut),len(l_signal_events)-1,array('d',l_signal_events))
h_BDTcut_SignalEvents = TH1F('h_BDTcut'," ; BDT classifier cut; Expected limit on signal strength (r)",len(l_BDT_cut)-1,array('d',l_BDT_cut))
h_BDTcut_minSigma_SignalEvents = TH1F('h_BDTcut_minSigma'," ; BDT classifier cut; Minus 1 \sigma expected limit on signal strength (r)",len(l_BDT_cut)-1,array('d',l_BDT_cut))
h_BDTcut_plusSigma_SignalEvents = TH1F('h_BDTcut_plusSigma'," ; BDT classifier cut; Plus 1 \sigma expected limit on signal strength (r)",len(l_BDT_cut)-1,array('d',l_BDT_cut))

filelist = glob.glob(os.path.join('DataCardsSexaq/', '*.card'))

for infile_name in sorted(filelist):
	infile = open(infile_name)
	BDT_cut = float(infile.readline()[1:])
	signal_rate = float(infile.readline()[1:])
	combineOut = subprocess.check_output(["combine", "-d", infile_name, "-M" ,"AsymptoticLimits","--run","blind"])
	limit = 9999.
	limit_minSigma = 9999.
	limit_plusSigma = 9999.
	for item in combineOut.split("\n"):
		if "Expected 50.0%" in item:
			limit = item.strip()	
		if 'Expected 16.0%' in item:
			limit_minSigma = item.strip()
		if 'Expected 84.0%' in item:
			limit_plusSigma = item.strip()	
	limit_short = float(limit[20:])
	limit_minSigma_short = float(limit_minSigma[20:])
	limit_plusSigma_short = float(limit_plusSigma[20:])
	print "file: ",infile_name, " BDT_cut: ", str(BDT_cut), " signal_rate: ", str(signal_rate), "  --->>> ", limit, " (", limit_short,")", '; min 1*sigma: ', limit_minSigma_short, ' ; plus 1*sigma', limit_plusSigma_short
	h2_BDTcut_SignalEvents.SetBinContent(h2_BDTcut_SignalEvents.GetXaxis().FindBin(BDT_cut),h2_BDTcut_SignalEvents.GetYaxis().FindBin(signal_rate),limit_short)
	h_BDTcut_SignalEvents.SetBinContent(h_BDTcut_SignalEvents.GetXaxis().FindBin(BDT_cut),limit_short)
	h_BDTcut_minSigma_SignalEvents.SetBinContent(h_BDTcut_minSigma_SignalEvents.GetXaxis().FindBin(BDT_cut),limit_minSigma_short)
	h_BDTcut_plusSigma_SignalEvents.SetBinContent(h_BDTcut_plusSigma_SignalEvents.GetXaxis().FindBin(BDT_cut),limit_plusSigma_short)

fOut = TFile.Open("Results/combineSexaqResults.root","RECREATE")
h2_BDTcut_SignalEvents.Write()
h_BDTcut_SignalEvents.Write()
h_BDTcut_minSigma_SignalEvents.Write()
h_BDTcut_plusSigma_SignalEvents.Write()
fOut.Close()
