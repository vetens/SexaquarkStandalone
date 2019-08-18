from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas


#fIn = TFile('/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/AnalyzerAllSteps/test/wihtMatchingOnHits/test_TrackMatchingOnHits.root', 'read')
fIn = TFile('/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/crmc/Sexaquark_trial10_scalingNAntiS/FlatTree/FlatTreeGEN_Sexaquark_trial10_scalingNAntiS.root', 'read')
tree = fIn.Get('FlatTreeProducerGEN/FlatTreeGENLevel') 

fOut = TFile('AnalyzerFlatTreeProducerGEN_trial10.root','RECREATE')

h_antiS_pt = TH1F('h_antiS_pt','; #bar{S} p_{t} (GeV); #Entries',200,0,20)
h_antiS_pz = TH1F('h_antiS_pz','; #bar{S} |p_{z}| (GeV); #Entries',100,0,100)
h_antiS_eta = TH1F('h_antiS_eta','; #bar{S} #eta ; #Entries',160,-8,8)
h_antiS_eta_cut_eta_smaller4 = TH1F('h_antiS_eta_cut_eta_smaller4','; #bar{S} #eta ; #Entries',160,-8,8)

h_antiS_eta_pt = TH2F('h_antiS_eta_pt',';#bar{S} #eta; #bar{S} p_{t} (GeV); #Entries',160,-8,8,200,0,20)
h_antiS_eta_pz = TH2F('h_antiS_eta_pz',';#bar{S} #eta; #bar{S} |p_{z}| (GeV); #Entries',160,-8,8,200,0,100)

nEntries = tree.GetEntries()
print 'Number of entries in the tree: ', nEntries
for i in range(0,nEntries):
	tree.GetEntry(i)
	h_antiS_pt.Fill(tree._S_pt[0])
	h_antiS_pz.Fill(tree._S_pz[0])
	h_antiS_eta.Fill(tree._S_eta[0])
	h_antiS_eta_pt.Fill(tree._S_eta[0],tree._S_pt[0])
	h_antiS_eta_pz.Fill(tree._S_eta[0],abs(tree._S_pz[0]))
	if abs(tree._S_eta[0]) < 4:
		h_antiS_eta_cut_eta_smaller4.Fill(tree._S_eta[0])

c_antiS_pt = TCanvas("c_antiS_pt","");
h_antiS_pt.DrawNormalized()
c_antiS_pt.Write()

c_antiS_pz = TCanvas("c_antiS_pz","");
h_antiS_pz.DrawNormalized()
c_antiS_pz.Write()

c_antiS_eta = TCanvas("c_antiS_eta","");
h_antiS_eta.DrawNormalized()
c_antiS_eta.Write()

c_antiS_eta_cut_eta_smaller4 = TCanvas("c_antiS_eta_cut_eta_smaller4","");
h_antiS_eta_cut_eta_smaller4.DrawNormalized()
c_antiS_eta_cut_eta_smaller4.Write()


c_antiS_eta_pt = TCanvas("c_antiS_eta_pt","")
h_antiS_eta_pt.DrawNormalized()
c_antiS_eta_pt.Write()

c_antiS_eta_pz = TCanvas("c_antiS_eta_pz","")
h_antiS_eta_pz.DrawNormalized()
c_antiS_eta_pz.Write()

fOut.Write()
fOut.Close()
