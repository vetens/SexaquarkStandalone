from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas

fIn1 = TFile('./AnalyzerAllStepsPlotsStep2WithPU.root', 'read')
fIn2 = TFile('./AnalyzerAllStepsPlotsStep2NoPU.root', 'read')

l_overlap = ["h_NonAntiSTrack_All_pt_cut_eta_Lxy_clone","h_NonAntiSTrack_All_lxy_cut_pt_eta_clone","h_NonAntiSTrack_All_lxy_cut_tight_pt_eta_clone","h_NonAntiSTrack_All_eta_cut_pt_lxy_clone",]

fOut = TFile('OverlapPlots.root','RECREATE')
for eff in l_overlap:
	teff1 = fIn1.Get(eff)
	teff2 = fIn2.Get(eff)
	
	c = TCanvas(eff, eff, 800, 600)
	teff1.Draw()
	teff2.SetLineColor(2)
	teff2.Draw("same")
	c.Write()



fOut.Write()
fOut.Close()


