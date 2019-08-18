from ROOT import *


fIn = TFile('/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaqWithPU2016NeutrinoGun_tryToFix_8_10072019_v1/combinedAnalyzedSkimmed_SexaqWithPU2016NeutrinoGun_tryToFix_8.root', 'read')
fOut = TFile('accPlots_Skimmed_WithPU2016NeutrinoGun_tryToFix_8.root','RECREATE')

accuracyPlots = []

#Non AntiS tracking efficiency:
parameters = ["pt","phi","eta","mass","lxy_interactionVertex","vz_interactionVertex"]
parameters_text = ["p_{t}(GEN #bar{S})-p_{t}(RECO #bar{S}) (GeV)","#phi(GEN #bar{S})-#phi(RECO #bar{S}) (rad)","#eta(GEN #bar{S})-#eta(RECO #bar{S}) (rad)","mass (GEN #bar{S})-(RECO #bar{S}) (GeV)","l_{0}(GEN #bar{S}) - l_{0}(RECO #bar{S}) (cm)","v_{z}(GEN #bar{S}) - v_{z}(RECO #bar{S}) (cm)"]
fit_range = [0.4,0.12,0.12,1,1,1]
for parameter in parameters:
#for parameter in ["pt","eta","phi","lxy","vz","dxy","pt_cut_eta_Lxy","lxy_cut_pt_eta","lxy_cut_pt_eta_dxy","eta_cut_pt_lxy"]:
	accuracyPlots.append(fIn.Get("AnalyzerGEN/GENRECO/GENRECO_AntiS/GENRECO_AntiS_accuracy/h_RECOAcc_AntiS_"+parameter).Clone())

gStyle.SetOptStat(0)
gStyle.SetOptFit()

color = 36

i = 0
for plots in accuracyPlots:

	print plots.GetName()
	c = TCanvas(plots.GetName(),"")

	plots.Fit("gaus","","",-fit_range[i],fit_range[i])
	print parameters[i]
  	plots.GetXaxis().SetTitle(parameters_text[i])
  	plots.GetXaxis().SetRangeUser(-0.8,0.8)
	plots.SetLineColor(color)
	plots.SetFillColorAlpha(color,0.5)
	plots.Draw("E1L")
#	plots.Write()
	c.Update()	
	c.Write()

	i = i+1

gStyle.SetOptStat(0)
gStyle.SetOptFit()
	
fOut.Write()
fOut.Close()
