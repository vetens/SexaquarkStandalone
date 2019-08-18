from ROOT import TFile, TH1F, TH2F, TEfficiency, TH1D, TH2D, TCanvas
import numpy as np

fIn1 = TFile('./analyzedFlatTreeV0s_test_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root', 'read') #MC: needs to be plots containing Ks and Lambdas which are not from the AntiS 
fIn2 = TFile('./analyzedFlatTreeV0s_t_SingleMuon_FlatTreeV0s.root', 'read') #Data: needs to be plots from ZeroBias sample with V0s reconstructed according to the new algorithm

#MC_dir = "AnalyzerGEN/GENRECO/GENRECO_KsNonAntiS/GENRECO_RECO_KsNonAntiS/RECO_Ks_extra_mass_cut/"
MC_dir = "AnalyzerRECO/RECO/RECO_Ks/RECO_Ks_extra_mass_cut/RECO_Ks_extra_mass_cut_kinregA/"
Data_dir = "AnalyzerRECO/RECO/RECO_Ks/RECO_Ks_extra_mass_cut/RECO_Ks_extra_mass_cut_kinregA/"

#l_overlap_MC = [
#MC_dir+"h_GENRECO_RECO_RECO_Ks_dxy",
#MC_dir+"h_GENRECO_RECO_RECO_Ks_dz",
#MC_dir+"h2_GENRECO_RECO_RECO_Ks_lxy_dz",
#MC_dir+"h2_GENRECO_RECO_RECO_Ks_lxy_vz"
#]

l_overlap = [
"h2_RECO_Ks_lxy_vz_kinregA",
"h2_RECO_Ks_pt_pz_kinregA",
"h2_RECO_Ks_dxy_dz_min_kinregA",
"h2_RECO_Ks_eta_phi_kinregA"
]


fOut = TFile('DividePlots.root','RECREATE')


def rebin2(h, ngx, ngy):

   #Rebin 2-d histogram h, grouping ngx bins together along X
   #and ngy bins together along Y
   #NB: this macro ignores histogram errors if defined
   
   #make a clone of h
   hold = h.Clone();
   hold.SetDirectory(0);

   nbinsx = hold.GetXaxis().GetNbins();
   nbinsy = hold.GetYaxis().GetNbins();
   xmin  = hold.GetXaxis().GetXmin();
   xmax  = hold.GetXaxis().GetXmax();
   ymin  = hold.GetYaxis().GetXmin();
   ymax  = hold.GetYaxis().GetXmax();
   nx = nbinsx/ngx;
   ny = nbinsy/ngy;
   h.SetBins (nx,xmin,xmax,ny,ymin,ymax);

   #loop on all bins to reset contents and errors
   cu = 0.;
   bx = 0.;
   by = 0.;
   ix = 0;
   iy = 0;
   ibin = 0;
   Bin = 0;
   binx = 0;
   biny = 0;
   for biny in range(1,nbinsy):
      for binx in range(1,nbinsx):
         ibin = h.GetBin(binx,biny);
         h.SetBinContent(ibin,0);
   
   #loop on all bins and refill
   for biny in range(1,nbinsy):
      by  = hold.GetYaxis().GetBinCenter(biny);
      iy  = h.GetYaxis().FindBin(by);
      for binx in range(1,nbinsx):
         bx  = hold.GetXaxis().GetBinCenter(binx);
         ix  = h.GetXaxis().FindBin(bx);
         Bin = hold.GetBin(binx,biny);
         ibin= h.GetBin(ix,iy);
         cu  = hold.GetBinContent(Bin);
         h.AddBinContent(ibin,cu);
   
   return h;



i = 0
for h in l_overlap:
	print i
	fOut.mkdir(str(i))
	fOut.cd(str(i))
	h1 = fIn1.Get(l_overlap[i]) #MC
	h2 = fIn2.Get(l_overlap[i]) #Data

	print l_overlap[i]
	print "MC number of entries in 2D histo: ", h1.GetEntries()
	print "Data number of entries in 2D histo: ", h2.GetEntries()

	for p in ['X','Y']:	
		
		h1_px = h1.ProjectionX()
		if(p =='Y'):
			h1_px = h1.ProjectionY()
		h1_px.Rebin(1)
		h1_px.SetName(p+"_proj_MC")
			

		h2_px = h2.ProjectionX()
		if(p == 'Y'):
			h2_px = h2.ProjectionY()
		h2_px.Rebin(1)
		h2_px.SetName(p+"_proj_Data")
		
	
		error_px = []
		for j in range(1,h1_px.GetNbinsX()):
			if(h1_px.GetBinContent(j) > 0 and h2_px.GetBinContent(j) > 0):
				error_px.append(h2_px.GetBinContent(j)/h1_px.GetBinContent(j)*np.sqrt(1/h1_px.GetBinContent(j)+1/h2_px.GetBinContent(j)))
			else:
				error_px.append(0)
		
		h_px_eff = h2_px.Clone()
		h_px_eff.Divide(h1_px)
		h_px_eff.SetName("p"+p+"_eff")
		for j in range(1,h_px_eff.GetNbinsX()):
			h_px_eff.SetBinError(j,error_px[j-1])
		h_px_eff.Write()

	#rebin the 2D histos and write them to file
	h1_rebin = rebin2(h1,15,15)
	h1_rebin.SetName("2D_MC_rebin")
	h1_rebin.Write()
	h2_rebin = rebin2(h2,15,15)
	h2_rebin.SetName("2D_Data_rebin")
        h2_rebin.Write()

	#Divide Data over MC and write the result to file
	h2_rebin_clone = h2_rebin.Clone()
	h2_rebin_clone.SetName("2D_Data_over_MC")
	h2_rebin_clone.Divide(h1_rebin)
	h2_rebin_clone.Write()

	
	#calculate the errors on the Data/MC
	twoD_Data_over_MC_error = TH2F("twoD_Data_over_MC_error","twoD_Data_over_MC_error",h1_rebin.GetNbinsX(),h1_rebin.GetXaxis().GetXmin(),h1_rebin.GetXaxis().GetXmax(),h1_rebin.GetNbinsY(),h1_rebin.GetYaxis().GetXmin(),h1_rebin.GetYaxis().GetXmax())
	for r in range(1,h1_rebin.GetNbinsX()):
		for q in range(1,h1_rebin.GetNbinsY()):	
			if(h1_rebin.GetBinContent(r,q) > 0 and h2_rebin.GetBinContent(r,q) > 0):
				error = (h2_rebin.GetBinContent(r,q)/h1_rebin.GetBinContent(r,q)*np.sqrt(1/h1_rebin.GetBinContent(r,q)+1/h2_rebin.GetBinContent(r,q)))
				twoD_Data_over_MC_error.SetBinContent(r,q,error) 
			else:
				twoD_Data_over_MC_error.SetBinContent(r,q,0)

	twoD_Data_over_MC_error.Write()

	i += 1

fOut.Write()
fOut.Close()


