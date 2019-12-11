from ROOT import TFile, TH1F, TH2F, TEfficiency 


#fIn = TFile('/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/AnalyzerAllSteps/test/wihtMatchingOnHits/test_TrackMatchingOnHits.root', 'read')
fIn = TFile('/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_9/src/TMVA/Step2/BDT/DiscrApplied_test_FlatTreeProducerData12062019.root', 'read')
fOut = TFile('AnalyzerFlatTreeData.root','RECREATE')


inTree = fIn.Get('FlatTree')

#plot the BDT variable versus the mass
h2_BDT_mass = TH2F('h2_BDT_mass','h2_BDT_mass; S mass (GeV); BDT variable;',200,0,20,200,-1,1)
h_BDT_Normalized = TH1F('h_BDT_Normalized','h_BDT_Normalized; BDT variable',200,-1,1)
for i in range(inTree.GetEntries()):
	inTree.GetEntry(i)
	h2_BDT_mass.Fill(inTree._S_mass[0],inTree.SexaqBDT)
	h_BDT_Normalized.Fill(inTree.SexaqBDT)

h_BDT_Normalized = h_BDT_Normalized.DrawNormalized("ALP")
h_BDT_Normalized.Write()

fOut.Write()
fOut.Close()
