import os
import ROOT

class TreeCloner(object):
	

	dirname = "BDT"

	os.mkdir(dirname)

	cwd = os.getcwd()

	iDir = '/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/CRAB_SimSexaq_Skimmed_completely_disabled_cosThetaXYCut_innerHitPosCut_25052019_v1/FlatTree/'

	for fileIn in os.listdir(iDir):
		print fileIn
		fileH  = ROOT.TFile.Open(iDir+fileIn)
		inTree = fileH.Get('FlatTreeProducer/FlatTree')

		fileOut = cwd+'/'+dirname+'/test.root'
		ofile   = ROOT.TFile(fileOut, 'recreate')
		print inTree.GetEntries()
		outTree = inTree.CloneTree()

		ofile.Write()
		ofile.Close()
