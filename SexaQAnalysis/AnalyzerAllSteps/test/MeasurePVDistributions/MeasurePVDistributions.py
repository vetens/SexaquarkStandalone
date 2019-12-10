#script to get the z distribution of PVs and the number of PV distributions (PU). You need these specific distributions for MC reweighing to data.

from ROOT import *

inputFilesList_name = '../FlatTreeProducerBDT/inputFilesLists/inputFiles_FlatTreeProducerBDT_SinglePhoton_Run2016H-07Aug17-v1_trialR.txt'
inputFile = open(inputFilesList_name,"r")

with open(inputFilesList_name) as fp:
	EDM_file = fp.readline()
	while EDM_file:
		print 'running on file: ', EDM_file[5:]
		f = TFile.Open('/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SinglePhoton/SinglePhoton_Run2016H-07Aug17-v1_trialR/190912_221721/0000/events_skimmed_19.root')
		f.cd()
		eventsTree = f.Get('Events') 
		print 'nentries: ', eventsTree.GetEntries()
		for event in eventsTree:
			#f.cd('Events/recoVertexs_offlinePrimaryVertices__RECO/obj/position_/fCoordinates')
			#pvTree = 
			print event.recoVertexs_offlinePrimaryVertices__RECO.obj.position_.fCoordinates.fZ

		EDM_file = fp.readline()
		
