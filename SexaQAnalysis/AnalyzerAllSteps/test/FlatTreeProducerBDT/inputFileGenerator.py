import re
import os

directories= [
#BTagCSV
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/BTagCSV/BTagCSV_Run2016G-07Aug17-v1_trialR/190923_081632/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/BTagCSV/BTagCSV_Run2016H-07Aug17-v1_trialR/190923_081802/0000",
#BTagMu
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/BTagMu/BTagMu_Run2016G-07Aug17-v1_trialR/190906_141052/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/BTagMu/BTagMu_Run2016H-07Aug17-v1_trialR/190906_142015/0000",
#Charmonium
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/Charmonium/Charmonium_Run2016G-07Aug17-v1_trialR/190906_142102/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/Charmonium/Charmonium_Run2016H-07Aug17-v1_trialR/190906_142149/0000",
#DisplacedJet
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/DisplacedJet/DisplacedJet_Run2016G-07Aug17-v1_trialR/190923_084417/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/DisplacedJet/DisplacedJet_Run2016H-07Aug17-v1_trialR/190923_084545/0000",
#DoubleEG
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/DoubleEG/DoubleEG_Run2016G-07Aug17-v1_trialR/190923_084715/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/DoubleEG/DoubleEG_Run2016H-07Aug17-v1_trialR/190923_084842/0000",
#DoubleMuon
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/DoubleMuon/DoubleMuon_Run2016G-07Aug17-v1_trialR/190908_140445/0000",
#DoubleMuonLowMass
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/DoubleMuonLowMass/DoubleMuonLowMass_Run2016G-07Aug17-v1_trialR/190923_082721/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/DoubleMuonLowMass/DoubleMuonLowMass_Run2016H-07Aug17-v1_trialR/190909_071419/0000",
#HTMHT
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/HTMHT/HTMHT_Run2016G-07Aug17-v1_trialR/191007_155452/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/HTMHT/HTMHT_Run2016H-07Aug17-v1_trialR/190910_071209/0000",
#JetHT
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/JetHT/JetHT_Run2016G-07Aug17-v1_trialR/190910_071252/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/JetHT/JetHT_Run2016H-07Aug17-v1_trialR/190910_071343/0000",
#MET
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/MET/MET_Run2016G-07Aug17-v1_trialR/190911_214956/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/MET/MET_Run2016H-07Aug17-v1_trialR/190911_215145/0000",
#MuOnia
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/MuOnia/MuOnia_Run2016G-07Aug17-v2_trialR/190911_215244/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/MuOnia/MuOnia_Run2016H-07Aug17-v1_trialR/190911_215535/0000",
#MuonEG
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/MuonEG/MuonEG_Run2016G-07Aug17-v1_trialR/190912_191107/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/MuonEG/MuonEG_Run2016H-07Aug17-v1_trialR/190917_210911/0000",
#SingleElectron
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SingleElectron/SingleElectron_Run2016G-07Aug17-v1_trialR/190803_145734/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SingleElectron/SingleElectron_Run2016H-07Aug17-v1_trialR/190923_082846/0000",
#SingleMuon
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SingleMuon/SingleMuon_Run2016G-07Aug17-v1_trialR/190730_070712/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SingleMuon/SingleMuon_Run2016H-07Aug17-v1_trialR/190731_114250/0000",
#SinglePhoton
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SinglePhoton/SinglePhoton_Run2016G-07Aug17-v1_trialR/190912_221547/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/SinglePhoton/SinglePhoton_Run2016H-07Aug17-v1_trialR/190912_221721/0000",
#Tau
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/Tau/Tau_Run2016G-07Aug17-v1_trialR/190912_221843/0000",
"/pnfs/iihe/cms/store/user/jdeclerc/data_Sexaq/trialR/Tau/Tau_Run2016H-07Aug17-v1_trialR/190912_222003/0000"
]

for directory in directories:
	directorySplit = re.split(r"/", directory) 
	print directorySplit[10]
	outputFileName = "inputFilesLists/inputFiles_FlatTreeProducerBDT_"+directorySplit[10]+".txt" 
	fOut  = open(outputFileName,"w")
	for filename in os.listdir(directory):
            if "events_skimmed" in filename:
		fOut.write("file:"+directory+"/"+filename + "\n")

	fOut.close()
