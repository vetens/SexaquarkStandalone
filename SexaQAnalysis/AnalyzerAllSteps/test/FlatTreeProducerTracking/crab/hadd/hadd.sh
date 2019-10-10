for i in 01 02 03 04 05 06 07 08 09 10 11 12
do
	hadd combined"$i".root @x"$i"
	#gfal-copy combined"$i".root srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/FlatTree_Skimmed/CRAB_SimSexaq_trial17/crab_FlatTreeProducerTracking_trial17_25092019_v1/190925_194120/combined_FlatTree_Tracking_Skimmed_trial17_part"$i".root
	rm combined"$i".root
done
