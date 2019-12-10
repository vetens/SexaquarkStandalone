#include <iostream>
#include <fstream>

void splithist(double ratio=0.2){
  TVirtualPad* pmain = TVirtualPad::Pad();
  
  double h = 1. - pmain->GetTopMargin() - pmain->GetBottomMargin();

  double xlow = 0.; //gStyle->GetPadLeftMargin();
  double xhigh = 1.; // - gStyle->GetPadRightMargin();

  double ytop = 1.; //- gStyle->GetPadTopMargin();
  double ybot = 0.; //gStyle->GetPadBottomMargin();
  double ymid = pmain->GetBottomMargin() + ratio * h; //ybot + ratio * (ytop-ybot);
  

  double yp1bot = ymid;
  double yp2top = ymid;

  TPad* p1 = new TPad(TString(pmain->GetName()) + "_1", pmain->GetTitle(),
		      xlow, yp1bot, xhigh, ytop);
  p1->SetNumber(1);
  TPad* p2 = new TPad(TString(pmain->GetName()) + "_2", pmain->GetTitle(), 
		      xlow, ybot, xhigh, yp2top);
  p2->SetNumber(2);
  p1->SetFillStyle(4000);
  p2->SetFillStyle(4000);

  double p1h = ytop - yp1bot;
  double p2h = yp2top - ybot;
  

/*
  p1->SetTopMargin(pmain->GetTopMargin()/p1h);
  p1->SetBottomMargin((ymid-yp1bot)/p1h);
  p1->SetLeftMargin(pmain->GetLeftMargin());
  p1->SetRightMargin(pmain->GetRightMargin());
  
  p2->SetTopMargin((ymid-yp2top)/p2h);
  p2->SetBottomMargin(pmain->GetBottomMargin()/p2h);
  p2->SetLeftMargin(pmain->GetLeftMargin());
  p2->SetRightMargin(pmain->GetRightMargin());
*/
    p1->SetTopMargin(0.11);
    p1->SetBottomMargin(0.);
    p1->SetRightMargin(0.02);
 
    p2->SetTopMargin(0.);
    p2->SetBottomMargin(0.3);
    p2->SetRightMargin(0.02);



  p2->Draw();
  p1->Draw();
  pmain->Modified();

  p1->cd();
}

void fixsplithist(TH1* htop, TH1* hbot){
  TVirtualPad* pmain = TVirtualPad::Pad()->GetCanvas();
  if(!pmain) return;
  TVirtualPad* p1 = pmain->cd(1);
  TVirtualPad* p2 = pmain->cd(2);

  if(!p1 || !p2) return;
  
  double scale = p1->GetHNDC() / p2->GetHNDC();

  double s = htop->GetYaxis()->GetLabelSize() * scale;
  double ss = htop->GetYaxis()->GetLabelSize();
  //hbot->GetYaxis()->SetLabelSize(s);
  hbot->GetYaxis()->SetLabelSize(0.05);
  //htop->GetYaxis()->SetLabelSize(ss);
  htop->GetYaxis()->SetLabelSize(.6);
  htop->GetYaxis()->SetTitleSize(.6);
  //s = htop->GetYaxis()->GetTitleSize() * scale;
  hbot->GetYaxis()->SetTitleSize(0.05);

   hbot->GetYaxis()->SetTitleOffset(.3);  
   htop->GetYaxis()->SetTitleOffset(.2);  

  s = htop->GetYaxis()->GetLabelOffset() * scale;

  s = htop->GetXaxis()->GetLabelSize() ;
  hbot->GetXaxis()->SetLabelSize(0.06);
  htop->GetXaxis()->SetLabelSize(0.);
  s = htop->GetXaxis()->GetTitleSize() * scale;
  hbot->GetXaxis()->SetTitleSize(0.06);

 

  s = htop->GetXaxis()->GetLabelOffset() * scale;
  hbot->GetXaxis()->SetLabelOffset(s);  


  hbot->GetYaxis()->SetNdivisions(5);
}
void colorIt(TH1D *hMidA, int kCyan ){
	hMidA->SetMarkerColor(kCyan);
	hMidA->SetLineColor(kCyan);
	hMidA->SetMarkerStyle(22);//24
    hMidA->SetMarkerSize(0.9);
//	hMidA->SetFillColor(kCyan);
//	hMidA->SetFillStyle(0);	
    hMidA->GetXaxis()->SetLabelSize(0.05);
	hMidA->GetXaxis()->SetTitleSize(0.1);
    hMidA->GetYaxis()->SetTitleSize(0.1);
	hMidA->GetYaxis()->SetTitleOffset(0.55);
	hMidA->GetYaxis()->SetLabelSize(0.12);
	hMidA->SetTitleFont(42, "XYZ");
    hMidA->SetLabelFont(42, "XYZ");
        }


void loop_macro_exc_inc(){
  TFile *f1   = TFile::Open("Results_test/output_DataToMC_RunG_H_with_dxy_dz_min_PV_reweighing_on_Ks_vz_Dz_min_PV.root");
  //for Lambda
  //TFile *f1   = TFile::Open("Results/output_DataToMC_RunG_H_with_dxy_dz_min_PV_cut_reweighing_on_Lambda_vz_Dz_min_PV.root");

   std::vector<const char *> v_file_name;
   v_file_name.push_back("Data_MC_plots/h_RECO_Ks_vz.dat");
   v_file_name.push_back("Data_MC_plots/h_RECO_Ks_lxy.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_pt.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_pt_tracks1.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_pt_tracks2.dat");
   v_file_name.push_back("Data_MC_plots/h_RECO_Ks_pt_tracks1_and_2.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_pz.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_pz_tracks1.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_pz_tracks2.dat");
   v_file_name.push_back("Data_MC_plots/h_RECO_Ks_pz_tracks1_and_2.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_dxy_PV.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_dz.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track1Track2_openingsAngle_ptCut1p2.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track1Track2_max_dxy_beamspot.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track1_dxy_beamspot.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track2_dxy_beamspot.dat");
   v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track1_and_2_dxy_beamspot.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track1Track2_max_dz_min_PV.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track1_dz_min_PV.dat");
   //v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track2_dz_min_PV.dat");
   v_file_name.push_back("Data_MC_plots/h_RECO_Ks_Track1_and_2_dz_min_PV.dat");
   
   std::vector<const char *> v_hist_name;
   v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_vz.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_lxy.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_pt.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_pt_tracks1.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_pt_tracks2.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_pt_tracks1_and_2.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_pz.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_pz_tracks1.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_pz_tracks2.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_pz_tracks1_and_2.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_dxy_PV.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_dz.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track1Track2_openingsAngle_ptCut1p2.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track1Track2_max_dxy_beamspot.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track1_dxy_beamspot.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track2_dxy_beamspot.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track1_and_2_dxy_beamspot.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track1Track2_max_dz_min_PV.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track1_dz_min_PV.pdf");
   //v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track2_dz_min_PV.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Ks_Track1_and_2_dz_min_PV.pdf");

   std::vector<const char *> v_hist_data;
   v_hist_data.push_back("Data_Histos/h_RECO_Ks_vz");
   v_hist_data.push_back("Data_Histos/h_RECO_Ks_lxy");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_pt");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_pt_tracks1");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_pt_tracks2");
   v_hist_data.push_back("Data_Histos/h_RECO_Ks_pt_tracks1_and_2");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_pz");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_pz_tracks1");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_pz_tracks2");
   v_hist_data.push_back("Data_Histos/h_RECO_Ks_pz_tracks1_and_2");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_dxy_PV");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_dz");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track1Track2_openingsAngle_ptCut1p2");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track1Track2_max_dxy_beamspot");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track1_dxy_beamspot");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track2_dxy_beamspot");
   v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track1_and_2_dxy_beamspot");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track1Track2_max_dz_min_PV");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track1_dz_min_PV");
   //v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track2_dz_min_PV");
   v_hist_data.push_back("Data_Histos/h_RECO_Ks_Track1_and_2_dz_min_PV");

   std::vector<const char *> v_hist_mc;
   v_hist_mc.push_back("MC_Histos/h_RECO_Ks_vz");
   v_hist_mc.push_back("MC_Histos/h_RECO_Ks_lxy");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_pt");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_pt_tracks1");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_pt_tracks2");
   v_hist_mc.push_back("MC_Histos/h_RECO_Ks_pt_tracks1_and_2");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_pz");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_pz_tracks1");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_pz_tracks2");
   v_hist_mc.push_back("MC_Histos/h_RECO_Ks_pz_tracks1_and_2");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_dxy_PV");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_dz");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track1Track2_openingsAngle_ptCut1p2");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track1Track2_max_dxy_beamspot");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track1_dxy_beamspot");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track2_dxy_beamspot");
   v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track1_and_2_dxy_beamspot");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track1Track2_max_dz_min_PV");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track1_dz_min_PV");
   //v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track2_dz_min_PV");
   v_hist_mc.push_back("MC_Histos/h_RECO_Ks_Track1_and_2_dz_min_PV");
  
	
   //for Lambda	
  /* std::vector<const char *> v_hist_name;
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_vz.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_lxy.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_pt.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_pt_tracks1.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_pt_tracks2.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_pz.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_pz_tracks1.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_pz_tracks2.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_dxy_PV.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_dz.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p2.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_Track1Track2_max_dxy_beamspot.pdf");
   v_hist_name.push_back("Data_MC_plots/h_RECO_Lambda_Track1Track2_max_dz_min_PV.pdf");

   std::vector<const char *> v_hist_data;
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_vz");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_lxy");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_pt");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_pt_tracks1");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_pt_tracks2");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_pz");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_pz_tracks1");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_pz_tracks2");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_dxy_PV");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_dz");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p2");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_Track1Track2_max_dxy_beamspot");
   v_hist_data.push_back("Data_Histos/h_RECO_Lambda_Track1Track2_max_dz_min_PV");

   std::vector<const char *> v_hist_mc;
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_vz");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_lxy");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_pt");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_pt_tracks1");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_pt_tracks2");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_pz");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_pz_tracks1");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_pz_tracks2");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_dxy_PV");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_dz");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_Track1Track2_openingsAngle_ptCut1p2");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_Track1Track2_max_dxy_beamspot");
   v_hist_mc.push_back("MC_Histos/h_RECO_Lambda_Track1Track2_max_dz_min_PV");
*/
   gStyle->SetOptStat(0);
   TFile *fOut = new TFile("final_plots_output_DataToMC_RunG_with_dxy_dz_min_PV_cut_reweighing_on_Ks_vz_Dz_min_PV.root","RECREATE");
   //TFile *fOut = new TFile("final_plots_output_DataToMC_RunG_H_with_dxy_dz_min_PV_cut_reweighing_on_Lambda_vz_Dz_min_PV.root","RECREATE");
   gStyle->SetOptStat(0);

   for(int i_h = 0; i_h < v_hist_data.size(); i_h++){

      std::cout << v_hist_data[i_h] << std::endl;

      TH1D *data     =(TH1D*)(f1->Get(v_hist_data[i_h]));
      TH1D *mc     =(TH1D*)(f1->Get(v_hist_mc[i_h]));
      data->Sumw2();
      mc->Sumw2();
	
      //scale the MC with a given factor
      //mc->Scale(1.20363184668);
      //for lambda
      //mc->Scale(1.7);

      gStyle->SetOptStat(0);

      TH1D *Mont= (TH1D*)data->Clone("Mont");

      Mont->Divide(data,mc,1.0,1.0);


      TCanvas *c1 = new TCanvas(v_hist_data[i_h],"", 700, 900);


      c1->cd(); 
      auto legend = new TLegend(0.75,0.75,0.9,0.9);
      colorIt(data,kBlue);  
      splithist(0.5);
      c1->cd(1);
      gPad->SetLogy(); 
      gPad->SetGridx(); 
      gPad->SetGridy();
      mc->GetYaxis()->SetLabelSize(0.05);
      data->GetYaxis()->SetLabelSize(0.05);
      mc->GetYaxis()->SetTitleOffset(1.);
      data->GetYaxis()->SetTitleOffset(1.);
      mc->GetYaxis()->SetTitleSize(0.05);
      data->GetYaxis()->SetTitleSize(0.05);
      mc->Draw("hhist");
      legend->AddEntry(mc,"MC","l");
      //for lambda
      //legend->AddEntry(mc,"1.72*MC","l");
      data->SetMarkerStyle(23);
      data->Draw("peXOCsames");
      legend->AddEntry(data,"data","p");
      legend->Draw();

      c1->cd(2);
      gPad->SetGridy();
      colorIt(Mont,kBlack);
      Mont->SetMinimum(0.5);
      Mont->SetMaximum(1.99);
      Mont->GetXaxis()->SetTitleOffset(1.2);
      Mont->GetYaxis()->SetTitleOffset(1.2);
      Mont->SetYTitle("Data/MC");
      Mont->Draw("peX0C");// rec/gen
      fixsplithist(data,Mont);

      //loop over all the bins and write the DATA/MC ratio to a file
       ofstream myfile;
       myfile.open(v_file_name[i_h]);
       myfile <<  "kinVariable" << "," << "DataOverMC" << "," << "DataOverMCError" << "," << "FiducialRegionRequirement" << "\n";     
     for(unsigned int i_bin = 1; i_bin <= Mont->GetXaxis()->GetNbins(); i_bin++){
	double kinVariable = Mont->GetXaxis()->GetBinCenter(i_bin);
	double DataOverMC = Mont->GetBinContent(i_bin);
	double DataOverMCError = Mont->GetBinError(i_bin);
	myfile << kinVariable << "," << DataOverMC << "," << DataOverMCError << "," << ( ( abs(DataOverMC-1) < 0.2 ) ) << "\n";
      }
      myfile.close();
      gStyle->SetOptStat(0);
      c1->Write();
      c1->SaveAs(v_hist_name[i_h],".pdf");
  }
  
  gStyle->SetOptStat(0);
  fOut->Close();

}

