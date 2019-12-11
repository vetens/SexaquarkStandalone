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
  hbot->GetYaxis()->SetLabelSize(0.06);
  //htop->GetYaxis()->SetLabelSize(ss);
  htop->GetYaxis()->SetLabelSize(.04);
  s = htop->GetYaxis()->GetTitleSize() * scale;
  hbot->GetYaxis()->SetTitleSize(0.06);

   hbot->GetYaxis()->SetTitleOffset(0.8);  
   htop->GetYaxis()->SetTitleOffset(1.1);  

  s = htop->GetYaxis()->GetLabelOffset() * scale;

  s = htop->GetXaxis()->GetLabelSize() * scale;
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
    hMidA->SetMarkerSize(0.8);
//	hMidA->SetFillColor(kCyan);
//	hMidA->SetFillStyle(0);	
    hMidA->GetXaxis()->SetLabelSize(0.05);
	hMidA->GetXaxis()->SetTitleSize(0.06);
    hMidA->GetYaxis()->SetTitleSize(0.05);
	hMidA->GetYaxis()->SetTitleOffset(0.55);
	hMidA->GetYaxis()->SetLabelSize(0.05);
	hMidA->SetTitleFont(42, "XYZ");
    hMidA->SetLabelFont(42, "XYZ");
        }


void loop_macro_exc_inc_dataSignalToBKGPartialUnblind(){
   std::vector<const char *> v_hist_name;
   v_hist_name.push_back("../DataSignal_DataBkg_ratios_plots/h_BDT_DataSBkgRef_to_DataAntiS.pdf");


   std::vector<const char *> v_hist_data_bkgRef;
   v_hist_data_bkgRef.push_back("h_BDT_ALL_Data");

 
   std::vector<const char *> v_hist_data_partialUnbl;
   v_hist_data_partialUnbl.push_back("h_BDT_ALL_Data");
  
   gStyle->SetOptStat(0);
   TFile *fOut = new TFile("DataSignal_DataBkg_ratios_plots/DataSignal_DataBkg_ratios.root","RECREATE");
   //TFile *fOut = new TFile("final_plots_output_DataToMC_RunG_H_with_dxy_dz_min_PV_cut_reweighing_on_Lambda_vz_Dz_min_PV.root","RECREATE");
   gStyle->SetOptStat(0);

   for(int i_h = 0; i_h < v_hist_data_bkgRef.size(); i_h++){

      TFile *f1   = TFile::Open("Results_OverLapBDTDistributions/Results_bkgReference_overlapCheckApplied/OverLapBDTDistributions.root");
      TFile *f2   = TFile::Open("Results_OverLapBDTDistributions/Results_partialUnblinging_overlapCheckApplied/OverLapBDTDistributions.root");
      TFile *f3   = TFile::Open("Results_OverLapBDTDistributions/Results_10%Unblinging_overlapCheckApplied/OverLapBDTDistributions.root");

      std::cout << v_hist_data_bkgRef[i_h] << std::endl;

      TH1D *data_bkgRef     =(TH1D*)(f1->Get(v_hist_data_bkgRef[i_h]));
      TH1D *data_partialUnbl     =(TH1D*)(f2->Get(v_hist_data_partialUnbl[i_h]));
      TH1D *data_10perCentUnbl     =(TH1D*)(f3->Get(v_hist_data_partialUnbl[i_h]));
      data_bkgRef->Sumw2();
      data_partialUnbl->Sumw2();
      data_10perCentUnbl->Sumw2();


      fOut->cd();

      data_bkgRef->Rebin(2.);
      data_partialUnbl->Rebin(2.);
      data_10perCentUnbl->Rebin(2.);

	
      //scale the MC with a given factor
      data_10perCentUnbl->Scale(10.); //have to scale with a factor 10, because I only use 10% of the data in the unblinding
      //for lambda
      //mc->Scale(1.7);

      gStyle->SetOptStat(0);

      TH1D *Mont= (TH1D*)data_bkgRef->Clone("Mont");

      Mont->Divide(data_partialUnbl,data_bkgRef,1.0,1.0);


      TCanvas *c1 = new TCanvas(v_hist_data_bkgRef[i_h],"", 700, 900);


      c1->cd(); 
      auto legend = new TLegend(0.65,0.7,0.99,0.99);
      colorIt(data_bkgRef,kRed);  
      splithist(0.3);
      c1->cd(1);
      gPad->SetLogy(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
	
      colorIt(data_10perCentUnbl,kBlue);

      data_partialUnbl->GetYaxis()->SetTitle("N_{ev}/0.2 BDT class.");
      data_10perCentUnbl->GetYaxis()->SetTitle("N_{ev}/0.2 BDT class.");

      data_bkgRef->SetMarkerStyle(23);
      data_bkgRef->Draw("peXOC");
      //data_partialUnbl->Draw("hhist");
      data_partialUnbl->Draw("peXOCsame");
      data_10perCentUnbl->Draw("peXOCsame");
      legend->AddEntry(data_bkgRef,"data S","ple");
      legend->AddEntry(data_partialUnbl,"data #bar{S} BDT class. < 0.1  ","ple");
      legend->AddEntry(data_10perCentUnbl,"data #bar{S} x 10 ","ple");

      data_partialUnbl->SetMarkerSize(0.8);
      legend->Draw();

      data_partialUnbl->GetYaxis()->SetTitle("N_{ev}/0.2 BDT class.");
      data_10perCentUnbl->GetYaxis()->SetTitle("N_{ev}/0.2 BDT class.");


      c1->cd(2);
      gPad->SetGridy();
      colorIt(Mont,kBlack);
      Mont->SetMinimum(0.5);
      Mont->SetMaximum(1.5);
      Mont->GetXaxis()->SetTitleOffset(1.2);
      Mont->GetYaxis()->SetTitleOffset(3.);
      Mont->SetYTitle("data #bar{S}/data S");
      Mont->Fit("pol0","","",-0.1,0.1);
     // Mont->Fit("pol1","","",-0.1,0.1);
     // TLine *line = new TLine(-.1,0.903,.1,0.903);
      Mont->Draw("peX0C");// rec/gen
      //line->Draw("same");
      fixsplithist(data_bkgRef,Mont);

      TH1D *Mont2= (TH1D*)data_bkgRef->Clone("Mont2");
      colorIt(Mont2,kBlack);
      Mont2->Divide(data_10perCentUnbl,data_bkgRef,1.0,1.0);
      Mont2->Draw("peX0Csame");

      gStyle->SetOptStat(0);
      c1->SetName(v_hist_name[i_h]);
      c1->Write();
      c1->SaveAs(v_hist_name[i_h],".pdf");
  }
  
  gStyle->SetOptStat(0);
  fOut->Close();

}

