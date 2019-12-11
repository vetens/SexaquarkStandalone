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
  hbot->GetYaxis()->SetLabelSize(0.04);
  //htop->GetYaxis()->SetLabelSize(ss);
  htop->GetYaxis()->SetLabelSize(.04);
  s = htop->GetYaxis()->GetTitleSize() * scale;
  hbot->GetYaxis()->SetTitleSize(0.04);

   hbot->GetYaxis()->SetTitleOffset(1.2);  
   htop->GetYaxis()->SetTitleOffset(.55);  

  s = htop->GetYaxis()->GetLabelOffset() * scale;

  s = htop->GetXaxis()->GetLabelSize() * scale;
  hbot->GetXaxis()->SetLabelSize(0.04);
  htop->GetXaxis()->SetLabelSize(0.);
  s = htop->GetXaxis()->GetTitleSize() * scale;
  hbot->GetXaxis()->SetTitleSize(0.04);

 

  s = htop->GetXaxis()->GetLabelOffset() * scale;
  hbot->GetXaxis()->SetLabelOffset(s);  


  hbot->GetYaxis()->SetNdivisions(5);
}
void colorIt(TH1D *hMidA, int kCyan ){
	hMidA->SetMarkerColor(kCyan);
	hMidA->SetLineColor(kCyan);
	hMidA->SetMarkerStyle(22);//24
    hMidA->SetMarkerSize(0.6);
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


void loop_macro_exc_inc_dataToMC_BKGReference(){
  TFile *f1   = TFile::Open("../Results_OverLapBDTDistributions/Results_bkgReference/OverLapBDTDistributions.root");
  //for Lambda
  //TFile *f1   = TFile::Open("Results/output_DataToMC_RunG_H_with_dxy_dz_min_PV_cut_reweighing_on_Lambda_vz_Dz_min_PV.root");

   std::vector<const char *> v_hist_name;
   v_hist_name.push_back("Data_MC_plots/h_BDT_DataS_to_MCAntiS.pdf");
   v_hist_name.push_back("Data_MC_plots/h_BDT_DataS_to_MCS.pdf");
   v_hist_name.push_back("Data_MC_plots/h_BDT_MCS_to_MCAntiS.pdf");

   std::vector<const char *> v_hist_data;
   v_hist_data.push_back("h_BDT_ALL_Data");
   v_hist_data.push_back("h_BDT_ALL_Data");
   v_hist_data.push_back("h_BDT_31_MC-AntiS-BKG-DYJets");

   std::vector<const char *> v_hist_mc;
   v_hist_mc.push_back("h_BDT_31_MC-AntiS-BKG-DYJets");
   v_hist_mc.push_back("h_BDT_32_MC-S-BKG-DYJets");
   v_hist_mc.push_back("h_BDT_32_MC-S-BKG-DYJets");
  
   gStyle->SetOptStat(0);
   TFile *fOut = new TFile("Data_MC_plots/Data_MC_ratios_BDTVariable.root","RECREATE");
   //TFile *fOut = new TFile("final_plots_output_DataToMC_RunG_H_with_dxy_dz_min_PV_cut_reweighing_on_Lambda_vz_Dz_min_PV.root","RECREATE");
   gStyle->SetOptStat(0);

   for(int i_h = 0; i_h < v_hist_data.size(); i_h++){

      std::cout << v_hist_data[i_h] << std::endl;

      TH1D *data     =(TH1D*)(f1->Get(v_hist_data[i_h]));
      TH1D *mc     =(TH1D*)(f1->Get(v_hist_mc[i_h]));
      data->Sumw2();
      mc->Sumw2();
	
      //scale the MC with a given factor
      double scale_factorMC = 1.;
      if(i_h == 0 || i_h ==1){ scale_factorMC = (data->GetEntries()/mc->GetEntries());}
      else{data->Rebin(4.);mc->Rebin(4.);}
	
      //for lambda
      //mc->Scale(1.7);

      gStyle->SetOptStat(0);

      TH1D *Mont= (TH1D*)data->Clone("Mont");

      Mont->Divide(data,mc,1.0,scale_factorMC);


      TCanvas *c1 = new TCanvas(v_hist_data[i_h],"", 700, 900);


      c1->cd(); 
      auto legend = new TLegend(0.6,0.8,0.99,0.99);
      colorIt(data,kBlack);  
      splithist(0.3);
      c1->cd(1);
      gPad->SetLogy(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
	
      data->SetMarkerStyle(23);
      data->Draw("peXOC");

      mc->Draw("hhistsames");
      TLegendEntry * entrydata;
      TLegendEntry * entrymc;
      if(i_h == 0 || i_h == 1 )entrydata = legend->AddEntry(data,"data","p");

      if(i_h == 0) entrymc = legend->AddEntry(mc,"MC S","l");
      if(i_h == 1) entrymc = legend->AddEntry(mc,"MC #bar{S}","l");
      if(i_h == 2) entrymc = legend->AddEntry(mc,"MC S","l");
	
      if(i_h == 2)entrydata = legend->AddEntry(data,"MC #bar{S}","p");
	
      entrydata->SetMarkerSize(2.);   
      entrymc->SetMarkerSize(2.);   

      legend->Draw();

    

      c1->cd(2);
      gPad->SetGridy();
      colorIt(Mont,kBlack);
      Mont->SetMinimum(0.5);
      Mont->SetMaximum(1.5);
      Mont->GetXaxis()->SetTitleOffset(1.2);
      Mont->GetYaxis()->SetTitleOffset(3.);
      Mont->SetYTitle("Data/MC scaled to data");
      if(i_h == 2) Mont->SetYTitle("MC #bar{S}/MC S");
      Mont->Draw("peX0C");// rec/gen
      fixsplithist(data,Mont);

      gStyle->SetOptStat(0);
      c1->SetName(v_hist_name[i_h]);
      c1->Write();
      c1->SaveAs(v_hist_name[i_h],".pdf");
  }
  
  gStyle->SetOptStat(0);
  fOut->Close();

}

