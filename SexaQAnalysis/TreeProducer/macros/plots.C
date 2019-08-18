#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>

#include "./classes/HQClass.h"

#include <boost/progress.hpp>
//#include "DataFormats/Math/interface/deltaR.h"



void bookHistos(std::map<string, TH1F*>& histos_1d, std::map<string, TH2F*>& histos_2d){

  histos_1d["track_pt"] = new TH1F("track_pt", "track_pt", 100, 0, 500);
  histos_1d["track_phi"] = new TH1F("track_phi", "track_phi", 200, -4, 4);

  histos_2d["track_phi_eta"] = new TH2F("track_phi_eta", "track_phi_eta", 200, -4, 4, 200, -4, 4);

}



void writeHistos(TFile*& outf, const std::map<string, TH1F*>& histos, const std::map<string, TH2F*>& histos2d){
        outf->cd();
        for (std::pair<std::string, TH1F*> element : histos) {
                std::string word = element.first;
                TH1F* histo = element.second;
                histo->Write();
        }

        for (std::pair<std::string, TH2F*> element2 : histos2d) {
                std::string word = element2.first;
                TH2F* histo = element2.second;
                histo->Write();
        }

}

void plots(){

    std::map<string, TH1F*> histos_1d;
    std::map<string, TH2F*> histos_2d;

    bookHistos(histos_1d, histos_2d);

    TChain *tree = new TChain("tree/HexaQAnalysis");
    tree->Add("../test/MC_tree.root");

    HQClass hqhand;
    hqhand.Init(tree);

    Int_t Nentries = hqhand.fChain->GetEntries();
    std::cout<<"Processing "<<Nentries<<" entries"<<std::endl;
    boost::progress_display show_progress( Nentries );

    for(Int_t entry = 0; entry < Nentries; ++entry){

        ++show_progress;

        hqhand.GetEntry(entry);



        //cout<<hqhand.nTrack<<endl;
        for(int trk=0; trk<hqhand.nTrack; trk++){
            //apply some cuts if needed
            //if(hqhand.track_pt->at(trk)<0.5) continue;

            histos_1d["track_pt"]->Fill(hqhand.track_pt->at(trk));
            histos_1d["track_phi"]->Fill(hqhand.track_phi->at(trk));

            histos_2d["track_phi_eta"]->Fill(hqhand.track_eta->at(trk), hqhand.track_phi->at(trk));

            //TLorentzVector proton;
            //TLorentzVector pi;
            //proton.SetPtEtaPhiM( hqhand.track_pt->at(trk) , hqhand.track_eta->at(trk), hqhand.track_phi->at(trk), XXX);
            //pi.SetPtEtaPhiM( hqhand.track_pt->at(trk) , hqhand.track_eta->at(trk), hqhand.track_phi->at(trk), YYY);
            //TLorentzVector sub = proton - pi;
            //TLorentzVector add = proton + pi;

        }

    }///for entries

    TString outfile = "./HQ_plots.root";
    TFile * fout = new TFile(outfile, "RECREATE");
    writeHistos(fout, histos_1d, histos_2d);
    fout->Close();
}
