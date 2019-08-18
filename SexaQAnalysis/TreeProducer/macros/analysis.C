#include <iostream>

#include "classes/HQClass.h"

using namespace std;


double covariance(vector<vector<double> >* &cov, int idx, int i, int j){
    int dim = 3;                   //vtx covariance is 3x3
    if(cov->size() == 16) dim = 4; //track covarinance is 4x4
    double covij = cov->at(idx).at(i*dim+j);
    return covij;
}


int analysis(){


    TChain *tree = new TChain("tree/HexaQAnalysis");
    tree->Add("../test/MC_tree.root");

    HQClass hqhand;
    hqhand.Init(tree);


    Long64_t nEntries = tree->GetEntries();
    cout<<nEntries<<endl;

    for (int i=0; i<nEntries; i++){
        hqhand.GetEntry(i);

        //cout<<hqhand.gen_mass->size()<<endl;
        //for(int genp=0; genp<hqhand.gen_mass->size(); genp++){
            //cout<<hqhand.gen_mass->at(genp)<<endl;
        //}

        cout<<hqhand.nTrack<<endl;
        for(int trk=0; trk<hqhand.nTrack; trk++){
            //cout<<hqhand.track_pt->at(trk)<<endl;
            for(int j=0; j<4; j++){
                for(int k=0; k<4; k++){
                    cout<<covariance(hqhand.track_covariance, trk, j, k )<<"\t";
                }
                cout<<endl;
            }
            cout<<endl;
        }


    }

    return 0;

}
