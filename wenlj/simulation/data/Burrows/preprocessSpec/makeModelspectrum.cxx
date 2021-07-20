#include <iostream>
#include <vector>
#include <unistd.h>//for getopt()

#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TVectorD.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"

#include "SNBurrowsIntegFcn.hh"

int main(int argc, char* argv[]){

    Int_t imode = -1;

    Int_t opt;
    while((opt = getopt(argc, argv, ":m:"))!= -1){
        switch(opt){
            case'm':
                imode = atoi(optarg);
                break;
            case'?':
                printf("Error: wrong command input option %s\n",optarg);
                break;
            case':':
                printf("Error: Forget to add parameter after option %s\n",optarg);
                break;
        }
    }
    std::cout << "Usage: -m Burrows model number" << std::endl;
    if(imode < 0){
        std::cout << "model number is not assigned!" << std::endl;
        exit(-1);
    }

    //neutrino types
    const Int_t ntype = 3;

    //prepare TGraphs
    Double_t Evmin = 0.;
    Double_t Evmax = 80;
    const Int_t nbin_Ev = 100;
    Double_t step_Ev = (Evmax-Evmin)/nbin_Ev;

    Double_t vecEv[nbin_Ev];
    Double_t spectrum[ntype][nbin_Ev];

    printf("Start processing Burrows model %d\n", imode);
    SNBurrowsIntegFcn sninteg(imode);
    for(Int_t ib=0; ib<nbin_Ev; ib++){
        vecEv[ib] = Evmin+step_Ev*(ib+0.5);
        printf("model %d,Nbin_Ev %d ==> %d bin\n", imode,nbin_Ev, ib);
        for(Int_t it=0; it<ntype; it++){
            spectrum[it][ib] = sninteg.getEnergySpectrum(vecEv[ib],it);
        }
    }
    printf("Finished Burrows model %d\n",imode);

    TGraph* grmodelNu[ntype];
    for(Int_t it=0; it<3; it++){
        grmodelNu[it] = new TGraph(nbin_Ev, vecEv,spectrum[it]);
        grmodelNu[it]->SetName(Form("grmod%dtype%d",imode, it));
    }
    TFile* fout = new TFile(Form("/junofs/users/lihl/neutrino/superN/SNsim/simulation/data/Burrows/energySpeFiles/energySpecBurrows_%d.root",imode),"RECREATE");
    fout->cd();

    for(Int_t it=0; it<ntype; it++){
        grmodelNu[it]->Write();
    }

    fout->Close();



}
