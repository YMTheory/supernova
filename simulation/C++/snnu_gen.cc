#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

#include "SNdetect.hh"
#include "SNsource.hh"
#include "SNnumGarchingSrc.hh"
#include "SNGarchingIntegFcn.hh"
#include "SNeffectLS.hh"

int main(int argc, char* argv[]) {
    Int_t imod = -1;//SN neutrino source model
    Int_t icha = -1;//nu interaction channel
    Double_t nuMass = 0.0;// nu mass
    Int_t MH = 0 ; // mass ordering
    Double_t Ethr = 0.2 ; // energy threshold MeV
    Double_t tmin = 0.0;
    Double_t tmax = 0.0;
    Double_t scale = 1;
    Int_t    dist  = 10;
    //==========getopt()===========//
    Int_t opt;
    while((opt = getopt(argc, argv, ":e:m:n:c:t:o:x:y:d:s:"))!= -1){
        switch(opt){
            case'm':
                imod = atoi(optarg);
                break;
            case'c':
                icha = atoi(optarg);
                break;
            case 'n':
                nuMass = atof(optarg);
                break;
            case 'o':
                MH = atoi(optarg);
            case 't':
                Ethr = atof(optarg);
            case 'x':
                tmin = atof(optarg);
            case 'y':
                tmax = atof(optarg);
            case 'd':
                dist = atoi(optarg);
            case 's':
                scale = atof(optarg);
            case'?':
                printf("Error: wrong command input option %s\n",optarg);
                break;
            case':':
                printf("Error: Forget to add parameter after option %s\n",optarg);
                break;
        }
    }

    std::cout << "Usage===>\n"
        << "-m SN burst neutrino model " << imod << "\n"
        << "-c LS channel " << icha<< "\n"
        << "-n neutrino mass " << nuMass<< "\n"
        << "-o mass ordering " << MH<< "\n"
        << "-t threshold " << Ethr << "\n"
        << "-x tmin " << tmin << "\n"
        << "-y tmax " << tmax << "\n"
        << std::endl;

    if(icha == -1){
        std::cout << "interaction channel; -c" << std::endl;
        exit(-1);
    }
    if(imod == -1){
        std::cout << "No SN neutrino model specified; -m" << std::endl;
        exit(-1);
    }

    std::cout << "Neutrio Mass = " << nuMass << " eV" << std::endl;

    TString MO[3] = {"NONE", "NO", "IO"};


    //=========================analyze===========================//
    //--------initial SNdetect--------//
    //Int_t icha = 2;
    SNdetect::channelName chaname;
    if(icha == 0)chaname = SNdetect::NuP;
    if(icha == 1)chaname = SNdetect::NuE;
    if(icha == 2)chaname = SNdetect::IBD;
    if(icha == 3)chaname = SNdetect::CEvNS;

    SNdetect* pdet = SNdetect::instance();
    pdet->setSrcModel(imod);
    pdet->getPointerSrc()->setSNDistance(dist);
    pdet->setChannel(chaname);
    Double_t erelres = 3e-2;
    pdet->getPointerEffectLS()->setERelativeRes(erelres);
    pdet->getPointerEffectLS()->setThresholdE(0.);
    pdet->initFCN(); 


    SNeffectLS* peffLS = pdet->getPointerEffectLS();

    if(0) {
    std::cout << "-----> fluence = " << pdet->getPointerSrc()->oneSNFluenceDetAtTime(0.02, 60, 1, 1) << "\n" 
              << "-----> total particle number = " << peffLS-> getNumberOfCarbon() << "\n"
              << "-----> differential cross section = " << peffLS->differentialXS(60, 0.02, 0) << " cm^2" << "\n"
              << "-----> getEventAboveEthrVisAtTime = " << pdet->getEventAboveEthrVisAtTime(0.02, 0.15, 1, 1) << "\n"
              << std::endl;
    }

    // time binning
    Double_t step_t = 0.001;
    Int_t nbin_time = ceil((tmax - tmin)/step_t);
    std::cout << "Time Binning: " <<  tmin << " to " << tmax << " per " << step_t << " s with total Nbins " << nbin_time << std::endl;
    
    // energy binning 

    Double_t Evmin = 0;//pdet->getPointerSrc()->getSNminEv();
    Double_t Evmax = 90; 
    Double_t step_Ev = 1;
    Int_t nbins_Ev = (Evmax-Evmin)/step_Ev;

    Double_t Evismin = Ethr;//0;
    Double_t Evismax;
    Double_t step_Evis;
    if(chaname == SNdetect::NuP){
        Evismax = 5.;
        step_Evis = 0.05;
    }
    if(chaname == SNdetect::NuE || chaname == SNdetect::IBD){
        Evismax = 80.;
        step_Evis = 0.2;
    }
    if(chaname == SNdetect::CEvNS){
        Evismax = 0.3;
        step_Evis = 0.001;
    }
    Int_t nbins_Evis = (Evismax-Evismin)/step_Evis;
    std::cout << "nbins_Evis " << nbins_Evis << std::endl;
    TString chaName[4] = {"pES","eES","IBD", "CEvNS"};



    // get the visible events time spectra above Ethr and save into histogram and root file.
    TString modelName;
    modelName = Form("Garching%d", imod);
    TString fn = Form("%s_PDF_%s_%dkpc_%s_%.2fMeV_%.3fs-%.3fs_scale%.3f_test.root",  modelName.Data(), MO[MH].Data(), dist, chaName[icha].Data(), Ethr, tmin, tmax, scale);
    std::cout << "output filename : " << fn << std::endl;
    TFile* f = new TFile(fn, "recreate");
    TH1D* h1 = new TH1D("h1", "visible energy spectrum rate", nbin_time, tmin*1000-0.5, tmax*1000-0.5);

    double factor = 1;
    double tshift = 0;
    for (int ipt=0; ipt<nbin_time; ipt++) {
        double t = tmin + ipt * step_t + tshift;
        std::cout << "Running bin " << ipt << " time " << t << std::endl;
        double flu = 0;
        for(int ii=0; ii<6; ii++) {
            //std::cout << "In main: " << t << " " << ii << " " << pdet->getEventAboveEthrVisAtTime(t, Ethr, ii, MH) << std::endl;
            flu += pdet->getEventAboveEthrVisAtTime(t, Ethr, ii, MH) * factor / 1000.;  // for 1 ms interval
        }
        h1->SetBinContent(ipt+1, flu*scale);
    }

    h1->Write();
    f->Close();


    //TString modelName;
    //modelName = Form("Garching%d", imod);
    //TString fn = Form("%s_PDF_%s_%s_%dkpc_%.3fs-%.3fs_scale%.3f_test2D.root",  modelName.Data(), chaName[icha].Data(), MO[MH].Data(), dist,  tmin, tmax, scale);
    //std::cout << "output filename : " << fn << std::endl;
    //TFile* f = new TFile(fn, "recreate");
    //TH2D* h1 = new TH2D("h1", "time-visible energy 2D hist", nbin_time, tmin*1000-0.5, tmax*1000-0.5, nbins_Evis, Evismin, Evismax);
    //double factor = 1;
    //double tshift = 0;
    //for (int ipt=0; ipt<nbin_time; ipt++) {
    //    double t = tmin + ipt * step_t + tshift;
    //    for (int iE=0; iE<nbins_Evis; iE++) {
    //        double Evis = Evismin + (iE + 0.5) * step_Evis;
    //        std::cout << "Running time " << t << ", visible energy " << Evis << std::endl;
    //        double flu = 0;
    //        for (int ii=0; ii<6; ii++) {
    //            flu += pdet->getEvisSpectrumAtTime(t, Evis, ii, MH) * factor / 1000;
    //        }
    //        std::cout << "=======> flu = " << flu << std::endl;
    //        h1->SetBinContent(ipt+1, iE+1, flu);
    //    }
    //}

    //h1->Write();
    //f->Close();

}
