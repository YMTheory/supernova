#include <iostream>
#include <cstdlib>
#include <unistd.h>//for getopt()

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

#include "SNdetect.hh"
#include "SNsource.hh"
#include "SNeffectLS.hh"

int main(int argc, char* argv[]){

    Int_t imod = -1;//SN neutrino source model
    Int_t icha = -1;//nu interaction channel
    //==========getopt()===========//
    Int_t opt;
    while((opt = getopt(argc, argv, ":e:m:c:t:"))!= -1){
        switch(opt){
            case'm':
                imod = atoi(optarg);
                break;
            case'c':
                icha = atoi(optarg);
                break;
            case'?':
                printf("Error: wrong command input option %s\n",optarg);
                break;
            case':':
                printf("Error: Forget to add parameter after option %s\n",optarg);
                break;
        }
    }

    std::cout << "Usage===>\n"
        << "-m SN burst neutrino model \n"
        << "-c LS channel \n"
        << std::endl;


    if(icha == -1){
        std::cout << "interaction channel; -c" << std::endl;
        exit(-1);
    }
    if(imod == -1){
        std::cout << "No SN neutrino model specified; -m" << std::endl;
        exit(-1);
    }

    //=========================analyze===========================//
    //--------initial SNdetect--------//
    //Int_t icha = 2;
    SNdetect::channelName chaname;
    if(icha == 0)chaname = SNdetect::NuP;
    if(icha == 1)chaname = SNdetect::NuE;
    if(icha == 2)chaname = SNdetect::IBD;
    SNdetect* pdet = SNdetect::instance();
    pdet->setSrcModel(imod);
    Double_t dist = 10;//kpc
    pdet->getPointerSrc()->setSNDistance(dist);
    pdet->setChannel(chaname);
    Double_t erelres = 3e-2;
    pdet->getPointerEffectLS()->setERelativeRes(erelres);
    Double_t Ethr = 0.2;
    pdet->getPointerEffectLS()->setThresholdE(Ethr);
    pdet->initFCN(); 

    //SN pregenitor mass
    Double_t mass=0;
    if(imod>10000){
        Int_t tomass = ((imod%10000)/10);
        mass = tomass/10.;
    }
    else{
        mass = imod/100;
    }



    //--------histogram------//
    Double_t tmin = -0.2;//s
    Double_t tmax = 20;//s
    Double_t t1_thr = 0.02;
    Double_t t2_thr = 0.6; 
    Double_t step_t0 = 0.002;
    Int_t npt0 = (t1_thr-tmin)/step_t0;
    Double_t step_t1 = 0.03;
    Int_t npt1 = (t2_thr-t1_thr)/step_t1;
    Double_t step_t2 = 0.1;
    Int_t npt2 = (tmax-t2_thr)/step_t2;

    Int_t nbin_time = npt0+npt1+npt2;
    std::cout << "nbin_time  " << nbin_time << std::endl;
    Double_t* binning_time = new Double_t[nbin_time+1];
    for(Int_t ip=0; ip<npt0; ip++){
        binning_time[ip] = tmin+step_t0*ip;
        //std::cout << binning_time[ip] << std::endl;
    }
    for(Int_t ip=npt0; ip<(npt0+npt1); ip++){
        binning_time[ip] = t1_thr+step_t1*(ip-npt0);
        //std::cout << binning_time[ip] << std::endl;
    }
    for(Int_t ip=(npt0+npt1); ip<(npt0+npt1+npt2); ip++){
        binning_time[ip] = t2_thr+step_t2*(ip-npt0-npt1);
        //std::cout << binning_time[ip] << std::endl;
    } 
    binning_time[nbin_time] = tmax;

    //**energy binning
    Double_t Evmin = Ethr;//pdet->getPointerSrc()->getSNminEv();
    Double_t Evmax = 80; 
    Double_t step_Ev = 1;
    Int_t nbins_Ev = (Evmax-Evmin)/step_Ev;

    Double_t Evismin = Ethr;//0;
    Double_t Evismax;
    Double_t step_Evis;
    if(chaname == SNdetect::NuP){
        Evismax = 5.;
        step_Evis = 0.01;
    }
    if(chaname == SNdetect::NuE || chaname == SNdetect::IBD){
        Evismax = 80.;
        step_Evis = 0.2;
    }

    Int_t nbins_Evis = (Evismax-Evismin)/step_Evis;
    std::cout << "nbins_Evis " << nbins_Evis << std::endl;
    TString chaName[3] = {"nup","nue","IBD"};

    //--------output file-----//
    TString path="/junofs/users/lihl/neutrino/timeSN/etSpec";
    TFile* fout = new TFile(Form("%s/rootfiles/evistSpec_mod%d_cha%d.root",path.Data(),imod,icha),"RECREATE");
    Double_t runtime[2];
    TTree* tinfo = new TTree("tinfo","model information");
    tinfo->Branch("model",   &imod,   "model/I");
    tinfo->Branch("channel", &icha,   "channel/I");
    tinfo->Branch("runtime", runtime, "runtime[2]/D"); 
    pdet->getPointerSrc()->getTimeRange(runtime[0],runtime[1],icha);
    tinfo->Fill();

    //----------2D distribution info---------//
    const Int_t nMH=3;
    TH2D* h2d_etSpec[nMH];
    for(Int_t imh=0; imh<nMH; imh++){
        h2d_etSpec[imh] = new TH2D(Form("hET_mod%d_cha%d_mh%d",imod,icha,imh),Form("spectrum of energy and time for model%d",imod),nbin_time,binning_time,nbins_Evis,Evismin, Evismax);
        h2d_etSpec[imh]->GetXaxis()->SetTitle("time [s]");
        h2d_etSpec[imh]->GetYaxis()->SetTitle("E_{vis} [MeV]");
        h2d_etSpec[imh]->GetZaxis()->SetTitle(Form("%s event",chaName[icha].Data()));
    }

    //channel event at [itime,ievis]
    std::cout << "Nbin_Time=" << nbin_time << "; Nbin_Evis=" << nbins_Evis << std::endl;
    for(Int_t ipt=0; ipt<nbin_time; ipt++){
        Double_t itwidth = binning_time[ipt+1]-binning_time[ipt];
        Double_t itime   = binning_time[ipt]+0.5*itwidth;
        for(Int_t ipe=0; ipe<nbins_Evis; ipe++){
            std::cout <<"[" << nbin_time << "," << nbins_Evis << "]====>"
                << "ibin_time " << ipt << "; ibin_evis " << ipe 
                << std::endl;

            Double_t ievis   = Evismin+step_Evis*(ipe+0.5);

            for(Int_t imh=0; imh<nMH; imh++){
                Double_t iflux=0;
                for(Int_t itype=0; itype<6; itype++){
                    iflux += pdet->getEvisSpectrumAtTime(itime,ievis,itype,imh);
                }
                h2d_etSpec[imh]->SetBinContent(ipt+1,ipe+1,iflux*itwidth*step_Evis);
            }

        }
    }

    fout->cd();
    tinfo->Write();
    for(Int_t imh=0; imh<3; imh++){
        h2d_etSpec[imh]->Write();
    }
    fout->Close();

    return 0;

}

