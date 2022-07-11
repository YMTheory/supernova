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
#include "SNBurrowsIntegFcn.hh"
#include "SNnumBurrowsSrc.hh"

int main(int argc, char* argv[]) {

    Int_t imod = -1;//SN neutrino source model
    Int_t icha = -1;//nu interaction channel
    Double_t nuMass = 0.0;// nu mass
    Int_t MH = 0 ; // mass ordering
    Double_t Ethr = 0.2 ; // energy threshold MeV
    //==========getopt()===========//
    Int_t opt;
    while((opt = getopt(argc, argv, ":e:m:n:c:t:o:"))!= -1){
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
    SNdetect* pdet = SNdetect::instance();
    pdet->setSrcModel(imod);
    Double_t dist = 10;//kpc
    pdet->getPointerSrc()->setSNDistance(dist);
    pdet->setChannel(chaname);
    Double_t erelres = 3e-2;
    pdet->getPointerEffectLS()->setERelativeRes(erelres);
    //Double_t Ethr = 0.;
    pdet->getPointerEffectLS()->setThresholdE(0.);
    pdet->initFCN(); 

    SNnumGarchingSrc* modelSrc = (SNnumGarchingSrc*)pdet->getPointerSrc();

    SNeffectLS* peffLS = pdet->getPointerEffectLS();
    
    // time binning
    Double_t tmin = 0.00;
    Double_t tmax = 0.060;
    Double_t step_t = 0.001;
    Int_t nbin_time = int((tmax - tmin)/step_t);
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

    Int_t nbins_Evis = (Evismax-Evismin)/step_Evis;
    std::cout << "nbins_Evis " << nbins_Evis << std::endl;
    TString chaName[3] = {"pES","eES","IBD"};



    //SNGarchingIntegFcn* pfcn = new SNGarchingIntegFcn(imod);
    //SNBurrowsIntegFcn* pfcn  = new SNBurrowsIntegFcn(imod);
    /*
    TFile* fout = new TFile("/junofs/users/miaoyu/supernova/production/PDFs/NeutrinoEnuTspec_Burrows2D.root", "recreate");
    TH2D* h2d_nue = new TH2D("h2d_nue", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nux = new TH2D("h2d_nux", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nuebar = new TH2D("h2d_nuebar", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nuxbar = new TH2D("h2d_nuxbar", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);

    for (Int_t ipt = 0; ipt<nbin_time; ipt++) {
        double timeTmp = tmin + step_t * (ipt + 0.5);
        for (Int_t i = 0; i < nbins_Ev; i++) {
            double EvTmp = Evmin + step_Ev * i;
            double nevt;
            nevt = pfcn->getEventAtTime(timeTmp, EvTmp, 0);
            h2d_nue->Fill(timeTmp, EvTmp, nevt);
            nevt = pfcn->getEventAtTime(timeTmp, EvTmp, 1);
            h2d_nuebar->Fill(timeTmp, EvTmp, nevt);
            nevt = pfcn->getEventAtTime(timeTmp, EvTmp, 2);
            h2d_nux->Fill(timeTmp, EvTmp, nevt);
            nevt = pfcn->getEventAtTime(timeTmp, EvTmp, 3);
            h2d_nuxbar->Fill(timeTmp, EvTmp, nevt);
        }
    }

    fout->cd();
    h2d_nue->Write();
    h2d_nuebar->Write();
    h2d_nux->Write();
    h2d_nuxbar->Write();
    fout->Close();
    */
    TString modelName;
    if (imod > 6000 and imod < 7000)  { imod -= 6500; modelName = Form("Burrows2D%d", imod); }
    else
        modelName = Form("Garching%d", imod);
    

    TFile* f = new TFile(Form("/junofs/users/miaoyu/supernova/production/PDFs/10kpc/%s_PDF_%s_10kpc_%s_%.2fMeV.root", modelName.Data(), MO[MH].Data(), chaName[icha].Data(), Ethr), "recreate");
    TH1D* h1 ;
    if (imod > 0 and imod < 100)
        h1 = new TH1D("h1", "visible energy spectrum rate", 60, -0.5, 60.5);
    else
        h1 = new TH1D("h1", "visible energy spectrum rate", 60, -30.5, 30.5);

    double factor;
    if (imod > 0 and imod < 100) {
        factor = 1e50;
    }
    else {
        factor = 1;
    }
    double tshift;
    if (imod > 0 and imod < 100) {
        tshift = 0;
    }
    else {
        tshift = -0.03;
    }
    for (int ipt=0; ipt<nbin_time; ipt++) {
        double t = tmin + ipt * step_t + tshift;
        double flu = 0;
        for(int ii=0; ii<6; ii++) {
            //std::cout << t << " " << ii << " " << pdet->getEventAboveEthrVisAtTime(t, Ethr, ii, MH) << std::endl;
            flu += pdet->getEventAboveEthrVisAtTime(t, Ethr, ii, MH) * factor / 1000.;  // for 1 ms interval
        }
        h1->SetBinContent(ipt+1, flu);
    }

    h1->Write();
    f->Close();

    
    //TFile* fout = new TFile("/junofs/users/miaoyu/supernova/production/PDFs/NeutrinoAverageET_Burrows2D.root", "recreate");
    //TH1D* h1d_nue = new TH1D("h1d_nue", "", nbin_time, tmin, tmax);
    //TH1D* h1d_nux = new TH1D("h1d_nux", "", nbin_time, tmin, tmax);
    //TH1D* h1d_nuebar = new TH1D("h1d_nuebar", "", nbin_time, tmin, tmax);
    //TH1D* h1d_nuxbar = new TH1D("h1d_nuxbar", "", nbin_time, tmin, tmax);


    ////SNGarchingIntegFcn* pgarfcn = new SNGarchingIntegFcn(11);
    //for (Int_t ipt = 0; ipt<nbin_time; ipt++) {
    //    double timeTmp = tmin + step_t * (ipt + 0.5);
    //    double lum = 0;

    //    lum = pfcn->getAverageET(timeTmp, 0) ;
    //    h1d_nue->SetBinContent(ipt+1, lum);
    //    lum = pfcn->getAverageET(timeTmp, 1) ;
    //    h1d_nuebar->SetBinContent(ipt+1, lum);
    //    lum = pfcn->getAverageET(timeTmp, 2) ;
    //    h1d_nux->SetBinContent(ipt+1, lum);
    //    lum = pfcn->getAverageET(timeTmp, 3) ;
    //    h1d_nuxbar->SetBinContent(ipt+1, lum);
    //}

    //fout->cd();
    //h1d_nue->Write();
    //h1d_nux->Write();
    //h1d_nuebar->Write();
    //h1d_nuxbar->Write();
    //fout->Close();
     

    /*
    
    TFile* fout = new TFile("/junofs/users/miaoyu/supernova/production/PDFs/NeutrinFlux_82503.root", "recreate");
    TH2D* h2d_nue_mh1 = new TH2D("h2d_nue_mh1", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nux_mh1 = new TH2D("h2d_nux_mh1", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nuebar_mh1 = new TH2D("h2d_nuebar_mh1", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nuxbar_mh1 = new TH2D("h2d_nuxbar_mh1", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nue_mh2 = new TH2D("h2d_nue_mh2", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nux_mh2 = new TH2D("h2d_nux_mh2", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nuebar_mh2 = new TH2D("h2d_nuebar_mh2", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);
    TH2D* h2d_nuxbar_mh2 = new TH2D("h2d_nuxbar_mh2", "", nbin_time, tmin, tmax, nbins_Ev, Evmin, Evmax);

    for(Int_t imh=1; imh<3; imh++){
        for (Int_t ipt = 0; ipt<nbin_time; ipt++) {
            // time loop
            std::cout << "MH " << imh << ", time bin " << ipt << std::endl;

                for (Int_t j=0; j<nbins_Ev; j++){
                    Double_t EvTmp = Evmin + step_Ev * (j+0.5); 
                    double timeTmp = tmin + step_t * (ipt + 0.5);

                    double fluence;
                    fluence = modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 0, 1);
                    h2d_nue_mh1->Fill(timeTmp, EvTmp, fluence);
                    fluence = modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 0, 2);
                    h2d_nue_mh2->Fill(timeTmp, EvTmp, fluence);
                    fluence = modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 1, 1);
                    h2d_nuebar_mh1->Fill(timeTmp, EvTmp, fluence);
                    fluence = modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 1, 2);
                    h2d_nuebar_mh2->Fill(timeTmp, EvTmp, fluence);
                    fluence = modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 2, 1) + modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 4, 1);
                    h2d_nux_mh1->Fill(timeTmp, EvTmp, fluence);
                    fluence = modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 2, 2) + modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 4, 2);
                    h2d_nux_mh2->Fill(timeTmp, EvTmp, fluence);
                    fluence = modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 3, 1) + modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 5, 1);
                    h2d_nuxbar_mh1->Fill(timeTmp, EvTmp, fluence);
                    fluence = modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 3, 2) + modelSrc->snFluenceDetAtTime(timeTmp, 0.0, EvTmp, 5, 2);
                    h2d_nuxbar_mh2->Fill(timeTmp, EvTmp, fluence);

            }
        }
    }

    fout->cd();
    h2d_nue_mh1->Write();
    h2d_nux_mh1->Write();
    h2d_nuebar_mh1->Write();
    h2d_nuxbar_mh1->Write();
    h2d_nue_mh2->Write();
    h2d_nux_mh2->Write();
    h2d_nuebar_mh2->Write();
    h2d_nuxbar_mh2->Write();
    fout->Close();
    */

    /* 
    //--------output file-----//
    TString path="/junofs/users/miaoyu/supernova/production/PDFs/10kpc";
    TFile* fout = new TFile(Form("%s/evistSpec_mod%d_cha%d_nuMass%.1f.root",path.Data(),imod,icha, nuMass),"RECREATE");
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
        h2d_etSpec[imh] = new TH2D(Form("hET_mod%d_cha%d_mh%d",imod,icha,imh),Form("spectrum of energy and time for model%d",imod),nbin_time, tmin, tmax, nbins_Evis,Evismin, Evismax);
        h2d_etSpec[imh]->GetXaxis()->SetTitle("time [s]");
        h2d_etSpec[imh]->GetYaxis()->SetTitle("E_{vis} [MeV]");
        h2d_etSpec[imh]->GetZaxis()->SetTitle(Form("%s event",chaName[icha].Data()));
    }


    // Modify PDFs due to neutrino masses
    
    Double_t Npar;
    if (chaname == SNdetect::NuP or chaname==SNdetect::IBD) {
        Npar = peffLS->getNumberOfProton();
    }
    if (chaname == SNdetect::NuE) {
        Npar = peffLS->getNumberOfProton() + 6*peffLS->getNumberOfCarbon();
    }

    for (int imh=1; imh<2; imh++) {
        double tmp_Evis[nbins_Evis];
        
        for (int ipt=0; ipt<nbin_time; ipt++) {
            double timeTmp = tmin + step_t * (ipt + 0.5);
            std::cout << "Mass ordering : " << imh << " , time bin number : " << ipt << " " << timeTmp << std::endl;
            for(int i=0; i<nbins_Evis; i++) { tmp_Evis[i] = 0.;}

            for(int k=0; k<nbins_Ev; k++) {
                double EvTmp   = Evmin + step_Ev * (k + 0.5);
                
                if (chaname == SNdetect::NuE or chaname == SNdetect::NuP) {
                    for (int j=0; j<nbins_Evis; j++) {
                        double EvisTmp = Evismin + step_Evis * (j + 0.5);
                        for(int ii=0; ii<6; ii++) {
                            double fluence = modelSrc->snFluenceDetAtTime(timeTmp, nuMass, EvTmp, ii, imh);
                            double T = pdet->getTFromEvis(EvisTmp);
                            double dxsdT   = peffLS->differentialXS(EvTmp, T, ii);
                            tmp_Evis[j] += Npar * fluence * dxsdT * step_Evis * step_Ev;
                        }
                    }
                }
            }

            if(chaname == SNdetect::IBD) {
                for (int j=0; j<nbins_Evis; j++) {
                    double EvisTmp = Evismin + step_Evis * (j + 0.5);
                    double EvTmp   = EvisTmp - 0.78; 
                    double tot_fluence = 0;
                    if (EvTmp <= 0) {
                        tot_fluence = 0;
                    } else {
                        double fluence = modelSrc->snFluenceDetAtTime(timeTmp, nuMass, EvTmp, 1, imh);
                        double totxs   = peffLS->totalXS(EvTmp, 1);
                        tot_fluence = Npar * fluence * totxs * step_Ev * step_Evis;
                    }
                    tmp_Evis[j] += tot_fluence;
                }
            }

            for(int jj=0; jj<nbins_Evis; jj++) {
                double EvisTmp = Evismin + step_Evis * (jj+0.5);
                h2d_etSpec[imh]->Fill(timeTmp, EvisTmp, tmp_Evis[jj]);
            }
        }
    }


    fout->cd();
    tinfo->Write();
    for(Int_t imh=1; imh<3; imh++){
        h2d_etSpec[imh]->Write();
    }
    fout->Close();
    */
    
    return 0;

}
