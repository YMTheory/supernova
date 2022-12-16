#include <iostream>
#include <string>
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSystem.h"

#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"

#include "SNeffectLS.hh"
#include "SNsource.hh"
#include "SNnumGarchingSrc.hh"
#include "SNdetect.hh"
#include "SNchannelNuP.hh"
#include "SNchannelNuE.hh"
#include "SNchannelIBD.hh"
#include "SNchannelCEvNS.hh"


SNdetect* SNdetect::pdet = NULL;
SNdetect::SNdetect():
    fEvmin(0),
    fEvmax(0),
    fTmin(0),
    fTmax(0){

        gRandom = new TRandom3(0);
        peffLS = NULL;
        psrc = NULL;
        file = NULL;
        iMCflavor=-1;

        f2rndEvT     = NULL;
        f1rndEv      = NULL;

    }

SNdetect::~SNdetect(){
    if(file) file->Close();
    deletFCN();
    if(psrc)   delete psrc;
    if(peffLS) delete peffLS;
}

SNdetect* SNdetect::instance(){
    if(!pdet){
        pdet = new SNdetect();
        printf("Object SNEdetect (%p) created as singleton...\n",pdet);
    }
    return pdet;
}

void SNdetect::initChannel(channelName cha){
    switch(cha){
        case NuP:
            {
                SNchannelNuP* chaNuP = new SNchannelNuP();
                peffLS = chaNuP->createChannel();
                fchannel = cha;

                if(file)file->Close();
                file = TFile::Open("/junofs/users/miaoyu/supernova/wenlj/simulation/data/quenching/quenchCurve.root");
                if(!file){
                    std::cout << "/junofs/users/miaoyu/supernova/wenlj/simulation/data/quenching/quenchCurve.root can't be opened!" << std::endl;
                    exit(-1);
                }
                grquench    = (TGraph*)file->Get("quenchCurve");
                grquenchInv = (TGraph*)file->Get("quenchCurveInv");

                std::cout << "nu-P channel has been created." << std::endl;
                break;
            }
        case NuE:
            {
                SNchannelNuE* chaNuE = new SNchannelNuE();
                peffLS = chaNuE->createChannel();
                fchannel = cha;
                std::cout << "nu-e channel has been created." << std::endl;
                break;
            }
        case IBD:
            {
                SNchannelIBD* chaIBD = new SNchannelIBD();
                peffLS = chaIBD->createChannel();
                fchannel = cha;
                std::cout << "IBD channel has been created." << std::endl;
                break;
            }
        case CEvNS:
            {
                SNchannelCEvNS* chaCEvNS = new SNchannelCEvNS();
                peffLS = chaCEvNS->createChannel();
                fchannel = cha;
                if(file)file->Close();
                file = TFile::Open("/junofs/users/miaoyu/supernova/wenlj/simulation/data/quenching/quenchCurve.root");
                if(!file){
                    std::cout << "/junofs/users/miaoyu/supernova/wenlj/simulation/data/quenching/quenchCurve.root can't be opened!" << std::endl;
                    exit(-1);
                }
                grquench    = (TGraph*)file->Get("quenchCurve");
                grquenchInv = (TGraph*)file->Get("quenchCurveInv");
                std::cout << "CEvNS channel has been created." << std::endl;
                break;
            }
        default:
            std::cout << "Error:  the channel is not included!" << std::endl;

    }
}

void SNdetect::setSrcModel(int imode){
    if(psrc) delete psrc;
    psrc = new SNnumGarchingSrc(imode);
    std::cout << "numerical Garching SN neutrino model " << imode << " is being used" << std::endl;

}
//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"



//---------------Random generator------------//
//-------------------FCN------------------//
//
//rnd generation

double fcnEvandTES(double* x, double* par){
    double Ev = x[0];
    double T  = x[1];
    int type = int(par[0]);
    int MH   = int(par[1]);

    SNdetect* pdet = SNdetect::instance();
    int ntype = 6;
    double xsweightFlu = 0;
    //all flavor mixing
    if(type==-1){
        for(int it=0; it<ntype; it++){
            double flu = pdet->getPointerSrc()->oneSNFluenceDet(Ev, it,MH);
            double xs  = pdet->getPointerEffectLS()->differentialXS(Ev, T, it);
            xsweightFlu += flu*xs;
        }
    }
    //single flavor
    else{
        double flu = pdet->getPointerSrc()->oneSNFluenceDet(Ev, type, MH);
        double xs  = pdet->getPointerEffectLS()->differentialXS(Ev, T, type);
        xsweightFlu = flu*xs;

    }
    return xsweightFlu;
}

double fcnFlatEvandTES(double* x, double* par){
    double Ev = x[0];
    double T  = x[1];
    int type = int(par[0]);
    int MH   = int(par[1]);

    SNdetect* pdet = SNdetect::instance();
    int ntype = 6;
    double xsweightFlu = 0;
    //all flavor mixing
    if(type==-1){
        for(int it=0; it<ntype; it++){
            //double flu = 1./pdet->getPointerEffectLS()->totalXS(Ev, it);
            double flu = 1;
            double xs  = pdet->getPointerEffectLS()->differentialXS(Ev, T, it);
            xsweightFlu += flu*xs;
        }
    }
    //single flavor
    else{
        double flu = 1;
        //double flu = 1./pdet->getPointerEffectLS()->totalXS(Ev, type);
        double xs  = pdet->getPointerEffectLS()->differentialXS(Ev, T, type);
        xsweightFlu = flu*xs;
    }
    return xsweightFlu;

}

double fcnEvCC(double* x, double* par){
    double Ev = x[0];
    int type  = int(par[0]);
    int MH    = int(par[1]);

    SNdetect* pdet = SNdetect::instance();
    double flu = pdet->getPointerSrc()->oneSNFluenceDet(Ev, type, MH);
    double xs  = pdet->getPointerEffectLS()->totalXS(Ev, type);

    return flu*xs;
}

double fcnEvCNC(double* x, double* par){
    double Ev = x[0];
    int MH = int(par[0]);

    SNdetect* pdet = SNdetect::instance();
    double xsweightFlu = 0;
    for(int it=0; it<6; it++){
        double flu = pdet->getPointerSrc()->oneSNFluenceDet(Ev, it, MH);
        double xs  = pdet->getPointerEffectLS()->totalXS(Ev, it);
        xsweightFlu += flu*xs;
    }
    return xsweightFlu;
}



//------------generating--------------//

void SNdetect::initFCN(){
    fEvmin = psrc->getSNminEv();
    fEvmax = psrc->getSNmaxEv();
    fEvmax = 90;    // concerns for Burrows2D dataset
    std::cout << "START of initFCN ..." << std::endl;
    // test line
    // std::cout << "initFCN: fEvmin, fEvmax " << fEvmin << ", " << fEvmax << std::endl;
    //
    fEthr  = peffLS->getThresholdE();
    switch(fchannel){
        case NuP:
            {
                //double mp = 938.;
                fTmin = 0.;
                //fTmax = 2.*fEvmax*fEvmax/mp;
                fTmax = getTmax();
                f2rndEvT = new TF2("f2rndEvT", fcnEvandTES, fEvmin, fEvmax, fTmin, fTmax, 2);
                //std::cout << fcnEvandTES << " " << fEvmin << " " << fEvmax << " " << fTmin << " " << fTmax << std::endl;
                f2rndEvT->SetNpx(200);
                f2rndEvT->SetNpy(200);

                f2rndFlatEvT = new TF2("f2rndFlatEvT", fcnFlatEvandTES, fEvmin, fEvmax, fTmin, fTmax,2);
                f2rndFlatEvT->SetNpx(200);
                f2rndFlatEvT->SetNpy(200);
                break;
            }
        case NuE:
            {
                double me = 0.511;
                fTmin = 0.;
                fTmax = fEvmax/(1+me/(2.*fEvmax));

                f2rndEvT = new TF2("f2rndEvT", fcnEvandTES, fEvmin, fEvmax, fTmin, fTmax, 2);
                f2rndEvT->SetNpx(200);
                f2rndEvT->SetNpy(200);

                f2rndFlatEvT = new TF2("f2rndFlatEvT", fcnFlatEvandTES, fEvmin, fEvmax, fTmin, fTmax,2);
                f2rndFlatEvT->SetNpx(200);
                f2rndFlatEvT->SetNpy(200);
                break;
            }
        case IBD:
            {
                fTmin = 0;
                fTmax = fEvmax-1.3;
                f1rndEv = new TF1("f1rndEv", fcnEvCC, fEvmin, fEvmax, 2);
                f1rndEv->SetNpx(200);
                break;
            }

    }
    std::cout << "END of initFCN ..." << std::endl;

}

void SNdetect::deletFCN(){
    if(fchannel == NuP || fchannel == NuE){
        delete f2rndEvT;
        delete f2rndFlatEvT;
    }
    if(fchannel == IBD){
        delete f1rndEv;
    }
}

void SNdetect::resetFCN(){
    deletFCN();
    initFCN();
}


//--------------Analytical result------------------//
//
//---------------------FCN---------------------//
double fcnAnaEv2T(double* x, double* par){
    double Ev = x[0];
    double T     = par[0];
    int type     = int(par[1]);
    int MH       = int(par[2]);

    SNdetect* pdet = SNdetect::instance();
    int ntype = 6;
    double totflu = 0;
    if(type==-1){
        for(int it=0; it<ntype; it++){
            double flu = pdet->getPointerSrc()->oneSNFluenceDet(Ev, it, MH);
            double dxs = pdet->getPointerEffectLS()->differentialXS(Ev, T, it);
            totflu += flu*dxs;
        }
    }
    else{
        double flu = pdet->getPointerSrc()->oneSNFluenceDet(Ev, type, MH);
        double dxs = pdet->getPointerEffectLS()->differentialXS(Ev, T, type);
        totflu = flu*dxs;

    }

    return totflu;
}






//-----------------------time distribution---------------//
//=============for a specific time============//
//FCN
double fcnAnaEv2TatTime(double* x, double* par){
    double Ev = x[0];
    double time = par[0];
    double T    = par[1];
    int    type = int(par[2]);
    int    MH   = int(par[3]);

    SNdetect* pdet = SNdetect::instance();
    int ntype = 6;
    double totflu = 0;
    if(type==-1){
        for(int it=0; it<ntype; it++){
            double flu = pdet->getPointerSrc()->oneSNFluenceDetAtTime(time, Ev, it, MH);
            double dxs = pdet->getPointerEffectLS()->differentialXS(Ev, T, it);
            totflu += flu*dxs;
        // test line
        //std::cout << "fcnAnaEv2TatTime flu, dxs: " << flu << ", " << dxs << std::endl;
        }
    }
    else{
        double flu = pdet->getPointerSrc()->oneSNFluenceDetAtTime(time, Ev, type, MH);
        double dxs = pdet->getPointerEffectLS()->differentialXS(Ev, T, type);
        totflu = flu*dxs;
        // test line
        //std::cout << "fcnAnaEv2TatTime flu, dxs: " << flu << ", " << dxs << std::endl;
        //
    }

    return totflu;
}

double fcnEventsVisatTime(double* x, double* par){
    double Evis = x[0];
    double time = par[0];
    int    type = int(par[1]);
    int    MH   = int(par[2]);

    SNdetect* pdet = SNdetect::instance();
    //std::cout << "NumEvis " << pdet->getEvisSpectrumAtTime(time, Evis, type, MH) << std::endl;
    return pdet->getEvisSpectrumAtTime(time, Evis, type, MH);
}


double SNdetect::getXSweightedEvSpectrumAtTime(double time, double Ev, int type, int MH){
    double Npar;
    //number of protons
    if(fchannel == NuP || fchannel == IBD){
        Npar = peffLS->getNumberOfProton();
    }   
    //number of electrons
    if(fchannel == NuE){
        //Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon();
        Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon() + 8*peffLS->getNumberOfOxygen() ;
    }
    if(fchannel == CEvNS) {
        Npar = peffLS->getNumberOfCarbon();
    }
    //number of carbon-12, 98.9% of carbons
    int ntype = 6;
    double totflu = 0;
    if(type==-1){
        for(Int_t it=0; it<ntype; it++){
            double flu = psrc->oneSNFluenceDetAtTime(time,Ev, it,MH);
            double xs  = peffLS->totalXS(Ev, it);
            totflu += flu*xs;
        }
    }
    else{
        double flu = psrc->oneSNFluenceDetAtTime(time, Ev, type, MH);
        //std::cout << "flu " << flu << std::endl;
        double xs  = peffLS->totalXS(Ev, type);
        //std::cout << time << " " << Ev << " " << Npar << " " << flu << " " << xs << std::endl;
        totflu = flu*xs;
    }
    return totflu*Npar;
}



double SNdetect::getTSpectrumAtTime(double time, double T, int type, int MH){
    double Evmin;
    double Npar;
    if(fchannel == NuP || fchannel == NuE || fchannel == CEvNS){
        if(fchannel == NuP){
            double mp = 938;//MeV
            Evmin = TMath::Sqrt(mp*T/2.);
            Npar  = peffLS->getNumberOfProton();
        }
        if(fchannel == NuE){
            double me = 0.51;//MeV
            Evmin = 0.5*(T+TMath::Sqrt(T*T+2*T*me));
            Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon();
            //Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon() + 8*peffLS->getNumberOfOxygen() ;
        }
        if(fchannel == CEvNS) {
            double mC12 = 11179.01;
            Evmin = TMath::Sqrt(mC12*T/2);
            Npar = peffLS->getNumberOfCarbon();
        }

        // test line
        //std::cout << "---> Evmin, fEvmax, Npar: " << Evmin << ", " << fEvmax << ", " << Npar << std::endl;
        //

        TF1 f("anaEv2TatTime",fcnAnaEv2TatTime, fEvmin, fEvmax, 4);
        f.SetParameter(0, time);
        f.SetParameter(1, T);
        f.SetParameter(2, type);
        f.SetParameter(3, MH);
        ROOT::Math::WrappedTF1 wf1(f);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(1e-3);
        double spect = Npar*ig.Integral(Evmin, fEvmax);
        // test line
        // std::cout << "spect: " << spect << std::endl;
        return spect;
    }
    if(fchannel==IBD){
        return getXSweightedEvSpectrumAtTime(time,getEvFromT(T), type, MH);
    }
    return 0;
}

double SNdetect::getEvisSpectrumAtTime(double time, double Evis, int type, int MH){
    if(fchannel==NuP or fchannel==CEvNS){
        double T = grquenchInv->Eval(Evis);
        double dTqdTp = (grquench->Eval(T+0.005)-grquench->Eval(T))/0.005; 
        return getTSpectrumAtTime(time, getTFromEvis(Evis),type, MH)/dTqdTp;
    }
    else{
        return getTSpectrumAtTime(time, getTFromEvis(Evis), type, MH);
    }

    return 0;
}


double SNdetect::getEventAboveEthrVisAtTime(double time, double Ethr, int type, int MH){
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                //fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
                double m_Tmax = 2 * fEvmax * fEvmax / mp;
                fEvismax = grquench->Eval(m_Tmax);
                break;
            }
        case NuE:
            {
                double me = 0.51;
                fEvismax = fEvmax/(1+me/(2.*fEvmax));
                break;
            }
        case IBD:
            {
                fEvismax = fEvmax-0.8;
                break;
            }
        case CEvNS:
            {
                double mC12 = 11179.01;
                double m_Tmax = 2 * fEvmax * fEvmax / mC12;
                fEvismax = grquench->Eval(m_Tmax);
                break;
            }
    }

    TF1 f("eveVisatTime", fcnEventsVisatTime, 0, fEvismax, 3);
    f.SetParameter(0, time);
    f.SetParameter(1, type);
    f.SetParameter(2, MH);
    
    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
    ig.SetFunction(wf1);
    ig.SetRelTolerance(1e-3);
    double integ = ig.Integral(Ethr, fEvismax);

    //int np = 500;
    //double *x = new double[np];
    //double *w = new double[np];
    //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
    //double integ = f.IntegralFast(np,x,w,Ethr,fEvismax);
    //delete[] x;
    //delete[] w;
    return integ;
    return 0;
}



double SNdetect::getEvisFromT(double T){
    switch(fchannel){
        case(NuP):
            return grquench->Eval(T);
        case(NuE):
            return T;
        case(IBD):
            return T+0.51;
        case(CEvNS):
            return grquench->Eval(T);
    }
    return 0;
}

double SNdetect::getTFromEvis(double Evis){
    switch(fchannel){
        case(NuP):
            return grquenchInv->Eval(Evis);
        case(NuE):
            return Evis;
        case(IBD):
            return Evis-0.51;
        case(CEvNS):
            return grquenchInv->Eval(Evis);
    }
    return 0;
}

double SNdetect::getTFromEv(double Ev){
    switch(fchannel){
        case(NuP):
            {
                std::cout << "NuP: T-Ev is not 1 to 1 related." << std::endl;
                return -1;
            }
        case(NuE):
            {    
                std::cout << "NuE: T-Ev is not 1 to 1 related." << std::endl;
                return -1;
            }
        case(IBD):
            return Ev-1.3;
        case(CEvNS):
            {
                std::cout << "CEvNS: T-Ev is not 1 to 1 related." << std::endl;
                return -1;
            }
    }
    return 0;
}
double SNdetect::getEvFromT(double T){
    switch(fchannel){
        case(NuP):
            {
                std::cout << "NuP: T-Ev is not 1 to 1 related." << std::endl;
                break;
            }
        case(NuE):
            {    
                std::cout << "NuE: T-Ev is not 1 to 1 related." << std::endl;
                break;
            }
        case(IBD):
            {
                return T+1.3;
                break;
            }
        case(CEvNS):
            {
                std::cout << "CEvNS: T-Ev is not 1 to 1 related." << std::endl;
                break;
            }
        default:
            std::cout << "Nothing..." << std::endl;
    }
    return -1;
}

double SNdetect::getTmax(){
    //for NuP
    double tmin = 0;
    double tmax = 15;
    TF1* f1 = new TF1("f1","TMath::Sqrt((1+x/2./938)*(938*x/2.))+x/2.",tmin,tmax);
    //TF1* f1 = new TF1("f1","TMath::Sqrt(938*x/2)",tmin,tmax);
    double Ev[200];
    double T[200];
    double step = (tmax-tmin)/200;
    for(int ip=0; ip<200; ip++){
        T[ip] = step*ip+tmin;
        Ev[ip] = f1->Eval(T[ip]);
    }
    TGraph* gr = new TGraph(200,Ev, T);
    double tthre = gr->Eval(fEvmax);

    delete f1;
    delete gr;
    return tthre;
}