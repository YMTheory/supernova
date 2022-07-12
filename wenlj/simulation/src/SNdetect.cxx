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
#include "SNanaSrc.hh"
#include "SNnumJapanSrc.hh"
#include "SNnumGarchingSrc.hh"
#include "SNnumPreGuoSrc.hh"
#include "SNnumPreKatoSrc.hh"
#include "SNnumBurrowsSrc.hh"
#include "SNnumNakazatoSrc.hh"
#include "SNdetect.hh"
#include "SNchannelNuP.hh"
#include "SNchannelNuE.hh"
#include "SNchannelIBD.hh"
#include "SNchannelNCC.hh"
#include "SNchannelBCC.hh"
#include "SNchannelCNC.hh"


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

                //quenching curve
                //std::string pathnow(gSystem->WorkingDirectory()); 
                //std::string nowstring = "/SNsim/simulation";
                //std::string targetstring = "/SNsim/simulation/data";

                //size_t last = pathnow.rfind(nowstring);
                //if(last != std::string::npos){
                //    size_t lenpath = pathnow.length();
                //    pathnow.erase(last, lenpath);
                //    pathnow.append(targetstring);
                //}
                //std::cout << pathnow << std::endl;
                //file = TFile::Open(Form("%s/quenchCurve.root", pathnow.data()));

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
        case NCC:
            {
                SNchannelNCC* chaNCC = new SNchannelNCC();
                peffLS = chaNCC->createChannel();
                fchannel = cha;
                std::cout << "NCC channel has been created." << std::endl;
                break;
            }
        case BCC:
            {
                SNchannelBCC* chaBCC = new SNchannelBCC();
                peffLS = chaBCC->createChannel();
                fchannel = cha;
                std::cout << "BCC channel has been created." << std::endl;
                break;
            }
        case CNC:
            {
                SNchannelCNC* chaCNC = new SNchannelCNC();
                peffLS = chaCNC->createChannel();
                fchannel = cha;
                std::cout << "CNC channel has been created." << std::endl;
                break;
            }
        default:
            std::cout << "Error:  the channel is not included!" << std::endl;

    }
}

void SNdetect::setSrcModel(int imode){
    if(psrc) delete psrc;
    if(imode == 0){
        psrc = new SNanaSrc();
        std::cout << "Analytical SN neutrino model is being used." << std::endl;
    }
    if(imode > 6000 && imode < 7999){
        std::cout << "numerical Burrows SN neutrino model " << imode << " is being used" << std::endl;
        psrc = new SNnumBurrowsSrc(imode);
    }
    if(imode == 11 or (imode > 70000 && imode < 100000)){
        psrc = new SNnumGarchingSrc(imode);
        std::cout << "numerical Garching SN neutrino model " << imode << " is being used" << std::endl;
    }
    //if(imode > 1000 && imode < 9999){
    //    psrc = new SNnumJapanSrc(imode);
    //    std::cout <<"numerical Japan SN neutrino model " << imode << " is being used" << std::endl;
    //}

    if(imode > 2000 and imode < 3000) {
        psrc = new SNnumNakazatoSrc(imode);
        std::cout << "numerical Nakazato SN neutrino model " << imode << " is being used" << std::endl;
    }

    //pre-SN models 
    if(imode>10000 && imode < 11000){
        psrc = new SNnumPreGuoSrc(imode);
        std::cout <<"numerical pre SN neutrino model from Guo " << imode << " is being used" << std::endl;

    }
    if(imode>20000 && imode < 21000){
        psrc = new SNnumPreKatoSrc(imode);
        std::cout <<"numerical pre SN neutrino model from Kato " << imode << " is being used" << std::endl;
    }


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
        case BCC:
            {
                fTmin = 0.;
                fTmax = fEvmax-13.88;
                f1rndEv = new TF1("f1rndEv", fcnEvCC, fEvmin, fEvmax, 2);
                f1rndEv->SetNpx(200);
                break;
            }
        case NCC:
            {
                fTmin = 0.;
                fTmax = fEvmax-16.827;
                f1rndEv = new TF1("f1rndEv", fcnEvCC, fEvmin, fEvmax, 2);
                f1rndEv->SetNpx(200);
                break;

            }
        case CNC:

            {
                fTmin = 15.;
                fTmax = 16.;
                f1rndEv = new TF1("f1rndEv", fcnEvCNC, fEvmin, fEvmax, 1);
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

//-------------Randomly generation of energy spectra-----------//
void SNdetect::setMCflavor(int type, int MH){
    iMCflavor = type;
    if(fchannel == CNC){
        f1rndEv->SetParameter(0,MH);
    }
    else if(fchannel == NuP || fchannel == NuE){
        f2rndEvT->SetParameter(0, type);
        f2rndEvT->SetParameter(1, MH);
        f2rndFlatEvT->SetParameter(0,type);
        f2rndFlatEvT->SetParameter(1,MH);
    }
    else{
        f1rndEv->SetParameter(0, type);
        f1rndEv->SetParameter(1, MH);
    }

}
void SNdetect::getRandomEvandT(double& Ev, double& T){
    switch(fchannel){
        case NuP:
            {
                f2rndEvT->GetRandom2(Ev, T);
                break;
            }
        case NuE:
            {
                f2rndEvT->GetRandom2(Ev, T);
                break;
            }
        case IBD:
            {
                if(iMCflavor==-1 || iMCflavor == 1){
                    Ev = f1rndEv->GetRandom();
                    T  = getTFromEv(Ev);
                    if(T<0) T=0;
                }
                else{
                    Ev = 0;
                    T  = 0;
                }
                break;
            }
        case NCC:
            {
                if(iMCflavor==-1 || iMCflavor == 0){
                    Ev = f1rndEv->GetRandom();
                    T  = getTFromEv(Ev);
                    if(T<0) T=0;
                }
                else{
                    Ev = 0;
                    T  = 0;
                }
                break;
            }

        case BCC:
            {
                if(iMCflavor==-1 || iMCflavor == 1){
                    Ev = f1rndEv->GetRandom();
                    T  = getTFromEv(Ev);
                    if(T<0) T=0;
                }
                else{
                    Ev = 0;
                    T  = 0;
                }
                break;

            }

        case CNC:
            {
                Ev = f1rndEv->GetRandom();
                T  = 15.11;
                break;
            }
    }    
}

void SNdetect::getRandomFlatEvandT(double& Ev, double& T){
    switch(fchannel){
        case NuP:
            {
                f2rndFlatEvT->GetRandom2(Ev, T);
                break;
            }
        case NuE:
            {
                f2rndEvT->GetRandom2(Ev, T);
                break;
            }
        case IBD:
            {
                if(iMCflavor==-1 || iMCflavor == 1){
                    Ev = f1rndEv->GetRandom();
                    T  = getTFromEv(Ev);
                    if(T<0) T=0;
                }
                else{
                    Ev = 0;
                    T  = 0;
                }
                break;
            }
        case NCC:
            {
                if(iMCflavor==-1 || iMCflavor == 0){
                    Ev = f1rndEv->GetRandom();
                    T  = getTFromEv(Ev);
                    if(T<0) T=0;
                }
                else{
                    Ev = 0;
                    T  = 0;
                }
                break;
            }

        case BCC:
            {
                if(iMCflavor==-1 || iMCflavor == 0){
                    Ev = f1rndEv->GetRandom();
                    T  = getTFromEv(Ev);
                    if(T<0) T=0;
                }
                else{
                    Ev = 0;
                    T  = 0;
                }
                break;
            }

        case CNC:
            {
                Ev = gRandom->Uniform(fEvmin, fEvmax);
                T  = 15.11;
                break;
            }
    }    

}


double SNdetect::getRandomEvis(double T){
    return getEvisFromT(T);
}

double SNdetect::getRandomEobs(double Evis){
    double sigE = peffLS->getEres(Evis);
    double Eobs = gRandom->Gaus(Evis, sigE);
    if(Eobs<0)Eobs=0;
    return Eobs;
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

double fcnAnaT2Eobs(double* x, double* par){
    double T    = x[0];
    double Eobs = par[0];
    int type    = int(par[1]);
    int MH      = int(par[2]);

    SNdetect* pdet = SNdetect::instance();
    double Evis = pdet->getEvisFromT(T);
    double prob = pdet->getPointerEffectLS()->getProbResGauss(Eobs, Evis);
    double flu  = pdet->getTSpectrum(T, type,MH);
    //if(isnan(flu))std::cout<< "T " << T << std::endl;

    return prob*flu;
}

//To calculate number of events
double fcnEventsVis(double* x, double* par){
    double Evis = x[0];
    int type = int(par[0]);
    int MH   = int(par[1]);

    SNdetect* pdet = SNdetect::instance();
    return pdet->getEvisSpectrum(Evis, type, MH);
}
double fcnEventsVisCNC(double* x, double* par){
    double Ev   = x[0];
    int type    = int(par[0]);
    int MH      = int(par[1]);

    SNdetect* pdet = SNdetect::instance();
    return pdet->getXSweightedEvSpectrum(Ev, type, MH);
}

double fcnEventsObs(double* x, double* par){
    double Eobs = x[0];
    int type = (int)par[0];
    int MH   = (int)par[1];

    SNdetect* pdet = SNdetect::instance();
    return pdet->getEobsSpectrum(Eobs, type, MH);
}





//-------------Analytical-------------//
double SNdetect::getXSweightedEvSpectrum(double Ev, int type, int MH){
    double Npar;

    //number of protons
    if(fchannel == NuP || fchannel == IBD){
        Npar = peffLS->getNumberOfProton();
    }   
    //number of electrons
    if(fchannel == NuE){
        Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon();
    }
    //number of carbon-12, 98.9% of carbons
    if(fchannel == NCC || fchannel == BCC || fchannel == CNC){
        Npar = 0.989*peffLS->getNumberOfCarbon();
    }
    int ntype = 6;
    double totflu = 0;
    if(type == -1){
        for(Int_t it=0; it<ntype; it++){
            double flu = psrc->oneSNFluenceDet(Ev, it, MH);
            double xs  = peffLS->totalXS(Ev, it);
            totflu += flu*xs;
        }
    }
    else{
        double flu = psrc->oneSNFluenceDet(Ev, type, MH);
        double xs  = peffLS->totalXS(Ev, type);
        totflu = flu*xs;
    }

    return totflu*Npar;
}

double SNdetect::getTSpectrum(double T, int type, int MH){
    double Evmin;
    double Npar;
    if(fchannel == NuP || fchannel == NuE){
        if(fchannel == NuP){
            double mp = 938;//MeV
            Evmin = TMath::Sqrt(mp*T/2.);
            Npar  = peffLS->getNumberOfProton();
        }
        if(fchannel == NuE){
            double me = 0.511;//MeV
            Evmin = 0.5*(T+TMath::Sqrt(T*T+2*T*me));
            Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon();
        }

        TF1 f("anaEv2T",fcnAnaEv2T, fEvmin, fEvmax, 3);
        f.SetParameter(0,T);
        f.SetParameter(1,type);
        f.SetParameter(2,MH);
        ROOT::Math::WrappedTF1 wf1(f);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(1e-3);
        double spect = Npar*ig.Integral(Evmin, fEvmax);
        return spect;
    }
    if(fchannel==IBD || fchannel==NCC || fchannel==BCC){
        return getXSweightedEvSpectrum(getEvFromT(T), type, MH);
    }
    if(fchannel == CNC){
        if((TMath::Abs(T-15.11))>1e-6) 
            return 0;
        else 
            return getEventAboveEthrVis(fEthr, type, MH);
    }


    return 0;
}

double SNdetect::getEvisSpectrum(double Evis, int type, int MH){
    if(fchannel == NuP){
        double T = grquenchInv->Eval(Evis);
        double dTqdTp = (grquench->Eval(T+0.005)-grquench->Eval(T))/0.005; 
        return getTSpectrum(T, type, MH)/dTqdTp;
    }
    else{
        return getTSpectrum(getTFromEvis(Evis), type,MH);
    }

    return 0;
}

double SNdetect::getEobsSpectrum(double Eobs, int type, int MH){
    if(fchannel==CNC){
        double prob = peffLS->getProbResGauss(Eobs, 15.11);
        return getEventAboveEthrVis(fEthr, type, MH)*prob;
    }
    else{
        TF1 f("anaEvisEobs", fcnAnaT2Eobs, fTmin, fTmax, 3);
        f.SetParameter(0, Eobs);
        f.SetParameter(1, type);
        f.SetParameter(2, MH);

        ROOT::Math::WrappedTF1 wf1(f);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(1e-2);
        double integ = ig.Integral(fTmin, fTmax);

        //int np = 500;
        //double *x = new double[np];
        //double *w = new double[np];
        //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-8);
        //double integ = f.IntegralFast(np,x,w,fTmin,fTmax);
        //delete[] x;
        //delete[] w;
        return integ;
    }
}

double SNdetect::getEventAboveEthrVis(double Ethr, int type, int MH){
    double fEvmax = psrc->getSNmaxEv();
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
                break;
            }
        case NuE:
            {
                double me = 0.511;
                fEvismax = fEvmax/(1+me/(2.*fEvmax));
                break;
            }
        case IBD:
            {
                fEvismax = fEvmax-0.8;
                break;
            }
        case NCC:
            {
                fEvismax = fEvmax-16.827;
                break;
            }
        case BCC:
            {
                fEvismax = fEvmax-13.88+0.511;
                break;
            }
        case CNC:
            {
                break;
            }
    }

    fEvismax = 80;
    if(fchannel==CNC){
        //for C12-NC, use Ev to estimate the total number of events
        TF1 f("eveVis", fcnEventsVisCNC, 0, fEvmax, 2);
        f.SetParameter(0, type);
        f.SetParameter(1, MH);
        int np = 700;
        double *x = new double[np];
        double *w = new double[np];
        f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
        double integ = f.IntegralFast(np,x,w,15.11,fEvmax);
        delete[] x;
        delete[] w;
        return integ;
    }

    else{
        TF1 f("eveVis", fcnEventsVis, 0, fEvismax, 2);
        std::cout << "fEvismax = " << fEvismax << std::endl;
        f.SetParameter(0, type);
        f.SetParameter(1, MH);
        //ROOT::Math::WrappedTF1 wf1(f);
        //ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
        //ig.SetFunction(wf1);
        //ig.SetRelTolerance(1e-3);

        //double integ = ig.Integral(Ethr, fEvismax);

        //return integ;
        int np = 700;
        double *x = new double[np];
        double *w = new double[np];
        f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
        double integ = f.IntegralFast(np,x,w,Ethr,fEvismax);
        delete[] x;
        delete[] w;
        return integ;
    }
    return 0;
}

double SNdetect::getEventAboveEthrObs(double Ethr, int type, int MH){
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
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
                fEvismax = fTmax+0.511;
                break;
            }
        case BCC:
            {
                fEvismax = fTmax+0.511;
                break;
            }
        case NCC:
            {
                fEvismax = fTmax;
                break;
            }
        case CNC:
            {
                fEvismax = fTmax;
                break;
            }
    }


    TF1 f("eveObs", fcnEventsObs, 0, fEvismax, 2);
    f.SetParameter(0, type);
    f.SetParameter(1, MH);
    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
    ig.SetFunction(wf1);
    ig.SetRelTolerance(1e-3);
    double integ = ig.Integral(Ethr, fEvismax);
    return integ;

    //int np = 500;
    //double *x = new double[np];
    //double *w = new double[np];
    //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
    //double integ = f.IntegralFast(np,x,w,Ethr,fEvismax);
    //delete[] x;
    //delete[] w;
    //return integ;
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

double fcnAnaT2EobsatTime(double* x, double* par){
    double T    = x[0];
    double time = par[0];
    double Eobs = par[1];
    int    type = int(par[2]);
    int    MH   = int(par[3]);

    SNdetect* pdet = SNdetect::instance();
    double Evis = pdet->getEvisFromT(T);
    double prob = pdet->getPointerEffectLS()->getProbResGauss(Eobs, Evis);
    double flu  = pdet->getTSpectrumAtTime(time,T,type,MH);

    return prob*flu;
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
double fcnEventsVisatTimeCNC(double* x, double* par){
    double Ev   = x[0];
    double time = par[0];
    int type    = int(par[1]);
    int MH      = int(par[2]);

    SNdetect* pdet = SNdetect::instance();
    return pdet->getXSweightedEvSpectrumAtTime(time, Ev, type, MH);
}
double fcnLumiVisatTime(double* x, double* par){
    double Evis = x[0];
    double time = par[0];
    int    type = int(par[1]);
    int      MH = int(par[2]);

    SNdetect* pdet = SNdetect::instance();
    return Evis*pdet->getEvisSpectrumAtTime(time, Evis, type, MH);
}

double fcnEventsObsatTime(double* x, double* par){
    double Eobs = x[0];
    double time = par[0];
    int    type = int(par[1]);
    int      MH = int(par[2]);

    //std::cout << "num test" << std::endl;
    SNdetect* pdet = SNdetect::instance();
    return pdet->getEobsSpectrumAtTime(time, Eobs, type, MH);
}

double fcnLumiObsatTime(double* x, double* par){
    double Eobs = x[0];
    double time = par[0];
    int    type = int(par[1]);
    int      MH = int(par[2]);

    std::cout << "lumi test" << std::endl;
    SNdetect* pdet = SNdetect::instance();
    return Eobs*pdet->getEobsSpectrumAtTime(time, Eobs, type, MH);
}
//======FCN


double SNdetect::getXSweightedEvSpectrumAtTime(double time, double Ev, int type, int MH){
    double Npar;
    //number of protons
    if(fchannel == NuP || fchannel == IBD){
        Npar = peffLS->getNumberOfProton();
    }   
    //number of electrons
    if(fchannel == NuE){
        Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon();
    }
    //number of carbon-12, 98.9% of carbons
    if(fchannel == NCC || fchannel == BCC || fchannel == CNC){
        Npar = 0.989*peffLS->getNumberOfCarbon();
    }
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
        totflu = flu*xs;
    }
    return totflu*Npar;
}

double SNdetect::getTSpectrumAtTime(double time, double T, int type, int MH){
    double Evmin;
    double Npar;
    if(fchannel == NuP || fchannel == NuE){
        if(fchannel == NuP){
            double mp = 938;//MeV
            Evmin = TMath::Sqrt(mp*T/2.);
            Npar  = peffLS->getNumberOfProton();
        }
        if(fchannel == NuE){
            double me = 0.51;//MeV
            Evmin = 0.5*(T+TMath::Sqrt(T*T+2*T*me));
            Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon();
        }

        // test line
        // std::cout << "Evmin, fEvmax, Npar: " << Evmin << ", " << fEvmax << ", " << Npar << std::endl;
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
    if(fchannel==IBD || fchannel==NCC || fchannel==BCC){
        return getXSweightedEvSpectrumAtTime(time,getEvFromT(T), type, MH);
    }
    if(fchannel == CNC){
        if(T!=15.11) 
            return 0;
        else 
            return getEventAboveEthrVisAtTime(time, 0.2, type, MH);
    }
    return 0;
}

double SNdetect::getEvisSpectrumAtTime(double time, double Evis, int type, int MH){
    if(fchannel==NuP){
        double T = grquenchInv->Eval(Evis);
        double dTqdTp = (grquench->Eval(T+0.005)-grquench->Eval(T))/0.005; 
        return getTSpectrumAtTime(time, getTFromEvis(Evis),type, MH)/dTqdTp;
    }
    else{
        return getTSpectrumAtTime(time, getTFromEvis(Evis), type, MH);
    }

    return 0;
}

double SNdetect::getEobsSpectrumAtTime(double time, double Eobs, int type, int MH){
    TF1 f("anaEvisEobsatTime", fcnAnaT2EobsatTime, fTmin, fTmax, 4);
    f.SetParameter(0, time);
    f.SetParameter(1, Eobs);
    f.SetParameter(2, type);
    f.SetParameter(3, MH);

    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
    ig.SetFunction(wf1);
    ig.SetRelTolerance(1e-3);
    double integ = ig.Integral(fTmin,fTmax);

    //int np = 500;
    //double *x = new double[np];
    //double *w = new double[np];
    //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
    //double integ = f.IntegralFast(np,x,w,fTmin,fTmax);
    //delete[] x;
    //delete[] w;

    return integ;

}

double SNdetect::getEventAboveEthrVisAtTime(double time, double Ethr, int type, int MH){
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
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
        case NCC:
            {
                fEvismax = fEvmax-16.827;
                break;
            }
        case BCC:
            {
                fEvismax = fEvmax-13.88+0.511;
                break;
            }
        case CNC:
            {
                break;
            }
    }

    if(fchannel==CNC){
        //for C12-NC, use Ev to estimate the total number of events
        TF1 f("eveVis", fcnEventsVisatTimeCNC, 0, fEvmax, 3);
        f.SetParameter(0, time);
        f.SetParameter(1, type);
        f.SetParameter(2, MH);
        int np = 500;
        double *x = new double[np];
        double *w = new double[np];
        f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
        double integ = f.IntegralFast(np,x,w,15.11,fEvmax);
        delete[] x;
        delete[] w;
        return integ;

    }
    else{
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
    }
    return 0;
}

double SNdetect::getEventAboveEthrObsAtTime(double time, double Ethr, int type, int MH){
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
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
        case BCC:
            {
                fEvismax = fTmax+0.511;
                break;
            }
        case NCC:
            {
                fEvismax = fTmax;
                break;
            }
        case CNC:
            {
                fEvismax = fTmax;
                break;
            }
    }

    TF1 f("eveObsatTime", fcnEventsObsatTime, 0, fEvismax, 3);
    f.SetParameter(0, time);
    f.SetParameter(1, type);
    f.SetParameter(2, MH);

    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
    ig.SetFunction(wf1);
    ig.SetRelTolerance(1e-3);
    double integ = ig.Integral(Ethr, fEvismax);
    return integ;


    //int np = 500;
    //double *x = new double[np];
    //double *w = new double[np];
    //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
    //double integ = f.IntegralFast(np,x,w,Ethr,fEvismax);
    //delete[] x;
    //delete[] w;
    //return integ;
}
double SNdetect::getLumiAboveEthrVisAtTime(double time, double Ethr, int type, int MH){
    double fEvmax = psrc->getSNmaxEv();
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
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
        case BCC:
            {
                fEvismax = fTmax+0.511;
                break;
            }
        case NCC:
            {
                fEvismax = fTmax;
                break;
            }
        case CNC:
            {
                fEvismax = fTmax;
                break;
            }
    }

    TF1 f("eveVisatTime", fcnLumiVisatTime, 0, fEvismax, 3);
    f.SetParameter(0, time);
    f.SetParameter(1, type);
    f.SetParameter(2, MH);
    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
    ig.SetFunction(wf1);
    ig.SetRelTolerance(1e-3);

    double integ = ig.Integral(Ethr, fEvismax);

    //int np = 1000;
    //double *x = new double[np];
    //double *w = new double[np];
    //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
    //double integ = f.IntegralFast(np,x,w,Ethr,fEvismax);
    //delete[] x;
    //delete[] w;
    return integ;
}

double SNdetect::getLumiAboveEthrObsAtTime(double time, double Ethr, int type,int MH){
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
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
        case BCC:
            {
                fEvismax = fTmax+0.511;
                break;
            }
        case NCC:
            {
                fEvismax = fTmax;
                break;
            }
        case CNC:
            {
                fEvismax = fTmax;
                break;
            }
    }

    TF1 f("eveObsatTime", fcnLumiObsatTime, 0, fEvismax, 3);
    f.SetParameter(0, time);
    f.SetParameter(1, type);
    f.SetParameter(2, MH);

    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
    ig.SetFunction(wf1);
    ig.SetRelTolerance(1e-3);
    double integ = ig.Integral(Ethr, fEvismax);
    return integ;

    //int np = 500;
    //double *x = new double[np];
    //double *w = new double[np];
    //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
    //double integ = f.IntegralFast(np,x,w,Ethr,fEvismax);
    //delete[] x;
    //delete[] w;
    //return integ;
}

//==================Time Range================//
//
//FCN
double fcnAnaEv2TTimeInterval(double* x, double* par){
    double Ev = x[0];
    double tmin = par[0];
    double tmax = par[1];
    double T    = par[2];
    int    type = int(par[3]);
    int    MH   = int(par[4]);

    SNdetect* pdet = SNdetect::instance();
    int ntype = 6;
    double totflu = 0;
    if(type==-1){
        for(int it=0; it<ntype; it++){
            double flu = pdet->getPointerSrc()->oneSNFluenceDetTimeInterval(Ev, tmin, tmax, it, MH);
            double dxs = pdet->getPointerEffectLS()->differentialXS(Ev, T, it);
            totflu += flu*dxs;
        }
    }
    else{
        double flu = pdet->getPointerSrc()->oneSNFluenceDetTimeInterval(Ev, tmin, tmax, type, MH);
        double dxs = pdet->getPointerEffectLS()->differentialXS(Ev, T, type);
        totflu = flu*dxs;
    }

    return totflu;
}

double fcnAnaT2EobsTimeInterval(double* x, double* par){
    double T    = x[0];
    double tmin = par[0];
    double tmax = par[1];
    double Eobs = par[2];
    int    type = int(par[3]);
    int    MH   = int(par[4]);

    SNdetect* pdet = SNdetect::instance();
    double Evis = pdet->getEvisFromT(T);
    double prob = pdet->getPointerEffectLS()->getProbResGauss(Eobs, Evis);
    double flu  = pdet->getTSpectrumTimeInterval(tmin, tmax, T,type,MH);

    return prob*flu;
}
double fcnEventsVisTimeInterval(double* x, double* par){
    double Evis = x[0];
    double tmin = par[0];
    double tmax = par[1];
    int    type = int(par[2]);
    int    MH   = int(par[3]);

    SNdetect* pdet = SNdetect::instance();
    //std::cout << "NumEvis " << pdet->getEvisSpectrumTimeInterval(tmin, tmax, Evis, type, MH) << std::endl;
    return pdet->getEvisSpectrumTimeInterval(tmin, tmax, Evis, type, MH);
}
double fcnEventsVisTimeIntervalCNC(double* x, double* par){
    double Ev   = x[0];
    double tmin = par[0];
    double tmax = par[1];
    int type    = int(par[2]);
    int MH      = int(par[3]);

    SNdetect* pdet = SNdetect::instance();
    return pdet->getXSweightedEvSpectrumTimeInterval(tmin,tmax, Ev, type, MH);
}
double fcnLumiVisTimeInterval(double* x, double* par){
    double Evis = x[0];
    double time = par[0];
    int    type = int(par[1]);
    int      MH = int(par[2]);

    SNdetect* pdet = SNdetect::instance();
    return Evis*pdet->getEvisSpectrumAtTime(time, Evis, type, MH);
}

double fcnEventsObsTimeInterval(double* x, double* par){
    double Eobs = x[0];
    double tmin = par[0];
    double tmax = par[1];
    int    type = int(par[2]);
    int      MH = int(par[3]);

    //std::cout << "num test" << std::endl;
    SNdetect* pdet = SNdetect::instance();
    return pdet->getEobsSpectrumTimeInterval(tmin, tmax, Eobs, type, MH);
}

double fcnLumiObsTimeInterval(double* x, double* par){
    double Eobs = x[0];
    double time = par[0];
    int    type = int(par[1]);
    int      MH = int(par[2]);

    std::cout << "lumi test" << std::endl;
    SNdetect* pdet = SNdetect::instance();
    return Eobs*pdet->getEobsSpectrumAtTime(time, Eobs, type, MH);
}
//======FCN


double SNdetect::getXSweightedEvSpectrumTimeInterval(double tmin, double tmax, double Ev, int type, int MH){
    double Npar;
    //number of protons
    if(fchannel == NuP || fchannel == IBD){
        Npar = peffLS->getNumberOfProton();
    }   
    //number of electrons
    if(fchannel == NuE){
        Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon();
    }
    //number of carbon-12, 98.9% of carbons
    if(fchannel == NCC || fchannel == BCC || fchannel == CNC){
        Npar = 0.989*peffLS->getNumberOfCarbon();
    }
    int ntype = 6;
    double totflu = 0;
    if(type==-1){
        for(Int_t it=0; it<ntype; it++){
            double flu = psrc->oneSNFluenceDetTimeInterval(Ev, tmin, tmax, it, MH);
            double xs  = peffLS->totalXS(Ev, it);
            totflu += flu*xs;
        }
    }
    else{
        double flu = psrc->oneSNFluenceDetTimeInterval(Ev, tmin, tmax, type, MH);
        double xs  = peffLS->totalXS(Ev, type);
        totflu = flu*xs;
    }
    return totflu*Npar;
}

double SNdetect::getTSpectrumTimeInterval(double tmin, double tmax, double T, int type, int MH){
    double Evmin;
    double Npar;
    if(fchannel == NuP || fchannel == NuE){
        if(fchannel == NuP){
            double mp = 938;//MeV
            Evmin = TMath::Sqrt(mp*T/2.);
            Npar  = peffLS->getNumberOfProton();
        }
        if(fchannel == NuE){
            double me = 0.51;//MeV
            Evmin = 0.5*(T+TMath::Sqrt(T*T+2*T*me));
            Npar = peffLS->getNumberOfProton()+6*peffLS->getNumberOfCarbon();
        }

        TF1 f("anaEv2TTimeInterval",fcnAnaEv2TTimeInterval, fEvmin, fEvmax, 5);
        f.SetParameter(0, tmin);
        f.SetParameter(1, tmax);
        f.SetParameter(2, T);
        f.SetParameter(3, type);
        f.SetParameter(4, MH);
        ROOT::Math::WrappedTF1 wf1(f);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(1e-3);
        double spect = Npar*ig.Integral(Evmin, fEvmax);
        return spect;
    }
    if(fchannel==IBD || fchannel==NCC || fchannel==BCC){
        return getXSweightedEvSpectrumTimeInterval(tmin, tmax, getEvFromT(T), type, MH);
    }
    if(fchannel == CNC){
        if(T!=15.11) 
            return 0;
        else 
            return getEventAboveEthrVisTimeInterval(tmin, tmax, 0.2, type, MH);
    }
    return 0;
}

double SNdetect::getEvisSpectrumTimeInterval(double tmin, double tmax, double Evis, int type, int MH){
    if(fchannel==NuP){
        double T = grquenchInv->Eval(Evis);
        double dTqdTp = (grquench->Eval(T+0.005)-grquench->Eval(T))/0.005; 
        return getTSpectrumTimeInterval(tmin, tmax, getTFromEvis(Evis),type, MH)/dTqdTp;
    }
    else{
        return getTSpectrumTimeInterval(tmin, tmax, getTFromEvis(Evis), type, MH);
    }

    return 0;
}

double SNdetect::getEobsSpectrumTimeInterval(double tmin, double tmax, double Eobs, int type, int MH){
    TF1 f("anaEvisEobsTimeInterval", fcnAnaT2EobsTimeInterval, fTmin, fTmax, 5);
    f.SetParameter(0, tmin);
    f.SetParameter(1, tmax);
    f.SetParameter(2, Eobs);
    f.SetParameter(3, type);
    f.SetParameter(4, MH);

    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
    ig.SetFunction(wf1);
    ig.SetRelTolerance(1e-3);
    double integ = ig.Integral(fTmin,fTmax);

    //int np = 500;
    //double *x = new double[np];
    //double *w = new double[np];
    //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
    //double integ = f.IntegralFast(np,x,w,fTmin,fTmax);
    //delete[] x;
    //delete[] w;

    return integ;

}

double SNdetect::getEventAboveEthrVisTimeInterval(double tmin, double tmax, double Ethr, int type, int MH){
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
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
        case NCC:
            {
                fEvismax = fEvmax-16.827;
                break;
            }
        case BCC:
            {
                fEvismax = fEvmax-13.88+0.511;
                break;
            }
        case CNC:
            {
                break;
            }
    }

    if(fchannel==CNC){
        //for C12-NC, use Ev to estimate the total number of events
        TF1 f("eveVis", fcnEventsVisTimeIntervalCNC, 0, fEvmax, 4);
        f.SetParameter(0, tmin);
        f.SetParameter(1, tmax);
        f.SetParameter(2, type);
        f.SetParameter(3, MH);
        int np = 500;
        double *x = new double[np];
        double *w = new double[np];
        f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
        double integ = f.IntegralFast(np,x,w,15.11,fEvmax);
        delete[] x;
        delete[] w;
        return integ;

    }
    else{
        TF1 f("eveVisTimeInterval", fcnEventsVisTimeInterval, 0, fEvismax, 4);
        f.SetParameter(0, tmin);
        f.SetParameter(1, tmax);
        f.SetParameter(2, type);
        f.SetParameter(3, MH);
        ROOT::Math::WrappedTF1 wf1(f);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(1e-3);

        double integ = ig.Integral(Ethr, fEvismax);

        //int np = 1000;
        //double *x = new double[np];
        //double *w = new double[np];
        //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
        //double integ = f.IntegralFast(np,x,w,Ethr,fEvismax);
        //delete[] x;
        //delete[] w;
        return integ;
    }
    return 0;
}

double SNdetect::getEventAboveEthrObsTimeInterval(double tmin, double tmax, double Ethr, int type, int MH){
    double fEvismax;
    switch(fchannel){
        case NuP:
            {
                double mp = 938.;
                fEvismax = grquench->Eval(2.*fEvmax*fEvmax/mp);
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
        case BCC:
            {
                fEvismax = fTmax+0.511;
                break;
            }
        case NCC:
            {
                fEvismax = fTmax;
                break;
            }
        case CNC:
            {
                fEvismax = fTmax;
                break;
            }
    }

    TF1 f("eveObsTimeInterval", fcnEventsObsTimeInterval, 0, fEvismax, 4);
    f.SetParameter(0, tmin);
    f.SetParameter(1, tmax);
    f.SetParameter(2, type);
    f.SetParameter(3, MH);

    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,ROOT::Math::Integration::kGAUSS61);
    ig.SetFunction(wf1);
    ig.SetRelTolerance(1e-3);
    double integ = ig.Integral(Ethr, fEvismax);
    return integ;


    //int np = 500;
    //double *x = new double[np];
    //double *w = new double[np];
    //f.CalcGaussLegendreSamplingPoints(np,x,w,1e-10);
    //double integ = f.IntegralFast(np,x,w,Ethr,fEvismax);
    //delete[] x;
    //delete[] w;
    //return integ;
}

//
//=====]



double SNdetect::getEvisFromT(double T){
    switch(fchannel){
        case(NuP):
            return grquench->Eval(T);
        case(NuE):
            return T;
        case(IBD):
            return T+0.51;
        case(BCC):
            return T+0.51;
        case(NCC):
            return T;
        case(CNC):
            return T;
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
        case(BCC):
            return Evis-0.51;
        case(NCC):
            return Evis;
        case(CNC):
            return Evis;
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
        case(BCC):
            return Ev-13.88;
        case(NCC):
            return Ev-16.827;
        case(CNC):
            {
                std::cout << "CNC: T is a fixed value, but not related with Ev." << std::endl;
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
        case(BCC):
            {
                return T+13.88;
                break;
            }
        case(NCC):
            {
                return T+16.827;
                break;
            }
        case(CNC):
            {
                std::cout << "CNC: T is a fixed value, but not related with Ev." << std::endl;
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
