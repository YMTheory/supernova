#include <iostream>

#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"

#include "SNsource.hh"
#include "SNnumBurrowsSrc.hh"
#include "SNBurrowsIntegFcn.hh"

SNnumBurrowsSrc::SNnumBurrowsSrc() : SNsource(){
    for(int it=0; it<3; it++){
        grspectrum[it] = NULL;
    }
}
SNnumBurrowsSrc::SNnumBurrowsSrc(int imode) : SNsource(){
    //readSpectrum(imode);
    pgarfcn = new SNBurrowsIntegFcn(imode);
    fcnTimeInterval = new TF1("fcnTimeInterval", pgarfcn, &SNBurrowsIntegFcn::fcnfluxtime, -1, 12, 2, "SNBurrowsIntegFcn", "fcnfluxtime");
}

SNnumBurrowsSrc::~SNnumBurrowsSrc(){
}

void SNnumBurrowsSrc::readSpectrum(int imode){
    std::cout << "*** Start to read spectrum ***" << std::endl;
    TString path="/junofs/users/miaoyu/supernova/wenlj/simulation/data/Burrows/energySpeFiles";
    TFile* fin = TFile::Open(Form("%s/energySpecBurrows_%d.root",path.Data(),imode));
    if(!fin){
        printf("%s/energySpecBurrows_%d.root can't be opened!!\n",path.Data(), imode);
        exit(-1);
    }
    for(int it=0; it<3; it++){
        grspectrum[it] = (TGraph*)fin->Get(Form("grmod%dtype%d",imode, it));
    }
    fin->Close();
}


double SNnumBurrowsSrc::oneSNFluenceDet(double E, int type){
    double index =1./(4*TMath::Pi()*TMath::Power(3.086e21,2))/(dist*dist);
    if(type < 2){
        //return index*grspectrum[type]->Eval(E,0,"S");
        return index*grspectrum[type]->Eval(E);
    }
    if(type >= 2){
        //return index*grspectrum[2]->Eval(E,0,"S");
        return index*grspectrum[2]->Eval(E);
    }
    return 0;
}

double SNnumBurrowsSrc::totalSNFluenceDet(double E){
    int ntype = 6;
    double tfluence = 0;
    for(int ii=0; ii<ntype; ii++){
        tfluence += oneSNFluenceDet(E,ii);
    }
    return tfluence;
}


double SNnumBurrowsSrc::oneSNFluenceDet(double E, int type, int MH){

    double fluence = 0;
    if(MH==0)return oneSNFluenceDet(E, type);
    else{
        double p=0;
        double pbar=0;
        if(MH==1){
            p    = 0.022;//sin^2theta13
            pbar = 0.687;//cos^2theta12cos^2theta13
        }
        if(MH==2){
            p    = 0.291;//sin^2theta12cos^2theta13
            pbar = 0.022;//sin^2theta13 
        }
        if(type==0) fluence = p*oneSNFluenceDet(E,0)+(1-p)*oneSNFluenceDet(E,2);//nue
        if(type==1) fluence = pbar*oneSNFluenceDet(E,1)+(1-pbar)*oneSNFluenceDet(E,3);//nue_bar
        if(type==2 || type==4) fluence = 0.5*(1-p)*oneSNFluenceDet(E,0)+0.5*(1+p)*oneSNFluenceDet(E,2);
        if(type==3 || type==5) fluence = 0.5*(1-pbar)*oneSNFluenceDet(E,1)+0.5*(1+pbar)*oneSNFluenceDet(E,3);
    }

    return fluence;
}

double SNnumBurrowsSrc::totalSNFluenceDet(double E, int MH){
    int ntype = 6;
    double tfluence = 0;
    for(int ii=0; ii<ntype; ii++){
        tfluence += oneSNFluenceDet(E,ii, MH);
    }
    return tfluence;
}



double SNnumBurrowsSrc::totalSNFluenceDetTimeIntegE(double time, int MH)
{
    int ntype = 6;
    double num = 0;
    for(int ii=0; ii<ntype; ii++){
        num += oneSNFluenceDetTimeIntegE(time,ii, MH);
    }
    return num;
}
double SNnumBurrowsSrc::oneSNFluenceDetTimeIntegE(double time, int type, int MH) 
{
    double index =1./(4*TMath::Pi()*TMath::Power(3.086e21,2))/(dist*dist);
    double fluence = 0;
    if(MH == 0){
        fluence = pgarfcn->getNumT(time, type);
    }
    else{
        double p=0;
        double pbar=0;
        if(MH==1){
            p    = 0.022;//sin^2theta13
            pbar = 0.687;//cos^2theta12cos^2theta13
        }
        if(MH==2){
            p    = 0.291;//sin^2theta12cos^2theta13
            pbar = 0.022;//sin^2theta13 
        }
        if(type == 0) fluence = p*pgarfcn->getNumT(time,0) + (1-p)*pgarfcn->getNumT(time,2);
        if(type == 1) fluence = pbar*pgarfcn->getNumT(time,1) + (1-pbar)*pgarfcn->getNumT(time,3);
        if(type == 2 || type == 4) fluence = 0.5*(1-p)*pgarfcn->getNumT(time,0) + 0.5*(1+p)*pgarfcn->getNumT(time,2);
        if(type == 3 || type == 5) fluence = 0.5*(1-pbar)*pgarfcn->getNumT(time,1) + 0.5*(1+pbar)*pgarfcn->getNumT(time,3);

    }
    return fluence*index;
}

double SNnumBurrowsSrc::totalSNFluenceDetAtTime(double time, double E, int MH)
{
    int ntype = 6;
    double flux = 0;
    for(int ii=0; ii<ntype; ii++){
        flux += oneSNFluenceDetAtTime(time,E,ii, MH);
    }
    return flux;
}

double SNnumBurrowsSrc::oneSNFluenceDetAtTime(double time, double E, int type, int MH)
{
    double index =1./(4*TMath::Pi()*TMath::Power(3.086e21,2))/(dist*dist);
    double fluence = 0;
    if(MH == 0){
        fluence = pgarfcn->getEventAtTime(time, E, type);
    }
    else{
        double p=0;
        double pbar=0;
        if(MH==1){
            p    = 0.022;//sin^2theta13
            pbar = 0.687;//cos^2theta12cos^2theta13
        }
        if(MH==2){
            p    = 0.291;//sin^2theta12cos^2theta13
            pbar = 0.022;//sin^2theta13 
        }
        if(type == 0) fluence = p*pgarfcn->getEventAtTime(time,E, 0) + (1-p)*pgarfcn->getEventAtTime(time, E, 2);
        if(type == 1) fluence = pbar*pgarfcn->getEventAtTime(time, E, 1) + (1-pbar)*pgarfcn->getEventAtTime(time, E, 3);
        if(type == 2 || type == 4) fluence = 0.5*(1-p)*pgarfcn->getEventAtTime(time, E, 0) + 0.5*(1+p)*pgarfcn->getEventAtTime(time, E, 2);
        if(type == 3 || type == 5) fluence = 0.5*(1-pbar)*pgarfcn->getEventAtTime(time, E, 1) + 0.5*(1+pbar)*pgarfcn->getEventAtTime(time, E, 3);

    }
    return fluence*index;
}

double SNnumBurrowsSrc::snFluenceDetAtTime(double& time, double nuMass, double E, int type, int MH)
{
    double index =1./(4*TMath::Pi()*TMath::Power(3.086e21,2))/(dist*dist);
    double fluence = 0;
    if(MH == 0){
        fluence = pgarfcn->getEventAtTime(time, E, type);
    }
    else{
        double p=0;
        double pbar=0;
        if(MH==1){
            p    = 0.022;//sin^2theta13
            pbar = 0.687;//cos^2theta12cos^2theta13
        }
        if(MH==2){
            p    = 0.291;//sin^2theta12cos^2theta13
            pbar = 0.022;//sin^2theta13 
        }
        if(type == 0) fluence = p*pgarfcn->getEventAtTime(time,E, 0) + (1-p)*pgarfcn->getEventAtTime(time, E, 2);
        if(type == 1) fluence = pbar*pgarfcn->getEventAtTime(time, E, 1) + (1-pbar)*pgarfcn->getEventAtTime(time, E, 3);
        if(type == 2 || type == 4) fluence = 0.5*(1-p)*pgarfcn->getEventAtTime(time, E, 0) + 0.5*(1+p)*pgarfcn->getEventAtTime(time, E, 2);
        if(type == 3 || type == 5) fluence = 0.5*(1-pbar)*pgarfcn->getEventAtTime(time, E, 1) + 0.5*(1+pbar)*pgarfcn->getEventAtTime(time, E, 3);
    }
    // modify the neutrino arrival time
    double DeltaT = 5.14e-3 * (nuMass*nuMass) * (100.0/E/E) * (dist/10); // unit: s
    time = time + DeltaT;

    return fluence*index;
}


#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
double SNnumBurrowsSrc::oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type)
{
    double num = 0;
    double index = 1./(4*TMath::Pi()*TMath::Power(3.086e21,2))/(dist*dist);
    int nbin_Ev = 100;
    double step_Ev = (fEvmax-fEvmin)/nbin_Ev;
    double* vecEv = new double[nbin_Ev];
    double* flux  = new double[nbin_Ev];
    for(int ib=0; ib<nbin_Ev; ib++){
        vecEv[ib] = fEvmin+step_Ev*(ib+0.5);
        fcnTimeInterval->SetParameter(0, vecEv[ib]);
        fcnTimeInterval->SetParameter(1, type);

        ROOT::Math::WrappedTF1 wf1(*fcnTimeInterval);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(1e-2);
        flux[ib] = ig.Integral(tfirst, tlast);
    }
    TGraph* grSpecTimeInterv = new TGraph(nbin_Ev, vecEv, flux);
    num = grSpecTimeInterv->Eval(E, 0, "S");
    delete[] vecEv;
    delete[] flux;
    delete grSpecTimeInterv;

    return num*index;
}
double SNnumBurrowsSrc::totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int MH)
{
    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetTimeInterval(E, tfirst, tlast, it, MH);
    }

    return totflu;
}


double SNnumBurrowsSrc::oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type, int MH)
{
    double fluence = 0;
    if(MH == 0){
        fluence = oneSNFluenceDetTimeInterval(E, tfirst, tlast, type);
    }
    else{
        double p=0;
        double pbar=0;
        if(MH==1){
            p    = 0.022;//sin^2theta13
            pbar = 0.687;//cos^2theta12cos^2theta13
        }
        if(MH==2){
            p    = 0.291;//sin^2theta12cos^2theta13
            pbar = 0.022;//sin^2theta13 
        }
        if(type == 0) fluence = p*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 0) + (1-p)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 2);
        if(type == 1) fluence = pbar*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 1) + (1-pbar)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 3);
        if(type == 2 || type == 4) fluence = 0.5*(1-p)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 0) + 0.5*(1+p)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 2);
        if(type == 3 || type == 5) fluence = 0.5*(1-pbar)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 1) + 0.5*(1+pbar)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 3);

    }
    return fluence;
}
double SNnumBurrowsSrc::totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast)
{
    int ntype = 6;
    double flux = 0;
    for(int ii=0; ii<ntype; ii++){
        flux += oneSNFluenceDetTimeInterval(E,tfirst,tlast,ii);
    }
    return flux;
}

double SNnumBurrowsSrc::oneSNLuminosityTime(double time, int type) {return pgarfcn->getLumT(time, type);}
double SNnumBurrowsSrc::oneSNAverageETime(double time, int type) {return pgarfcn->getAverageET(time, type);}

void SNnumBurrowsSrc::getTimeRange(double& tmin, double& tmax, int type)
{
    pgarfcn->getTimeLimits(tmin, tmax, type);
}



