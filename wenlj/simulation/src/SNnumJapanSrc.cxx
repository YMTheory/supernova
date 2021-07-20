#include <iostream>
#include <fstream>
#include <vector>
#include <array>

#include "TVectorD.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TSpline.h"
#include "TF1.h"

#include "SNsource.hh"
#include "SNnumJapanSrc.hh"
#include "SNJapanIntegFcn.hh"

SNnumJapanSrc::SNnumJapanSrc() : SNsource(){
    pjapInteg = NULL;
    f1TE = NULL;
}

SNnumJapanSrc::SNnumJapanSrc(int imode) : SNsource(){
    pjapInteg = new SNJapanIntegFcn(imode);
    f1TE = new TF1("f1TE",pjapInteg, &SNJapanIntegFcn::fluxTimeE,0,400,3,"SNJapanIntegFcn","fluxTimeE");
}

SNnumJapanSrc::~SNnumJapanSrc(){
    if(pjapInteg)delete pjapInteg;
    if(f1TE)delete f1TE;
}



double SNnumJapanSrc::oneSNFluenceDet(double E, int type){
    double index = 1./(4*TMath::Pi()*TMath::Power(dist*3.096e21,2));//cm^-2
    double flux = pjapInteg->getFluxE(E, type);
    return flux*index;
}

double SNnumJapanSrc::oneSNFluenceDet(double E, int type, int MH){
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

double SNnumJapanSrc::totalSNFluenceDet(double E){
    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDet(E,it);
    }

    return totflu;
}

double SNnumJapanSrc::totalSNFluenceDet(double E, int MH){
    int ntype = 6;
    double tfluence=0;
    for(int ii=0; ii<ntype; ii++){
        tfluence += oneSNFluenceDet(E,ii, MH);
    }
    return tfluence;
}


double SNnumJapanSrc::oneSNFluenceDetTimeIntegE(double time, int type, int MH){
    double index =1./(4*TMath::Pi()*TMath::Power(3.086e21,2))/(dist*dist);
    double fluence = 0;
    if(MH == 0){
        fluence = pjapInteg->getNumT(time, type);
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
        if(type == 0) fluence = p*pjapInteg->getNumT(time,0) + (1-p)*pjapInteg->getNumT(time,2);
        if(type == 1) fluence = pbar*pjapInteg->getNumT(time,1) + (1-pbar)*pjapInteg->getNumT(time,3);
        if(type == 2 || type == 4) fluence = 0.5*(1-p)*pjapInteg->getNumT(time,0) + 0.5*(1+p)*pjapInteg->getNumT(time,2);
        if(type == 3 || type == 5) fluence = 0.5*(1-pbar)*pjapInteg->getNumT(time,1) + 0.5*(1+pbar)*pjapInteg->getNumT(time,3);

    }
    return fluence*index;
}

double SNnumJapanSrc::totalSNFluenceDetTimeIntegE(double time, int MH){
    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetTimeIntegE(time,it, MH);
    }

    return totflu;
}

double SNnumJapanSrc::oneSNFluenceDetAtTime(double time, double E, int type, int MH){
    double index =1./(4*TMath::Pi()*TMath::Power(3.086e21,2))/(dist*dist);
    double fluence = 0;
    if(MH == 0){
        fluence = pjapInteg->getEventAtTime(time, E, type);
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
        if(type == 0) fluence = p*pjapInteg->getEventAtTime(time,E, 0) + (1-p)*pjapInteg->getEventAtTime(time, E, 2);
        if(type == 1) fluence = pbar*pjapInteg->getEventAtTime(time, E, 1) + (1-pbar)*pjapInteg->getEventAtTime(time, E, 3);
        if(type == 2 || type == 4) fluence = 0.5*(1-p)*pjapInteg->getEventAtTime(time, E, 0) + 0.5*(1+p)*pjapInteg->getEventAtTime(time, E, 2);
        if(type == 3 || type == 5) fluence = 0.5*(1-pbar)*pjapInteg->getEventAtTime(time, E, 1) + 0.5*(1+pbar)*pjapInteg->getEventAtTime(time, E, 3);

    }
    return fluence*index;
}

double SNnumJapanSrc::totalSNFluenceDetAtTime(double time, double E, int MH){
    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetAtTime(time, E, it, MH);
    }
    return totflu;

}

//===>Time Interval
double SNnumJapanSrc::oneSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int type){
    double index = 1./(4*TMath::Pi()*TMath::Power(dist*3.096e21,2));//cm^-2
    f1TE->SetParameters(tfirst,tlast,type);
    double fluence = f1TE->Eval(E);
    return fluence*index;
}
double SNnumJapanSrc::totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast){

    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetTimeInterval(E, tfirst, tlast, it);
    }

    return totflu;
}

double SNnumJapanSrc::oneSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int type, int MH){
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

double SNnumJapanSrc::totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int MH){

    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetTimeInterval(E, tfirst, tlast, it, MH);
    }

    return totflu;
}
double SNnumJapanSrc::oneSNLuminosityTime(double time, int type){
    double lumi = pjapInteg->getLumT(time, type);
    return lumi;
}

double SNnumJapanSrc::oneSNAverageETime(double time, int type){
    return pjapInteg->getAverageET(time, type);
}


void SNnumJapanSrc::getTimeRange(double& tmin, double& tmax, int type){
    pjapInteg->getTimeLimits(tmin, tmax, type);
}
