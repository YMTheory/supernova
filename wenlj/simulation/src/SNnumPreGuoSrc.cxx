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
#include "SNnumPreGuoSrc.hh"
#include "SNpreGuoIntegFcn.hh"

SNnumPreGuoSrc::SNnumPreGuoSrc() : SNsource(){
    ppreInteg = NULL;
    f1TE = NULL;
    setSNmaxEv(40);
    setSNminEv(0);
}

SNnumPreGuoSrc::SNnumPreGuoSrc(int imode) : SNsource(){
    ppreInteg = new SNpreGuoIntegFcn(imode);
    f1TE = new TF1("f1TE",ppreInteg, &SNpreGuoIntegFcn::fluxTimeE,0,40,3,"SNpreGuoIntegFcn","fluxTimeE");
}

SNnumPreGuoSrc::~SNnumPreGuoSrc(){
    if(ppreInteg)delete ppreInteg;
    if(f1TE)delete f1TE;
}



double SNnumPreGuoSrc::oneSNFluenceDet(double E, int type){
    double index = 1./(4*TMath::Pi()*TMath::Power(dist*3.086e21,2));//cm^-2
    double flux = ppreInteg->getFluxE(E, type);
    return flux*index;
}

//PDG 2018
double SNnumPreGuoSrc::oneSNFluenceDet(double E, int type, int MH){
    double fluence = 0;
    double sin2theta12 = 0.297;
    double sin2theta13 = 0.022;
    double oscp, oscpbar;
    if(MH==0)return oneSNFluenceDet(E, type);
    if(MH==1){
        //NH  
        oscp = sin2theta13;
        oscpbar = (1-sin2theta12)*(1-sin2theta13);
    }
    if(MH==2){
        //IH
        oscp = sin2theta12*(1-sin2theta13);
        oscpbar = sin2theta13;
    }
    if(type==0) fluence = oscp*oneSNFluenceDet(E,0)+(1-oscp)*oneSNFluenceDet(E,2);
    if(type==1) fluence = oscpbar*oneSNFluenceDet(E,1)+(1-oscpbar)*oneSNFluenceDet(E,3);
    if(type==2 || type==4) fluence = 0.5*(1-oscp)*oneSNFluenceDet(E,0)+0.5*(1+oscp)*oneSNFluenceDet(E,2);
    if(type==3 || type==5) fluence = 0.5*(1-oscpbar)*oneSNFluenceDet(E,1)+0.5*(1+oscpbar)*oneSNFluenceDet(E,3);
    return fluence;
}

double SNnumPreGuoSrc::totalSNFluenceDet(double E){
    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDet(E,it);
    }

    return totflu;
}

double SNnumPreGuoSrc::totalSNFluenceDet(double E, int MH){
    int ntype = 6;
    double tfluence=0;
    for(int ii=0; ii<ntype; ii++){
        tfluence += oneSNFluenceDet(E,ii, MH);
    }
    return tfluence;
}


double SNnumPreGuoSrc::oneSNFluenceDetTimeIntegE(double time, int type, int MH){
    double index = 1./(4*TMath::Pi()*TMath::Power(dist*3.086e21,2));//cm^-2
    if(MH == 0){
        double flux = ppreInteg->getNumT(time, type);
        return index*flux;
    }
    double sin2theta12 = 0.297;
    double sin2theta13 = 0.022;
    double oscp, oscpbar;
    double fluence = 0;
    if(MH==1){
        //NH  
        oscp = sin2theta13;
        oscpbar = (1-sin2theta12)*(1-sin2theta13);
    }
    if(MH==2){
        //IH
        oscp = sin2theta12*(1-sin2theta13);
        oscpbar = sin2theta13;
    }
    if(type==0)fluence = oscp*ppreInteg->getNumT(time,0)+(1-oscp)*ppreInteg->getNumT(time,2);
    if(type==1)fluence = oscpbar*ppreInteg->getNumT(time,1)+(1-oscpbar)*ppreInteg->getNumT(time,3);
    if(type==2||type==4)fluence = 0.5*(1-oscp)*ppreInteg->getNumT(time,0)+0.5*(1+oscp)*ppreInteg->getNumT(time,2);
    if(type==3||type==5)fluence = 0.5*(1-oscpbar)*ppreInteg->getNumT(time,1)+0.5*(1+oscpbar)*ppreInteg->getNumT(time,3);

    return fluence*index;
}

double SNnumPreGuoSrc::totalSNFluenceDetTimeIntegE(double time, int MH){
    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetTimeIntegE(time,it, MH);
    }

    return totflu;
}
double SNnumPreGuoSrc::oneSNFluenceDetAtTime(double time, double E, int type, int MH){
    double index = 1./(4*TMath::Pi()*TMath::Power(dist*3.086e21,2));//cm^-2
    if(MH == 0){
        double flux = ppreInteg->getEventAtTime(time,E, type);
        return index*flux;
    }
    double sin2theta12 = 0.297;
    double sin2theta13 = 0.022;
    double oscp = 0;
    double oscpbar = 0;
    double fluence = 0;
    if(MH==1){
        //NH  
        oscp = sin2theta13;
        oscpbar = (1-sin2theta12)*(1-sin2theta13);
    }
    else if(MH==2){
        //IH
        oscp = sin2theta12*(1-sin2theta13);
        oscpbar = sin2theta13;
    }
    if(type==0)fluence = oscp*ppreInteg->getEventAtTime(time, E ,0)+(1-oscp)*ppreInteg->getEventAtTime(time, E, 2);
    if(type==1)fluence = oscpbar*ppreInteg->getEventAtTime(time,E,1)+(1-oscpbar)*ppreInteg->getEventAtTime(time, E, 3);
    if(type==2||type==4)fluence = 0.5*(1-oscp)*ppreInteg->getEventAtTime(time, E, 0)+0.5*(1+oscp)*ppreInteg->getEventAtTime(time, E, 2);
    if(type==3||type==5)fluence = 0.5*(1-oscpbar)*ppreInteg->getEventAtTime(time, E, 1)+0.5*(1+oscpbar)*ppreInteg->getEventAtTime(time, E, 3);

    return fluence*index;
}
double SNnumPreGuoSrc::totalSNFluenceDetAtTime(double time, double E, int MH){
    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetAtTime(time, E, it, MH);
    }
    return totflu;

}
double SNnumPreGuoSrc::oneSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int type){
    double index = 1./(4*TMath::Pi()*TMath::Power(dist*3.086e21,2));//cm^-2
    f1TE->SetParameters(tfirst,tlast,type);
    double flux = f1TE->Eval(E);
    return flux*index;
}
double SNnumPreGuoSrc::totalSNFluenceDetTimeInterval(double E, double tfiirst, double tlast){

    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetTimeInterval(E,tfiirst, tlast,it);
    }

    return totflu;
}
double SNnumPreGuoSrc::oneSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int type, int MH){
    double fluence = 0;
    double sin2theta12 = 0.297;
    double sin2theta13 = 0.022;
    double oscp, oscpbar;
    if(MH==0)return oneSNFluenceDetTimeInterval(E, tfirst, tlast, type);
    if(MH==1){
        //NH  
        oscp = sin2theta13;
        oscpbar = (1-sin2theta12)*(1-sin2theta13);
    }
    else if(MH==2){
        //IH
        oscp = sin2theta12*(1-sin2theta13);
        oscpbar = sin2theta13;
    }
    if(type==0) fluence = oscp*oneSNFluenceDetTimeInterval(E,tfirst, tlast, 0)+(1-oscp)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 2);
    if(type==1) fluence = oscpbar*oneSNFluenceDetTimeInterval(E,tfirst, tlast,1)+(1-oscpbar)*oneSNFluenceDetTimeInterval(E, tfirst, tlast,3);
    if(type==2 || type==4) fluence = 0.5*(1-oscp)*oneSNFluenceDetTimeInterval(E,tfirst, tlast,0)+0.5*(1+oscp)*oneSNFluenceDetTimeInterval(E,tfirst, tlast,2);
    if(type==3 || type==5) fluence = 0.5*(1-oscpbar)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 1)+0.5*(1+oscpbar)*oneSNFluenceDetTimeInterval(E, tfirst, tlast, 3);
    return fluence;
}
double SNnumPreGuoSrc::totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int MH){

    int ntype = 6;
    double totflu = 0;
    for(int it=0; it<ntype; it++){
        totflu += oneSNFluenceDetTimeInterval(E, tfirst, tlast, it, MH);
    }

    return totflu;
}
double SNnumPreGuoSrc::oneSNLuminosityTime(double time, int type){
    double lumi = ppreInteg->getLumT(time, type);
    return lumi;
}

double SNnumPreGuoSrc::oneSNAverageETime(double time, int type){
    return ppreInteg->getAverageET(time, type);
}


void SNnumPreGuoSrc::getTimeRange(double& tmin, double& tmax, int type){
    ppreInteg->getTimeLimits(tmin, tmax, type);
}
