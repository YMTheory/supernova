#include "SNanaSrc.hh"
#include "TMath.h"


SNanaSrc::SNanaSrc():SNsource(){
}

SNanaSrc::~SNanaSrc(){

}

double SNanaSrc::oneSNFluenceDet(double E, int type){
    double gammaA[6];
    int ntype = 6;
    double Eatmp[6] = {0};
    for(int it=0; it<ntype; it++){
        Eatmp[it] = totE/6/1e52;
        gammaA[it] = 3;

    }
    double factor = 6.921e14/(4*TMath::Pi()*dist*dist);
    double ifluence = (Eatmp[type]/averagE[type])*(TMath::Power(1+gammaA[type],1+gammaA[type])/TMath::Gamma(1+gammaA[type]))*(TMath::Power(E,gammaA[type])/TMath::Power(averagE[type],1+gammaA[type]))*TMath::Exp(-(1+gammaA[type])*E/averagE[type]);

    return ifluence*factor;
}

//mass hierarchy: 0 no osc; 
//                1 NH; 2 IH;(MSW effect)
double SNanaSrc::oneSNFluenceDet(double E, int type, int MH){
    double fluence = 0;
    double sin2theta12 = 0.297;
    double cos2theta12 = 1-sin2theta12;
    if(MH==0)return oneSNFluenceDet(E, type);
    else{
        if(MH==1){
            if(type==0) fluence = oneSNFluenceDet(E,2);
            if(type==1) fluence = cos2theta12*oneSNFluenceDet(E,1)+sin2theta12*oneSNFluenceDet(E,2);
            if(type>=2) fluence = 0.25*oneSNFluenceDet(E,0)+0.25*sin2theta12*oneSNFluenceDet(E,1)+0.25*(2+cos2theta12)*oneSNFluenceDet(E,2);
        }
        if(MH==2){
            if(type==0) fluence = sin2theta12*oneSNFluenceDet(E,0)+cos2theta12*oneSNFluenceDet(E,2);
            if(type==1) fluence = oneSNFluenceDet(E,2);
            if(type>=2) fluence = 0.25*cos2theta12*oneSNFluenceDet(E,0)+0.25*oneSNFluenceDet(E,1)+0.25*(2+sin2theta12)*oneSNFluenceDet(E,2);
        }
    }
    return fluence;
}

double SNanaSrc::totalSNFluenceDet(double E){
    //equipartitioned in all 6 neutrino flavors
    int ntype = 6;
    double tfluence = 0;
    for(int ii=0; ii<ntype; ii++){
        tfluence += oneSNFluenceDet(E,ii);
    }

    return tfluence;
}
double SNanaSrc::totalSNFluenceDet(double E, int MH){
    int ntype = 6;
    double tfluence=0;
    for(int ii=0; ii<ntype; ii++){
        tfluence += oneSNFluenceDet(E,ii, MH);
    }
    return tfluence;

}

double SNanaSrc::oneSNFluenceDetTimeIntegE(double time, int type, int MH){
    return -1;
}
double SNanaSrc::totalSNFluenceDetTimeIntegE(double time, int MH){
    return -1;
}
double SNanaSrc::oneSNFluenceDetAtTime(double time, double E, int type, int MH){
    return -1;
}
double SNanaSrc::totalSNFluenceDetAtTime(double time, double E, int MH){
    return -1;
}

double SNanaSrc::oneSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int type){
    return -1;
}
double SNanaSrc::totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast){
    return -1;
}
double SNanaSrc::oneSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int type, int MH){
    return -1;
}
double SNanaSrc::totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int MH){
    return -1;
}

double SNanaSrc::oneSNLuminosityTime(double time, int type){
    return -1;
}
double SNanaSrc::oneSNAverageETime(double time, int type){
    return -1;
}

void SNanaSrc::getTimeRange(double& tmin, double& tmax, int type){
}
