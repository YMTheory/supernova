#include "SNnueLS.hh"
#include "TMath.h"


SNnueLS::SNnueLS():SNeffectLS(){
}

SNnueLS::~SNnueLS(){

}

double SNnueLS::differentialXS(double E, double T, int type){
    //type: 0,1->nu_e,   anti_nu_e
    //      2,3->nu_mu,  anti_nu_mu
    //      4,5->nu_tau, anti_nu_tau

    double me = 0.51;//electron mass MeV
    //double Tmax = E/(1+me/2./E);//maximum T for E
    //if(T <0 || T>Tmax)return 0;
    if(E < (T/2+TMath::Sqrt(T*(T+me))/2))return 0;

    double diffxs;
    double sin2_thetaW = 0.23;
    double epsi_p;//epsilon+
    double epsi_m;//epsilon-
    double indexG = 1.732e-44;//2*G_u^2*m_e/pi
    if(type==2 || type==4){
        //nu_mu, nu_tau
        epsi_p = -1*sin2_thetaW;
        epsi_m = 0.5-sin2_thetaW;
        diffxs = indexG*(epsi_m*epsi_m+epsi_p*epsi_p*TMath::Power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E));
        return diffxs;
    }
    if(type==3 || type==5){
        //anti_nu_mu, anti_nu_tau
        epsi_p = -1*sin2_thetaW;
        epsi_m = 0.5-sin2_thetaW;
        diffxs = indexG*(epsi_p*epsi_p+epsi_m*epsi_m*TMath::Power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E));
        return diffxs;
    }
    if(type == 0){
        //nu_e
        epsi_p = -1*sin2_thetaW;
        epsi_m = -0.5-sin2_thetaW;
        diffxs = indexG*(epsi_m*epsi_m+epsi_p*epsi_p*TMath::Power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E));
        return diffxs;
    }
    if(type == 1){
        //anti_nu_e
        epsi_p = -1*sin2_thetaW;
        epsi_m = -0.5-sin2_thetaW;
        diffxs = indexG*(epsi_p*epsi_p+epsi_m*epsi_m*TMath::Power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E));

        return diffxs;
    }

    return 0;
}

double SNnueLS::totalXS(double E, int type){

    //type: 0,1->nu_e,   anti_nu_e
    //      2,3->nu_mu,  anti_nu_mu
    //      4,5->nu_tau, anti_nu_tau

    double totalxs;
    double sin2_thetaW = 0.23;
    double indexG = 1.732e-44;//2*G_u^2*m_e/pi
    if(type==2 || type==4){
        //nu_mu, nu_tau
        totalxs = (indexG/4.*E)*(1.-4.*sin2_thetaW+16./3.*TMath::Power(sin2_thetaW,2));
        return totalxs;
    }
    if(type==3 || type==5){
        //anti_nu_mu, anti_nu_tau
        totalxs = (indexG/4.*E)*(1./3.-4./3.*sin2_thetaW+16./3.*TMath::Power(sin2_thetaW,2));
        return totalxs;
    }
    if(type == 0){
        //nu_e
        totalxs = (indexG/4.*E)*(1.+4.*sin2_thetaW+16./3.*TMath::Power(sin2_thetaW,2));
        return totalxs;
    }
    if(type == 1){
        //anti_nu_e
        totalxs = (indexG/4.*E)*(1./3.+4./3.*sin2_thetaW+16./3.*TMath::Power(sin2_thetaW,2));
        return totalxs;

    }
    return 0;
}

