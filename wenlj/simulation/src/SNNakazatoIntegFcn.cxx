#include <iostream>
#include <fstream>
#include <vector>

#include "TVectorD.h"
#include "TMath.h"
#include "TGraph.h"
#include "SNNakazatoIntegFcn.hh"

SNNakazatoIntegFcn::SNNakazatoIntegFcn(){
    for(int it=0; it<3; it++){
        grLuminosity[it] = NULL;
        grAlpha[it]      = NULL;
        grAverageE[it]   = NULL;
        timeMin[it] = -999.;
        timeMax[it] = -999.;
    }
}
SNNakazatoIntegFcn::SNNakazatoIntegFcn(int imode){
    for(int it=0; it<3; it++){
        grLuminosity[it] = NULL;
        grAlpha[it]      = NULL;
        grAverageE[it]   = NULL;
        timeMin[it] = -999.;
        timeMax[it] = -999.;
    }
    readFluxGraph(imode);
}
SNNakazatoIntegFcn::~SNNakazatoIntegFcn(){
    for(Int_t it=0; it<3; it++){
        if(grLuminosity[it])delete grLuminosity[it];
        if(grAverageE[it])delete grAverageE[it];
        if(grAlpha[it])delete grAlpha[it];
    }
}

void SNNakazatoIntegFcn::readFluxGraph(int imode){

    std::cout << "SNNakazatoIntegFcn readFluxGraph" << std::endl;
    //read in data
    TString path="/junofs/users/miaoyu/supernova/models/Nakazato_2013/";
    std::ifstream fSNmod_antinue;
    std::ifstream fSNmod_nue;
    std::ifstream fSNmod_nux;
    TString filename = Form("%s/neutrino_signal_nubar_e", path.Data());
    fSNmod_antinue.open(filename);
    filename = Form("%s/neutrino_signal_nu_e", path.Data());
    fSNmod_nue.open(filename);
    filename = Form("%s/neutrino_signal_nu_x", path.Data());
    fSNmod_nux.open(filename);

    if((!fSNmod_antinue) || (!fSNmod_nue) || (!fSNmod_nux)){
        std::cout << "Error: " << imode << " files of Nakazato models can't be opened!" << std::endl;
        exit(-1);
    }

    double time, luminosity, aveE, aveE2, alpha;

    //nu_e
    std::vector<double> t_nue;
    std::vector<double> lumin_nue;//foe/s=10^51 erg/s
    std::vector<double> averagE_nue;//MeV; 1erg = 6.24151*10^5MeV
    std::vector<double> alpha_nue;
    //std::vector<double> averagE2_nue;
    while(fSNmod_nue >> time >> luminosity >> aveE >> alpha){
        t_nue.push_back(time);
        vecTime.push_back(time);//to be accessed by the user
        lumin_nue.push_back(luminosity/1e51);
        averagE_nue.push_back(aveE);
        //averagE2_nue.push_back(aveE2);
        alpha_nue.push_back(alpha);
    }
    int npnt_nue = t_nue.size();
    TVectorD vecTime_nue(npnt_nue);
    TVectorD vecLumin_nue(npnt_nue);
    TVectorD vecAlpha_nue(npnt_nue);
    TVectorD vecAverageE_nue(npnt_nue);
    for(int ip=0; ip<npnt_nue; ip++){
        vecTime_nue[ip] = t_nue[ip];
        vecLumin_nue[ip] = lumin_nue[ip];
        vecAlpha_nue[ip] = alpha_nue[ip]; //1./(averagE2_nue[ip]/TMath::Power(averagE_nue[ip],2)-1)-1;
        vecAverageE_nue[ip] = averagE_nue[ip];
    }

    //anti_nu_e
    std::vector<double> t_antinue;
    std::vector<double> lumin_antinue;
    std::vector<double> averagE_antinue;
    std::vector<double> averagE2_antinue;
    while(fSNmod_antinue >> time >> luminosity >> aveE >> aveE2){
        t_antinue.push_back(time);
        lumin_antinue.push_back(luminosity/1e51);
        averagE_antinue.push_back(aveE);
        averagE2_antinue.push_back(aveE2);
    }
    int npnt_antinue = t_antinue.size();
    TVectorD vecTime_antinue(npnt_antinue);
    TVectorD vecLumin_antinue(npnt_antinue);
    TVectorD vecAlpha_antinue(npnt_antinue);
    TVectorD vecAverageE_antinue(npnt_antinue);
    for(int ip=0; ip<npnt_antinue; ip++){
        vecTime_antinue[ip] = t_antinue[ip];
        vecLumin_antinue[ip] = lumin_antinue[ip];
        vecAlpha_antinue[ip] = 1./(averagE2_antinue[ip]/TMath::Power(averagE_antinue[ip],2)-1)-1;
        vecAverageE_antinue[ip] = averagE_antinue[ip];
    }



    //nu_x
    std::vector<double> t_nux;
    std::vector<double> lumin_nux;
    std::vector<double> averagE_nux;
    std::vector<double> averagE2_nux;
    while(fSNmod_nux >> time >> luminosity >> aveE >> aveE2){
        t_nux.push_back(time);
        lumin_nux.push_back(luminosity/1e51);
        averagE_nux.push_back(aveE);
        averagE2_nux.push_back(aveE2);
    }
    int npnt_nux = t_nux.size();
    TVectorD vecTime_nux(npnt_nux);
    TVectorD vecLumin_nux(npnt_nux);
    TVectorD vecAlpha_nux(npnt_nux);
    TVectorD vecAverageE_nux(npnt_nux);
    for(int ip=0; ip<npnt_nux; ip++){
        vecTime_nux[ip] = t_nux[ip];
        vecLumin_nux[ip] = lumin_nux[ip];
        vecAlpha_nux[ip] = 1./(averagE2_nux[ip]/TMath::Power(averagE_nux[ip],2)-1)-1;
        vecAverageE_nux[ip] = averagE_nux[ip];
    }

    //parameter vs time graphs
    grLuminosity[0] = new TGraph(vecTime_nue, vecLumin_nue);
    grLuminosity[1] = new TGraph(vecTime_antinue, vecLumin_antinue);
    grLuminosity[2] = new TGraph(vecTime_nux, vecLumin_nux);

    grAlpha[0] = new TGraph(vecTime_nue, vecAlpha_nue);
    grAlpha[1] = new TGraph(vecTime_antinue, vecAlpha_antinue);
    grAlpha[2] = new TGraph(vecTime_nux, vecAlpha_nux);

    grAverageE[0] = new TGraph(vecTime_nue, vecAverageE_nue);
    grAverageE[1] = new TGraph(vecTime_antinue, vecAverageE_antinue);
    grAverageE[2] = new TGraph(vecTime_nux, vecAverageE_nux);


    //time limits
    timeMin[0] = t_nue[0];     timeMax[0] = t_nue[npnt_nue-1];
    timeMin[1] = t_antinue[0]; timeMax[1] = t_antinue[npnt_antinue-1];
    timeMin[2] = t_nux[0];     timeMax[2] = t_nux[npnt_nux-1];
    std::cout << timeMin[0]  << ";  " << timeMax[0] << ";"
        << timeMin[1]  << ";  " << timeMax[1] << ";"
        << timeMin[2]  << ";  " << timeMax[2] << ";"
        <<std::endl;
}


double SNNakazatoIntegFcn::fcnfluxtime(double* x, double* par){
    double time = x[0];
    double Ev   = par[0];
    int type = int(par[1]);

    double luminosity, A, alpha;
    if(type < 2){
        luminosity = grLuminosity[type]->Eval(time);
        A = grAverageE[type]->Eval(time);
        alpha = grAlpha[type]->Eval(time);
    }
    if(type >= 2){
        luminosity = grLuminosity[2]->Eval(time);
        A = grAverageE[2]->Eval(time);
        alpha = grAlpha[2]->Eval(time);
    }
   
    double index = 6.24151e56;
    double flux = (luminosity/A)*(TMath::Power(Ev,alpha)/TMath::Gamma(1+alpha))*TMath::Power((alpha+1)/A,alpha+1)*TMath::Exp(-(alpha+1)*Ev/A);
    return index*flux;
}

double SNNakazatoIntegFcn::getNumT(double time, int type){
    double num=0;
    if(type < 2 ){
        if(time < timeMin[type])return 0;
        if(time > timeMax[type])return 0;
        num = grLuminosity[type]->Eval(time)/(grAverageE[type]->Eval(time));
    }
    if(type >= 2){
        if(time < timeMin[2])return 0;
        if(time > timeMax[2])return 0;
        num = grLuminosity[2]->Eval(time)/(grAverageE[2]->Eval(time));
    }
    double index = 6.24151e56;
    return index*num ;
}

double SNNakazatoIntegFcn::getLumT(double time, int type){
    double lumi=0;
    if(type < 2 ) {
        if(time < timeMin[type])return 0;
        if(time > timeMax[type])return 0;
        lumi = grLuminosity[type]->Eval(time);
    }
    if(type >= 2) {
        if(time < timeMin[2])return 0;
        if(time > timeMax[2])return 0;
        lumi = grLuminosity[2]->Eval(time);
    }
    return 1e51*lumi;
}

double SNNakazatoIntegFcn::getEventAtTime(double time, double E, int type){
    double luminosity, A, alpha;
    if(type < 2){
        if(time < timeMin[type]) return 0;
        if(time > timeMax[type]) return 0;
        luminosity = grLuminosity[type]->Eval(time);
        A = grAverageE[type]->Eval(time);
        alpha = grAlpha[type]->Eval(time);
        if (A == 0 ) return 0 ;
        // test line
        //std::cout << "getEventAtTime: time, luminosity, A, alpha: " 
        //         << time << ", " << luminosity << ", " << A << ", " << alpha << std::endl;
        //
    }
    if(type >= 2){
        if(time < timeMin[2])return 0;
        if(time > timeMax[2])return 0;
        luminosity = grLuminosity[2]->Eval(time);
        A = grAverageE[2]->Eval(time);
        alpha = grAlpha[2]->Eval(time);
        if (A == 0 ) return 0 ;
    }
   
    double index = 6.24151e56;
    double flux = (luminosity/A)*(TMath::Power(E,alpha)/TMath::Gamma(1+alpha))*TMath::Power((alpha+1)/A,alpha+1)*TMath::Exp(-(alpha+1)*E/A);
    return index*flux;
}

double SNNakazatoIntegFcn::getAverageET(double time, int type){
    double aveE = 0;
    //std::cout << "time " << time <<"; type" << type << std::endl;
    if(type < 2) {
        if(timeMin[type]<time && timeMax[type]>time)aveE = grAverageE[type]->Eval(time);
        else aveE = 0;
    }
    else {
        if(timeMin[2] < time && timeMax[2]>time)aveE = grAverageE[2]->Eval(time);
        else aveE = 0;
    }
    
    return aveE;
}

