#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "TVectorD.h"
#include "TMath.h"
#include "TGraph.h"

#include "SNBurrowsIntegFcn.hh"

SNBurrowsIntegFcn::SNBurrowsIntegFcn(){
    for(int it=0;it<3;it++){
        timeMin[it] = -999.;
        timeMax[it] = -999.;
        grLuminosity[it] = NULL;
        grAverageE[it] = NULL;
        grAlpha[it]    = NULL;
    }
    for (int ii=0; ii<80; ii++) {
        for (int jj=0; jj<100; jj++) {
            etspec_nua[ii][jj] = 0.;
            etspec_nue[ii][jj] = 0.;
            etspec_nux[ii][jj] = 0.;
        }
    }
}

SNBurrowsIntegFcn::SNBurrowsIntegFcn(int imode){
    for(int it=0;it<3;it++){
        timeMin[it] = -999.;
        timeMax[it] = -999.;
        grLuminosity[it] = NULL;
        grAverageE[it] = NULL;
        grAlpha[it]    = NULL;
    }
    for (int ii=0; ii<80; ii++) {
        for (int jj=0; jj<100; jj++) {
            etspec_nua[ii][jj] = 0.;
            etspec_nue[ii][jj] = 0.;
            etspec_nux[ii][jj] = 0.;
        }
    }
    //readFluxGraph(imode);
    readNumFlux(imode);
}

SNBurrowsIntegFcn::~SNBurrowsIntegFcn(){

}

void SNBurrowsIntegFcn::readNumFlux(int imass) {

    imass = imass - 6500;
    std::ifstream fSNmod_antinue;
    std::ifstream fSNmod_nue;
    std::ifstream fSNmod_nux;
    TString filename;
    filename = Form("/junofs/users/miaoyu/supernova/models/Fornax_2021/nuebar_%dM.dat", imass);
    std::cout << filename << std::endl;
    fSNmod_antinue.open(filename);

    filename = Form("/junofs/users/miaoyu/supernova/models/Fornax_2021/nux_%dM.dat", imass);
    std::cout << filename << std::endl;
    fSNmod_nux.open(filename);
    
    filename = Form("/junofs/users/miaoyu/supernova/models/Fornax_2021/nue_%dM.dat", imass);
    std::cout << filename << std::endl;
    fSNmod_nue.open(filename);

    if( (!fSNmod_nue)){
        std::cout << "Error: " << imass << " mass nue files of Burrows models can't be opened!" << std::endl;
        exit(-1);
    }

    if((!fSNmod_antinue) ){
        std::cout << "Error: " << imass << " mass antinue files of Burrows models can't be opened!" << std::endl;
        exit(-1);
    }
    if((!fSNmod_nux)){
        std::cout << "Error: " << imass << " mass nux files of Burrows models can't be opened!" << std::endl;
        exit(-1);
    }



    double time;
    std::string line;
    int iline = 0; 
    while(getline(fSNmod_antinue, line)) {
        std::istringstream ss(line);
        double s; 
        int ipar = 0;
        while(ss >> s) {
            if(ipar > 0) {
                etspec_nua[iline][ipar-1] = s / 1e50;
            }
            ipar++;
        }
        iline++;
    }
    iline = 0;
    while(getline(fSNmod_nue, line)) {
        std::istringstream ss(line);
        double s; 
        int ipar = 0;
        while(ss >> s) {
            if(ipar > 0) {
                etspec_nue[iline][ipar-1] = s / 1e50;
            }
            ipar++;
        }
        iline++;
    }
    iline = 0;
    while(getline(fSNmod_nux, line)) {
        std::istringstream ss(line);
        double s; 
        int ipar = 0;
        while(ss >> s) {
            if (ipar == 0) {
                vecTime.push_back(s);
            }
            if(ipar > 0) {
                etspec_nux[iline][ipar-1] = s / 1e50;
            }
            ipar++;
        }
        iline++;
    }
    std::cout << ">>>>> readNumFlux successfully!!!" << std::endl;

    timeMin[0] = vecTime[0]; timeMax[0] = vecTime[79];
    timeMin[1] = vecTime[0]; timeMax[1] = vecTime[79];
    timeMin[2] = vecTime[0]; timeMax[2] = vecTime[79];
    std::cout << timeMin[0]  << ";  " << timeMax[0] << ";"
        << timeMin[1]  << ";  " << timeMax[1] << ";"
        << timeMin[2]  << ";  " << timeMax[2] << ";"
        <<std::endl;
}

void SNBurrowsIntegFcn::readFluxGraph(int imode){
    //TString path = "/junofs/users/miaoyu/supernova/production/nuFlux/Fornax_2021/";
    TString path="/junofs/users/miaoyu/supernova/wenlj/YB_fromSZ/";
    std::ifstream fSNmod_antinue;
    std::ifstream fSNmod_nue;
    std::ifstream fSNmod_nux;
    fSNmod_antinue.open("/junofs/users/miaoyu/supernova/wenlj/YB_fromSZ/na");
    fSNmod_nue.open("/junofs/users/miaoyu/supernova/wenlj/YB_fromSZ/ne");
    fSNmod_nux.open("/junofs/users/miaoyu/supernova/wenlj/YB_fromSZ/nx");

    if((!fSNmod_antinue) || (!fSNmod_nue) || (!fSNmod_nux)){
        std::cout << "Error: " << imode << " files of Burrows models can't be opened!" << std::endl;
        exit(-1);
    }

    double time, luminosity, aveE, aveE2;
    //nu_e
    std::vector<double> t_nue;
    std::vector<double> lumin_nue;//foe/s=10^51 erg/s
    std::vector<double> averagE_nue;//MeV; 1erg = 6.24151*10^5MeV
    std::vector<double> averagE2_nue;
    while(fSNmod_nue >> time >> luminosity >> aveE >> aveE2){
        t_nue.push_back(time);
        vecTime.push_back(time);//to be accessed by the user
        lumin_nue.push_back(luminosity);
        averagE_nue.push_back(aveE);
        averagE2_nue.push_back(aveE2);
    }
    int npnt_nue = t_nue.size();
    TVectorD vecTime_nue(npnt_nue);
    TVectorD vecLumin_nue(npnt_nue);
    TVectorD vecAlpha_nue(npnt_nue);
    TVectorD vecAverageE_nue(npnt_nue);
    for(int ip=0; ip<npnt_nue; ip++){
        vecTime_nue[ip] = t_nue[ip];
        vecLumin_nue[ip] = lumin_nue[ip];
        vecAlpha_nue[ip] = 1./(averagE2_nue[ip]/TMath::Power(averagE_nue[ip],2)-1)-1;
        vecAverageE_nue[ip] = averagE_nue[ip];
    }


    //anti_nu_e
    std::vector<double> t_antinue;
    std::vector<double> lumin_antinue;
    std::vector<double> averagE_antinue;
    std::vector<double> averagE2_antinue;
    while(fSNmod_antinue >> time >> luminosity >> aveE >> aveE2){
        t_antinue.push_back(time);
        lumin_antinue.push_back(luminosity);
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
        lumin_nux.push_back(luminosity);
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

double SNBurrowsIntegFcn::fcnfluxtime(double* x, double* par) {
    double time = x[0];
    double Ev   = par[0];
    int type    = int(par[1]);

    double luminosity, A, alpha;
    if (type<2) {
        luminosity = grLuminosity[type]->Eval(time);
        A = grAverageE[type]->Eval(time);
        alpha = grAlpha[type]->Eval(time);
    }
    if (type >=2) {
        luminosity = grLuminosity[2]->Eval(time);
        A = grAverageE[2]->Eval(time);
        alpha = grAlpha[2]->Eval(time);
    }

    double index = 6.24151e56;
    double flux = (luminosity/A)*(TMath::Power(Ev,alpha)/TMath::Gamma(1+alpha))*TMath::Power((alpha+1)/A,alpha+1)*TMath::Exp(-(alpha+1)*Ev/A);
    return index*flux;
}

double SNBurrowsIntegFcn::getNumT(double time, int type) {
    if (type > 2) type = 2;
    if(time < timeMin[type])return 0;
    if(time > timeMax[type])return 0;
    int binlow = int(time/0.001);
    int binhig = binlow + 1;
    if (binlow > 80 or binhig > 80) {
        binlow = 79;
        binhig = 80;
        std::cout << "Error !!! ROI has exceeded the Burrows PDF range!!!!" << std::endl;
    }
    double index = 624151.;
    double lowVal = 0;
    double higVal = 0; 
    for( int ibin=0; ibin<100; ibin++) {
        if (type == 0)
            lowVal += etspec_nue[binlow][ibin] ;
        if (type == 1)
            lowVal += etspec_nua[binlow][ibin] ;
        if (type >= 2)
            lowVal += etspec_nux[binlow][ibin] ;
    }
    for( int ibin=0; ibin<100; ibin++) {
        if (type == 0)
            higVal += etspec_nue[binhig][ibin] ;
        if (type == 1)
            higVal += etspec_nua[binhig][ibin] ;
        if (type >= 2)
            higVal += etspec_nux[binhig][ibin] ;
    }

    double r = (time - binlow * 0.001) / 0.001;
    double num = r * higVal + (1 - r) * lowVal;

    return num;
}

double SNBurrowsIntegFcn::getLumT(double time, int type)
{
    double lumi=0;
    if (type > 2) type = 2;
    if(time < timeMin[type])return 0;
    if(time > timeMax[type])return 0;
    int binlow = int(time/0.001);
    int binhig = binlow + 1;
    if (binlow > 79 or binhig > 79) {
        binlow = 78;
        binhig = 79;
        std::cout << "Error !!! ROI has exceeded the Burrows PDF range!!!!" << std::endl;
    }
    double index = 624151;
    double lowVal = 0;
    double higVal = 0; 
    for( int ibin=0; ibin<100; ibin++) {
        if (type == 0)
            lowVal += etspec_nue[binlow][ibin] * (ibin * 1) / index;
        if (type == 1)
            lowVal += etspec_nua[binlow][ibin] * (ibin * 1) / index;
        if (type >= 2)
            lowVal += etspec_nux[binlow][ibin] * (ibin * 1) / index ;
    }
    for( int ibin=0; ibin<100; ibin++) {
        if (type == 0)
            higVal += etspec_nue[binhig][ibin] * (ibin * 1) / index;
        if (type == 1)
            higVal += etspec_nua[binhig][ibin] * (ibin * 1) / index;
        if (type >= 2)
            higVal += etspec_nux[binhig][ibin] *(ibin * 1) / index ;
    }

    double r = (time - binlow * 0.001) / 0.001;
    lumi = r * higVal + (1 - r) * lowVal;
    
    return lumi;
}


double SNBurrowsIntegFcn::getEventAtTime(double time, double E, int type){
    if (type > 2) type = 2;
    if(time < timeMin[type]) return 0;
    if(time > timeMax[type]) return 0;
    
    int Ebinlow = int(E);
    int Ebinhig = Ebinlow + 1;
    if (Ebinlow >= 100 or Ebinhig >= 100) {
        std::cout << "Error !!! ROI has exceeded the Burrows PDF energy range !!!!!" << std::endl;
        Ebinlow = 98;
        Ebinhig = 99;
    }
    int binlow = int(time/0.001);
    int binhig = binlow + 1;
    if (binlow > 79 or binhig > 79) {
        binlow = 78;
        binhig = 79;
        std::cout << "Error !!! ROI has exceeded the Burrows PDF time range!!!!" << std::endl;
    }

    double E2 = Ebinhig * 1;   // MeV
    double E1 = Ebinlow * 1;   // MeV
    double t1 = binlow * 0.001;   // s
    double t2 = binhig * 0.001;   // s
    double w11 = (t2 - time) * (E2 - E) / (t2 - t1) / (E2 - E1);
    double w12 = (t2 - time) * (E - E1) / (t2 - t1) / (E2 - E1) ;
    double w21 = (time - t1) * (E2 - E) / (t2 - t1) / (E2 - E1);
    double w22 = (time - t1) * (E - E1) / (t2 - t1) / (E2 - E1);

    double val = 0;
    if (type == 0)  {
        double val = w11 * etspec_nue[binlow][Ebinlow] + w12 * etspec_nue[binlow][Ebinhig] + w21 * etspec_nue[binhig][Ebinlow] + w22 * etspec_nue[binhig][Ebinhig];
        return val;
    }
    else if (type == 1) {
        double val = w11 * etspec_nua[binlow][Ebinlow] + w12 * etspec_nua[binlow][Ebinhig] + w21 * etspec_nua[binhig][Ebinlow] + w22 * etspec_nua[binhig][Ebinhig];
        return val;
    }
    else if (type == 2) {
        double val = (w11 * etspec_nux[binlow][Ebinlow] + w12 * etspec_nux[binlow][Ebinhig] + w21 * etspec_nux[binhig][Ebinlow] + w22 * etspec_nux[binhig][Ebinhig]) ;
        return val;
    }
}

double SNBurrowsIntegFcn::getAverageET(double time, int type){
    double aveE = 0;
    if (type > 2) type = 2;
    if(time < timeMin[type]) return 0;
    if(time > timeMax[type]) return 0;

    int binlow = int(time/0.001);
    int binhig = binlow + 1;
    if (binlow > 79 or binhig > 79) {
        binlow = 78;
        binhig = 79;
        std::cout << "Error !!! ROI has exceeded the Burrows PDF time range!!!!" << std::endl;
    }

    double lowVal = 0; double lowN = 0;
    double higVal = 0; double higN = 0;

    for(int i=0; i<100; i++) {
        double E = i;
        if (type == 0) {
            lowVal += etspec_nue[binlow][i] * E;
            lowN   += etspec_nue[binlow][i];
        }
        if (type == 1) {
            lowVal += etspec_nua[binlow][i] * E;
            lowN   += etspec_nua[binlow][i];
        }
        if (type == 2) {
            lowVal += etspec_nux[binlow][i] * E;
            lowN   += etspec_nux[binlow][i];
        }
    }

    for(int i=0; i<100; i++) {
        double E = i;
        if (type == 0) {
            higVal += etspec_nue[binhig][i] * E;
            higN   += etspec_nue[binhig][i];
        }
        if (type == 1) {
            higVal += etspec_nua[binhig][i] * E;
            higN   += etspec_nua[binhig][i];
        }
        if (type == 2) {
            higVal += etspec_nux[binhig][i] * E;
            higN   += etspec_nux[binhig][i];
        }
    }

    if(0) {
        std::cout << time << " " 
                  << lowN << " " << lowVal  << " "
                  << higN << " " << higVal  << " "
                  << std::endl;
    }

    lowVal /= lowN;
    higVal /= higN;
    
    aveE = (time - binlow*0.001)/0.001 * higVal + (binhig*0.001 - time)/0.001 * lowVal;
    return aveE;
}

