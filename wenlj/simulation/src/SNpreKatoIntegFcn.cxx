#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <array>

#include "TGraph.h"
#include "TVectorD.h"
#include "TSpline.h"
#include "TMath.h"
#include "SNpreKatoIntegFcn.hh"

SNpreKatoIntegFcn::SNpreKatoIntegFcn(){
    initialize();
}

SNpreKatoIntegFcn::SNpreKatoIntegFcn(int imode){
    initialize();
    readFluxGraph(imode);
}

void SNpreKatoIntegFcn::initialize(){
    for(int it=0; it<4; it++){
        timeMin[it] = -999;
        timeMax[it] = -999;
    }
    for(int it=0; it<4; it++){
        for(int ip=0; ip<501; ip++){
            grBinNum_T[it][ip]     = NULL;
        }
        grNum_T[it] = NULL ;
        grLum_T[it] = NULL;
        grNum_E[it] = NULL;
    }
}

SNpreKatoIntegFcn::~SNpreKatoIntegFcn(){
}


void SNpreKatoIntegFcn::readFluxGraph(int imode){
    TString filePath;
    if(imode == 20150) filePath = "/junofs/users/miaoyu/supernova/wenlj/simulation/data/preSN/kato/15";
    if(imode == 20120) filePath = "/junofs/users/miaoyu/supernova/wenlj/simulation/data/preSN/kato/12";

    const int npnt_E = 501;//energy bin
    //-------reading nue file--------//
    readFile(Form("%s/total_nue",filePath.Data()),0);
    readFile(Form("%s/total_nueb",filePath.Data()),1);
    readFile(Form("%s/total_nux",filePath.Data()),2);
    readFile(Form("%s/total_nuxb",filePath.Data()),3);

    double** binNum = new double*[4];
    for(int itype=0; itype<4; itype++){
        binNum[itype] = new double[npnt_E];
        for(int ibin=0; ibin<npnt_E; ibin++){
            //binNum[itype][ibin] = integIbinTgraph(grBinNum_T[itype][ibin], timeMin[itype], timeMax[itype]);
            double factor = 3600*24;
            binNum[itype][ibin] = integIbinTgraph(grBinNum_T[itype][ibin], -1*factor,-1e-6*factor);
        }
        grNum_E[itype] = new TGraph(npnt_E, binE_E, binNum[itype]);
    }

    
}

void SNpreKatoIntegFcn::readFile(TString path, int itype){


    std::ifstream fSNmod;
    //-----light curve-----//
    fSNmod.open(Form("%s/lightcurve_all.dat",path.Data()));
    if(!fSNmod){
        std::cout << "Error to open the model file!!!" << std::endl;
        exit(-1);
    } 
    
    vecTime.clear();
    std::vector<int> vecStep;
    std::vector<double> vecNumLumTot;
    std::vector<double> vecELumTot; 
    double factorTime = -1;
    double itime;
    double inumLum;
    double ieLum;
    int istep;
    while(fSNmod>>itime){
        fSNmod >> inumLum >> ieLum >> istep;
        vecTime.push_back(itime*factorTime);
        vecStep.push_back(istep);
        vecNumLumTot.push_back(inumLum);
        vecELumTot.push_back(ieLum);
    }
    fSNmod.close();
    fSNmod.clear();

    long int nsteps = vecStep.size();
    timeMin[itype] = vecTime[0];
    timeMax[itype] = vecTime[nsteps-1];
    grNum_T[itype] = new TGraph(nsteps);
    grLum_T[itype] = new TGraph(nsteps);
    for(int ip=0; ip<nsteps; ip++){
        grNum_T[itype]->SetPoint(ip, vecTime[ip], vecNumLumTot[ip]);
        grLum_T[itype]->SetPoint(ip, vecTime[ip], vecELumTot[ip]);
    }
    //-------energy spectrum for each time------//
    const int npnt_E = 501;//energy bin
    for(int ibin=0; ibin<npnt_E; ibin++){
        grBinNum_T[itype][ibin] = new TGraph(nsteps);
    }
    for(int ip=0; ip<nsteps; ip++){
        fSNmod.open(Form("%s/spe_all0%d.dat", path.Data(),vecStep[ip]));
        double energy, flux;
        int icout = 0;
        while(fSNmod>>energy>>flux){
            binE_E[icout] = energy;
            grBinNum_T[itype][icout]->SetPoint(ip,vecTime[ip],flux);
            icout++;
        }
        fSNmod.close();
        fSNmod.clear();
    }

}



double SNpreKatoIntegFcn::integIbinTgraph(TGraph* igraph, double first, double last){
    int nstep = 1000;
    double step = (last-first)/nstep;
    double integ = 0;
    for(int ist=0; ist<nstep; ist++){
        double time = first+step*ist;
        integ += step*igraph->Eval(time);
    }

    return integ;
}

//FCN for SN src 
double SNpreKatoIntegFcn::fluxTimeE(double* x, double* par){
    double Ev = x[0];
    double tfirst = par[0];
    double tlast  = par[1];
    int itype = (int)(par[2]);

    if(itype==2 || itype==4)itype=2;
    if(itype==3 || itype==5)itype=3;
    double flux = 0;
    int npnt_E = 501;
    double* NumE = new double[npnt_E];
    for(int ib=0; ib<npnt_E; ib++){
        NumE[ib] = integIbinTgraph(grBinNum_T[itype][ib],tfirst,tlast);
        //std::cout <<"ib" << ib <<"; NumE " << NumE[ib] << std::endl;
    }
    TGraph gr(npnt_E,binE_E, NumE);
    flux = gr.Eval(Ev);
    if(flux<0)flux=0;
    delete[] NumE;
    return flux;
}


double SNpreKatoIntegFcn::getFluxE(double E, int type){
    if(type==2 || type==4)type=2;
    if(type==3 || type==5)type=3;
    double flux = grNum_E[type]->Eval(E); 
    if(flux<0)flux=0;
    return flux;
}

double SNpreKatoIntegFcn::getNumT(double T, int type){
    if(type==2 || type==4)type=2;
    if(type==3 || type==5)type=3;
    if(T < timeMin[type] || T>timeMax[type]){
        //std::cout << "Time range Error.." << std::endl;
        return 0;
    }
    //if(T>timeMax[type])return 0;
    double flux = grNum_T[type]->Eval(T); 
    if(flux<0)
    {
        flux=0;
    }
    return flux;
}
double SNpreKatoIntegFcn::getLumT(double T, int type){
    if(type==2 || type==4)type=2;
    if(type==3 || type==5)type=3;
    if(T < timeMin[type] || T>timeMax[type])return 0;
    //if(T>timeMax[type]){
    //    std::cout << "T " << T << "; maxT " << timeMax[type] << std::endl;
    //    return 0;
    //}
    double lumi = grLum_T[type]->Eval(T);
    if(lumi<0)lumi=0;
    double factor = 1.60218e-6;//MeV to erg
    return factor*lumi;
}

double SNpreKatoIntegFcn::getEventAtTime(double time, double E, int type){
    if(type==2 || type==4)type=2;
    if(type==3 || type==5)type=3;
    double flux = 0; 
    int npnt_E = 501;
    double* NumE = new double[npnt_E];
    if(time < timeMin[type] || time>timeMax[type])return 0;
    //if(time>timeMax[type])return 0;
    for(int ib=0; ib<npnt_E; ib++){
        NumE[ib] = grBinNum_T[type][ib]->Eval(time);
    }
    TGraph* gr = new TGraph(npnt_E, binE_E, NumE);
    flux = gr->Eval(E);
    if(flux<0)flux=0;
    delete[] NumE;
    delete gr;
    return flux;
}

double SNpreKatoIntegFcn::getAverageET(double time, int type){
    if(type==2 || type==4)type=2;
    if(type==3 || type==5)type=3;
    double aveE = 0;
    if(time < timeMin[type] || time > timeMax[type])return 0;
    //if(time > timeMax[type])return 0;
    aveE = (grLum_T[type]->Eval(time))/(grNum_T[type]->Eval(time));
    return aveE;
}
