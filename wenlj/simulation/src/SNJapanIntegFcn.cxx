#include <iostream>
#include <fstream>
#include <vector>
#include <array>

#include "TGraph.h"
#include "TVectorD.h"
#include "TSpline.h"
#include "SNJapanIntegFcn.hh"

SNJapanIntegFcn::SNJapanIntegFcn(){
    nbin_E = 0;
    binE_E = NULL;
    for(int it=0; it<3; it++){
        timeMin[it] = -999;
        timeMax[it] = -999;
    }
    grNum_nux_T     = NULL;
    grNum_antinue_T = NULL;
    grNum_nux_T     = NULL;
    grLum_nux_T     = NULL;
    grLum_antinue_T = NULL;
    grLum_nux_T     = NULL;
    grNum_nux_E     = NULL;
    grNum_antinue_E = NULL;
    grNum_nux_E     = NULL;

}

SNJapanIntegFcn::SNJapanIntegFcn(int imode){
    nbin_E = 0;
    binE_E = NULL;
    for(int it=0; it<3; it++){
        timeMin[it] = -999;
        timeMax[it] = -999;
    }
    readFluxGraph(imode);
}

SNJapanIntegFcn::~SNJapanIntegFcn(){
}


void SNJapanIntegFcn::readFluxGraph(int imode){
    TString path;
    if(imode<6000) {
        path = Form("/junofs/users/miaoyu/supernova/wenlj/simulation/data/Japan/intpdata/intp%d.data",imode);
        nbin_E = 20;
    }
    else if(imode== 9001){
        path = "/junofs/users/miaoyu/supernova/wenlj/simulation/data/Japan/cooling/spectob147S.data";
        nbin_E = 25;
    }
    std::ifstream fSNmod;
    fSNmod.open(path);
    if(!fSNmod){
        std::cout << "Error: /junofs/users/miaoyu/supernova/wenlj/simulation/data/Japan/..." << imode <<".data can't be opened!" << std::endl;
        exit(-1);
    }

    std::vector<double> vecTime;
    std::vector<std::vector<double>> vecNum_nue;
    std::vector<std::vector<double>> vecNum_antinue;
    std::vector<std::vector<double>> vecNum_nux;

    std::vector<std::vector<double>> vecLum_nue;
    std::vector<std::vector<double>> vecLum_antinue;
    std::vector<std::vector<double>> vecLum_nux;
    double* binE = new double[nbin_E+1];
    double time;

    //reading from file
    while(fSNmod >> time){
        vecTime.push_back(time);
        std::vector<double> arrNum_nue;
        std::vector<double> arrNum_antinue;
        std::vector<double> arrNum_nux;

        std::vector<double> arrLum_nue;
        std::vector<double> arrLum_antinue;
        std::vector<double> arrLum_nux;
       

        for(int ib=0; ib<nbin_E; ib++){
            double iLum[3] = {0};
            double iNum[3] = {0};
            fSNmod >> binE[ib] >> binE[ib+1];
            fSNmod >> iNum[0] >> iNum[1] >> iNum[2];
            fSNmod >> iLum[0] >> iLum[1] >> iLum[2];

            arrNum_nue.push_back(iNum[0]);
            arrNum_antinue.push_back(iNum[1]);
            arrNum_nux.push_back(iNum[2]);
            
            arrLum_nue.push_back(iLum[0]);
            arrLum_antinue.push_back(iLum[1]);
            arrLum_nux.push_back(iLum[2]);
        }

        vecNum_nue.push_back(arrNum_nue);
        vecNum_antinue.push_back(arrNum_antinue);
        vecNum_nux.push_back(arrNum_nux);

        vecLum_nue.push_back(arrLum_nue);
        vecLum_antinue.push_back(arrLum_antinue);
        vecLum_nux.push_back(arrLum_nux);
    }
    //flux distribution
    //number and luminosity time distribution of each energy bin 
    int ntbin = vecTime.size();
    for(int ity=0; ity<3; ity++){
        timeMin[ity] = vecTime[0];
        timeMax[ity] = vecTime[ntbin-1];
    }
    double* binTime = new double[ntbin];
    double** ibinNum_nue_T     = new double*[nbin_E];
    double** ibinNum_antinue_T = new double*[nbin_E];
    double** ibinNum_nux_T     = new double*[nbin_E];
    
    for(int ib=0;ib<nbin_E;ib++){
        ibinNum_nue_T[ib]     = new double[ntbin];
        ibinNum_antinue_T[ib] = new double[ntbin];
        ibinNum_nux_T[ib]     = new double[ntbin];
    }
    
    double* timeNum_nue     = new double[ntbin];
    double* timeNum_antinue = new double[ntbin];
    double* timeNum_nux     = new double[ntbin];
   
    double* timeLum_nue     = new double[ntbin];
    double* timeLum_antinue = new double[ntbin];
    double* timeLum_nux     = new double[ntbin];
    for(int it=0; it<ntbin; it++){
        binTime[it] = vecTime[it];
        timeNum_nue[it]     = 0;
        timeNum_antinue[it] = 0;
        timeNum_nux[it]     = 0;
        timeLum_nue[it]     = 0;
        timeLum_antinue[it] = 0;
        timeLum_nux[it]     = 0;
        for(int ib=0; ib<nbin_E; ib++){
            ibinNum_nue_T[ib][it]     = vecNum_nue[it].at(ib);
            ibinNum_antinue_T[ib][it] = vecNum_antinue[it].at(ib);
            ibinNum_nux_T[ib][it]     = vecNum_nux[it].at(ib);

            timeNum_nue[it]     += (binE[ib+1]-binE[ib])*vecNum_nue[it].at(ib);
            timeNum_antinue[it] += (binE[ib+1]-binE[ib])*vecNum_antinue[it].at(ib);
            timeNum_nux[it]     += (binE[ib+1]-binE[ib])*vecNum_nux[it].at(ib);
            
            timeLum_nue[it]     += (binE[ib+1]-binE[ib])*vecLum_nue[it].at(ib);
            timeLum_antinue[it] += (binE[ib+1]-binE[ib])*vecLum_antinue[it].at(ib);
            timeLum_nux[it]     += (binE[ib+1]-binE[ib])*vecLum_nux[it].at(ib);
        }
    }
    
    for(int ib=0; ib<nbin_E; ib++){
        TGraph* igrNue  = new TGraph(ntbin,binTime, ibinNum_nue_T[ib]);
        TGraph* igrNueb = new TGraph(ntbin,binTime, ibinNum_antinue_T[ib]);
        TGraph* igrNux  = new TGraph(ntbin,binTime, ibinNum_nux_T[ib]);
        grBinNum_nue_T.push_back(igrNue);
        grBinNum_antinue_T.push_back(igrNueb);
        grBinNum_nux_T.push_back(igrNux);
    }

    //time distribution for all bins
    //number
    grNum_nue_T     = new TGraph(ntbin, binTime, timeNum_nue);
    grNum_antinue_T = new TGraph(ntbin, binTime, timeNum_antinue);
    grNum_nux_T     = new TGraph(ntbin, binTime, timeNum_nux);

    grLum_nue_T     = new TGraph(ntbin, binTime, timeLum_nue);
    grLum_antinue_T = new TGraph(ntbin, binTime, timeLum_antinue);
    grLum_nux_T     = new TGraph(ntbin, binTime, timeLum_nux);
    //number distribution by integrating all the time
    double* binNum_nue     = new double[nbin_E];
    double* binNum_antinue = new double[nbin_E];
    double* binNum_nux     = new double[nbin_E];
    binE_E = new double[nbin_E];
    for(int ib=0; ib<nbin_E; ib++){
        binE_E[ib] = binE[ib+1];
        binNum_nue[ib]     = integIbinTgraph(grBinNum_nue_T[ib], timeMin[0],timeMax[0]);
        binNum_antinue[ib] = integIbinTgraph(grBinNum_antinue_T[ib], timeMin[1],timeMax[1]);
        binNum_nux[ib]     = integIbinTgraph(grBinNum_nux_T[ib], timeMin[2],timeMax[2]);
    }

    grNum_nue_E     = new TGraph(nbin_E, binE_E, binNum_nue);
    grNum_antinue_E = new TGraph(nbin_E, binE_E, binNum_antinue);
    grNum_nux_E     = new TGraph(nbin_E, binE_E, binNum_nux);

}


double SNJapanIntegFcn::integIbinTgraph(TGraph* igraph, double first, double last){
    int nstep = 1000;
    double step = (last-first)/nstep;
    double integ = 0;
    for(int ist=0; ist<nstep; ist++){
        double time = first+step*ist;
        integ += step*igraph->Eval(time,0,"S");
    }

    return integ;
}

//FCN for SN src 
double SNJapanIntegFcn::fluxTimeE(double* x, double* par){
    double Ev = x[0];
    double tfirst = par[0];
    double tlast  = par[1];
    int itype = (int)(par[2]);

    double flux = 0;
    double* NumE = new double[nbin_E];

    if( itype==0 ){
        for(int ib=0; ib<nbin_E; ib++){
            NumE[ib] = integIbinTgraph(grBinNum_nue_T[ib],tfirst,tlast);
        }
        TGraph gr(nbin_E,binE_E, NumE);
        flux = gr.Eval(Ev,0,"S");
    }
    if(itype == 1){
        for(int ib=0; ib<nbin_E; ib++){
            NumE[ib] = integIbinTgraph(grBinNum_antinue_T[ib],tfirst,tlast);
        }
        TGraph gr(nbin_E,binE_E, NumE);
        flux = gr.Eval(Ev,0,"S");
    }
    if(itype > 1){
        for(int ib=0; ib<nbin_E; ib++){
            NumE[ib] = integIbinTgraph(grBinNum_nux_T[ib],tfirst,tlast);
        }
        TGraph gr(nbin_E,binE_E, NumE);
        flux = gr.Eval(Ev,0,"S");
    }
    if(flux<0)flux=0;
    delete[] NumE;
    return flux;
}


double SNJapanIntegFcn::getFluxE(double E, int type){
    double flux = 0;
    if(type == 0) flux = grNum_nue_E->Eval(E,0,"S"); 
    if(type == 1) flux = grNum_antinue_E->Eval(E,0,"S"); 
    if(type > 1)  flux = grNum_nux_E->Eval(E,0,"S");
    if(flux<0)flux=0;
    return flux;
}

double SNJapanIntegFcn::getNumT(double T, int type){
    double flux = 0;
    if(type == 0) {
        if(T < timeMin[0] || T>timeMax[0])return 0;
        flux = grNum_nue_T->Eval(T); 
    }
    if(type == 1){
        if(T < timeMin[1] || T>timeMax[1])return 0;
        flux = grNum_antinue_T->Eval(T); 
    }
    if(type > 1){
        if(T < timeMin[2] || T>timeMax[2])return 0;
        flux = grNum_nux_T->Eval(T);
    }
    if(flux<0)flux=0;
    return flux;
}
double SNJapanIntegFcn::getLumT(double T, int type){
    double lumi = 0;
    if(type == 0){
        if(T < timeMin[0] || T>timeMax[0])return 0;
        lumi = grLum_nue_T->Eval(T); 
    }
    if(type == 1){
        if(T < timeMin[1] || T>timeMax[1])return 0;
        lumi = grLum_antinue_T->Eval(T); 
    }
    if(type > 1){
        if(T < timeMin[2] || T>timeMax[2])return 0;
        lumi = grLum_nux_T->Eval(T);
    }
    if(lumi<0)lumi=0;
    return lumi;
}

double SNJapanIntegFcn::getEventAtTime(double time, double E, int type){
    double flux = 0; 
    double* NumE = new double[nbin_E];
    if(type == 0){
        if(time < timeMin[0] || time>timeMax[0]){
            delete[] NumE;
            return 0;
        }
        for(int ib=0; ib<nbin_E; ib++){
            NumE[ib] = grBinNum_nue_T[ib]->Eval(time,0,"S");
        }
    }
    if(type == 1){
        if(time < timeMin[1] || time>timeMax[1]){
            delete[] NumE;
            return 0;
        }
        for(int ib=0; ib<nbin_E; ib++){
            NumE[ib] = grBinNum_antinue_T[ib]->Eval(time, 0, "S");
        }
    }
    if(type > 1){
        if(time < timeMin[2] || time>timeMax[2]){
            delete[] NumE;
            return 0;
        }
        for(int ib = 0; ib<nbin_E; ib++){
            NumE[ib] = grBinNum_nux_T[ib]->Eval(time, 0, "S");
        }
    }
    TGraph gr(nbin_E, binE_E, NumE);
    flux = gr.Eval(E, 0, "S");
    if(flux<0)flux=0;
    delete[] NumE;
    return flux;
}

double SNJapanIntegFcn::getAverageET(double time, double type){
    double aveE = 0;
    if(type == 0) {
        if(time < timeMin[0] || time > timeMax[0])return 0;
        aveE = grLum_nue_T->Eval(time)/(grNum_nue_T->Eval(time));
    }
    if(type == 1){
        if(time < timeMin[1] || time>timeMax[1])return 0;
        aveE = grLum_antinue_T->Eval(time)/(grNum_antinue_T->Eval(time));
    }
    if(type > 1){
        if(time < timeMin[2] || time>timeMax[2])return 0;
        aveE = grLum_nux_T->Eval(time)/(grNum_nux_T->Eval(time));
    }
    if(aveE<0)aveE=0;
    double index = 6.24151e5;
    return index*aveE;
}
