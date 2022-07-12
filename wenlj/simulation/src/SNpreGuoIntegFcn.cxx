#include <iostream>
#include <fstream>
#include <vector>
#include <array>

#include "TGraph.h"
#include "TVectorD.h"
#include "TSpline.h"
#include "TMath.h"
#include "SNpreGuoIntegFcn.hh"

SNpreGuoIntegFcn::SNpreGuoIntegFcn(){
    initialize();
}

SNpreGuoIntegFcn::SNpreGuoIntegFcn(int imode){
    initialize();
    readFluxGraph(imode);
}

void SNpreGuoIntegFcn::initialize(){
    for(int it=0; it<4; it++){
        timeMin[it] = -999;
        timeMax[it] = -999;
    }
    for(int ip=0; ip<80; ip++){
        grBinNum_nue_T[ip]     = NULL;
        grBinNum_antinue_T[ip] = NULL;
        grBinNum_nux_T[ip]     = NULL;
    }

    grNum_nux_T     = NULL;
    grNum_antinue_T = NULL;
    grNum_nux_T     = NULL;
    grNum_antinux_T = NULL;
    grLum_nux_T     = NULL;
    grLum_antinue_T = NULL;
    grLum_nux_T     = NULL;
    grLum_antinux_T = NULL;
    grNum_nux_E     = NULL;
    grNum_antinue_E = NULL;
    grNum_nux_E     = NULL;
    grNum_antinux_E = NULL;

}

SNpreGuoIntegFcn::~SNpreGuoIntegFcn(){
}


void SNpreGuoIntegFcn::readFluxGraph(int imode){
    std::ifstream fSNmod;
    int tag = (imode%10000)/10;
    fSNmod.open(Form("/junofs/users/miaoyu/supernova/wenlj/simulation/data/preSN/guo/S%d_flux_tE.dat",tag));
    if(!fSNmod){
        std::cout << "/junofs/users/miaoyu/supernova/wenlj/simulation//data/preSN/guo " << imode <<".data can't be opened!" << std::endl;
        exit(-1);
    }

    const int npnt_E = 80;//energy bin
    std::vector<std::array<double,npnt_E> > vecNum_nue;
    std::vector<std::array<double,npnt_E> > vecNum_antinue;
    std::vector<std::array<double,npnt_E> > vecNum_nux;
    std::vector<std::array<double,npnt_E> > vecNum_antinux;

    //reading from file
    double time;
    double factorTime = -3600*24;
    while(fSNmod>>time){
        //[/day]
        vecTime.push_back(time*factorTime);
        
        std::array<double, npnt_E> arrNum_nue;
        std::array<double, npnt_E> arrNum_antinue;
        std::array<double, npnt_E> arrNum_nux;
        std::array<double, npnt_E> arrNum_antinux;

        for(int ip=0; ip<npnt_E; ip++){
            //[/MeV/s/cm^2]
            fSNmod >> arrNum_nue[ip];
            //std::cout << arrNum_nue[ip] << std::endl;
        }
        fSNmod >> time;
        for(int ip=0; ip<npnt_E; ip++){
            //[/MeV/s/cm^2]
            fSNmod >> arrNum_antinue[ip];
        }
        fSNmod >> time;
        for(int ip=0; ip<npnt_E; ip++){
            //[/MeV/s/cm^2]
            fSNmod >> arrNum_nux[ip];
        }
        fSNmod >> time;
        for(int ip=0; ip<npnt_E; ip++){
            //[/MeV/s/cm^2]
            fSNmod >> arrNum_antinux[ip];
        }

        vecNum_nue.push_back(arrNum_nue);
        vecNum_antinue.push_back(arrNum_antinue);
        vecNum_nux.push_back(arrNum_nux);
        vecNum_antinux.push_back(arrNum_antinux);
    }

    //flux distribution
    //number and luminosity time distribution of each energy bin 
    int ntbin = vecTime.size();
    for(int ity=0; ity<4; ity++){
        timeMin[ity] = vecTime[0];
        timeMax[ity] = vecTime[ntbin-1];
        //std::cout << "Tmin " << timeMin[ity] << "; Tmax " << timeMax[ity] << std::endl;
    }
    double* binTime = new double[ntbin];
    double** ibinNum_nue_T     = new double*[npnt_E];
    double** ibinNum_antinue_T = new double*[npnt_E];
    double** ibinNum_nux_T     = new double*[npnt_E];
    double** ibinNum_antinux_T = new double*[npnt_E];

    for(int ib=0;ib<npnt_E;ib++){
        ibinNum_nue_T[ib]     = new double[ntbin];
        ibinNum_antinue_T[ib] = new double[ntbin];
        ibinNum_nux_T[ib]     = new double[ntbin];
        ibinNum_antinux_T[ib]     = new double[ntbin];
    }
    //for time integral
    double* timeNum_nue     = new double[ntbin];
    double* timeNum_antinue = new double[ntbin];
    double* timeNum_nux     = new double[ntbin];
    double* timeNum_antinux = new double[ntbin]; 
    double* timeLum_nue     = new double[ntbin];
    double* timeLum_antinue = new double[ntbin];
    double* timeLum_nux     = new double[ntbin];
    double* timeLum_antinux = new double[ntbin];


    // 1kpc = 3.086e21
    //from [/MeV/s/cm^2] to [/MeV/s] 
    double factorNum = 4*TMath::Pi()*TMath::Power((0.2*3.086e21),2);

    //from MeV to erg
    double factorLum = 1.60218e-6;
    for(int it=0; it<ntbin; it++){
        binTime[it] = vecTime[it];
        timeNum_nue[it]     = 0;
        timeNum_antinue[it] = 0;
        timeNum_nux[it]     = 0;
        timeNum_antinux[it] = 0;
        timeLum_nue[it]     = 0;
        timeLum_antinue[it] = 0;
        timeLum_nux[it]     = 0;
        timeLum_antinux[it] = 0;
        for(int ib=0; ib<npnt_E; ib++){
            ibinNum_nue_T[ib][it]     = factorNum*vecNum_nue[it].at(ib);
            ibinNum_antinue_T[ib][it] = factorNum*vecNum_antinue[it].at(ib);
            ibinNum_nux_T[ib][it]     = factorNum*vecNum_nux[it].at(ib);
            ibinNum_antinux_T[ib][it] = factorNum*vecNum_antinux[it].at(ib);

            if(ib<npnt_E-1){
                double deltaE = TMath::Power(10,0.02*(ib+1))-TMath::Power(10,0.02*ib);
                double iE = 0.5*(TMath::Power(10,0.02*(ib+1))+TMath::Power(10,0.02*ib));
                timeNum_nue[it]     += factorNum*deltaE*vecNum_nue[it].at(ib);
                timeNum_antinue[it] += factorNum*deltaE*vecNum_antinue[it].at(ib);
                timeNum_nux[it]     += factorNum*deltaE*vecNum_nux[it].at(ib); 
                timeNum_antinux[it] += factorNum*deltaE*vecNum_antinux[it].at(ib); 
               
                //Luminosity [erg/s]
                timeLum_nue[it]     += factorNum*factorLum*deltaE*iE*vecNum_nue[it].at(ib);
                timeLum_antinue[it] += factorNum*factorLum*deltaE*iE*vecNum_antinue[it].at(ib);
                timeLum_nux[it]     += factorNum*factorLum*deltaE*iE*vecNum_nux[it].at(ib);
                timeLum_nux[it]     += factorNum*factorLum*deltaE*iE*vecNum_antinux[it].at(ib);
            }
            //std::cout << "TimeNum " << timeNum_nue[it] << std::endl;
        }
    }

    //===========time distribution for number and luminosity===========//
    grNum_nue_T     = new TGraph(ntbin, binTime, timeNum_nue);
    grNum_antinue_T = new TGraph(ntbin, binTime, timeNum_antinue);
    grNum_nux_T     = new TGraph(ntbin, binTime, timeNum_nux);
    grNum_antinux_T = new TGraph(ntbin, binTime, timeNum_nux);

    grLum_nue_T     = new TGraph(ntbin, binTime, timeLum_nue);
    grLum_antinue_T = new TGraph(ntbin, binTime, timeLum_antinue);
    grLum_nux_T     = new TGraph(ntbin, binTime, timeLum_nux);
    grLum_antinux_T = new TGraph(ntbin, binTime, timeLum_nux);
    
    //time distribution for each energy bin
    for(int ib=0; ib<npnt_E; ib++){
        grBinNum_nue_T[ib]     = new TGraph(ntbin,binTime, ibinNum_nue_T[ib]);
        grBinNum_antinue_T[ib] = new TGraph(ntbin,binTime, ibinNum_antinue_T[ib]);
        grBinNum_nux_T[ib]     = new TGraph(ntbin,binTime, ibinNum_nux_T[ib]);
        grBinNum_antinux_T[ib] = new TGraph(ntbin,binTime, ibinNum_antinux_T[ib]);
    }
    //number distribution by integrating all the time
    double* binNum_nue     = new double[npnt_E];
    double* binNum_antinue = new double[npnt_E];
    double* binNum_nux     = new double[npnt_E];
    double* binNum_antinux = new double[npnt_E];
    for(int ib=0; ib<npnt_E; ib++){
        binE_E[ib] = TMath::Power(10,0.02*ib);
        //binNum_nue[ib]     = integIbinTgraph(grBinNum_nue_T[ib], timeMin[0],timeMax[0]);
        //binNum_antinue[ib] = integIbinTgraph(grBinNum_antinue_T[ib], timeMin[1],timeMax[1]);
        //binNum_nux[ib]     = integIbinTgraph(grBinNum_nux_T[ib], timeMin[2],timeMax[2]);
        //binNum_antinux[ib] = integIbinTgraph(grBinNum_antinux_T[ib], timeMin[3],timeMax[3]);
        double factor = 3600*24;
        binNum_nue[ib]     = integIbinTgraph(grBinNum_nue_T[ib], -1*factor,-1e-6*factor);
        binNum_antinue[ib] = integIbinTgraph(grBinNum_antinue_T[ib], -1*factor,-1e-6*factor);
        binNum_nux[ib]     = integIbinTgraph(grBinNum_nux_T[ib], -1*factor,-1e-6*factor);
        binNum_antinux[ib] = integIbinTgraph(grBinNum_antinux_T[ib], -1*factor,-1e-6*factor);
    }

    //=======influence distribution=======//
    grNum_nue_E     = new TGraph(npnt_E, binE_E, binNum_nue);
    grNum_antinue_E = new TGraph(npnt_E, binE_E, binNum_antinue);
    grNum_nux_E     = new TGraph(npnt_E, binE_E, binNum_nux);
    grNum_antinux_E = new TGraph(npnt_E, binE_E, binNum_antinux);

}


double SNpreGuoIntegFcn::integIbinTgraph(TGraph* igraph, double first, double last){
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
double SNpreGuoIntegFcn::fluxTimeE(double* x, double* par){
    double Ev = x[0];
    double tfirst = par[0];
    double tlast  = par[1];
    int itype = (int)(par[2]);

    double flux = 0;
    int npnt_E = 80;
    double* NumE = new double[npnt_E];
    if( itype==0 ){
        for(int ib=0; ib<npnt_E; ib++){
            NumE[ib] = integIbinTgraph(grBinNum_nue_T[ib],tfirst,tlast);
            //std::cout <<"ib" << ib <<"; NumE " << NumE[ib] << std::endl;
        }
        TGraph gr(npnt_E,binE_E, NumE);
        flux = gr.Eval(Ev);
    }
    if(itype == 1){
        for(int ib=0; ib<npnt_E; ib++){
            NumE[ib] = integIbinTgraph(grBinNum_antinue_T[ib],tfirst,tlast);
        }
        TGraph gr(npnt_E,binE_E, NumE);
        flux = gr.Eval(Ev);
    }
    if(itype ==2 || itype==4){
        for(int ib=0; ib<npnt_E; ib++){
            NumE[ib] = integIbinTgraph(grBinNum_nux_T[ib],tfirst,tlast);
        }
        TGraph gr(npnt_E,binE_E, NumE);
        flux = gr.Eval(Ev);
    }
    
    if(itype ==3 || itype==5){
        for(int ib=0; ib<npnt_E; ib++){
            NumE[ib] = integIbinTgraph(grBinNum_antinux_T[ib],tfirst,tlast);
        }
        TGraph gr(npnt_E,binE_E, NumE);
        //flux = gr.Eval(Ev,0,"S");
        flux = gr.Eval(Ev);
    }

    delete[] NumE;
    return flux;
}


double SNpreGuoIntegFcn::getFluxE(double E, int type){
    double flux = 0;
    if(type == 0) flux = grNum_nue_E->Eval(E); 
    if(type == 1) flux = grNum_antinue_E->Eval(E); 
    if(type ==2 || type==4)  flux = grNum_nux_E->Eval(E);
    if(type ==3 || type==5)  flux = grNum_antinux_E->Eval(E);
    return flux;
}

double SNpreGuoIntegFcn::getNumT(double T, int type){
    double flux = 0;
    //if(T < timeMin[0] || T>timeMax[0])return 0;
    if(T>timeMax[0])return 0;
    if(type == 0) {
        flux = grNum_nue_T->Eval(T); 
    }
    if(type == 1){
        flux = grNum_antinue_T->Eval(T); 
    }
    if(type == 2 || type==4){
        flux = grNum_nux_T->Eval(T);
    }
    if(type == 3 || type==5){
        flux = grNum_antinux_T->Eval(T);
    }
    return flux;
}
double SNpreGuoIntegFcn::getLumT(double T, int type){
    double lumi = 0;
    //if(T < timeMin[0] || T>timeMax[0])return 0;
    if(T>timeMax[0])return 0;
    if(type == 0){
        lumi = grLum_nue_T->Eval(T); 
    }
    if(type == 1){
        lumi = grLum_antinue_T->Eval(T); 
    }
    if(type ==2 || type==4){
        lumi = grLum_nux_T->Eval(T);
    }
    
    if(type ==3 || type==5){
        lumi = grLum_antinux_T->Eval(T);
    }
    return lumi;
}

double SNpreGuoIntegFcn::getEventAtTime(double time, double E, int type){
    double flux = 0; 
    int npnt_E = 80;
    double* NumE = new double[npnt_E];
    //if(time < timeMin[0] || time>timeMax[0])return 0;
    if(time>timeMax[0])return 0;
    if(type == 0){
        for(int ib=0; ib<npnt_E; ib++){
            NumE[ib] = grBinNum_nue_T[ib]->Eval(time);
        }
        TGraph gr(npnt_E, binE_E, NumE);
        flux = gr.Eval(E);
    }
    if(type == 1){
        for(int ib=0; ib<npnt_E; ib++){
            NumE[ib] = grBinNum_antinue_T[ib]->Eval(time, 0);
        }
        TGraph gr(npnt_E, binE_E, NumE);
        flux = gr.Eval(E);
    }
    if(type == 2 || type == 4){
        for(int ib = 0; ib<npnt_E; ib++){
            NumE[ib] = grBinNum_nux_T[ib]->Eval(time);
        }
        TGraph gr(npnt_E, binE_E, NumE);
        flux = gr.Eval(E);
    }
    
    if(type == 3 || type == 5){
        for(int ib = 0; ib<npnt_E; ib++){
            NumE[ib] = grBinNum_antinux_T[ib]->Eval(time);
        }
        TGraph gr(npnt_E, binE_E, NumE);
        flux = gr.Eval(E);
    }
    if(flux<0)flux=0;
    delete[] NumE;
    return flux;
}

double SNpreGuoIntegFcn::getAverageET(double time, int type){
    double aveE = 0;
    //if(time < timeMin[0] || time > timeMax[0])return 0;
    if(time > timeMax[0])return 0;
    //nue
    if(type == 0) {
        aveE = grLum_nue_T->Eval(time)/(grNum_nue_T->Eval(time));
    }

    //antinue
    if(type == 1){
        aveE = grLum_antinue_T->Eval(time)/(grNum_antinue_T->Eval(time));
    }
    
    //nux
    if(type ==2 || type==4){
        aveE = grLum_nux_T->Eval(time)/(grNum_nux_T->Eval(time));
    }
    
    //antinux
    if(type ==3 || type==5){
        aveE = grLum_antinux_T->Eval(time)/(grNum_antinux_T->Eval(time));
    }
    double index = 6.24151e5;
    return index*aveE;
}
