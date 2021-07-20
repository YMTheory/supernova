#include "SNcncLS.hh"
#include "TGraph.h"
#include "TMath.h"

SNcncLS::SNcncLS():SNeffectLS(){
}

SNcncLS::~SNcncLS(){

}

double SNcncLS::differentialXS(double E, double T, int type){
    //double E_thr = 1.806;//MeV
    return 0;
}

double SNcncLS::totalXS(double E, int type){

    double inputenu[17] = {16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100};
    double sigmacnu[17] = {0.010, 0.106, 0.302, 0.599, 0.994, 1.49, 2.07, 2.74, 4.78, 7.26, 10.1, 13.1, 19.5, 25.4, 30.2, 33.7, 35.8};
    double sigmacnubar[17] = {0.0095, 0.099, 0.279, 0.547, 0.896, 1.32, 1.82, 2.38, 4.03, 5.95, 8.03, 10.2, 14.4, 17.9, 20.7, 22.5, 23.6};
    if(E<15.11) return 0;
    if(type%2==0){  
        TGraph* grcnc = new TGraph(17,inputenu, sigmacnu);
        double totxs = grcnc->Eval(E, 0, "S")*1e-42;
        delete grcnc;
        return totxs;
    } 
    if(type%2==1){  
        TGraph* grcnc = new TGraph(17,inputenu, sigmacnubar);
        double totxs = grcnc->Eval(E, 0, "S")*1e-42;
        delete grcnc;
        return totxs;
    }
    return 0;
}
