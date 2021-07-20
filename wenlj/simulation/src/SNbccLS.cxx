#include "SNbccLS.hh"
#include "TGraph.h"
#include "TMath.h"

SNbccLS::SNbccLS():SNeffectLS(){
}

SNbccLS::~SNbccLS(){

}

double SNbccLS::differentialXS(double E, double T, int type){
    //double E_thr = 1.806;//MeV
    return 0;
}

double SNbccLS::totalXS(double E, int type){
    if(type == 1){
        if(E<14.39)return 0;
        double inputenu[17] = {16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100};
        double sigmabcc[17] = {0.086, 0.327, 0.711, 1.23, 1.87, 2.62, 3.48, 4.42, 7.10, 10.1, 13.2, 16.4, 22.2, 27.0, 30.5, 32.8, 34.2};
        TGraph* grbcc = new TGraph(17,inputenu, sigmabcc);
        double totxs = grbcc->Eval(E, 0, "S")*1e-42;
        delete grbcc;
        return totxs; 
    }
    return 0;
}
