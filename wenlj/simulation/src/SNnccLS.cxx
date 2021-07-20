#include "SNnccLS.hh"
#include "TGraph.h"
#include "TMath.h"

SNnccLS::SNnccLS():SNeffectLS(){
}

SNnccLS::~SNnccLS(){

}

double SNnccLS::differentialXS(double E, double T, int type){
    //double E_thr = 1.806;//MeV
    return 0;
}

double SNnccLS::totalXS(double E, int type){
    if(type == 0){
        if(E<17.34) return 0;
        double inputenu[17] = {16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100};
        double sigmancc[17] = {0, 0.036, 0.287, 0.772, 1.49, 2.44,3.62, 5.03, 9.47, 15.1, 21.8, 29.2, 45.2, 60.8, 74.2, 84.2, 90.6};
        TGraph* grncc = new TGraph(17,inputenu, sigmancc);
        double totxs = grncc->Eval(E, 0, "S")*1e-42;
        delete grncc;
        return totxs; 
    }
    return 0;
}
