#include "SNnupLS.hh"
#include "TMath.h"

SNnupLS::SNnupLS():SNeffectLS(){
}
SNnupLS::~SNnupLS(){

}

double SNnupLS::differentialXS(double E, double T, int type){
    double mp = 938;//MeV
    //double Emin = TMath::Sqrt(mp*T/2);//with approximation
    double Emin = TMath::Sqrt((1+T/2./mp)*(mp*T/2))+T/2.;
    if(E<Emin)return 0;
    double diffxs = (1+466*T/(E*E))*4.83*1e-42;

    return diffxs;

}

double SNnupLS::totalXS(double E, int type){
    //double mp = 938;
    //double totxs = (2./mp+233*4/(mp*mp))*E*E*4.83*1e-42;
    return 0;
}


