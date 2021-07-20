#include "SNeffectLS.hh"
#include "TMath.h"

SNeffectLS::SNeffectLS():
    fEres(3e-2){
    double mLS = 20.;//kt
    double fracP = 12e-2;
    double fracC = 88e-2;
    double Na = 6.02e23;
    double molarH = 1.;//g/mol
    double molarC = 12.;
    fNumP = mLS*1e9*fracP*Na/molarH;
    fNumC = mLS*1e9*fracC*Na/molarC;
    fEthr = 0.2;//MeV
}
SNeffectLS::~SNeffectLS(){

}
double SNeffectLS::getEres(double Evis){
    return fEres*TMath::Sqrt(Evis);
}

double SNeffectLS::getProbResGauss(double Eobs,double Evis){
    double sigmaE = getEres(Evis);
    double prob   = TMath::Gaus(Eobs, Evis, sigmaE, 1);

    return prob;
}
