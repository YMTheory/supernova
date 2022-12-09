#include "SNeffectLS.hh"
#include "TMath.h"

SNeffectLS::SNeffectLS():
    // JUNO Configurations:
    fEres(3e-2){
    double mLS = 20.;//kt
    double fracP = 12e-2;
    double fracC = 88e-2;
    double Na = 6.02e23;
    double molarH = 1.;//g/mol
    double molarC = 12.;
    fNumP = mLS*1e9*fracP*Na/molarH;
    fNumC = mLS*1e9*fracC*Na/molarC;
    fNumO = 0;
    fEthr = 0.2;//MeV


    // LENA Configurations
    /// double mLS = 50.;//kt
    /// double fracP = 12e-2;
    /// double fracC = 88e-2;
    /// double Na = 6.02e23;
    /// double molarH = 1.;//g/mol
    /// double molarC = 12.;
    /// fNumP = mLS*1e9*fracP*Na/molarH;
    /// fNumC = mLS*1e9*fracC*Na/molarC;
    /// fEthr = 0.25;//MeV
    
    // THEIA-25 Configurations
    //double mLS = 25.;//kt
    //double fracH2O = 0.9;
    //double fracLAB = 0.1;
    //double Na = 6.02e23;
    //double fracP0 = 12.5e-2;
    //double fracO = 1 - fracH2O;
    //double fracP = 12e-2;
    //double fracC = 88e-2;
    //double molarH = 1.;//g/mol
    //double molarC = 12.;
    //double molarO = 16.;
    //fNumP = mLS*1e9*(fracH2O*fracP0+fracLAB*fracP)*Na/molarH;
    //fNumC = mLS*1e9*fracC*fracLAB*Na/molarC;
    //fNumO = mLS*1e9*fracH2O*fracO*Na/molarO;
    //fEthr = 0.60;//MeV
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
