#include "SNeffectLS.hh"
#include "TMath.h"
#include "TFile.h"

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


    TFile* f = new TFile("/junofs/users/miaoyu/supernova/simulation/C++/detector/energy_response_beta.root", "read");
    std::cout << "+++++> Load energy nonlinearity and energy resolution for electrons and positrons here !" << std::endl;
    gElecNonl = (TGraph*)f->Get("nonle1");
    gElecNonlInv = (TGraph*)f->Get("nonle1inv");
    gElecRes  = (TGraph*)f->Get("rese1");
    gPosiNonl = (TGraph*)f->Get("nonle2");
    gPosiNonlInv = (TGraph*)f->Get("nonle2inv");
    gPosiRes  = (TGraph*)f->Get("rese2");
    // test line:
    /// std::cout << "Electron energy nonlinearity test " << "\n"
    ///           << " Edep = " << 1 << " MeV, Evis/Edep = " << gElecNonl->Eval(1) << "\n"
    ///           << "Electron inverse nonlinearity test " << "\n"
    ///           << "Evis = " << 1 << " MeV, Edep = " << gElecNonlInv->Eval(1) << " MeV \n"
    ///           << "Electron energy resoltuion test " << "\n"
    ///           << " Evis = " << 1 << " MeV, sigma/Evis = " << gElecRes->Eval(1) 
    ///           << std::endl;
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

double SNeffectLS::getElecNonl(double E) {
    return gElecNonl->Eval(E);
}

double SNeffectLS::getElecRes(double Eobs, double Evis) {
    double sigma = Evis * gElecRes->Eval(Evis) ;
    double prob = TMath::Gaus(Eobs, Evis, sigma, 1);
    return prob;
}

double SNeffectLS::getPosiNonl(double E) {
    return gPosiNonl->Eval(E);
}

double SNeffectLS::getPosiRes(double Eobs, double Evis) {
    double sigma = Evis * gPosiRes->Eval(Evis) ;
    double prob = TMath::Gaus(Eobs, Evis, sigma, 1);
    if (prob < 0) {
        std::cout << "Here the probability<0???"
                  << ", Evis=" << Evis 
                  << ", Eobs=" << Eobs
                  << ", sigma=" << sigma
                  << ", probability=" << prob
                  << std::endl;
    }
    return prob;
}

double SNeffectLS::getElecTfromEvis(double Evis) {
    return gElecNonlInv->Eval(Evis);
}

double SNeffectLS::getPosiTfromEvis(double Evis) {
    return gPosiNonlInv->Eval(Evis);
}

