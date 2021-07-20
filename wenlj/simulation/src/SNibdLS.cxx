#include "SNibdLS.hh"
#include "TMath.h"

SNibdLS::SNibdLS():SNeffectLS(){
}

SNibdLS::~SNibdLS(){

}

double SNibdLS::differentialXS(double E, double T, int type){
    //double E_thr = 1.806;//MeV
    return 0;
}

double SNibdLS::totalXS(double E, int type){
    if(type == 1){
        double E_thr = 1.81;//MeV
        if(E < E_thr)return 0;

        double mn = 939.57;//neutron mass
        double mp = 938.27;//proton mass
        double deltam_np = mn-mp;//MeV
        double me = 0.51;//MeV 
        double Ee = E-deltam_np;//MeV
        double Pe = TMath::Sqrt(Ee*Ee-me*me);
        double totxs = 0.0952*(Ee*Pe)*1e-42;//cm^2

        return totxs; 
    }
    return 0;
}
