#include "SNnupLS.hh"
#include "TMath.h"

SNnupLS::SNnupLS():SNeffectLS(){
}
SNnupLS::~SNnupLS(){

}

double SNnupLS::differentialXS_old(double E, double T, int type){
    double mp = 938;//MeV
    //double Emin = TMath::Sqrt(mp*T/2);//with approximation
    double Emin = TMath::Sqrt((1+T/2./mp)*(mp*T/2))+T/2.;
    if(E<Emin)return 0;
    double diffxs = (1+466*T/(E*E))*4.83*1e-42;

    return diffxs;

}

double SNnupLS::totalXS_old(double E, int type){
    double mp = 938;
    double totxs = (2./mp+233*4/(mp*mp))*E*E*4.83*1e-42;
    return totxs;
    //return 0;
}


double SNnupLS::differentialXS(double E, double T, int type) {
    //nu-P scattering parameter
    double sw2 = 0.23126;
    double hbarc = 1.973269718e-11;	//Mev * cm
    double gf = 1.1663787e-11;	//Mev^2
    double mp = 938.272046;	//MeV
    double mvv = 840.0;	//MeV
    double kappap = 1.792847356;
    double kappan = -1.9130427;
    //double mav = 1032.0;	//MeV
    double ga = -1.267;
    double alpha = 1.0 - 2.0 * sw2;
    double beta = -1.0;
    double gamma = - 2.0 / 3.0 * sw2;
    double constt = mp * pow(gf * hbarc, 2) / M_PI;
    //Input
    double eps = E / mp;
    double tau = T / 2.0 / mp;
    if(eps <= 0.0 || tau <= 0.0 || tau > pow(eps, 2) / (1 + 2 * eps))
    {
        return 0.0;
    }

    double qf3v = 0.5 * (kappap - kappan) / ((1 + tau) * pow(1 + 4 * tau * pow(mp / mvv, 2), 2));
    double qf0v = 1.5 * (kappap + kappan) / ((1 + tau) * pow(1 + 4 * tau * pow(mp / mvv, 2), 2));
    double qg3v = 0.5 * (kappap - kappan + 1) / pow(1 + 4 * tau * pow(mp / mvv, 2), 2);
    double qg0v = 1.5 * (kappap + kappan + 1) / pow(1 + 4 * tau * pow(mp / mvv, 2), 2);

    double qf2 = alpha * qf3v + gamma * qf0v;
    double qf1 = alpha * qg3v + gamma * qg0v - qf2;

    double qga3 = 0.5 * ga / pow(1 + 4 * tau * pow(mp / mvv, 2), 2);
    double qg = qga3 * beta;

    double funca = 4 * tau * qf1 * qf2 
        + (1 - tau) * (tau * pow(qf2, 2) - pow(qf1, 2)) + (1 + tau) * pow(qg, 2);
    double funcb = 4 * qg * (qf1 + qf2);
    double funcc = pow(qf1, 2) + tau * pow(qf2, 2) + pow(qg, 2);

    if (type == 0 or type == 2 or type == 4)
        // neutrino
        return constt * (pow(eps - tau, 2) * funcc 
            +  (eps - tau) * tau * funcb 
            + tau * funca) / pow(eps, 2);
    else if (type == 1 or type ==3 or type == 5)
        // anti-neutrino
        return constt * (pow(eps - tau, 2) * funcc 
            - (eps - tau) * tau * funcb 
            + tau * funca) / pow(eps, 2);

}

double SNnupLS::totalXS(double E, int type) {
    return 0.0;
}



