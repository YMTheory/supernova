#include "SNCEvNSLS.hh"
#include "TMath.h"

SNCEvNSLS::SNCEvNSLS() : SNeffectLS() {}

SNCEvNSLS::~SNCEvNSLS() {}

double SNCEvNSLS::differentialXS(double E, double T, int type) {
    //J. Phys. G, Nucl. Part. Phys. 39, 095204
    //return differential cross section;
    //Ev unit: MeV
    //XS unit: cm2/MeV
    
    if (T == 0 || E == 0)
        return 0;
    
    double GF = 1.1664e-5;
    double NA = 6.022e23;
    double cmtoeVminus1 = 51000;
    double eVminus1tocm = 1./cmtoeVminus1;
    double gtoeV = 5.61e32;
    double mC12 = 12 / NA * gtoeV;
    double q = TMath::Sqrt(T*T + 2*mC12*1e-6*T);
    int Z = 6;
    int N = 6;

    // form factor
    q = q * 1e6; // eV
    int A = Z + N;
    double R = 1.2 * TMath::Power(A, 1./3.) * 1e-13; // cm
    double s = 1e-13; // cm
    double R0 = TMath::Sqrt(R*R - 5*s*s);
    double x = q * R0 * cmtoeVminus1;
    double j1 = TMath::Sin(x) / x / x - TMath::Cos(x)/x;
    double formFactor = 3 * j1 / x * TMath::Exp(-0.5*(q*s*cmtoeVminus1) * (q*s*cmtoeVminus1));

    double xs = GF*GF * (mC12*1e-9) / (8*TMath::Pi()) * ((4*0.2325-1)*Z+N)*((4*0.2325-1)*Z+N) * (1 + (1-T/E)*(1-T/E) - mC12*1e-6*T/E/E) * formFactor*formFactor * 1e-18 * eVminus1tocm*eVminus1tocm * 1e-3;

    if (xs < 0)
        return 0;
    else
        return xs;

}

double SNCEvNSLS::totalXS(double E, int type) {
    return 0;

}


