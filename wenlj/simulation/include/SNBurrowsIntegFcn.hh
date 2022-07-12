#ifndef _SNBurrowsIntegFcn_hh
#define _SNBurrowsIntegFcn_hh

#include <vector>
class TGraph;
class SNBurrowsIntegFcn{
    public:
        SNBurrowsIntegFcn();
        explicit SNBurrowsIntegFcn( int imode);
        ~SNBurrowsIntegFcn();
 
        void getTimeLimits(double& tmin, double& tmax, int type){
            if(type < 2){
                tmin = timeMin[type];
                tmax = timeMax[type];
            }
            if(type >= 2){
                tmin = timeMin[2];
                tmax = timeMax[2];
            }
        }
        
        std::vector<double> getTimeVector() {return vecTime;}
        double fcnfluxtime(double* x, double* par);
        double getNumT(double time, int type);
        double getLumT(double time, int type);
        double getEventAtTime(double time, double E, int type);
        double getAverageET(double time, int type);
        double getEventAtTimeNum(double time, double E, int type);
    private:
        double timeMin[3];
        double timeMax[3];
        std::vector<double> vecTime;

        TGraph*   grLuminosity[3];//unnormalized energy spectrum
        TGraph*   grAverageE[3];//time evolution for average energy
        TGraph*   grAlpha[3];
        void readFluxGraph(int imode);
        void readNumFlux(int imass);

        double etspec_nua[80][100];
        double etspec_nue[80][100];
        double etspec_nux[80][100];
};

#endif
