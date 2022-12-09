#ifndef _SNGarchingIntegFcn_hh
#define _SNGarchingIntegFcn_hh
//to integral of time for Garching model
#include <vector>
class TGraph;

class SNGarchingIntegFcn{
    public:
    SNGarchingIntegFcn();
    //For Garching Model, naming with four digits
    //8+(integ(solar mass))+(0, c 1, o 2, co 3) LS EoS
    //9+(integ(solar mass))+(0, c 1, o 2, co 3) Shen EoS
    explicit SNGarchingIntegFcn(int imode);
    ~SNGarchingIntegFcn();

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
    std::vector<double> getTimeVector(){ return vecTime;}
    double fcnfluxtime(double* x, double* par);//time distribution at one energy
    double getNumT(double time, int type);
    double getLumT(double time, int type);
    double getEventAtTime(double time, double E, int type);
    double getAverageET(double time, int type);
    private:
    double timeMin[3];
    double timeMax[3];
    std::vector<double> vecTime;
    TGraph* grLuminosity[3];
    TGraph* grAlpha[3];
    TGraph* grAverageE[3];
    void readFluxGraph(int imode);

};

#endif
