#ifndef _SNpreKatoIntegFcn_hh
#define _SNpreKatoIntegFcn_hh
#include <vector>

class TGraph;
class TString;

class SNpreKatoIntegFcn{
    public:
        SNpreKatoIntegFcn();
        explicit SNpreKatoIntegFcn(int imode);
        ~SNpreKatoIntegFcn();

        void getTimeLimits(double& tmin, double& tmax, int type){
            if(type==2 || type==3)type=2;
            if(type==4 || type==5)type=3;
            tmin = timeMin[type];
            tmax = timeMax[type];
        }

        void readFile(TString path, int itype);
        
        std::vector<double> getTimeVector(){ return vecTime;}
        double getFluxE(double E, int type);
        double getNumT(double T, int type);
        double getLumT(double T, int type);
        double getEventAtTime(double time, double E, int type);
        double getAverageET(double time, int type);
        double fluxTimeE(double* x, double* par);

    private:
        void initialize();
        void readFluxGraph(int imode);
        double integIbinTgraph(TGraph* igraph, double first, double last);//to integrate time distribution
        std::vector<double> vecTime;
        double timeMin[4];
        double timeMax[4];

        double binE_E[501];
        //time distribution of each energy bin
        //0 nu_e; 1 nu_e_bar; 2 nu_x; 3 nu_x_bar
        TGraph* grBinNum_T[4][501];

        //time distribution for all bins
        TGraph* grNum_T[4];
        TGraph* grLum_T[4];

        //energy distribution integrating all the time
        TGraph* grNum_E[4];
};

#endif
