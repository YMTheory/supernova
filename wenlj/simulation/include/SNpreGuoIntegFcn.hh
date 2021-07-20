#ifndef _SNpreGuoIntegFcn_hh
#define _SNpreGuoIntegFcn_hh
#include <vector>
class TGraph;

class SNpreGuoIntegFcn{
    public:
        SNpreGuoIntegFcn();
        explicit SNpreGuoIntegFcn(int imode);
        ~SNpreGuoIntegFcn();

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

        double binE_E[80];
        //time distribution of each energy bin
        TGraph* grBinNum_nue_T[80];
        TGraph* grBinNum_antinue_T[80];
        TGraph* grBinNum_nux_T[80];
        TGraph* grBinNum_antinux_T[80];

        //time distribution for all bins
        TGraph* grNum_nue_T;
        TGraph* grNum_antinue_T;
        TGraph* grNum_nux_T;
        TGraph* grNum_antinux_T;

        TGraph* grLum_nue_T;
        TGraph* grLum_antinue_T;
        TGraph* grLum_nux_T;
        TGraph* grLum_antinux_T;
        //energy distribution integrating all the time
        TGraph* grNum_nue_E;
        TGraph* grNum_antinue_E;
        TGraph* grNum_nux_E;
        TGraph* grNum_antinux_E;
};

#endif
