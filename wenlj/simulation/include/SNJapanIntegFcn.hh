#ifndef _SNJapanIntegFcn_hh
#define _SNJapanIntegFcn_hh
#include <vector>
class TGraph;

class SNJapanIntegFcn{
    public:
        SNJapanIntegFcn();
        explicit SNJapanIntegFcn(int imode);
        ~SNJapanIntegFcn();

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

        double getFluxE(double E, int type);
        double getNumT(double T, int type);
        double getLumT(double T, int type);
        double getEventAtTime(double time, double E, int type);
        double getAverageET(double time, double type);
        double fluxTimeE(double* x, double* par);

    private:
        void readFluxGraph(int imode);
        double integIbinTgraph(TGraph* igraph, double first, double last);//to integrate time distribution
        double timeMin[3];
        double timeMax[3];
    
        int nbin_E;
        double* binE_E;//energy points
        //time distribution of each energy bin
        std::vector<TGraph*> grBinNum_nue_T;
        std::vector<TGraph*> grBinNum_antinue_T;
        std::vector<TGraph*> grBinNum_nux_T;

        //time distribution for all bins
        TGraph* grNum_nue_T;
        TGraph* grNum_antinue_T;
        TGraph* grNum_nux_T;

        TGraph* grLum_nue_T;
        TGraph* grLum_antinue_T;
        TGraph* grLum_nux_T;
        //energy distribution integrating all the time
        TGraph* grNum_nue_E;
        TGraph* grNum_antinue_E;
        TGraph* grNum_nux_E;
};

#endif
