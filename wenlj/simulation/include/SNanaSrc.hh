#ifndef _SNanaSrc_hh
#define _SNanaSrc_hh

//-----------------------------------------------------------------
//SN neutrino fluence over a burst time deltaT
//
//Created by Huiling, 2016-12-22
//------------------------------------------------------------------

#include "SNsource.hh"

class SNanaSrc: public SNsource{

    public:
        SNanaSrc();
        ~SNanaSrc();

        //Input:
        //-----E: the energy of all neutrino flavor
        //-----Ea: burst energy for one neutrino flavor
        //-----type: specify the kind of neutrino,i.e. averagE[type]
        //Output:
        //-----unit: cm^-2*Mev^-1
        //type: 0,1->nu_e,   anti_nu_e
        //      2,3->nu_mu,  anti_nu_mu
        //      4,5->nu_tau, anti_nu_tau
        virtual double oneSNFluenceDet(double E, int type);
        //mass hierarchy: 0 no osc; 1 NH; 2 IH;
        virtual double oneSNFluenceDet(double E, int type, int MH);
        
        //Input:
        //-----E: the energy of the neutrino
        //Output:
        //-----unit: cm^-2*Mev^-1
        //mass hierarchy: 0 no osc; 1 NH; 2 IH;
        virtual double totalSNFluenceDet(double E);
        virtual double totalSNFluenceDet(double E, int MH);
        virtual double totalSNFluenceDetTimeIntegE(double time, int MH);
        virtual double oneSNFluenceDetTimeIntegE(double time, int type, int MH); 
         
        
        virtual double oneSNFluenceDetAtTime(double time, double E, int type, int MH);
        virtual double totalSNFluenceDetAtTime(double time, double E, int MH);

        virtual double oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type);
        virtual double totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast);
        virtual double oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type, int MH);
        virtual double totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int MH); 
        
        virtual double oneSNLuminosityTime(double time, int type);
        virtual double oneSNAverageETime(double time, int type);
        
        virtual void getTimeRange(double& tmin, double& tmax, int type);
        
    public:
        ClassDef(SNanaSrc,1);
};


#endif
