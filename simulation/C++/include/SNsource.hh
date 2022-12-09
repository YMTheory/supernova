#ifndef _SNsource_hh
#define _SNsource_hh

//-----------------------------------------------------------------
//SN neutrino fluence over a burst time deltaT
//
//Created by Huiling, 2016-12-22
//------------------------------------------------------------------

#include "TObject.h"

class SNsource: public TObject{

    public:
        SNsource();
        virtual ~SNsource();

        //Input:
        //-----E: the energy of all neutrino flavor
        //-----Ea: burst energy for one neutrino flavor
        //-----type: specify the kind of neutrino,i.e. averagE[type]
        //Output:
        //-----unit: cm^-2*Mev^-1
        //type: 0,1->nu_e,   anti_nu_e
        //      2,3->nu_mu,  anti_nu_mu
        //      4,5->nu_tau, anti_nu_tau
        virtual double oneSNFluenceDet(double E, int type) = 0;
        virtual double oneSNFluenceDet(double E, int type, int MH) = 0;
        
        //Input:
        //-----E: the energy of the neutrino
        //Output:
        //-----unit: cm^-2*Mev^-1
        //mass hierarchy: 0 no osc; 1 NH; 2 IH;
        virtual double totalSNFluenceDet(double E) = 0;
        virtual double totalSNFluenceDet(double E, int MH) = 0;
        virtual double totalSNFluenceDetTimeIntegE(double time, int MH) = 0;
        virtual double oneSNFluenceDetTimeIntegE(double time, int type, int MH) = 0; 

        virtual double totalSNFluenceDetAtTime(double time, double E, int MH) = 0;
        virtual double oneSNFluenceDetAtTime(double time, double E, int type, int MH) = 0;
        
        virtual double oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type, int MH) = 0;
        virtual double totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int MH) = 0;
        virtual double oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type) = 0;
        virtual double totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast) = 0;
        
        virtual double oneSNLuminosityTime(double time, int type) = 0;
        virtual double oneSNAverageETime(double time, int type) = 0;
        
        virtual void getTimeRange(double& tmin, double& tmax, int type) = 0;


        void setSNTotalEnergy(double E)  { totE = E; }
        void setSNDistance(double d)     { dist = d; }
        void setAverageEnergy(double E[]) { 
            averagE[0] = E[0];
            averagE[1] = E[1];
            averagE[2] = E[2];
            averagE[3] = E[3];
            averagE[4] = E[4];
            averagE[5] = E[5];
            
        }
        void setSNmaxEv(double Emax) { fEvmax = Emax; }
        void setSNminEv(double Emin) { fEvmin = Emin; }
        
        
        double  getSNmaxEv()      { return fEvmax; }  
        double  getSNminEv()      { return fEvmin; }  
        double  getSNTotalEnergy(){ return totE; }
        double  getSNDistance()   { return dist; }
        double* getAverageE()     { return averagE; }
        double* getEnergyForEachFlavor() {
            for(int it=0; it<6; it++){
                Ea[it] = totE/6;
            }
            return Ea;
        };
        

        double totE;//SN burst total energy over time deltaT; unit:1erg
        double dist;//supernova source distance; unit:1kpc
        double fEvmin;
        double fEvmax;
        
        double Ea[6];//enenrgy for each flavor neutrino
        double averagE[6];//average energy for six flavor neutrinos; unit: MeV
        // the sequence of the 6 neutrino flavors: nu_e, anti_nu_e,nu_mu, anti_nu_mu, nu_tau, anti_nu_tau

    //public:
    //    ClassDef(SNsource,1);
};


#endif
