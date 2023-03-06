#ifndef _SNeffectLS_hh
#define _SNeffectLS_hh

//-------------------------------------------------------------------
//SN effect of Liquid Scintillator abstract class
//
//Created by Huiling, 2016-12-21
//-------------------------------------------------------------------

#include <iostream>
#include "TObject.h"
#include "TGraph.h"


class SNeffectLS : public TObject{
    public:
        SNeffectLS();
        virtual ~SNeffectLS();

        //Input:
        //-----E: neutrino energy
        //-----T:recoiled kinetic energy
        //Output:
        //differential cross section 
        //unit: cm^2MeV^-1

        //type: 0,1->nu_e,   anti_nu_e
        //      2,3->nu_mu,  anti_nu_mu
        //      4,5->nu_tau, anti_nu_tau
        virtual double differentialXS(double E, double T, int type) = 0;
        virtual double totalXS(double E, int type) = 0;


        void setERelativeRes(double res)  { fEres = res;}
        void setNumberOfProton(double np) { fNumP = np; }
        void setNumberOfCarbon(double nc) { fNumC = nc; }
        void setNumberOfOxygen(double no) { fNumO = no; }
        void setThresholdE(double E)      { fEthr = E;  }

        double getNumberOfProton()  { return fNumP; }
        double getNumberOfCarbon()  { return fNumC; }
        double getNumberOfOxygen()  { return fNumO; }
        double getERelativeResPar() { return fEres; }
        double getThresholdE()      { return fEthr; }
        double getProbResGauss(double Eobs, double Evis);
        double getEres(double Evis);
        double getElecNonl(double E);
        double getElecTfromEvis(double T);
        double getElecRes(double Eobs, double Evis);
        double getPosiNonl(double E);
        double getPosiRes(double Eobs, double Evis);
        double getPosiTfromEvis(double T);

    protected:
        double fEres;// relative energy resolution
        double fNumP;//number of free protons in LS
        double fNumC;//number of free protons in LS
        double fEthr;//threshold energy 
        double fNumO;//number of free oxygen in LS

        TGraph* gElecNonl;
        TGraph* gElecNonlInv;
        TGraph* gElecRes;
        TGraph* gPosiNonl;
        TGraph* gPosiRes;
        TGraph* gPosiNonlInv;

    //public:
    //    ClassDef(SNeffectLS,1);

};
#endif

