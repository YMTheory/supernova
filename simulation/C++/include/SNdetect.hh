#ifndef _SNdetect_hh
#define _SNdetect_hh

#include "TObject.h"

class TRandom3;
class SNeffectLS;
class SNsource;
class TF1;
class TF2;
class TF3;
class TFile;
class TGraph;

class SNdetect : public TObject{
    private:   
        SNdetect();
    public:
       enum channelName {NuP, NuE, IBD};
       
       ~SNdetect();

       static SNdetect* instance();

       void initChannel(channelName cha);
       void setChannel(channelName cha){ initChannel(cha);}
       void setSrcModel(int imode);
       channelName getChannel(){ return fchannel; }

       //pointer
       SNsource* getPointerSrc(){ return psrc; }
       SNeffectLS* getPointerEffectLS(){ return peffLS;}

       //toyMC information
       void initFCN();//fcn for random generation
       void deletFCN();//delete pointers when exiting 
       void resetFCN();//after change parameters in source or LS, should reset fcn
       void setMCflavor(int type, int MH);
       double getEvisFromT(double T);
       double getTFromEvis(double Evis);
       double getTFromEv(double Ev);
       double getEvFromT(double T);
       double getTmax();
      

       //-----------------------------------------------------//
       //---------time-denpendent spectra distribution--------//
       //number of events
       //different mass hierarchy of neutrinos, 0 no osc; 1 NH; 2 IH
       double getXSweightedEvSpectrumAtTime(double time, double Ev,int type, int MH);
       double getTSpectrumAtTime(double time, double T, int type, int MH);
       double getEvisSpectrumAtTime(double time, double Evis, int type, int MH);
       double getEventAboveEthrVisAtTime(double time,double Ethr, int type, int MH);//number of events in LS w/o energy resolution at a specific time
      
    private:
       static SNdetect* pdet ;
       SNeffectLS*  peffLS;
       SNsource*    psrc; 
       channelName fchannel;//specific channel in this detection
       int iMCflavor;
       double fEvmin;
       double fEvmax;
       double fTmin;
       double fTmax;
       double fEthr;
       //for quenching of Nu-P ES channel
       TFile* file;
       TGraph* grquench;
       TGraph* grquenchInv;

       TF2* f2rndEvT;
       TF2* f2rndFlatEvT;
       TF1* f1rndEv;

};



#endif
