#ifndef _SNnueLS_hh
#define _SNnueLS_hh

//----------------------------------------------------------------------
//class for "
//SN neutrino-electron elastic scattering cross section 
//"
//Created by Huiling, 2016-12-21
//----------------------------------------------------------------------


#include "SNeffectLS.hh"

class SNnueLS : public SNeffectLS{
    public:
        SNnueLS();
        ~SNnueLS();

        virtual double differentialXS(double E, double T, int type);
        virtual double totalXS(double E, int type);
        
    public:
        ClassDef(SNnueLS, 1);
};



#endif
