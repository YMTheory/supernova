#ifndef _SNnupLS_hh
#define _SNnupLS_hh

//----------------------------------------------------------------------
//class for "
//SN neutrino-proton elastic scattering in LS
//"
//Created by Huiling, 2016-12-21
//-----------------------------------------------------------------------


#include "SNeffectLS.hh"

class SNnupLS : public SNeffectLS{
    public: 
        SNnupLS();
        ~SNnupLS();

        virtual double differentialXS(double E, double T, int type);
        virtual double totalXS(double E, int type);


    public:
        ClassDef(SNnupLS, 1);
}; 

#endif

