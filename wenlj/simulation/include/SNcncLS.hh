#ifndef _SNcncLS_hh
#define _SNcncLS_hh

//-----------------------------------------------------------------------
//class for "
//SN anti_nu_e-proton inverse beta decay in LS
//"
//Created by Huiling, 2016-12-21
//-----------------------------------------------------------------------


#include "SNeffectLS.hh"

class SNcncLS : public SNeffectLS{
    public:
        SNcncLS();
        ~SNcncLS();

        virtual double differentialXS(double E, double T, int type);
        virtual double totalXS(double E, int type);

        public:
        ClassDef(SNcncLS, 1);
};


#endif
