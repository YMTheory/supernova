#ifndef _SNbccLS_hh
#define _SNbccLS_hh

//-----------------------------------------------------------------------
//class for "
//SN anti_nu_e-proton inverse beta decay in LS
//"
//Created by Huiling, 2016-12-21
//-----------------------------------------------------------------------


#include "SNeffectLS.hh"

class SNbccLS : public SNeffectLS{
    public:
        SNbccLS();
        ~SNbccLS();

        virtual double differentialXS(double E, double T, int type);
        virtual double totalXS(double E, int type);

        public:
        ClassDef(SNbccLS, 1);
};


#endif
