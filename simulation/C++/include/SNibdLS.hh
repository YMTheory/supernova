#ifndef _SNibdLS_hh
#define _SNibdLS_hh

//-----------------------------------------------------------------------
//class for "
//SN anti_nu_e-proton inverse beta decay in LS
//"
//Created by Huiling, 2016-12-21
//-----------------------------------------------------------------------


#include "SNeffectLS.hh"

class SNibdLS : public SNeffectLS{
    public:
        SNibdLS();
        ~SNibdLS();

        virtual double differentialXS(double E, double T, int type);
        virtual double totalXS(double E, int type);

    //public:
    //    ClassDef(SNibdLS, 1);
};


#endif
