#ifndef SNCEvNSLS_hh
#define SNCEvNSLS_hh

#include "SNeffectLS.hh"

class SNCEvNSLS : public SNeffectLS {

    public:
        SNCEvNSLS();
        ~SNCEvNSLS();

        virtual double differentialXS(double E, double T, int type);
        virtual double totalXS(double E, int type);

};

#endif


