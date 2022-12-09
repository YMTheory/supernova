#ifndef _SNchannelIBD_hh
#define _SNchannelIBD_hh

#include "SNchannels.hh"
#include "SNibdLS.hh"

class SNchannelIBD : public SNchannels{
    public:
        SNchannelIBD();
        virtual ~SNchannelIBD();

        virtual SNeffectLS* createChannel();

    //public:
    //    ClassDef(SNchannelIBD, 1);
};




#endif
