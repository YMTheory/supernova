#ifndef _SNchannelBCC_hh
#define _SNchannelBCC_hh

#include "SNchannels.hh"
#include "SNbccLS.hh"

class SNchannelBCC : public SNchannels{
    public:
        SNchannelBCC();
        virtual ~SNchannelBCC();

        virtual SNeffectLS* createChannel();

    public:
        ClassDef(SNchannelBCC, 1);
};




#endif
