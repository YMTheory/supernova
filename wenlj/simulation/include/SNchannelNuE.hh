#ifndef _SNchannelNuE_hh
#define _SNchannelNuE_hh

#include "SNchannels.hh"
#include "SNnueLS.hh"
class SNchannelNuE : public SNchannels{
    public:
        SNchannelNuE();
        virtual ~SNchannelNuE();

        virtual SNeffectLS* createChannel();

    public:
        ClassDef(SNchannelNuE, 1);
};




#endif
