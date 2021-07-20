#ifndef _SNchannelNuP_hh
#define _SNchannelNuP_hh

#include "SNchannels.hh"
#include "SNnupLS.hh"
class SNchannelNuP : public SNchannels{
    public:
        SNchannelNuP();
        virtual ~SNchannelNuP();

        virtual SNeffectLS* createChannel();

    public:
        ClassDef(SNchannelNuP, 1);
};




#endif
