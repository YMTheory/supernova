#ifndef _SNchannelNCC_hh
#define _SNchannelNCC_hh

#include "SNchannels.hh"
#include "SNnccLS.hh"

class SNchannelNCC : public SNchannels{
    public:
        SNchannelNCC();
        virtual ~SNchannelNCC();

        virtual SNeffectLS* createChannel();

    public:
        ClassDef(SNchannelNCC, 1);
};




#endif
