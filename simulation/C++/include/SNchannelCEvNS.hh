#ifndef SNchannelCEvNS_hh
#define SNchannelCEvNS_hh

#include "SNchannels.hh"
#include "SNCEvNSLS.hh"

class SNchannelCEvNS : public SNchannels {
    public:
        SNchannelCEvNS();
        virtual ~SNchannelCEvNS();

        virtual SNeffectLS* createChannel();

};

#endif

