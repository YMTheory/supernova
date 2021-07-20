#ifndef _SNchannelCNC_hh
#define _SNchannelCNC_hh

#include "SNchannels.hh"
#include "SNcncLS.hh"

class SNchannelCNC : public SNchannels{
    public:
        SNchannelCNC();
        virtual ~SNchannelCNC();

        virtual SNeffectLS* createChannel();

    public:
        ClassDef(SNchannelCNC, 1);
};




#endif
