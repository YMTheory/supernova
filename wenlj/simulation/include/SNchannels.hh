#ifndef _SNchannels_hh
#define _SNchannels_hh

#include "TObject.h"
#include "SNeffectLS.hh"
class SNchannels : public TObject{
    public:
        SNchannels();
        virtual ~SNchannels();

        virtual SNeffectLS* createChannel() = 0;

    public:
        ClassDef(SNchannels, 1);
};



#endif
