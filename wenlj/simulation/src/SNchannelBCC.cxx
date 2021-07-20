#include "SNchannelBCC.hh"


SNchannelBCC::SNchannelBCC():SNchannels(){

}

SNchannelBCC::~SNchannelBCC(){

}

SNeffectLS* SNchannelBCC::createChannel(){
    return new SNbccLS();
}

