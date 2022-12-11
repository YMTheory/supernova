#include "SNchannelCEvNS.hh"

SNchannelCEvNS::SNchannelCEvNS():SNchannels() {}

SNchannelCEvNS::~SNchannelCEvNS() {}

SNeffectLS* SNchannelCEvNS::createChannel(){
    return new SNCEvNSLS();
}
