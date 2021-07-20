#include "SNchannelIBD.hh"


SNchannelIBD::SNchannelIBD():SNchannels(){

}

SNchannelIBD::~SNchannelIBD(){

}

SNeffectLS* SNchannelIBD::createChannel(){
    return new SNibdLS();
}

