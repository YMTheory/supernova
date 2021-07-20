#include "SNchannelCNC.hh"


SNchannelCNC::SNchannelCNC():SNchannels(){

}

SNchannelCNC::~SNchannelCNC(){

}

SNeffectLS* SNchannelCNC::createChannel(){
    return new SNcncLS();
}

