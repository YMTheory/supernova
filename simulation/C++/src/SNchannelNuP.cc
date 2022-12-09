#include "SNchannelNuP.hh"


SNchannelNuP::SNchannelNuP():SNchannels(){

}

SNchannelNuP::~SNchannelNuP(){

}

SNeffectLS* SNchannelNuP::createChannel(){
    return new SNnupLS();
}

