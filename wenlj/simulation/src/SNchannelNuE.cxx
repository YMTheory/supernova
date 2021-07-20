#include "SNchannelNuE.hh"


SNchannelNuE::SNchannelNuE():SNchannels(){

}

SNchannelNuE::~SNchannelNuE(){

}

SNeffectLS* SNchannelNuE::createChannel(){
    return new SNnueLS();
}

