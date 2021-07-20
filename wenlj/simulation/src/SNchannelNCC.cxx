#include "SNchannelNCC.hh"


SNchannelNCC::SNchannelNCC():SNchannels(){

}

SNchannelNCC::~SNchannelNCC(){

}

SNeffectLS* SNchannelNCC::createChannel(){
    return new SNnccLS();
}

