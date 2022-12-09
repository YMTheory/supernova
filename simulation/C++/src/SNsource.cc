#include "SNsource.hh"

SNsource::SNsource():
    totE(30e52),dist(10),
    fEvmin(0.),fEvmax(100){
   
    for(int ii=0; ii<6; ii++){    
        Ea[ii] = 0;
        averagE[ii] = 12;
    }
}

SNsource::~SNsource(){

}




