#!/bin/bash
source SNsl7.sh 

path=$SNsim/simulation/examples/etSpec
runfile=${path}/etSpec2dVis

${runfile} -m $1 -c $2 
