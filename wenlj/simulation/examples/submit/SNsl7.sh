#!/bin/bash

###JUNO###
source /cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J20v2r0-Pre0/setup.sh

###SNsim#####
export LD_LIBRARY_PATH=$SNsim/simulation/lib:$LD_LIBRARY_PATH

