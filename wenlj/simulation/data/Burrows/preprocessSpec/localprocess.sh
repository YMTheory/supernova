#!/bin/bash

path=/junofs/users/lihl/neutrino/superN/SNsim/simulation/data/Burrows/preprocessSpec
for imod in {6120,6150,6200,6250,7120,7150,7200,7250}
do
    ${path}/makeModelspectrum -m ${imod}
done
