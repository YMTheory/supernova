## Models
In simulation/data, you can find all the numerical models collected from different groups, covering SN, preSN and long cooling simualtions.

## To build SNsim in your local directory:

1. set the ROOT environment $ROOTSYS
    e.g use JUNO software environment or your own ROOT

2. replace the path in simulation/src/ to your path:
cd your-simulation-path
sed -i "s/\/junofs\/users\/lihl\/neutrino\/superN\/SNsim\/simulation/your-simulation-path/g" `grep /junofs/users/lihl/neutrino/superN/SNsim/simulation -rl src/`

3. make
    In JUNO software environment, please 'make' twice. Makefile should be update later.

## To load libSNsim.so in ROOT interactive mode, please add the following lines into your rootlogon.C
    gROOT->ProcessLine(Form(".include %s/include",path.data()));
    gSystem->Load(Form("%s/lib/libSNsim.so", path.data()));

## Find examples in simulation/examples/


