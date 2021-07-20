#!/bin/bash

###------Cooling models---------//
#find ./Cooling -name "*nu*" | xargs sed -i '/^[@#]/d'
#path=./Cooling/LS220
#files=$(find ${path} -name \*.\*)
#for iname in ${files}
#do
#    inameCurr=$(echo ${iname} | awk -F"/" '{print $4}')
#    echo ${inameCurr}
#    inamePath=$(echo ${iname} | awk -F"/" '{print $1"/"$2"/"$3"/"}')
#    echo ${inamePath}
#    inameNew=$(echo ${inameCurr} | sed -e 's/$/0/g' -e 's/^s/8/g' -e 's/\.//g' -e 's/co0/3/g' -e 's/c0/1/g' -e 's/o0/2/g' -e 's/r0/4/g')
#    echo ${inameNew}
#    mv ${iname} ./${inameNew}
#done
#
#
#path=./Cooling/Shen
#files=$(find ${path} -name \*.\*)
#for iname in ${files}
#do
#    inameCurr=$(echo ${iname} | awk -F"/" '{print $4}')
#    echo ${inameCurr}
#    inamePath=$(echo ${iname} | awk -F"/" '{print $1"/"$2"/"$3"/"}')
#    echo ${inamePath}
#    inameNew=$(echo ${inameCurr} | sed -e 's/$/0/g' -e 's/^s/9/g' -e 's/\.//g' -e 's/co0/3/g' -e 's/c0/1/g' -e 's/o0/2/g' -e 's/r0/4/g')
#    echo ${inameNew}
#    mv ${iname} ./${inameNew}
#done



##-----------Accretion models-----------//
find ./Accretion -name "*nu*" | xargs sed -i '/^[@#]/d'
path=./Accretion/Shen
files=$(find ${path} -name \*.\*)
for iname in ${files}
do
    inameCurr=$(echo ${iname} | awk -F"/" '{print $4}')
    echo ${inameCurr}
    inamePath=$(echo ${iname} | awk -F"/" '{print $1"/"$2"/"$3"/"}')
    echo ${inamePath}
    inameNew=$(echo ${inameCurr} | sed -e 's/$/0/g' -e 's/^s/9/g' -e 's/\.//g' -e 's/co0/3/g' -e 's/c0/1/g' -e 's/o0/2/g' -e 's/r0/4/g')
    echo ${inameNew}
    if [ ! -d ${inameNew} ]; then
        echo "Oh, new !!"
        mv ${iname} ./${inameNew}
    fi
done

path=./Accretion/LS220
files=$(find ${path} -name \*.\*)
for iname in ${files}
do
    inameCurr=$(echo ${iname} | awk -F"/" '{print $4}')
    echo ${inameCurr}
    inamePath=$(echo ${iname} | awk -F"/" '{print $1"/"$2"/"$3"/"}')
    echo ${inamePath}
    inameNew=$(echo ${inameCurr} | sed -e 's/$/0/g' -e 's/^s/8/g' -e 's/\.//g' -e 's/co0/3/g' -e 's/c0/1/g' -e 's/o0/2/g' -e 's/r0/4/g')
    echo ${inameNew}
    if [ ! -d ${inameNew} ]; then
        echo "Oh, new !!"
        mv ${iname} ./${inameNew}
    fi
done

path=./Accretion/LS180
files=$(find ${path} -name \*.\*)
for iname in ${files}
do
    inameCurr=$(echo ${iname} | awk -F"/" '{print $4}')
    echo ${inameCurr}
    inamePath=$(echo ${iname} | awk -F"/" '{print $1"/"$2"/"$3"/"}')
    echo ${inamePath}
    inameNew=$(echo ${inameCurr} | sed -e 's/$/0/g' -e 's/^s/7/g' -e 's/\.//g' -e 's/co0/3/g' -e 's/c0/1/g' -e 's/o0/2/g' -e 's/r0/4/g')
    echo ${inameNew}
    if [ ! -d ${inameNew} ]; then
        echo "Oh, new !!"
        mv ${iname} ./${inameNew}
    fi
done

