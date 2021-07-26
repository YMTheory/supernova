#!/bin/bash
export SNSIMDIR=/junofs/users/miaoyu/supernova/wenlj
export OutPath=${SNSIMDIR}/dataset
export SubPath=${SNSIMDIR}/job
export LogPath=${SNSIMDIR}/log
dist=10

function create_job()
{
SH_NAME=${SubPath}/run-${1}_mh${2}_${dist}kpc_group${3}.sh

if [ -f $SH_NAME ]; then
    rm $SH_NAME
fi

#k=`expr $i/2`
k=`echo ${1} | awk '{printf("%g",$1/10.)}'`
echo $k

cat>$SH_NAME <<EOF
#!/bin/bash
cd $SNSIMDIR
#(time python generatePDFs.py 82503 2 1 $k ${dist}.0)
(time python generateDataSet.py 82503 1 0 $k ${2} 10.0 0.1 ${3})
EOF

chmod +x $SH_NAME
hep_sub $SH_NAME -o $LogPath/log-${1}_mh${2}_${dist}kpc_group${3}.out -e $LogPath/log-${1}_mh${2}_${dist}kpc_group${3}.err -g juno
}

for i in {0..20}
do
    for j in {0..2}
    do
        create_job $i $j 1

    done
done
