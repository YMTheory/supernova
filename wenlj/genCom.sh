#!/bin/bash
export SNSIMDIR=/junofs/users/wenlj/juno/supernova/workStation/SNsim
export OutPath=${SNSIMDIR}/dataset
export SubPath=${SNSIMDIR}/job
export LogPath=${SNSIMDIR}/log
dist=3

function create_job()
{
SH_NAME=${SubPath}/run-${1}_${dist}kpc_group${2}.sh

if [ -f $SH_NAME ]; then
    rm $SH_NAME
fi

#k=`expr $i/2`
k=`echo ${1} | awk '{printf("%g",$1/10.)}'`
echo $k

cat>$SH_NAME <<EOF
#!/bin/bash
cd $SNSIMDIR
(time python generatePDFs.py 82503 2 1 $k ${dist}.0)
EOF
#(time python generateDataSet.py 82503 1 0 $k 10.0 0.1 ${2})

chmod +x $SH_NAME
hep_sub $SH_NAME -o $LogPath/log-${1}_${dist}kpc_group${2}.out -e $LogPath/log-${1}_${dist}kpc_group${2}.err -g juno
}

for i in {0..20}
do
    create_job $i 1
done
