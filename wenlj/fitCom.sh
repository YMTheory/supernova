#!/bin/bash
export SNSIMDIR=/junofs/users/miaoyu/supernova/wenlj
export OutPath=${SNSIMDIR}/dataset
export SubPath=${SNSIMDIR}/job
export LogPath=${SNSIMDIR}/log

dist=10

function create_job()
{
SH_NAME=${SubPath}/fit-data${1}_mh${2}-pdf${3}_mh${4}-group${5}_${dist}kpc.sh

if [ -f $SH_NAME ]; then
    rm $SH_NAME
fi

#k=`expr $i/2`
k=`echo ${1} | awk '{printf("%g",$1/10.)}'`
j=`echo ${3} | awk '{printf("%g",$1/10.)}'`
echo $k, $j

cat>$SH_NAME <<EOF
#!/bin/bash
cd $SNSIMDIR
(time python nllFit.py 82503 1 0 $k ${2} $j ${4} ${dist}.0 0.1 ${5})
EOF

chmod +x $SH_NAME
hep_sub $SH_NAME -o $LogPath/logFit-data${1}_mh${2}-pdf${3}_mh${4}-group${5}_${dist}kpc.out -e $LogPath/logFit-data${1}_mh${2}-pdf${3}_mh${4}-group${5}_${dist}kpc.err -g juno
}

for i in {0..0}
do
    for m in {0..20}
    do
        for n in {1..1}
        do
           create_job $i 1 $m 1 $n
           create_job $i 1 $m 2 $n
           create_job $i 2 $m 1 $n
           create_job $i 2 $m 2 $n
        done
    done
done

#for i in $(seq 1 4); do
#    for j in $(seq 1 9); do
#        echo $i $j
#    done
#done
