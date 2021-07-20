#!/bin/bash
export SNSIMDIR=/junofs/users/wenlj/juno/supernova/workStation/SNsim
export OutPath=${SNSIMDIR}/dataset
export SubPath=${SNSIMDIR}/job
export LogPath=${SNSIMDIR}/log

function create_job()
{
SH_NAME=${SubPath}/runMH2-${1}_2kpc.sh

if [ -f $SH_NAME ]; then
    rm $SH_NAME
fi

#k=`expr $i/2`
k=`echo ${1} | awk '{printf("%g",$1/10.)}'`
echo $k

cat>$SH_NAME <<EOF
#!/bin/bash
cd $SNSIMDIR
(time python generatePDFsMH2.py 82503 1 0 $k 2.0)
EOF

chmod +x $SH_NAME
hep_sub $SH_NAME -o $LogPath/logMH2-${1}_2kpc.out -e $LogPath/logMH2-${1}_2kpc.err -g juno
}

for i in {0..20}
do
    create_job $i
done
