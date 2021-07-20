#!/bin/bash
export JobName=NewModel_RefOff_AbsOff
# export JobName=default
export TopPath=/junofs/users/wangyg/Software/juno-dev/J20v2r0-Pre1
export JobPath=/junofs/users/wangyg/PMT/AngleResponse/junoImpact/fullMC/juno
export OutPath=${JobPath}/output/pmt_v2/${JobName}
export SubPath=${OutPath}/sub
export LogPath=${OutPath}/log

function create_job()
{
SH_NAME=${SubPath}/run-$1-50.sh

if [ -f $SH_NAME ]; then
    rm $SH_NAME
fi

cat>$SH_NAME <<EOF
#!/bin/bash
source $TopPath/bashrc
cd $JobPath
(time python tut_detsim.py --evtmax 50 --enable-pmt-optical-model --seed $1 --no-gdml --output $OutPath/evt-$1-50.root --user-output $OutPath/user-evt-$1-50.root gun --particles e- --momentums 1.0 --momentums-interp KineticEnergy)
EOF

chmod +x $SH_NAME
hep_sub $SH_NAME -o $LogPath/log-${1}.out -e $LogPath/log-${1}.err -g juno
}

for i in {10000..10399}
do
    create_job $i
done
