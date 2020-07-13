#!/bin/bash


if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    echo "./gen.sh \"startRun-EndRun\" "
    echo "example: ./gen \"20825-20827\" "
    exit
fi

#!/bin/bash
string=$1
read num1 num2 <<<${string//[^0-9]/ }
runArray=()
bashdir="/home/newdriver/pyQuant/prex_replayed/rootfile/"
while [ $num2 -gt $num1 ]
do
    if test -f "${bashdir}/prexLHRS_${num1}_-1.root"; then
        echo "${bashdir}/prexLHRS_${num1}_-1.root"
        runArray+=($num1)
    fi
    
    if test -f "${bashdir}/prexRHRS_${num1}_-1.root"; then
        runArray+=($num1)
    fi
    num1=$(($num1+1))
done

for run in ${runArray[@]}
do
    echo " >working on run $run"
    analyzer  -b -q .L rootlogon.C 'TargetCheck.C('${run}')'
done