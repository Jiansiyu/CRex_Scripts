#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
    echo "./gen.sh \"runIDs\" nfoil"
    echo "example: ./gen \"20825 20826 20827\" 1"
    exit
fi
LIST=$1
nfoil=$2

((nfoil=nfoil -1))
# Number of Event in each Cut area
N=100000
nCol=13
echo $LIST $nfoil
for run in $LIST 
  do
    echo " >working on run $run"
    # generate the files used for store the cut data
    echo $run
    : >prexRHRS_$run.SieveCut.cut
    : >prexRHRS_$run.SieveCut.ncut
    for i in $(seq 0 ${nfoil})
       do 
            #loop on the col ID, because all the raw in a col will store in a same file
            for j in  $(seq 0 ${nCol})
            do
               name=$(echo -e "prexRHRS_${run}_-1.root.SieveCut.$i"_"$j.cut") #col
               echo "   >>>> Looking for file $name"
               if [[ -e "$name" ]]; then
                 while read l
                 do 
                    echo -e "$l" >>prexRHRS_${run}.SieveCut.cut
                    if [[ "$l" == "fEvtHdr.fRun==0" ]];then
                        echo -e "0" >>prexRHRS_${run}.SieveCut.ncut 
                    else
                        echo -e "$N" >>prexRHRS_${run}.SieveCut.ncut
                    fi
                 done < $name
                else
               for k in $(seq 0 6) #holes
                do
                    echo -e "fEvtHdr.fRun==0" >>prexRHRS_$run.SieveCut.cut
                    echo -e "0" >>prexRHRS_$run.SieveCut.ncut
                done
               fi
               echo -e  "" >>prexRHRS_$run.SieveCut.cut 
               echo -e  "" >>prexRHRS_$run.SieveCut.ncut
            done
            echo -e  "" >>prexRHRS_$run.SieveCut.cut
            echo -e  "" >>prexRHRS_$run.SieveCut.ncut
       done
done
