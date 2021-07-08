#!/bin/bash
nodes=8
ppn=40
let nmpi=$nodes*$ppn
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=${ppn}]"
#BSUB -R "select[type=any]"
#BSUB -n ${nmpi}
#BSUB -x
#BSUB -q excl
#BSUB -W 30
#---------------------------------------
ulimit -s 10240
/opt/ibm/spectrum_mpi/bin/mpirun --tag-output --bind-to none -np $nmpi ./variability -t 10 -i  5000 -x 10000 -n 50
EOF
#---------------------------------------
bsub -core_isolation y  <batch.job
#bsub  < batch.job
