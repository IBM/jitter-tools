#!/bin/bash 
nodes=8
ppn=40
let nmpi=$nodes*$ppn
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -nnodes ${nodes}
#BSUB -q excl 
#BSUB -W 15
#---------------------------------
ulimit -s 10240
/opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun --rs_per_host ${ppn} --nrs ${nmpi} ./variability -t 10 -i  500 -x 10000 -n 50
EOF
#---------------------------------------
bsub -core_isolation y -m c699c031 <batch.job
