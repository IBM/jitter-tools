#!/bin/bash
nodes=2
ppn=6
let cores_per_rank=42/$ppn
let nmpi=$nodes*$ppn
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o $(pwd)/output.out
#BSUB -e $(pwd)/output.err
#BSUB -nnodes ${nodes}
##BSUB -P VEN101
#BSUB -alloc_flags smt4
#BSUB -q excl_6gpus
#BSUB -W 15
#---------------------------------------
ulimit -s 10240
export OMP_NUM_THREADS=1
export SAVE_LIST=0
#export NVML_LOG=yes
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5
/opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun -D CUDA_VISIBLE_DEVICES \
       --tasks_per_rs ${ppn} --rs_per_host 1 --nrs ${nodes} --cpu_per_rs 42 --gpu_per_rs 6  \
      -d plane:${ppn} --bind=proportional-packed:${cores_per_rank}  $(pwd)/jitter-bench-cu -m lutexp_cuda -p 1000 -c 5000 
# -c 5 -t 300 -x 100 -n 50
EOF
#---------------------------------------
bsub < batch.job
