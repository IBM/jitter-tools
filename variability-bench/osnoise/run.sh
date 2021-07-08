#!/bin/bash
nodes=4500
ppn=40
let nmpi=$nodes*$ppn
tstamp=`date +%m_%d_%H_%M_%S`
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -nnodes ${nodes}
#BSUB -P VEN201
##BSUB -alloc_flags "cpublink"
#BSUB -alloc_flags "smt4 cpublink"
#BSUB -q tested
#BSUB -W 40
#---------------------------------------
ulimit -s 10240
export OMP_NUM_THREADS=1
export SAVE_LIST=0
export RANKS_PER_NODE=${ppn}
export PAMI_IBV_ENABLE_DCT=1
export PAMI_PMIX_DATACACHE=1

#export HPM_EVENT_LIST=PM_RUN_CYC,PM_RUN_INST_CMPL,PM_L1_DCACHE_RELOADED_ALL,PM_FLUSH
export HPM_GROUP=1
export OMPI_LD_PRELOAD_POSTPEND=/ccs/home/walkup/mpitrace/spectrum_mpi/libmpihpm.so

export BIND_THREADS=yes
export SYSTEM_CORES=2

#  --smpiargs="-MXM"

  jsrun  -X 1 --progress ${tstamp}.progress    \
   --tasks_per_rs ${ppn} --rs_per_host 1 --nrs ${nodes} --cpu_per_rs 42 --gpu_per_rs 6 -d plane:${ppn} \
   --bind proportional-packed:1  ./no_aslr ./osnoise -c 5 -t 400 -x 100 -n 31

#jsrun --tasks_per_rs ${ppn} --rs_per_host 1 --nrs ${nodes} --cpu_per_rs 42 --gpu_per_rs 6 -d plane:${ppn}  ./h4.sh ./scalar
#jsrun --tasks_per_rs ${ppn} --rs_per_host 1 --nrs ${nodes} --cpu_per_rs 42 --gpu_per_rs 6 -d plane:${ppn}  ./h4.sh ./simple
#jsrun --tasks_per_rs ${ppn} --rs_per_host 1 --nrs ${nodes} --cpu_per_rs 42 --gpu_per_rs 6 -d plane:${ppn}  ./h4.sh ./root
EOF
#---------------------------------------
 bsub  batch.job
