#!/bin/bash
nodes=1
ppn=4
let nmpi=$nodes*$ppn
let cores_per_rank=40/$ppn
let cores=40*${nodes}
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
##BSUB -nnodes ${nodes}
#BSUB -csm y 
#BSUB -R "1*{select[LN]} + ${cores}*{select[CN&&(hname==sierra2837)&&(type==any)]span[ptile=42]}"
#BSUB -G guests
#BSUB -alloc_flags "smt1 cpublink"
#BSUB -core_isolation 2
#BSUB -J osnoise
#BSUB -q pdebug
#BSUB -W 30
#---------------------------------------
jsrun  --tasks_per_rs ${ppn} --rs_per_host 1 --nrs ${nodes} --cpu_per_rs 40 --gpu_per_rs 4  \
      -d plane:${ppn} --bind=proportional-packed:${cores_per_rank}  ./osnoise -d -c 5 -t 40 -x 100 -n 31

EOF
#---------------------------------------
 bsub  <batch.job
