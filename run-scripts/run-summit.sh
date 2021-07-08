#!/bin/bash

if [[ $# -lt 6 ]]; then
    echo "usage: ${0} dir nodes ppn smt iso app"
    exit 1
fi

dir=${1};   shift
nodes=${1}; shift
ppn=${1};   shift
smt=${1};   shift
iso=${1};   shift
app=${*}

mkdir -p ${dir}
let nmpi=${nodes}*${ppn}

bsub << EOF

#BSUB -nnodes ${nodes}
#BSUB -o ${dir}/job%J.out
#BSUB -q tested
#BSUB -P VEN201
#BSUB -U IBM_MLX
#BSUB -W 60
#BSUB -J SUMMIT(${nodes})
##BSUB -alloc_flags "smt${smt}"
#BSUB -alloc_flags "smt${smt} cpublink"

jsrun\
    --rs_per_host ${ppn}\
    --nrs ${nmpi}\
    ${app}

EOF

