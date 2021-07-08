#!/bin/bash

if [[ $# -lt 4 ]]; then
    echo "usage: ${0} nodes ppn iso app"
    exit 1
fi

nodes=${1}; shift
ppn=${1};   shift
iso=${1};   shift
app=${*}

let nmpi=${nodes}*${ppn}

bsub << EOF

#BSUB -nnodes ${nodes}
#BSUB -o job%J.out
#BSUB -q pibm
#BSUB -G guests
#BSUB -W 60
#BSUB -core_isolation ${iso}
#BSUB -alloc_flags "cpublink"

jsrun\
    --rs_per_host ${ppn}\
    --nrs ${nmpi}\
    --bind packed:2\
    ${app}

EOF
