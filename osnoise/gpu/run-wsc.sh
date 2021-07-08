#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

TRACEFS=${SDIR}/../../run-scripts/tracefs
TRACEDIR=/tmp

nod=""
ppn=""
iso=""
app=""
dir="."
trc="off"

function usage {
  echo "usage: ${0} REQUIRED [OPTIONAL] -- ./osnoise PARAMS"
  echo "required:"
  echo "  --nod NOD   number of nodes"
  echo "  --ppn PPN   processes per node (GPUs per node)"
  echo "  --iso ISO   number of jitter isolation cores"
  echo "  --    APP   osnoise binary and parameters"
  echo "optional:"
  echo "  --dir DIR   base output directory (default: ${dir})"
  echo "  --trc TRC   tracing option (default: ${trc})"
  exit 1
}

while ((${#})); do
  if   [[ ${1} == "-h" || ${1} == "--help" ]]; then usage
  elif [[ ${1} == "-n" || ${1} == "--nod"  ]]; then shift; nod=${1}; shift
  elif [[ ${1} == "-p" || ${1} == "--ppn"  ]]; then shift; ppn=${1}; shift
  elif [[ ${1} == "-i" || ${1} == "--iso"  ]]; then shift; iso=${1}; shift
  elif [[ ${1} == "-x" || ${1} == "--"     ]]; then shift; app=${*}; break
  elif [[ ${1} == "-d" || ${1} == "--dir"  ]]; then shift; dir=${1}; shift
  elif [[ ${1} == "-t" || ${1} == "--trc"  ]]; then shift; trc=${1}; shift
  else
    echo "error: invalid parameter: ${1}"
    usage
  fi
done

function check {
  param=${1}
  if [[ -z ${!param} ]]; then
    echo "error: missing required parameter: ${param}"
    usage
  fi
}

check "nod"
check "ppn"
check "iso"
check "app"

nmpi=$(( nod * ppn ))
cores=$(( 44 - 2 * iso ))
cores_per_rank=$(( cores / ppn ))
cores_total=$(( nod * cores ))

alloc_flags="smt1 cpublink "
if [[ ${trc} != "off" ]]; then
  alloc_flags+="${trc} "
fi

mkdir -p ${dir}
dir=$(readlink -f ${dir})

bsub << EOF

#BSUB -o ${dir}/job%J.out
#BSUB -nnodes ${nod}
#BSUB -core_isolation ${iso}
#BSUB -alloc_flags "${alloc_flags}"
#BSUB -J osn(${nod})

function echoval {
  local cmd=\${*}
  echo "\${cmd}" | tr -s " "
  eval "\${cmd}"
}

echo
echo ---------------------------------------------------------------------------
echo SETUP
echo ---------------------------------------------------------------------------
echo

run="/opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun \
  --tasks_per_rs ${ppn} \
  --rs_per_host 1 \
  --nrs ${nod} \
  --cpu_per_rs ${cores} \
  --gpu_per_rs ${ppn} \
  --launch_distribution plane:${ppn} \
  --bind proportional-packed:${cores_per_rank} \
  -D CUDA_VISIBLE_DEVICES \
  ${app}"

job=\${LSB_JOBID}
hosts=\$(tail -n +2 \${LSB_DJOB_HOSTFILE} | sort -u)
hostfile="${dir}/job\${job}.hosts"
echo \${hosts} | sed "s/ /\n/g" > \${hostfile}

echo nod:   ${nod}
echo ppn:   ${ppn}
echo iso:   ${iso}
echo app:   ${app}
echo trc:   ${trc}
echo dir:   ${dir}
echo run:   \${run} | tr -s " "
echo job:   \${job}
echo hosts: \${hosts}

echo
echo ---------------------------------------------------------------------------
echo TRACE START
echo ---------------------------------------------------------------------------
echo

if [[ ${trc} == "off" ]]; then
  echo "tracing: off -- nothing to do"
else
  echo "tracing: start"
  for host in \${hosts}; do
    cmd="${TRACEFS}/tracefs-start.sh ${trc}"
    echoval "ssh \${host} \"\${cmd}\" &"
  done
  wait
  echo "tracing: started"
fi

echo
echo ---------------------------------------------------------------------------
echo RUN
echo ---------------------------------------------------------------------------
echo

\${run}

echo
echo ---------------------------------------------------------------------------
echo TRACE STOP
echo ---------------------------------------------------------------------------
echo

if [[ ${trc} == "off" ]]; then
  echo "tracing: off -- nothing to do"
else
  echo "tracing: stop"
  for host in \${hosts}; do
    cmd="${TRACEFS}/tracefs-stop.sh"
    echoval "ssh \${host} \"\${cmd}\" &"
  done
  wait
  echo "tracing: stopped"
fi

echo
echo ---------------------------------------------------------------------------
echo TRACE COLLECT
echo ---------------------------------------------------------------------------
echo

if [[ ${trc} == "off" ]]; then
  echo "tracing: off -- nothing to do"
else
  echo "tracing: collect"
  pdsh -w ^\${hostfile} \
    "${TRACEFS}/tracefs-collect.sh \${job} ${dir} ${TRACEDIR}"
  echo "tracing: collect done"
  echo
  echo "tracing: move"
  pdsh -w ^\${hostfile} \
    "scp ${TRACEDIR}/job\${job}* $(hostname):${dir}/; rm ${TRACEDIR}/job\${job}*"
  echo "tracing: move done"
  echo
  echo "tracing: clean"
  for host in \${hosts}; do
    cmd="${TRACEFS}/tracefs-reset.sh"
    echoval "sudo ssh \${host} \"\${cmd}\" &"
  done
  wait
  echo "tracing: clean done"
fi

EOF

