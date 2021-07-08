#!/bin/bash

#
#   IBM Corporation (C) 2018
#

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

jsrun=/opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun
mpirun=${MPI_ROOT}/bin/mpirun
smt=1

dir=""
nod=""
ppn=""
app=""

iso="2"
trc="off"

function usage {
    echo "usage: $0 REQUIRED [OPTIONAL] -- APP..."
    echo "required:"
    echo "  -d, --dir DIR       directory to write outputs"
    echo "  -n, --nod NOD       number of compute nodes"
    echo "  -p, --ppn PPN       processes per node"
    echo "  -x, --    APP       application to run on compute nodes"
    echo "optional:"
    echo "  -i, --iso ISO       core isolation (0-6) (default: ${iso})"
    echo "  -t, --trc TRC       enable tracing (trace[f][p])"
    echo "  "
    echo
    exit 1
}

while ((${#})); do
    if   [[ ${1} == "-h" || ${1} == "--help" ]]; then usage
    elif [[ ${1} == "-d" || ${1} == "--dir"  ]]; then dir=${2}; shift; shift
    elif [[ ${1} == "-n" || ${1} == "--nod"  ]]; then nod=${2}; shift; shift
    elif [[ ${1} == "-p" || ${1} == "--ppn"  ]]; then ppn=${2}; shift; shift
    elif [[ ${1} == "-i" || ${1} == "--iso"  ]]; then iso=${2}; shift; shift
    elif [[ ${1} == "-t" || ${1} == "--trc"  ]]; then trc=${2}; shift; shift
    elif [[ ${1} == "-x" || ${1} == "--"     ]]; then shift; app=${*}; break
    else
        echo "[run] error: invalid option: '${1}'"
        usage
    fi
done

if [[ -z ${dir} ||
      -z ${nod} ||
      -z ${ppn} ||
      -z ${app} ]]; then
    echo "[run] error: missing required parameters"
    usage
fi

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

mkdir -p ${dir}
dir=$(readlink -f ${dir})

lsf=""
lsf+="-core_isolation ${iso} "

run=""
run+="${jsrun} "
run+="--rs_per_host ${ppn} "
run+="--nrs $((${nod}*${ppn})) "
run+="${app}"

# ------------------------------------------------------------------------------
# Launch
# ------------------------------------------------------------------------------

bsub << EOF

#BSUB -nnodes ${nod}
#BSUB ${lsf}
#BSUB -q pibm
#BSUB -G guests
#BSUB -o ${dir}/job%J.out
#BSUB -J RUN(${nod})
#BSUB -alloc_flags "tracefp"

function echoval {
    cmd=\${*}
    echo -n "[run] $ "
    echo \${cmd}
    eval \${cmd}
}

echo
echo ============================================================
echo SETUP
echo ============================================================
echo

echo
echo ----------------------------------
echo Job
echo ----------------------------------
echo

if [[ -z \${LSB_JOBID} ]]; then
    echo "[run] error: variable LSB_JOBID is not set"
    exit 1
fi

job=\${LSB_JOBID}
echo \${job}

echo
echo ----------------------------------
echo Hosts
echo ----------------------------------
echo

if [[ -z \${LSB_DJOB_HOSTFILE} ]]; then
    echo "[run] error: variable LSB_DJOB_HOSTFILE is not set"
    exit 1
fi

if [[ ! -f \${LSB_DJOB_HOSTFILE} ]]; then
    echo "[run] error: could not find hostfile: '\${LSB_DJOB_HOSTFILE}'"
    exit 1
fi

hosts=\$(cat \${LSB_DJOB_HOSTFILE} | tail --lines=+2 | uniq)
hostfile=${dir}/job\${job}.hostfile
for host in \${hosts}; do
    echo "\${host} slots=${ppn}" >> \${hostfile}
    echo "\${host}"
done

echo
echo ----------------------------------
echo Tracing
echo ----------------------------------
echo

if [[ ${trc} != "off" ]]; then
    echo "start tracing..."
    export JT_JOB=\${job}
    export JT_DIR=${dir}
    export JT_ON=1

    for host in \${hosts}; do
        cmd="${SDIR}/tracefs/tracefs-start.sh ${trc}"
        echoval "ssh \${host} \"\${cmd}\" &"
    done
    wait
    echo "tracing started"
fi

echo
echo ============================================================
echo RUN
echo ============================================================
echo

echoval ${run}

echo
echo ============================================================
echo CLEANUP
echo ============================================================
echo

echo 
echo ----------------------------------
echo Tracing
echo ----------------------------------
echo

if [[ ${trc} != "off" ]]; then
    echo "stop tracing..."
    for host in \${hosts}; do
        cmd="${SDIR}/tracefs/tracefs-stop.sh"
        echoval "ssh \${host} \"\${cmd}\" &"
    done
    wait
    echo "tracing stopped"
    echo

    #tracefs_root="/sys/kernel/debug/tracing/instances"
    #for host in \${hosts}; do
    #    for instance in ftrace tracep; do
    #        #echoval "ssh \${host} cp \${tracefs_root}/\${instance}/trace /tmp"
    #        echoval "ssh \${host} cp \${tracefs_root}/\${instance}/per_cpu/cpu0/trace /tmp"
    #        echoval "scp \${host}:/tmp/trace $(hostname):${dir}/job\${job}-trace-\${host}.\${instance}"
    #        echoval "ssh \${host} rm /tmp/trace"
    #        echo
    #    done
    #done
    #echo "traces collected"
    #echo

    echo "collect traces..."
    echoval "${SDIR}/tracefs/tracefs-collect.sh \${job} ${dir}"
    echo "traces collected"
    echo

    echo "reset tracing..."
    for host in \${hosts}; do
        cmd="${SDIR}/tracefs/tracefs-reset.sh"
        echoval "ssh \${host} \"\${cmd}\" &"
    done
    wait
    echo "tracing reset" 
    echo
fi

EOF

