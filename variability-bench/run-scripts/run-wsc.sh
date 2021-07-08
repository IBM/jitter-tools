#!/bin/bash

# ------------------------------------------------
# Script's path
# ------------------------------------------------

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# ------------------------------------------------
# Definitions
# ------------------------------------------------

jsrun="/opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun"
mpirun="/opt/ibm/spectrum_mpi/bin/mpirun"
cpuset="${SDIR}/cpuset-wsc.sh"
user="${USER}"
tracedir="/tmp/traces"

# ------------------------------------------------
# Arguments
# ------------------------------------------------

# Required paramters
dir=""
nod=""
app=""

# Optional parameters (default values)
ppn="auto"
smt="1"
iso="1"
gpu="4"
lau="jsrun"
mpi="pami"
mpc="libcoll"
omp="1"
que="excl"
mac="off"
trc="off"
trf="auto"
jid="RUN"
djt="off"
win="0:1000000"
gpf="off"

# Print usage information and exit
function usage {
    echo "usage: ${0} -d DIR -n NOD [OPTIONAL] -x APP"
    echo "required:"
    echo "  -d DIR  output directory"
    echo "  -n NOD  number of nodes"
    echo "  -x APP  application executable and arguments"
    echo "optional:"
    echo "  -p PPN  processes per node (default: ${ppn})"
    echo "  -s SMT  SMT configuration (default: ${smt})"
    echo "  -i ISO  isolation option (off|1|2|...) (default: ${iso})"
    echo "  -g GPU  GPUs per node (default: ${gpu})"
    echo "  -l LAU  launcher (jsrun|mpirun) (default: ${lau})"
    echo "  -k MPI  MPI interconnect (pami|mxm) (default: ${mpi})"
    echo "  -c MPC  MPI collectives (libcoll|hcoll) (default: ${mpc})"
    echo "  -M OMP  OpenMP threads (default: ${omp})"
    echo "  -q QUE  LSF queue (default: ${que})"
    echo "  -m MAC  select machines to run (default: ${mac})"
    echo "            off (disabled)"
    echo "            comma-separated list of hosts"
    echo "  -t TRC  tracing (default: ${trc})"
    echo "            off             disabled"
    echo "            perf            perf (mmap => file) [REQUIRES SUDO]"
    echo "            trace[f][p]  tracefs (circular buffer)"
    echo "                  │  └── tracepoints (Perf events)"
    echo "                  └───── ftracing (kernel functions)"
    echo "  -f TRF  tracing CPU filter (default: ${trf})"
    echo "            auto       automatically select CPUs based on SMT and iso"
    echo "            off        disable filter (trace all CPUs)"
    echo "            CPU...     comma-separated list of CPUs"
    echo "  -j JID  job identification (default: ${jid})"
    echo "  -D DJT  dejitter allocation cgroup in prolog (default: ${djt})"
    echo "  -W WIN  jwind <window>:<period> in prolog (default: ${win})"
    echo "  -G GPF  use GPFS option (default: ${gpf})"
    echo
    exit 1
}

# Parse command line arguments
while ((${#})); do
    if   [[ ${1} == "-h" || ${1} == "--help" ]]; then usage
    elif [[ ${1} == "-d" || ${1} == "--dir"  ]]; then dir=${2}; shift; shift
    elif [[ ${1} == "-n" || ${1} == "--nod"  ]]; then nod=${2}; shift; shift
    elif [[ ${1} == "-p" || ${1} == "--ppn"  ]]; then ppn=${2}; shift; shift
    elif [[ ${1} == "-s" || ${1} == "--smt"  ]]; then smt=${2}; shift; shift
    elif [[ ${1} == "-i" || ${1} == "--iso"  ]]; then iso=${2}; shift; shift
    elif [[ ${1} == "-l" || ${1} == "--lau"  ]]; then lau=${2}; shift; shift
    elif [[ ${1} == "-k" || ${1} == "--mpi"  ]]; then mpi=${2}; shift; shift
    elif [[ ${1} == "-c" || ${1} == "--mpc"  ]]; then mpc=${2}; shift; shift
    elif [[ ${1} == "-M" || ${1} == "--omp"  ]]; then omp=${2}; shift; shift
    elif [[ ${1} == "-q" || ${1} == "--que"  ]]; then que=${2}; shift; shift
    elif [[ ${1} == "-m" || ${1} == "--mac"  ]]; then mac=${2}; shift; shift
    elif [[ ${1} == "-t" || ${1} == "--trc"  ]]; then trc=${2}; shift; shift
    elif [[ ${1} == "-f" || ${1} == "--trf"  ]]; then trf=${2}; shift; shift
    elif [[ ${1} == "-j" || ${1} == "--jid"  ]]; then jid=${2}; shift; shift
    elif [[ ${1} == "-D" || ${1} == "--djt"  ]]; then djt=${2}; shift; shift
    elif [[ ${1} == "-W" || ${1} == "--win"  ]]; then win=${2}; shift; shift
    elif [[ ${1} == "-G" || ${1} == "--gpf"  ]]; then gpf=${2}; shift; shift
    elif [[ ${1} == "-x" || ${1} == "--app"  ]]; then shift; app=${*}; break
    else
        echo "error: invalid option: '$1'"
        usage
    fi
done

# Check required flags
if [[ -z ${nod} || -z ${dir} || -z ${app} ]]; then
    echo "error: missing required arguments"
    usage
fi

# ------------------------------------------------
# Setup
# ------------------------------------------------

# Prepare output directory
mkdir -p ${dir}
dir=$(readlink -f ${dir})

# Prepare LSF options
lsf=""
aflags=""

# Handle isolation option; increment LSF options if necessary
if   [[ ${iso} == "off" || ${iso} == 0 ]]; then
    if [[ ${ppn} == "auto" ]]; then
        ppn=$((44*${smt}))
    fi
else
    lsf+=" -core_isolation ${iso}"
    if [[ ${ppn} == "auto" ]]; then
        ppn=$(((44-${iso}*2)*${smt}))
    fi
fi

# Handle SMT option
aflags+="smt${smt} "

# Handle MPI interconnect option
mpiargs=""
if   [[ ${mpi} == "pami" ]]; then
    true
elif [[ ${mpi} == "mxm" ]]; then
    mpiargs+="-MXM "
else
    echo "error: invalid MPI interconnect option: '${mpi}'"
    usage
fi

# Handle MPI collectives option
if   [[ ${mpc} == "libcoll" ]]; then
    true
elif [[ ${mpc} == "hcoll" ]]; then
    mpiargs+="-HCOLL "
else
    echo "error: invalid MPI collectives option: '${mpc}'"
    usage
fi

# Handle launcher option and prepare run command
if   [[ ${lau} == "jsrun" ]]; then
    run="${jsrun} "
    run+="--nrs ${nod} "
    run+="--tasks_per_rs ${ppn} "
    run+="--gpu_per_rs ${gpu} "
    run+="--rs_per_host 1 "
    run+="--launch_distribution packed "
    run+="--smpiargs=\\\"${mpiargs}\\\" "
elif [[ ${lau} == "mpirun" ]]; then
    run="${mpirun} "
    run+="-x USING_MPIRUN -x SMT -x RANKS_PER_NODE -x ISOLATION -x OMP_NUM_THREADS "
    run+="-x JTF_JOB -x JTF_DIR "
    run+="--hostfile \\\${hostfile} "
    run+="${mpiargs} "
else
    echo "error: invalid launcher: '${lau}'"
    usage
fi
run+="${cpuset} ${app}"

# Handle machine selection flag
if [[ ${mac} == "off" ]]; then
    mac="##BSUB -m"
else
    mac="#BSUB -m \"c699launch01 $(echo ${mac} | tr ',' ' ')\""
fi

# Check tracing option
if [[ ${trc} != "off"      &&\
      ${trc} != "perf"     &&\
      ${trc} != "tracef"   &&\
      ${trc} != "tracep"   &&\
      ${trc} != "tracefp"  ]]; then
    echo "error: invalid tracing option: '${trc}'"
    usage
elif [[ ${trc} != "off" && ${trc} != "perf" ]]; then
    aflags+="${trc} "
fi

# Handle automatic tracing filter options
if [[ ${trf} == "auto" || ${trf} == "off" ]]; then
    sockets=2
    cores=22
    hwts=4
    if [[ ${iso} == "off" ]]; then  iso_num=0;
    else                            iso_num=${iso}; fi
    if [[ ${trf} == "auto" ]]; then last_core=$((${cores}-1-${iso_num}));
    else                            last_core=$((${cores}-1)); fi

    trf=""
    for socket in $(seq 0 $((${sockets}-1))); do
        for core in $(seq 0 $((${cores}-1))); do
            if [[ ${core} -gt ${last_core} ]]; then continue; fi
            for hwt in $(seq 1 ${smt}); do
                cpu=$(((${socket}*${cores}+${core})*${hwts}+${hwt}-1))
                trf+="${cpu} "
            done
        done
    done

    trf=$(echo ${trf} | tr -s " " ",")
fi

if [[ ${djt} != "off" ]]; then
    aflags+="dejitter "
fi

if [[ ${win} != "0:1000000" ]]; then
    aflags+="jwind${win/:/x} "
fi

# Finish LSF options
lsf+=" -alloc_flags \"${aflags}\""

# ------------------------------------------------
# Launch
# ------------------------------------------------

#echo "dir=${dir}"
#echo "nod=${nod}"
#echo "ppn=${ppn}"
#echo "smt=${smt}"
#echo "iso=${iso}"
#echo "lau=${lau}"
#echo "mpi=${mpi}"
#echo "mpc=${mpc}"
#echo "omp=${omp}"
#echo "que=${que}"
#echo "mac=${mac}"
#echo "trc=${trc}"
#echo "trf=${trf}"
#echo "jid=${jid}"
#echo "djt=${djt}"
#echo "win=${win}"
#echo "lsf=${lsf}"

bsub << EOF

    #BSUB -nnodes ${nod}
    #BSUB -o ${dir}/job%J.out
    #BSUB -q ${que}
    #BSUB -J ${jid}(${nod})
    ${mac}
    #BSUB ${lsf}

    # ----------------------------------
    # functions
    # ----------------------------------

    function echosep
    {
        echo
        printf "%80s\n" | tr " " "="
        echo
    }

    function echoval
    {
        cmd=\${*}
        echo "\${cmd}" | tr -s " "
        eval "\${cmd}"
    }

    # ----------------------------------
    # internal variables
    # ----------------------------------

    job=\${LSB_JOBID}
    hosts=\$(echo \${LSB_MCPU_HOSTS} | sed "s/ /\n/g" | grep c699c | sort -u)
    outhosts="${dir}/job\${job}.hosts"
    echo \${hosts} | sed "s/ /\n/g" > \${outhosts}

    env > ${dir}/job\${job}.env

    # ----------------------------------
    # mpirun hostfile
    # ----------------------------------

    hostfile="${dir}/job\${job}.hostfile"
    for host in \${hosts}; do
        echo "\${host} slots=${ppn}" >> \${hostfile}
    done

    if [[ "${lau}" == "mpirun" ]]; then
        echosep
        echo "mpirun setup"

        export USING_MPIRUN=1
        for host in \${hosts}; do
            echoval "sudo ssh \$host \"chmod a+rw /sys/fs/cgroup/cpuset/tasks\""
        done
    fi

    #export FORCE_CPUSET=1

    # ----------------------------------
    # export environment variables
    # ----------------------------------

    export SMT=${smt}
    export RANKS_PER_NODE=${ppn}
    export ISOLATION=${iso}
    export OMP_NUM_THREADS=${omp}
    export OMP_PLACES=threads

    export JTF_JOB=\${job}
    export JTF_DIR="${dir}"

    # ----------------------------------
    # collect pre-run data
    # ----------------------------------

    #echosep
    #echo "collect pre-run data"
    #
    #for host in \${hosts}; do
    #    echoval "ssh \${host} \"${SDIR}/collect-data.sh \
    #        &> ${dir}/job\${job}-\${host}-pre\" &"
    #done
    #wait

    # ----------------------------------
    # Start GPFS option
    # ----------------------------------

    if [[ ${gpf} != "off" ]]; then
        echosep
        echo "Moving mmfsd"
        for host in \${hosts}; do
             echoval "sudo ssh \${host} \"${SDIR}/../cpuset-scripts/cpuset_move \
                /sys/fs/cgroup/cpuset/csm_system /sys/fs/cgroup/cpuset mmfsd\" &"
            #echoval "sudo ssh \${host} \"${SDIR}/../cpuset-scripts/cpuset_move \
            #    /sys/fs/cgroup/cpuset/csm_system                               \
            #    /sys/fs/cgroup/cpuset/allocation_\${CSM_ALLOCATION_ID} mmfsd\" &"
        done
    fi
    wait

    # ----------------------------------
    # start tracing
    # ----------------------------------

    echosep
    if   [[ ${trc} == "off" ]]; then
        echo "tracing off: nothing to do"

    elif [[ ${trc} == "perf" ]]; then
        echo "perf tracing: start"

        for host in \${hosts}; do
            perfdat="/tmp/job\${job}-perf-\${host}.dat"
            perfout="/tmp/job\${job}-perf-\${host}.out"
            cmd="sudo perf record -e sched:*,irq:*,workqueue:*,powerpc:*
                --mmap-pages 2G --cpu ${trf} -k CLOCK_MONOTONIC -a
                -o \${perfdat} &> \${perfout}"
            echoval "ssh \${host} \"\${cmd}\" &"
        done
        echoval "sleep 30"

    else
        echo "tracefs tracing: start"
        for host in \${hosts}; do
            cmd="${SDIR}/tracefs/tracefs-start.sh ${trc}"
            echoval "ssh \${host} \"\${cmd}\" &"
        done
        wait
        echo "tracefs tracing: all set"
    fi

    # ----------------------------------
    # run application
    # ----------------------------------

    echosep
    echoval "${run}"

    # ----------------------------------
    # stop tracing
    # ----------------------------------

    echosep
    if   [[ ${trc} == "off" ]]; then
        echo "tracing off: nothing to do"

    elif [[ ${trc} == "perf" ]]; then
        echo "perf tracing: stop"
        for host in \${hosts}; do
            echoval "sudo ssh \${host} \"killall perf\" &"
        done
        wait
        echoval "sleep 60"
        echo "perf tracing: stopped"

    else
        echo "tracefs tracing: stop"
        for host in \${hosts}; do
            cmd="${SDIR}/tracefs/tracefs-stop.sh"
            echoval "ssh \${host} \"\${cmd}\" &"
        done
        wait
        echo "tracefs tracing: stopped"
    fi

    # ----------------------------------
    # copy and clean traces
    # ----------------------------------

    echosep
    if   [[ ${trc} == "off" ]]; then
        echo "tracing off: nothing to copy"

    elif [[ ${trc} == "perf" ]]; then
        echo "perf tracing: copy and clean"
        for host in \${hosts}; do
            echoval "sudo scp \${host}:/tmp/job\${job}-perf-\${host}.* \
                $(hostname):${tracedir}"
            echoval "sudo ssh \${host} rm /tmp/job\${job}-perf-\${host}.*"
            echoval "sudo ssh $(hostname) \
                chown ${user}:${user} ${tracedir}/job\${job}-perf-\${host}.*"
        done
        echo "perf tracing: done"

    else
        echo "tracefs tracing: collect"
        pdsh -w ^${dir}/job\${job}.hosts "     \
            ${SDIR}/tracefs/tracefs-collect.sh \
            \${job} ${dir} ${tracedir}"
        echo "tracefs tracing: collect done"

        echo
        echo "tracefs tracing: move"
        pdsh -w ^${dir}/job\${job}.hosts "                   \
            scp ${tracedir}/job\${job}* $(hostname):${dir}/; \
            rm  ${tracedir}/job\${job}*"
        echo "tracefs tracing: move done"

        echo
        echo "tracefs tracing: clean"
        for host in \${hosts}; do
            cmd="${SDIR}/tracefs/tracefs-reset.sh"
            echoval "ssh \${host} \"\${cmd}\" &"
        done
        wait
        echo "tracefs tracing: clean done"
    fi

    # ----------------------------------
    # collect post-run data
    # ----------------------------------

    #echosep
    #echo "collect post-run data"
    #
    #for host in \${hosts}; do
    #    echoval "ssh \${host} \"${SDIR}/collect-data.sh
    #        &> ${dir}/job\${job}-\${host}-post\" &"
    #done
    #wait

    # ----------------------------------
    # stop GPFS option
    # ----------------------------------

    echosep
    echo "Check GPFS processes"
    echo
    for host in \${hosts}; do
        echoval "ssh \${host} \"ps -eo comm,cgroup | grep mmfsd\" &"
    done
    wait

    if [[ ${gpf} != "off" ]]; then
        echosep
        echo "Reset GPFS config"
        echo
        for host in \${hosts}; do
             echoval "sudo ssh \${host} \"${SDIR}/../cpuset-scripts/cpuset_move \
                /sys/fs/cgroup/cpuset /sys/fs/cgroup/cpuset/csm_system mmfsd\" &"
            #echoval "sudo ssh \${host} \"${SDIR}/../cpuset-scripts/cpuset_move \
            #    /sys/fs/cgroup/cpuset/allocation_\${CSM_ALLOCATION_ID}         \
            #    /sys/fs/cgroup/cpuset/csm_system mmfsd\" &"
        done
    fi
    wait

    echosep
    echo "Check GPFS processes"
    echo
    for host in \${hosts}; do
        echoval "ssh \${host} \"ps -eo comm,cgroup | grep mmfsd\" &"
    done
    wait

EOF

