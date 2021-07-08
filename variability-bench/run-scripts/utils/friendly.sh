#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

function friendly_lim
{
    echo 1
}

function friendly_get
{
    a=$(bjobs 2> /dev/null | grep JITTER | wc -l)
    sleep 1
    b=$(bjobs 2> /dev/null | grep JITTER | wc -l)
    if [[ ${a} -eq ${b} ]]; then
        echo ${a}
    else
        return 1
    fi
}

function friendly_lst
{
    cmdlst[0]="${SDIR}/../run-wsc.sh -d /gpfs/wscgpfs01/nmimura/work/jitter-web/exps/32-mpirun-pami-libcoll -n 32 -l mpirun                 -j JITTER -x ../jitter-bench/jitter-bench -p 60000 -c 5000 -n -e"
    cmdlst[1]="${SDIR}/../run-wsc.sh -d /gpfs/wscgpfs01/nmimura/work/jitter-web/exps/32-jsrun-pami-libcoll  -n 32 -l  jsrun                 -j JITTER -x ../jitter-bench/jitter-bench -p 60000 -c 5000 -n -e"
    cmdlst[2]="${SDIR}/../run-wsc.sh -d /gpfs/wscgpfs01/nmimura/work/jitter-web/exps/32-mpirun-pami-hcoll   -n 32 -l mpirun        -c hcoll -j JITTER -x ../jitter-bench/jitter-bench -p 60000 -c 5000 -n -e"
    cmdlst[3]="${SDIR}/../run-wsc.sh -d /gpfs/wscgpfs01/nmimura/work/jitter-web/exps/32-jsrun-pami-hcoll    -n 32 -l  jsrun        -c hcoll -j JITTER -x ../jitter-bench/jitter-bench -p 60000 -c 5000 -n -e"
    cmdlst[4]="${SDIR}/../run-wsc.sh -d /gpfs/wscgpfs01/nmimura/work/jitter-web/exps/32-mpirun-mxm-libcoll  -n 32 -l mpirun -k mxm          -j JITTER -x ../jitter-bench/jitter-bench -p 60000 -c 5000 -n -e"
    cmdlst[5]="${SDIR}/../run-wsc.sh -d /gpfs/wscgpfs01/nmimura/work/jitter-web/exps/32-jsrun-mxm-libcoll   -n 32 -l  jsrun -k mxm          -j JITTER -x ../jitter-bench/jitter-bench -p 60000 -c 5000 -n -e"
    cmdlst[6]="${SDIR}/../run-wsc.sh -d /gpfs/wscgpfs01/nmimura/work/jitter-web/exps/32-mpirun-mxm-hcoll    -n 32 -l mpirun -k mxm -c hcoll -j JITTER -x ../jitter-bench/jitter-bench -p 60000 -c 5000 -n -e"
    cmdlst[7]="${SDIR}/../run-wsc.sh -d /gpfs/wscgpfs01/nmimura/work/jitter-web/exps/32-jsrun-mxm-hcoll     -n 32 -l  jsrun -k mxm -c hcoll -j JITTER -x ../jitter-bench/jitter-bench -p 60000 -c 5000 -n -e"
}

friendly_lst
idx=0

while :; do
    lim=$(friendly_lim)
    if [[ $? -ne 0 ]]; then continue; fi
    num=$(friendly_get)
    if [[ $? -ne 0 ]]; then continue; fi

    echo "$(date +'%Y-%m-%d %H:%M:%S') --> ${num}/${lim}"

    if [[ ${num} -lt ${lim} ]]; then
        cmd=${cmdlst[${idx}]}
        idx=$(((${idx}+1)%(${#cmdlst[@]})))
        echo "RUN  ${cmd}"
        eval ${cmd}
    fi

    sleep 1
done

