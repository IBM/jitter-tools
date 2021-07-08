#!/bin/bash

#
#   IBM Corporation (C) 2018
#   Nelson Mimura -- nmimura@us.ibm.com
#
#   Start tracing using tracefs
#

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
# Variables and arguments
# ------------------------------------------------

ROOT=/sys/kernel/debug/tracing
FTRACE=${ROOT}/instances/ftrace
TRACEP=${ROOT}/instances/tracep
EVENTS="sched:* irq:* workqueue:* powerpc:*"

function usage {
    echo "usage: ${0} TRC"
    echo "  TRC: tracing option passed to run-script"
    echo
    exit 1
}

trc=$1
ftrace="0"
tracep="0"
if [[ ${trc} == *f* ]]; then
    ftrace=1
fi
if [[ ${trc} == *p* ]]; then
    tracep=1
fi

# ------------------------------------------------
# Completetly reset tracing
# ------------------------------------------------

#${SDIR}/tracefs-reset.sh

# ------------------------------------------------
# Start tracing
# ------------------------------------------------

if [[ ${ftrace} == 1 ]]; then
    #mkdir -p ${FTRACE}
    echo function > ${FTRACE}/current_tracer
    echo          > ${FTRACE}/set_event
    echo mono     > ${FTRACE}/trace_clock
    echo 100000   > ${FTRACE}/buffer_size_kb
    echo 1        > ${FTRACE}/tracing_on
fi

if [[ ${tracep} == 1 ]]; then
    #mkdir -p ${TRACEP}
    echo nop       > ${TRACEP}/current_tracer
    echo ${EVENTS} > ${TRACEP}/set_event
    echo mono      > ${TRACEP}/trace_clock
    echo 100000    > ${TRACEP}/buffer_size_kb
    echo 1         > ${TRACEP}/tracing_on
fi

