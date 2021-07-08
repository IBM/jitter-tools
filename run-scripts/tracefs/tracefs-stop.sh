#!/bin/bash

#
#   IBM Corporation (C) 2018
#   Nelson Mimura -- nmimura@us.ibm.com
#
#   Stop tracefs tracing
#

ROOT=/sys/kernel/debug/tracing
FTRACE=${ROOT}/instances/ftrace
TRACEP=${ROOT}/instances/tracep

#echo 0 > ${ROOT}/tracing_on

if [[ -d ${FTRACE} ]]; then 
    echo 0 > ${FTRACE}/tracing_on
fi

if [[ -d ${TRACEP} ]]; then
    echo 0 > ${TRACEP}/tracing_on
fi

