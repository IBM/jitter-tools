#!/bin/bash

#
#   IBM Corporation (C) 2018
#   Nelson Mimura -- nmimura@us.ibm.com
#
#   Reset tracefs configuration
#

ROOT=/sys/kernel/debug/tracing
FTRACE=${ROOT}/instances/ftrace
TRACEP=${ROOT}/instances/tracep

#echo 0   > ${ROOT}/tracing_on
#echo nop > ${ROOT}/current_tracer
#echo     > ${ROOT}/set_event
#echo     > ${ROOT}/trace
#echo 100 > ${ROOT}/buffer_size_kb

if [[ -d ${FTRACE} ]]; then
    echo 0   > ${FTRACE}/tracing_on
    echo nop > ${FTRACE}/current_tracer
    echo     > ${FTRACE}/set_event
    echo     > ${FTRACE}/trace
    echo 100 > ${FTRACE}/buffer_size_kb
    #rmdir ${FTRACE}
fi

if [[ -d ${TRACEP} ]]; then
    echo 0   > ${TRACEP}/tracing_on
    echo nop > ${TRACEP}/current_tracer
    echo     > ${TRACEP}/set_event
    echo     > ${TRACEP}/trace
    echo 100 > ${TRACEP}/buffer_size_kb
    #rmdir ${TRACEP}
fi

