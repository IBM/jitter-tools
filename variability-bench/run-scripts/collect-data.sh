#!/bin/bash

function box
{
    title=$1
    printf "%80s\n" | tr " " "="
    echo ${title}
    printf "%80s\n" | tr " " "="
}

function test_versions
{
    echo
    box "VERSIONS"
    echo

    /opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun --version | head -n 1
    /opt/ibm/spectrum_mpi/bin/mpirun --version | head -n 1
}

function test_env
{
    echo
    box "ENVIRONMENT"

    echo
    env
}

function test_ib_counters
{
    echo
    box "INFINIBAND COUNTERS"

    root="/sys/class/infiniband"
    for dev in $(ls ${root}); do
        echo
        for dir in counters hw_counters; do
            path="${root}/${dev}/ports/1/${dir}"
            for counter in $(ls ${path}); do
                counter_path="${path}/${counter}"
                printf "%-80s %s\n" "${counter_path}" "$(cat ${counter_path})"
            done
        done
    done
}

test_versions
test_env
test_ib_counters

