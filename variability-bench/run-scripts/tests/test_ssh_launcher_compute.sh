#!/bin/bash

if [[ $# -lt 1 ]]; then
    echo "usage: $0 hosts..."
    echo "  hosts...    list of hosts to test"
    echo
    exit 1
fi

hosts=${*}

function test_host {
    host=${1}
    rc=0

    ssh ${host} hostname &> /dev/null
    if [[ $? -ne 0 ]]; then
        echo "[test_ssh_launcher_compute] failed: '${host}'" 1>&2
        rc=1
    fi

    echo ${rc}
}

error=0
for host in ${hosts}; do
    let error=${error}+$(test_host ${host} &)
done
wait

if [[ ${error} -eq 0 ]]; then
    echo "[test_ssh_launcher_compute] success!"
fi

exit ${error}

