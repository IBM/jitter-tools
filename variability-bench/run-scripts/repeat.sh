#!/bin/bash

if [[ $# -lt 2 ]]; then
    echo "usage: $0 reps cmd"
    echo "  reps: repetitions"
    echo "  cmd: command to run"
    exit 1
fi

reps=$1; shift
cmd=$*

for i in $(seq 1 ${reps}); do
    echo ${cmd}
    eval ${cmd}
done

