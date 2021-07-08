#!/bin/bash

tmp=".watch"

function log {
    cmd=${*}
    eval ${cmd} &>> ${tmp}
}

function calls {
    date
    echo
    bjobs
}

function show {
    clear
    cat ${tmp}
    > ${tmp}
}

while :; do
    log calls
    show
    sleep 1
done

