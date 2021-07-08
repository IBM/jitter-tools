#!/bin/bash

rank=$PMIX_RANK
cpu=$(cat /proc/self/status | grep Cpus_allowed_list | awk '{print $2}' | cut -d '-' -f1 | cut -d ',' -f1)
cpus=$(cat /proc/self/status | grep Cpus_allowed_list | awk '{print $2}')
host=$HOSTNAME
printf '[check] rank=%03d cpu=%03d host=%s cpus=%s\n' "$rank" "$cpu" "$host" "$cpus"

