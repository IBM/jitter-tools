#!/bin/bash

#
#   IBM Corporation (C) 2018
#   Bryan Rosenburg -- rosnbrg@us.ibm.com
#   Nelson Mimura -- nmimura@us.ibm.com
#
#   Set up cpu pinning for WSC experiments
#

# ------------------------------------------------
# Read application to be executed
# ------------------------------------------------

app=$*

# ------------------------------------------------
# Check for input env vars
# ------------------------------------------------

if [ -z "$PMIX_RANK" ]; then
    echo "[cpuset-wsc] PMIX_RANK not set; stop"
    exit 1
fi

if [ -z "$SMT" ]; then
    echo "[cpuset-wsc] SMT not set; stop"
    exit 1
fi

if [ -z "$RANKS_PER_NODE" ]; then
    echo "[cpuset-wsc] RANKS_PER_NODE not set; stop"
    exit 1
fi

if [ -z "$ISOLATION" ]; then
    echo "[cpuset-wsc] ISOLATION not set; stop"
    exit 1
fi

# ------------------------------------------------
# Determine local rank
# ------------------------------------------------

world_rank=$PMIX_RANK
local_rank=$(expr $world_rank % $RANKS_PER_NODE)
export LOCAL_RANK=$local_rank

# ------------------------------------------------
# Base CPU parameters
# ------------------------------------------------

sockets=2
smt_max=4
cores_per_socket=22

# ------------------------------------------------
# Choose a cpu
# This algorithm distributes ranks uniformly
# across available cpus.
# ------------------------------------------------
if [[ "$SMT" != "1" && "$SMT" != "4" ]]; then 
    echo "[cpuset-wsc] error: invalid smt value: $SMT"
    exit 1
fi
iso=$ISOLATION
if   [[ ${iso} == "off" || ${iso} == "n" ]]; then
    iso=0
elif [[ ${iso} == "on"  || ${iso} == "y" ]]; then
    iso=1
fi
avail_threads_per_socket=$(( ($cores_per_socket - $iso) * $SMT ))
avail_threads=$(( $sockets * $avail_threads_per_socket ))
logical_thread=$(( ($local_rank * $avail_threads) / $RANKS_PER_NODE ))
if [[ $(( $logical_thread >= $avail_threads_per_socket )) != 0 ]]; then
    # Skip the core-isolation "gap" between sockets
    logical_thread=$(( $logical_thread + ($iso * $SMT) ))
fi
#
# @nmimura 07/12/19: change this to reflect the new core
# allocation for system activities; before the system cores 
# were the last ones from each socket, now we use the first ones.
cpu=$(( ($logical_thread * $smt_max) / $SMT ))
cpu=$(( $cpu + $iso * 4  ))

# ------------------------------------------------
# Force cpuset (for mpirun setup or faulty nodes)
# ------------------------------------------------

if [[ "${USING_MPIRUN}" == 1 || "${FORCE_CPUSET}" == 1 ]]; then
    echo $$ > /sys/fs/cgroup/cpuset/tasks
fi

# ------------------------------------------------
# Run command
# ------------------------------------------------

cmd="/usr/bin/taskset -c $cpu $app"
#printf "[cpuset-wsc] HOST=%s WORLD=%03d LOCAL=%03d SMT=%d CPU=%03d CPUSET=%s\n" "${HOSTNAME}" "$world_rank" "$local_rank" "$SMT" "$cpu" "$(cat /proc/self/cpuset)"
eval $cmd

if [[ $? -ne 0 ]]; then
    printf "[cpuset-wsc] error: HOST=%s WORLD=%03d LOCAL=%03d SMT=%d CPU=%03d CPUSET=%s\n" "${HOSTNAME}" "$world_rank" "$local_rank" "$SMT" "$cpu" "$(cat /proc/self/cpuset)"
fi

