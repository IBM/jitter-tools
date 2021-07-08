#!/bin/bash

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
# Echo and execute command
# ------------------------------------------------

function echoval
{
  cmd=${*}
  echo ${cmd} | tr -s " "
  eval ${cmd}
}

# ------------------------------------------------
# Process arguments
# ------------------------------------------------

ROOT="/sys/kernel/debug/tracing"
FTRACE="${ROOT}/instances/ftrace"
TRACEP="${ROOT}/instances/tracep"
FANOUT=4

if [[ ${#} != 3 ]]; then
  echo "usage: ${0} job dir trd"
  echo "  job   job id"
  echo "  dir   job directory"
  echo "  trd   trace directory"
  exit 1
fi

job="${1}"
dir="${2}"
trd="${3}"

# ------------------------------------------------
# Filtered collect function
# ------------------------------------------------

function filtered_collect
{
  if [[ ${#} -ne 5 ]]; then
    echo "error: filtered_collect(): missing parameters"
    exit 1
  fi

  local hst=${1}
  local cpu=${2}
  local trd=${3}
  local job=${4}
  local cnt=${5}

  cpu_pad=$(printf "%03d" "${cpu}")

  if [[ -d ${FTRACE} ]]; then
    ftrace="${trd}/job${job}-${HOSTNAME}-cpu${cpu_pad}.ftrace"
    echoval "cp ${FTRACE}/per_cpu/cpu${cpu}/trace ${ftrace}.raw"
    echoval "${SDIR}/filter/tracefs-filter ${jtf} ${cnt} ${hst} ${cpu} \
      ${ftrace}.raw ${ftrace}"
    rm "${ftrace}.raw"
  fi

  if [[ -d ${TRACEP} ]]; then
    tracep="${trd}/job${job}-${HOSTNAME}-cpu${cpu_pad}.tracep"
    echoval "cp ${TRACEP}/per_cpu/cpu${cpu}/trace ${tracep}.raw"
    echoval "${SDIR}/filter/tracefs-filter ${jtf} ${cnt} ${hst} ${cpu} \
      ${tracep}.raw ${tracep}"
    rm "${tracep}.raw"
  fi
}

# ------------------------------------------------
# Copy traces from tracefs and filter if possible
# ------------------------------------------------

mkdir -p ${trd}
jtf="${dir}/job${job}.jtf"

# ----------------------------
# Filtering enabled
# ----------------------------

if [[ -f ${jtf} ]]; then

  # Get number of JTF entries
  cnt=$(cat ${jtf} | wc -l)

  # Generate list of cpus
  cpus=$(cat ${jtf} | grep ${HOSTNAME} | tr -s " " | cut -d " " -f 3 | uniq)  

  # Copy traces per cpu and run filter program
  counter=0
  for cpu in ${cpus}; do
    filtered_collect ${HOSTNAME} ${cpu} ${trd} ${job} ${cnt} &

    counter=$(( ${counter} + 1))
    if [[ $(( ${counter} % ${FANOUT} )) -eq 0 ]]; then
      wait
    fi
  done

# ----------------------------
# Filtering disabled
# ----------------------------

else

  # Copy consolidated trace file
  if [[ -d ${FTRACE} ]]; then
    echoval "cp ${FTRACE}/trace ${trd}/job${job}-${HOSTNAME}.ftrace"
  fi
  if [[ -d ${TRACEP} ]]; then
    echoval "cp ${TRACEP}/trace ${trd}/job${job}-${HOSTNAME}.tracep"
  fi

fi

