#!/bin/bash

GPFS=/gpfs/alpinetds/ven201/scratch/
stat ${GPFS} &> /dev/null
echo "gpfs_check ${HOSTNAME} ${?}"

