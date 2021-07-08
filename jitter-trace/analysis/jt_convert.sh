#!/bin/bash

# ==================================== #
# IBM Corporation (C) 2017             #
# Nelson Mimura Gonzalez               #
# nmimura@us.ibm.com                   #
# ==================================== #

for perfdat in "$@"; do

    perftml=${perfdat/".dat"/".tml"} # Textual timeline.
    perfper=${perfdat/".dat"/".per"} # Errors during conversion.

    cmd="perf script -i ${perfdat} -v --ns > ${perftml} 2> ${perfper} &"
    echo "[convert] ${cmd}"
    eval ${cmd}

done

wait

