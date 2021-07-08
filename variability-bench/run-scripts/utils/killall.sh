#!/bin/bash

for job in $(bjobs | grep "PEND\|PSUSP" | cut -d " " -f 1); do
    bkill ${job}
done

