#!/bin/bash

if [ $# -eq 0 ]; then
  echo "Usage: ./script *.rpm"
  exit 0
fi

rpm -e dysc
rpm -i $1
systemctl daemon-reload

echo "Install csm?"
read x
if [ "x${x}" = "xy" ]; then
  echo "Installing csm!"
  cp csm/dyncgroup_prolog /shared/lsf-csm/csm_cn/prologs/WSC/dyncgroup_prolog
  cp csm/dyncgroup_epilog /shared/lsf-csm/csm_cn/prologs/WSC/dyncgroup_epilog
fi


# /opt/ibm/csm/prologs/privileged_prolog
#
#                elif "DYNCGROUP" == fparse.group(1):
#                    cmd = ['/usr/bin/python', '%s/dyncgroup_prolog' % (JOB_SCRIPT_DIR_WSC), feature.upper()]
#                    runcmd (cmd)
#

# /opt/ibm/csm/prologs/privileged_epilog
#
#                elif "DYNCGROUP" == fparse.group(1):
#                    cmd = ['/usr/bin/python', '%s/dyncgroup_epilog' % (JOB_SCRIPT_DIR_WSC), feature.upper()]
#                    runcmd (cmd)

