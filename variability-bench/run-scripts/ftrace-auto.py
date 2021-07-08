#!/usr/bin/python

import os
import sys

sdir   = os.path.dirname(os.path.realpath(__file__))
ftproc = '%s/../jitter-trace/ftrace/jt_ftproc.py' % (sdir)
th     = '0.001'

if len(sys.argv[1:]) < 2:
    print 'usage: %s ftdir jtf' % (sys.argv[0])
    print '  ftdir  directory containing traces'
    print '  jtf    path to job jtf file'
    exit(1)

ftdir = sys.argv[1]

for match in sys.argv[2:]:
    if 'jtf' not in match or '.pro' in match:
        continue

    f = open(match, 'r')
    f.readline()
    t0 = f.readline().strip()
    t1 = f.readline().strip()
    f.close()
 
    cmd  = '%s ' % (ftproc)
    cmd += '-f %s/%s.ftrace ' % (ftdir, match.split('/')[-1])
    cmd += '-s %s -e %s -t %s ' % (t0, t1, th)
    cmd += '-g -q -u -r '
    cmd += '> %s.pro' % (match)

    print cmd
    os.system(cmd)

