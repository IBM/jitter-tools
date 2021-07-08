#!/usr/bin/python

import os
import sys

# ------------------------------------------------------------------------------

class var: pass

#
#   Program parameters
#

var.dir  = None
var.sort = 'n'

#
#   Program variables
#

var.exps = []

# ------------------------------------------------------------------------------

def main(args):
    parse(args)
    check()
    process()
    show()

# ------------------------------------------------------------------------------

#
#   Parse command line arguments
#

def parse(args):
    try:
        var.dir = parse_dir(args[0]);
    except IndexError:
        print >> sys.stderr, 'error: missing arguments'
        usage()

#
#   Print usage information and exit
#

def usage():
    print >> sys.stderr, 'usage: %s DIR' % (sys.argv[0])
    print >> sys.stderr, '  DIR     directory containing outputs'
    exit(1)

#
#   Parse input directory option
#

def parse_dir(arg):
    if not os.path.exists(arg):
        print 'error: could not find directory: "%s"' % (arg)
        usage
    return arg

# 
#   Parse sort option
#   

def parse_sort(arg):
    return arg

#
#   Check program parameters
#

def check():
    if any(i is None for i in [var.dir]):
        print >> sys.stderr, 'error: missing required arguments'
        usage()

# ------------------------------------------------------------------------------

#
#   Class to hold experiment (output) information
#

class Exp:
    def __init__(self):
        self.np      = '-' 
        self.ppn     = '-' 
        self.nrs     = '-' 
        self.smt     = '-' 
        self.mxm     = '-'
        self.p       = '-' 
        self.c       = '-' 
        self.barrier = '-' 
        self.extra   = '-' 
        self.raw     = '-' 
        self.trace   = '-' 
        self.ftrace  = '-' 
        self.elapsed = '-' 
        self.slow    = '-' 
        self.t_comp  = '-' 
        self.t_jit   = '-' 
        self.t_comm  = '-' 
        self.cj_max  = '-' 
        self.cj_var  = '-' 
        self.ph_max  = '-' 
        self.ph_var  = '-' 
        self.when    = '-'
        self.done    = False

#
#   Process input directory
#

def process():
    for i in os.listdir(var.dir):
        if i.startswith('job') and i.endswith('.out') and 'jt' not in i:
            path = os.path.join(var.dir, i)
            process_file(path)

#
#   Process experiment file
#

def process_file(path):
    f   = open(path, 'r')
    exp = Exp()
    exp.job = path.split('job')[1].split('.out')[0]
    for l in f:
        l = l.strip()
        if l.startswith('/opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun'):
            exp.np      = int(l.split('--np ')[1].split(' ')[0])
            exp.ppn     = int(l.split('--rs_per_host ')[1].split(' ')[0])
            exp.nrs     = int(l.split('--nrs ')[1].split(' ')[0])
            exp.n       = exp.nrs / exp.ppn
            exp.smt     = exp.np  / exp.nrs
            exp.p       = int(l.split('-p ')[1].split(' ')[0])
            exp.c       = int(l.split('-c ')[1].split(' ')[0])
            exp.barrier = '-r'  in l.split(' ')
            exp.extra   = '-e'  in l.split(' ')
            exp.raw     = '-d'  in l.split(' ')
            exp.ftrace  = '-t'  in l.split(' ')
            exp.mxm     = 'mxm' in l
        elif l.startswith('Runtime elapsed'):
            exp.elapsed = l.strip().split(' ')[-2]
        elif l.startswith('Slowdown'):
            exp.slow    = '%.4f' % (float(l.strip().split(' ')[-1][:-1]))
        elif l.endswith('cpu time is compute'):
            exp.t_comp  = l.strip().split(' ')[0]
        elif l.endswith('cpu time is jitter'):
            exp.t_jit   = l.strip().split(' ')[0]
        elif l.endswith('cpu time is communication'):
            exp.t_comm  = l.strip().split(' ')[0]
        elif l.startswith('Job was executed on'):
            exp.when    = l.strip().split('at')[1]
    f.close()

    f = open(path, 'r')
    try:
        data       = f.read()
        cj         = data.split('Compute and Jitter histogram')[1]
        exp.cj_max = cj.split('max: ')[1].split(' ')[0]
        ph         = data.split('Phase histogram')[1]
        exp.ph_max = ph.split('max: ')[1].split(' ')[0]
        exp.done   = True
    except:
        pass
    f.close()

    var.exps.append(exp)

# ------------------------------------------------------------------------------

def show():
    fmt  = '%6s %2d %6d %6d %4d %4d %3d'
    fmt += '%7d %6d %3d %1d %1d %1d '
    fmt += '%10s %10s %10s %10s %10s '
    fmt += '%12s %12s '
    fmt += '%28s'
    header = fmt.replace('d', 's') % (\
            'job', 'ok', 'n', 'np', 'ppn', 'smt', 'mxm',
            'p', 'c', 'r', 'e', 'd', 't',
            'elapsed', 'slowdown', 'comp', 'jit', 'comm',
            'c+j max', 'phase max',
            'executed at')
    print len(header) * '-'
    print header
    print len(header) * '-'
    for exp in sorted(var.exps, 
            key=lambda x: (x.p, x.c, x.n, x.slow), 
            reverse=True):
        if exp.done:
            print fmt % (\
                    exp.job, exp.done, exp.n, exp.np, exp.ppn, exp.smt, exp.mxm,
                    exp.p, exp.c, exp.barrier, exp.extra, exp.raw, exp.ftrace, 
                    exp.elapsed, exp.slow, exp.t_comp, exp.t_jit, exp.t_comm,
                    exp.cj_max, exp.ph_max, 
                    exp.when)
    print len(header) * '-'
    for exp in var.exps:
        if not exp.done:
            print fmt % (\
                    exp.job, exp.done, exp.n, exp.np, exp.ppn, exp.smt, exp.mxm,
                    exp.p, exp.c, exp.barrier, exp.extra, exp.raw, exp.ftrace,
                    exp.elapsed, exp.slow, exp.t_comp, exp.t_jit, exp.t_comm,
                    exp.cj_max, exp.ph_max, 
                    exp.when)

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1:])

