#!/usr/bin/python

import numpy as np
import os
import sys

try:
    dir = sys.argv[1]
except IndexError:
    print 'usage: %s dir' % (sys.argv[0])
    print '  dir    directory containing .pro files'
    exit(1)

pros = []

for item in os.listdir(dir):
    if item.endswith('.pro'):
        pros.append(os.path.join(dir, item))

def process_stats(stats):

    THRESH   = 0.000001
    fsum_evs = {} # filtered events
    fsum_all = {} # filtered sum
    fsum_ufe = {} # unfiltered events
    fsum_ufs = {} # unfiltered sum

    for pro in pros:
        f = open(pro, 'r')
        for l in f:
            if stats in l:
                break
        for l in f:
            if l.startswith('-'):
                break
        for l in f:
            if not l.strip():
                break
            split = filter(None, l.split(' '))
            fname = split[0]
            fsum  = float(split[2])
            if fname not in fsum_ufe:
                fsum_ufe[fname] = []
                fsum_ufs[fname] = 0.0
            fsum_ufe[fname].append(fsum)
            fsum_ufs[fname] += fsum
            if fsum < THRESH:
                continue
            if fname not in fsum_evs:
                fsum_evs[fname] = []
                fsum_all[fname] = 0.0
            fsum_evs[fname].append(fsum)
            fsum_all[fname] += fsum
        f.close()

    HLINE = 110 * '='

    print
    print stats
    print HLINE
    print
    print '%-30s %16s %16s %8s %8s %8s %8s' % (
            'function', 
            'sum (s)',
            'avg(s)', 
            'count',
            '>100us',
            '>1ms',
            'traces')
    print 100 * '-'
    for fname in sorted(fsum_all, key=fsum_all.get, reverse=True):
        print '%-30s %16.6f %16.6f %8d %8d %8d %8d' % (
                fname[:30],
                fsum_all[fname], 
                np.mean(fsum_evs[fname]),
                len(fsum_evs[fname]), 
                len([i for i in fsum_evs[fname] if 0.0001 < i and i < 0.001]),
                len([i for i in fsum_evs[fname] if 0.001 < i]),
                len(pros))

    print
    print 'HISTOGRAMS OF', stats
    print HLINE
    for fname in sorted(fsum_all, key=fsum_all.get, reverse=True):
        print
        print fname
        print '  min=%.6f max=%.6f avg(unf)=%.6f count(unf)=%d' % (
                min(fsum_ufe[fname]), 
                max(fsum_ufe[fname]), 
                np.mean(fsum_ufe[fname]),
                len(fsum_ufe[fname]))
        hist, bins = np.histogram(
                fsum_ufe[fname],
                bins=[1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1E0, float('inf')])
        for i in range(len(hist)):
            print '  %10.6f %10.6f -- %d' % (\
                    bins[i], 
                    bins[i+1],
                    hist[i])

process_stats('GLOBAL FUNCTION STATISTICS')
process_stats('GLOBAL ROOT FUNCTION STATISTICS')

