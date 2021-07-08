#!/usr/bin/python

import numpy as np
import os
import sys

# ------------------------------------------------------------------------------

DETAIL_NONE  = 0
DETAIL_PLAIN = 1
DETAIL_GRAPH = 2

# ------------------------------------------------------------------------------

class var: pass

var.file       = None
var.start      = None
var.end        = None
var.thresh     = None

var.top        = 3
var.detail     = DETAIL_NONE
var.quiet      = False
var.superquiet = False
var.funstat    = False
var.rootstat   = False

var.last       = None  # last event read from trace
var.first      = True  # first event from trace
var.started    = False # jitter started
var.app        = []    # app events
var.jit        = []    # jit events
var.jit_accum  = 0     # jit time accumulator
var.jit_start  = []    # jit start timestamps
var.jit_end    = []    # jit end timestamps
var.jit_seq    = []    # jit current function sequence
var.jit_seqs   = []    # jit sequences
var.funstats   = {}    # global function statistics
var.rootstats  = {}    # global root function statistics

# ------------------------------------------------------------------------------

def main(args):
    setup(args)
    process()
    show()

def setup(args):
    parse_args(args)
    check_var()
    show_var()

def process():
    process_start()
    process_file()
    process_end()

def show():
    show_app()
    show_jit()

# ------------------------------------------------------------------------------

def parse_args(args):
    try:
        i = 0
        while i < len(args):
            if   args[i] in ['-h', '--help']:
                usage()
            elif args[i] in ['-f', '--file']:
                var.file     = parse_file(args[i+1]); i += 2
            elif args[i] in ['-s', '--start']:
                var.start    = parse_time(args[i+1]); i += 2
            elif args[i] in ['-e', '--end']:
                var.end      = parse_time(args[i+1]); i += 2
            elif args[i] in ['-t', '--thresh']:
                var.thresh   = parse_time(args[i+1]); i += 2
            elif args[i] in ['-p', '--top']:
                var.top      = int(args[i+1]);        i += 2
            elif args[i] in ['-d', '--detailed']:
                var.detail   = DETAIL_PLAIN;          i += 1
            elif args[i] in ['-g', '--graph']:
                var.detail   = DETAIL_GRAPH;          i += 1
            elif args[i] in ['-q', '--quiet']:
                var.quiet    = True;                  i += 1
            elif args[i] in ['-u', '--funstat']:
                var.funstat  = True;                  i += 1
            elif args[i] in ['-r', '--rootstat']:
                var.rootstat = True;                  i += 1
            else:
                print >> sys.stderr, 'error: invalid argument: "%s"' % (args[i])
                usage()
    except IndexError:
        print >> sys.stderr, 'error: missing argument value: "%s"' % (args[i])
        usage()

def usage():
    print >> sys.stderr
    print >> sys.stderr, 'usage: %s REQUIRED [OPTIONAL]' % (sys.argv[0])
    print >> sys.stderr
    print >> sys.stderr, 'required:'
    print >> sys.stderr, '  -f FILE    path to trace file'
    print >> sys.stderr, '  -s TIME    timestamp to start analysis'
    print >> sys.stderr, '  -e TIME    timestamp to end analysis'
    print >> sys.stderr, '  -t TIME    threshold to identify aplication time'
    print >> sys.stderr, 'optional:'
    print >> sys.stderr, '  -p TOP      number of top jitter events to show'
    print >> sys.stderr, '  -d          show information for top jitter events'
    print >> sys.stderr, '  -g          create call graphs'
    print >> sys.stderr, '  -q          suppress graphs (function stats only)'
    print >> sys.stderr, '  -u          collect all function statistics'
    print >> sys.stderr
    exit(1)

def parse_file(arg):
    if not os.path.exists(arg):
        print >> sys.stderr, 'error: could not find trace file: "%s"' % (arg)
        print >> sys.stderr, 'check path and file does exist'
        exit(1)
    return arg

def parse_time(arg):
    try:
        split = arg.split('.')
        sec   = split[0]
        nsec  = split[1] + (9 - len(split[1])) * '0'
        return int('%s%s' % (sec, nsec))
    except ValueError:
        print >> sys.stderr, 'error: invalid timestamp: "%s"' % (arg)
        print >> sys.stderr, 'invalid characters (e.g., letters)?'
        exit(1)
    except IndexError:
        print >> sys.stderr, 'error: malformed timestamp: "%s"' % (arg)
        print >> sys.stderr, 'missing fractional part?'
        exit(1)

def check_var():
    if any(i is None for i in [var.file,
                               var.start,
                               var.end,
                               var.thresh]):
        print >> sys.stderr, 'error: missing required arguments'
        usage()

def show_var():
    print
    print '#jt_ftproc -- process ftrace output'
    print '==================================='
    print
    print 'file:      %s' % (var.file)
    print 'start:     %s' % (time_i2s(var.start))
    print 'end:       %s' % (time_i2s(var.end))
    print 'threshold: %s' % (time_i2s(var.thresh))
    print 'print top: %d' % (var.top)
    print 'detail:    %s' % (detail_i2s(var.detail))

def detail_i2s(detail):
    if detail == DETAIL_NONE:
        return 'none'
    if detail == DETAIL_PLAIN:
        return 'plain'
    if detail == DETAIL_GRAPH:
        return 'graph'

# ------------------------------------------------------------------------------

class Event:
    def __init__(self, comm, pid, time, function, parent):
        self.comm     = comm
        self.pid      = pid
        self.time     = time
        self.function = function
        self.parent   = parent

# ------------------------------------------------------------------------------

def process_file():
    fd = open(var.file, 'r')
    for line in fd:
        if not line.startswith('#'):
            process_line(line)
    fd.close()

def process_line(line):
    event = parse_line(line)
    if var.start <= event.time and event.time <= var.end:
        process_event(event)

def parse_line(line):
    split    = line.strip().split('[')
    sub      = split[0].split('-')
    comm     = '-'.join(sub[:-1])
    pid      = sub[-1]
    split    = split[1].split(']')[1].split(' ')
    time     = parse_time(split[2][:-1])
    function = split[3]
    parent   = split[4][2:]
    return Event(comm, pid, time, function, parent)

def process_event(event):
    if var.last:
        process_delta(event)
    var.last = event

def process_delta(event):
    delta = event.time - var.last.time
    if delta > var.thresh or var.first:
        process_app(event, delta)
    else:
        process_jit(event, delta)

# ------------------------------------------------------------------------------

def process_app(event, delta):
    if not var.first:
        jit_end()
    var.first = False
    app(delta)
    jit_start(event)

def process_jit(event, delta):
    jit_inc(event, delta)

# ------------------------------------------------------------------------------

def jit_start(event):
    var.jit_start.append(event.time)
    var.jit_seq.append((event.time, event.function, event.parent))
    var.started = True

def jit_inc(event, delta):
    var.jit_accum += delta
    var.jit_seq.append((event.time, event.function, event.parent))

def jit_end():
    if var.jit_accum:
        var.jit.append(var.jit_accum)
        var.jit_accum = 0
        var.jit_end.append(var.last.time)
        var.jit_seqs.append(tuple(var.jit_seq))
        del var.jit_seq[:]
    var.started = False

def app(delta):
    if (delta):
        var.app.append(delta)

# ------------------------------------------------------------------------------

LIMIT_PID   = 0
LIMIT_BEGIN = 'BEGIN'
LIMIT_END   = 'END'

def process_start():
    process_event(Event(
                LIMIT_BEGIN, 
                LIMIT_PID, 
                var.start, 
                LIMIT_BEGIN, 
                LIMIT_BEGIN))

def process_end():
    process_event(Event(
                LIMIT_END,   
                LIMIT_PID, 
                var.end,   
                LIMIT_END,   
                LIMIT_END))

    if var.started:
        jit_end()

# ------------------------------------------------------------------------------

def time_statistics(list):
    print 'Total: %s' % (time_i2s(sum(list)))
    print 'Number of events: %d' % (len(list))
    print 'Min: %s' % (time_i2s(min(list)))
    print 'Max: %s' % (time_i2s(max(list)))
    print 'Avg: %s' % (time_i2s(np.mean(list)))
    print 'Std: %s' % (time_i2s(np.std(list)))
    print
    print 'Histogram:'
    time_histo(list)

def time_i2s(time):
    ts = '%010d' % (time)
    s  = int(ts[:-9])
    us = int(ts[-9:]) / 1000
    return '%d.%06d' % (s, us)

def time_histo(list):
    hist, bins = np.histogram(list)
    for i in range(len(hist)):
        print '  %s %s: %s' % (time_i2s(bins[i]), time_i2s(bins[i+1]), hist[i])

# ------------------------------------------------------------------------------

HLINE_A = 100 * '='
HLINE_B =  90 * '-'
HLINE_C = 135 * '-'

def show_app():
    print
    print 'APPLICATION STATISTICS'
    print HLINE_A
    print

    print 'First (forced): %s' % (time_i2s(var.app[0]))
    if len(var.app) > 1:
        time_statistics(var.app[1:])

def show_jit():
    print
    print 'JITTER STATISTICS'
    print HLINE_A
    print

    time_statistics(var.jit)

    if var.top or var.funstat:
        show_jit_top()
    if var.funstat:
        show_jit_funstat()
    if var.rootstat:
        show_jit_rootstat()
    print

def show_jit_top():
    print
    print 'TOP JITTER EVENTS'
    print HLINE_A
    print

    count = 0
    for i in reversed(np.argsort(var.jit)):
        count += 1
        if count > var.top:
            if not var.funstat:
                break
            else: 
                var.quiet      = True
                var.superquiet = True

        if var.detail and not var.superquiet:
            print HLINE_B
        if not var.superquiet:
            print '#%-6d EVENT: duration = %s | start = %s | end = %s' % (\
                    count,
                    time_i2s(var.jit[i]),
                    time_i2s(var.jit_start[i]),
                    time_i2s(var.jit_end[i]))
        if var.detail and not var.superquiet:
            print HLINE_B

        if var.detail == DETAIL_PLAIN:
            show_jit_detail_plain(i, count)
        if var.detail == DETAIL_GRAPH or var.funstat:
            show_jit_detail_graph(i, count)

def show_jit_detail_plain(i, count):
    functions = []
    for time, function, parent in var.jit_seqs[i]:
        if parent not in functions:
            functions.append(parent)
        if function not in functions:
            functions.append(function)

    if not var.superquiet:
        print ', '.join(functions)
        print

TAB = 4 * ' '

def show_jit_detail_graph(i, count):
    level  = 0
    stack  = []
    fu_evs = {}

    # Print call graph and gather function statistics
    for time, function, parent in var.jit_seqs[i]:
        while stack and parent != stack[-1][1]:
            pop_time, pop_function = stack.pop()
            dts = ''
            if pop_time:
                dt  = time - pop_time
                dts = time_i2s(dt)
                if pop_function not in fu_evs:
                    fu_evs[pop_function] = []
                if pop_function not in var.funstats:
                    var.funstats[pop_function] = []
                fu_evs[pop_function].append(dt)
                var.funstats[pop_function].append(dt)
                # root statistics
                if len(stack) == 1:
                    if pop_function not in var.rootstats:
                        var.rootstats[pop_function] = []
                    var.rootstats[pop_function].append(dt)
            level -= 1
            if not var.quiet:
                print '%20s | %20s | %s} # %s'    % (\
                        time_i2s(time),
                        dts,
                        level * TAB,
                        pop_function)
        if stack:
            stack.append((time, function))
            if not var.quiet:
                print '%20s | %20s | %s%s {' % (\
                        time_i2s(time),
                        '',
                        level * TAB,
                        function)
            level += 1
        else:
            stack.append((None, parent))
            if not var.quiet:
                print '%20s | %20s | %s%s {' % (\
                        '-',
                        '',
                        level * TAB,
                        parent)
            level += 1
            stack.append((time, function))
            if not var.quiet:
                print '%20s | %20s | %s%s {' % (\
                        time_i2s(time),
                        '',
                        level * TAB,
                        function)
            level += 1

    if var.superquiet or var.detail != DETAIL_GRAPH:
        return

    # Prepare sorting
    sorter_evs = []
    sorter_sum = []
    zero_list  = []
    for ev in fu_evs:
        tsum = sum(fu_evs[ev])
        if tsum > 0:
            sorter_evs.append(ev)
            sorter_sum.append(tsum)
        else:
            zero_list.append(ev)

    # Print functions
    print
    print '%-40s %10s %20s %20s %20s %20s' % (\
            '#%-6d FUNCTIONS' % (count),
            'count',
            'sum',
            'min',
            'max',
            'avg')
    print HLINE_C

    for i in reversed(np.argsort(sorter_sum)):
        ev = sorter_evs[i]
        print '%-40s %10d %20s %20s %20s %20s' % (\
                ev,
                len(fu_evs[ev]),
                time_i2s(sum(fu_evs[ev])),
                time_i2s(min(fu_evs[ev])),
                time_i2s(max(fu_evs[ev])),
                time_i2s(np.mean(fu_evs[ev])))
    print
    print 'zero list:', ','.join(zero_list)
    print

def show_jit_funstat():
    print
    print 'GLOBAL FUNCTION STATISTICS'
    print HLINE_A
    print

    sorter        = {}
    zero_sum_list = []
    zero_avg_list = []
    for ev in var.funstats:
        tsum = sum(var.funstats[ev])
        tavg = np.mean(var.funstats[ev])
        if   tsum < 1000:
            zero_sum_list.append(ev)
        elif tavg < 1000:
            zero_avg_list.append(ev)
        else:
            sorter[ev] = tsum

    print '%-40s %10s %20s %20s %20s %20s' % (\
            'function',
            'count',
            'sum',
            'min',
            'max',
            'avg')
    print HLINE_C

    for ev in sorted(sorter, key=sorter.get, reverse=True):
        print '%-40s %10d %20s %20s %20s %20s' % (\
                ev,
                len(var.funstats[ev]),
                time_i2s(sum(var.funstats[ev])),
                time_i2s(min(var.funstats[ev])),
                time_i2s(max(var.funstats[ev])),
                time_i2s(np.mean(var.funstats[ev])))
    print
    print 'zero sum list:', ','.join(zero_sum_list)
    print
    print 'zero avg list:', ','.join(zero_avg_list)
    print

def show_jit_rootstat():
    print                                       
    print 'GLOBAL ROOT FUNCTION STATISTICS'          
    print HLINE_A                               
    print                                           

    sorter = {}
    for ev in var.rootstats:
        tsum = sum(var.rootstats[ev])
        tavg = np.mean(var.rootstats[ev])
        sorter[ev] = tsum

    print '%-40s %10s %20s %20s %20s %20s' % (\
            'function',
            'count',
            'sum',
            'min',
            'max',
            'avg')
    print HLINE_C

    for ev in sorted(sorter, key=sorter.get, reverse=True):
        print '%-40s %10d %20s %20s %20s %20s' % (\
                ev,
                len(var.rootstats[ev]),
                time_i2s(sum(var.rootstats[ev])),
                time_i2s(min(var.rootstats[ev])),
                time_i2s(max(var.rootstats[ev])),
                time_i2s(np.mean(var.rootstats[ev])))


# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1:])

