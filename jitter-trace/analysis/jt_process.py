#!/usr/bin/python -u

import numpy as np
import sys
import time

# ==============================================================================
# global variables
# ==============================================================================

class var: pass

var.tml = None
var.app = None

var.sft = False
var.rcu = False
var.tim = False
var.ppc = False

var.beg = '0.0'
var.end = 'inf'
var.was = False

# ==============================================================================
# main function
# ==============================================================================

def main(args):

    parse_args(args)
    d = process_tml(var.tml, var.app)
    show(d, var.app)

# ==============================================================================
# parse command line arguments
# ==============================================================================

def parse_args(args):

    i = 0
    while i < len(args):
        if args[i] == '-h':
            usage()
        elif args[i] == '-f':
            var.tml = args[i+1]
            i += 2
        elif args[i] == '-a':
            var.app = args[i+1]
            i += 2
        elif args[i] == '-s':
            var.sft = True
            i += 1
        elif args[i] == '-r':
            var.rcu = True
            i += 1
        elif args[i] == '-t':
            var.tim = True
            i += 1
        elif args[i] == '-p':
            var.ppc = True
            i += 1
        elif args[i] == '-b':
            var.beg = args[i+1]
            i += 2
        elif args[i] == '-e':
            var.end = args[i+1]
            i += 2
        else:
            usage()

    if any(i is None for i in [var.tml, var.app]):
        print 'error: missing arguments'
        usage()

    handler_list = sorted(handlers.keys())
    if not var.sft:
        for handler in handler_list:
            if handler.startswith('irq:softirq_'):
                del handlers[handler]
    if not var.rcu:
        for handler in handler_list:
            if handler.startswith('rcu:'):
                del handlers[handler]
    if not var.tim:
        for handler in handler_list:
            if handler.startswith('timer:'):
                del handlers[handler]
    if not var.ppc:
        for handler in handler_list:
            if handler.startswith('powerpc:'):
                del handlers[handler]

# ==============================================================================
# print usage information and exit
# ==============================================================================

def usage():

    print 'usage: %s -f FILE -a APP [-s] [-r] [-t] [-b BEGIN] [-e END]' % (sys.argv[0])
    print 'required:'
    print '  -f FILE    timeline (.tml) file'
    print '  -a APP     application (must match process name)'
    print 'optional:'
    print '  -s         enable softirq analysis (default: False)'
    print '  -r         enable rcu analysis (default: False)'
    print '  -t         enable timer analysis (default: False)'
    print '  -b         timestamp to begin analysis'
    print '  -e         timestamp to end analysis'
    print
    exit(1)

# ==============================================================================
# data structure
# ==============================================================================

class Data:

    def __init__(self):

        self.cpus    = [] # list of cpus
        self.active  = {} # active process          [cpu]
        self.stack   = {} # secondary stack         [cpu]
        self.watch   = {} # watch map               [cpu]
        self.app_rt  = {} # app runtime             [cpu]
        self.app_jit = {} # app jitter              [cpu]
        self.off_jit = {} # off jitter              [cpu][off]
        self.rt      = {} # app sched runtime       [cpu]
        self.rtj     = {} # off sched runtime       [off]

# ==============================================================================
# group process names
# ==============================================================================

proc_aliases = [
    'kworker',
    'ksoftirqd',
    'swapper',
    'watchdog',
    'migration'
]

def proc_group(proc):

    for alias in proc_aliases:
        if ('%s/' % (alias)) in proc:
            return alias

    return proc

# ==============================================================================
# determine if time is inside enabled interval
# ==============================================================================

def inside(t):

    ts, tns = time_convert(t)
    bs, bns = time_convert(var.beg)

    if (ts > bs) or (ts == bs and tns >= bns):
        if var.end == 'inf':
            return True
        else:
            es, ens = time_convert(var.end)
            if (ts < es) or (ts == es and tns <= ens):
                return True

    return False

# ==============================================================================
# convert timestamp
# ==============================================================================

def time_convert(time):

    sec, nsec = time.split('.')
    sec  = int(sec)
    nsec = int(nsec.ljust(9, '0'))

    return sec, nsec

# ==============================================================================
# calculate time difference
# ==============================================================================

US_PER_S = 1000000
NS_PER_S = 1000000000

def time_diff(start, end):

    ts_sec, ts_nsec = time_convert(start)
    te_sec, te_nsec = time_convert(end)

    dt_sec  = te_sec  - ts_sec
    dt_nsec = te_nsec - ts_nsec

    if dt_nsec < 0:
        dt_sec  -= 1
        dt_nsec += NS_PER_S

    try:
        diff = float('%d.%09d' % (dt_sec, dt_nsec))
    except ValueError, e:
        print e
        print '[debug] time_diff: start:', start
        print '[debug] time_diff: end:  ', end
        exit(1)
    return diff

# ==============================================================================
# create histogram
# ==============================================================================

def histogram(data):

    if len(data) > 1:
        return np.histogram(data)
    else:
        return np.histogram(data, range=(0.0, data[0]))

# ==============================================================================
# process timeline file
# ==============================================================================

STEP_TOTAL = 1000.0
SPACE      = 20 * ' '

def process_tml(tml, app):

    # ----------------------
    # prepare data structure
    # ----------------------

    d = Data()

    # ---------
    # open file
    # ---------

    try:
        f = open(tml, 'r')
    except IOError:
        print 'error: could not open "%s"' % (tml)
        exit(1)

    # ----------------------------
    # count lines to show progress
    # ----------------------------

    print >> sys.stderr,\
        'counting lines...',
    start = time.time()
    lines = 0
    for l in f:
        lines += 1
    end   = time.time()
    print >> sys.stderr,\
        '%d (%d seconds)' % (lines, end - start)

    # ---------------------------------------
    # calculate progress step and start timer
    # ---------------------------------------

    step_line  = lines / int(STEP_TOTAL)
    start      = time.time()

    # -----------------------------
    # reset reader and process file
    # -----------------------------

    f.seek(0)
    line = 0
    step = 0

    for l in f:
        if process_line(l, app, d):
            break
        line += 1

        # -------------
        # show progress
        # -------------

        if line % step_line == 0:
            step += 1
            print >> sys.stderr,\
                'processing %.1f%%%s\r' % (100.0 * step / STEP_TOTAL, SPACE),

    # ---------------------
    # stop timer and return
    # ---------------------

    end = time.time()
    print >> sys.stderr,\
        'processing done (%d seconds)' % (end - start)

    return d

# ==============================================================================
# process line
# ==============================================================================

def process_line(l, app, d):

    # ---------------------------------------------
    # parse and add cpu to all pertinent structures
    # ---------------------------------------------

    split = l.strip().split('[')
    cpu   = int(split[1].split(']')[0])

    # -----------------------------------
    # add cpu to all pertinent structures
    # -----------------------------------

    if cpu not in d.cpus:
        d.cpus.append(cpu)
        d.active[cpu]  = None
        d.stack[cpu]   = []
        d.watch[cpu]   = False
        d.app_rt[cpu]  = []
        d.app_jit[cpu] = []
        d.off_jit[cpu] = {}
        d.rt[cpu]      = 0

    # ------------------------
    # process supported events
    # ------------------------

    sub2  = filter(None, '['.join(split[1:]).split(' '))
    type  = sub2[2][:-1]

    if type in handlers:
        sub1  = filter(None, split[0].split(' '))
        pid   = sub1[-1]
        proc  = proc_group(' '.join(sub1[:-1]))
        time  = sub2[1][:-1]
        args  = ' '.join(sub2[3:])

        # -----------------------------------------
        # process event if within analysis interval
        # -----------------------------------------

        if inside(time):
            if not var.was:
                print >> sys.stderr, 'inside analysis interval: ', var.beg, time
                for cpu in d.cpus:
                    do_switch(None, var.app, d, var.app, '0', cpu, var.beg, None,
                              'swapper', '0', 'S', 'swapper:0',
                              var.app, '0', '%s:0' % (var.app))
                var.was = True
            handlers[type](l, app, d, proc, pid, cpu, time, args)
        else:
            if var.was:
                print >> sys.stderr, 'outside analysis interval:', var.end, time
                for cpu in d.cpus:
                    do_switch(None, var.app, d, 'swapper', '0', cpu, var.end, None,
                              var.app, '0', 'S', '%s:0' % (var.app),
                              'swapper', '0', 'swapper:0')
                var.was = False
                return True # stop

# ==============================================================================
# handle runtime
# ==============================================================================

def handle_runtime(l, app, d, proc, pid, cpu, time, args):

    runtime = int(args.split('runtime=')[1].split(' ')[0])

    if proc == app:
        d.rt[cpu] += runtime
    else:
        if proc not in d.rtj:
            d.rtj[proc] = 0.0
        d.rtj[proc] += runtime

# ==============================================================================
# handle switch
# ==============================================================================

def handle_switch(l, app, d, proc, pid, cpu, time, args):

    split = args.split(' ==> ')

    swo       = split[0].split(':')
    swo_proc  = proc_group(':'.join(swo[:-1]))
    swo_pid   = swo[-1].split(' ')[0]
    swo_state = swo[-1].split(' ')[-1]
    swo_id    = '%s:%s' % (swo_proc, swo_pid)

    swi       = split[1].split(':')
    swi_proc  = proc_group(':'.join(swi[:-1]))
    swi_pid   = swi[-1].split(' ')[0]
    swi_id    = '%s:%s' % (swi_proc, swi_pid)

    do_switch(l, app, d, proc, pid, cpu, time, args,
              swo_proc, swo_pid, swo_state, swo_id,
              swi_proc, swi_pid, swi_id)

    # reset stack
    del d.stack[cpu][:]

# ==============================================================================
# implement switch
# ==============================================================================

def do_switch(l, app, d, proc, pid, cpu, time, args,
              swo_proc, swo_pid, swo_state, swo_id,
              swi_proc, swi_pid, swi_id):

    # handle application swout
    # --------------------------------------------------------------------------

    if swo_proc == app:

        # try to compute runtime (if active process matches)
        if d.active[cpu]:
            active_proc, active_pid, active_id, active_time = d.active[cpu]
            if swo_id == active_id:
                d.app_rt[cpu].append(time_diff(active_time, time))

        # update watch map
        if swo_state in ['R', 'D']:
            d.watch[cpu] = (swo_id, time)

    # handle offender swout
    # --------------------------------------------------------------------------

    else:

        # try to compute offender perspective jitter
        if d.active[cpu] and d.watch[cpu]:
            active_proc, active_pid, active_id, active_time = d.active[cpu]
            if swo_id == active_id:
                if swo_proc not in d.off_jit[cpu]:
                    d.off_jit[cpu][swo_proc] = []
                d.off_jit[cpu][swo_proc].append(time_diff(active_time, time))


    # handle application swin
    # --------------------------------------------------------------------------

    if swi_proc == app:

        # set active process
        d.active[cpu] = (swi_proc, swi_pid, swi_id, time)

        # try to compute application perspective jitter, then reset watch map
        if d.watch[cpu]:
            watch_id, watch_time = d.watch[cpu]
            if swi_id == watch_id:
                d.app_jit[cpu].append(time_diff(watch_time, time))
        d.watch[cpu] = None

    # handle offender swin
    # --------------------------------------------------------------------------

    else:

        # set active process
        d.active[cpu] = (swi_proc, swi_pid, swi_id, time)

# ==============================================================================
# handle stack event start
# ==============================================================================

def handle_stack_start(l, app, d, proc, pid, cpu, time, args, name, id):

    # treat events that are processed on behalf of other processes as stacked
    # events; stacked events, such as interrupts and workqueues, temporarily
    # steal time from the underlying process; this implementation handles
    # stacked events as a regular switch with specific values

    if not d.active[cpu]:
        return

    swo_proc, swo_pid, swo_id, _ = d.active[cpu]
    d.stack[cpu].append((swo_proc, swo_pid, swo_id))

    do_switch(l, app, d, proc, pid, cpu, time, args,
              swo_proc, swo_pid, 'R', swo_id,
              name, pid, id)

# ==============================================================================
# handle stack event end
# ==============================================================================

def handle_stack_end(l, app, d, proc, pid, cpu, time, args, id):

    if not d.active[cpu] or not d.stack[cpu]:
        return

    swo_proc, swo_pid, swo_id, _ = d.active[cpu]
    swi_proc, swi_pid, swi_id    = d.stack[cpu].pop()

    do_switch(l, app, d, proc, pid, cpu, time, args,
              swo_proc, swo_pid, 'S', swo_id,
              swi_proc, swi_pid, swi_id)

# ==============================================================================
# stack handlers
# ==============================================================================

def handle_wq_start(l, app, d, proc, pid, cpu, time, args):

    wq_str = args.split('struct ')[1].split(':')[0]
    wq_fun = args.split('function ')[1]
    name   = '%s#wq:%s' % (proc, wq_fun)
    id     = '%s:%s#wq%s' % (proc, pid, wq_str)

    handle_stack_start(l, app, d, proc, pid, cpu, time, args, name, id)

# ==============================================================================

def handle_wq_end(l, app, d, proc, pid, cpu, time, args):

    wq_str = args.split('struct ')[1].split(':')[0]
    id     = '%s:%s#wq%s' % (proc, pid, wq_str)

    handle_stack_end(l, app, d, proc, pid, cpu, time, args, id)

# ==============================================================================

def handle_irq_start(l, app, d, proc, pid, cpu, time, args):

    irq_code = args.split('irq=')[1].split(' name=')[0]
    irq_name = args.split('name=')[1]
    name     = '%s#irq:%s=%s' % (proc, irq_code, irq_name)
    id       = '%s:%s#irq:%s' % (proc, pid, irq_code)

    handle_stack_start(l, app, d, proc, pid, cpu, time, args, name, id)

# ==============================================================================

def handle_irq_end(l, app, d, proc, pid, cpu, time, args):

    irq_code = args.split('irq=')[1].split(' ret=')[0]
    id       = '%s:%s#irq:%s' % (proc, pid, irq_code)

    handle_stack_end(l, app, d, proc, pid, cpu, time, args, id)

# ==============================================================================

def handle_softirq_start(l, app, d, proc, pid, cpu, time, args):

    softirq_vec = args.split('vec=')[1].split(' ')[0]
    softirq_act = args.split('action=')[1].split(']')[0]

    name = '%s#softirq:%s=%s' % (proc, softirq_vec, softirq_act)
    id   = '%s:%s#softirq:%s' % (proc, pid, softirq_vec)

    handle_stack_start(l, app, d, proc, pid, cpu, time, args, name, id)

# ==============================================================================

def handle_softirq_end(l, app, d, proc, pid, cpu, time, args):

    softirq_vec = args.split('vec=')[1].split(' ')[0]
    softirq_act = args.split('action=')[1].split(']')[0]
    id          = '%s:%s#softirq:%s' % (proc, pid, softirq_vec)

    handle_stack_end(l, app, d, proc, pid, cpu, time, args, id)

# ==============================================================================

def handle_rcu(l, app, d, proc, pid, cpu, time, args):

    name = '%s#rcu'    % (proc)
    id   = '%s:%s#rcu' % (proc, pid)

    if 'Start' in args:
        handle_stack_start(l, app, d, proc, pid, cpu, time, args, name, id)
    elif 'End' in args:
        handle_stack_end(l, app, d, proc, pid, cpu, time, args, id)
    else:
        print >> sys.stderr,\
            'unhandled rcu event: time=%s cpu=%03d, proc=%s' % (time, cpu, proc)

# ==============================================================================

def handle_timer_expire_start(l, app, d, proc, pid, cpu, time, args):

    timer = args.split('timer=')[1].split(' ')[0]
    name  = '%s#timer'       % (proc)
    id    = '%s:%s#timer:%s' % (proc, pid, timer)

    handle_stack_start(l, app, d, proc, pid, cpu, time, args, name, id)

# ==============================================================================

def handle_timer_expire_end(l, app, d, proc, pid, cpu, time, args):

    timer = args.split('timer=')[1].split(' ')[0]
    id    = '%s:%s#timer:%s' % (proc, pid, timer)

    handle_stack_end(l, app, d, proc, pid, cpu, time, args, id)

# ==============================================================================

def handle_hrtimer_expire_start(l, app, d, proc, pid, cpu, time, args):

    hrtimer = args.split('hrtimer=')[1].split(' ')[0]
    name    = '%s#hrtimer'       % (proc)
    id      = '%s:%s#hrtimer:%s' % (proc, pid, hrtimer)

    handle_stack_start(l, app, d, proc, pid, cpu, time, args, name, id)

# ==============================================================================

def handle_hrtimer_expire_end(l, app, d, proc, pid, cpu, time, args):

    hrtimer = args.split('hrtimer=')[1].split(' ')[0]
    id      = '%s:%s#hrtimer:%s' % (proc, pid, hrtimer)

    handle_stack_end(l, app, d, proc, pid, cpu, time, args, id)

# ==============================================================================

def handle_ppc_timer_start(l, app, d, proc, pid, cpu, time, args):

    addr = args.split('pt_regs=')[1].split(' ')[0]
    name = '%s#ppctimer'       % (proc)
    id   = '%s:%s#ppctimer:%s' % (proc, pid, addr)

    handle_stack_start(l, app, d, proc, pid, cpu, time, args, name, id)

# ==============================================================================

def handle_ppc_timer_end(l, app, d, proc, pid, cpu, time, args):

    addr = args.split('pt_regs=')[1].split(' ')[0]
    id   = '%s:%s#ppctimer:%s' % (proc, pid, addr)

    handle_stack_end(l, app, d, proc, pid, cpu, time, args, id)

# ==============================================================================
# handlers
# ==============================================================================

handlers = {
    'sched:sched_switch'                : handle_switch,
    'sched:sched_stat_runtime'          : handle_runtime,
    'workqueue:workqueue_execute_start' : handle_wq_start,
    'workqueue:workqueue_execute_end'   : handle_wq_end,
    'irq:irq_handler_entry'             : handle_irq_start,
    'irq:irq_handler_exit'              : handle_irq_end,
    'irq:softirq_entry'                 : handle_softirq_start,
    'irq:softirq_exit'                  : handle_softirq_end,
    'rcu:rcu_utilization'               : handle_rcu,
    'timer:timer_expire_entry'          : handle_timer_expire_start,
    'timer:timer_expire_exit'           : handle_timer_expire_end,
    'timer:hrtimer_expire_entry'        : handle_hrtimer_expire_start,
    'timer:hrtimer_expire_exit'         : handle_hrtimer_expire_end,
    'powerpc:timer_interrupt_entry'     : handle_ppc_timer_start,
    'powerpc:timer_interrupt_exit'      : handle_ppc_timer_end,
}

# ==============================================================================
# show results
# ==============================================================================

HLINE = 180 * '='

def show(d, app):

    show_app_rt_summ(d, app)
    show_app_rt_histo(d, app)
    show_app_jit_summ(d, app)
    show_app_jit_histo(d, app)
    show_off_jit_summ(d, app)
    show_off_jit_histo(d, app)
    show_rt(d, app)
    show_rtj(d, app)

# ==============================================================================
# show application runtime summary per cpu
# ==============================================================================

def show_app_rt_summ(d, app):

    print
    print 'application runtime (summary per cpu)'
    print HLINE
    print

    print '%12s %12s %12s %12s %12s %12s %12s %12s' % (\
            'cpu',
            '#events',
            'min(s)',
            'max(s)',
            'avg(s)',
            'std(s)',
            'var(%)',
            'sum(s)')
    print 103 * '-'

    for cpu in sorted(d.cpus):

        if d.app_rt[cpu]:
            _len = len(d.app_rt[cpu])
            _min = min(d.app_rt[cpu])
            _max = max(d.app_rt[cpu])
            _avg = np.mean(d.app_rt[cpu])
            _std = np.std(d.app_rt[cpu])
            _var = 100.0 * _std / _avg
            _sum = sum(d.app_rt[cpu])

            print '%12s %12d %12.6f %12.6f %12.6f %12.6f %12.2f %12.6f' % (\
                    '%03d' % (cpu),
                    _len,
                    _min,
                    _max,
                    _avg,
                    _std,
                    _var,
                    _sum)

        else:
            print '%12s %12s %12s %12s %12s %12s %12s %12s' % (\
                    '%03d' % (cpu),
                    '-',
                    '-',
                    '-',
                    '-',
                    '-',
                    '-',
                    '-')

# ==============================================================================
# show application runtime histogram per cpu
# ==============================================================================

def show_app_rt_histo(d, app):

    print
    print 'application runtime (histogram per cpu)'
    print HLINE

    for cpu in sorted(d.cpus):
        if d.app_rt[cpu]:
            print
            print 'cpu %03d' % (cpu)
            print 38 * '-'

            hist, edges= histogram(d.app_rt[cpu])
            for i in range(len(hist)):
                print '%12.6f %12.6f %12d' % (\
                        edges[i],
                        edges[i+1],
                        hist[i])

# ==============================================================================
# show application perspective jitter summary per cpu
# ==============================================================================

def show_app_jit_summ(d, app):

    print
    print 'application perspective jitter (summary per cpu)'
    print HLINE
    print

    print '%12s %12s %12s %12s %12s %12s %12s %12s' % (\
            'cpu',
            '#events',
            'min(us)',
            'max(us)',
            'avg(us)',
            'std(us)',
            'var(%)',
            'sum(us)')
    print 103 * '-'

    for cpu in sorted(d.cpus):

        if d.app_jit[cpu]:
            _len = len(d.app_jit[cpu])
            _min = min(d.app_jit[cpu])      * US_PER_S
            _max = max(d.app_jit[cpu])      * US_PER_S
            _avg = np.mean(d.app_jit[cpu])  * US_PER_S
            _std = np.std(d.app_jit[cpu])   * US_PER_S
            _var = 100.0 * _std / _avg
            _sum = sum(d.app_jit[cpu])      * US_PER_S

            print '%12s %12d %12.3f %12.3f %12.3f %12.3f %12.2f %12.3f' % (\
                '%03d' % (cpu),
                _len,
                _min,
                _max,
                _avg,
                _std,
                _var,
                _sum)

        else:
            print '%12s %12s %12s %12s %12s %12s %12s %12s' % (\
                    '%03d' % (cpu),
                    '-',
                    '-',
                    '-',
                    '-',
                    '-',
                    '-',
                    '-')

# ==============================================================================
# show application perspective jitter histogram per cpu
# ==============================================================================

def show_app_jit_histo(d, app):

    print
    print 'application perspective jitter (histogram per cpu)'
    print HLINE

    for cpu in sorted(d.cpus):
        if d.app_jit[cpu]:
            print
            print 'cpu %03d' % (cpu)
            print 38 * '-'

            hist, edges = histogram(d.app_jit[cpu])
            for i in range(len(hist)):
                print '%12.6f %12.6f %12d' % (\
                        edges[i],
                        edges[i+1],
                        hist[i])

# ==============================================================================
# show offender perspective jitter summary per cpu
# ==============================================================================

def show_off_jit_summ(d, app):

    print
    print 'offender perspective jitter (summary per cpu)'
    print HLINE

    for cpu in sorted(d.cpus):
        if d.off_jit[cpu]:
            print
            print '%12s %50s %12s %12s %12s %12s %12s %12s %12s' % (\
                    '%03d' % cpu,
                    'offender',
                    '#events',
                    'min(us)',
                    'max(us)',
                    'avg(us)',
                    'std(us)',
                    'var(%)',
                    'sum(us)')
            print 154 * '-'

            for off in sorted(d.off_jit[cpu]):
                jitters = d.off_jit[cpu][off]
                _len = len(jitters)
                _min = min(jitters)         * US_PER_S
                _max = max(jitters)         * US_PER_S
                _avg = np.mean(jitters)     * US_PER_S
                _std = np.std(jitters)      * US_PER_S
                _var = 100.0 * _std / _avg
                _sum = sum(jitters)         * US_PER_S

                print '%12s %50s %12d %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f' % (\
                        '',
                        off,
                        _len,
                        _min,
                        _max,
                        _avg,
                        _std,
                        _var,
                        _sum)

# ==============================================================================
# show offender perspective jitter histogram per cpu
# ==============================================================================

def show_off_jit_histo(d, app):

    print
    print 'offender perspective jitter (histogram per cpu)'
    print HLINE

    for cpu in sorted(d.cpus):
        if d.off_jit[cpu]:
            print
            print 'cpu %03d' % (cpu)
            print 38 * '-'

            for off in sorted(d.off_jit[cpu]):
                print
                print off
                print 'occurrences=%d min=%.3fusec max=%.3fusec sum=%.3fusec' % (\
                        len(d.off_jit[cpu][off]),
                        min(d.off_jit[cpu][off]),
                        max(d.off_jit[cpu][off]),
                        sum(d.off_jit[cpu][off]))

                hist, edges = histogram(d.off_jit[cpu][off])
                for i in range(len(hist)):
                    print '%12.3f %12.3f %13d' % (\
                            edges[i]   * US_PER_S,
                            edges[i+1] * US_PER_S,
                            hist[i])

# ==============================================================================
# show application sched runtime
# ==============================================================================

def show_rt(d, app):

    print
    print 'application sched runtime (per cpu)'
    print HLINE
    print

    print '%12s %s' % ('cpu', 'runtime(s)')
    for cpu in sorted(d.cpus):
        print '%12s %12.9f' % ('%03d' % cpu, float(d.rt[cpu]) / float(NS_PER_S))

# ==============================================================================
# show offender sched runtime
# ==============================================================================

def show_rtj(d, app):

    print
    print 'offender sched runtime'
    print HLINE
    print

    print '%50s %12s' % (\
            'offender',
            'runtime(s)')
    print 63 * '-'

    for off in sorted(d.rtj):
        print '%50s %12.9f' % (\
                off,
                float(d.rtj[off]) / float(NS_PER_S))
    print '%50s %12.9f' % (\
            'TOTAL',
            sum(map(float, d.rtj.values())) / float(NS_PER_S))

# ==============================================================================
# call main
# ==============================================================================

if __name__ == '__main__':
    main(sys.argv[1:])

