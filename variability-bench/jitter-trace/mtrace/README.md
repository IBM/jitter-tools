# mtrace: Communication pattern tracer

This library "overrides" the MPI blocking calls and traces them. It's faster
than general solutions as ltrace or mpitrace but doesn't provide as much
information, used by default writes some equivalent phase lengths to use as
`-c` parameter in jitter-bench.

## Build

build it with:

`make`

or, if you also want the dump of all the traced events:

`make full_log`

## Use

Make it so that this library is loaded before the standard mpi one:

`export LD_PRELOAD=$PWD/mtrace.so`

## Output

a file named 'eq.{nranks}.log' with:

- Equavilent phase length for non-propagating-noise communication events and
  number of iterations;
- Equivalent phase length for propagating-noise communication events and number
  of iterations;
- Equivalent phase length for other communication events and number of
  iterations;
- Total time spent by the application.
