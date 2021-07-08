# Compile for IBM POWER
CC = gcc
MPICC = OMPI_CC=$(CC) mpicc -std=c99 -Ofast -W -Wall
LUTEXP_CC = xlc_r -g -O3 -mcpu=power9
LUTEXP_CU_CC = nvcc -ccbin mpicc -O3 -gencode arch=compute_70,code=sm_70
export TMRAND_CC = gcc -std=c99 -c -Ofast -mcpu=power8

# Compile for Cray Intel
#MPICC = cc -O3
#LUTEXP_CC = cc -O3
#export TMRAND_CC = cc -c -O3



all:    jitter-bench readjbdump

jitter-bench: jitter-bench.c lutexp.o tmrand/libtmrand.a
	$(MPICC) $^ -lm -o $@

jitter-bench-cu: lutexp-cu.o jitter-bench.c lutexp.o tmrand/libtmrand.a
	$(MPICC) -c jitter-bench.c -o jitter-bench.o
	OMPI_CC=g++ mpic++ jitter-bench.o lutexp.o lutexp-cu.o tmrand/libtmrand.a -lm -L/usr/local/cuda-9.2/targets/ppc64le-linux/lib/ -lcudart -o $@

jitter-bench-cu-run: jitter-bench-cu
	sh jitter-bench-cu-run.sh

readjbdump: readjbdump.c
	$(CC) readjbdump.c -o readjbdump -std=gnu11 -Wall -Wextra

# ------------------
# Utils
# ------------------

lutexp.o: lutexp.c lutexp.h
	$(LUTEXP_CC) -c $< -o $@

lutexp-cu.o: lutexp.cu
	$(LUTEXP_CU_CC) -c $< -o $@

tmrand/libtmrand.a:
	make -C tmrand/

clean:
	rm -f *.o jitter-bench jitter-bench-cu
	make -C tmrand/ clean
