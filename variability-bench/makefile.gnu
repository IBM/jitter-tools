RAND  = tmrand/libtmrand.a
CC    = /usr/bin/gcc
MPICC = OMPI_CC=$(CC) /opt/ibm/spectrum_mpi/bin/mpicc
MAKE  = /usr/bin/make
RM    = /usr/bin/rm -f

variability : variability.c lutexp.o build-tmrand
	$(MPICC) -g -O2 variability.c -o variability-bench lutexp.o $(RAND) -lm

lutexp.o : lutexp.c
	$(CC) -c -g -O3 lutexp.c

build-tmrand:
	$(MAKE) -C tmrand/

clean :
	$(RM) variability-bench lutexp.o
	$(MAKE) -C tmrand/ clean
