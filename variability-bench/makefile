RAND = tmrand/libtmrand.a

variability : variability.c lutexp.o
	mpicc -g -O2 variability.c -o variability-bench lutexp.o $(RAND) -lm

lutexp.o : lutexp.c
	xlc_r -c -g -O3 -mcpu=power9 lutexp.c

clean :
	rm -f variability-bench lutexp.o
