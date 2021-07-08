RAND = tmrand/libtmrand.a
BIND = bindthreads.o -lnuma

osnoise : osnoise.c lutexp.o bindthreads.o
	mpicc -g -O2 osnoise.c -o osnoise lutexp.o $(RAND) $(BIND) -lm

lutexp.o : lutexp.c
	gcc -c -g -O3 -mcpu=power8 lutexp.c

bindthreads.o : bindthreads.c
	mpicc -c -g -O2 bindthreads.c

clean :
	rm -f osnoise lutexp.o
