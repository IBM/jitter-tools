RAND = ../tmrand/libtmrand.a

osnoise : libtmrand osnoise.o
	nvcc -ccbin mpicc -g -O2 osnoise.o -o osnoise  $(RAND) -lm  -lcudart_static -lculibos -lnvidia-ml -lrt -ldl

libtmrand:
	make -C ../tmrand

osnoise.o : osnoise.cu
	nvcc -ccbin mpicc -c -O3 -I. -gencode arch=compute_70,code=sm_70 osnoise.cu

clean :
	rm -f osnoise lutexp.o
