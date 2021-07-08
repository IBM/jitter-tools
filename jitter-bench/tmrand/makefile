libtmrand.a : tmrand.o tinymt64.o
	ar rv libtmrand.a tmrand.o tinymt64.o

clean :
	rm -f libtmrand.a tmrand.o tinymt64.o

main.x : main.c libtmrand.a
	$(TMRAND_CC) -std=c99 main.c libtmrand.a -o main.x


tmrand.o : tmrand.c tinymt64.h
	$(TMRAND_CC) tmrand.c
 
tinymt64.o : tinymt64.c tinymt64.h
	$(TMRAND_CC) tinymt64.c
