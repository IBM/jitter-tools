libtmrand.a : tmrand.o tinymt64.o
	ar rv libtmrand.a tmrand.o tinymt64.o

clean :
	rm -f libtmrand.a tmrand.o tinymt64.o main.x

main.x : main.c libtmrand.a
	gcc -g -O2 -std=c99 main.c libtmrand.a -o main.x

tmrand.o : tmrand.c tinymt64.h
	gcc -c -Ofast -mcpu=power8 -std=c99 tmrand.c
 
tinymt64.o : tinymt64.c tinymt64.h
	gcc -c -Ofast -mcpu=power8 -std=c99 tinymt64.c
