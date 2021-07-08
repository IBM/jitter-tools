all:    no_aslr

no_aslr: no_aslr.o
	gcc -o $@ $^

no_aslr.o: no_aslr.c
	gcc -c $^

clean:
	rm -f no_aslr *.o
