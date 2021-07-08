CC = mpicc -Wall -g -shared -O3 -fPIC -mcpu=power9

all: mtrace.so

full_log: CFLAGS += -Dfull_log	
full_log: mtrace

mtrace.so: mtrace.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f *.so
