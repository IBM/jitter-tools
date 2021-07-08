all:    jwind_cgroup

jwind_cgroup: jwind_cgroup.o
	gcc -o $@ $^ -lrt

jwind_cgroup.o: jwind_cgroup.c
	gcc -O2 -c $^
