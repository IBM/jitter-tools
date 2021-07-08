/* Copyright IBM corp.
 * author: Bryan Rosenburg
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/personality.h>

int main(int argc, char *argv[], char *envp[])
{
    if (argc < 2) {
	printf("Usage: %s <program> [ <arg1> <arg2> ... ]\n", argv[0]);
	exit(-1);
    }
    int persona;
    persona = personality(ADDR_NO_RANDOMIZE);
    //printf("old persona: 0x%08x\n", persona);
    int rc = execve(argv[1], &argv[1], envp);
    fprintf(stderr, "execve(\"%s\", ...) returned %d, errno %d (%s)\n",
	    argv[1], rc, errno, strerror(errno));
    return 0;
}
