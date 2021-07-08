#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <cstring>
#include <fstream>
#include <time.h>

#define THOUSAND        1000
#define MILLION         1000000
#define BILLION         1000000000

#define USEC_PER_SEC    MILLION  // one sec  = million  usec
#define NSEC_PER_SEC    BILLION  // one sec  = billion  nsec
#define NSEC_PER_MSEC   MILLION  // one msec = million  nsec
#define NSEC_PER_USEC   THOUSAND // one usec = thousand nsec

class stat_t;

struct config_t {
  unsigned int conf_cpus; // configured cpus
  unsigned int conf_sockets;
  unsigned int smt_mode;

  bool verbose;

  struct timespec poll_delay;
  std::string root_cpuset;
  std::string root_cpuacct;
  std::ofstream log;

  stat_t *cpu_stats;
};

#endif /* #ifndef _CONFIG_H_ */
