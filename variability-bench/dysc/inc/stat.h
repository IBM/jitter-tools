#ifndef _STAT_H_
#define _STAT_H_

#include <stdio.h>
#include <sched.h>
#include <vector>
#include <list>
#include "sched.h"

struct cpu_util_t {
  unsigned long long tdiff;
  unsigned int cpu;
};

static bool cmp_sort(const cpu_util_t& first, const cpu_util_t& second)
{
  return (first.tdiff < second.tdiff);
}

class stat_t {
private:
  FILE *fd;
  std::vector< unsigned long long > lasttimes;

public:
  explicit stat_t(unsigned cpuc) :
    fd(fopen("/proc/stat", "r")),
    lasttimes(cpuc)
  {
    if (!fd) {
      std::cerr << "ERR: couldn't open '/proc/stat'" << std::endl;
      return;
    }

    fseek(fd, 0, SEEK_SET);

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    // first line is a don't care -> summary
    getline(&line, &len, fd);

    unsigned t_cpu = 0;
    while ((read = getline(&line, &len, fd)) != -1
           && t_cpu < lasttimes.size()) {
      char *t_line = line;
      t_line += 3; // cpu
      t_line += t_cpu / 10; // NUM
      ++t_line; // empty space
      unsigned long long utime = strtoull(t_line, &t_line, 10);
      ++t_line;
      unsigned long long ntime = strtoull(t_line, &t_line, 10);
      ++t_line;
      unsigned long long stime = strtoull(t_line, &t_line, 10);
      this->lasttimes[ t_cpu ] = utime + ntime + stime;
      ++t_cpu;
    }
    free(line);
  }
  ~stat_t(void)
  {
    if (this->fd)
      fclose(this->fd);
  }

  cpuset_t lowest_utilization(cpuset_t& avail_cpu_mask,
                              unsigned cpu_req_amount)
  {
    cpuset_t mask(lasttimes.size());
    if (!this->fd)
      return mask;

    fseek(this->fd, 0, SEEK_SET);

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    // first line is a don't care -> summary
    getline(&line, &len, this->fd);

    std::list< cpu_util_t > g;
    unsigned t_cpu = 0;
    while ((read = getline(&line, &len, this->fd)) != -1
           && t_cpu < lasttimes.size()) {
      char *t_line = line;
      t_line += 3; // cpu
      t_line += t_cpu / 10; // NUM
      ++t_line; // empty space
      unsigned long long utime = strtoull(t_line, &t_line, 10);
      ++t_line;
      unsigned long long ntime = strtoull(t_line, &t_line, 10);
      ++t_line;
      unsigned long long stime = strtoull(t_line, &t_line, 10);

      unsigned long long cmp_diff = (utime + ntime + stime) - this->lasttimes[ t_cpu ];
      this->lasttimes[ t_cpu ] = utime + ntime + stime;

      if (avail_cpu_mask.isset_bit(t_cpu))
        g.push_back({ cmp_diff, t_cpu });

      ++t_cpu;
    }
    free(line);

    g.sort(cmp_sort);

    for (unsigned c = 0; c < cpu_req_amount
         && !g.empty(); ++c) {
      mask.set_bit(g.front().cpu);
      g.pop_front();
    }

    return mask;
  }
};

#endif /* #ifndef _STAT_H_ */
