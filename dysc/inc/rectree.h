#ifndef _RECTREE_H_
#define _RECTREE_H_

#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <ctime>
#include <dirent.h>
#include <fcntl.h>
#include <stdio_ext.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "sched.h"
#include "config.h"
#include "stat.h"
#include "policy.h"


static inline void time_get(
  struct timespec *t)
{
  clock_gettime(CLOCK_MONOTONIC, t);
}

static inline void time_diff(
  struct timespec *dt,
  struct timespec *t1,
  struct timespec *t0)
{
  dt->tv_sec = t1->tv_sec - t0->tv_sec;
  dt->tv_nsec = t1->tv_nsec - t0->tv_nsec;
  if (t1->tv_nsec < t0->tv_nsec) {
    dt->tv_sec -= 1;
    dt->tv_nsec += NSEC_PER_SEC;
  }
}


static const char * folder_iter(const char *path)
{
  const char *str_it = path;
  while (*str_it != '\0') {
    if (*str_it == '/'
        && !(str_it == path
             || *(str_it - 1) == '\\')) {
      return str_it;
    }
    ++str_it;
  }
  return nullptr;
}


class cgroup_leave;

class cgroup_abstract_t {
protected:
// --- global config
  config_t& cfg;

// --- paths and open decriptors
  std::string root_cpuset;
  int fd_cpuset;
  std::string root_cpuacct;
  int fd_cpuacct;

// --- parent and children
  std::unordered_map< std::string, cgroup_leave* > children;

// --- counters
  unsigned long long old_cpu_usage_ns;

// --- policies and own mask
  cpuset_t cur_mask;

  int check_fd_cpuacct(void)
  {
    if (__builtin_expect(this->fd_cpuacct < 0, 0))
      this->fd_cpuacct = open((this->root_cpuacct + "/cpuacct.usage").c_str(), O_RDONLY);
    return this->fd_cpuacct;
  }

  unsigned long long read_cpu_usage(int fd)
  {
    unsigned long long cpuacct_usage = 0;
    if (fd < 0)
      return cpuacct_usage;

    lseek(fd, 0, SEEK_SET);

    // in case of 64bit systems: ULLONG_MAX requires 21 digits
    // add +1 to avoid second while loop on 64bit systems
    // consequence: 1 read call!
    char buf[ 21 + 1 ];
    int len;
    do {
      len = read(fd, buf, sizeof(buf));
      for (char *c = buf; c < buf + len; ++c) {
        switch (*c) {
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
          cpuacct_usage *= 10;
          cpuacct_usage += (*c - '0');
          break;
        case '\0':
        case '\n':
          return cpuacct_usage;
        default: // shouldn't happen
          std::cerr << "ERR:" << __FUNCTION__ << ":" << __LINE__ << ": '" << *c << "' occured!" << std::endl;
          return ~0x0ull;
        }
      }
    } while (len == sizeof(buf));

    return cpuacct_usage;
  }

  void display(unsigned int ttotal_ns, unsigned int ton_ns, cpuset_maskset_t *own_policy)
  {
    printf("(%6.2f%%) %3d %c",
           ((double)ton_ns / (double)ttotal_ns) * 100.0,
           this->cur_mask.get_corec(),
           own_policy ? '*' : ' ');
    this->cur_mask.show_bitmask();
    printf(" (");
    if (own_policy) {
      own_policy->get_min().show_bitmask();
    }
    else {
      for (unsigned i = 0; i < (this->cfg.conf_cpus > 8 ? 8 : this->cfg.conf_cpus); ++i)
        printf("-");
      printf(",");
      for (unsigned i = 0; i < (this->cfg.conf_cpus > 8 ? 8 : this->cfg.conf_cpus); ++i)
        printf("-");
    }
    printf(") %s\n",
           this->root_cpuset.c_str() + this->cfg.root_cpuset.size());
  }

public:
  cgroup_abstract_t(config_t& cfg,
                    std::string& root_cpuset,
                    std::string& root_cpuacct) :
    cfg(cfg),
    root_cpuset(root_cpuset),
    fd_cpuset(open((root_cpuset + "/cpuset.cpus").c_str(), O_RDWR)),
    root_cpuacct(root_cpuacct),
    // reports the total CPU time (in nanoseconds) consumed by all tasks in this cgroup (including tasks lower in the hierarchy).
    fd_cpuacct(open((root_cpuacct + "/cpuacct.usage").c_str(), O_RDONLY)),
    old_cpu_usage_ns(this->read_cpu_usage(fd_cpuacct)),
    cur_mask(cfg.conf_cpus, fd_cpuset)
  {
    if (this->fd_cpuset < 0)
      std::cout << "WARN: couldn't open '" << root_cpuset << "/cpuset.cpus'." << std::endl;
    else
      std::cout << "INFO: opened '" << root_cpuset << "/cpuset.cpus'." << std::endl;

    if (this->fd_cpuacct < 0)
      std::cout << "WARN: couldn't open '" << root_cpuacct << "/cpuacct.usage'." << std::endl;
    else
      std::cout << "INFO: opened '" << root_cpuacct << "/cpuacct.usage'." << std::endl;
  }
  virtual ~cgroup_abstract_t(void);

  // will be overwritten by leaves!
  virtual int enable_cpus(cpuset_t mask) = 0;
};



class cgroup_leave : public cgroup_abstract_t {
private:
// --- parent and children
  cgroup_abstract_t& parent;

// --- policies and own mask
  struct timespec decr_migration_cost;
  struct timespec incr_migration_cost;
  cpuset_maskset_t *own_policy;
  policy *phandle;

  int check_fd_cpuset(void)
  {
    if (__builtin_expect(this->fd_cpuset < 0, 0)) {
      this->fd_cpuset = open((this->root_cpuset + "/cpuset.cpus").c_str(), O_RDWR);
      // since fd_cpuset was not yet open, we need to set the initial state here!
      lseek(this->fd_cpuset, 0, SEEK_SET);
      if (this->own_policy
          && this->fd_cpuset >= 0)
        this->own_policy->set_min_bitmask(this->fd_cpuset);
    }
    return this->fd_cpuset;
  }

  int set_cpu_mask(cpuset_t& cpu_mask)
  {
    if (this->check_fd_cpuset() >= 0) {
      lseek(this->fd_cpuset, 0, SEEK_SET);
      //unsigned ctr = 100; // we try it 100 times!
      do {
        int ret = cpu_mask.get_bitmask(this->fd_cpuset);
        if (ret)
          return ret;

        lseek(this->fd_cpuset, 0, SEEK_SET);
        cpuset_t tmp_mask(this->cfg.conf_cpus, this->fd_cpuset);
        if (tmp_mask == cpu_mask)
          return 0;
      } while (1);
      std::cout << "WARN: couldn't set mask!" << std::endl;
    }
    return -1;
  }

public:
  cgroup_leave(config_t& cfg,
               std::string& root_cpuset,
               std::string& root_cpuacct,
               cgroup_abstract_t& parent) :
    cgroup_abstract_t(cfg, root_cpuset, root_cpuacct),
    parent(parent),
    own_policy(nullptr)
  {
    this->decr_migration_cost.tv_sec = 0;
    this->decr_migration_cost.tv_nsec = 0;
    this->incr_migration_cost.tv_sec = 0;
    this->incr_migration_cost.tv_nsec = 0;
  }

  virtual ~cgroup_leave(void)
  {
    // write initial mask back (only if we really read it once)
    if (this->own_policy
        && this->fd_cpuset > 0)
      this->set_cpu_mask(this->own_policy->get_min());

    if(this->phandle)
      delete this->phandle;

    close(this->fd_cpuset);
    close(this->fd_cpuacct);

    std::cout << "INFO: closed cgroup '" << root_cpuset << "'\n";
    std::cout << "INFO: closed cgroup '" << root_cpuacct << "'";
    std::cout << std::endl;
  }

  void set_policy(cpuset_maskset_t *own_policy)
  {
    this->own_policy = own_policy;
    this->own_policy->set_min_bitmask(cur_mask);
    this->phandle = new policy(this->cfg,
                               this->decr_migration_cost,
                               this->incr_migration_cost);
  }

  int enable_cpus(cpuset_t mask)
  {
    int ret = 0;
    if (! this->parent.enable_cpus(mask)
        && mask != this->cur_mask) {
      auto tmp_mask = this->cur_mask | mask;

      struct timespec bgn, end;
      time_get(&bgn);
      if (!(ret = this->set_cpu_mask(tmp_mask)))
        this->cur_mask = tmp_mask;
      time_get(&end);
      time_diff(&incr_migration_cost, &end, &bgn);

      if (this->cfg.verbose && ! ret) {
        std::cout << "INFO: migration cost (increase to " << this->cur_mask.get_corec() << "): ";
        std::cout << this->incr_migration_cost.tv_sec << " sec. ";
        std::cout << this->incr_migration_cost.tv_nsec << " nsec." << std::endl;
      }
    }
    return ret;
  }

  int disable_cpus(cpuset_t mask)
  {
    int ret = 0;
    auto tmp = (this->cur_mask & mask);
    for (auto it = children.begin(); it != children.end(); ++it)
      ret |= it->second->disable_cpus(tmp);

    if (!ret
        && mask != this->cur_mask) {
      auto tmp_mask = tmp | this->own_policy->get_min(); // bug! -> own_policy!

      struct timespec bgn, end;
      time_get(&bgn);
      if ( ! (ret = this->set_cpu_mask(tmp_mask)) )
        this->cur_mask = tmp_mask;
      time_get(&end);
      time_diff(&decr_migration_cost, &end, &bgn);

      if (this->cfg.verbose && ! ret) {
        std::cout << "INFO: migration cost (decrease to " << this->cur_mask.get_corec() << "): ";
        std::cout << this->decr_migration_cost.tv_sec << " sec. ";
        std::cout << this->decr_migration_cost.tv_nsec << " nsec." << std::endl;
      }
    }

    return ret;
  }

  void update(unsigned long long ttotal_ns, struct timespec& tcur)
  {
    for (auto c_it = this->children.begin(); c_it != this->children.end(); ++c_it)
      c_it->second->update(ttotal_ns, tcur);

    if (this->phandle
        && this->check_fd_cpuset() > 0
        && this->check_fd_cpuacct() > 0) {

      lseek(this->fd_cpuset, 0, SEEK_SET);
      this->cur_mask.set_bitmask(this->fd_cpuset);
      unsigned long long ton_ns = this->read_cpu_usage(this->fd_cpuacct);

      // get the current diff
      ton_ns -= this->old_cpu_usage_ns;
      this->old_cpu_usage_ns += ton_ns;

      // after children are handled -> take care about ours!
      // -> children always need to be done first!
      unsigned cur_cpu_count = this->cur_mask.get_corec();
      auto pol_min = this->own_policy->get_min();
      unsigned vout = this->phandle->handle_policy(ttotal_ns, ton_ns,
                                                   pol_min.get_corec(),
                                                   this->own_policy->get_opt_cpu_count(),
                                                   cur_cpu_count);

      if (vout > cur_cpu_count) {
        cpuset_t t(this->cfg.conf_cpus, true);
        t ^= this->cur_mask;
        auto leftover = this->cfg.cpu_stats->lowest_utilization(t, vout - cur_cpu_count);

        this->enable_cpus(leftover);
        this->phandle->feedback_incr(vout);
        
        if (this->cfg.log.is_open())
          this->cfg.log << tcur.tv_sec << ";" << tcur.tv_nsec << ";" << this->cur_mask.get_corec() << "\n";
      }
      else if (vout < cur_cpu_count) {
        auto leftover = this->cur_mask ^ pol_min;
        leftover = leftover.mask_cores(vout - pol_min.get_corec());

        this->disable_cpus(leftover);
        this->phandle->feedback_decr(vout);

        if (this->cfg.log.is_open())
          this->cfg.log << tcur.tv_sec << ";" << tcur.tv_nsec << ";" << this->cur_mask.get_corec() << "\n";
      }
    }
  }

  void append(const char *path, cpuset_maskset_t *own_policy)
  {
    const char *bgn_str_it = path;
    const char *end_str_it = folder_iter(bgn_str_it);

    std::string fname;
    if (end_str_it != nullptr)
      fname.assign(bgn_str_it, end_str_it - bgn_str_it);
    else
      fname = bgn_str_it;

    cgroup_leave *child;
    auto children_it = this->children.find(fname);
    if (children_it != this->children.end())
      child = children_it->second;
    else {
      std::string cpuset_abs_fname = this->root_cpuset + "/" + fname;
      std::string cpuacct_abs_fname = this->root_cpuacct + "/" + fname;
      child = this->children[ fname ] = new cgroup_leave(this->cfg,
                                                         cpuset_abs_fname,
                                                         cpuacct_abs_fname,
                                                         *this);
    }

    const char *next_path = path + fname.size() + 1; // +1 as of /
    if (next_path < path + strlen(path))
      child->append(next_path, own_policy);
    else
      child->set_policy(own_policy);
  }

  void verify_policy(cpuset_maskset_t *parent_policy)
  {
    if (!this->own_policy)
      this->own_policy = parent_policy;

    for (auto c_it = this->children.begin(); c_it != this->children.end(); ++c_it)
      c_it->second->verify_policy(this->own_policy);

    if (this->own_policy) {
      if (!this->set_cpu_mask(this->own_policy->get_min()))
        this->cur_mask = this->own_policy->get_min();
      else
        std::cerr << "ERR: setting minimum was not successful!" << std::endl;
    }
  }
};



class cgroup_root : public cgroup_abstract_t {
private:
  struct timespec tcur;
  uint64_t daemon_avg_sum_ns;
  unsigned daemon_avg_ctr;

public:
  cgroup_root(config_t& cfg,
              std::string& root_cpuset,
              std::string& root_cpuacct) :
    cgroup_abstract_t(cfg, root_cpuset, root_cpuacct),
    daemon_avg_sum_ns(0),
    daemon_avg_ctr(0)
  {
    time_get(&tcur);
  }
  virtual ~cgroup_root(void)
  {
    std::cout << "INFO: closed cgroup '" << root_cpuset << "'" << std::endl;
    close(this->fd_cpuset);
    std::cout << "INFO: closed cgroup '" << root_cpuacct << "'" << std::endl;
    close(this->fd_cpuacct);
  }

  int enable_cpus(cpuset_t mask)
  {
    return 0;
  }

  void update(void)
  {
    struct timespec tdiff, tcur;
    time_get(&tcur);
    time_diff(&tdiff, &tcur, &this->tcur);
    this->tcur = tcur;
    unsigned long long ttotal_ns = tdiff.tv_sec * NSEC_PER_SEC + tdiff.tv_nsec;

    for (auto c_it = this->children.begin(); c_it != this->children.end(); ++c_it)
      c_it->second->update(ttotal_ns, tcur);

    uint64_t delay_ns = this->cfg.poll_delay.tv_sec * NSEC_PER_SEC + this->cfg.poll_delay.tv_nsec;
    this->daemon_avg_sum_ns += ttotal_ns - delay_ns;
    if ((++this->daemon_avg_ctr * delay_ns) >= 60E9) { // each minute
      std::cout << "INFO: daemon avg. time: ";
      std::cout << (this->daemon_avg_sum_ns / this->daemon_avg_ctr) << " ns" << std::endl;
      this->daemon_avg_sum_ns = 0;
      this->daemon_avg_ctr = 0;
    }
  }

  void append(const char *path, cpuset_maskset_t *own_policy)
  {
    const char *bgn_str_it = path;
    const char *end_str_it = folder_iter(bgn_str_it);

    std::string fname;
    if (end_str_it != nullptr)
      fname.assign(bgn_str_it, end_str_it - bgn_str_it);
    else
      fname = bgn_str_it;

    cgroup_leave *child;
    auto children_it = this->children.find(fname);
    if (children_it != this->children.end())
      child = children_it->second;
    else {
      std::string cpuset_abs_fname = this->root_cpuset + "/" + fname;
      std::string cpuacct_abs_fname = this->root_cpuacct + "/" + fname;
      child = this->children[ fname ] = new cgroup_leave(this->cfg,
                                                         cpuset_abs_fname,
                                                         cpuacct_abs_fname,
                                                         *this);
    }

    const char *next_path = path + fname.size() + 1; // +1 as of /
    if (next_path < path + strlen(path))
      child->append(next_path, own_policy);
    else
      child->set_policy(own_policy);
  }

  void verify_policy(void)
  {
    for (auto c_it = this->children.begin(); c_it != this->children.end(); ++c_it)
      c_it->second->verify_policy(nullptr);
  }
};

static inline cgroup_root * tree_build(config_t& cfg,
                                       std::unordered_map< std::string, class cpuset_maskset_t* >& mask_policies)
{
  cgroup_root *root = new cgroup_root(cfg, cfg.root_cpuset, cfg.root_cpuacct);

  for (auto it = mask_policies.begin(); it != mask_policies.end(); ++it)
    root->append((char*)it->first.c_str() + 1, it->second);  // +1 as of /

  // sets policies to parents if none existent
  root->verify_policy();

  return root;
}

#endif /* #ifndef _RECTREE_H_ */

