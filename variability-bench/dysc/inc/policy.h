#ifndef _POLICY_H_
#define _POLICY_H_

#include <cstdint>
#include <vector>
#include <time.h>
#include "config.h"
#include "sched.h"



class policy {
private:
  config_t& cfg;
  struct timespec& decr_migration_cost;
  struct timespec& incr_migration_cost;

  std::vector< unsigned > history;
  std::vector< unsigned >::iterator h_cur;

  
  unsigned speculate_next(void)
  {
    unsigned* vals = history.data();
    int elem = (h_cur - this->history.begin());
    
    unsigned cnt = 5;
    
    unsigned x_avg = 0;
    unsigned y_avg = 0;
    for(unsigned i = 1; i <= cnt; ++i) {
      x_avg += i;
      y_avg += vals[ (elem + history.size() - i) % history.size() ];
    }
    x_avg /= cnt;
    y_avg /= cnt;
    
    unsigned nenner = 0;
    unsigned zaehler = 0;
    for(unsigned x = 1; x <= cnt; ++x) {
      nenner += (vals[ (elem + history.size() - x) % history.size() ] - y_avg) * (x - x_avg);
      zaehler += x * x;
    }
    float alpha_1 = (float)nenner / (float)zaehler;

    float alpha_0 = y_avg - alpha_1 * x_avg;

    //std::cout << "  m = " << alpha_1 << ", c = " << alpha_0 << std::endl;   
    int next = (alpha_0 + alpha_1 * (cnt + 1));

    return (next > 0) ? next : 0;
  }
  
  void add_val(unsigned val)
  {
    *this->h_cur = val;
    ++this->h_cur;
    if( this->h_cur == this->history.end() )
      this->h_cur = this->history.begin();
  }

// --- time it takes to modify cgroup
#if 0
  double max_prop_runtime;
  unsigned int daemon_ticks_since_last_sched;

  bool shall_schedule(void)
  {
    uint64_t daemon_delay_ns = this->cfg.poll_delay.tv_sec * NSEC_PER_SEC + this->cfg.poll_delay.tv_nsec;
    uint64_t sched_delay_ns = this->migration_cost.tv_sec * NSEC_PER_SEC + this->migration_cost.tv_nsec;

    return (daemon_delay_ns * this->daemon_ticks_since_last_sched * this->max_prop_runtime >= sched_delay_ns);
  }
#endif

public:
  policy(config_t& cfg,
         struct timespec& decr_migration_cost,
         struct timespec& incr_migration_cost) :
    cfg(cfg),
    decr_migration_cost(decr_migration_cost),
    incr_migration_cost(incr_migration_cost),
    history(10, 0.0),
    h_cur(history.begin())
#if 0
    , max_prop_runtime(0.10),
    daemon_ticks_since_last_sched(0)
#endif
  {
  }
  ~policy(void) {
  
  }

  unsigned handle_policy(unsigned long long ttotal_ns,
                         unsigned long long ton_ns,
                         unsigned def_smt_c,
                         unsigned opt_smt_c,
                         unsigned cur_smt_c)
  {
#if 0
    ++this->daemon_ticks_since_last_sched;
    if (this->shall_schedule())
#endif
    {

      unsigned max_smt_c = opt_smt_c + def_smt_c;

      // ton_ns / ttotal_ns is the percentage related to all cores in the system (100% means 1 core, 200% means 2 cores, ...)
      // hence: (ton_ns/ttotal_ns) / current cores -> provides scale between 0 and 100% (still related to all cores)
      unsigned vout = ceil((((float)ton_ns * (float)max_smt_c) / ((float)ttotal_ns * (float)def_smt_c)));

      if (vout < def_smt_c)
        vout = def_smt_c;
      else if (vout > max_smt_c)
        vout = max_smt_c;

#if 0
      if(vout != cur_smt_c)
        this->daemon_ticks_since_last_sched = 0;
#endif
      
      add_val(vout);
      
      unsigned val = this->speculate_next();
     // std::cout << "spec: " << val << " inserting: " << vout << "(" << max_smt_c << ")" << std::endl;
      if(val < def_smt_c)
        val = def_smt_c;
      if(val > max_smt_c) {
        val = max_smt_c;
      }
      
      return val;
    }
#if 0
    return cur_smt_c;
#endif
  }
  
  void feedback_decr(unsigned vout)
  {
    uint64_t daemon_delay_ns = this->cfg.poll_delay.tv_sec * NSEC_PER_SEC + this->cfg.poll_delay.tv_nsec;
    uint64_t sched_delay_ns = this->decr_migration_cost.tv_sec * NSEC_PER_SEC + this->decr_migration_cost.tv_nsec;
    
    unsigned steps = ceil((float)sched_delay_ns / (float)daemon_delay_ns);
    //std::cout << "decr: " << ((float)sched_delay_ns / (float)daemon_delay_ns) << " " << steps << std::endl;
    for(unsigned i = 0; i < steps; ++i)
      this->add_val(vout);
  }
  
  void feedback_incr(unsigned vout)
  {
    uint64_t daemon_delay_ns = this->cfg.poll_delay.tv_sec * NSEC_PER_SEC + this->cfg.poll_delay.tv_nsec;
    uint64_t sched_delay_ns = this->incr_migration_cost.tv_sec * NSEC_PER_SEC + this->incr_migration_cost.tv_nsec;

    unsigned steps = ceil((float)sched_delay_ns / (float)daemon_delay_ns);
    //std::cout << "incr: " << ((float)sched_delay_ns / (float)daemon_delay_ns) << " " << steps << std::endl;
    for(unsigned i = 0; i < steps; ++i)
      this->add_val(vout);
  }
};

#endif /* #ifndef _POLICY_H_ */

