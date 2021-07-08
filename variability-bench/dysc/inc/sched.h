#ifndef _SCHED_H_
#define _SCHED_H_

#include <iostream>
#include <iomanip> //setw()
#include <cstdint>
#include <cstring>
#include <cassert>
#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif /* #ifndef _GNU_SOURCE */
#include <sched.h>
#include "config.h"

class cpuset_t {
private:
  unsigned int conf_cpus;
  cpu_set_t cpuset;

  void __bitmask_to_string(std::string& str, unsigned bitbgn, unsigned bitend)
  {
    if (str.size())
      str += ",";
    str += std::to_string(bitbgn);
    if (bitbgn != bitend) {
      str += "-";
      str += std::to_string(bitend);
    }
  }

public:
  ~cpuset_t(void)
  {
  }
  cpuset_t(unsigned int conf_cpus, const char *mask) :
    conf_cpus(conf_cpus)
  {
    if (mask)
      this->set_bitmask(mask);
    else
      CPU_ZERO(&this->cpuset);
  }

  cpuset_t(unsigned int conf_cpus, bool ones) :
    conf_cpus(conf_cpus)
  {
    if (ones)
      for (unsigned c = 0; c < conf_cpus; ++c)
        CPU_SET(c, &this->cpuset);
    else
      CPU_ZERO(&this->cpuset);
  }

  cpuset_t(unsigned int conf_cpus, int fd) :
    conf_cpus(conf_cpus)
  {
    this->set_bitmask(fd);
  }

  // used in operator overloading
  explicit cpuset_t(unsigned int conf_cpus) :
    conf_cpus(conf_cpus)
  {
    CPU_ZERO(&this->cpuset);
  }

  int isset_bit(int c)
  {
    return CPU_ISSET(c, &this->cpuset);
  }

  void set_bit(int c)
  {
    CPU_SET(c, &this->cpuset);
  }

  // end is included
  void set_bitrange(unsigned bgn, unsigned end)
  {
    for (unsigned c = bgn; c <= end; ++c)
      CPU_SET(c, &this->cpuset);
  }

  int set_bitmask(const char *mask)
  {
    unsigned int bgn = 0, end = 0;
    unsigned int *cur = &bgn;
    const char *it = mask;
    CPU_ZERO(&this->cpuset);

    if (!mask)
      return -1;

    while (1) {
      switch (*it) {
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
        (*cur) *= 10;
        (*cur) += (*it - '0');
        break;
      case '-':
        cur = &end;
        break;
      case ',':
        if (cur == &bgn)
          this->set_bit(bgn);
        else {
          this->set_bitrange(bgn, end);
          end = 0;
          cur = &bgn;
        }
        bgn = 0;
        break;
      case '\0':
      case '\n':
        goto done;
      default: // shouldn't happen
        std::cerr << "ERR:" << __FUNCTION__ << ":" << __LINE__ << ": '" << mask << "' occured!" << std::endl;
        return -1;
      }
      ++it;
    }
 done:
    if (cur == &bgn)
      this->set_bit(bgn);
    else
      this->set_bitrange(bgn, end);

    return 0;
  }

  int set_bitmask(int fd)
  {
    unsigned int bgn = 0, end = 0;
    unsigned int *cur = &bgn;
    CPU_ZERO(&this->cpuset);

    if (fd < 0)
      return -1;

    // created a buffer just to mitigate a lot of read calls
    char buf[ 256 ];
    int len;
    do {
      len = read(fd, buf, sizeof(buf));
      for (char *c = buf; c < buf + len; ++c) {
        switch (*c) {
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
          (*cur) *= 10;
          (*cur) += (*c - '0');
          break;
        case '-':
          cur = &end;
          break;
        case ',':
          if (cur == &bgn)
            this->set_bit(bgn);
          else {
            this->set_bitrange(bgn, end);
            end = 0;
            cur = &bgn;
          }
          bgn = 0;
          break;
        case '\0':
        case '\n':
          goto done;
        default: // shouldn't happen
          std::cerr << "ERR:" << __FUNCTION__ << ":" << __LINE__ << ": '" << *c << "' occured!" << std::endl;
          return -1;
        }
      }
    } while (len == sizeof(buf));
 done:
    if (cur == &bgn)
      this->set_bit(bgn);
    else
      this->set_bitrange(bgn, end);

    return 0;
  }

  cpu_set_t * get_bitmask(void)
  {
    return &this->cpuset;
  }

  int get_bitmask(int fd)
  {
    if (fd < 0)
      return -1;

    std::string tmp;
    this->bitmask_to_string(tmp);

    int ret = write(fd, tmp.c_str(), tmp.size());
    if (ret < (int)tmp.size()) {
      std::cerr << "WARN: couldn't write '" << tmp << "' to fd '" << fd << "'" << std::endl;
      return -1;
    }
    return 0;
  }

  void bitmask_to_string(std::string& str)
  {
    unsigned pos = 0;
    bool isvalid = false;
    for (unsigned c = 0; c < this->conf_cpus; ++c) {
      if (CPU_ISSET(c, &this->cpuset)) {
        if (!isvalid) {
          pos = c;
          isvalid = true;
        }
      }
      else if (isvalid) {
        __bitmask_to_string(str, pos, c - 1);
        isvalid = false;
      }
    }
    if (isvalid)
      __bitmask_to_string(str, pos, this->conf_cpus - 1);
  }

  cpuset_t operator^(cpuset_t& rhs)
  {
    cpuset_t next(this->conf_cpus);
    CPU_XOR(next.get_bitmask(), rhs.get_bitmask(), &this->cpuset);
    return next;
  }
  cpuset_t operator^=(cpuset_t& rhs)
  {
    CPU_XOR(&this->cpuset, rhs.get_bitmask(), &this->cpuset);
    return *this;
  }
  bool operator!=(cpuset_t& rhs)
  {
    return !CPU_EQUAL(&this->cpuset, rhs.get_bitmask());
  }
  bool operator==(cpuset_t& rhs)
  {
    return CPU_EQUAL(&this->cpuset, rhs.get_bitmask());
  }
  cpuset_t operator&(cpuset_t& rhs)
  {
    cpuset_t next(this->conf_cpus);
    CPU_AND(next.get_bitmask(), rhs.get_bitmask(), &this->cpuset);
    return next;
  }
  cpuset_t operator|(cpuset_t& rhs)
  {
    cpuset_t next(this->conf_cpus);
    CPU_OR(next.get_bitmask(), rhs.get_bitmask(), &this->cpuset);
    return next;
  }
  cpuset_t operator|=(cpuset_t& rhs)
  {
    CPU_OR(&this->cpuset, rhs.get_bitmask(), &this->cpuset);
    return *this;
  }

  cpuset_t mask_cores(unsigned amount)
  {
    cpuset_t next(this->conf_cpus);
    for (unsigned c = 0; c < this->conf_cpus && amount; ++c)
      if (CPU_ISSET(c, &this->cpuset)) {
        next.set_bit(c);
        --amount;
      }
    return next;
  }

  unsigned get_corec(void)
  {
    return CPU_COUNT(&this->cpuset);
  }

  void show_bitmask(void)
  {
    for (unsigned int c = 0; c < this->conf_cpus; c += 8) {
      uint8_t cur = 0x0ull;
      for (unsigned int i = 0; i < 8; ++i) {
        cur |= ((CPU_ISSET(c + i, &this->cpuset) & 0x1) << (7 - i));
      }
      std::cout << std::hex << std::setfill('0') << std::setw(2) << (unsigned)cur << std::dec;
    }
  }
};

struct cpuset_maskset_t {
private:
  cpuset_t min_mask;
  cpuset_t smt_mask;
  unsigned int opt_cpus_per_socket;
  config_t& cfg;

public:
  explicit cpuset_maskset_t(config_t& cfg) :
    min_mask(cfg.conf_cpus),
    smt_mask(cfg.conf_cpus),
    opt_cpus_per_socket(0),
    cfg(cfg)
  {
    for (unsigned i = 0; i < cfg.conf_cpus; i += 4)
      for (unsigned j = 0; j < cfg.smt_mode; ++j)
        this->smt_mask.set_bit(i + j);
  }
  ~cpuset_maskset_t(void)
  {
  }

  int set_min_bitmask(cpuset_t& cpuset)
  {
    this->min_mask = cpuset;

    std::cout << "INFO: min: ";
    this->min_mask.show_bitmask();
    std::string str = "";
    min_mask.bitmask_to_string(str);
    std::cout << "---> " << str << std::endl;

    this->set_max_bitmask(this->opt_cpus_per_socket);
    return 0;
  }

  int set_min_bitmask(int fd)
  {
    std::cout << "INFO: min: ";
    this->min_mask.show_bitmask();
    std::string str = "";
    min_mask.bitmask_to_string(str);
    std::cout << "---> " << str << "\n";

    int ret = this->min_mask.set_bitmask(fd);
    if (!ret)
      this->set_max_bitmask(this->opt_cpus_per_socket);
    return ret;
  }

  void set_max_bitmask(unsigned opt_cpus_per_socket)
  {
    this->opt_cpus_per_socket = opt_cpus_per_socket;
    
    std::cout << "INFO: max: additional " << opt_cpus_per_socket << " cpu(s) per socket\n";
  }

  unsigned get_opt_cpu_count(void)
  {
    return this->opt_cpus_per_socket * this->cfg.conf_sockets;
  }

  cpuset_t& get_min(void)
  {
    return this->min_mask;
  }
};

#endif /* #ifndef _SCHED_H_ */
