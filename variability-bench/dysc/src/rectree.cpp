#include "../inc/rectree.h"

cgroup_abstract_t::~cgroup_abstract_t(void)
{
  for (auto it = this->children.begin(); it != this->children.end(); ++it)
    delete it->second;
  if (this->fd_cpuset > 0)
    close(this->fd_cpuset);
  if (this->fd_cpuacct > 0)
    close(this->fd_cpuacct);
}
