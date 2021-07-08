#include <cstring>
#include <iostream>
#include <unordered_set>
#include <regex.h> // using old POSIX regex, as C++ regex is experimental on WSC
#include <getopt.h>
#include <unistd.h>
#include <signal.h>
#include <libgen.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include "../inc/rectree.h"
#include "../inc/config.h"

bool folder_exists(const char *str)
{
  struct stat buf;
  return ! ((stat(str, &buf) != 0)
            || ! (buf.st_mode & S_IFDIR));
}

volatile bool poll_daemon;
volatile int irpt_id;
volatile int irpt_param;
void (*prev_hdl_sigterm)(int);
void (*prev_hdl_sigint)(int);
void hdl_sigterm(int param)
{
  irpt_id = SIGTERM;
  irpt_param = param;
  poll_daemon = false;
}
void hdl_sigint(int param)
{
  irpt_id = SIGINT;
  irpt_param = param;
  poll_daemon = false;
}

void usage(char *exename, config_t& cfg)
{
  std::cout << "Usage: " << exename << " [OPTION]...\n"
    "\n"
    "  -D              run binary headless as daemon.\n"
    "  -i TIME_IN_US   polling time of daemon.\n"
    "                  default: " << (cfg.poll_delay.tv_sec * NSEC_PER_SEC + cfg.poll_delay.tv_nsec) << "us\n"
    "  -s PATH         path to root cpuset directory.\n"
    "                  default: '" << cfg.root_cpuset << "'\n"
    "  -a PATH         path to root cpuacct directory.\n"
    "                  default: '" << cfg.root_cpuacct << "'\n"
    "  -c FILE         sysconfig file such as /etc/sysconfig/dysc.\n"
    "  -o FILE         write a trace that captures the cgroup changes.\n"
    "  -v              verbose mode prints the migration cost.\n"
    "  -h              print help.\n";
  std::cout << std::flush;
}


int
mkpath(const char *s, mode_t mode){
    char *q, *r = NULL, *path = NULL, *up = NULL;
    int rv;

    rv = -1;
    if (strcmp(s, ".") == 0 || strcmp(s, "/") == 0)
        return (0);

    if ((path = strdup(s)) == NULL)
        exit(1);
 
    if ((q = strdup(s)) == NULL)
        exit(1);

    if ((r = dirname(q)) == NULL)
        goto out;
    
    if ((up = strdup(r)) == NULL)
        exit(1);

    if ((mkpath(up, mode) == -1) && (errno != EEXIST))
        goto out;

    if ((mkdir(path, mode) == -1) && (errno != EEXIST))
        rv = -1;
    else
        rv = 0;

out:
    if (up != NULL)
        free(up);
    free(q);
    free(path);
    return (rv);
}



int parse_cfgfile(std::ifstream& cfg_stream, config_t& cfg,
                  std::unordered_map< std::string, class cpuset_maskset_t* >& mask_policies)
{
  regex_t regex_smt_mode;
  const char *r_smt_mode = " *SMT_MODE *= *([1,2,4])";
  if (regcomp(&regex_smt_mode, r_smt_mode, REG_EXTENDED)) {
    std::cerr << "ERR: couldn't generate regex '" << r_smt_mode << "'" << std::endl;
    return -1;
  }

#if 0
  regex_t regex_cgroup;
  const char *r_cgroup = " *DYSC__CGROUP *= *(.*)";
  if (regcomp(&regex_cgroup, r_cgroup, REG_EXTENDED)) {
    std::cerr << "ERR: couldn't generate regex '" << r_cgroup << "'" << std::endl;
    return -1;
  }
#endif

  regex_t regex_cgroup_cpus;
  const char *r_cgroup_cpus = " *DYSC__GROUP_OPT_CPUS *= *([0-9]*)";
  if (regcomp(&regex_cgroup_cpus, r_cgroup_cpus, REG_EXTENDED)) {
    std::cerr << "ERR: couldn't generate regex '" << r_cgroup_cpus << "'" << std::endl;
    return -1;
  }

  regex_t regex_fileoutput;
  const char *r_fileoutput = " *FILEOUTPUT *= *(.*)";
  if (regcomp(&regex_fileoutput, r_fileoutput, REG_EXTENDED)) {
    std::cerr << "ERR: couldn't generate regex '" << r_fileoutput << "'" << std::endl;
    return -1;
  }

  std::string str;
  class cpuset_maskset_t *cur_cgroup = mask_policies[ "/dysc_mmfsd_iso" ] = new cpuset_maskset_t(cfg);
  while (getline(cfg_stream, str)) {
    regmatch_t pmatch[2]; // no. 1 is entire, no. 2 is group

    if (regexec(&regex_fileoutput, str.c_str(), 2, pmatch, REG_NOTBOL) != REG_NOMATCH) {
      std::string f(str.begin() + pmatch[1].rm_so, str.begin() + pmatch[1].rm_eo);
      std::cout << "INFO: cfg: output file:            '" << f << "'" << std::endl;
      if(f[0] != '\0') {
        int offset = f.length() - 1;

        while(f[offset] != '/'
              && offset != 0)
          --offset;

        std::string fldr(f.c_str(), f.c_str() + offset);
        std::cout << "INFO: creating folder:  " << fldr << std::endl;
        std::cout << "INFO: opening log file: " << f << std::endl;
        mkpath(fldr.c_str(), ACCESSPERMS);
        cfg.log.open(f.c_str());
      }
    }
    else if (regexec(&regex_smt_mode, str.c_str(), 2, pmatch, REG_NOTBOL) != REG_NOMATCH) {
      cfg.smt_mode = *(str.c_str() + pmatch[1].rm_so) - '0';
      std::cout << "INFO: cfg: smt mode:               " << cfg.smt_mode << std::endl;
    }
#if 0
    else if (regexec(&regex_cgroup, str.c_str(), 2, pmatch, REG_NOTBOL) != REG_NOMATCH) {
      std::string f(str.begin() + pmatch[1].rm_so, str.begin() + pmatch[1].rm_eo);
      std::cout << "INFO: cfg: policy for cgroup:      " << f << std::endl;
      cur_cgroup = mask_policies[ f ] = new cpuset_maskset_t(cfg);
    }
#endif
    else if (cur_cgroup) {
      if (regexec(&regex_cgroup_cpus, str.c_str(), 2, pmatch, REG_NOTBOL) != REG_NOMATCH) {
        std::string f(str.begin() + pmatch[1].rm_so, str.begin() + pmatch[1].rm_eo);
        // per group one cpu ..
        int opt_cpus_per_socket = atoi(f.c_str());
        std::cout << "INFO: cfg: > opt. cpus per socket: " << opt_cpus_per_socket << std::endl;
        cur_cgroup->set_max_bitmask(opt_cpus_per_socket);
        cur_cgroup = nullptr;
      }
    }
  }

  regfree(&regex_smt_mode);
#if 0
  regfree(&regex_cgroup);
#endif
  regfree(&regex_cgroup_cpus);
  regfree(&regex_fileoutput);
  return 0;
}



FILE* w_fh_acct;
FILE* r_fh_acct;
FILE* fh_set;
char *line = NULL;
size_t len = 0;

long fix_ctr = 0;

long fix_ctr_waiting = 100;
void fix_cpuacct(config_t& cfg, std::string dyncgroup_home)
{
  ++fix_ctr;
  if(fix_ctr > fix_ctr_waiting) // all 10 seconds~
    fix_ctr = 0;
  else
    return;

  if( ! r_fh_acct)
    r_fh_acct = fopen((cfg.root_cpuacct + dyncgroup_home + std::string("/tasks")).c_str(), "r");

  if( ! w_fh_acct)
    w_fh_acct = fopen((cfg.root_cpuacct + dyncgroup_home + std::string("/tasks")).c_str(), "w");

  if( ! fh_set)
    fh_set  = fopen((cfg.root_cpuset + dyncgroup_home + std::string("/tasks")).c_str(), "r");

  if( r_fh_acct && w_fh_acct && fh_set )
  {
    std::unordered_set<long long unsigned> all_pids;
    ssize_t nread;
    fseek(fh_set, 0, SEEK_SET);
    while((nread = getline(&line, &len, fh_set)) != -1) {
      long long unsigned pid = strtoull(line, NULL, 10);
      all_pids.insert(pid);
    }
    std::cout << "BUGFIX: pids cpuset '" << dyncgroup_home << "': " << all_pids.size() << std::endl;

    fseek(r_fh_acct, 0, SEEK_SET);
    while((nread == getline(&line, &len, r_fh_acct)) != -1) {
      long long unsigned pid = strtoull(line, NULL, 10);

      // further hack
      if(all_pids.find(pid) == all_pids.end())
        break;
      all_pids.erase(pid);
    }

    std::cout << "BUGFIX: pids cpuacct '" << dyncgroup_home << "'missing: " << all_pids.size() << std::endl;
    if(all_pids.size() < 20 )
      fix_ctr_waiting = 10000;
    for(auto it = all_pids.begin(); it != all_pids.end(); ++it)
    {
      std::string val = std::to_string(*it);
      try {
        fseek(w_fh_acct, 0, SEEK_SET);
        fputs(val.c_str(), w_fh_acct);
      }
      catch(const std::exception&) {}
    }
  }
}

void fix_cpuacct_close(void)
{
  if( r_fh_acct )
    fclose( r_fh_acct );
  if( w_fh_acct )
    fclose( w_fh_acct );
  if( fh_set )
    fclose( fh_set );
}




///////////////////////////////////////////////////////
/* checks if the string is purely an integer
 * we can do it with `strtol' also
 */
int check_if_number (char *str)
{
  int i;
  for (i=0; str[i] != '\0'; ++i)
    if (!isdigit (str[i]))
      return 0;
  return 1;
}
 
#define MAX_BUF 1024
#define PID_LIST_BLOCK 32
 
int *pidof (char *pname)
{
  DIR *dirp;
  FILE *fp;
  struct dirent *entry;
  int *pidlist, pidlist_index = 0, pidlist_realloc_count = 1;
  char path[MAX_BUF], read_buf[MAX_BUF];
 
  dirp = opendir ("/proc/");
  if (dirp == NULL)
  {
    perror ("Fail");
    return NULL;
  }
 
  pidlist = (int*) malloc (sizeof (int) * PID_LIST_BLOCK);
  if (pidlist == NULL)
  {
    return NULL;
  }
 
  while ((entry = readdir (dirp)) != NULL)
  {
    if (check_if_number (entry->d_name))
    {
      strcpy (path, "/proc/");
      strcat (path, entry->d_name);
      strcat (path, "/comm");
 
      /* A file may not exist, it may have been removed.
       * dut to termination of the process. Actually we need to
       * make sure the error is actually file does not exist to
       * be accurate.
       */
      fp = fopen (path, "r");
      if (fp != NULL)
      {
        fscanf (fp, "%s", read_buf);
        if (strcmp (read_buf, pname) == 0)
        {
          /* add to list and expand list if needed */
          pidlist[pidlist_index++] = atoi (entry->d_name);
          if (pidlist_index == PID_LIST_BLOCK * pidlist_realloc_count)
          {
            pidlist_realloc_count++;
            pidlist = (int*) realloc (pidlist, sizeof (int) * PID_LIST_BLOCK * pidlist_realloc_count); //Error check todo
            if (pidlist == NULL)
            {
              return NULL;
            }
          }
        }
        fclose (fp);
      }
    }
  }
 
  closedir (dirp);
  pidlist[pidlist_index] = -1; /* indicates end of list */
  return pidlist;
}
///////////////////////////////////////////////////////

int cgtaskset( std::string& cgname_taskfile, std::unordered_set< pid_t >& pids )
{
  int ret = 0;
  FILE* fh_cgset_name = fopen(cgname_taskfile.c_str(), "w");
  if( fh_cgset_name ) {
    for(auto it = pids.begin(); it != pids.end(); ++it) {
      std::string val = std::to_string(*it);
      try {
        fseek(fh_cgset_name, 0, SEEK_SET);
        if(fputs(val.c_str(), fh_cgset_name) == EOF) {
          ret = -1;
          break;
        }
      }
      catch(const std::exception&) {}
    }
    fclose( fh_cgset_name );
  }
  else
    ret = -1;

  return ret;
}

int cgcopy( std::string& from, std::string& to )
{
  int ret = 0;
  FILE* fh_from = fopen(from.c_str(), "r");
  if( fh_from ) {
    size_t len = 0;
    char *line = NULL;
    if (getline(&line, &len, fh_from) != -1) {
      FILE *fh_to = fopen( to.c_str(), "w" );
      if( fh_to ) {
        if( fputs( line, fh_to ) == EOF )
          ret = -1;
        fclose( fh_to );
      }
      free( line );
    }
    else
      ret = -1;
    fclose( fh_from );
  }
  else {
    std::cerr << "ERR: couln't open " << from << std::endl;
    return -1;
  }
  return ret;
}


int mmfsd_cgroup_init( config_t& cfg, const char* cgname, std::string& system_cgroup, std::unordered_set< pid_t >& csm_system_tasks )
{
  std::string s_cgset_name = cfg.root_cpuset;
  if( *(cfg.root_cpuset.end() - 1) != '/' )
    s_cgset_name += "/";
  s_cgset_name += cgname;

  std::string s_cgacct_name = cfg.root_cpuacct;
  if( *(cfg.root_cpuacct.end() - 1) != '/' )
    s_cgacct_name += "/";
  s_cgacct_name += cgname;
  
  // create dedicated cpuset and cpuacct
  if( mkpath( s_cgset_name.c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IXOTH | S_ISVTX  ) == -1 ) {
    std::cerr << "ERR: cgroup set '" << s_cgset_name << "' does not exist and can't be created" << std::endl;
    return -1;
  }
  if( mkpath( s_cgacct_name.c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IXOTH | S_ISVTX  ) == -1 ) {
    std::cerr << "ERR: cgroup acct '" << s_cgacct_name << "' does not exist and can't be created" << std::endl;
    return -1;
  }

  std::string s_cgset_csm_system_dir = cfg.root_cpuset;
  s_cgset_csm_system_dir += system_cgroup;
  
  // copy initial cpuset.cpus and cpuset.mems of /csm_system
  std::string s_cgset_csm_system_cpuset = s_cgset_csm_system_dir + "/cpuset.cpus";
  std::string s_cgset_name_cpuset = s_cgset_name + "/cpuset.cpus";
  if( cgcopy(s_cgset_csm_system_cpuset, s_cgset_name_cpuset) == -1 ) {
    std::cerr << "ERR: couln't get cpuset.cpus of intial cgroup" << std::endl;
    return -1;
  }

  std::string s_cgset_csm_system_mems = s_cgset_csm_system_dir + "/cpuset.mems";
  std::string s_cgset_name_mems = s_cgset_name + "/cpuset.mems";
  if( cgcopy(s_cgset_csm_system_mems, s_cgset_name_mems) == -1 ) {
    std::cerr << "ERR: couln't get cpuset.mems of intial cgroup" << std::endl;
    return -1;
  }
  
  // find all of the mmfsd_tasks
  std::unordered_set< pid_t > all_mmfsd_tasks;
  int* mmfsd_pids = pidof("mmfsd");
  for( int i = 0; mmfsd_pids[i] != -1; ++i ) {
    std::string s_mmfsd_tasks = "/proc/" + std::to_string(mmfsd_pids[i]) + "/task";
    all_mmfsd_tasks.insert(mmfsd_pids[i]);

    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (s_mmfsd_tasks.c_str())) != NULL) {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
        all_mmfsd_tasks.insert(atoi(ent->d_name));
      closedir (dir);
    }
  }
  free(mmfsd_pids);

  // before moving find the tasks / pids that are in /csm_system cpuset
  std::string s_csm_system_tasks = s_cgset_csm_system_dir + "/tasks";
  FILE* fh_csm_system_tasks = fopen(s_csm_system_tasks.c_str(), "r");
  if( fh_csm_system_tasks ) {
    char *line = NULL;
    size_t len = 0;
    while( getline( &line, &len, fh_csm_system_tasks ) != -1 ) {
      pid_t pid = atoi( line );
      if( all_mmfsd_tasks.find( pid ) != all_mmfsd_tasks.end() )
        csm_system_tasks.insert( pid );
    }
    if(line)
      free(line);
    
    fclose( fh_csm_system_tasks );
  }

  // move all tasks / pids into the new cgroup
  std::string s_cgset_taskfile = s_cgset_name + "/tasks";
  if(cgtaskset( s_cgset_taskfile, all_mmfsd_tasks ) == -1)
    std::cerr << "ERR: couldn't move pids to cgroup set" << std::endl;

  std::string s_cgacct_taskfile = s_cgacct_name + "/tasks";
  if(cgtaskset( s_cgacct_taskfile, all_mmfsd_tasks ) == -1)
    std::cerr << "ERR: couldn't move pids to cgroup acct" << std::endl;

  return 0;
}

void mmfsd_cgroup_finish( config_t& cfg, const char* cgname, std::string& system_cgroup, std::unordered_set< pid_t >& csm_system_tasks )
{
  std::string s_cgset_name = cfg.root_cpuset;
  if( *(cfg.root_cpuset.end() - 1) != '/' )
    s_cgset_name += "/";
  s_cgset_name += cgname;

  // cpuset
  // move all pids back to the root cgroup and the individual ones to /csm_system
  std::string s_cgset_tasks = s_cgset_name + "/tasks";
  FILE* fh_cgset_tasks = fopen(s_cgset_tasks.c_str(), "r");
  if( fh_cgset_tasks ) {
    char *line = NULL;
    size_t len = 0;
    std::unordered_set< pid_t > root_tasks;
    while( getline( &line, &len, fh_cgset_tasks ) != -1 ) {
      pid_t pid = atoi(line);
      if( csm_system_tasks.find( pid ) == csm_system_tasks.end() )
        root_tasks.insert( pid );
    }
    if(line)
      free(line);
    fclose( fh_cgset_tasks );
    
    std::string s_cgset_root_tasks = cfg.root_cpuset;
    s_cgset_root_tasks += "/tasks";
    if(cgtaskset( s_cgset_root_tasks, root_tasks ) == -1)
      std::cerr << "ERR: couldn't move pids back to /csm_system cgroup" << std::endl;
  }

  std::string s_csm_system_root_tasks = cfg.root_cpuset;
  s_csm_system_root_tasks += system_cgroup + "/tasks";
  if(cgtaskset( s_csm_system_root_tasks, csm_system_tasks ) == -1)
    std::cerr << "ERR: couldn't move pids back to /csm_system cgroup" << std::endl;


  std::string s_cgacct_name = cfg.root_cpuacct;
  if( *(cfg.root_cpuacct.end() - 1) != '/' )
    s_cgacct_name += "/";
  s_cgacct_name += cgname;

  // cpuacct
  std::string s_cgacct_tasks = s_cgacct_name + "/tasks";
  FILE* fh_cgacct_tasks = fopen(s_cgacct_tasks.c_str(), "r");
  if( fh_cgacct_tasks ) {
    char *line = NULL;
    size_t len = 0;
    std::unordered_set< pid_t > root_tasks;
    while( getline( &line, &len, fh_cgacct_tasks ) != -1 )
      root_tasks.insert( atoi(line) );
    if(line)
      free(line);
    fclose( fh_cgacct_tasks );

    std::string s_cgacct_root_tasks = cfg.root_cpuacct;
    s_cgacct_root_tasks += "/tasks";
    if(cgtaskset( s_cgacct_root_tasks, root_tasks ) == -1)
      std::cerr << "ERR: couldn't move pids back to /csm_system cgroup" << std::endl;
  }

  // delete cgroup folders
  rmdir( s_cgset_name.c_str() );
  rmdir( s_cgacct_name.c_str() );
}







//  cgroup -> prolog .... epilog -> cgroup
int main(int argc, char *argv[])
{
  config_t cfg;
  cfg.conf_cpus = sysconf(_SC_NPROCESSORS_CONF);
  cfg.conf_sockets = 2;
  cfg.smt_mode = 4;
  cfg.root_cpuset = "/sys/fs/cgroup/cpuset";
  cfg.root_cpuacct = "/sys/fs/cgroup/cpuacct";
  cfg.poll_delay.tv_sec = 1; // 1 sec
  cfg.poll_delay.tv_nsec = 0;
  cfg.cpu_stats = new stat_t(cfg.conf_cpus);
  cfg.verbose = false;

  long clk_tck = sysconf(_SC_CLK_TCK); // clock ticks per second
  if (clk_tck != -1) {
    unsigned long long ns = (NSEC_PER_SEC / clk_tck);
    cfg.poll_delay.tv_sec = ns / NSEC_PER_SEC;
    cfg.poll_delay.tv_nsec = ns % (unsigned long long)NSEC_PER_SEC; // mask
  }


  std::string system_cgroup("/");
  std::unordered_map< std::string, class cpuset_maskset_t* > mask_policies;
  int c;
  while ((c = getopt(argc, argv, "hi:s:a:o:p:c:v")) != -1)
    switch (c) {
    case 'i':
      cfg.poll_delay.tv_sec = atoi(optarg) / USEC_PER_SEC;
      cfg.poll_delay.tv_nsec = atoi(optarg) % (unsigned long long)NSEC_PER_USEC;
      break;
    case 's':
      cfg.root_cpuset = optarg;
      break;
    case 'a':
      cfg.root_cpuacct = optarg;
      break;
    case 'o':
      std::cout << "INFO: opening log file: " << optarg << std::endl;
      cfg.log.open(optarg);
      break;
    case 'v':
      cfg.verbose = true;
      break;
    case 'c':
    {
      std::ifstream cfg_stream(optarg);
      if (cfg_stream.is_open()) {
        parse_cfgfile(cfg_stream, cfg, mask_policies);
        cfg_stream.close();
      }
      break;
    }
    case 'p':
    {
      if (optarg[0] == '/')
        system_cgroup = optarg;
      else
        system_cgroup += std::string(optarg);

      // if /system cpuset exists, we should move ourself in there.
      if (folder_exists((std::string(cfg.root_cpuset) + system_cgroup).c_str())) {
        std::string tasks = std::string(cfg.root_cpuset) + system_cgroup + std::string("/tasks");
        int fd = open(tasks.c_str(), O_WRONLY);
        if (fd >= 0) {
          pid_t pid = syscall(SYS_getpid);
          std::string str = std::to_string(pid);
          std::cout << "INFO: appending own daemon pid " << str << " to cgroup '" << tasks << "'" << std::endl;
          write(fd, (char*)str.c_str(), str.size());
          close(fd);
        }
        else {
          std::cerr << "ERR: couldn't open '" << tasks << "': " << fd << std::endl;
          return -1;
        }
      }
      else {
        std::cout << "WARN: no " << system_cgroup << " cgroup defined! Can't move to it." << std::endl;
        return -1; // default
      }
      break;
    }
    case 'h':
    default:
      usage(argv[0], cfg);
      return -1;
    }

  if( ! strcmp(system_cgroup.c_str(), "/") ) {
    std::cout << "ERR: no system cgroup set!" << std::endl;
    return -1;
  }


  if (getuid()) {
    std::cerr << "ERR: daemon needs to run with root privileges!" << std::endl;
    return -1;
  }

  std::cout << "INFO: starting daemon '" << argv[0] << "'" << std::endl;
  std::cout << "INFO: polling delay is set to : ";
  const char *unit[] = {
    "us", "ms", "sec"
  };
  unsigned unit_id = 0;
  double t_poll_us = cfg.poll_delay.tv_sec * USEC_PER_SEC + cfg.poll_delay.tv_nsec / NSEC_PER_USEC;
  for (; unit_id < sizeof(unit) / sizeof(unit[0]); ) {
    if (t_poll_us >= 1000) {
      t_poll_us /= 1000;
      ++unit_id;
    }
    else
      break;
  }
  std::cout << t_poll_us << " " << unit[unit_id] << std::endl;

  poll_daemon = true;
  irpt_id = -1;
  prev_hdl_sigterm = signal(SIGTERM, hdl_sigterm);
  prev_hdl_sigint = signal(SIGINT, hdl_sigint);
  
  
  std::unordered_set< pid_t > csm_system_pids;
  const char *cgname = "/dysc_mmfsd_iso";
  if( mmfsd_cgroup_init( cfg, cgname, system_cgroup, csm_system_pids ) == -1 ) {
    std::cerr << "ERR: couldn't set up dedicated mmfsd cgroup" << std::endl;
    return -1;
  }

  cgroup_root *n = tree_build(cfg, mask_policies);
  do {
    if (!n) {
      std::cerr << "ERR: could not obtain next cgroup root state!" << std::endl;
      break;
    }
    n->update();
    // make absolutely sure to get our prefered delay time
    auto poll_delay = cfg.poll_delay; // will be overwritten
    while (clock_nanosleep(CLOCK_MONOTONIC, 0, &poll_delay, &poll_delay)) ;

    fix_cpuacct(cfg, cgname);

  } while (poll_daemon);

  if (cfg.log.is_open())
    cfg.log.close();

  // tree needs to be killed before policies!
  delete n;
  for (auto it = mask_policies.begin(); it != mask_policies.end(); ++it)
    delete it->second;

  switch (irpt_id) {
  case -1:
    break;
  case SIGTERM:
    prev_hdl_sigterm(irpt_param);
    break;
  case SIGINT:
    //prev_hdl_sigint(irpt_param);
    break;
  }
  
  mmfsd_cgroup_finish( cfg, cgname, system_cgroup, csm_system_pids );

  fix_cpuacct_close();

  std::cout << "INFO: exiting daemon" << std::endl;
  return 0;
}
