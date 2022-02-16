#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <limits.h>

#include "checkCPUtime.h"


int checkCPUtime() {
  tms myTimes; /* Variable qui va contenir les temps */
  long       ticks;
  double     sec;
  times( &myTimes );  /* Demande des valeurs */
  ticks = myTimes.tms_utime;
  sec = (double)ticks/(double)sysconf(_SC_CLK_TCK);
  return (int)sec;
}
