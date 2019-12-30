#include <sys/time.h>    
void gettimeofday_(int *tim, int *tmz)
{
  struct timeval tv;
  struct timezone tz;

    gettimeofday(&tv, &tz);
    tim[0] = tv.tv_sec;
    tim[1] = tv.tv_usec;  
    tmz[0] = tz.tz_minuteswest;
    tmz[1] = tz.tz_dsttime;

}

