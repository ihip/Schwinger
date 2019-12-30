#include <sys/time.h>    
void gethostid_(int *hhh)
{
    hhh[0] = gethostid();
}

