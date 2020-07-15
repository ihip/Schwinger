#include <unistd.h>

void gethostid_(int *hid)
{
    hid[0] = gethostid();
}

