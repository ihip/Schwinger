#include <unistd.h>

int gethostname_(char* name, int* len) {
  int ierr, i;

  for(i = 0; i < *len; i++) name[i] = 0;

  ierr = gethostname(name, *len);

  return ierr;
}
