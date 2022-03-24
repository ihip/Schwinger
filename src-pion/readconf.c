// hip // 29. Nov 2007 //

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

readconf_(char *fname, double *u)
{
  FILE *fin;
  char filename[64];
  char read_buffer1[8], result1[8];
  char read_buffer2[8], result2[8];
  char read_buffer3[8], result3[8];
  char read_buffer4[8], result4[8];
  double dr1, di1, dr2, di2, a;
  int x, y, k, i;

  // copy fname to filename
  i = 0;
  while((fname[i] != 32) && (i < 63)) { 
    filename[i] = fname[i]; i++;
  }
  if(i == 63) {
    printf("readconf: file name to long!\n");
    exit(0);
  }
  filename[i] = 0;

  printf("readconf: filename = %s\n", filename);

  fin = fopen(filename, "r");
  if(fin == NULL) {
    printf("readconf: no such file!\n");
    exit(0);
  }
  printf("Reading in %s\n", filename);

  for(x = 0; x < NTIME; x++) 
    for(y = 0; y < NSPACE; y++) 
  {
    fread(&read_buffer1, sizeof(double), 1, fin);
    fread(&read_buffer2, sizeof(double), 1, fin);
    fread(&read_buffer3, sizeof(double), 1, fin);
    fread(&read_buffer4, sizeof(double), 1, fin);

    for(i = 0; i < sizeof(double); i++) {
      result1[sizeof(double) - 1 - i] = read_buffer1[i];
      result2[sizeof(double) - 1 - i] = read_buffer2[i];
      result3[sizeof(double) - 1 - i] = read_buffer3[i];
      result4[sizeof(double) - 1 - i] = read_buffer4[i];
     }
    dr1 = *((double *)result2);
    di1 = *((double *)result1);
    dr2 = *((double *)result4);
    di2 = *((double *)result3);

    k = NTIME * y + x;
    u[2 * k] = dr1; u[2 * k + 1] = di1;
    //u[512 + 2 * k] = dr2; u[512 + 2 * k + 1] = di2;
    u[NTIME * NSPACE * 2 + 2 * k] = dr2;
    u[NTIME * NSPACE * 2 + 2 * k + 1] = di2;

    //a = sqrt(dr1 * dr1 + di1 * di1);
    //printf("%3d / %18.15g + I * %18.15g / %18.15g\n", k, dr1, di1, a);
    //a = sqrt(dr2 * dr2 + di2 * di2);
    //printf("%3d / %18.15g + I * %18.15g / %18.15g\n", k, dr2, di2, a);
  } //for k

  fclose(fin);
}
