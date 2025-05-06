#include <stdio.h>
#include <stdlib.h>

int main() {
    char input_file_name[64];
    FILE *infile;
    int nread, usize;

    int marker;
    int ntime, nspace;
    double beta, akap, eps;
    int iter;
    double *u;

    int it, ix, is;
    double u1real, u1imag, u2real, u2imag;

    printf("Input file name: ");
    scanf("%s", input_file_name);
    printf("\n");

// open binary file to read
    infile = fopen(input_file_name, "rb");

// every record starts with a marker: 4-byte integer with the length of the record

// first record in FORTRAN written as: write(7) NTIME, NSPACE
    nread = fread(&marker, sizeof(int), 1, infile); // marker should be 8 (two 4-byte integers) 
    nread = fread(&ntime, sizeof(int), 1, infile);  // 4-byte integer ntime
    printf("ntime =  %d\n", ntime);
    nread = fread(&nspace, sizeof(int), 1, infile); // 4-byte integer nspace
    printf("nspace = %d\n", nspace);
    nread = fread(&marker, sizeof(int), 1, infile); // closing marker should be 8

// second record: write(7) beta, akap, eps
    nread = fread(&marker, sizeof(int), 1, infile); // marker should be 24 (three 8-byte reals) 
    nread = fread(&beta, sizeof(double), 1, infile);
    printf("beta = %20.15lf\n", beta);
    nread = fread(&akap, sizeof(double), 1, infile);
    printf("akap = %20.15f\n", akap);
    nread = fread(&eps, sizeof(double), 1, infile);
    printf("eps = %20.15f\n", eps);
    nread = fread(&marker, sizeof(int), 1, infile); // closing marker should be 24

// third record: write(7) ntherm + it * mstep
    nread = fread(&marker, sizeof(int), 1, infile); // marker should be 4 (4-byte integer) 
    nread = fread(&iter, sizeof(int), 1, infile);
    printf("iter = %d\n", iter);
    nread = fread(&marker, sizeof(int), 1, infile); // closing marker should be 4

// allocate memory for gauge fields depending on lattice size
    usize = ntime * nspace * 4 * sizeof(double);
    u = (double *)malloc(usize);

// fourth record contains gauge fields: write(7) u
    nread = fread(&marker, sizeof(int), 1, infile);  // marker should be usize
    printf("\ngauge field marker: expected = %d > found = %d\n\n", usize, marker);
    nread = fread(u, sizeof(double), ntime * nspace * 4, infile);
    nread = fread(&marker, sizeof(int), 1, infile);  // closing marker should be usize

// close input file
    fclose(infile);

// write gauge fields from array u
// in the memory there is a block of ntime * nspace gauge fields with index 1: u(is, 1)
// and then follows another block of ntime * nspace gauge fields with index 2: u(is, 2)
    for(it = 1; it <= ntime; it++)
      for(ix = 1; ix <= nspace; ix++) {
        is = (it - 1) * nspace + (ix - 1);
        u1real = u[is * 2];
        u1imag = u[is * 2 + 1];
        printf("%8d %8d %22.17lf %22.17lf\n", it, ix, u1real, u1imag);
        u2real = u[(ntime * nspace + is) * 2];
        u2imag = u[(ntime * nspace + is) * 2 + 1];
        printf("%8d %8d %22.17lf %22.17lf\n", it, ix, u2real, u2imag);    
      }

} // main