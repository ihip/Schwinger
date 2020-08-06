########################################################################
# Hip, 24 Nov 07; last modified: 03 Aug 20 (Linux - Ubuntu)            #
########################################################################

NTIME = 16
NSPACE = 16

########################################################################
# Linux version:
FC = gfortran
CC = gcc
OPTIMIZATION = -O
LIB = /usr/lib/x86_64-linux-gnu/
LAPACK = $(LIB)lapack/liblapack.so.3 $(LIB)blas/libblas.so.3
########################################################################

FFLAGS = -c $(OPTIMIZATION) -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)
CFLAGS = -c -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)

eig: main_eig.o hf.o itopch.o lapackeig.o scalp.o\
     overlap.o save_eig.o Makefile
	$(FC) $(OPTIMIZATION) main_eig.o hf.o itopch.o scalp.o\
          lapackeig.o overlap.o save_eig.o\
          $(MYLIB) $(LAPACK) -o eig$(NTIME).x

main_eig.o: main_eig.F Makefile
	$(FC) $(FFLAGS) main_eig.F

hf.o: hf.F Makefile
	$(FC) $(FFLAGS) hf.F

itopch.o: itopch.F Makefile
	$(FC) $(FFLAGS) itopch.F

lapackeig.o: lapackeig.F Makefile
	$(FC) $(FFLAGS) lapackeig.F

scalp.o: scalp.f Makefile
	$(FC) $(FFLAGS) scalp.f

overlap.o: overlap.F Makefile
	$(FC) $(FFLAGS) overlap.F

save_eig.o: save_eig.f Makefile
	$(FC) $(FFLAGS) save_eig.f

########################################################################

clean:
	rm *.o

########################################################################
