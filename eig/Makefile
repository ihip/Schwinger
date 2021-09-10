########################################################################
# ihip, 24 Nov 07; last modified: 10 Aug 20 (Linux - Ubuntu)           #
########################################################################

NTIME = 20
NSPACE = 20

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

cfeig: main_cfeig.o hf.o mk_index.o itopch.o lapackeig.o scalp.o\
	 normalize.o overlap.o save_eig.o Makefile
	$(FC) $(OPTIMIZATION) main_cfeig.o hf.o mk_index.o itopch.o\
           scalp.o normalize.o lapackeig.o overlap.o save_eig.o\
          $(MYLIB) $(LAPACK) -o cfeig$(NSPACE)x$(NTIME).x
		  
eig: main_eig.o hf.o mk_index.o itopch.o lapackeig.o scalp.o\
	 normalize.o overlap.o save_eig.o Makefile
	$(FC) $(OPTIMIZATION) main_eig.o hf.o mk_index.o itopch.o\
           scalp.o normalize.o lapackeig.o overlap.o save_eig.o\
          $(MYLIB) $(LAPACK) -o eig$(NSPACE)x$(NTIME).x

main_eig.o: main_eig.F Makefile
	$(FC) $(FFLAGS) main_eig.F

main_cfeig.o: main_cfeig.F Makefile
	$(FC) $(FFLAGS) main_cfeig.F

hf.o: hf.F Makefile
	$(FC) $(FFLAGS) hf.F

mk_index.o: mk_index.F Makefile
	$(FC) $(FFLAGS) mk_index.F
	
itopch.o: itopch.F Makefile
	$(FC) $(FFLAGS) itopch.F

lapackeig.o: lapackeig.F Makefile
	$(FC) $(FFLAGS) lapackeig.F

scalp.o: scalp.f Makefile
	$(FC) $(FFLAGS) scalp.f
	
normalize.o: normalize.f Makefile
	$(FC) $(FFLAGS) normalize.f

overlap.o: overlap.F Makefile
	$(FC) $(FFLAGS) overlap.F

save_eig.o: save_eig.f Makefile
	$(FC) $(FFLAGS) save_eig.f

########################################################################

clean:
	rm *.o

########################################################################
