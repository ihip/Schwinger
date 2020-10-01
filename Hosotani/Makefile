########################################################################
# ihip, 01 Oct 20; last modified: 01 Oct 20 (Linux - Ubuntu)           #
########################################################################
FC = gfortran
OPTIMIZATION = -O
LIB = /usr/lib/x86_64-linux-gnu/
LAPACK = $(LIB)lapack/liblapack.so.3 $(LIB)blas/libblas.so.3
########################################################################

FFLAGS = -c $(OPTIMIZATION)

hop: hop.o Makefile
	$(FC) $(OPTIMIZATION) hop.o $(LAPACK) -o hop.x

hop.o: hop.F Makefile
	$(FC) $(FFLAGS) hop.F

Hosotani: Hosotani.o Makefile
	$(FC) $(OPTIMIZATION) Hosotani.o $(LAPACK) -o Hosotani.x

Hosotani.o: Hosotani.F Makefile
	$(FC) $(FFLAGS) Hosotani.F

########################################################################

clean:
	rm *.o

########################################################################
