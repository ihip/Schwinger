######################################################################
# Makefile for Linux
# modified by ihip 2019-12-30
######################################################################

F77 = gfortran
FFLAGS = -cpp

.f:
	$(F77) $(FFLAGS) $< -o $@.x ../mylib.a

all: test_cputime test_datim test_gethostid test_idate test_ranf\
	test_rcarry test_rgauss_new test_rtim test_gethostname\
	test_gettimeofday

clean:
	rm -f *.x *.o
