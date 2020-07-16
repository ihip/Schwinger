######################################################################
# Makefile for Linux
# modified by ihip 2020-07-15
######################################################################

F77 = gfortran
FFLAGS = -cpp

all:
	$(F77) $(FFLAGS) gauge_met_c.f gauge_met_c_test.f\
	  ../cbl-mylib/mylib.a -o u1_met.x

clean:
	rm -f *.x *.o
