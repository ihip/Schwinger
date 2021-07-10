######################################################################
# Makefile for Linux
# modified by ihip 2020-07-15 / 2021-07-10
######################################################################

F77 = gfortran
FFLAGS = -cpp

all:
	$(F77) $(FFLAGS) gauge_met_c.f gauge_met_c_chiT.F\
	  ../../cbl-mylib/mylib.a -o u1_chiT.x

clean:
	rm -f *.x *.o
