######################################################################
# Makefile for Linux
# modified by ihip 2019-12-30
######################################################################
# for production of mylib.a
######################################################################
CC = gcc
CFLAGS =
F77 = gfortran
FFLAGS = -cpp

.f:
	$(F77) $(FFLAGS) -c $<

.c:
	$(CC) $(CFLAGS) -c $<

mylib:	cbl_cputime cbl_datim cbl_header cbl_len cbl_open\
	cbl_random cbl_ranf cbl_rcarry cbl_rgauss cbl_rtim\
	gettimeofday_ gethostid_
	ar rc mylib.a *.o
	rm -f *.o
	@echo Now move mylib.a to the final directory

all:	mylib

clean:
	rm -f *.x *.o *.bak
