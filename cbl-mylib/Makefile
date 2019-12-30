######################################################################
# Makefile V.1.0
# Tue Sep 20 11:11:26 MDT 1994
# (cbl) - modified for 64 bits (ivh, 28 Jan 97)
######################################################################
# for production of mylib64.a
######################################################################
CC = cc
CFLAGS = -mips4 -64
F77 =  /usr/bin/f77
FFLAGS = -O3 -mips4 -64

.f:
	$(F77) $(FFLAGS) -c $<

.c:
	$(CC) $(CFLAGS) -c $<

mylib:	cbl_cputime cbl_datim cbl_header cbl_len cbl_open\
	cbl_random cbl_ranf cbl_rcarry cbl_rtim gettimeofday_
	ar rc mylib64.a *.o
	rm -f *.o
	echo Now move mylib64.a to the final directory

all:	mylib

clean:
	rm -f *.x *.o *.bak
