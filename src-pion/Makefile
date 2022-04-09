########################################################################
# Hip, 17 Jun 08; last modified: 2022-03-24 (x86_64 Ubuntu 20.04 LTS)  #
########################################################################

NTIME = 20
NSPACE = 20

########################################################################
# Linux version:
FC = gfortran
CC = cc
OPTIMIZATION = -O
LIB = /usr/lib/x86_64-linux-gnu/
LAPACK = $(LIB)liblapack.so $(LIB)libblas.so
########################################################################

FFLAGS = -c $(OPTIMIZATION) -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)
CFLAGS = -c -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)

inv: main_inv.o hf.o itopch.o scalp.o readconf.o pbp.o imassp.o\
     neuberger.o invert.o wdet.o Makefile
	$(FC) $(OPTIMIZATION) main_inv.o hf.o itopch.o scalp.o\
          readconf.o neuberger.o invert.o pbp.o imassp.o wdet.o\
          $(MYLIB) $(LAPACK) -o dsp_L$(NSPACE)x$(NTIME).x

main_inv.o: main_inv.F Makefile
	$(FC) $(FFLAGS) main_inv.F

invert.o: invert.F Makefile
	$(FC) $(FFLAGS) invert.F

pbp.o: pbp.F Makefile
	$(FC) $(FFLAGS) pbp.F

imassp.o: imassp.F Makefile
	$(FC) $(FFLAGS) imassp.F

hf.o: hf.F Makefile
	$(FC) $(FFLAGS) hf.F

itopch.o: itopch.F Makefile
	$(FC) $(FFLAGS) itopch.F

scalp.o: scalp.F Makefile
	$(FC) $(FFLAGS) scalp.F

neuberger.o: neuberger.F Makefile
	$(FC) $(FFLAGS) neuberger.F

readconf.o: readconf.F Makefile
	$(FC) $(FFLAGS) readconf.F

#readconf.o: readconf.c Makefile
#	$(CC) $(CFLAGS) readconf.c

wdet.o: wdet.F Makefile
	$(FC) $(FFLAGS) wdet.F

########################################################################

clean:
	rm *.o

########################################################################
