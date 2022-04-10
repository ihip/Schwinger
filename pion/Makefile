########################################################################
# Hip, 17 Jun 08; last modified: 2022-04-10 (x86_64 Ubuntu 20.04 LTS)  #
########################################################################

NTIME = 20
NSPACE = 20

########################################################################
# Linux version:
FC = gfortran
CC = cc
OPTIMIZATION = -O
# WSL (Windows Subsystem for Linux) Ubuntu
ifeq ($(WSL_DISTRO_NAME),Ubuntu)
	LIB = /usr/lib/x86_64-linux-gnu/
endif
# cluster Isabella (SRCE Zagreb)
ifeq ($(SGE_CLUSTER_NAME),isabella)
	LIB = /usr/lib64/
endif
LAPACK = $(LIB)liblapack.so.3
########################################################################

FFLAGS = -c $(OPTIMIZATION) -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)
CFLAGS = -c -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)

dsp: main_dsp.o hf.o itopch.o scalp.o readconf.o pbp.o imassp.o\
     neuberger.o invert.o wdet.o Makefile
	$(FC) $(OPTIMIZATION) main_dsp.o hf.o itopch.o scalp.o\
          readconf.o neuberger.o invert.o pbp.o imassp.o wdet.o\
          $(MYLIB) $(LAPACK) -o dsp_L$(NSPACE)x$(NTIME).x
	rm *.o

main_dsp.o: main_dsp.F Makefile
	$(FC) $(FFLAGS) main_dsp.F

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
