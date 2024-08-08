########################################################################
# Hip, 17 Jun 08; last modified: 2023-01-08 (x86_64 Ubuntu 20.04 LTS)  #
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
# cluster at UNAM (Mexico City) - 
ifeq ($(SITE_NAME),ICN-UNAM)
	LIB = /software/icn/lapack/3.9/
endif
LAPACK = $(LIB)liblapack.so.3
MYLIB = ../../cbl-mylib/mylib.a
########################################################################

FFLAGS = -c $(OPTIMIZATION) -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)
CFLAGS = -c -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)

dsp_cbWf: main_dsp_cbWf.o itopch.o scalp.o readconf.o pbp.o imassp.o\
     invert.o wdet.o save_old_header.o mvec_cbWf.o Makefile
	$(FC) $(OPTIMIZATION) main_dsp_cbWf.o itopch.o scalp.o\
          readconf.o invert.o pbp.o imassp.o\
		  wdet.o save_old_header.o mvec_cbWf.o\
          $(MYLIB) $(LAPACK) -o dsp_cbWf_L$(NSPACE)x$(NTIME).x
	rm *.o

main_dsp_cbWf.o: main_dsp_cbWf.F Makefile
	$(FC) $(FFLAGS) main_dsp_cbWf.F

mvec_cbWf.o: mvec_cbWf.F Makefile
	$(FC) $(FFLAGS) mvec_cbWf.F

invert.o: ../invert.F Makefile
	$(FC) $(FFLAGS) ../invert.F

pbp.o: ../pbp.F Makefile
	$(FC) $(FFLAGS) ../pbp.F

imassp.o: ../imassp.F Makefile
	$(FC) $(FFLAGS) ../imassp.F

itopch.o: ../itopch.F Makefile
	$(FC) $(FFLAGS) ../itopch.F

scalp.o: ../scalp.F Makefile
	$(FC) $(FFLAGS) ../scalp.F

readconf.o: ../readconf.F Makefile
	$(FC) $(FFLAGS) ../readconf.F

#readconf.o: readconf.c Makefile
#	$(CC) $(CFLAGS) readconf.c

wdet.o: wdet.F Makefile
	$(FC) $(FFLAGS) wdet.F

save_old_header.o: ../save_old_header.F Makefile
	$(FC) $(FFLAGS) ../save_old_header.F

########################################################################

clean:
	rm *.o

########################################################################
