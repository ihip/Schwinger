########################################################################
# ihip, 24 Nov 07; last modified: 10 Aug 20 (Linux - Ubuntu)           #
########################################################################

NTIME = 16
NSPACE = 16

########################################################################
# Linux version:
FC = gfortran
CC = gcc
OPTIMIZATION = -O
########################################################################

FFLAGS = -c $(OPTIMIZATION) -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)
CFLAGS = -c -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)

main_c2mv: main_c2mv.o mk_index.o mvec.o mvec_Wilson.o hf.o Makefile

	$(FC) $(OPTIMIZATION) main_c2mv.o mk_index.o mvec.o mvec_Wilson.o\
	  hf.o -o c2mv$(NSPACE)x$(NTIME).x

main_eig.o: main_c2mv.F Makefile
	$(FC) $(FFLAGS) main_c2mv.F

mk_index.o: mk_index.F Makefile
	$(FC) $(FFLAGS) mk_index.F

mvec.o: mvec.F Makefile
	$(FC) $(FFLAGS) mvec.F
	
mvec_Wilson.o: mvec_Wilson.F Makefile
	$(FC) $(FFLAGS) mvec_Wilson.F
	
hf.o: hf.F Makefile
	$(FC) $(FFLAGS) hf.F

########################################################################

clean:
	rm *.o

########################################################################
