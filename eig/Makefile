########################################################################
# ihip, 24 Nov 07; last modified: 2023-01-07                           #
########################################################################

NTIME =  20
NSPACE = 20

########################################################################
# Linux version:
FC = gfortran
CC = gcc
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
########################################################################

FFLAGS = -c $(OPTIMIZATION) -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)
CFLAGS = -c -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)

eig: main_eig.o hf.o mk_index.o itopch.o lapackeig.o scalp.o\
	 normalize.o overlap.o save_eig.o Makefile
	$(FC) $(OPTIMIZATION) main_eig.o hf.o mk_index.o itopch.o\
           scalp.o normalize.o lapackeig.o overlap.o save_eig.o\
          $(MYLIB) $(LAPACK) -o eig$(NSPACE)x$(NTIME).x

# eigenvalues of hypercube overlap fermions (HO)
cfeig: main_cfeig.o hf.o mk_index.o itopch.o lapackeig.o scalp.o\
	 normalize.o overlap.o save_eig.o Makefile
	$(FC) $(OPTIMIZATION) main_cfeig.o hf.o mk_index.o itopch.o\
           scalp.o normalize.o lapackeig.o overlap.o save_eig.o\
          $(MYLIB) $(LAPACK) -o cfeig$(NSPACE)x$(NTIME).x

# eigenvalues of central-branch Wilson fermions (cbWf)
cbWf_cfeig: main_cbWf_cfeig.o mvec_cbWf.o mk_index.o itopch.o lapackeig.o scalp.o\
	 normalize.o save_eig.o Makefile
	$(FC) $(OPTIMIZATION) main_cbWf_cfeig.o mvec_cbWf.o mk_index.o itopch.o\
           scalp.o normalize.o lapackeig.o save_eig.o\
          $(MYLIB) $(LAPACK) -o cbWf_cfeig$(NSPACE)x$(NTIME).x
		  
sigma: sigma.o Makefile
	$(FC) $(OPTIMIZATION) sigma.o -o sigma.x	

sigma.o: sigma.F Makefile
	$(FC) $(FFLAGS) sigma.F

main_eig.o: main_eig.F Makefile
	$(FC) $(FFLAGS) main_eig.F

main_cfeig.o: main_cfeig.F Makefile
	$(FC) $(FFLAGS) main_cfeig.F

main_cbWf_cfeig.o: main_cbWf_cfeig.F Makefile
	$(FC) $(FFLAGS) main_cbWf_cfeig.F

mvec_cbWf.o: mvec_cbWf.F Makefile
	$(FC) $(FFLAGS) mvec_cbWf.F

hf.o: hf.F Makefile
	$(FC) $(FFLAGS) hf.F

mk_index.o: mk_index.F Makefile
	$(FC) $(FFLAGS) mk_index.F
	
itopch.o: itopch.F Makefile
	$(FC) $(FFLAGS) itopch.F

lapackeig.o: lapackeig.F Makefile
	$(FC) $(FFLAGS) lapackeig.F

scalp.o: scalp.f Makefile
	$(FC) $(FFLAGS) scalp.f
	
normalize.o: normalize.f Makefile
	$(FC) $(FFLAGS) normalize.f

overlap.o: overlap.F Makefile
	$(FC) $(FFLAGS) overlap.F

save_eig.o: save_eig.f Makefile
	$(FC) $(FFLAGS) save_eig.f

########################################################################

clean:
	rm *.o

########################################################################
