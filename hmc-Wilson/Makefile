########################################################################
# Makefile for hmc ##### ihip, 28 Jul 20 - Last modified: 10 Aug 20 ####
########################################################################
# You can edit the following parameters according to your needs:
#
# currently available integration routines: fullstep
INTEGRATE = fullstep
#
# currently available inversion routines for M^+ M:
# - cgr
INVERT = cgr
# - cgr_b (WARNING: cgr_b is identical with cgr!)
INVERTB = cgr_b
#
# currently available inversion routines for M:
# - cgf1, bicg1, bicgg5_1
INVERT1 = cgf1
# - bicg1_b (WARNING: bicg1_b is identical with bicg1!)
INVERT1B = bicg1_b
#
# inversion accuracy
ACCURACY = 0.2e-07
ACCURACY1 = 0.2e-07
#
# lattice size:
# it should be carefully checked when NTIME != NSPACE
NTIME = 16
NSPACE = 16
#
########################################################################
F77 = gfortran
OPTIMIZATION = -O
FFLAGS = -c -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)\
         -DACCURACY=$(ACCURACY) -DACCURACY1=$(ACCURACY1)
MYLIB = ../cbl-mylib/mylib.a
########################################################################

hmc: main.o hmc.o hmc_update.o create_u.o scalp.o norm.o unitarize.o\
	  mk_index.o sgauge.o itopch.o $(INTEGRATE).o\
	  $(INVERT).o $(INVERT1).o $(INVERTB).o $(INVERT1B).o\
	  inv.o save.o invert_excp.o mvec.o dsferm.o Makefile
	  
	$(F77) $(OPTIMIZATION)\
	  main.o hmc.o hmc_update.o create_u.o scalp.o norm.o unitarize.o\
	  mk_index.o sgauge.o itopch.o $(INTEGRATE).o\
	  $(INVERT).o $(INVERT1).o $(INVERTB).o $(INVERT1B).o\
	  inv.o save.o invert_excp.o mvec.o dsferm.o\
	  $(MYLIB) -o hmc$(NSPACE)x$(NTIME).x

########################################################################
	  
main.o: main.F Makefile
	$(F77) $(FFLAGS) main.F
	
hmc.o: hmc.F Makefile
	$(F77) $(FFLAGS) hmc.F

hmc_update.o: hmc_update.F Makefile
	$(F77) $(FFLAGS) hmc_update.F

fullstep.o: fullstep.F Makefile
	$(F77) $(FFLAGS) fullstep.F
 
create_u.o: create_u.F Makefile
	$(F77) $(FFLAGS) create_u.F
 
inv.o: inv.F Makefile
	$(F77) $(FFLAGS) inv.F

save.o: save.f Makefile
	$(F77) $(FFLAGS) save.f
 
invert_excp.o: invert_excp.F Makefile
	$(F77) $(FFLAGS) invert_excp.F

### standard functions #################################################

scalp.o: scalp.f Makefile
	$(F77) $(FFLAGS) scalp.f

norm.o: norm.f Makefile
	$(F77) $(FFLAGS) norm.f

### "standard" subroutines #############################################

unitarize.o: unitarize.f Makefile
	$(F77) $(FFLAGS) unitarize.f

mk_index.o: mk_index.F Makefile
	$(F77) $(FFLAGS) mk_index.F

sgauge.o: sgauge.F Makefile
	$(F77) $(FFLAGS) sgauge.F

itopch.o: itopch.F Makefile
	$(F77) $(FFLAGS) itopch.F

### fermion specific subroutines #######################################

mvec.o: mvec.F Makefile
	$(F77) $(FFLAGS) mvec.F

dsferm.o: dsferm.F Makefile
	$(F77) $(FFLAGS) dsferm.F

### invert subroutines #################################################

cgr.o: cgr.F Makefile
	$(F77) $(FFLAGS) cgr.F

cgr_b.o: cgr_b.F Makefile
	$(F77) $(FFLAGS) cgr_b.F

### invert1 subroutines ################################################

cgf1.o: cgf1.F Makefile
	$(F77) $(FFLAGS) cgf1.F

bicg1.o: bicg1.F Makefile
	$(F77) $(FFLAGS) bicg1.F

bicg1_b.o: bicg1_b.F Makefile
	$(F77) $(FFLAGS) bicg1_b.F

bicgg5_1.o: bicgg5_1.F Makefile
	$(F77) $(FFLAGS) bicgg5_1.F

########################################################################

clean:
	rm *.o

########################################################################
