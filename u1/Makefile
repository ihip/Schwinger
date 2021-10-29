########################################################################
# Makefile for metropolis ##### Hip, 2021-07-12                    #####
########################################################################
# You can edit the following parameters according to your needs:
#
# lattice size:
# [Hip: for metropolis it is not clear if NTIME != NSPACE is ok...]
NTIME = 8
NSPACE = 64
#
########################################################################
F77 = gfortran
OPTIMIZATION =
FFLAGS = -c -cpp -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)
MYLIB = ../cbl-mylib/mylib.a
SRC = ../src
VPATH = $(SRC)
########################################################################

metropolis: main_met.o metropolis.o save.o mk_index.o create_u.o\
		sgauge.o itopch.o sintopch.o v_dummy.o Makefile
	$(F77) $(OPTIMIZATION) main_met.o metropolis.o mk_index.o\
		create_u.o sgauge.o save.o itopch.o sintopch.o v_dummy.o\
		$(MYLIB) -o met_L$(NSPACE)x$(NTIME).x

main_met.o: main_met.F Makefile
	$(F77) $(FFLAGS) main_met.F

metropolis.o: metropolis.F Makefile
	$(F77) $(FFLAGS) metropolis.F

create_u.o: create_u.F Makefile
	$(F77) $(FFLAGS) create_u.F

save.o: save.f Makefile
	$(F77) $(FFLAGS) save.f

v_dummy.o: v_dummy.f Makefile
	$(F77) $(FFLAGS) v_dummy.f

### subroutines from ../src ############################################

mk_index.o: mk_index.F Makefile
	$(F77) $(FFLAGS) $(SRC)/mk_index.F

sgauge.o: sgauge.F Makefile
	$(F77) $(FFLAGS) $(SRC)/sgauge.F

itopch.o: itopch.F Makefile
	$(F77) $(FFLAGS) $(SRC)/itopch.F

sintopch.o: sintopch.F Makefile
	$(F77) $(FFLAGS) $(SRC)/sintopch.F

########################################################################

clean:
	rm *.o

########################################################################
