########################################################################
# Makefile for metropolis ##### Hip, 2021-07-12                    #####
########################################################################
# You can edit the following parameters according to your needs:
#
# lattice size:
# [Hip: for metropolis it is not clear if NTIME != NSPACE is ok...]
NTIME = 16
NSPACE = 16
#
########################################################################
F77 = gfortran
OPTIMIZATION =
FFLAGS = -c -cpp -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)
MYLIB = ../cbl-mylib/mylib.a
########################################################################

metropolis: main_met.o metropolis.o save.o hmc.o itopch.o sintopch.o\
        v_dummy.o Makefile
	$(F77) $(OPTIMIZATION) main_met.o metropolis.o save.o hmc.o\
	    itopch.o sintopch.o v_dummy.o $(MYLIB) -o met$(NTIME).x

main_met.o: main_met.f Makefile
	$(F77) $(FFLAGS) main_met.f

metropolis.o: metropolis.f Makefile
	$(F77) $(FFLAGS) metropolis.f

save.o: save.f Makefile
	$(F77) $(FFLAGS) save.f

hmc.o: hmc.f Makefile
	$(F77) $(FFLAGS) hmc.f

itopch.o: itopch.f Makefile
	$(F77) $(FFLAGS) itopch.f

sintopch.o: sintopch.f Makefile
	$(F77) $(FFLAGS) sintopch.f

v_dummy.o: v_dummy.f Makefile
	$(F77) $(FFLAGS) v_dummy.f

########################################################################

clean:
	rm *.o

########################################################################
