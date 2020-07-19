########################################################################
# Makefile for hmc ##### hip, 16 Jul 20 - Last modified: 16 Jul 20 #####
########################################################################
# You can edit the following parameters according to your needs:
#
# currently available integration routines: fullstep
INTEGRATE = fullstep
#
# currently available inversion routines for M^+ M:
# - cgr(b), eo_bicgg5, eo_bicgs, eo_bicgs2
INVERT = eo_bicgg5
INVERTB = cgr_b
#
# currently available inversion routines for M:
# - bicg1(b), bicgg5_1, cgf1, eo1_bicgg5, eo1_bicgs, eo1_bicgs2
INVERT1 = eo1_bicgg5
INVERT1B = bicg1_b
#
# inversion accuracy
ACCURACY = 0.2e-07
ACCURACY1 = 0.2e-07
#
# lattice size:
# !if you use eo-preconditioning it should be NTIME = NSPACE!
NTIME = 16
NSPACE = 16
#
########################################################################
F77 = gfortran
OPTIMIZATION =
FFLAGS = -c -cpp -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)\
         -DACCURACY=$(ACCURACY) -DACCURACY1=$(ACCURACY1)
MYLIB = ../cbl-mylib/mylib.a
#LIB = /usr/local/lib/
#NAG = $(LIB)libnag.a $(LIB)libnagblas.a
#GL = -lfgl -lgl
########################################################################

hmc: main.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o $(INVERT1).o eo.o eo1.o\
	  itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile
	$(F77) $(OPTIMIZATION) main.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o\
	  $(INVERT1).o eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o\
	  op_ubudbd.o inv.o mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  $(MYLIB) -o hmc$(NTIME).x

hmc_pp: main_pp.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o $(INVERT1).o eo.o eo1.o\
	  itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile
	$(F77) $(OPTIMIZATION) main_pp.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o\
	  $(INVERT1).o eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o\
	  op_ubudbd.o inv.o mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  $(MYLIB) -o hmc_pp$(NTIME).x

# for check of dispersion relations:
hmc_dsp: main_dsp.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o $(INVERT1).o eo.o eo1.o\
	  itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o imassp.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile
	$(F77) $(OPTIMIZATION) main_dsp.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o\
	  $(INVERT1).o eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o imassp.o\
	  op_ubudbd.o inv.o mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  $(MYLIB) -o hmc_dsp$(NTIME).x

hmc_w: main_w.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o $(INVERT1).o eo.o eo1fix.o\
	  itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  wave.o fix.o Makefile
	$(F77) $(OPTIMIZATION) main_w.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o\
	  $(INVERT1).o eo.o eo1fix.o itopch.o eigen.o escalp.o init_gamma.o\
	  op_ubudbd.o inv.o mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  wave.o fix.o $(MYLIB) -o hmc_w$(NTIME).x

naghmc: main_nag.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o $(INVERT1).o eo.o eo1.o\
	  itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile \
	  nageig.o
	$(F77) $(OPTIMIZATION) main_nag.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o\
	  $(INVERT1).o eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o\
	  op_ubudbd.o inv.o mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  nageig.o $(NAG) $(MYLIB) -o nag_hmc$(NTIME).x

dhmc: main_d.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o $(INVERT1).o eo.o eo1.o\
	  itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile \
	  nageig.o
	$(F77) $(OPTIMIZATION) main_d.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o\
	  $(INVERT1).o eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o\
	  op_ubudbd.o inv.o mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  nageig.o $(NAG) $(MYLIB) -o dhmc$(NTIME).x

thmc: main_t.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o $(INVERT1).o eo.o eo1.o\
	  itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile \
	  nagevect.o
	$(F77) $(OPTIMIZATION) main_t.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o\
	  $(INVERT1).o eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o\
	  op_ubudbd.o inv.o mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  nagevect.o $(NAG) $(MYLIB) -o thmc$(NTIME).x

ghmc: main.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o $(INVERT1).o eo.o eo1.o\
	  itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile \
	  nageig.o graf.o
	$(F77) $(OPTIMIZATION) main.o hmc.o hmc_update.o $(INTEGRATE).o $(INVERT).o\
	  $(INVERT1).o eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o\
	  op_ubudbd.o inv.o mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  nageig.o graf.o $(NAG) $(GL) $(MYLIB) -o ghmc$(NTIME).x

quenched: main_q.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o eo.o\
	  eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile
	$(F77) $(OPTIMIZATION) main_q.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o\
	  eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
	  mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  $(MYLIB) -o q$(NTIME).x

metropolis: main_met.o metropolis.o save.o hmc.o itopch.o v_dummy.o Makefile
	$(F77) $(OPTIMIZATION) main_met.o metropolis.o save.o hmc.o itopch.o v_dummy.o\
          $(MYLIB) -o met$(NTIME).x

# quenched + eigenvalues in Hasenfratz projection
qhas: main_has.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o eo.o\
	  eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o\
          save.o invert_excp.o has.o $(INVERTB).o $(INVERT1B).o Makefile
	$(F77) $(OPTIMIZATION) main_has.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o\
	  eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o\
	  save.o invert_excp.o has.o $(INVERTB).o $(INVERT1B).o\
	  nageig.o $(NAG) $(MYLIB) -o qhas$(NTIME).x

# quenched + eigenvalues in Hasenfratz projection
qhas0: main_has0.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o eo.o\
	  eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o\
          save.o invert_excp.o has.o $(INVERTB).o $(INVERT1B).o Makefile
	$(F77) $(OPTIMIZATION) main_has0.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o\
	  eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o\
	  save.o invert_excp.o has.o $(INVERTB).o $(INVERT1B).o\
	  nageig.o $(NAG) $(MYLIB) -o qhas0$(NTIME).x

# preparing quenched gauge configurations
pqgc: main_pqgc.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o eo.o\
	  eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o\
          save.o invert_excp.o has.o $(INVERTB).o $(INVERT1B).o Makefile
	$(F77) $(OPTIMIZATION) main_pqgc.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o\
	  eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o\
	  save.o invert_excp.o has.o $(INVERTB).o $(INVERT1B).o\
	  nageig.o $(NAG) $(MYLIB) -o pqgc$(NTIME).x

# quenched mass & PCAC measurements (01 Jul 98):
# identical with normal main, but calls q_update (+ dependencies)
qm: main_qm.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o eo.o\
	  eo1.o itopch.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile \
	  nageig.o
	$(F77) $(OPTIMIZATION) main_qm.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o\
	  eo.o eo1.o itopch.o escalp.o init_gamma.o op_ubudbd.o inv.o\
	  save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  nageig.o $(NAG) $(MYLIB) -o qm$(NTIME).x

gquenched: main_q.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o eo.o\
	  eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
          mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o Makefile \
	  nageig.o graf.o
	$(F77) $(OPTIMIZATION) main_q.o hmc.o q_update.o quenched.o $(INVERT).o $(INVERT1).o\
	  eo.o eo1.o itopch.o eigen.o escalp.o init_gamma.o op_ubudbd.o inv.o\
	  mass.o save.o invert_excp.o $(INVERTB).o $(INVERT1B).o\
	  nageig.o graf.o $(NAG) $(GL) $(MYLIB) -o gq$(NTIME).x

main_met.o: main_met.f Makefile
	$(F77) $(FFLAGS) main_met.f

main_pqgc.o: main_pqgc.f Makefile
	$(F77) $(FFLAGS) main_pqgc.f

main_has.o: main_has.f Makefile
	$(F77) $(FFLAGS) main_has.f

main_has0.o: main_has0.f Makefile
	$(F77) $(FFLAGS) main_has0.f

has.o: has.f Makefile
	$(F77) $(FFLAGS) has.f

main_qm.o: main_qm.f Makefile
	$(F77) $(FFLAGS) main_qm.f

main.o: main.f Makefile
	$(F77) $(FFLAGS) main.f

main_w.o: main_w.f Makefile
	$(F77) $(FFLAGS) main_w.f

main_nag.o: main_nag.f Makefile
	$(F77) $(FFLAGS) main_nag.f

main_q.o: main_q.f Makefile
	$(F77) $(FFLAGS) main_q.f

main_t.o: main_t.f Makefile
	$(F77) $(FFLAGS) main_t.f

main_d.o: main_d.f Makefile
	$(F77) $(FFLAGS) main_d.f

main_dsp.o: main_dsp.f Makefile
	$(F77) $(FFLAGS) main_dsp.f

hmc.o: hmc.f Makefile
	$(F77) $(FFLAGS) hmc.f

hmc_update.o: hmc_update.f Makefile
	$(F77) $(FFLAGS) hmc_update.f

q_update.o: q_update.f Makefile
	$(F77) $(FFLAGS) q_update.f

eo.o: eo.f Makefile
	$(F77) $(FFLAGS) eo.f

eo1.o: eo1.f Makefile
	$(F77) $(FFLAGS) eo1.f

eo1fix.o: eo1fix.f Makefile
	$(F77) $(FFLAGS) eo1fix.f

fullstep.o: fullstep.f Makefile
	$(F77) $(FFLAGS) fullstep.f

quenched.o: quenched.f Makefile
	$(F77) $(FFLAGS) quenched.f

itopch.o: itopch.f Makefile
	$(F77) $(FFLAGS) itopch.f

eigen.o: eigen.f Makefile
	$(F77) $(FFLAGS) eigen.f

escalp.o: escalp.f Makefile
	$(F77) $(FFLAGS) escalp.f

init_gamma.o: init_gamma.f Makefile
	$(F77) $(FFLAGS) init_gamma.f

op_ubudbd.o: op_ubudbd.f Makefile
	$(F77) $(FFLAGS) op_ubudbd.f

inv.o: inv.f Makefile
	$(F77) $(FFLAGS) inv.f

wave.o: wave.f Makefile
	$(F77) $(FFLAGS) wave.f

fix.o: fix.f Makefile
	$(F77) $(FFLAGS) fix.f

mass.o: mass.f Makefile
	$(F77) $(FFLAGS) mass.f

save.o: save.f Makefile
	$(F77) $(FFLAGS) save.f

invert_excp.o: invert_excp.f Makefile
	$(F77) $(FFLAGS) invert_excp.f

nageig.o: nageig.f Makefile
	$(F77) $(FFLAGS) nageig.f

nagevect.o: nagevect.f Makefile
	$(F77) $(FFLAGS) nagevect.f

graf.o: graf.f
	$(F77) $(FFLAGS) graf.f

imassp.o: imassp.f Makefile
	$(F77) $(FFLAGS) imassp.f

metropolis.o: metropolis.f Makefile
	$(F77) $(FFLAGS) metropolis.f

v_dummy.o: v_dummy.f Makefile
	$(F77) $(FFLAGS) v_dummy.f

### invert functions ###################################################

cgr.o: cgr.f Makefile
	$(F77) $(FFLAGS) cgr.f

cgr_b.o: cgr_b.f Makefile
	$(F77) $(FFLAGS) cgr_b.f

eo_bicgg5.o: eo_bicgg5.f Makefile
	$(F77) $(FFLAGS) eo_bicgg5.f

eo_bicgs.o: eo_bicgs.f Makefile
	$(F77) $(FFLAGS) eo_bicgs.f

eo_bicgs2.o: eo_bicgs2.f Makefile
	$(F77) $(FFLAGS) eo_bicgs2.f

### invert1 functions ##################################################

bicg1.o: bicg1.f Makefile
	$(F77) $(FFLAGS) bicg1.f

bicg1_b.o: bicg1_b.f Makefile
	$(F77) $(FFLAGS) bicg1_b.f

bicgg5_1.o: bicgg5_1.f Makefile
	$(F77) $(FFLAGS) bicgg5_1.f

cgf1.o: cgf1.f Makefile
	$(F77) $(FFLAGS) cgf1.f

eo1_bicgg5.o: eo1_bicgg5.f Makefile
	$(F77) $(FFLAGS) eo1_bicgg5.f

eo1_bicgs.o: eo1_bicgs.f Makefile
	$(F77) $(FFLAGS) eo1_bicgs.f

eo1_bicgs2.o: eo1_bicgs2.f Makefile
	$(F77) $(FFLAGS) eo1_bicgs2.f

########################################################################

clean:
	rm *.o

########################################################################
