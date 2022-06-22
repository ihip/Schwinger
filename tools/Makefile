########################################################################
# Makefile for tools #### ivh, 21 Mar 97 - Last modified: 14 Jul 21 ####
########################################################################
#
# maximal lattice size in time dimension (recommended 64):
MAX_NTIME = 64
#
# maximal lattice size in space dimension (recommended 64):
MAX_NSPACE = 64
#
# maximal number of computables (recommended 30):
MAX_NCOMP = 30
#
# maximal number of mass computables (recommended 20):
MAX_MCOMP = 20
#
# maximal number of measurements (recommended 10000):
MAX_NMEAS = 10000
#
########################################################################
#
F77 = gfortran
FFLAGS = -cpp -DMAX_NCOMP=$(MAX_NCOMP) -DMAX_MCOMP=$(MAX_MCOMP) \
         -DMAX_NMEAS=$(MAX_NMEAS) -DMAX_NTIME=$(MAX_NTIME) \
	 -DMAX_NSPACE=$(MAX_NTIME)
#
########################################################################

datan: datan.f load.o jackknife.o jackknife_s.o
	$(F77) $(FFLAGS) datan.f load.o jackknife.o jackknife_s.o -o datan.x

########################################################################

massan: massan.f load.o jackknife.o bisec.o mjack.o weinfit.o mpjack.o
	$(F77) $(FFLAGS) massan.f load.o jackknife.o bisec.o mpjack.o\
	  mjack.o weinfit.o -o massan.x

########################################################################

m: m.f load.o jackknife.o bisec.o mjack.o weinfit.o
	$(F77) $(FFLAGS) m.f load.o jackknife.o bisec.o \
	  mjack.o weinfit.o -o m.x

########################################################################

wavan: wavan.f load.o jackknife.o
	$(F77) $(FFLAGS) wavan.f load.o jackknife.o -o wavan.x

########################################################################

collect: collect.f load.o jackknife.o jackknife_s.o
	$(F77) $(FFLAGS) collect.f load.o jackknife.o jackknife_s.o \
	  -o collect.x

########################################################################

gcollect: gcollect.f load.o jackknife.o jackknife_s.o
	$(F77) $(FFLAGS) gcollect.f load.o jackknife.o jackknife_s.o \
	  -o gcollect.x

########################################################################

pcaccoll: pcaccoll.f load.o jackknife.o
	$(F77) $(FFLAGS) pcaccoll.f load.o jackknife.o -o pcaccoll.x

########################################################################

pcaccollq: pcaccollq.f load.o jackknife.o
	$(F77) $(FFLAGS) pcaccollq.f load.o jackknife.o -o pcaccollq.x

########################################################################

masscoll: masscoll.f load.o jackknife.o mjack.o weinfit.o bisec.o
	$(F77) $(FFLAGS) masscoll.f load.o jackknife.o mjack.o \
	  weinfit.o bisec.o -o masscoll.x

########################################################################

maspcoll: maspcoll.f load.o djackknife.o mpjack.o weinfit.o bisec.o
	$(F77) $(FFLAGS) maspcoll.f load.o djackknife.o mpjack.o \
	  weinfit.o bisec.o -o maspcoll.x
	make clean

########################################################################

checksigma: checksigma.f load.o djackknife.o mpjack.o weinfit.o bisec.o
	$(F77) $(FFLAGS) checksigma.f load.o djackknife.o mpjack.o \
	  weinfit.o bisec.o -o checksigma.x
	  
########################################################################

nkcoll: nkcoll.f icjack.o load.o jackknife.o
	$(F77) $(FFLAGS) nkcoll.f icjack.o load.o jackknife.o -o nkcoll.x

########################################################################

nkhcoll: nkhcoll.f load.o jackknife.o
	$(F77) $(FFLAGS) nkhcoll.f load.o jackknife.o -o nkhcoll.x

########################################################################

nrealchk: nrealchk.f load.o jackknife.o
	$(F77) $(FFLAGS) nrealchk.f jackknife.o load.o -o nrealchk.x

########################################################################

jktop: jktop.f load.o jackknife.o
	$(F77) $(FFLAGS) jktop.f jackknife.o load.o -o jktop.x

########################################################################

pcacfit: pcacfit.f fit.o
	$(F77) $(FFLAGS) pcacfit.f fit.o -o pcacfit.x

########################################################################

pcacfitx: pcacfitx.f fitexy.o
	$(F77) $(FFLAGS) pcacfitx.f fitexy.o -o pcacfitx.x

########################################################################

join_data: join_data.f load.o
	$(F77) $(FFLAGS) join_data.f load.o -o join_data.x

########################################################################

join_masp: join_masp.f load.o
	$(F77) $(FFLAGS) join_masp.f load.o -o join_masp.x

########################################################################

join_neig: join_neig.f load.o
	$(F77) $(FFLAGS) join_neig.f -o join_neig.x

########################################################################

test_bisec: test_bisec.f bisec.o
	$(F77) $(FFLAGS) test_bisec.f bisec.o -o test_bisec.x

########################################################################

retfit: retfit.f fit.o
	$(F77) $(FFLAGS) retfit.f fit.o -o retfit.x

########################################################################

jackknife.o: jackknife.f
	$(F77) -c $(FFLAGS) jackknife.f

djackknife.o: djackknife.f
	$(F77) -c $(FFLAGS) djackknife.f

jackknife_s.o: jackknife_s.f
	$(F77) -c $(FFLAGS) jackknife_s.f

load.o: load.f
	$(F77) -c $(FFLAGS) load.f

fit.o: fit.f
	$(F77) -c $(FFLAGS) fit.f

fitexy.o: fitexy.f
	$(F77) -c $(FFLAGS) fitexy.f

bisec.o: bisec.f
	$(F77) -c $(FFLAGS) bisec.f

mjack.o: mjack.f
	$(F77) -c $(FFLAGS) mjack.f

mpjack.o: mpjack.f
	$(F77) -c $(FFLAGS) mpjack.f

icjack.o: icjack.f
	$(F77) -c $(FFLAGS) icjack.f

weinfit.o: weinfit.f
	$(F77) -c $(FFLAGS) weinfit.f

########################################################################

clean:
	rm *.o

########################################################################
