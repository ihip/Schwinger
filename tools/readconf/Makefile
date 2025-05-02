NTIME = 16
NSPACE = 16

F77 = gfortran
FFLAGS = -cpp -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)

readconf: readconf.F
	$(F77) $(FFLAGS) readconf.F -o readconf.x

clean:
	rm *.o

