NTIME = 16
NSPACE = 16

F77 = gfortran
FFLAGS = -cpp -DNTIME=$(NTIME) -DNSPACE=$(NSPACE)

CC = gcc

readconf: readconf.F
	$(F77) $(FFLAGS) readconf.F -o readconf.x

creadconf: creadconf.c
	$(CC) creadconf.c -o creadconf.x

clean:
	rm *.o

