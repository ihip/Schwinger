c	>>> sigma >>> Ivan Hip >>> 19.01.2008.	

# define MAX_LIST_LEN 10000

	character*64 listname, eigfile, logname
	real*8 re, im, m
	complex*16 lambda, sum, pro
	complex*16 sigma_av
	complex*16 sigma(MAX_LIST_LEN)
	real*8 det(MAX_LIST_LEN)
	real*8 dsum, det_av

	write(*, *)
	write(*, *) 'Sigma / 2021-08-27'
	write(*, *) '=================='
	write(*, *)

c       >>> read listname
	write(*, '(1x, ''List name of .eig files: '', $)')
	read(*, '(a)') listname

	open(1, file = listname, form = 'formatted', status = 'old')

c	>>> open log file
c	logname = listname(1:lnblnk(listname))//'.sig'
c	open(3, file = logname, form = 'formatted', status = 'unknown')

c	>>> read lattice size
	write(*, *)
	write(*, '(1x, ''npace = ntime = '', $)')
	read(*, *) nspace
	nsites = nspace * nspace
	neigen = nsites * 2
	write(*, *) '--> ', neigen, ' eigenvalues expected'
	write(*, *)

c	>>> read number of flavors
	write(*, '(1x, ''N_f = '', $)')
	read(*, *) nf
	write(*, *)

c	>>> read mass parameter
c	write(*, '(1x, ''mass: '', $)')
c	read(*, *) m
c	write(*, *)

c	  >>> loop over m
      m = 0.15d0

22    ifile = 0
c       >>> main loop
c11	read(1, '(a, i, i)', end = 99) eigfile, itop, izero
11    read(1, *, end = 99) eigfile
	  ifile = ifile + 1

c	eigfile = eigfile(1:lnblnk(eigfile))//'.eig'
c	write(*, '(a, $)') eigfile

	open(2, file = eigfile, form = 'formatted', status = 'old')

	sum = dcmplx(0.0, 0.0)
	pro = dcmplx(1.0d0, 0.0)
	do i = 1, neigen
	  read(2, *) re, im
c	  if(dabs(re) .lt. 1.0d-12) re = 0.0d0
c	  if(dabs(im) .lt. 1.0d-12) im = 0.0d0
	  lambda = re + dcmplx(0.0, 1.0d0) * im
c	  write(*, *) lambda
	  sum = sum + (dcmplx(1.0, 0.0) - lambda / dcmplx(2.0, 0.0)) /
     &      ((dcmplx(1.0, 0.0) - m / dcmplx(2.0, 0.0)) * lambda + m)
	  pro = pro * 
     &      ((dcmplx(1.0, 0.0) - m / dcmplx(2.0, 0.0)) * lambda + m)
	end do
	det(ifile) = dreal(pro) / 1.0d50
	sigma(ifile) = (sum / dcmplx(nsites * 1.0, 0.0)) * det(ifile)**nf
	
	
c	write(*, *) ifile, dreal(sigma(ifile)), det(ifile) * 1.0d50
c	write(3, '(a46, f4.2, i4, i3, $)') eigfile, m
c	write(3, '(2g23.15)') dreal(sigma), dimag(sigma)

	close(2)
	go to 11

c	  >>> this is the end
99    continue
c	  close(3)

c	  >>> number of input files = number of computed sigma's
	  nfile = ifile

	  sum = dcmplx(0.0d0, 0.0d0)
	  dsum = 0.0d0
	  do i = 1, nfile
		sum = sum + sigma(i)
		dsum = dsum + det(i)**nf
	  end do
c	  >>> write avarage sigma
	  sigma_av = sum / dble(nfile)
	  det_av = dsum / dble(nfile)
	  write(*, '(1x, f6.3, $)') m
	  write(*, *) dreal(sigma_av) / det_av

c	  >>> end of loop over m
      m = m - 0.001d0
	  rewind 1
	  if(m .gt. 0) goto 22

	  close(1)
	  write(*, *)
	  write(*, *) '... done!'
	  write(*, *)

	  end

