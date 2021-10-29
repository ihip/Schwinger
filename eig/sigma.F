c	>>> sigma >>> Ivan Hip >>> 2008-01-19 >>> 2021-10-29	

#define MAX_LIST_LEN 10000
#define MAX_EIGEN 1000

	character*64 listname, eigfile, logname
	real*8 re, im, m
	complex*16 lambda(MAX_EIGEN), sum, pro
	complex*16 sigma(MAX_LIST_LEN)
	real*8 det(MAX_LIST_LEN)
	real*8 dsum, sigma_av
	real*8 aa, var

	write(*, *)
	write(*, *) 'Sigma / 2021-10-29'
	write(*, *) '=================='
	write(*, *)

c       >>> read listname
	write(*, '(1x, ''List name of .eig files: '', $)')
	read(*, '(a)') listname

	open(1, file = listname, form = 'formatted', status = 'old')

c	>>> open log file
c	logname = listname(1:lnblnk(listname))//'.sig'
c	open(3, file = logname, form = 'formatted', status = 'unknown')

c	>>> read number of flavors
	write(*, '(1x, ''N_f = '', $)')
	read(*, *) nf
	write(*, *)

c	>>> read number of jackknife blocks
	write(*, '(1x, ''jkblocks = '', $)')
	read(*, *) jkblocks
	write(*, *)

c	>>> read mass parameter
c	write(*, '(1x, ''mass: '', $)')
c	read(*, *) m
c	write(*, *)

c	  >>> loop over m
      m = 0.15d0

22    ifile = 0
c     >>> loop over all configurations

11    read(1, *, end = 99) eigfile
	  ifile = ifile + 1

c	  >>> read eigenvalues from eigfile
	  open(2, file = eigfile, form = 'formatted', status = 'old')
	  neigen = 0
7	  read(2, *, end = 77) re, im
	  	neigen = neigen + 1
c	  	if(dabs(re) .lt. 1.0d-12) re = 0.0d0
c	  	if(dabs(im) .lt. 1.0d-12) im = 0.0d0
	  	lambda(neigen) = re + dcmplx(0, 1) * im
	  	goto 7
77	  continue
	  close(2)
	  nsites = neigen / 2

c	  >>> compute sigma
	  sum = dcmplx(0.0, 0.0)
	  pro = dcmplx(1.0d0, 0.0)
      do i = 1, neigen
	    sum = sum + (dcmplx(1, 0) - lambda(i) / dcmplx(2, 0)) /
     & ((dcmplx(1, 0) - m / dcmplx(2, 0)) * lambda(i) + m)
	    pro = pro * 
     & ((dcmplx(1, 0) - m / dcmplx(2, 0)) * lambda(i) + m)
	  end do
	  det(ifile) = dreal(pro) / 1.0d50
	  sigma(ifile) = (sum / dcmplx(nsites, 0)) * det(ifile)**nf

	  go to 11

c	  >>> this is the end
99    continue
c	  close(3)

c	  >>> number of input files = number of computed sigmas
	  nconf = ifile

	  sum = dcmplx(0.0d0, 0.0d0)
	  dsum = 0.0d0
	  do i = 1, nconf
		sum = sum + sigma(i)
		dsum = dsum + det(i)**nf
	  end do
c	  >>> write average sigma
	  sigma_av = dreal(sum) / dsum

c	  >>> compute njack sigmas
	  call jack_sigma(nconf, jkblocks, sigma, det, nf, aa, var)
c	  print *, 'aa = ', aa
c	  print *, 'sqrt(var) = ', dsqrt(var)

	  write(*, '(1x, f6.3, $)') m
	  write(*, *) sigma_av, dsqrt(var)

c	  >>> end of loop over m
      m = m - 0.001d0
	  rewind 1
	  if(m .gt. 0) goto 22

	  close(1)
	  write(*, *)
	  write(*, *) '... done!'
	  write(*, *)

	  end


C     ******************************************************************
      subroutine jack_sigma(nconf, jkblocks, sigma, det, nf, aa, var)
C     ******************************************************************
C     **                   **                                         **
C     ** JACK_SIGMA        **   I. Hip, 2021-10-29                    **
C     **                   **   Last modified: 2021-10-29             **
C     ******************************************************************
C     IN nmeas               - number of configurations
C     IN jkblocks            - number of jackknife subsamples to be used
C     IN (complex*16) sigma  - sigma for each configuration
C     IN (real*8) det        - determinant for each configuration
C     IN nf                  - number of flavors
C     OUT aa (real*8)        - average of jkblocks
C     OUT var (real*8)       - variance (result of jackknife)
C     ******************************************************************
	parameter(max_jkblocks=100)
	complex*16 sigma(*)
	real*8 det(*)
	complex*16 sum
	real*8 dsum
	real*8 a(max_jkblocks), sum_a, sum_a2
	real*8 aa, aa2, var

	if(jkblocks .gt. max_jkblocks) then
	  jkblocks = max_jkblocks
	  write(*, *) 'Jackknife WARNING: number of blocks to big, ',
     &      'set to ', max_jkblocks, '.'
	end if

	if(mod(nconf, jkblocks) .ne. 0) then
	  write(*, *) 'mod(nconf, jkblocks) != 0'
	  write(*, *) 'Jackknife ERROR: not all conf. included'
	  stop
	end if

C   >>> prepare jkblocks different ensembles
	npart = nconf / jkblocks

	do i = 1, jkblocks
	  sum = dcmplx(0.0d0, 0.0d0)
	  dsum = 0.0d0
	  do j = 1, nconf
	    if((j .lt. (i - 1) * npart + 1) .or. (j .gt. i * npart)) then
		  sum = sum + sigma(j)
		  dsum = dsum + det(j)**nf
		end if
	  end do
	  a(i) = dreal(sum) / dsum
c	  print *, i, a(i)
	end do

	sum_a = 0.0d0
	sum_a2 = 0.0d0
	do i = 1, jkblocks
	  sum_a = sum_a + a(i)
	  sum_a2 = sum_a2 + a(i)**2
	end do
	aa = sum_a / dble(jkblocks)
	aa2 = sum_a2 / dble(jkblocks)
	var = dble(jkblocks - 1) * (aa2 - aa**2)
	
	return
	end

