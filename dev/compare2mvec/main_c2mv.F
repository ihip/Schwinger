      program compare2mvec

c	  ihip - 10 Aug 20 - last modified: 10 Aug 20

c	  >>> necessary initializations for 2D U(1)
	  parameter (ndrc=2,ngrp=1,ndim=2)
	  parameter (nsite=NTIME*NSPACE)
	  parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	  parameter (nall=2*ndrc*ngrp*nsite)

C     >>> index matrix - be sure to call mk_index to initialize !!!
      integer ind(4, nsite)

c     >>> gauge fields
	  complex*16 u(nsite, ndim)

c     >>> some necessary definitions
	  real*8 beta, akap, eps

c     >>> gauge file name
	  character*64 input_file

c     >>> vectors
	  complex*16 a(ndrc, nsite), m1xa(ndrc, nsite), m2xa(ndrc, nsite)

	  complex*16 diff

C     ------------------------------------------------------------------
	  common /gauge_fields/ u
	  common /index_arrays/ ind
	  common /kappa/ akap
C     ------------------------------------------------------------------

C     >>> initialize index arrays
	  call mk_index()
C	  write(*, *) 'mk_index() ... passed!'
	  call make_ihf()
C	  write(*, *) 'make_ihf() ... passed!'

	  write(*, *)
	  write(*,'(1x,''Compare2mvec (ihip, 10 Aug 20)'')')
	  write(*, *)   '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	  write(*, *)
	  write(*, '(1x, ''Lattice:   '', i3, '' x '', i3)') NSPACE, NTIME
	  write(*, *)

C     >>> read .conf file name
	  write(*, '(1x, ''Enter .conf file name: '', $)')
	  read(*, '(a)') input_file
	  write(*, *)

C	  >>> read .conf file
	  open(1, file = input_file, form = 'unformatted', status = 'old')
	    read(1) ntime, nspace
	    read(1) beta, akap, eps
	    read(1) it
	    read(1) u
	  close(1)
	
	  write(*, '(1x, ''Lattice: '', i3, '' x '', i3)') nspace, ntime
	  if(nspace .ne. NSPACE) stop 'nspace different from NSPACE!'
	  if(ntime .ne. NTIME) stop 'ntime different from NTIME!'

	  print *, 'beta = ', beta	  
	  print *, 'akap = ', akap
	  print *, 'eps  = ', eps
	  print *, 'it   = ', it
	  print *

	  write(*, '(1x, ''Enter n: '', $)')
	  read(*, *) n
	  write(*, *)

c	  do i = 1, nsite
c	  	u(i, 1) = dble(i)
c	  	u(i, 2) = dble(i)
c	  end do

	  do i = 1, nsite
	  	a(1, i) = dcmplx(0.0d0, 0.0d0)
		a(2, i) = dcmplx(0.0d0, 0.0d0)
		if (i .eq. n) then
		  a(1, i) = dcmplx(1.0d0, 0.0d0)
		  a(2, i) = dcmplx(0.0d0, 0.0d0)
		end if
	  end do

C	  >>> multiplication with the 1st matrix
	  akap = 1.0
	  call mvec(a, m1xa)

C	  >>> multiplication with the 2nd matrix
C	  call init_hf('GWRoptimal_ws   ')
C	  call mvec_dhf(a, m2xa)
	  akap = 1.0
	  call mvec_Wilson(a, m2xa)

	  do i = 1, nsite
	    do j = 1, ndrc
	  	  write(*, *) i, j, m1xa(j, i), m2xa(j, i)
		end do
	  end do
	  write(*, *)

	  do i = 1, nsite
	    do j = 1, ndrc
	      diff = m2xa(j, i) - m1xa(j, i)
	  	  if(cdabs(diff) .gt. 0.0d0) write(*, *) i, j, diff
		end do
	  end do

	  end
