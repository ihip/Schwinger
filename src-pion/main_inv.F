c >>> riješiti učitavanje nu i nuf
c >>> pokrpati sve Makefileove

	program compute_inverse_of_PHF
c       hip - 17 Jun 08 - last modified: 26 Jun 08

c       >>> necessary initializations for 2D U(1)
	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)

c       >>> needed for condensates
	real*8 pbp, pbg5p, uudd, ug5udg5d
	real*8 op_pbp_e, op_pbg5p_e, op_ubudbd_e, op_ubg5udbg5d_e

c	>>> DSP: needed for dispersion relations check
	real*8 dsp(NTIME - 1, 0:NSPACE / 2, 4)
	real*8 con(0:NSPACE / 2, 2)

c       >>> gauge fields
	complex*16 u(nsite, ndim)

c	>>> mass
	real*8 mass, massscale

c       >>> some necessary definitions
	real*8 beta, akap, eps

c       >>> determinant - needed for compatibility with massan
	real*8 detr, edetr
	parameter (detr=0.0,edetr=0.0)

c       >>> gauge file name
	character*64 listname, fname, maspname, logname

c       >>> 'dense' working array
	complex*16 a(nferm, nferm)

c	>>> shift of the hy-fermions before Neubergerization
	real*8 mu_shift

	real*8 rtol

c	>>> needed for tests
	complex*16 sum

C       ------------------------------------------------------------------
	common /gauge_fields/ u
	common /kappa/ akap
        common /dummy_main/ a, g 
C       ------------------------------------------------------------------

C       >>> initialize index arrays
	call mk_index()
	call make_ihf()

	write(*, *)
	write(*,'(1x,''Compute inverse of PHF (hip, 26 Jun 08)'')')
	write(*, *)   '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)
	write(*, '(1x, ''Lattice:   '', i2, '' x '', i2)') NTIME, NSPACE
	write(*, *)

C       >>> read gauge field name and parameters
	write(*, '(1x, ''List of gauge configurations: '', $)')
	read(*, '(a)') listname

	open(1, file = listname, form = 'formatted', status = 'old')
	nfiles = 0 

c	>>> read mass
	write(*, '(1x, ''mass: '', $)')
	read(*, *) mass

c	>>> open log file
	logname = listname(1:lnblnk(listname))//'.lgm'
        write(*, *) 'logname: ', logname
	open(2, file = logname, form = 'formatted', status = 'unknown')

c	>>> open .masp file
	write(*, '(1x, ''.masp name: '', $)')
	read(*, *) maspname
	open(3, file = maspname, form = 'unformatted', status = 'unknown')
c	>>> write header
	nmeas = 0
	beta = 5.0
	write(3) NSPACE, NTIME, nmeas, beta, mass	

c	>>> set mu_shift = 1.0d0
	mu_shift = 1.0d0

c	write(*, '(1x, ''rtol: '', $)')
c	read(*, *) rtol
	rtol = 0.000001

	write(*, '(1x, ''Top. charge (0-selected, 1-all): '', $)')
	read(*, *) iall
	if(iall .eq. 0) then
	  write(*, '(1x, ''  Which |nu| to collect: '', $)')
	  read(*, *) inu
	end if
	write(*, *)

c       >>> main loop
11	read(1, *, end = 99) fname, nu, nuf
	nfiles = nfiles + 1
	write(*, *) 'Processing: ', fname

	call readconf(fname, u)

c	do i = 1, nsite
c	  do j = 1, ndim
c	    write(*, *) i, j, u(i, j)
c	  end do
c	end do

	if((iall .eq. 1) .or. (nuf .eq. inu)) then

	i0 = i0 + 1
	write(*, *) 'meas. no. ### ', i0
        write(*, *) 'nuf = ', nuf

c       >>> fill in of appropriate Dirac operator

c	  >>> fill-in the hypercube fermions
	  call init_hf('GWRoptimal_ws   ')
	  call fill_dhf(a)
c	  >>> appropriate shift is necessary:
	  do i = 1, nferm
	    a(i, i) = a(i, i) - mu_shift
	  end do
c	  >>> perform "neubergerization"
	  call neuberger(a)

c	>>> scale according to mass parameter
c	>>> D(m) = (1 - m/2)D(0) + m
        massscale = 1.0d0 - mass / 2.0d0
        do i = 1, nferm
	  do j = 1, nferm
	     a(j, i) = massscale  * a(j, i)
	  end do
	end do
c	>>> mass-shift
	do i = 1, nferm
	  a(i, i) = a(i, i) + mass
	end do

c       >>> compute inverse (a becomes inverse!)
	write(*, *) 'inversion'
	call invert(a)

c	>>> compute pbp
c	sum = dcmplx(0.0d0, 0.0d0)
c	do i = 1, nferm
c	  sum = sum + a(i, i) 
c	end do
c	write(*, *) 'pbp = ', sum

c ------------- taken from np.F (version: 02 Feb 1999) ----------------
c       >>> compute condensates
	pbp = op_pbp_e(a)
	pbg5 = op_pbg5p_e(a)
	uudd = op_ubudbd_e(a)
	ug5udg5d = op_ubg5udbg5d_e(a)

c       >>> compute dispersion relations
	call massdsp(a, dsp, con)

c       >>> write to a .masp file
	write(3) nu, nuf, rtol
	write(3) detr, edetr, pbp, pbg5p, uudd, ug5udg5d
	write(3) (((dsp(j, ip, ic), j = 1, NTIME - 1),
     &    ip = 0, NSPACE / 2), ic = 1, 4)  
	write(3) ((con(ip, ic), ip = 0, NSPACE / 2), ic = 1, 2)  
c ---------------------------------------------------------------------

c	>>> write entry to log-file
	write(2, *) fname(1:lnblnk(fname)), mass, nu, nuf, pbp, pbg5

	write(*, *) 

	end if

	goto 11

99	write(*, *) 'nfiles = ', nfiles

c	>>> close .masp-file
	close(3)

c       >>> close log-file
        close(2)
c	>>> close list-file
	close(1)

	end




