	program compute_inverse_of_cbWf
c   >>> Hip / 17 Jun 08 / last modified: 2023-01-08

c   >>> necessary initializations for 2D U(1)
	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)

c	>>> needed for old header
	parameter (mcomp=1)
	character*64 c_mcomputable(mcomp)
	character*64 c_source, c_direction, c_mfuture

c   >>> needed for tail
	real*4 tcpu

c   >>> needed for condensates
	real*8 pbp, pbg5p, uudd, ug5udg5d
	real*8 op_pbp_e, op_pbg5p_e, op_ubudbd_e, op_ubg5udbg5d_e

c	>>> DSP: needed for dispersion relations check
	real*8 dsp(NTIME - 1, 0:NSPACE / 2, 4)
	real*8 con(0:NSPACE / 2, 2)

c   >>> gauge fields
	complex*16 u(nsite, ndim)

c   >>> some necessary definitions
	real*8 beta, akap, eps, rtol 

c   >>> determinant
	real*8 detr, edetr

c   >>> file names handling
	character*64 listname, outname, cfname, maspname
	character*4 c_nspace, c_ntime
	character*10 c_size
	character*6 c_beta, c_fmass

c   >>> 'dense' working array
	complex*16 a(nferm, nferm)

c	>>> point 'source' used for fill in
	complex*16 s(nferm)

c	>>> chiral condensate
	real*8 sigma

C   ------------------------------------------------------------------
	common /gauge_fields/ u
	common /kappa/ akap
	common /dummy_main/ a, g 
C   ------------------------------------------------------------------

C   >>> initialize index arrays
	call mk_index()

C	>>> set tolerance for identification of real eigenvalues
	rtol = 1.0d-12

	write(*, *)
	write(*,'(1x,''Compute inverse of cbWf / Hip, 2023-01-08'')')
	write(*, *)   '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)
	write(*, '(1x, ''Lattice:   '', i2, '' x '', i2)') NSPACE, NTIME
	write(*, *)

C   >>> read file with list of gauge configurations
	write(*, '(1x, ''List of gauge configurations: '', $)')
	read(*, '(a)') listname

	open(1, file = listname, form = 'formatted', status = 'old')
c	>>> count configurations
	nmeas = 0
22	read(1, '(a)', end = 33) cfname
	nmeas = nmeas + 1
	goto 22

33	close(1)
	write(*, *) 'nmeas = ', nmeas
	write(*, *)	

	write(*, '(1x, ''Top. charge (0-selected, 1-all): '', $)')
	read(*, *) iall
	if(iall .eq. 0) then
	  write(*, '(1x, ''  Which |nu| to collect: '', $)')
	  read(*, *) inu
	end if
	write(*, *)

c	>>> loop over all configurations
	nfiles = 0
	open(1, file = listname, form = 'formatted', status = 'old')	
11	read(1, '(a)', end = 99) cfname
	nfiles = nfiles + 1
	write(*, *) 'Processing: ', cfname

	call readconf(cfname, beta, u)

c	>>> open .masp file and write header
	if(nfiles .eq. 1) then
C   	>>> prepare the file name L[nspace]x[ntime]_b[beta]m[fmass].data
		write(c_nspace, '(i4)') NSPACE + 1000
		write(c_ntime, '(i4)') NTIME + 1000
		c_size = 'L'//c_nspace(2:4)//'x'//c_ntime(2:4)//'_b'
		ival = int((beta + 100.0000001d0) * 1000.0d0)
		write(c_beta, '(i6)') ival
		ival = int((fmass + 10.00000001d0) * 10000.0d0)
		write(c_fmass, '(i6)') ival
		maspname = c_size//c_beta(2:6)//'m'//c_fmass(2:6)//'.masp'
		open(3, file = maspname, form = 'unformatted',
     &		status = 'unknown')
c		write(3) NSPACE, NTIME, nmeas, beta, fmass
		call save_header(3, 0, 'dsp-cbW', NSPACE, NTIME, beta, fmass,
     &		0.0, 0, 0, nmeas, 0, 0, 0)

C       .mass specific part - saving desc. of mass computables
		c_source = 'point'
		c_direction = 'time'
		write(3) mcomp, c_source, c_direction, c_mfuture
		write(3) (c_mcomputable(i), i = 1, mcomp)

		write(*, *) 'beta = ', beta
		write(*, *) '.masp file created: ', maspname
		write(*, *)

c	    >>> start measuring CPU time
		call cbl_cputime(0, tcpu)
	end if

	call wdet(cfname, fmass, rtol, nuf, edetr, sigma)

c	>>> topological charge (geometrical definition)
	nu = itopch()

	if((iall .eq. 1) .or. (nuf .eq. inu)) then

	write(*, *) 'meas. no. ### ', nfiles
	write(*, *) 'nu | nuf = ', nu, ' | ', nuf
	if(iabs(nu) .ne. nuf) write(*, *) 'ALERT: nu <> nuf !!!'

c   >>> fill in of appropriate Dirac operator
	do i = 1, nferm
	  s(i) = dcmplx(0.0d0, 0.0d0)
	end do
	do k = 1, nferm
	  s(k) = dcmplx(1.0d0, 0.0d0)
c     >>> central-branch Wilson fermions
	  call mvec_cbWf(s, a(1, k))
	  s(k) = dcmplx(0.0d0, 0.0d0)
	end do

c   >>> compute inverse (a becomes inverse!)
	write(*, *) 'inversion'
	call invert(a)

c ------------- taken from np.F (version: 02 Feb 1999) ----------------
c   >>> compute condensates
	pbp = op_pbp_e(a)
	write(*, *) 'pbp = ', pbp
	pbg5 = op_pbg5p_e(a)
	uudd = op_ubudbd_e(a)
	ug5udg5d = op_ubg5udbg5d_e(a)

c   >>> compute dispersion relations
	call massdsp(a, dsp, con)

c   >>> write to a .masp file
	write(3) nu, nuf, rtol
	write(3) sigma, edetr, pbp, pbg5p, uudd, ug5udg5d
	write(3) (((dsp(j, ip, ic), j = 1, NTIME - 1),
     &    ip = 0, NSPACE / 2), ic = 1, 4)  
	write(3) ((con(ip, ic), ip = 0, NSPACE / 2), ic = 1, 2)  
c ---------------------------------------------------------------------

	write(*, *) '---------------------------------------------------'

	end if

	goto 11

99	write(*, *) 'nfiles = ', nfiles
	write(*, *)

c	>>> close .masp-file
	call cbl_cputime(1, tcpu)
	call save_tail(3, 0, 0, 0, 0, tcpu)
	close(3)

c	>>> close list-file
	close(1)

	end
