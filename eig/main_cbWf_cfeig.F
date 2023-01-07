      program compute_eigenvalues
c		of central-branch Wilson fermions (cbWf)
c	  ihip - 24 Nov 07 - last modified: 07 Jan 23

c	  >>> necessary initializations for 2D U(1)
	  parameter (ndrc=2,ngrp=1,ndim=2)
	  parameter (nsite=NTIME*NSPACE)
	  parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	  parameter (nall=2*ndrc*ngrp*nsite)

c	  >>> maximal number of real eigenvalues
	  parameter (maxevect=64)

c     >>> gauge fields
	  complex*16 u(nsite, ndim)

c     >>> some necessary definitions
	  real*8 beta, akap, eps
		
c     >>> gauge file name
	  character*64 input_file, list_name

c     >>> 'dense' working array
	  complex*16 a(nferm, nferm)

c     >>> auxiliary vector
	  complex*16 s(nferm)

c     >>> EVs
	  complex*16 eig(nferm)
	  complex*16 reig(maxevect), vr(nferm, maxevect)
	  real*8 rtol, edetr, edeti

C     ------------------------------------------------------------------
	  common /gauge_fields/ u
	  common /kappa/ akap
      common /dummy_main/ a, g 
C     ------------------------------------------------------------------

C     >>> initialize index arrays
	  call mk_index()

C	  >>> set tolerance for real eigenvalues
	  rtol = 1.0d-12

	  write(*, *)
	  write(*,'(1x,''cbWf_cfeig - eigenvalues (Hip, 07 Jan 23)'')')
	  write(*, *)   '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	  write(*, *)
	  write(*, '(1x, ''Lattice:   '', i3, '' x '', i3)') NSPACE, NTIME
	  write(*, *)
	  write(*, '(1x, ''rtol set to: '', d10.3)') rtol
	  write(*, *)

C     >>> read .conf list file name
	  write(*, '(1x, ''Enter list of .conf files: '', $)')
	  read(*, '(a)') list_name
	  write(*, *)

	  open(2, file = list_name, form = 'formatted', status = 'old')

c	  >>> loop through the files in the list
10    read(2, '(a)', end = 99) input_file

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

c	  >>> compute topological charge
	  nu = itopch()
	  print *, 'nu   = ', nu
	  write(*, *)
	  write(*, *)

c	  >>> lattice Dirac operator
	  write(*, *) '### cbWf operator ###'

c	  >>> fill-in the central-branch Wilson fermions
	  do i = 1, nferm
	    s(i) = dcmplx(0.0d0, 0.0d0)
	  end do
	  do k = 1, nferm
	    s(k) = dcmplx(1.0d0, 0.0d0)
	    call mvec_cbWf(s, a(1, k))
	    s(k) = dcmplx(0.0d0, 0.0d0)
	  end do

c     >>> compute the eigenvalues of a (WARNING: a is overwritten!)
	  write(*, *) 'computing eigenvalues...'
	  write(*, *)

	  nreal = 0
	  nuf = 0
	  call lapackeig(a, eig, rtol, nreal, reig, vr, maxevect,
     &    edetr, edeti, nuf)

	  print *, 'nreal = ', nreal	  
	  print *, 'nuf   = ', nuf
	  print *, 'edetr = ', edetr
	  print *, 'edeti = ', edeti
	  write(*, *)

c	  >>> write eigenvalues to file
	  call save_eig(nferm, eig, input_file, 'CB')
	  goto 10

99	  continue
	  end

