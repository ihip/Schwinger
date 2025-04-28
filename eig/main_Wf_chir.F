      program compute_chiralities

c	  ihip - 22 Aug 24 - last modified: 22 Aug 24

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
	  character*64 input_file

c     >>> 'dense' working array
	  complex*16 a(nferm, nferm)

c     >>> EVs
	  complex*16 eig(nferm)
	  complex*16 reig(maxevect), vr(nferm, maxevect)
	  real*8 rtol, edetr, edeti

c     >>> unit vector needed for dense fill-in
      complex*16 s(nferm)

c     >>> chiralities (for real modes only)
      real*8 chir(maxevect)

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
	  write(*,'(1x,''Compute chiralities (Hip, 22 Aug 24)'')')
	  write(*, *)   '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	  write(*, *)
	  write(*, '(1x, ''Lattice:   '', i3, '' x '', i3)') NSPACE, NTIME
	  write(*, *)
	  write(*, '(1x, ''rtol set to: '', d10.3)') rtol
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

c	  do i = 1, nsite
c	    do j = 1, ndim
c	      write(*, *) i, j, u(i, j)
c	    end do
c	  end do

c	  >>> compute topological charge
	  nu = itopch()
	  print *, 'nu   = ', nu
	  write(*, *)
	  write(*, *)

c     >>> read kappa
	  write(*, '(1x, ''Enter kappa: '', $)')
      read(*, *) akap
      write(*, *) 'kappa = ', akap

c	  >>> fill-in the dense matrix
      do i = 1, nferm
        s(i) = dcmplx(0.0d0, 0.0d0)
      end do
      do k = 1, nferm
        s(k) = dcmplx(1.0d0, 0.0d0)
        call mvec(s, a(1, k))
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

c     >>> compute chiralities of real modes
      do i = 1, nreal
        chir(i) = 0.0d0
        do k = 1, nferm - 1, 2
          chir(i) = chir(i) + dconjg(vr(k, i)) * vr(k, i) -
     &      dconjg(vr(k + 1, i)) * vr(k + 1, i)
        end do
        write(*, *) dreal(reig(i)), dimag(reig(i)), chir(i)
      end do
	  write(*, *)
      
c	  >>> write eigenvalues
	  call save_eig(nferm, eig, input_file, 'Wf')

	  end

