C     ******************************************************************
      subroutine v_invert1(ch)
      character*32 ch
      write(ch, '(a, e7.1, a)') 'BiCGgamma5 (', ACCURACY1, ') v2'
      end
C     ******************************************************************
      integer*4 function invert1(x, b)
C     ******************************************************************
C     **                       **                                     **
C     **  INVERT1 (BiCGgamma5) **    Program by I. Hip, 11 Oct 95     **
C     **  v2                   **    Last modified: 15 Feb 97         **
C     **                       **                                     **
C     ******************************************************************
C     Based on pseudo-code from: P. de Forcrand "Progress on lattice QCD
C     algorithms", hep-lat/9509082
C     ------------------------------------------------------------------
C     Program compute x which satisfies Q x = b
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input variables:   vector x, b
C     ------------------------------------------------------------------
C     Output variables:  vector x
C     ------------------------------------------------------------------
C     Remarks:- index arrays have to be initialized (call mk_index(ind))
C             - subroutines MVEC and MTVEC should be available
c             - maybe the vectors r, p, u, t should be put into a 
c               common for dummy storage
C     -----------------------------------------------------------------
C     On exit x contains Q^(-1) b   (i.e. the inverse of Q)
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (accuracy=ACCURACY1,maxits=1000,max_limit=100)
C     ------------------------------------------------------------------
      real*8 x(nall), b(nall)
      real*8 r(nall), p(nall)
      real*8 u(nall)
      real*8 delta, delta_old, ddelta, omega, sum, epsilon
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      epsilon=nall*accuracy**2
      limit = 0
      niter = 0
      
      call mvec(x, r)
      delta = 0.0
      do i = 1, nall, 4
	j = i + 1
	k = j + 1
	l = k + 1
	r(i) = b(i) - r(i)
	r(j) = b(j) - r(j)
	r(k) = b(k) - r(k)
	r(l) = b(l) - r(l)
	p(i) = r(i)
	p(j) = r(j)
	p(k) = r(k)
	p(l) = r(l)
	delta = delta + r(i)*r(i) + r(j)*r(j) - r(k)*r(k) - r(l)*r(l)
      end do
C     -----------------------------------------------------------------
C     ITERATION:
C     -----------------------------------------------------------------
C     Step (1)         omega = delta / < p | gamma5 Q | p>              
C     -----------------------------------------------------------------
10    call mvec(p, u) 
      
      sum = 0.0
      do i = 1, nall, 4
	j = i + 1
	k = j + 1
	l = k + 1
	sum = sum + p(i)*u(i) + p(j)*u(j) - p(k)*u(k) - p(l)*u(l) 
      end do

      if(sum .eq. 0.0) then
	write(*, *) '!!! STOP in BiCGg5_1: <p|g5 Q^+|p> = 0.0 !!!'
	go to 13
      end if
      omega = delta / sum

C     -----------------------------------------------------------------
C     Step (2) + (3) + (4)
C     -----------------------------------------------------------------
      delta_old = delta
      delta = 0.0
      do i = 1, nall, 4
	j = i + 1
	k = j + 1
	l = k + 1
	x(i) = x(i) + omega * p(i)
	x(j) = x(j) + omega * p(j)
	x(k) = x(k) + omega * p(k)
	x(l) = x(l) + omega * p(l)
	r(i) = r(i) - omega * u(i)
	r(j) = r(j) - omega * u(j)
	r(k) = r(k) - omega * u(k)
	r(l) = r(l) - omega * u(l)
	delta = delta + r(i)*r(i) + r(j)*r(j) - r(k)*r(k) - r(l)*r(l)
      end do

C      write(*, *) niter, '. delta = ', delta

      if(dabs(delta) .lt. epsilon) goto 99

C     >>> instability exit:
      if(dabs(delta) .ge. dabs(delta_old)) then
	limit = limit + 1
	if(limit .gt. max_limit) then
	  invert1 = -2
	  return
	end if
      else
	limit = 0
      end if 

C     -----------------------------------------------------------------
C     Step (5)
C     -----------------------------------------------------------------
      ddelta = delta / delta_old
      
      do i = 1, nall
	p(i) = r(i) + ddelta * p(i)
      end do  
      
      niter = niter + 1
      if(niter .lt. maxits) goto 10 

C     >>> maxits is reached, niter is set to -1 which means:
C     >>> "no convergence achieved"
13    invert1 = -1
      return

99    continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     keep the next line during the test phase:
C      write(*, *) 'BiCGg5_1: niter = ', niter, '   <r|g5|r> = ', delta
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     -----------------------------------------------------------------
C     On exit x contains Q^(-1) b   (i.e. the inverse of Q)
C     -----------------------------------------------------------------
      invert1 = niter
      return
      end


C     ******************************************************************
      integer*4 function invert1t(x, b)
C     ******************************************************************
C     **                       **                                     **
C     **  INVERT1 (BiCGgamma5) **    Program by I. Hip, 11 Oct 95     **
C     **  v2                   **    Last modified: 15 Feb 97         **
C     **                       **                                     **
C     ******************************************************************
C     Based on pseudo-code from: P. de Forcrand "Progress on lattice QCD
C     algorithms", hep-lat/9509082
C     ------------------------------------------------------------------
C     Program compute x which satisfies (Q^+) x = b
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input variables:   vector x, b
C     ------------------------------------------------------------------
C     Output variables:  vector x
C     ------------------------------------------------------------------
C     Remarks:- index arrays have to be initialized (call mk_index(ind))
C             - subroutines MVEC and MTVEC should be available
c             - maybe the vectors r, p, u, t should be put into a 
c               common for dummy storage
C     -----------------------------------------------------------------
C     On exit x contains (Q^+)^(-1) b   (i.e. the inverse of Q^+)
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (accuracy=ACCURACY1,maxits=1000,max_limit=100)
C     ------------------------------------------------------------------
      real*8 x(nall), b(nall)
      real*8 r(nall), p(nall)
      real*8 u(nall)
      real*8 delta, delta_old, ddelta, omega, sum, epsilon
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      epsilon=nall*accuracy**2
      limit = 0
      niter = 0
      
      call mtvec(x, r)
      delta = 0.0
      do i = 1, nall, 4
	j = i + 1
	k = j + 1
	l = k + 1
	r(i) = b(i) - r(i)
	r(j) = b(j) - r(j)
	r(k) = b(k) - r(k)
	r(l) = b(l) - r(l)
	p(i) = r(i)
	p(j) = r(j)
	p(k) = r(k)
	p(l) = r(l)
	delta = delta + r(i)*r(i) + r(j)*r(j) - r(k)*r(k) - r(l)*r(l)
      end do
C     -----------------------------------------------------------------
C     ITERATION:
C     -----------------------------------------------------------------
C     Step (1)         omega = delta / < p | gamma5 Q^+ | p>              
C     -----------------------------------------------------------------
10    call mtvec(p, u) 
      
      sum = 0.0
      do i = 1, nall, 4
	j = i + 1
	k = j + 1
	l = k + 1
	sum = sum + p(i)*u(i) + p(j)*u(j) - p(k)*u(k) - p(l)*u(l) 
      end do

      if(sum .eq. 0.0) then
	write(*, *) '!!! STOP in BiCGg5_1t: <p|g5 Q^+|p> = 0.0 !!!'
	go to 13
      end if
      omega = delta / sum

C     -----------------------------------------------------------------
C     Step (2) + (3) + (4)
C     -----------------------------------------------------------------
      delta_old = delta
      delta = 0.0
      do i = 1, nall, 4
	j = i + 1
	k = j + 1
	l = k + 1
	x(i) = x(i) + omega * p(i)
	x(j) = x(j) + omega * p(j)
	x(k) = x(k) + omega * p(k)
	x(l) = x(l) + omega * p(l)
	r(i) = r(i) - omega * u(i)
	r(j) = r(j) - omega * u(j)
	r(k) = r(k) - omega * u(k)
	r(l) = r(l) - omega * u(l)
	delta = delta + r(i)*r(i) + r(j)*r(j) - r(k)*r(k) - r(l)*r(l)
      end do

C      write(*, *) niter, '. delta = ', delta

      if(dabs(delta) .lt. epsilon) goto 99

C     >>> instability exit:
      if(dabs(delta) .ge. dabs(delta_old)) then
	limit = limit + 1
	if(limit .gt. max_limit) then
	  invert1t = -2
	  return
	end if
      else
	limit = 0
      end if 

C     -----------------------------------------------------------------
C     Step (5)
C     -----------------------------------------------------------------
      ddelta = delta / delta_old
      
      do i = 1, nall
	p(i) = r(i) + ddelta * p(i)
      end do  
      
      niter = niter + 1
      if(niter .lt. maxits) goto 10 

C     >>> maxits is reached, niter is set to -1 which means:
C     >>> "no convergence achieved"
13    invert1t = -1
      return

99    continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     keep the next line during the test phase:
C      write(*, *) 'BiCGg5_1t: niter = ', niter, '   <r|g5|r> = ', delta
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     -----------------------------------------------------------------
C     On exit x contains (Q^+)^(-1) b   (i.e. the inverse of Q^+)
C     -----------------------------------------------------------------
      invert1t = niter
      return
      end

