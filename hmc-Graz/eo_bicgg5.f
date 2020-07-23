C     EVEN-ODD preconditioning assumed

C     ******************************************************************
      subroutine v_invert(ch)
      character*32 ch
      write(ch, '(a, e7.1, a)') 'eo_BiCGg5 (', ACCURACY, ') v3'
      end
C     ******************************************************************
      function invert1eo(xx, bb)
C     - inversion of even-odd preconditioned matrix
C     ******************************************************************
C     **                     **                                       **
C     ** INVERT (BiCGgamma5) **    Program by I. Hip, 11 Oct 95       **
C     **                     **    Last modified: 13 Mar 97           **
C     ******************************************************************
C     This implementation of BiCGgamma5 algorithm is based on paper by
C     P. de Forcrand: "Progress on lat. QCD algorithms", hep-lat/9509082
C     ------------------------------------------------------------------
C     Program compute x which satisfies Q x = b
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input var: vector b (x vector is internally initialized with b)
C     ------------------------------------------------------------------
C     Output variables: vector x
C     ------------------------------------------------------------------
C     Remarks:- index arrays have to be initialized (call mk_index(ind))
C             - subroutines MVEC and MTVEC should be available
c             - maybe the vectors u, r, p, s should be put into a
c               common for dummy storage
C     -----------------------------------------------------------------
C     On exit x contains Q^(-1) b   (i.e. the inverse of Q)
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (nh=nall/2)
C     >>> only for bicgg5: correction on acc
      parameter (acc=ACCURACY*0.05,maxits=1000,max_limit=100)
C     ------------------------------------------------------------------
      real*8 xx(nall), bb(nall)
      real*8 x(nh), b(nh)
      real*8 r(nh), p(nh)
      real*8 u(nh), t(nh)
      real*8 delta, delta_old, ddelta, omega, sum, epsilon, bnorm2
C     ------------------------------------------------------------------
      save x
      data x /nh * 0.0/
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      limit = 0
      niter = 0

      call mdeo(bb(nh + 1), b)
      bnorm2 = 0.0d0
      do i = 1, nh
	b(i) = b(i) + bb(i)
        bnorm2 = bnorm2 + b(i)**2
      end do
      epsilon = bnorm2 * acc**2

C     initial guess x(0) is set to vector b
C     (do not remove blindly - this is assumed few lines later...)
c      do i = 1, nh
c	x(i) = b(i)
c      end do

      call mdoe(x, t)
      call mdeo(t, r)
      delta = 0.0
      do i = 1, nh, 4
	j = i + 1
	k = j + 1
	l = k + 1
C     the following 4 lines are not necessary if x = b
       r(i) = b(i) - x(i) + r(i)
       r(j) = b(j) - x(j) + r(j)
       r(k) = b(k) - x(k) + r(k)
       r(l) = b(l) - x(l) + r(l)
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
10    call mdoe(p, t)
      call mdeo(t, u)

      sum = 0.0
      do i = 1, nh, 4
	j = i + 1
	k = j + 1
	l = k + 1
C these 4 lines are part of multiplication
	u(i) = p(i) - u(i)
	u(j) = p(j) - u(j)
	u(k) = p(k) - u(k)
	u(l) = p(l) - u(l)
	sum = sum + p(i)*u(i) + p(j)*u(j) - p(k)*u(k) - p(l)*u(l)
      end do

      omega = delta / sum

C     -----------------------------------------------------------------
C     Step (2) + (3) + (4)
C     -----------------------------------------------------------------
      delta_old = delta
      delta = 0.0
      do i = 1, nh, 4
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
          do i = 1, nh
            x(i) = 1.0d0
          end do
          invert1eo = -2
          return
        end if
      else
        limit = 0
      end if 

C     -----------------------------------------------------------------
C     Step (5)
C     -----------------------------------------------------------------
      ddelta = delta / delta_old
      
      do i = 1, nh
	p(i) = r(i) + ddelta * p(i)
      end do

      niter = niter + 1
      if(niter .lt. maxits) goto 10 

C     >>> maxits is reached, niter is set to -1 which means:
C     >>> "no convergence achieved"
C     >>> return without changing starting guess xx      
13    do i = 1, nh
        x(i) = 1.0d0
      end do
      invert1eo = -1
      return

99    continue

      call mdoe(x, t)
      do i = 1, nh
	xx(i) = x(i)
	xx(nh + i) = bb(nh + i) + t(i)
      end do

C      write(*, *) 'BiCGg5 : niter = ', niter, '  <r|g5|r> = ', delta 

      invert1eo = niter
      return
      end

C     ******************************************************************
      function invert1teo(xx, bb)
C     - inversion of even-odd preconditioned matrix
C     ******************************************************************
C     **                     **                                       **
C     ** INVERT (BiCGgamma5) **    Program by I. Hip, 11 Oct 95       **
C     **                     **    Last modified: 13 Mar 97           **
C     **                     **                                       **
C     ******************************************************************
C     This implementation of BiCGgamma5 algorithm is based on paper by
C     P. de Forcrand: "Progress on lat. QCD algorithms", hep-lat/9509082
C     ------------------------------------------------------------------
C     Program compute x which satisfies Q x = b
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input var: vector b (x vector is internally initialized with b)
C     ------------------------------------------------------------------
C     Output variables: vector x
C     ------------------------------------------------------------------
C     Remarks:- index arrays have to be initialized (call mk_index(ind))
C             - subroutines MVEC and MTVEC should be available
c             - maybe the vectors u, r, p, s should be put into a
c               common for dummy storage
C     -----------------------------------------------------------------
C     On exit x contains Q^(-1) b   (i.e. the inverse of Q)
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (nh=nall/2)
C     >>> only for bicgg5: correction on acc
      parameter (acc=ACCURACY*0.05,maxits=1000,max_limit=100)
C     ------------------------------------------------------------------
      real*8 xx(nall), bb(nall)
      real*8 x(nh), b(nh)
      real*8 r(nh), p(nh)
      real*8 u(nh), t(nh)
      real*8 delta, delta_old, ddelta, omega, sum, epsilon, bnorm2
C     ------------------------------------------------------------------
      save x
      data x /nh * 0.0/
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      limit = 0
      niter = 0

      call mtdeo(bb(nh + 1), b)
      bnorm2 = 0.0d0
      do i = 1, nh
	b(i) = b(i) + bb(i)
        bnorm2 = bnorm2 + b(i)**2
      end do
      epsilon = bnorm2 * acc**2

C     initial guess x(0) is set to vector b
C     (do not remove blindly - this is assumed few lines later...)
c      do i = 1, nh
c	x(i) = b(i)
c      end do

      call mtdoe(x, t)
      call mtdeo(t, r)
      delta = 0.0
      do i = 1, nh, 4
	j = i + 1
	k = j + 1
	l = k + 1
C     the following 4 lines are not necessary if x = b
       r(i) = b(i) - x(i) + r(i)
       r(j) = b(j) - x(j) + r(j)
       r(k) = b(k) - x(k) + r(k)
       r(l) = b(l) - x(l) + r(l)
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
10    call mtdoe(p, t)
      call mtdeo(t, u)

      sum = 0.0
      do i = 1, nh, 4
	j = i + 1
	k = j + 1
	l = k + 1
C these 4 lines are part of multiplication
	u(i) = p(i) - u(i)
	u(j) = p(j) - u(j)
	u(k) = p(k) - u(k)
	u(l) = p(l) - u(l)
	sum = sum + p(i)*u(i) + p(j)*u(j) - p(k)*u(k) - p(l)*u(l) 
      end do

      omega = delta / sum

C     -----------------------------------------------------------------
C     Step (2) + (3) + (4)
C     -----------------------------------------------------------------
      delta_old = delta
      delta = 0.0
      do i = 1, nh, 4
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
          do i = 1, nh
            x(i) = 1.0d0
          end do
          invert1teo = -2
          return
        end if
      else
        limit = 0
      end if 

C     -----------------------------------------------------------------
C     Step (5)
C     -----------------------------------------------------------------
      ddelta = delta / delta_old
      
      do i = 1, nh
	p(i) = r(i) + ddelta * p(i)
      end do  
      
      niter = niter + 1
      if(niter .lt. maxits) goto 10 

C     >>> maxits is reached, niter is set to -1 which means:
C     >>> "no convergence achieved"      
C     >>> return without changing starting guess xx      
13    do i = 1, nh
        x(i) = 1.0d0
      end do
      invert1teo = -1
      return

99    continue

      call mtdoe(x, t)
      do i = 1, nh
	xx(i) = x(i)
	xx(nh + i) = bb(nh + i) + t(i)
      end do

C      write(*, *) 'BiCGg5T: niter = ', niter, '  <r|g5|r> = ', delta 
      
      invert1teo = niter
      end

