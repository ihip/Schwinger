C     ******************************************************************
      subroutine v_invert1(ch)
      character*32 ch
      write(ch, '(a, e7.1, a)') 'CGF1 (', ACCURACY1, ') v2'
      end
C     ******************************************************************
      integer*4 function invert1(x, b)
C     ******************************************************************
C     **                        **                                    **
C     **  INVERT1 (CGForcrand)  **    Program by I. Hip, 11 Oct 95    **
C     **  v2                    **    Last modified: 13 Feb 97        **
C     **                        **                                    **
C     ******************************************************************
C     Seems to be the simplest version of conjugate gradient: taken from
C     P. de Forcrand: "Progress on lat. QCD algorithms", hep-lat/9509082
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
      real*8 x(nall), b(nall), bb(nall)
      real*8 r(nall), p(nall)
      real*8 u(nall), t(nall)
      real*8 delta, delta_old, ddelta, omega, sum, epsilon
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      epsilon=nall*accuracy**2
      limit = 0
      niter = 0
      
C     >>> we are solving: M^+ M x = M^+ b = bb
      call mtvec(b, bb)

      call mvec(x, u)
      call mtvec(u, r)
      delta = 0.0
      do i = 1, nall
	r(i) = bb(i) - r(i)
	p(i) = r(i)
	delta = delta + r(i) * r(i)
      end do
C     -----------------------------------------------------------------
C     ITERATION:
C     -----------------------------------------------------------------
C     Step (1)         omega = delta / < p | gamma5 Q | p>              
C     -----------------------------------------------------------------
10    call mvec(p, t)
      call mtvec(t, u) 
      
      sum = 0.0
      do i = 1, nall
	sum = sum + p(i) * u(i) 
      end do

      omega = delta / sum

C     -----------------------------------------------------------------
C     Step (2) + (3) + (4)
C     -----------------------------------------------------------------
      delta_old = delta
      delta = 0.0
      do i = 1, nall
	x(i) = x(i) + omega * p(i)
	r(i) = r(i) - omega * u(i)
	delta = delta + r(i) * r(i)
      end do

c      write(*, *) niter, '. delta = ', delta

      if(delta .lt. epsilon) goto 99

      niter = niter + 1
      if(niter .gt. maxits) then
	invert1 = -1
	return
      endif

C     >>> instability exit:
      if(delta .ge. delta_old) then
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

99    continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     keep the next lines during the test phase:
c      print *,'CGF1: niter = ', niter
c      print *,'CGF1: ||r||^2 = ', delta
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     -----------------------------------------------------------------
C     On exit x contains Q^(-1) b   (i.e. the inverse of Q)
C     -----------------------------------------------------------------
      invert1 = niter
      end

