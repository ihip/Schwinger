C     ******************************************************************
      subroutine v_invert(ch)
      character*32 ch
      write(ch, '(a, e7.1, a)') 'CGR - classic (', ACCURACY, ') v3'
      end
C     ******************************************************************
      integer*4 function invert(x, b)
C     ******************************************************************
C     **                    **                                        **
C     **  CGR               **    Program by      I. Hip, 26.11.1992  **
C     **                    **    Last amendment (cbl)    31.11.1992  **
C     **  v3                **    Last modified: 13 Mar 97            **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program compute x which satisfies M^+ M x = b
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:   vector x, b
C                        ATTENTION: x is start value for the iterations
C                        in order to allow improved start
C     ------------------------------------------------------------------
C     Output variables:  vector x, niter
C     ------------------------------------------------------------------
C     Remarks:- index arrays have to be initialized (call mk_index)
C             - subroutines MVEC and MTVEC should be available
C             - maybe the vectors u, r, p, s should be put into a
C               common for dummy storage
C     cf. C.B. Lang report - page 6 - based on M. Creutz's program
C     (version without 10th step of (4a) different...)
C     -----------------------------------------------------------------
C     On exit x contains (M^+ M)^(-1) b   (i.e. the inverse of M^+M)
C     and     u contains (M^+)^(-1) b     (i.e. the inverse of M^+)
C     i.e. one could add u to the (output) parameter list
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (acc=ACCURACY,maxits=1000,max_limit=100)
C     ------------------------------------------------------------------
      real*8 x(nall), b(nall)
      real*8 u(nall), r(nall), p(nall), s(nall)
      real*8 alpha, alpha_old, sum, beta, sigma, epsilon, bnorm2
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      limit = 0
      niter = 0
      call mvec(x, u)
      call mtvec(u, r)
      alpha = 0.0d0
      bnorm2 = 0.0d0
      do i = 1, nall
        bnorm2 = bnorm2 + b(i)**2
	r(i) = r(i) - b(i)
	p(i) = r(i)
	alpha = alpha + r(i)*r(i)
      end do
      epsilon = bnorm2 * acc**2

C     -----------------------------------------------------------------
C     ITERATION:
C     -----------------------------------------------------------------
C     Step (1)                        s_k = M * p_k
C     -----------------------------------------------------------------
10    call mvec(p, s)
C     -----------------------------------------------------------------
C     Step (2)                        sigma = alpha_k / (s_k, s_k)
C     -----------------------------------------------------------------
      sum =  0.
      do i = 1, nall
	sum = sum + s(i)**2
      end do
      sigma = alpha / sum
C     -----------------------------------------------------------------
C     Step (3) + (4a)                 x_k = x_(k-1) - sigma * p_(k-1)
C                                     u_k = u_(k-1) - sigma * s_k
C     -----------------------------------------------------------------
C     Creutz suggests every 10th step:  u_k = M * x_k
C     (instead of  u_k = u_(k-1) - sigma * s_k )
C     but we find no improvement for that version
C     -----------------------------------------------------------------
      do i = 1, nall
	x(i) = x(i) - sigma * p(i)
	u(i) = u(i) - sigma * s(i)
      end do
C     -----------------------------------------------------------------
C                                     Here  u contains (M^+)^(-1) b
C                                     (i.e. the inverse of M^+)
C     -----------------------------------------------------------------
C     -----------------------------------------------------------------
C     Step (4b) + (5)                 r_k = M^+ * u_k - b
C                                     alpha_k = (r_k, r_k)
C     -----------------------------------------------------------------
      call mtvec(u, r)
      alpha_old = alpha
      alpha = 0.
      do i = 1, nall
	r(i) = r(i) - b(i)
	alpha = alpha + r(i)**2
      end do
C     -----------------------------------------------------------------
C     Step (5), exit tests:
C     -----------------------------------------------------------------
C     exit if
C            Number of iterations exceed maxits
C     or     alpha < nall*acc**2
C                (alpha is residual vector squared)
C     or     after max_limit successive non-decreases of res. vector
C     -----------------------------------------------------------------
C
C     Accuracy obtained exit:
C
      if(alpha .lt. epsilon) go to 99
C
C     No convergence achieved:
C
      niter = niter + 1
      if(niter .gt. maxits) then
	invert = -1
	return
      end if
C
C     Instability exit:
C
      if(alpha .ge. alpha_old) then
	limit = limit + 1
	if(limit .gt. max_limit) then
	  invert = -2
	  return
	end if
      else
	limit = 0
      end if
C     -----------------------------------------------------------------
C     Step (6)                        beta = alpha_k / alpha_(k-1)
C     -----------------------------------------------------------------
      beta = alpha / alpha_old
C     -----------------------------------------------------------------
C     Step (7)                        p_k = r_k + beta * p_(k-1)
C     -----------------------------------------------------------------
      do i = 1, nall
	p(i) = r(i) + beta * p(i)
      end do

      goto 10

99    continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     keep the next line during the test phase:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      print *,'CGR: niter = ',niter,'  alpha = ',alpha
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     -----------------------------------------------------------------
C     On exit x contains (M^+ M)^(-1) b   (i.e. the inverse of M^+M)
C     and     u contains (M^+)^(-1) b     (i.e. the inverse of M^+)
C     -----------------------------------------------------------------
      invert = niter
      return
      end
