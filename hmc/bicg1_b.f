C     ******************************************************************
      subroutine v_invert1_b(ch)
      character*32 ch
      write(ch, '(a, e7.1, a)') 'BiCG1 (', ACCURACY1, ') v2'
      end
C     ******************************************************************
      integer*4 function invert1_b(x, b)
C     ******************************************************************
C     **                    **                                        **
C     **  INVERT1 (BiCG)    **    Program by I. Hip, 05 Sep 95        **
C     **  v2                **    Last modified: 15 Feb 97            **
C     **                    **                                        **
C     ******************************************************************
C     This implementation of BiCG algorithm is based on Templates
C     ------------------------------------------------------------------
C     Program compute x which satisfies Q^+ Q x = b
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
C     On exit x contains (Q^+ Q)^(-1) b   (i.e. the inverse of Q^+Q)
C     ****************************************************************** 
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (accuracy=ACCURACY1,maxits=1000,max_limit=100)
C     ------------------------------------------------------------------
      real*8 x(nall), b(nall)
      real*8 r1(nall), r2(nall), p1(nall), p2(nall), q1(nall), q2(nall)
      real*8 alpha, ro, ro_old, beta
      real*8 sum, epsilon, res, res_old
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      epsilon = nall * accuracy**2
      limit = 0
      niter = 0

C     initial guess x(0) is set to vector b
      do i = 1, nall
	x(i) = b(i)
      end do
	
      call mvec(x, r1)
      ro = 0.0
      do i = 1, nall
	r1(i) = b(i) - r1(i)
	r2(i) = r1(i)
	ro = ro + r1(i) * r2(i)
	p1(i) = r1(i)
	p2(i) = r1(i)
      end do
C     -----------------------------------------------------------------
C     ITERATION:
C     -----------------------------------------------------------------
C     Step (1)          q1 = Q^+ Q p1
C                       q2 = Q Q^+ p2
C     -----------------------------------------------------------------
10    call mvec(p1, q1)
      call mtvec(p2, q2)
C     -----------------------------------------------------------------
C     Step (2)          alpha = ro / (p2 . q1) 
C     -----------------------------------------------------------------
      sum = 0.0
      do i = 1, nall
	sum = sum + p2(i) * q1(i)
      end do
      alpha = ro / sum
C     -----------------------------------------------------------------
C     Step (3) +        x = x + alpha p1
C     Step (4)          r1 = r1 - alpha q1
C                       r2 = r2 - alpha q2
C     + res = r1 . r1 which is needed to test convergence      
C     -----------------------------------------------------------------
      res_old = res
      res = 0.0
      do i = 1, nall
	x(i) = x(i) + alpha * p1(i)
	r1(i) = r1(i) - alpha * q1(i)
	r2(i) = r2(i) - alpha * q2(i)
	res = res + r1(i) * r1(i)
      end do  
c      write(*, *) niter, res
C     -----------------------------------------------------------------
C     Step (5)          exit test
C     -----------------------------------------------------------------
C     Accuracy atained exit: res < nall * accuracy**2 (res = r . r)
C     -----------------------------------------------------------------
      if(res .lt. epsilon) goto 99
C     -----------------------------------------------------------------
C     Emergency exit: Number of iterations exceed MAXITS
C     -----------------------------------------------------------------
      niter = niter + 1
      if(niter .gt. maxits) then
	invert1_b = -1
	return
      endif

C     >>> instability exit:
      if(res .ge. res_old) then
	limit = limit + 1
	if(limit .gt. max_limit) then
	  invert1_b = -2
	  return
	end if
      else
	limit = 0
      end if 

C     -----------------------------------------------------------------
C     Step (6)          ro_old = ro; ro = r1 . r2; beta = ro / ro_old          
C     -----------------------------------------------------------------
      ro_old = ro
      ro = 0.0
      do i = 1, nall
	ro = ro + r1(i) * r2(i)
      end do
      beta = ro / ro_old
C     -----------------------------------------------------------------
C     Step (7)          p1 = r1 + beta p1          
C                       p2 = r2 + beta p2
C     -----------------------------------------------------------------
      do i = 1, nall
	p1(i) = r1(i) + beta * p1(i)
	p2(i) = r2(i) + beta * p2(i)
      end do
      goto 10
 
99    continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     keep the next lines during the test phase:
C      write(*, *) 'BiCG: niter = ', niter, '   res = ', res
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      invert1_b = niter
      return
      end


C     ******************************************************************
      integer*4 function invert1t_b(x, b)
C     ******************************************************************
C     **                    **                                        **
C     **  INVERT1 (BiCG)    **    Program by I. Hip, 05 Sep 95        **
C     **  v2                **    Last modified: 15 Feb 97            **
C     **                    **                                        **
C     ******************************************************************
C     This implementation of BiCG algorithm is based on Templates
C     ------------------------------------------------------------------
C     Program compute x which satisfies Q^+ Q x = b
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
C     On exit x contains (Q^+ Q)^(-1) b   (i.e. the inverse of Q^+Q)
C     ****************************************************************** 
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (accuracy=ACCURACY1,maxits=1000,max_limit=100)
C     ------------------------------------------------------------------
      real*8 x(nall), b(nall)
      real*8 r1(nall), r2(nall), p1(nall), p2(nall), q1(nall), q2(nall)
      real*8 alpha, ro, ro_old, beta
      real*8 sum, epsilon, res, res_old
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      epsilon = nall * accuracy**2
      limit = 0
      niter = 0

C     initial guess x(0) is set to vector b
      do i = 1, nall
	x(i) = b(i)
      end do
	
      call mtvec(x, r1)
      ro = 0.0
      do i = 1, nall
	r1(i) = b(i) - r1(i)
	r2(i) = r1(i)
	ro = ro + r1(i) * r2(i)
	p1(i) = r1(i)
	p2(i) = r1(i)
      end do
C     -----------------------------------------------------------------
C     ITERATION:
C     -----------------------------------------------------------------
C     Step (1)          q1 = Q^+ Q p1
C                       q2 = Q Q^+ p2
C     -----------------------------------------------------------------
10    call mtvec(p1, q1)
      call mvec(p2, q2)
C     -----------------------------------------------------------------
C     Step (2)          alpha = ro / (p2 . q1) 
C     -----------------------------------------------------------------
      sum = 0.0
      do i = 1, nall
	sum = sum + p2(i) * q1(i)
      end do
      alpha = ro / sum
C     -----------------------------------------------------------------
C     Step (3) +        x = x + alpha p1
C     Step (4)          r1 = r1 - alpha q1
C                       r2 = r2 - alpha q2
C     + res = r1 . r1 which is needed to test convergence      
C     -----------------------------------------------------------------
      res_old = res
      res = 0.0
      do i = 1, nall
	x(i) = x(i) + alpha * p1(i)
	r1(i) = r1(i) - alpha * q1(i)
	r2(i) = r2(i) - alpha * q2(i)
	res = res + r1(i) * r1(i)
      end do  
c      write(*, *) niter, res
C     -----------------------------------------------------------------
C     Step (5)          exit test
C     -----------------------------------------------------------------
C     Accuracy atained exit: res < nall * accuracy**2 (res = r . r)
C     -----------------------------------------------------------------
      if(res .lt. epsilon) goto 99
C     -----------------------------------------------------------------
C     Emergency exit: Number of iterations exceed MAXITS
C     -----------------------------------------------------------------
      niter = niter + 1
      if(niter .gt. maxits) then
	invert1t_b = -1
	return
      endif

C     >>> instability exit:
      if(res .ge. res_old) then
	limit = limit + 1
	if(limit .gt. max_limit) then
	  invert1t_b = -2
	  return
	end if
      else
	limit = 0
      end if 

C     -----------------------------------------------------------------
C     Step (6)          ro_old = ro; ro = r1 . r2; beta = ro / ro_old          
C     -----------------------------------------------------------------
      ro_old = ro
      ro = 0.0
      do i = 1, nall
	ro = ro + r1(i) * r2(i)
      end do
      beta = ro / ro_old
C     -----------------------------------------------------------------
C     Step (7)          p1 = r1 + beta p1          
C                       p2 = r2 + beta p2
C     -----------------------------------------------------------------
      do i = 1, nall
	p1(i) = r1(i) + beta * p1(i)
	p2(i) = r2(i) + beta * p2(i)
      end do
      goto 10
 
99    continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     keep the next lines during the test phase:
C      write(*, *) 'BiCG: niter = ', niter, '   res = ', res
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      invert1t_b = niter
      return
      end


