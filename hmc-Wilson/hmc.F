C     ******************************************************************
      subroutine adjust_eps(eps,iacc,iter)
C     ******************************************************************
C     **                    **                                        **
C     **  ADJUST_EPS        **  Program by H. Gausterer, 16-01-1993   **
C     **                    **  (slightly modified by I. Hip, Nov 95) **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program to adjust eps (updating stepsize) to have a certain
C     value (accept)
C     ------------------------------------------------------------------
C     Part of package:   HMC2DFU1
C     ------------------------------------------------------------------
C     Input parameters:  eps, iacc, iter
C     ------------------------------------------------------------------
C     Output parameters: eps
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter(length = 30, accept = 0.8)
C     ------------------------------------------------------------------
      real*8 eps
      integer iaccept(length)
C     >>> iaccept should be static variable
      save iaccept
C     ------------------------------------------------------------------
      it = mod(iter-1,length) + 1
      iaccept(it) = iacc
C     ------------------------------------------------------------------
C     eps is changed after n = length updates
C     ------------------------------------------------------------------
      if(it.eq.length) then
C     ------------------------------------------------------------------
C     calculate the acceptance rate
C     ------------------------------------------------------------------
      amean = 0.
      do ii = 1,length
	    amean = amean + float(iaccept(ii))
      enddo
      amean = amean/length
C      write(*, *) iter, amean
C     ------------------------------------------------------------------
C     formula to perform an underrelaxed change for eps
C     ------------------------------------------------------------------
C      expo = (1. - log(amean)/log(accept))
      expo = -2.265 *log((amean-1.36)/(0.8-1.36))
      eps = eps * 2**expo
C
      endif
C
      end


C     ******************************************************************
      subroutine dsgauge(beta, u, dsg)
C     ******************************************************************
C     **                    **                                        **
C     **  DSGAUGE           **  Based on GET_DSB by C.B. Lang         **
C     **                    **  (Modified by I. Hip, Nov 95)          **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for determination of  d S_G(U) / d A
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:   beta, u(1:nsite, 1:ndim)
C     ------------------------------------------------------------------
C     Output variables:  dsg(1:nsite, 1:ndim)
C     ------------------------------------------------------------------
C     Remarks:   d S_G(U) / d A = iU d S / d U
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 u(nsite, ndim)
      real*8 beta, dsg(nsite, ndim)
      integer ind(4, nsite)
C     ------------------------------------------------------------------
      common /index_arrays/ind
C     ------------------------------------------------------------------
      do i = 1, nsite
	i1=ind(1,i)
	i2=ind(2,i)
	i3=ind(3,i)
	i4=ind(4,i)
	i5=ind(1,i4)
	i6=ind(2,i3)
       dsg(i,1) = beta *
     &     dimag(u(i,1)* (u(i1,2) * dconjg(u(i2,1) * u(i,2))
     &                   +u(i4,2) * dconjg(u(i5,2) * u(i4,1))))

       dsg(i,2) = beta *
     &     dimag(u(i,2)* (u(i3,1) * dconjg(u(i6,1) * u(i3,2))
     &                   +u(i2,1) * dconjg(u(i1,2) * u(i,1))))
      end do
      end


C     ******************************************************************
      subroutine gauss(x, n, ran)
C     ******************************************************************
C     **                    **                                        **
C     **  GAUSS             **   Based on NORM_RAN1 by H. Gausterer   **
C     **                    **     (Modified by I. Hip, Nov 95)       **
C     ******************************************************************
C     Generate initial conjugate impulses: all directions are equally
C     likely, but "lengths" vary with Gaussian distribution
C     mean value = 0, variance = 1
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:  n
C     ------------------------------------------------------------------
C     Output variables: x(1:n) (should be < x > = 0 and < x^2 > = 1)
C     ******************************************************************
      real*8 x(n)
      real ran(n)
      data twopi/6.28318530718/
C     ---------------------------------------------------------
      nh = n / 2
      call cbl_rc_log(ran, nh, 0.0)
      call cbl_rcarry(ran(nh + 1), nh, 0.0)
      do i = 1, nh
	a = sqrt(-2.0 * ran(i))
	b = twopi * ran(nh + i)
	x(i) = a * cos(b)
	x(nh + i) = a * sin(b)
      end do
      end


C     ******************************************************************
      subroutine cgauss(x, n, ran)
C     ******************************************************************
C     **                    **                                        **
C     **  CGAUSS            **   Based on NORM_RAN2 by H. Gausterer   **
C     **                    **     (Modified by I. Hip, Nov 95)       **
C     ******************************************************************
C     Generate complex Gaussian vector with mean value = 0, var = 0.5
C     ( < x1**2 + x2**2 > = 1 )
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:  n
C     ------------------------------------------------------------------
C     Output variables: x(1:2, 1:n)  ( < x > = 0, < (x^+, x) > = 1 )
C     ******************************************************************
      real*8 x(2, n)
      real ran(n, 2)
      data twopi/6.28318530718/
C     ---------------------------------------------------------
      call cbl_rc_log(ran(1, 1), n, 0.0)
      call cbl_rcarry(ran(1, 2), n, 0.0)
      do i = 1, n
	a = sqrt(-ran(i,1))
	b = twopi * ran(i,2)
	x(1, i) = a * cos(b)
	x(2, i) = a * sin(b)
      end do
      end
