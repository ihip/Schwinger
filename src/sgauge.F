C     ******************************************************************
      real*8 function sgauge(beta, u)
C     ******************************************************************
C     **                    **                                        **
C     **  SGAUGE            **  Based on MPLAQ by C.B. Lang           **
C     **                    **  (Modified by I. Hip, Nov 95)          **
C     **                    **  Last modified: ivh, 02 Apr 97         **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for determination of gauge action (sum over plaquettes)
C
C     S_G = \beta \sum_x ( 1 - Re( U_P(x) ) )
C
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input  variables:  beta, u (gauge field)
C     ------------------------------------------------------------------
C     Output variables:  sgauge(beta, u) gauge action
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 u(nsite, ndim)
      real*8 beta, sum
      integer ind(4, nsite)
C     ------------------------------------------------------------------
C         i2 --1--  .
C         |         |
C         2         2
C         |         |
C         i  --1-- i1
C
C     ------------------------------------------------------------------
      common /index_arrays/ind
C     ------------------------------------------------------------------
      sum = 0.0d0
      do i = 1, nsite
	  i1 = ind(1, i)
	  i2 = ind(2, i)
	  sum = sum + 1.0d0 - dreal(u(i,1)*u(i1,2)*dconjg(u(i2,1)*u(i,2)))
      end do
      sgauge = beta * sum
      end
      