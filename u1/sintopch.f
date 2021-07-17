C     ******************************************************************
      real*8 function sintopch()
C     ******************************************************************
C     **                    **                                        **
C     **  SINTOPCH          **  Function by I. Hip, 2021-07-12        **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for determination of topological charge
C     !!! ALTERNATIVE DEFINITION with sin() - see Bardeen et al. (1998)
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input  variables:  u (gauge field) via common block
C     ------------------------------------------------------------------
C     Output: topological charge
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 u(nsite, ndim)
      integer ind(4, nsite)
      complex*16 z
      real*8 pi, getphi, tpc
C     ------------------------------------------------------------------
      common /gauge_fields/u
      common /index_arrays/ind
C     ------------------------------------------------------------------
      getphi(z) = datan2(dimag(z), dreal(z))
      pi = datan2(0.0d0, -1.0d0)

      tpc = 0.0d0
      do i = 1, nsite
        i1 = ind(1, i)
        i2 = ind(2, i)
        chi = (getphi(u(i, 1)) + getphi(u(i1, 2)) -
     &        getphi(u(i2, 1)) - getphi(u(i, 2)))
        tpc = tpc + sin(chi)
      end do
      sintopch = tpc / (2.0d0 * pi)
      end
