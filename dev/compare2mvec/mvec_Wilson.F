      subroutine mvec_Wilson(a, b)
c       IN complex*16 a(ndrc, nsite)
c       OUT complex*16 a(ndrc, nsite)
c     /* ihip, 11 Aug 20, last modified: 11 Aug 20 */

c     ------------------------------------------------------------------
c     sigma_1 = ( 0, 1; 1, 0)
c     sigma_2 = ( 0,-i; i, 0)
c     sigma_3 = ( 1, 0; 0,-1)
c     ------------------------------------------------------------------

      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)

      complex*16 a(ndrc, nsite), b(ndrc, nsite)

      complex*16 u(nsite, ndim)
      integer ind(4, nsite)
      real*8 akap
c     ------------------------------------------------------------------
      common /gauge_fields/u
      common /index_arrays/ind
      common /kappa/akap 
c     ------------------------------------------------------------------

      complex*16 ci, cu

c     >>> complex i
      ci = dcmplx(0.0d0, 1.0d0)

      do i = 1, nsite
	
c       >>> diagonal term
        do j = 1, ndrc
	      b(j, i) = 1.0d0 * a(j, i)
	    end do

c       >>> term x + 1 (contains sigma_1 matrix)
c       b(i) = u(i, 1) * (1 + sigma_1) * a(ind(1, i))
        i1 = ind(1, i)
        b(1, i) = b(1, i) - akap * u(i, 1) * (a(1, i1) + a(2, i1))
        b(2, i) = b(2, i) - akap * u(i, 1) * (a(1, i1) + a(2, i1))

c       >>> term x - 1 (contains sigma_1 matrix)
c       b(i) = u(ind(3, i), 1)^* * (1 - sigma_1) * a(ind(3, i))
        i3 = ind(3, i)
        cu = dconjg(u(i3, 1))
        b(1, i) = b(1, i) - akap * (cu * (a(1, i3) - a(2, i3)))
        b(2, i) = b(2, i) - akap * (cu * (-a(1, i3) + a(2, i3)))        
        
c       >>> term t + 1 (contains sigma_2 matrix)
c       b(i) = u(i, 2) * (1 + sigma_2) * a(ind(2, i))
        i2 = ind(2, i)
        b(1, i) = b(1, i) - akap * u(i, 2) * (a(1, i2) - ci * a(2, i2))
        b(2, i) = b(2, i) - akap * u(i, 2) * (ci * a(1, i2) + a(2, i2))

c       >>> term t - 1 (contains sigma_2 matrix)
c       b(i) = u(ind(4, i), 2)^* * (1 - sigma_1) * a(ind(4, i))
        i4 = ind(4, i)
        cu = dconjg(u(i4, 2))
        b(1, i) = b(1, i) - akap * (cu * (a(1, i4) + ci *
     &    a(2, i4)))
        b(2, i) = b(2, i) - akap * (cu * (-ci * a(1, i4) +
     &    a(2, i4)))

      end do

      end
