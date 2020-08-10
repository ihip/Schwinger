C     ******************************************************************
      subroutine mvec(a,b)
C     ******************************************************************
C     **                    **                                        **
C     **  MVEC              **    Program by C.B. Lang, 20-12-1992    **
C     **                    **    Last Amendment (cbl), 12-01-1993    **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for multiplication b := M a
C     d=2
C     gauge group: U(1)
C     Wilson fermions (r=1)
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:   vector a
C     ------------------------------------------------------------------
C     Output variables:  vector b
C     ------------------------------------------------------------------
C     Remarks: index arrays have to be initialized (call mk_index)
C              eventual antiperiodic b.c. should be already
C              applied on the u-fields
C     ------------------------------------------------------------------
C     sigma_1 = ( 0, 1; 1, 0)
C     sigma_2 = ( 0,-i; i, 0)
C     sigma_3 = ( 1, 0; 0,-1)
C     ------------------------------------------------------------------
C     index = (1:nall) = (1:ndrc,1:ngrp,1:NTIME,1:NSPACE)
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      real*8 a(2,ndrc,nsite),b(2,ndrc,nsite),u(2,nsite,ndim)
      real*8 xr,xi,c1r,c1i,c2r,c2i,c3r,c3i,c4r,c4i
      real*8 akap
      integer ind(4,nsite)
C     ------------------------------------------------------------------
      common /gauge_fields/u
      common /index_arrays/ind
      common /kappa/akap 
C     ------------------------------------------------------------------
      do i=1,nsite

	    i1 = ind(1,i)
	    i2 = ind(2,i)
	    i3 = ind(3,i)
	    i4 = ind(4,i)

	    xr  = a(1,1,i1) + a(1,2,i1)
	    xi  = a(2,1,i1) + a(2,2,i1)
	    c1r = u(1,i,1)*xr - u(2,i,1)*xi
	    c1i = u(1,i,1)*xi + u(2,i,1)*xr

	    xr  = a(1,1,i2) + a(2,2,i2)
	    xi  = a(2,1,i2) - a(1,2,i2)
	    c2r = u(1,i,2)*xr - u(2,i,2)*xi
	    c2i = u(1,i,2)*xi + u(2,i,2)*xr

	    xr  = a(1,1,i3) - a(1,2,i3)
	    xi  = a(2,1,i3) - a(2,2,i3)
	    c3r = u(1,i3,1)*xr + u(2,i3,1)*xi
	    c3i = u(1,i3,1)*xi - u(2,i3,1)*xr

	    xr  = a(1,1,i4) - a(2,2,i4)
	    xi  = a(2,1,i4) + a(1,2,i4)
	    c4r = u(1,i4,2)*xr + u(2,i4,2)*xi
	    c4i = u(1,i4,2)*xi - u(2,i4,2)*xr

	    b(1,1,i) = a(1,1,i) - akap*( c1r + c2r + c3r + c4r)
	    b(2,1,i) = a(2,1,i) - akap*( c1i + c2i + c3i + c4i)
	    b(1,2,i) = a(1,2,i) - akap*( c1r - c3r - c2i + c4i)
	    b(2,2,i) = a(2,2,i) - akap*( c1i - c3i + c2r - c4r)

      enddo
      end

C     ******************************************************************
      subroutine mtvec(a,b)
C     ******************************************************************
C     **                    **                                        **
C     **  MTVEC             **    Program by C.B. Lang, 20-12-1992    **
C     **                    **    Last Amendment (cbl), 12-01-1993    **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for multiplication b := (M^+) a
C     d=2
C     gauge group: U(1)
C     Wilson fermions (r=1)
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:   vector a
C     ------------------------------------------------------------------
C     Output variables:  vector b
C     ------------------------------------------------------------------
C     Remarks: index arrays have to be initialized (call mk_index)
C     ------------------------------------------------------------------
C     sigma_1 = ( 0, 1; 1, 0)
C     sigma_2 = ( 0,-i; i, 0)
C     sigma_3 = ( 1, 0; 0,-1)
C     ------------------------------------------------------------------
C     index = (1:nall) = (1:ndrc,1:ngrp,1:NTIME,1:NSPACE)
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      real*8 a(2,ndrc,nsite),b(2,ndrc,nsite),u(2,nsite,ndim)
      real*8 xr,xi,c1r,c1i,c2r,c2i,c3r,c3i,c4r,c4i
      real*8 akap
      integer ind(4,nsite)
C     ------------------------------------------------------------------
      common /gauge_fields/u
      common /index_arrays/ind
      common /kappa/akap 
C     ------------------------------------------------------------------
      do i=1,nsite

        i1 =ind(1,i)
        i2 =ind(2,i)
        i3 =ind(3,i)
        i4 =ind(4,i)

        xr  = a(1,1,i1) - a(1,2,i1)
        xi  = a(2,1,i1) - a(2,2,i1)
        c1r = u(1,i,1)*xr - u(2,i,1)*xi
        c1i = u(1,i,1)*xi + u(2,i,1)*xr

        xr  = a(1,1,i2) - a(2,2,i2)
        xi  = a(2,1,i2) + a(1,2,i2)
        c2r = u(1,i,2)*xr - u(2,i,2)*xi
        c2i = u(1,i,2)*xi + u(2,i,2)*xr

        xr  = a(1,1,i3) + a(1,2,i3)
        xi  = a(2,1,i3) + a(2,2,i3)
        c3r = u(1,i3,1)*xr + u(2,i3,1)*xi
        c3i = u(1,i3,1)*xi - u(2,i3,1)*xr

        xr  = a(1,1,i4) + a(2,2,i4)
        xi  = a(2,1,i4) - a(1,2,i4)
        c4r = u(1,i4,2)*xr + u(2,i4,2)*xi
        c4i = u(1,i4,2)*xi - u(2,i4,2)*xr

        b(1,1,i) = a(1,1,i) - akap*( c1r + c2r + c3r + c4r)
        b(2,1,i) = a(2,1,i) - akap*( c1i + c2i + c3i + c4i)
        b(1,2,i) = a(1,2,i) - akap*(-c1r + c3r + c2i - c4i)
        b(2,2,i) = a(2,2,i) - akap*(-c1i + c3i - c2r + c4r)

      enddo
      end