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
C     ******************************************************************
C     Standard subroutines
C     ******************************************************************
C     ******************************************************************

      real*8 function norm(v, length)
      real*8 v(*), sum
      sum = 0.0
      do i = 1, length
	sum = sum + v(i)**2
      end do
      norm = sum
      end

      complex*16 function scalp(a, b, length)
      complex*16 a(*), b(*), sum
      sum = dcmplx(0.0, 0.0)
      do i = 1, length
	sum = sum + dconjg(a(i)) * b(i)
      end do
      scalp = sum
      end

C     ******************************************************************
      subroutine create_u(istart, u, input_file)
C     ******************************************************************
C     **                    **                                        **
C     **  CREATE_U          **  Subroutine by I. Hip, Nov 95          **
C     **                    **  Last modified: 18 Feb 96              **
C     ******************************************************************
C     Create initial gauge fields u
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input  variables:  istart (Start: 0 = cold, 1 = hot, 2 = file)
C                        input_file (file with saved configuration)
C     ------------------------------------------------------------------
C     Output variables:  u (gauge field)
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 u(nsite, ndim)
      character*32 input_file
      real rand(nsite, ndim)
      real*8 pi2, alpha
      real*8 beta, akap, eps
C     ------------------------------------------------------------------
     
C     >>> cold start
      if(istart .eq. 0) then
	do nd = 1, ndim
	  do ns = 1, nsite
	    u(ns, nd) = dcmplx(1.0d0, 0.0d0)
	  end do
	end do
	return
      end if             

C     >>> hot start
      if(istart .eq. 1) then
	pi2 = 2.0d0 * datan2(0.0d0, -1.0d0)
	call cbl_rcarry(rand, ngaug, 0.0)
	do nd = 1, ndim
	  do ns = 1, nsite
	    alpha = pi2 * rand(ns, nd)
	    u(ns, nd) = dcmplx(dcos(alpha), dsin(alpha))
	  end do
	end do
	return
      end if

C     >>> start from saved configuration
      if(istart .eq. 2) then
	open(1, file = input_file, form = 'unformatted', status = 'old')
	read(1) ntime, nspace
        read(1) beta, akap, eps
	read(1) it
	read(1) u
	close(1)
	return
      endif
      
      stop 'Invalid istart parameter'
      end

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

C     ******************************************************************
      subroutine mk_index(ind)
C     ******************************************************************
C     **                    **                                        **
C     **  MK_INDEX          **    Program by C.B. Lang, 02-03-1993    **
C     **                    **    Debugged by I. Hip, 26 Nov 96       **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for construction of d=2 index vectors
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Output variables:   ind(1:4,1:nsite)
C     ------------------------------------------------------------------
C     Remarks: index (1:nsite) = (1:NSPACE,1:NTIME)
C                            i = (it-1)*NSPACE+ix
C    ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      integer ind(4,nsite)
C     ------------------------------------------------------------------
C     first the n.n. index arrays
C     ------------------------------------------------------------------
      i=0
      do it=1,NTIME
	itp=mod(it,NTIME)+1
	itm=mod(it+NTIME-2,NTIME)+1
	do ix=1,NSPACE
	  i=i+1
	  ixp=mod(ix,NSPACE)+1
	  ixm=mod(ix+NSPACE-2,NSPACE)+1
	  ind(1,i)=(it-1)*NSPACE+ixp
	  ind(2,i)=(itp-1)*NSPACE+ix
	  ind(3,i)=(it-1)*NSPACE+ixm
	  ind(4,i)=(itm-1)*NSPACE+ix
	enddo
      enddo
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
      subroutine dsferm(akap, u, a, dsf)
C     ******************************************************************
C     **                    **                                        **
C     **  DSFERM            **  Based on GET_DSF by C.B. Lang         **
C     **                    **  (Modified by I. Hip, Nov 95)          **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for determination of  d S_F(U) / d A
C     Wilson fermions (r = 1)
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:  akap, u, a = (M^+ M)^(-1) phi
C     ------------------------------------------------------------------
C     Output variables:  dsf(1:nsite,1:ndim)
C     ------------------------------------------------------------------
C     Remarks:   d S_F(U) / d U =
C                - (a, [(iU dM^+/dU)M + M^+(iU dM/dU)] a) =
C                    (* where  (a = (M^+ M)^(-1) phi) *)
C                = - 2 Re (a, (i U dM^+/dU)  M a) = (* b = M a *)
C                = - 2 Re (a, (i U dM^+/dU)  b)
C     ------------------------------------------------------------------
C     sigma_1 = ( 0, 1; 1, 0)
C     sigma_2 = ( 0,-i; i, 0)
C     sigma_3 = ( 1, 0; 0,-1)
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      real*8 u(2,nsite,ndim), a(2,ndrc,nsite), b(2,ndrc,nsite)
      real*8 dsf(nsite,ndim), akap, akap2
      integer ind(4,nsite)
C     ------------------------------------------------------------------
      common /index_arrays/ind
C     ------------------------------------------------------------------
      call mvec(a, b)
      akap2 = 2.0d0 * akap

      do i = 1, nsite
	i1=ind(1,i)
	i2=ind(2,i)
	dsf(i,1) = -akap2 * (
     &   ((a(1,1,i) -a(1,2,i) )*(b(1,1,i1)-b(1,2,i1))
     &   +(a(2,1,i) -a(2,2,i) )*(b(2,1,i1)-b(2,2,i1))
     &   +(a(1,1,i1)+a(1,2,i1))*(b(1,1,i) +b(1,2,i) )
     &   +(a(2,1,i1)+a(2,2,i1))*(b(2,1,i) +b(2,2,i) ))*u(2,i,1)
     &   +((a(1,1,i) -a(1,2,i) )*(b(2,1,i1)-b(2,2,i1))
     &   -(a(2,1,i) -a(2,2,i) )*(b(1,1,i1)-b(1,2,i1))
     &   -(a(1,1,i1)+a(1,2,i1))*(b(2,1,i) +b(2,2,i) )
     &   +(a(2,1,i1)+a(2,2,i1))*(b(1,1,i) +b(1,2,i) ))*u(1,i,1) )
	dsf(i,2) = -akap2 * (
     &   ((a(1,1,i) -a(2,2,i) )*(b(1,1,i2)-b(2,2,i2))
     &   +(a(1,2,i) +a(2,1,i) )*(b(2,1,i2)+b(1,2,i2))
     &   +(a(1,1,i2)+a(2,2,i2))*(b(1,1,i) +b(2,2,i) )
     &   -(a(1,2,i2)-a(2,1,i2))*(b(2,1,i) -b(1,2,i) ))*u(2,i,2)
     &   +((a(1,1,i) -a(2,2,i) )*(b(2,1,i2)+b(1,2,i2))
     &   -(a(1,2,i) +a(2,1,i) )*(b(1,1,i2)-b(2,2,i2))
     &   -(a(1,1,i2)+a(2,2,i2))*(b(2,1,i) -b(1,2,i) )
     &   -(a(1,2,i2)-a(2,1,i2))*(b(1,1,i) +b(2,2,i) ))*u(1,i,2) )
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

C     ******************************************************************
      subroutine unitarize(u, n)
C     ******************************************************************
C     **                    **                                        **
C     **  UNITARIZE         **    Program by C.B. Lang, 02-03-1993    **
C     **                    **     (Modified by I. Hip, Nov 1995)     **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program unitarizes u-fields
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:   u - gauge field; n - length (usually = ngaug)
C     ------------------------------------------------------------------
C     Output variables:  u
C     ------------------------------------------------------------------
C     ******************************************************************
      complex*16 u(*)
C     ------------------------------------------------------------------
      do i = 1, n
	u(i) = u(i) / cdabs(u(i))
      end do
      end

C     ******************************************************************
C     ******************************************************************
C     End of standard subroutines
C     ******************************************************************
C     ******************************************************************
