C !!! CLOVER COMMENTED OUT !!! - Nov 2007

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Collection of routines which implement hypercube fermions:           C
C                                                                      C
C init_hf - initialization of parameters                               C
C                                                                      C
C make_ihf - initialization of necessary index array                   C
C                                                                      C
C mvec_dhf - implementation of the matrix-vector multiplication        C
C                                                                      C
C mhvec_dhf - hermitian conjugate matrix times vector                  C
C                                                                      C
C fill_dhf - construct the dense Dirac operator matrix                 C
C                                                                      C
C fill_dhc - construct the dense matrix for 1st order corrected hf     C
C                                                                      C
C hip, 12. Dec '98, ITP Graz; last modified: 07. Feb '99               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     ******************************************************************
      subroutine init_hf(ferm_name)
C     ******************************************************************
C     IN character*16 ferm_name - fermion action name
C     OUT (through COMMON!) real*8 lambda_1, lambda_2, rho_2, w_l, w_r -
C       parameters which define hypercube fermions
C     ------------------------------------------------------------------
C     WARNING: be careful to call init_hf with string of length 16
C              (put additional blanks if necessary)
C     ------------------------------------------------------------------
C     hip, 06. Jan '99, ITP Graz; last modified: 08. Jan '99
C     ------------------------------------------------------------------
      character*16 ferm_name
      real*8 lambda_1, lambda_2, rho_2, w_l, w_r
C     ------------------------------------------------------------------
      common /ihf_param/ lambda_1, lambda_2, rho_2, w_l, w_r
C     ------------------------------------------------------------------

      write(*, *)
      if(ferm_name .eq. 'TruncPerfect') then 
        lambda_1 = -0.24477247651781d0
        lambda_2 = -0.12761376173876d0
        rho_2 =     9.53057717269d-02
        w_l =       0.20913880d0
        w_r =       0.05330606d0
	write(*, *) 'init_hf: TruncPerfect'
      else if(ferm_name .eq. 'TruncPerfect_ws') then 
        lambda_1 = -0.24477247651781d0
        lambda_2 = -0.12761376173876d0
        rho_2 =     9.53057717269d-02
        w_l =       0.0d0
        w_r =       0.0d0
	write(*, *) 'init_hf: TruncPerfect_ws'
      else if(ferm_name .eq. 'GWRoptimal') then
        lambda_1 = -0.24771369d0
        lambda_2 = -0.12501304d0
        rho_2 =     0.09708390d0
        w_l =       0.21137856d0
        w_r =       0.04989817d0
	write(*, *) 'init_hf: GWRoptimal'
      else if(ferm_name .eq. 'GWRoptimal_ws') then
        lambda_1 = -0.24771369d0
        lambda_2 = -0.12501304d0
        rho_2 =     0.09708390d0
        w_l =       0.0d0
        w_r =       0.0d0
	write(*, *) 'init_hf: GWRoptimal_ws'
      else if(ferm_name .eq. 'TestFermion1') then
        lambda_1 = -0.25d0
        lambda_2 = -0.125d0
        rho_2 =     0.083d0
        w_l =       0.21504235d0
        w_r =       0.04035802d0
	write(*, *) 'init_hf: TestFermion1'
      else if(ferm_name .eq. 'TestFermion1_ws') then
        lambda_1 = -0.25d0
        lambda_2 = -0.125d0
        rho_2 =     0.083d0
        w_l =       0.0d0
        w_r =       0.0d0
	write(*, *) 'init_hf: TestFermion1_ws'
      else
        write(*, *) ferm_name, 'not yet implemented...'
        stop
      end if
      write(*, *)

      return
      end

C     ******************************************************************
      subroutine set_staples(a, b)
C     ******************************************************************
      real*8 a, b
      real*8 lambda_1, lambda_2, rho_2, w_l, w_r
C     ------------------------------------------------------------------
      common /ihf_param/ lambda_1, lambda_2, rho_2, w_l, w_r
C     ------------------------------------------------------------------

      w_l = a
      w_r = b

      return
      end


C     ******************************************************************
      subroutine make_ihf()
C     ******************************************************************
C     - index table for hypercube fermions (8 nearest neighbours)
C
C         6    2    5
C
C         3         1
C
C         7    4    8
C
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      integer*4 ihf(8, nsite)
C     ------------------------------------------------------------------
      common /ihf_index/ ihf
C     ------------------------------------------------------------------
      
      i = 0
      do it = 1, NTIME
	itp = mod(it, NTIME) + 1
	itm = mod(it + NTIME - 2, NTIME) + 1
	do ix = 1, NSPACE
	  i = i + 1
	  ixp = mod(ix, NSPACE) + 1
	  ixm = mod(ix + NSPACE - 2, NSPACE) + 1
	  ihf(1, i) = (it - 1) * NSPACE + ixp
	  ihf(2, i) = (itp - 1) * NSPACE + ix
	  ihf(3, i) = (it - 1) * NSPACE + ixm
	  ihf(4, i) = (itm - 1) * NSPACE + ix
	  ihf(5, i) = (itp - 1) * NSPACE + ixp
	  ihf(6, i) = (itp - 1) * NSPACE + ixm
	  ihf(7, i) = (itm - 1) * NSPACE + ixm
	  ihf(8, i) = (itm - 1) * NSPACE + ixp
	enddo
      enddo
      end


C     ******************************************************************
      subroutine mvec_dhf(a, b)
C     ******************************************************************
C     b = D_hf a, where D_hf is hypercube fermion operator (Bietenholz)
C     (this routine has been checked on 12. Dec '98 by WB & IH)
C     >>> coefficients from letter dated on 27. Oct '98 <<<
C     >>> new parametrization from letter dated on 18. Dec '98 <<<  
C     ------------------------------------------------------------------
C     !! DON'T FORGET TO INITIALIZE PARAMETERS BY CALLING init_hf()   !!
C     !! DON'T FORGET TO INITIALIZE INDEX TABLE BY CALLING make_ihf() !!
C     ------------------------------------------------------------------
C     IN complex*16 a(ndrc, nsite) == a(nferm) - vector to be
C       multiplied by the Dirac operator matrix
C     OUT complex*16 b(ndrc, nsite) == b(nferm) - result: b = D_hf a
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm = ndrc*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      integer*4 ihf(8, nsite)
      complex*16 u(nsite, ndim)

      real*8 lambda_1, lambda_2, rho_2, w_l, w_r
C     ------------------------------------------------------------------
      common /gauge_fields/ u
      common /ihf_index/ ihf
      common /ihf_param/ lambda_1, lambda_2, rho_2, w_l, w_r
C     ------------------------------------------------------------------
      complex*16 a(ndrc, nsite), b(ndrc, nsite)
      real*8 lambda, rho, f
      complex*16 ic
      real*8 rho_1

      rho_1 = 0.5d0 - 2.0d0 * rho_2

      ic = dcmplx(0.0d0, 1.0d0)

      do i = 1, nsite

C       >>> diagonal term
C	lambda = 1.48954495
	lambda = -4.0d0 * (lambda_1 + lambda_2)

C	write(*, *) 'diag: lambda = ', lambda
	
	do j = 1, ndrc
	  b(j, i) = lambda * a(j, i)
	end do

C       >>> 1-link paths
C	lambda = -0.14238963
C	rho = 0.27640390
        lambda = lambda_1 * (1.0d0 - 2.0d0 * w_l)
	rho = rho_1 * (1.0d0 - 2.0d0 * w_r)

C	write(*, *) '1-link: lambda = ', lambda
C	write(*, *) '        rho =    ', rho

	do j = 1, ndrc
	  b(j, i) = b(j, i) + lambda * (
     &      u(i, 1) * a(j, ihf(1, i)) + 
     &      u(i, 2) * a(j, ihf(2, i)) +
     &      dconjg(u(ihf(3, i), 1)) * a(j, ihf(3, i)) +
     &      dconjg(u(ihf(4, i), 2)) * a(j, ihf(4, i)) )
	end do

	b(1, i) = b(1, i) + rho * (
     &    u(i, 1) * a(2, ihf(1, i)) + 
     &    u(i, 2) * (-ic) * a(2, ihf(2, i)) -
     &    dconjg(u(ihf(3, i), 1)) * a(2, ihf(3, i)) -
     &    dconjg(u(ihf(4, i), 2)) * (-ic) * a(2, ihf(4, i)) )
      
	b(2, i) = b(2, i) + rho * (
     &    u(i, 1) * a(1, ihf(1, i)) + 
     &    u(i, 2) * ic * a(1, ihf(2, i)) -
     &    dconjg(u(ihf(3, i), 1)) * a(1, ihf(3, i)) -
     &    dconjg(u(ihf(4, i), 2)) * ic * a(1, ihf(4, i)) )
      
C       >>> 3-link paths corresponding to 1-link paths
C	lambda = -0.05119142
C	rho = 0.01649228
	lambda = lambda_1 * w_l
	rho = rho_1 * w_r

C	write(*, *) '3-link: lambda = ', lambda
C	write(*, *) '        rho =    ', rho

	do j = 1, ndrc
	  b(j, i) = b(j, i) + lambda * (
     &      (u(i, 2) * u(ihf(2, i), 1) * dconjg(u(ihf(1, i), 2)) +
     &      dconjg(u(ihf(4, i), 2)) * u(ihf(4, i), 1) * u(ihf(8, i), 2)) 
     &      * a(j, ihf(1, i)) +
     &      (u(i, 1) * u(ihf(1, i), 2) * dconjg(u(ihf(2, i), 1)) +
     &      dconjg(u(ihf(3, i), 1)) * u(ihf(3, i), 2) * u(ihf(6, i), 1))
     &      * a(j, ihf(2, i)) +
     &      (u(i, 2) * dconjg(u(ihf(6, i), 1)) * dconjg(u(ihf(3, i), 2))
     &      + dconjg(u(ihf(4, i), 2)) * dconjg(u(ihf(7, i), 1)) * 
     &      u(ihf(7, i), 2)) * a(j, ihf(3, i)) +
     &      (u(i, 1) * dconjg(u(ihf(8, i), 2)) * dconjg(u(ihf(4, i), 1))
     &      + dconjg(u(ihf(3, i), 1)) * dconjg(u(ihf(7, i), 2)) * 
     &      u(ihf(7, i), 1)) * a(j, ihf(4, i)) )
	end do

	b(1, i) = b(1, i) + rho * (
     &    (u(i, 2) * u(ihf(2, i), 1) * dconjg(u(ihf(1, i), 2)) +
     &    dconjg(u(ihf(4, i), 2)) * u(ihf(4, i), 1) * u(ihf(8, i), 2)) 
     &    * a(2, ihf(1, i)) +
     &    (u(i, 1) * u(ihf(1, i), 2) * dconjg(u(ihf(2, i), 1)) +
     &    dconjg(u(ihf(3, i), 1)) * u(ihf(3, i), 2) * u(ihf(6, i), 1))
     &    * (-ic) * a(2, ihf(2, i)) +
     &    (u(i, 2) * dconjg(u(ihf(6, i), 1)) * dconjg(u(ihf(3, i), 2))
     &    + dconjg(u(ihf(4, i), 2)) * dconjg(u(ihf(7, i), 1)) * 
     &    u(ihf(7, i), 2)) * (-1) * a(2, ihf(3, i)) +
     &    (u(i, 1) * dconjg(u(ihf(8, i), 2)) * dconjg(u(ihf(4, i), 1))
     &    + dconjg(u(ihf(3, i), 1)) * dconjg(u(ihf(7, i), 2)) * 
     &    u(ihf(7, i), 1)) * ic * a(2, ihf(4, i)) )

	b(2, i) = b(2, i) + rho * (
     &    (u(i, 2) * u(ihf(2, i), 1) * dconjg(u(ihf(1, i), 2)) +
     &    dconjg(u(ihf(4, i), 2)) * u(ihf(4, i), 1) * u(ihf(8, i), 2)) 
     &    * a(1, ihf(1, i)) +
     &    (u(i, 1) * u(ihf(1, i), 2) * dconjg(u(ihf(2, i), 1)) +
     &    dconjg(u(ihf(3, i), 1)) * u(ihf(3, i), 2) * u(ihf(6, i), 1))
     &    * ic * a(1, ihf(2, i)) +
     &    (u(i, 2) * dconjg(u(ihf(6, i), 1)) * dconjg(u(ihf(3, i), 2))
     &    + dconjg(u(ihf(4, i), 2)) * dconjg(u(ihf(7, i), 1)) * 
     &    u(ihf(7, i), 2)) * (-1) * a(1, ihf(3, i)) +
     &    (u(i, 1) * dconjg(u(ihf(8, i), 2)) * dconjg(u(ihf(4, i), 1))
     &    + dconjg(u(ihf(3, i), 1)) * dconjg(u(ihf(7, i), 2)) * 
     &    u(ihf(7, i), 1)) * (-ic) * a(1, ihf(4, i)) )


C       >>> 2-link paths to the nearest diagonal neighbours
C	lambda = -0.06380688
C	rho = 0.04765289

	lambda = lambda_2 / 2.0d0
	rho = rho_2 / 2.0d0

C	write(*, *) '2-link: lambda = ', lambda
C	write(*, *) '        rho =    ', rho

	do j = 1, ndrc
	  b(j, i) = b(j, i) + lambda * (
     &      (u(i, 1) * u(ihf(1, i), 2) + u(i, 2) * u(ihf(2, i), 1)) *
     &      a(j, ihf(5, i)) +
     &      (u(i, 2) * dconjg(u(ihf(6, i), 1)) + 
     &      dconjg(u(ihf(3, i), 1)) * u(ihf(3, i), 2)) * 
     &      a(j, ihf(6, i)) +
     &      (dconjg(u(ihf(3, i), 1)) * dconjg(u(ihf(7, i), 2)) +
     &      dconjg(u(ihf(4, i), 2)) * dconjg(u(ihf(7, i), 1))) *
     &      a(j, ihf(7, i)) +
     &      (u(i, 1) * dconjg(u(ihf(8, i), 2)) +
     &      dconjg(u(ihf(4, i), 2)) * u(ihf(4, i), 1)) *
     &      a(j, ihf(8, i)) )
	end do

	b(1, i) = b(1, i) + rho * (
     &    (u(i, 1) * u(ihf(1, i), 2) + u(i, 2) * u(ihf(2, i), 1)) *
     &    (1 - ic) * a(2, ihf(5, i)) +
     &    (u(i, 2) * dconjg(u(ihf(6, i), 1)) + 
     &    dconjg(u(ihf(3, i), 1)) * u(ihf(3, i), 2)) * 
     &    (-1 - ic) * a(2, ihf(6, i)) +
     &    (dconjg(u(ihf(3, i), 1)) * dconjg(u(ihf(7, i), 2)) +
     &    dconjg(u(ihf(4, i), 2)) * dconjg(u(ihf(7, i), 1))) *
     &    (-1 + ic) * a(2, ihf(7, i)) +
     &    (u(i, 1) * dconjg(u(ihf(8, i), 2)) +
     &    dconjg(u(ihf(4, i), 2)) * u(ihf(4, i), 1)) *
     &    (1 + ic) * a(2, ihf(8, i)) )

	b(2, i) = b(2, i) + rho * (
     &    (u(i, 1) * u(ihf(1, i), 2) + u(i, 2) * u(ihf(2, i), 1)) *
     &    (1 + ic) * a(1, ihf(5, i)) +
     &    (u(i, 2) * dconjg(u(ihf(6, i), 1)) + 
     &    dconjg(u(ihf(3, i), 1)) * u(ihf(3, i), 2)) * 
     &    (-1 + ic) * a(1, ihf(6, i)) +
     &    (dconjg(u(ihf(3, i), 1)) * dconjg(u(ihf(7, i), 2)) +
     &    dconjg(u(ihf(4, i), 2)) * dconjg(u(ihf(7, i), 1))) *
     &    (-1 - ic) * a(1, ihf(7, i)) +
     &    (u(i, 1) * dconjg(u(ihf(8, i), 2)) +
     &    dconjg(u(ihf(4, i), 2)) * u(ihf(4, i), 1)) *
     &    (1 - ic) * a(1, ihf(8, i)) )

C !!! CLOVER COMMENTED OUT !!!
C	>>> clover term
c	f = dimag(
c     &    u(i, 1) * u(ihf(1, i), 2) *
c     &    dconjg(u(ihf(2, i), 1)) * dconjg(u(i, 2)) +
c     &    u(ihf(3, i), 1) * u(i, 2) *
c     &    dconjg(u(ihf(6, i), 1)) * dconjg(u(ihf(3, i), 2)) +
c     &    u(ihf(7, i), 1) * u(ihf(4, i), 2) *
c     &    dconjg(u(ihf(3, i), 1)) * dconjg(u(ihf(7, i), 2)) +
c     &    u(ihf(4, i), 1) * u(ihf(8, i), 2) *
c     &    dconjg(u(i, 1)) * dconjg(u(ihf(4, i), 2)) ) / 4.0d0

c	b(1, i) = b(1, i) + f * a(1, i) / 2.0d0
c       b(2, i) = b(2, i) - f * a(2, i) / 2.0d0

      end do        

      end


C     ******************************************************************
      subroutine mhvec_dhf(a, b)
C     ******************************************************************
C     b = (D_hf)^+ a, where (D_hf)^+ is hermitian conjugate of D_hf
C     - the gamma_5 hermiticity relation is used: D^+ = g_5 D g_5
C     ------------------------------------------------------------------
C     !! DON'T FORGET TO INITIALIZE INDEX TABLE BY CALLING make_ihf() !!
C     ------------------------------------------------------------------
C     IN complex*16 a(ndrc, nsite) == a(nferm) - vector to be
C       multiplied by the Dirac operator matrix
C     OUT complex*16 b(ndrc, nsite) == b(nferm) - result: b = (D_hf)^+ a
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm = ndrc*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      integer*4 ihf(8, nsite)
      complex*16 u(nsite, ndim)
C     ------------------------------------------------------------------
      common /gauge_fields/ u
      common /ihf_index/ ihf
C     ------------------------------------------------------------------
      complex*16 a(ndrc, nsite), b(ndrc, nsite)
C     >>> tmp is used to keep a unchanged
      complex*16 tmp(ndrc, nsite)

C     >>> multiply with gamma_5
      do i = 1, nsite
	tmp(1, i) = a(1, i)
        tmp(2, i) = -a(2, i)
      end do

C     >>> multiply with D_hf
      call mvec_dhf(tmp, b)

C     >>> multiply with g5 - just "repair" every second component of b
      do i = 1, nsite
C       b(1, i) = b(1, i) <<< this is not necessary :)
        b(2, i) = -b(2, i)
      end do
        
      end 


C     ******************************************************************
      subroutine mvec_dha(mu, a, b)
C     ******************************************************************
C     ------------------------------------------------------------------
C     !! DON'T FORGET TO INITIALIZE INDEX TABLE BY CALLING make_ihf() !!
C     ------------------------------------------------------------------
C     IN complex*16 a(ndrc, nsite) == a(nferm) - vector to be
C       multiplied by the Dirac operator matrix
C     OUT complex*16 b(ndrc, nsite) == b(nferm) - result: b = (D_hf)^+ a
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm = ndrc*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      integer*4 ihf(8, nsite)
      complex*16 u(nsite, ndim)
C     ------------------------------------------------------------------
      common /gauge_fields/ u
      common /ihf_index/ ihf
C     ------------------------------------------------------------------
      real*8 mu
      complex*16 a(nferm), b(nferm)

C     >>> multiply with D_hf
      call mvec_dhf(a, b)

C     >>> compute \mu - b
      do i = 1, nferm
	b(i) = mu * a(i) - b(i)
      end do
        
      end 


C     ******************************************************************
      subroutine mhvec_dha(mu, a, b)
C     ******************************************************************
C     ------------------------------------------------------------------
C     !! DON'T FORGET TO INITIALIZE INDEX TABLE BY CALLING make_ihf() !!
C     ------------------------------------------------------------------
C     IN complex*16 a(ndrc, nsite) == a(nferm) - vector to be
C       multiplied by the Dirac operator matrix
C     OUT complex*16 b(ndrc, nsite) == b(nferm) - result: b = (D_hf)^+ a
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm = ndrc*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      integer*4 ihf(8, nsite)
      complex*16 u(nsite, ndim)
C     ------------------------------------------------------------------
      common /gauge_fields/ u
      common /ihf_index/ ihf
C     ------------------------------------------------------------------
      real*8 mu
      complex*16 a(nferm), b(nferm)

C     >>> multiply with D_hf
      call mhvec_dhf(a, b)

C     >>> compute \mu - b
      do i = 1, nferm
	b(i) = mu * a(i) - b(i)
      end do
        
      end 


C     ******************************************************************
      subroutine fill_dhf(a)
C     ******************************************************************
C     OUT complex*16 a(nferm, nferm) - Dirac operator matrix
C     ------------------------------------------------------------------
	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)
C       ----------------------------------------------------------------
	complex*16 a(nferm, nferm)
	complex*16 s(nferm)
C       ----------------------------------------------------------------

c       >>> fill in 'dense' hfpercube fermions
c       (this can be done faster...)
	do i = 1, nferm
	  s(i) = dcmplx(0.0d0, 0.0d0)
	end do
	do k = 1, nferm
	  s(k) = dcmplx(1.0d0, 0.0d0)
	  call mvec_dhf(s, a(1, k))
	  s(k) = dcmplx(0.0d0, 0.0d0)
	end do

	return
	end


C     ******************************************************************
      subroutine fill_dhc1(mu, a)
C     ******************************************************************
C     dhc: corrected hypercube fermions: expression for correction
C          taken from draft paper version 03 Feb 99
C     ------------------------------------------------------------------
C     hip, Jan 99; last modified: 04 Feb 99
C     ------------------------------------------------------------------
C     IN real*8 mu - shift of the origin; radius of the projected circle
C     OUT complex*16 a(nferm, nferm) - Dirac operator matrix
C     ------------------------------------------------------------------
	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)
C       ----------------------------------------------------------------
	real*8 mu, inv_mu2
	complex*16 a(nferm, nferm)
	complex*16 s(nferm)
	complex*16 d(nferm), dhd(nferm), dd(nferm)
C       ----------------------------------------------------------------

	write(*, *) 'mu = ', mu

	inv_mu2 = 0.5d0 / mu**2

c       >>> fill in 'dense' corrected hypercube fermions
	do i = 1, nferm
	  s(i) = dcmplx(0.0d0, 0.0d0)
	end do
	do k = 1, nferm
	  s(k) = dcmplx(1.0d0, 0.0d0)

c	  >>> compute d = D s and dh = D^+ s
	  call mvec_dha(mu, s, d)

c	  >>> compute dhd = D^+ D
	  call mhvec_dha(mu, d, dhd)

c	  >>> compute d = (3 - D^+ D / \mu^2) s / 2
	  do i = 1, nferm
	    d(i) = 1.5d0 * s(i) - inv_mu2 * dhd(i)
	  end do

c	  >>> compute dd = D d
	  call mvec_dha(mu, d, dd)

c	  >>> finally, a = \mu s - dd
	  do i = 1, nferm
	    a(i, k) = mu * s(i) - dd(i)
	  end do

	  s(k) = dcmplx(0.0d0, 0.0d0)
	end do

	return
	end


C     ******************************************************************
      subroutine fill_dhc2(mu, a)
C     ******************************************************************
C     dhc: corrected hypercube fermions: expression for correction
C          taken from draft paper version 03 Feb 99
C     ------------------------------------------------------------------
C     hip, Jan 99; last modified: 04 Feb 99
C     ------------------------------------------------------------------
C     IN real*8 mu - shift of the origin; radius of the projected circle
C     OUT complex*16 a(nferm, nferm) - Dirac operator matrix
C     ------------------------------------------------------------------
	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)
C       ----------------------------------------------------------------
	real*8 mu
	real*8 mu102, mu34
	complex*16 a(nferm, nferm)
	complex*16 s(nferm)
	complex*16 d(nferm), dhd(nferm), y(nferm)
C       ----------------------------------------------------------------

	mu102 = 10.0d0 / mu**2
	mu34 = 3.0d0 / mu**4

c       >>> fill in 'dense' corrected hypercube fermions
	do i = 1, nferm
	  s(i) = dcmplx(0.0d0, 0.0d0)
	end do
	do k = 1, nferm
	  s(k) = dcmplx(1.0d0, 0.0d0)

c	  >>> compute d = D s
	  call mvec_dha(mu, s, d)

c	  >>> compute dhd = D^+ D
	  call mhvec_dha(mu, d, dhd)

c	  >>> compute y = (10 / mu^2 - 3 dhd / \mu^4)
	  do i = 1, nferm
	    y(i) = mu102 * s(i) - mu34 * dhd(i)
	  end do

c	  >>> compute dhd = D^+ D y
	  call mvec_dha(mu, y, d)
	  call mhvec_dha(mu, d, dhd)

c	  >>> compute y = (15 - dhd) / 8
	  do i = 1, nferm
	    y(i) = 1.875d0 * s(i) - 0.125d0 * dhd(i)
	  end do
	  
c	  >>> compute d = A y
	  call mvec_dha(mu, y, d)

c	  >>> finally, a = \mu s - d
	  do i = 1, nferm
	    a(i, k) = mu * s(i) - d(i)
	  end do

	  s(k) = dcmplx(0.0d0, 0.0d0)
	end do

	return
	end


C     ******************************************************************
      subroutine fill_dhc_old(a)
C     ******************************************************************
C     dhc: corrected hypercube fermions: expression for correction
C          taken from 13 Dec 98 mail by W. Bietenholz
C     ------------------------------------------------------------------
C     OUT complex*16 a(nferm, nferm) - Dirac operator matrix
C     ------------------------------------------------------------------
	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)
C       ----------------------------------------------------------------
	complex*16 a(nferm, nferm)
	complex*16 s(nferm)
	complex*16 d(nferm), dh(nferm), dd(nferm)
        complex*16 dhd(nferm), ddh(nferm), ddhd(nferm)
C       ----------------------------------------------------------------

c       >>> fill in 'dense' corrected hypercube fermions
c       (this can be done faster...)
	do i = 1, nferm
	  s(i) = dcmplx(0.0d0, 0.0d0)
	end do
	do k = 1, nferm
	  s(k) = dcmplx(1.0d0, 0.0d0)

c	  >>> compute d = D a and dh = D^+ a
	  call mvec_dhf(s, d)
	  call mhvec_dhf(s, dh)

c	  >>> compute dhd = D^+ D a and ddh = D D^+ a
	  call mhvec_dhf(d, dhd)
	  call mvec_dhf(dh, ddh)

c	  >>> compute dd = D^2 a and ddhd = D D^+ D a
	  call mvec_dhf(d, dd)
	  call mvec_dhf(dhd, ddhd)

c	  >>> compute (3 * D + D^+ - D^+ D - D^2 - D D^+ + D D^+ D) / 2
	  do i = 1, nferm
	    a(i, k) = 0.5d0 * 
     &        (d(i) - dh(i) + dhd(i) + dd(i) + ddh(i) - ddhd(i))
	  end do

	  s(k) = dcmplx(0.0d0, 0.0d0)
	end do

	return
	end



