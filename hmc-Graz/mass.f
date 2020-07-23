C This subroutine measures propagation between time-slices only!!!

C     ******************************************************************
      integer*4 function mass(s0, s1, s2, s3)
C     ******************************************************************
C     **                    **                                        **
C     **  MASS              **  Ivan Hip, 13 Apr 96                   **
C     **                    **  Last modified: 18 Sep 96              **
C     ******************************************************************
C     Compute correlations between density currents
C     ------------------------------------------------------------------
C     Part of package:   HMC2DU1
C     ------------------------------------------------------------------
C     Input parameters:  
C     ------------------------------------------------------------------
C     Output parameters: 
C     ------------------------------------------------------------------
C     WARNING: it is assumed that ndrc = 2
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 q(nferm, ndrc), s(nferm)
      real*8 s0(NTIME - 1), s1(NTIME - 1), s2(NTIME - 1), s3(NTIME - 1)
      real*8 absdg, absnd, mixdg, mixnd
C     ------------------------------------------------------------------
      save q, qt
      data q/nall*1.0d0/

      do it = 1, NTIME - 1
        s0(it) = 0.0d0
        s1(it) = 0.0d0
        s2(it) = 0.0d0
        s3(it) = 0.0d0
      end do

      do ix = 1, NSPACE

	call msource(ix, 1, s)
        niter1 = invert1(q(1, 1), s)
        if(niter1 .eq. 0) then
          mass = 0
          return
        end if

	call msource(ix, 2, s)
        niter2 = invert1(q(1, 2), s)
        if(niter2 .eq. 0) then
          mass = 0
          return
        end if
 
	do it = 1, NTIME - 1
          do j = 1, NSPACE
            k = ndrc * ((it * NSPACE) + j) - 1
            absdg = dreal(
     & dconjg(q(k, 1)) * q(k, 1) + dconjg(q(k + 1, 2)) * q(k + 1, 2))
            absnd = dreal(
     & dconjg(q(k, 2)) * q(k, 2) + dconjg(q(k + 1, 1)) * q(k + 1, 1))
            mixdg = dreal(
     & dconjg(q(k + 1, 2)) * q(k, 1) + dconjg(q(k, 1)) * q(k + 1, 2))
            mixnd = dreal(
     & dconjg(q(k, 2)) * q(k + 1, 1) + dconjg(q(k + 1, 1)) * q(k, 2))
c	  >>> scalar density
	    s0(it) = s0(it) + (absdg - absnd)
c	  >>> j_1 current
	    s1(it) = s1(it) + (mixdg - mixnd)
c	  >>> j_2 current
	    s2(it) = s2(it) + (mixdg + mixnd)
c	  >>> pseudoscalar density
	    s3(it) = s3(it) + (absdg + absnd)
          end do
        end do
      end do

      mass = niter1 + niter2
      return
      end


      subroutine msource(ix, k, s)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 s(nferm)

      do i = 1, nferm
        s(i) = 0.0d0
      end do
      s((ix - 1) * ndrc + k) = 1.0d0

      return
      end

