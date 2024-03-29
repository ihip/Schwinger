C     ******************************************************************
      integer*4 function itopch()
C     ******************************************************************
C     **                    **                                        **
C     **  ITOPCH            **  Function by I. Hip, Jan 96            **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for determination of topological charge
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
     &        getphi(u(i2, 1)) - getphi(u(i, 2))) / pi
        chi = dmod(chi + 5.0d0, 2.0d0) - 1.0d0
        tpc = tpc + chi
      end do
      itopch = idnint(- tpc / 2.0d0)
      end


C     ******************************************************************
      integer*4 function itopch2(itp, itm)
C     ******************************************************************
C     **                    **                                        **
C     **  ITOPCH2           **  Function by I. Hip, 22 May 97         **
C     **                    **                                        **
C     ******************************************************************
C     Function returns topological charge of gauge field configuration u
C     ------------------------------------------------------------------
C     OUT (integer*4) itp  - positive top. charge
C     OUT (integer*4) itm  - negative top. charge
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

      itp = 0
      itm = 0

      tpc = 0.0d0
      do i = 1, nsite
        i1 = ind(1, i)
        i2 = ind(2, i)
        chi = (getphi(u(i, 1)) + getphi(u(i1, 2)) -
     &        getphi(u(i2, 1)) - getphi(u(i, 2))) / pi
	if((chi .gt. 1) .and. (chi .lt. 3)) itp = itp + 1 
	if(chi .gt. 3) itp = itp + 2
	if((-chi .gt. 1) .and. (-chi .lt. 3)) itm = itm - 1
	if(-chi .gt. 3) itm = itm - 2
        chi = dmod(chi + 5.0d0, 2.0d0) - 1.0d0
        tpc = tpc + chi
      end do
	write(*, *) itp, itm, itp + itm, idnint(- tpc / 2.0d0) 
      itopch2 = idnint(- tpc / 2.0d0)
      return
      end


C     ******************************************************************
	subroutine mk_index()
C     ******************************************************************
C     **                  **                                          **
C     **  MK_INDEX        **    Program by C.B. Lang, 02-03-1993      **
C     **                  **    Debugged by I. Hip, 26 Nov 96         **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for construction of d=2 index vectors
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Output variables:   ind(1:4,1:nsite)
C     ------------------------------------------------------------------
C     Remarks: index (1:nsite) = (1:NSPACE,1:NTIME)
C                                   i = (it-1)*NSPACE+ix
C    ------------------------------------------------------------------
C     ******************************************************************
	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
	integer ind(4,nsite)
C     ------------------------------------------------------------------
	common /index_arrays/ind
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

