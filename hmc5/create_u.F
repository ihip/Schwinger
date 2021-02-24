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
        open(1, file = input_file, form = 'unformatted',
     & status = 'old')
	    read(1) ntime, nspace
        read(1) beta, akap, eps
	    read(1) it
	    read(1) u
	    close(1)
	    return
      endif
      
      stop 'Invalid istart parameter'
      end
