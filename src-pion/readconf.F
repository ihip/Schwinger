C     ******************************************************************
      subroutine readconf(fname, beta, u)
C     ******************************************************************
C     **                    **                                        **
C     **  READCONF          **  Subroutine by I. Hip, 2022-03-24      **
C     **                    **  Last modified:                        **
C     ******************************************************************
C     Read gauge fields u
C     ------------------------------------------------------------------
C     Input  variables:  fname (file with saved configuration)
C     ------------------------------------------------------------------
C     Output variables:  beta, u (gauge field)
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 u(nsite, ndim)
      character*32 fname
      real*8 beta, akap, eps
C     ------------------------------------------------------------------    

C     >>> read a gauge configuration
      open(11, file = fname, form = 'unformatted', status = 'old')
	  read(11) ntime, nspace
      read(11) beta, akap, eps
	  read(11) it
	  read(11) u
	  close(11)
      
      if((ntime .ne. NTIME) .or. (nspace .ne. NSPACE)) 
     &  stop 'Inappropriate lattice size!'
      return
      end
