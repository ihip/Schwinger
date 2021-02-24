C v3 - ivh, 24 Jan 97 - Last modified: 21 Feb 97

#define SAVE_EXCP .true.

      integer*4 function invert_excp(x, b)

      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)

      real*8 x(nall), b(nall)
      real*8 xx(nall)

      common/excpcomm/ nexcp, nexcp_b, nexcp1, nexcp1_b

C     >>> copy initial guess
      do i = 1, nall
        xx(i) = x(i)
      end do

C     >>> call first invert
      niter = invert(xx, b)

c     >>> if no convergence, call invert_b
      if(niter .lt. 0) then
        nexcp = nexcp + 1
        do i = 1, nall
          xx(i) = x(i)
        end do
        niter = invert_b(xx, b)
        if(niter .lt. 0) then
          nexcp_b = nexcp_b + 1
          if(niter .eq. -1) niter = -3
          if(niter .eq. -2) niter = -4
          if(SAVE_EXCP) call save_gfield(nexcp_b, 0)
          invert_excp = niter
          return
        end if
      end if

C     >>> return the result (and use it as initial guess)
      do i = 1, nall
        x(i) = xx(i)
      end do

      invert_excp = niter
      return
      end


      subroutine save_gfield(nexcp_b, iflag1)
C     >>> if iflag1 = 1 then .exc1 file is created
      
      parameter (ndrc = 2,ngrp = 1,ndim = 2)
      parameter (nsite = NTIME * NSPACE)

      character*32 fname
      character*6 c_beta, c_kappa
      character*2 c_nexcp_b
      complex*16 u(nsite, ndim)
      real*8 akap, beta, dummyeps

      common/kappa/akap
      common/constants/beta
      common/gauge_fields/ u

C     >>> if there are to many excp-configurations, there is no more
C         sense to save all of them (first 99 would be enough!)
      if(nexcp_b .gt. 99) return

C     >>> prepare the file name
      ival = int((beta + 100.0000001d0) * 1000.0d0)
      write(c_beta, '(i6)') ival
      ival = int((akap + 1.000000001d0) * 100000.0d0)
      write(c_kappa, '(i6)') ival
      write(c_nexcp_b, '(i2)') nexcp_b
      if(iflag1 .eq. 1) then
        fname = 'b'//c_beta(2:6)//'k'//c_kappa(2:6)//c_nexcp_b//'.exc1'
      else
        fname = 'b'//c_beta(2:6)//'k'//c_kappa(2:6)//c_nexcp_b//'.excp'
      end if
      
C     >>> save the gauge configuration
      open(13, file = fname, form = 'unformatted', status = 'new')
      write(13) NTIME, NSPACE
C     >>> save dummyeps for compatibility reasons
      dummyeps = 0.0d0
      write(13) beta, akap, dummyeps
      write(13) it
      write(13) u
      close(13)

      return
      end

       
