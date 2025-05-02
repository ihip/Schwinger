      program readconf
C     reads binary .conf file and save its content as text file

      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
C     ------------------------------------------------------------------
      complex*16 u(nsite, ndim)
      character*32 input_file_name, output_file_name
      real*8 beta, akap, eps
      real*8 sgauge
C     ------------------------------------------------------------------
C     >>> file name preparation
      character*4 c_nspace, c_ntime
      character*10 c_size
      character*6 c_beta, c_kappa

C     reads name of .conf file from keyboard
      write(*, *)
      write(*, *) 'readconf / Hip / 2025-05-02'
      write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*, *)
      write(*, '(1x, ''.conf file name = '', $)')
      read(*, *) input_file_name

C     reads .conf file
      open(1, file = input_file_name, form = 'unformatted',
     & status = 'old')
	  read(1) ntime, nspace
      read(1) beta, akap, eps
	  read(1) it
	  read(1) u
	  close(1)

C     writes header information on screen
      write(*, *)
      write(*, *) 'ntime =  ', ntime
      if(ntime .ne. NTIME) stop "ntime <> NTIME"
      write(*, *) 'nspace = ', nspace
      if(ntime .ne. NSPACE) stop "nspace <> NSPACE"
      write(*, *) 'beta = ', beta            
      write(*, *) 'akap = ', akap
      write(*, *) 'eps = ', eps
      write(*, *) 'it = ', it
      write(*, *)

C     computes gauge action as a kind of checksum
      call mk_index()
      write(*, *) 'sgauge = ', sgauge(beta, u) / (ntime * nspace)
      write(*, *)

C     prepares output file name
      write(c_nspace, '(i4)') nspace + 1000
      write(c_ntime, '(i4)') ntime + 1000
      c_size = 'L'//c_nspace(2:4)//'x'//c_ntime(2:4)//'_b'
      ival = int((beta + 100.0000001d0) * 1000.0d0)
      write(c_beta, '(i6)') ival
      ival = int((akap + 1.000000001d0) * 100000.0d0)
      write(c_kappa, '(i6)') ival
      output_file_name =
     & c_size//c_beta(2:6)//'_k'//c_kappa(2:6)//'.ctxt'
      write(*, *) 'Gauge fields written to: ', output_file_name
      write(*, *)

C     writes gauge configuration as a text file
      open(2, file = output_file_name, form = 'formatted',
     & status = 'unknown')
      do it = 1, ntime
        do ix = 1, nspace
C         computes site index
          is = (it - 1) * nspace + ix
C         writes two links
          write(2, *) it, ix, dreal(u(is, 1)), dimag(u(is, 1))
          write(2, *) it, ix, dreal(u(is, 2)), dimag(u(is, 2))
        enddo
      enddo
      close(2)
      
      end


C >>> needed only to compute gauge action:

C     ******************************************************************
      subroutine mk_index()
C     ******************************************************************
C     **                  **    Program by C.B. Lang, 02-03-1993      **
C     **  MK_INDEX        **    Debugged by I. Hip, 26 Nov 96         **
C     **                  **    Last modified by ihip, 11 Aug 20      **
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
C     ------------------------------------------------------------------
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
