C v3: 07 Jan 97, ivh - Last modified: 11 Aug 98, hip

C !!! vidit koje sve vrijednosti vratiti kao parametre !!!
      subroutine load_header(nf, nspace, ntime, nmeas, beta, akap)
C     ------------------------------------------------------------------
      character*8 c_ident
      real*8 beta, akap, eps
      character*64 c_line
      character*32 c_integrate, c_invert, c_invert_b, c_invert1,
     &  c_invert1_b, c_future, hostname
      character*8  c_date, c_time
      integer*4 ifut(8)
C     ------------------------------------------------------------------

      read(nf) c_ident, c_line,
     &  c_integrate, c_invert, c_invert_b, c_invert1,
     &  c_invert1_b, c_future, c_date, c_time, hostname,
     &  nspace, ntime, beta, akap, eps, nsteps, ntherm, nmeas,
     &  mstep, istart, iseed, ifut

C !!! ovo josh treba doraditi !!!
      call header('    S c h w i n g e r   ', c_date, c_time)
      write(*, 1001) nspace, ntime, c_integrate, c_invert,
     &  c_invert_b,
     &  beta, akap, eps, nsteps, ntherm, nmeas, mstep, istart, iseed
1001  format(' Schwinger Model on a ',i3,' x',i3,' Lattice',/,
     &       1x, 13('******'),/,
     &       ' Integration: ', a32, /,
     &       ' Inv:  ', a32, ' (', a32, ')', /,
     &       1x, 13('******'),/,
     &       ' beta                                = ',f10.5,/,
     &       ' akap                                = ',f10.5,/,
     &       ' eps                                 = ',f10.5,/,
     &       ' Length of Traj.                     = ',i10,/,
     &       ' Thermalization                      = ',i10,/,
     &       ' Measurements                        = ',i10,/,
     &       ' Measurement step                    = ',i10,/,
     &       ' Start (0 = COLD, 1 = HOT, 2 = FILE) = ',i10,/,
     &       ' Seed                                = ',i10,/,
     &       1x, 13('******'))
      write(*, *) 'Inv1: ', c_invert1, ' (', c_invert1_b, ')'

      write(*, *) 'Host: ', hostname

      return
      end


      subroutine load_dheader(nf, ncomp, c_computable)
      character*64 c_computable(MAX_NCOMP)

      read(nf) ncomp
      read(nf) (c_computable(i), i = 1, ncomp)

      return
      end


      subroutine load_mheader(nf, mcomp)
      character*64 c_source, c_direction, c_mfuture
      character*64 c_mcomputable(MAX_MCOMP)

      read(nf) mcomp, c_source, c_direction, c_mfuture
      read(nf) (c_mcomputable(i), i = 1, mcomp)

      write(*, *) 'Source:       ', c_source
      write(*, *) 'Direction:    ', c_direction
      write(*, '(1x, 13(''******''))')

      return
      end


      subroutine load_tail(nf)
      integer*4 nexcp, nexcp_b, nexcp1, nexcp1_b, itfut(4)
      real*4 tcpu

      read(nf) nexcp, nexcp_b, nexcp1, nexcp1_b, itfut, tcpu
      write(*, *)
      write(*, *) 'nexcp  = ', nexcp, ' (', nexcp_b, ')'
      write(*, *) 'nexcp1 = ', nexcp1, ' (', nexcp1_b, ')'
      write(*, *) 'Used CPU time: ', tcpu, ' sec'
      write(*, *)

      return
      end


C     ******************************************************************
      subroutine header(string, date, time)
C     ******************************************************************
C     **                   **                                         **
C     **  HEADER           **   Modified cbl_header, I. Hip, Mar 96   **
C     **                   **                                         **
C     ******************************************************************        
C     String length should be <= 42, else it is truncated
C     ------------------------------------------------------------------
C     ******************************************************************
      character*(*) string
      character*8 date,time
      character*78 line1,line2,line3

      do i=1,78
	line1(i:i)='*'
	line2(i:i)=' '
      enddo

      line2(1:2)='**'
      line2(77:78)='**'

      line3=line2

      ile=len(string)
      if(ile.gt.42)ile=42
      ite1=(45-ile)/2
      line2(2+ite1:1+ite1+ile)=string(1:ile)
      line2(53:60)=date
      line2(62:62)='/'
      line2(64:71)=time
      print *,line1
      print *,line3
      print *,line2
      print *,line3
      print *,line1
      end



