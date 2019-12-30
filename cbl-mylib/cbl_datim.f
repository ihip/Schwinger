C     ******************************************************************
      subroutine cbl_datim(date,time)
C     ******************************************************************
C     **                    **                                        **
C     **  C B L _ D A T I M **    PROGRAM BY C.B. LANG,Dec. 1993      **
C     **                    **                                        **
C     ******************************************************************
c
c     subroutine gives date 'dd-mm-yy' and time 'hh:mm:ss'
c     as two character*8 strings
c
C     ******************************************************************
      integer dati(3),tim(3)
      character*8 date,time
      call idate(dati(2),dati(1),dati(3))
      call itime(tim) 
      write(date,'(i2,1h-,i2,1h-,i2)')dati(1),dati(2),dati(3)
      write(time,'(i2,1h:,i2,1h:,i2)')tim(1),tim(2),tim(3)
      end
