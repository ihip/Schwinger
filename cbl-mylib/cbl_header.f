C     ******************************************************************
      subroutine cbl_header(string)
C     ******************************************************************
C     **                   **                                         **
C     **  CBL_HEADER       **     Implemented by C.B. Lang, Dec. 1993 **
C     **                   **                                         **
C     ******************************************************************	
C     ------------------------------------------------------------------
C     Prints a standardized header for the calling program.
C     Usage: call cbl_header('My Program Title')  
C     Stringlength should be <= 42, else it is truncated
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
      call cbl_datim(date,time)
      line2(53:60)=date
      line2(62:62)='/'
      line2(64:71)=time  
      print *,line1
      print *,line3
      print *,line2
      print *,line3
      print *,line1
      end
