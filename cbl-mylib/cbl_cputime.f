C     ******************************************************************	
      subroutine cbl_cputime(iflag,cputime) 
C     ******************************************************************
C     **                   **                                         **
C     **  CBL_CPUTIME      **     Implemented by C.B. Lang, Dec. 1993 **
C     **                   **                                         **
C     ******************************************************************	
C     ------------------------------------------------------------------
c     user CPU time in seconds 
c     iflag = 0 gives elapsed cputime since beginning of process 
c     iflag = 1 gives elapsed cputime since last call to cputime
C     ------------------------------------------------------------------
C     ******************************************************************	
      real tarray(2)
      save sect0
      data sect0/0.0/
c
      total=etime(tarray)
c      
      if ( iflag .eq. 0 ) then
        cputime = tarray(1)
      else  
        cputime = tarray(1) - sect0
      endif
      sect0 = tarray(1)   
c
      end
