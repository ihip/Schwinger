C     ******************************************************************
      integer function cbl_len(st)
C     ******************************************************************
C     **                    **                                        **
C     **      cbl_len       **   Program by C.B. Lang, DATE: 01/18/94 **  
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for determination of actual length of a string (like
C     len(string) but not counting trailing blanks
C     ------------------------------------------------------------------
C     Output variables:
C     cbl_len
C     ------------------------------------------------------------------
C     Remarks:
c     don not forget to declare cbl_len intereger!
C     ------------------------------------------------------------------
C     Compile with:
C     f77 cbl_len.f
C     ******************************************************************
      character*(*) st
      le=len(st)
      do i=le,1,-1
        j=i
        if(st(j:j).ne.' ')goto 1
      enddo
      j=0
1     cbl_len=j
999   end
