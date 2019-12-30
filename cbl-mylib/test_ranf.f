      real ran(21007)
      sum=0
      call cbl_cputime(0,t)
      call cbl_rcarry(ran,21007,0.)
c      do i=1,21007
c      sum=sum+ran(i)
c      enddo
c        sum=sum+cbl_ranf()
c      enddo
      call cbl_cputime(1,td)
      print *,sum,td/21007
      end
C     ******************************************************************
      function cbl_ranf()
C     ******************************************************************
C     **                    **                                        **
C     **      cbl_ranf      **   Program by C.B. Lang, DATE: 01/14/94 **  
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for producing just one random number                               
C     (this is used instead of the usual bad cbl_random)
C     ------------------------------------------------------------------
C     Output variables:
C     cbl_ranf
C     ------------------------------------------------------------------
C     Remarks:
C     functions gets block of random numbers from an efficient generator
C     and the returns them one by one. This should only be used, if
C     one really needs just single random numbers...
C
C     CPU per call on R4000 (IRIX5.1.2): 1.904 * 10**(-6) sec
C     ------------------------------------------------------------------
C     Compile with:
C     f77 cbl_ranf.f
C     ******************************************************************
      parameter (len=1000)
      real ran(len)
      save ran,ilen
      data ilen/0/
      if(ilen.eq.0)then
        call cbl_rcarry(ran,len,0.)
        ilen=1
      endif
      cbl_ranf=ran(ilen)
      ilen=mod(ilen+1,len)
      end
