      function cbl_random(dummy)
C     ------------------------------------------------------------------
C     ******************************************************************
C     **                    **                                        **
C     **    R A N D O M     **    PROGRAM BY C.B. LANG,FEB 1989       **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Function produces real random numbers 0<=r<1 with the linear
C     congruential method. The multiplier is due to
C     D. Knuth, The Art of Computer Programming, Vol.2
C     ch. 3.3.4 (Lavaux and Janssens) and passes the spectral test
C     nicely.
C     ------------------------------------------------------------------
C     PROBLEMS
C     This type of random number generator has notoriuosly small periods.
C     For the following seeds one finds:   
C          seed         0        1       2       3     123     124
C        period   4194304  2097152 4194304  524288  131072 4194304
C                 = 2**24    2**23   2**24   2**21   2**17   2**24
C     
C     Thus even seeds are recommended. For better random numbers use
C     vecran.
C     ------------------------------------------------------------------
C     Calling time is 8.5 microsecs on a APOLLO DN4000 with FPX
C     Calling time is 1.8 microsecs on a APOLLO DN10000
C     ------------------------------------------------------------------
      integer*4 iseed, iout 
      common/rancom/iseed
      iseed=iseed*1664525
C
C                         mod 2**24
C
C      random= float(and(iseed,2**24-1))/float(2**24)
C     i=and(iseed,2**24-1)
c
      i=and(iseed,16777215)
      a=float(i)
      cbl_random=a*0.59604644775390625e-07
      end 
      subroutine cbl_ranget(iout)
C     ------------------------------------------------------------------
C     Delivers the current seed-value
C     ------------------------------------------------------------------
      integer*4 iout,iseed
      common/rancom/iseed
      data iseed/314159269/
      iout=iseed
      end
      subroutine cbl_ranset(iout)
C     ------------------------------------------------------------------
C     Overrides the current seed value
C     ------------------------------------------------------------------
      integer*4 iout,iseed
      common/rancom/iseed
      iseed=iout+314159269
      end

C     ------------------------------------------------------------------
      subroutine cbl_rantim(iseed)
C     ------------------------------------------------------------------
C     Initialize seed with time
C     ------------------------------------------------------------------
C     Get the 48-bit long time of day from the system and turn it
C     into a 32-bit long integer.
C     ------------------------------------------------------------------
      integer clock(3)
      call itime(clock)
      i = clock(1)
      j = clock(2)*128
      k = clock(3)*128*128
      iseed = i+j+k
      if((iseed.lt.0).or.(iseed.gt.900000000))then
        print *,'Error in rtim, iseed set to 54217137'
        iseed=54217137
      endif
      end


