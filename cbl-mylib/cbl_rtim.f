C     ******************************************************************	
      subroutine cbl_rtim(iseed)
C     ******************************************************************
C     **                   **                                         **
C     **    CBL_RTIM       **     Implemented by C.B. Lang, Jan. 1993 **
C     **                   **     Modified by I. Hip, Jul. 2020       **
C     ******************************************************************	
C     Initialize seed with time
C     ******************************************************************	
      integer itm(2),izone(2),id
      integer*8 lseed
c
c     add the hostid
c
      call gethostid(id)
C       print *, 'id =          ', id
      lseed=mod(iabs(id),900000000)   
C       print *, 'iseed 1 =     ', lseed
      call gettimeofday(itm,izone)
C       print *, 'itm(1) =      ', itm(1)
C       print *, 'itm(2) =      ', itm(2)
C       print *, 'iseed 2 =     ', lseed+itm(1)+itm(2)
C       print *, 'iseed 2 mod = ', mod(lseed+itm(1)+itm(2),900000000)
      iseed = mod(lseed+itm(1)+itm(2),900000000)
C       print *, 'iseed 3 =     ', iseed
      if((iseed.lt.0).or.(iseed.gt.900000000))then
        print *,'Error in cbl_rtim, iseed set to 54217137'
        iseed=54217137
      endif
      end
    