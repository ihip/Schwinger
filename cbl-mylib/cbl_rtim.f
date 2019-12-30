C     ******************************************************************	
      subroutine cbl_rtim(iseed)
C     ******************************************************************
C     **                   **                                         **
C     **    CBL_RTIM       **     Implemented by C.B. Lang, Jan. 1993 **
C     **                   **                                         **
C     ******************************************************************	
C     Initialize seed with time
C     ******************************************************************	
      integer itm(2),izone(2),gethostid
c
c     add the hostid
c
      do i=1,10
        iseed=mod(iabs(gethostid()),900000000)
      enddo    

      call gettimeofday(itm,izone)
      iseed = mod(iseed+itm(1)+itm(2),900000000)
      if((iseed.lt.0).or.(iseed.gt.900000000))then
        print *,'Error in cbl_rtim, iseed set to 54217137'
        iseed=54217137
      endif
      end

