      subroutine cbl_open(iunit,afile,aform,astat)
      character*(*) afile,aform,astat
      open(unit=iunit,file=afile,form=aform,err=10,
     &recl=80,status=astat)
      goto 50
10    print *,'Error in my_open'
50    end
