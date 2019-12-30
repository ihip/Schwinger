      character*8 date,time
      call cbl_cputime(0,c) 
      print *,c
      call cbl_datim(date,time )
      print *,date
      print *,time
      sum=0
      do i=1,1000000
         sum=sum+sin(float(i))
      enddo
      print *,"dummy:",sum
      call cbl_cputime(1,c) 
      print *,c
      end
