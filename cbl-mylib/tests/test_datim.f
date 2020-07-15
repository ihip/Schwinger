      program test_datim
C ihip / 2019-12-30 /		
      character*8 date, time
      
      write(*, *) 'Test datim'
      call cbl_datim(date, time)
      write(*, *) date
      write(*, *) time
      
      end
