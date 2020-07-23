      program test_gettimeofday
C ihip / 2020-07-23 /		
      integer itm(2), izone(2)
      
      write(*, *) 'Test gettimeofday'
      call gettimeofday(itm, izone)
      write(*, *) 'itm(1) = ', itm(1)
      write(*, *) 'itm(2) = ', itm(2)

      end

