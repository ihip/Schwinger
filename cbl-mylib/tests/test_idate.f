      program test_idate
C ihip / 2019-12-30 /		
      integer id(3)
      
      write(*, *) 'Test idate'
      call idate(id)
      write(*, *) id(1), id(2), id(3)

      end

