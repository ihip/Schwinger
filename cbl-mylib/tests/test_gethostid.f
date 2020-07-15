      program test_gethostid
C ihip / 2020-01-03 /		
      integer id
      
      write(*, *) 'Test gethostid'
      call gethostid(id)
      write(*, *) id

      end

