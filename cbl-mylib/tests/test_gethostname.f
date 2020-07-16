      program test_gethostname
C ihip / 2020-07-16 /
      character*32 hostname

      ierr = gethostname(hostname, 32)

      print *, 'ierr:     ', ierr
      print *, 'hostname: ', hostname
      
      stop
      end
