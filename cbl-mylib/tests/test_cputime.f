      program test_cputime
C ihip / 2020-01-03 /		
      real cputime
      real*8 x, sum
      
      write(*, *) 'Test cputime'
      write(*, *)
      write(*, *) 'Number of iterations: '
      read(*, *) niter
      write(*, *)
      write(*, *) 'Computing...'

      call cbl_cputime(0, cputime)
      sum = 0.0d0
      do i = 1, niter
        sum = sum + 1.0d0 / dble(i)**2
      end do
      call cbl_cputime(1, cputime)

      write(*, *)
      write(*, *) dsqrt(6.0d0 * sum)
      write(*, *)
      write(*, *) 'CPU time: ', cputime
      end
