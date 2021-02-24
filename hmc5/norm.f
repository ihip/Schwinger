      real*8 function norm(n, v)
c     modified by ihip, 10 Aug 20
      integer n
      real*8 v(*), sum
      
      sum = 0.0d0

      do i = 1, n
	    sum = sum + v(i)**2
      end do

      norm = sum
      return
      end
      