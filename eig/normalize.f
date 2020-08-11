      subroutine normalize(n, v)
c     ivh, 09 Nov 97, v2.0
c     >>> normalizes the complex vector v
      complex*16 v(*)

      complex*16 scalp
      real*8 norm

      norm = dsqrt(dreal(scalp(n, v, v)))

      do i = 1, n
        v(i) = v(i) / norm
      end do

      return
      end

