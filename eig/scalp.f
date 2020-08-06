      complex*16 function scalp(n, v1, v2)
c     ivh, 09 Nov 97, v2.0
c     >>> computes scalar product of complex vectors v1 and v2
      complex*16 v1(*), v2(*), sum

      sum = 0.0d0
      do i = 1, n
        sum = sum + dconjg(v1(i)) * v2(i)
      end do

      scalp = sum
      return
      end


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

