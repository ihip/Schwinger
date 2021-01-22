      complex*16 function scalp(n, v1, v2)
c     ivh, 09 Nov 97, v2.0
c     >>> computes scalar product of complex vectors v1 and v2
      integer n
      complex*16 v1(*), v2(*), sum

      sum = dcmplx(0.0, 0.0)
      do i = 1, n
        sum = sum + dconjg(v1(i)) * v2(i)
      end do

      scalp = sum
      return
      end


