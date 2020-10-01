      program hop
c   harmonic oscillator with periodic boundary

      parameter(ngrid = 100)

c   size of LAPACK working array
      parameter(lwork = 34 * ngrid)

c   matrix
      real*8 m(ngrid, ngrid)

c   eigenvalues
      real*8 ev(ngrid)

c   LAPACK working arrays
      real*8 work(lwork)

c   initialize matrix m with zeroes
      do i = 1, ngrid
        do j = 1, ngrid
          m(i, j) = 0.0d0
        end do
      end do

c   terms on the diagonal
      do i = 1, ngrid
        m(i, i) = -2
      end do

c   terms above and bellow diagonal
      do i = 1, ngrid - 1
        m(i + 1, i) = 1
        m(i, i + 1) = 1
      end do

c   periodic boundaries
      m(ngrid, 1) = 1
      m(1, ngrid) = 1

c   call LAPACK
      call DSYEV('N', 'U', ngrid, m, ngrid, ev, work, lwork, info)

      print *, 'info = ', info
      print *, 'optimal lwork = ', work(1)   
      print *

c   write eigenvalues
      do i = 1, ngrid
        print *, i, ev(i)
      end do

      end
