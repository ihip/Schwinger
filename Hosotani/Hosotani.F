      program Hosotani
c   Hosotani equation (\delta alpha = 0, \kappa = 1)

      parameter(nmax = 10000)

c   size of LAPACK working array
      parameter(lwork = 34 * nmax)

c   matrix
      real*8 m(nmax, nmax)

c   eigenvalues
      real*8 ev(nmax)

c   LAPACK working arrays
      real*8 work(lwork)

c   grid spacing and grid spacing squared
      real*8 df, df2

      print *
      print *, '"Hosotani" equation'
      print *, '~~~~~~~~~~~~~~~~~~~'
      print *
      print *, 'nmax = ', nmax
      print *
      write(*, '(1x, ''ngrid (<= nmax) = '', $)')
      read(*, *) ngrid
      if (ngrid .gt. nmax) stop 'ngrid > nmax'

c   df = 2 * pi / n
      df = 8.0d0 * datan(1.0d0) / ngrid
      df2 = df**2
      print *
      print *, 'df = ', df
      print *

c   initialize matrix m with zeroes
      do i = 1, ngrid
        do j = 1, ngrid
          m(i, j) = 0.0d0
        end do
      end do

c   terms on the diagonal
      do i = 1, ngrid
        m(i, i) = 2 / df2 - dcos((i - 1) * df)
      end do

c   terms above and bellow diagonal
      do i = 1, ngrid - 1
        m(i + 1, i) = -1 / df2
        m(i, i + 1) = -1 / df2
      end do

c   periodic boundaries
      m(ngrid, 1) = -1 / df2
      m(1, ngrid) = -1 / df2

c   call LAPACK
      call DSYEV('N', 'U', ngrid, m, nmax, ev, work, lwork, info)

      print *, 'info = ', info
      print *, 'optimal lwork = ', work(1)   
      print *

c   write 10 smallest eigenvalues
      do i = 1, 10
        print *, i, ev(i)
      end do
      print *

      end
