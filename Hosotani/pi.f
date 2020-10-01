      program pi

      real spi2
      real*8 dpi2

      print *, 'Two pi'
      print *

      spi2 = atan(1.0) * 8
      print *, 'real (single precision)'
      print *, spi2
      print *

      dpi2 = datan(1.0d0) * 8
      print *, 'real*8 (double precision)'
      print *, dpi2
      print *

      print *, '20 digits (Sagemath)'
      print *, '  6.2831853071795864769'
      print *

      end
     