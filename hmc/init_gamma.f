      subroutine init_gamma()
C       Ivan Hip, 30. Mar '96 (based on code by HCH & RET)

      parameter (ndrc=2,ngrp=1,ndim=2)

      complex*16 gamma1(ndrc, ndrc), gamma2(ndrc, ndrc), 
     &  gamma5(ndrc, ndrc), gamma1x5(ndrc, ndrc), gamma2x5(ndrc, ndrc) 

      common /gamma/ gamma1, gamma2, gamma5, gamma1x5, gamma2x5

      gamma1(1, 1) = dcmplx( 0.0, 0.0)
      gamma1(1, 2) = dcmplx( 1.0, 0.0)
      gamma1(2, 1) = dcmplx( 1.0, 0.0)
      gamma1(2, 2) = dcmplx( 0.0, 0.0)      

      gamma2(1, 1) = dcmplx( 0.0, 0.0)
      gamma2(1, 2) = dcmplx( 0.0,-1.0)
      gamma2(2, 1) = dcmplx( 0.0, 1.0)
      gamma2(2, 2) = dcmplx( 0.0, 0.0)      

      gamma5(1, 1) = dcmplx( 1.0, 0.0)
      gamma5(1, 2) = dcmplx( 0.0, 0.0)
      gamma5(2, 1) = dcmplx( 0.0, 0.0)
      gamma5(2, 2) = dcmplx(-1.0, 0.0)      

      gamma1x5(1, 1) = dcmplx( 0.0, 0.0)
      gamma1x5(1, 2) = dcmplx(-1.0, 0.0)
      gamma1x5(2, 1) = dcmplx( 1.0, 0.0)
      gamma1x5(2, 2) = dcmplx( 0.0, 0.0)
      
      gamma2x5(1, 1) = dcmplx(0.0, 0.0)
      gamma2x5(1, 2) = dcmplx(0.0, 1.0)
      gamma2x5(2, 1) = dcmplx(0.0, 1.0)
      gamma2x5(2, 2) = dcmplx(0.0, 0.0)     

      return
      end


