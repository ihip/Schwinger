      real*8 function op_pbp_e(g)
C     ------------------------------------------------------------------
C     computes chiral condensate (ivh, 16 May 96); mod. 03 Aug 98
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm
	sum = sum + g(i, i) 
      end do
      op_pbp_e = dreal(sum)

      return
      end

      real*8 function op_pbg5p_e(g)
C     ------------------------------------------------------------------
C     computes \bar{\Psi} \gamma_5 \Psi (ivh, 18 May 96)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
	sum = sum + (g(i, i) - g(i + 1, i + 1))
      end do
      op_pbg5p_e = dreal(sum)

      return
      end

      real*8 function op_ubudbd_e(g)
C     ------------------------------------------------------------------
C     computes $ \bar{u} u \bar{d} d $ (ivh, 29 May 96)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
	sum = sum + g(i, i)**2 + 2.0d0 * g(i, i) * g(i + 1, i + 1) +
     &    g(i + 1, i + 1)**2
      end do
      op_ubudbd_e = dreal(sum)

      return
      end


      real*8 function op_ubg5udbg5d_e(g)
C     ------------------------------------------------------------------
C     computes $ \bar{u} u \bar{d} d $ (ivh, 16 Jun 97)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
	sum = sum + g(i, i)**2 - 2.0d0 * g(i, i) * g(i + 1, i + 1) +
     &    g(i + 1, i + 1)**2
      end do
      op_ubg5udbg5d_e = dreal(sum)

      return
      end
