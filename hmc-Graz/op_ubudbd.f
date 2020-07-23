      complex*16 function op_ubudbd(xi, chi)
C code by ret, modified by ivh, 30. Mar '96
C computes ubudbd using noisy estimator

      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)

      complex*16 xi(nferm), chi(nferm)
      complex*16 xihat(nferm), phihat(nferm), chihat(nferm)
      complex*16 fourscalp
C     >>> dimension of rando is max(2*nferm, ngaug)
      real rando(2*nferm)

C     >>> we need some Gaussian "noise" xihat
C     >>> (xihat should be complex and < (xihat, xihat) > = 1)
      call cgauss(xihat, nferm, rando)

C     >>>  phihat = (M^+) xihat
      call mtvec(xihat, phihat)

      niter = invert(chihat, phihat)

      op_ubudbd = fourscalp(xi, chi, xihat, chihat, nsite, ndrc)

      end


      complex*16 function fourscalp(a, b, c, d, nsite, ndrc)
C code by ret, 23 May 96
C sum over all sites (x = 1, nsite; 11, 22 are Dirac indices):
C   (M^{-1}_{xx11})^{2} + 2*M^{-1}_{xx11}M^{-1}_{xx22} +
C   (M^{-1}_{xx22)^{2}

      integer nsite, ndrc
      complex*16 a(*), b(*), c(*), d(*), sum

      sum = dcmplx(0.0, 0.0)
      do i = 0, nsite - 1
        do j = 1, ndrc
         do k=1,ndrc

          sum = sum + dconjg(a(i * ndrc + j)) *
     & dconjg(c(i * ndrc + k)) * b(i * ndrc + j) * d(i * ndrc + k)

           enddo
        end do
      end do

      fourscalp = sum
      end

