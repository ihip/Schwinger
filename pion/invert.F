	subroutine invert(a)

	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)
	parameter (lwork=64*nferm)

	complex*16 a(nferm, nferm)
	complex*16 work(lwork)
	integer*4 ipiv(nferm)

c       >>> compute the LU factorization (NAG: f07arf)
	info = 0
c	call cbl_cputime(0, tcpu1)
	call zgetrf(nferm, nferm, a, nferm, ipiv, info)
c	call cbl_cputime(1, tcpu1)
c       write(*, *) 'info: ', info
c       write(*, *) 'f07arf (LU factorization): ', tcpu1

c       >>> compute inverse (out of LU) (NAG: f07awf)
	info = 0
c	call cbl_cputime(0, tcpu1)
	call zgetri(nferm, a, nferm, ipiv, work, lwork, info)
c	call cbl_cputime(1, tcpu1)
c       write(*, *) 'info: ', info
c       write(*, *) 'f07awf (matrix inversion): ', tcpu1
c       write(*, *)

	end

