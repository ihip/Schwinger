
	subroutine overlap(a)
c       compute overlap projection using epsilon(h)
c	INOUT complex*16 a(nferm, nferm) -
c         on INPUT: matrix which should be 'projected on the circle'
c           WARNING: it should obey g5-hermiticity property!!!
c         on OUTPUT: resulting matrix
c       hip, 27 Jul 98, last modified: 02 Feb 99

c     >>> necessary initializations for 2D U(1)
	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)

c       >>> size of work array for NAG eigenvalue subroutine
	parameter (lwork=64*nferm)

c       >>> 'dense' output matrix
	complex*16 a(nferm, nferm)

c       >>> auxiliary 'dense' matrix
	complex*16 b(nferm, nferm)

c       >>> vector with eigenvalues
	real*8 eig(nferm)

c       >>> necessary for LAPACK call
	real*8 rwork(3 * nferm)
	complex*16 work(lwork)

	complex*16 s(nferm), sum

c	>>> multiply with g5 to get hermitian matrix
	do j = 2, nferm, 2
	  do i = 1, nferm
	    a(j, i) = -a(j, i)	
	  end do
	end do
	
c       (diagonal entries should be exact real)
        do i = 1, nferm
          if(dimag(a(i, i)) .ne. 0.0d0) then
            write(*, *) 'im = ', dimag(a(i, i))
            a(i, i) = dcmplx(dreal(a(i, i)), 0.0d0)
	  end if
        end do

c	>>> diagonalization of herm. matrix (NAG: f02haf)
c	>>> WARNING: work and rwork interchanged compared to NAG
	ifail = 0
c	call cbl_cputime(0, tcpu1)
	call zheev('V', 'U', nferm, a, nferm, eig,
     &    work, lwork, rwork, ifail)
c	call cbl_cputime(1, tcpu1)
c       write(*, *) 'ifail = ', ifail
c       write(*, *) 'f02haf: ', tcpu1

c       >>> take signum of eigenvalues
	do i = 1, nferm
	  eig(i) = dsign(1.0d0, eig(i))
	end do

c       >>> reconstruct e(D) using eigenvectors
	do j = 1, nferm
	  do i = 1, nferm
	    s(i) = a(j, i) * eig(i)
	  end do
	  do m = 1, nferm
	    sum = dcmplx(0.0d0, 0.0d0)
	    do i = 1, nferm
		sum = sum + s(i) * dconjg(a(m, i))
	    end do
	    b(j, m) = sum
	  end do
	end do

c	>>> multiply by g5 to get A / |g5 A|
	do j = 1, nferm
	  do i = 1, nferm, 2
	    a(i, j) = b(i, j)
	    a(i + 1, j) = -b(i + 1, j)
	  end do
	end do

c       >>> add unity
	do m = 1, nferm
	  a(m, m) = a(m, m) + dcmplx(1.0d0, 0.0d0)
	end do

	end


