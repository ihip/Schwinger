	subroutine lapackeig(q, eig, rtol, nreal, reig,
     &    vr, maxevect, detr, deti, nuf)
c     INOUT complex*16 q(nferm, nferm) - on entry: q matrix
c       WARNING: on exit, q is overwritten and destroyed!!!
c     OUT complex*16 eig(nferm) - eigenvalues of q
c     IN real*8 rtol - tolerance for identification of real EVs
c     OUT integer*4 nreal - number of identified real EVs
c     OUT complex*16 reig(maxevect) - identified real EVs
c     OUT vr(nferm, maxevect) - right eigenvectors of q which
c       belong to the real EVs in reig
c     IN maxevect - maximal number of evectors which can be stored
c     OUT real*8 detr, deti - real and imaginary part of the 
c       determinant computed using EVs
c     OUT integer*4 nuf - fermionic nu, defined via zero modes
c     >>> hip, 05 Aug 98; last modified: 11 Aug 98

	parameter (ndrc=2,ngrp=1,ndim=2)
	parameter (nsite=NTIME*NSPACE)
	parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
	parameter (nall=2*ndrc*ngrp*nsite)

      complex*16 q(nferm, nferm)
      complex*16 eig(nferm), reig(maxevect), vr(nferm, maxevect)

	real*8 rtol, detr, deti
	complex*16 det

c     >>> needed for LAPACK subroutines
      complex*16 hq(nferm, nferm)
      complex*16 work(nferm), tau(nferm), vl
      real*8 rwork(2*nferm)
      logical select(nferm)
      integer ifaill, ifailr(nferm)

c     >>> reduction to Hessenberg form
      info = 0
c      call cbl_cputime(0, tcpu)
      call zgehrd(nferm, 1, nferm, q, nferm, tau, hq, nferm**2, info)
c      call cbl_cputime(1, tcpu)
c      write(1, *) 'zgehrd: info = ', info
c      write(1, *) 'zgehrd: hq(1, 1) = ', hq(1, 1)
c      write(1, *) 'zgehrd: tcpu = ', tcpu
c      write(1, *)

c     >>> save the contents of q (Hessenberg form + details of
c         unitary transformation) for future use
      do j = 1, nferm
	do i = 1, nferm
	  hq(i, j) = q(i, j)
	end do
      end do

c     >>> computation of (all) eigenvalues (contents of q is destroyed)
      info = 0
c      call cbl_cputime(0, tcpu)
      call zhseqr('E', 'N', nferm, 1, nferm, q, nferm, eig,
     &  vl, 1, work, nferm, info)
c      call cbl_cputime(1, tcpu)
c      write(1, *) 'zhseqr: info = ', info
c      write(1, *) 'zhseqr: tcpu = ', tcpu
c      write(1, *)

c     >>> select the real eigenvalues
      irealeig = 0
      inuf = 0
      det = dcmplx(1.0d0, 0.0d0)
      do i = 1, nferm
	det = det * eig(i)
	if(dabs(dimag(eig(i))) .lt. rtol) then
	  irealeig = irealeig + 1
	  if(dreal(eig(i)) .lt. 0.5d0) inuf = inuf + 1
	  if(irealeig .le. maxevect) then
	    reig(irealeig) = eig(i)
	    select(i) = .true.
	  else
	    select(i) = .false.
	  end if
	else
	  select(i) = .false.
	end if
      end do

c     >>> computation of the selected right evects of Hessenberg matrix
      msel = 0
      info = 0
c      call cbl_cputime(0, tcpu)
      call zhsein('R', 'Q', 'N', select, nferm, hq, nferm, eig, vl, 1,
     &  vr, nferm, maxevect, msel, q, rwork, ifaill, ifailr, info)
c      call cbl_cputime(1, tcpu)
c      write(1, *) 'zhsein: msel = ', msel
c      write(1, *) 'zhsein: info = ', info
c      write(1, *) 'zhsein: tcpu = ', tcpu
c      write(1, *)

c     >>> transformation of evecs to that of the original matrix
      info = 0
c      call cbl_cputime(0, tcpu)
      call zunmhr('L', 'N', nferm, msel, 1, nferm, hq, nferm, tau,
     &  vr, nferm, q, nferm**2, info)
c      call cbl_cputime(1, tcpu)
c      write(1, *) 'zunmhr: info = ', info
c      write(1, *) 'zunmhr: tcpu = ', tcpu
c      write(1, *)

c     >>> normalization of eigenvectors
      do i = 1, msel
	call normalize(nferm, vr(1, i))
      end do

      detr = dreal(det)
      deti = dimag(det)
      nreal = irealeig
      nuf = inuf
      
      return
      end

