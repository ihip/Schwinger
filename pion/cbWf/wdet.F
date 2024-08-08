c maximal lattice size 64 * 64 = 4096 (* 2 spinor components)
#define MAX_EIGEN 8192

	subroutine wdet(cfname, fmass, rtol, nuf, edetr, sigma)
c   IN character*64 cfname - gauge configuration file name (.conf)
c   IN real*8 fmass - fermion mass
c   IN real*8 rtol - tolerance for identification of real EVs
c   OUT integer*4 nuf - fermionic nu, defined via zero modes
c   OUT real*8 edetr - normalized determinant
c	OUT real*8 sigma - chiral condensate
c   >>> read appropriate EVs and computes normalized determinant
c       which is used for reweighting
c   >>> Hip, 2022-04-03; last modified: 2022-04-03

	real*8 re, im
	character*64 cfname, eigfile
	real*8 fmass, rtol, edetr, sigma
	complex*16 lambda(MAX_EIGEN), sum, pro

c   >>> get eigfile from cfname (cut .conf and append .eig)
	eigfile = cfname(1:lnblnk(cfname) - 5)//'-CB.eig'
	write(*, *) 'WGET: eigfile = ', eigfile

c	>>> read eigenvalues from eigfile
	open(12, file = eigfile, form = 'formatted', status = 'old')
	neigen = 0
	nuf = 0
7	read(12, *, end = 77) re, im
	  neigen = neigen + 1
	  lambda(neigen) = re + dcmplx(0, 1) * im
c	  >>> count real modes, i.e. fermionic top. charge
	  if((dabs(im) .lt. rtol) .and. (re .lt. 0.5d0)) nuf = nuf + 1
	goto 7
	
77	continue
	close(12)
    
	  nsites = neigen / 2

c	  >>> compute sigma
	  sum = dcmplx(0.0, 0.0)
	  pro = dcmplx(1.0d0, 0.0)
      do i = 1, neigen
	    sum = sum + (dcmplx(1, 0) - lambda(i) / dcmplx(2, 0)) /
     & ((dcmplx(1, 0) - fmass / dcmplx(2, 0)) * lambda(i) + fmass)
c	  >>> WARNING: normalization of the determinant
	    pro = pro * 0.5d0 *
     & ((dcmplx(1, 0) - fmass / dcmplx(2, 0)) * lambda(i) + fmass)
	  end do
	  edetr = dreal(pro)
c	  >>> check the magnitude of the determinant
	  write(*, *) 'edetr = ', edetr
	  sigma = dreal(sum) / dble(nsites)
	  write(*, *) 'sigma = ', sigma

	return
	end
