c Subroutines for computing eigenvalues of M^+ M
c Ivan Hip, 02 Mar 96, Varazhdin
c Last modified: 04 Apr 97, Graz 
c
c - three smallest eigenvalues are computed and returned to main
c   program in single precision
c
c WARNING - in subroutine eigen there are two adjustable parameters 
c which are needed in subroutine bisec; these are lbd and ubd which
c means lower and upper bound of the interval in which the eigenvalues 
c are searched for; to find required eigenvalues they should be 
c properly tuned - lbd = 0.0d0 is reasonable if the M^+ M is positive
c definite (but for check I put lbd = -1.0d-2), and ubd is now set
c to 1.0d-1; to set this bound higher means no error, but makes
c program more inefficient because more eigenvalues will be computed
c (which makes no sense if we need the three smallest ones)
c
c TO DO:
c - to check against standard eigenvalue algorithms 
c - to implement automatical adjusting of lbd and ubd parameters
c - to compute eigenvalues of D, M and/or Q matrices...
c
c----------------------------------------------------------------------
c Description of routines:
c
c   integer*4 eigen(eig) - uses no input parameters, and gives (up to)
c     three smallest eigenvalues in single precision as output; integer
c     number which is returned gives the number of returned eigenvalues
c     in eig (up to 3); if this number is smaller than 3, that means
c     that lbd and ubd parameters should be tuned (or some more serious
c     problem...);
c     - only this function should be called from main program!
c
c   lanczos(niter, a, b) - computes tridiagonal Lanczos matrix using
c     Lanczos recursion - niter is input parameter and defines 
c     dimension of the tridiagonal matrix; a and b are diagonal and
c     subdiagonal of the tridiagonal Lanczos matrix
c
c   bisec(...) - subroutine bisec computes eigenvalues of the
c     tridiagonal Lanczos matrix by method of bisection and recognizes
c     which of them are spurious - for detailed description of 
c     parameters see the comments in the subroutine
c  
c   find_tkmax(...) - function which is needed with bisec
c  
c----------------------------------------------------------------------
      
C     ******************************************************************
      integer*4 function eigen(eig)
C     ******************************************************************
C     **                    **                                        **
C     **  Eigen             **    Program by I. Hip, 02. Mar '96      **
C     **                    **    Last modified: 23. May '96                                    **
C     ******************************************************************
C     Computes 3 smallest eigenvalues of squared Wilson Dirac operator
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:   
C     ------------------------------------------------------------------
C     Output variables: eig - array of 3 smallest eigenvalues in single 
C                             precision
C     ------------------------------------------------------------------
C     USE lanczos, bisec, find_tkmax
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      real*4 eig(3)
      real*8 a(nall), b(nall + 1), e(nall)
      real*8 b2(nall + 1), tkmax, vb(2 * nall), find_tkmax, lbd, ubd
      integer*4 mp(nall)

      call lanczos(nall, a, b)
      
c WARNING: parameters lbd and ubd should be properly tuned!      
c (* old values which were used for 16^2 lattice:
c      lbd = 0.0d0
c      ubd = 2.0d-2 *)
c (* old values which were used for 8^2 lattice:
c      lbd = -1.0d2
c      ubd = 1.0d-1
      lbd = 0.0d0
      ubd = 1.5d0

      tkmax = find_tkmax(a, b, nall)
      call bisec(a, b, b2, lbd, ubd, nall, tkmax, vb, ndis, e, mp, ic)

      j = 1
      do i = 1, ndis
	if(mp(i) .gt. 0) then
	  eig(j) = e(i)
	  if(j .eq. 3) go to 99
	  j = j + 1
	end if
      end do

c less than 3 ev-s found in given interval - missing eig entries
c are filled with non-physical value 
      do i = j, 3
        eig(i) = -13.0d13
      end do

99    eigen = j
      return
      end


C     ******************************************************************
      subroutine lanczos(niter, a, b)
C     ******************************************************************
C     **                    **                                        **
C     **  Lanczos           **    Program by I. Hip, 08 Feb 96        **
C     **                    **    Last modified: 04 Apr 97            **
C     ******************************************************************
C     Computes tridiagonal matrix using Lanczos recursion
C     ------------------------------------------------------------------
C     Part of package:   hmc2dfu1
C     ------------------------------------------------------------------
C     Input variables:   niter - number of iterations
C     ------------------------------------------------------------------
C     Output variables:  a - diagonal of tridiagonal matrix
C                        b - subdiagonal of tridiagonal matrix
C     ------------------------------------------------------------------
C     USES mvec, mtvec - matrix multiplication routines
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      real*8 a(nall), b(nall)
      real*8 u(nall), v(nall), t(nall), tt(nall)
      real*8 alpha, beta, tmp
C     ------------------------------------------------------------------
C     INITIALIZATION:
C     -----------------------------------------------------------------
      iter = 1

      do i = 1, nall
	u(i) = 1.0d0 / dsqrt(dble(nall))
	v(i) = 0.0d0
      end do

C     -----------------------------------------------------------------
C     ITERATION:
C     -----------------------------------------------------------------
10    call mvec(u, t)
      call mtvec(t, tt)

      alpha =  0.0d0
      do i = 1, nall
	v(i) = v(i) + tt(i)
	alpha = alpha + v(i) * u(i)
      end do

      beta = 0.0d0
      do i = 1, nall
	v(i) = v(i) - alpha * u(i)
	beta = beta + v(i) * v(i)
      end do
      beta = dsqrt(beta)

      a(iter) = alpha
c      write(*, *) iter, alpha, beta
      iter = iter + 1
      if(iter .gt. niter) go to 99
      b(iter) = beta

      do i = 1, nall
	tmp = u(i)
	u(i) = v(i) / beta
	v(i) = - beta * tmp
      end do

      go to 10

99    continue

      end


C     ******************************************************************
      subroutine bisec(alpha, beta, beta2, lbd, ubd, mev, tkmax, vb,
     &  ndis, vs, mp, ic)
C     ******************************************************************
C     **                    **                                        **
C     **  Bisec             **    Program by I. Hip, Feb '96          **
C     **                    **                                        **
C     ******************************************************************
      double precision alpha(*), beta(*), beta2(*), lbd, ubd, tkmax,
     &  vb(*), vs(*)
      integer mev, ndis, mp(*), ic
C     ******************************************************************
c     IN:   alpha - diagonal of 3-diagonal matrix
c     IN:   beta - subdiagonal of 3-diagonal matrix
c     TMP:  beta2 - beta ** 2 stored in beta2
c     IN:   lbd - lower boundary of given interval
c     IN:   ubd - upper boundary of given interval
c     IN:   mev - part of T matrix to be diagonalized (e.g. mev = 20,
c             first 20 alphas and betas used to compute eigenvalues)
c     IN:   tkmax - biggest element in 3-diagonal matrix (abs. value)
c     TMP:  vb - upper and lower bounds in Sturm sequencing
c     OUT:  ndis - how many distinct eigenvalues computed
c     OUT:  vs - eigenvalues
c     OUT:  mp - information about eigenvalues which are returned in vs:
c             mp(i) = 0  -> i-th ev is spurious
c             mp(i) = 1  -> i-th ev is good
c             mp(i) = n (n > 1) -> n-multiple but good
c     OUT:  ic - total number of Sturms used (whatever that means...)
C     ******************************************************************
c     This function is slightly modified func. bisec which is published
c     in Cullum and Willoughby: Lanczos algorithms for large symmetric
c     eigenvalues computations, Vol.2. Programs, Birkhauser, Boston 1985
c
c     COMMENTS ON DIFFERENCES:
c     - without nint - that means we compute in only one interval:
c       (lbd, upd]
c     - without eps - defined in subroutine
c     - instead of ttol, tmax is used as a parameter - it could be
c       computed using function find_tkmax
c     - ic on input should set maximal number of sturms, but that makes
c       not to much sense to limit
C     ******************************************************************

	double precision eps, ttol, betam, temp, ep0, ep1, lb, ub,
     &    x1, yu, yv, xl, xu, xs, x0, ept
	integer mp1, na, md, ng, ict, isturm, nev, mt

c eps - 2 * machine precision (mach. prec. on SGI: 1.2D-16)
	eps = 2.0d0 * 1.2d-16
	ttol = eps * tkmax
c        write(*, *) 'tkmax = ', tkmax, ' ttol = ', ttol

	ndis = 0
	ic = 0
	mp1 = mev + 1
	betam = beta(mp1)
	beta(mp1) = 0.0
	do i = 1, mp1
	  beta2(i) = beta(i) * beta(i)
	end do
	
	temp = dfloat(mev + 1000)
	ep0 = temp * ttol
	ep1 = dsqrt(temp) * ttol
c        write(*, *) 'ep0 = ', ep0
	
	lb = lbd
	ub = ubd

	na = 0
	md = 0
	ng = 0
c ict is total Sturm count on (lb, ub)
	ict = 0

c start of T-eigenvalues calculations
	x1 = ub
	isturm = 1
	go to 330

50      na = nev
	x1 = lb
	isturm = 2
	go to 330

60      mt = nev
	ict = ict + 2
c mt is number of ev-s in given interval
c        write(*, *) mt, ' ev-s found in interval, ', na, ' ev-s > upb'

c estimate number of sturms
c       iest = 30 * mt

	if(mt .eq. 0) then
	  write(*, *) 'No eigenvalues in given interval.'
	  go to 430
	end if

	do i = 1, mt
	  vb(i) = lb
	  mti = mt + i
	  vb(mti) = ub
	end do

	k = mt
150     continue
	ico = 0
	xl = vb(k)
	mtk = mt + k
	xu = vb(mtk)

	isturm = 3
	x1 = xu
	ico = ico + 1
	go to 330

160     nu = nev
	isturm = 4
170     continue
	x1 = (xl + xu) * 0.5d+00
	xs = dabs(xl) + dabs(xu)
	x0 = xu - xl
	ept = eps * xs + ep1

	if(x0 .le. ept) go to 230

	ico = ico + 1
	go to 330

180     continue
	if(nev .lt. k) go to 190
	xl = x1
	go to 170

190     continue
	xu = x1
	nu = nev
	if(nev .ne. 0) then
	  do i = 1, nev
	    vb(i) = dmax1(x1, vb(i))
	  end do
	end if
	nev1 = nev + 1

	do ii = nev1, k
	  i = mt + ii
	  vb(i) = dmin1(x1, vb(i))
	end do
	go to 170

230     continue
	ndis = ndis + 1
	md = md + 1
	vs(ndis) = x1

	jsturm = 1
	x1 = xl - ep0
	go to 370

240     kl = kev
	jl = jev

	jsturm = 2
	ico = ico + 2
	x1 = xu + ep0
	go to 370

250     ju = jev
	ku = kev

c is this simple T-eigenvalue?
	if((kl - ku - 1) .eq. 0) go to 290
c if multiple than is also good (?)
	if(ku .eq. nu) go to 280

260     continue
	isturm = 5
	x1 = x1 + ep0
	ico = ico + 1
	go to 330

270     kne = ku - nev
	ku = nev
	if(kne .ne. 0) go to 260
	
c compute ev multiplicity: mp(..)
280     mpev = kl - ku
	knew = ku
	go to 300

290     continue
	mpev = 1
	if(ju .lt. jl) mpev = 0
	knew = k - 1

300     k = knew
	mp(ndis) = mpev
	if(mpev .ge. 1) ng = ng + 1
	ict = ict + ico

	if(k .le. 0) go to 410

	do i = 1, knew
	  vb(i) = dmax1(x1, vb(i))
	end do
	go to 150

330     nev = -na
	yu = 1.0d+00

	do i = 1, mev
	  if(yu .ne. 0.0) then
	    yv = beta2(i) / yu
	  else
	    yv = beta(i) / eps
	  end if
	  yu = x1 - alpha(i) - yv
	  if(yu .lt. 0.0) nev = nev + 1
	end do

	go to (50, 60, 160, 180, 270), isturm

370     kev = -na
	yu = 1.0d+00

	do ii = 1, mev
	  i = mp1 - ii
	  if(yu .ne. 0.0) then
	    yv = beta2(i + 1) / yu
	  else
	    yv = beta(i + 1) / eps
	  end if
	  yu = x1 - alpha(i) - yv
	  jev = 0
	  if(yu .lt. 0.0) then
	    kev = kev + 1
	    jev = 1
	  end if
	end do
	jev = kev - jev
	go to (240, 250) jsturm

410     ic = ict + ic

430     beta(mp1) = betam
	return
	end


C     ******************************************************************
      double precision function find_tkmax(alpha, beta, mev)
C     ******************************************************************
	  double precision alpha(*), beta(*)
	  integer mev

	double precision tkmax

	tkmax = 0.0
	do i = 1, mev
	  if(dabs(alpha(i)) .gt. tkmax) tkmax = alpha(i)
	  if(beta(i) .gt. tkmax) tkmax = beta(i)
	end do

	find_tkmax = tkmax
	end
