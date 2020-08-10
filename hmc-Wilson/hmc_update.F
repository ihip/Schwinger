C     ******************************************************************
      integer*4 function update(eps, nsteps, iacc, s_final, h_delta)
C     ******************************************************************
C     **                    **                                        **
C     **  UPDATE            **  Based on UPDATE by H. Gausterer       **
C     **  (HMC)             **  (Modified by I. Hip, Nov 95)          **
C     **                    **  Last modified: ivh, 08 May 97         **
C     ******************************************************************
C     Updating the gauge fields with the HMC algorithm
C     ------------------------------------------------------------------
C     IN (real*8) eps        - epsilon - integration step
C     IN (integer*4) nsteps  - number of steps in one trajectory
C     OUT (integer*4) iacc   - result of Metropolis test, new
C                              configuration is: 1 accepted, 0 rejected
C     OUT (real*8) s_final   - final action
C     OUT (real*8) h_delta   - energy change
C     ------------------------------------------------------------------
C     COMMON IN (real*8) beta                   - \beta
C     COMMON INOUT (complex*16) u(nsite, ndim)  - gauge fields
C     ------------------------------------------------------------------
C     Remark: This program calls GAUSS, CGAUSS, MTVEC, NORM, SGAUGE, 
C       INTEGRATE, SCALP, CBL_RC_LOG
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      real*8 eps, s_final, h_delta
      real*8 p(nsite, ndim), sgauge, norm
      real*8 h_initial, h_final, s_initial
      complex*16 xi(nferm), phi(nferm), chi(nferm), scalp
C     >>> dimension of rando is max(2*nferm, ngaug)
      real rando(2*nferm)
      real*8 beta 
      complex*16 u(nsite, ndim), ukeep(nsite, ndim)
C     ------------------------------------------------------------------
      common/constants/ beta
      common/gauge_fields/ u
      save ukeep
      data chi/nferm*(0.0d0, 0.0d0)/
C     ------------------------------------------------------------------
C     >>> keep the old gauge fields in ukeep
      do nd = 1, ndim
	do ns = 1, nsite
	  ukeep(ns, nd) = u(ns, nd)
	end do
      end do
	  
C     >>> initial conjugate impulses p should have Gaussian distribution 
C     >>> (p's are conjugate to vector potential A and they are real!)
      call gauss(p, ngaug, rando)

C     >>> to construct pseudofermion we need some Gaussian "noise" xi
C     >>> (xi should be complex and < (xi, xi) > = 1)
      call cgauss(xi, nferm, rando)

C     >>> pseudofermion: phi = (M^+) xi
      call mtvec(xi, phi)

C     >>> compute: h_initial = (p, p) / 2 + S_G(beta, U) + (xi, xi)
      s_initial = sgauge(beta, u) + norm(nall, xi)
      h_initial = 0.5d0 * norm(ngaug, p) + s_initial

      niter = integrate(eps, nsteps, p, xi, phi, chi)
C     >>> if invert didn't converge, reject the new configuration and
C     >>> restore the original one
      if(niter .lt. 0) then
	do nd = 1, ndim
	  do ns = 1, nsite
	    u(ns, nd) = ukeep(ns, nd)
	  end do
	end do
	iacc = 0
C       >>> restore original s
        s_final = s_initial
        h_delta = 0.0d0
	hmc_update = niter
        return
      end if

C     >>> compute: h_final = (p, p) / 2 + S_G(beta, U) + (phi, chi)
      s_final = sgauge(beta, u) + dreal(scalp(nferm, phi, chi))        
      h_final = 0.5d0 * norm(ngaug, p) + s_final

C     >>> METROPOLIS TEST <<<
      h_delta = h_final - h_initial
c	write(*, *) 'h_delta = ', h_delta
      call cbl_rc_log(ran, 1, 0.0)
      if(h_delta .gt. -ran) then
C     >>> new configuration is rejected - restore the original one
	do nd = 1, ndim
	  do ns = 1, nsite
	    u(ns, nd) = ukeep(ns, nd)
	  end do
	end do
C       >>> restore original s
        s_final = s_initial
	iacc = 0
      else
	iacc = 1
      endif

      update = niter
      return
      end
