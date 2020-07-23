C     ******************************************************************
      character*32 function v_integrate()
      v_integrate = 'full step'
      end
C     ******************************************************************
      integer*4 function integrate(eps, nsteps, p, xi, phi, chi)
C     ******************************************************************
C     **                    **                                        **
C     **  INTEGRATE         **  Based on Full step by H. Gausterer    **
C     **  (Full step)       **  (Modified by I. Hip, Nov 95)          **
C     **                    **  Last modified: 25 Jan 97              **
C     ******************************************************************
C     Integrate "equations of motion" for gauge fields
C     ------------------------------------------------------------------
C     Part of package:   HMC2DU1
C     ------------------------------------------------------------------
C     Input parameters:  eps, nsteps, p, xi, phi, chi
C     ------------------------------------------------------------------
C     Output parameters: p, phi, chi
C     ------------------------------------------------------------------
C     Output variables:  u(1:nsite, 1:ndim)
C     ------------------------------------------------------------------
C     Remark: This program calls INVERT, DSGAUGE, DSFERM, (SCALP)
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 u(nsite, ndim)
      complex*16 xi(nferm), phi(nferm), chi(nferm)
      complex*16 scalp, escalp, op_ubudbd
      real*8 p(nsite,ndim)
      real*8 dsg(nsite,ndim), dsf(nsite,ndim), dsa(nsite,ndim)
      real*8 beta, akap
      real*8 eps, sqrteps, hsqrteps, alpha
      integer*4 niter, max_niter
C     ------------------------------------------------------------------
      common /kappa/akap 
      common /constants/beta
      common /gauge_fields/u
C     ------------------------------------------------------------------
      max_niter = 1

      sqrteps = dsqrt(eps)
      hsqrteps = 0.5d0 * sqrteps
      
      call dsgauge(beta, u, dsg)
      niter = invert_excp(chi, phi)
      if(niter .lt. 0) then
        integrate = niter
        return
      end if

      call dsferm(akap, u, chi, dsf)

      do nd = 1,ndim
	do ns = 1,nsite
	  dsa(ns,nd) = dsg(ns,nd) + dsf(ns,nd)
	end do
      end do

      do 10 nst = 1, nsteps

	do nd = 1,ndim
	  do ns = 1,nsite
	    alpha = - 0.5d0 * dsa(ns,nd) * eps + p(ns,nd) * sqrteps
	    u(ns,nd) = u(ns,nd) * dcmplx(dcos(alpha), dsin(alpha))
	  end do
	end do

	call dsgauge(beta, u, dsg)
	niter = invert_excp(chi, phi)
	if(niter .lt. 0) then
	  integrate = niter
          return
        else
	  if(niter .gt. max_niter) max_niter = niter
	end if

	call dsferm(akap, u, chi, dsf)

	do nd = 1, ndim
	  do ns = 1, nsite
	    p(ns,nd) = p(ns,nd) - hsqrteps * dsa(ns,nd)
	    dsa(ns,nd) = dsg(ns,nd) + dsf(ns,nd)
	    p(ns,nd) = p(ns,nd) - hsqrteps * dsa(ns,nd)
	  end do
	end do

10    continue
      integrate = max_niter      
      end

