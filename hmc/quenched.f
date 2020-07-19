C     ******************************************************************
      character*32 function v_integrate()
      v_integrate = 'quenched'
      end
C     ******************************************************************
      subroutine qintegrate(eps, nsteps, p)
C     ******************************************************************
C     **                        **                                    **
C     **  QINTEGRATE            **  I. Hip, 08 May 97                 **
C     **  (Quenched full step)  **  Last modified: 08 May 97          **
C     **                        **                                    **
C     ******************************************************************
C     Integrate "equations of motion" for gauge fields (QUENCHED!)
C     ------------------------------------------------------------------
C     IN (real*8) eps        - epsilon - integration step
C     IN (integer*4) nsteps  - number of steps in one trajectory
C     IN (real*8) p          - conjugate impulses
C     ------------------------------------------------------------------
C     COMMON IN (real*8) beta                   - \beta
C     COMMON INOUT (complex*16) u(nsite, ndim)  - gauge fields
C     ------------------------------------------------------------------
C     Remark: This program calls INVERT, DSGAUGE, DSFERM, (SCALP)
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      real*8 eps, p(nsite, ndim)
      real*8 dsg(nsite, ndim), dsg_new(nsite, ndim)
      real*8 sqrteps, hsqrteps, alpha
      real*8 beta
      complex*16 u(nsite, ndim)
C     ------------------------------------------------------------------
      common /constants/beta
      common /gauge_fields/u
C     ------------------------------------------------------------------
      sqrteps = dsqrt(eps)
      hsqrteps = 0.5d0 * sqrteps
      
      call dsgauge(beta, u, dsg)

      do nst = 1, nsteps

	do nd = 1, ndim
	  do ns = 1, nsite
	    alpha = - 0.5d0 * dsg(ns, nd) * eps + p(ns, nd) * sqrteps
	    u(ns, nd) = u(ns, nd) * dcmplx(dcos(alpha), dsin(alpha))
	  end do
	end do

	call dsgauge(beta, u, dsg_new)

	do nd = 1, ndim
	  do ns = 1, nsite
	    p(ns, nd) = p(ns, nd) - hsqrteps * dsg(ns, nd)
	    dsg(ns, nd) = dsg_new(ns, nd)
	    p(ns, nd) = p(ns, nd) - hsqrteps * dsg(ns, nd)
	  end do
	end do

      end do

      return
      end

