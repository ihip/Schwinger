C     ******************************************************************
      integer*4 function update(eps, nsteps, iacc, s_final, h_delta)
C     ******************************************************************
C     **                    **                                        **
C     **  UPDATE            **  I.Hip, 08 May 97                      **
C     **  (quenched)        **  Last modified: ivh, 08 May 97         **
C     **                    **                                        **
C     ******************************************************************
C     Updating the gauge fields with the HMC algorithm (QUENCHED!)
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
C     Remark: This program calls GAUSS, NORM, SGAUGE, QINTEGRATE,
C       CBL_RC_LOG
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      real*8 eps, s_final, h_delta
      real*8 p(nsite,ndim), sgauge, norm
      real*8 h_initial, h_final, s_initial
      real*4 rando(ngaug)
      real*8 beta
      complex*16 u(nsite, ndim), ukeep(nsite, ndim)
C     ------------------------------------------------------------------
      common/constants/ beta
      common/gauge_fields/ u
      save ukeep
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

C     >>> compute: h_initial = (p, p) / 2 + S_G(beta, U)
      s_initial = sgauge(beta, u)
      h_initial = 0.5d0 * norm(p, ngaug) + s_initial

      call qintegrate(eps, nsteps, p)

C     >>> compute: h_final = (p, p) / 2 + S_G(beta, U)
      s_final = sgauge(beta, u)        
      h_final = 0.5d0 * norm(p, ngaug) + s_final

C     >>> METROPOLIS TEST <<<
      h_delta = h_final - h_initial
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

c     >>> for compatibility reasons returns niter - always 1!
      update = 1
      return
      end
