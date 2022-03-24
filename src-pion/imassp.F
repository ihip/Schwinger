      subroutine massdsp(g, dsp, con)      
C     ------------------------------------------------------------------
C	IN complex*16 g(nferm, nferm) - quark propagator matrix
C	OUT real*8 dsp(ntime - 1, 0:nspace / 2, 4) - correlations
C	  needed for dispersion relations
C	OUT real*8 con(0:nspace / 2, 2)
C     ------------------------------------------------------------------
C     03 Nov 96, hip - Last modified: 11 Aug 98, hip
C     ------------------------------------------------------------------
C     The following correlations between time-slices are computed:
C      dsp(index)
C       1 : Tr(sigma_1 G_{xy} sigma_1 G_{yx}) - sigma_1 pion
C       2 : Tr(sigma_1 G_{xx}) Tr(sigma_1 G_{yy}) - vacuum fluctuation
C       3 : Tr(sigma_3 G_{xy} sigma_3 G_{yx}) = Tr(G_{yx}^+ G_{yx})
C	4 : Tr(sigma_3 G_{xx}) Tr(sigma_3 G_{yy})
C      con(index)
C	1 : Tr(sigma_1 G_{xx})
C	2 : Tr(sigma_3 G_{xx})
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (pi2=2*3.1415926535897932)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm)
      real*8 absdg, absnd, mixdg, mixnd, phase

	real*8 dsp(NTIME - 1, 0:NSPACE / 2, 4)
	real*8 con(0:NSPACE / 2, 2)

c	>>> needed for \eta
	real*8 x1, x3, y1, y3

C     ------------------------------------------------------------------
c      common/invers/ g
C     ------------------------------------------------------------------

C     >>> initialization, everything is set to 0
	do ic = 1, 4
	  do ip = 0, NSPACE / 2
            do it = 1, NTIME - 1
	      dsp(it, ip, ic) = 0.0d0
	    end do
	  end do
	end do

c     >>> loop through all time-distances
      do it = 1, NTIME - 1

c	>>> this is not mistake - only last con should be saved
	do ic = 1, 2
          do ip = 0, NSPACE / 2
	    con(ip, ic) = 0.0d0
	  end do
	end do

c     >>> loop through all slices
      do it0 = 0, NTIME - 1

c       >>> loop through the source slice
        do ix = 1, NSPACE

c         >>> loop through the sink slice
	  do j = 1, NSPACE

	    kx = ndrc * (it0 * NSPACE + ix) - 1
	    ky = ndrc * ((mod(it0 + it, NTIME) * NSPACE) + j) - 1
	    idx = ix - j

	    absdg = dreal( dconjg(g(ky, kx)) * g(ky, kx) +
     & dconjg(g(ky + 1, kx + 1)) * g(ky + 1, kx + 1))
	    absnd = dreal( dconjg(g(ky, kx + 1)) * g(ky, kx + 1) +
     & dconjg(g(ky + 1, kx)) * g(ky + 1, kx))
	    mixdg = dreal( dconjg(g(ky + 1, kx + 1)) * g(ky, kx) +
     & dconjg(g(ky, kx)) * g(ky + 1, kx + 1))
	    mixnd = dreal( dconjg(g(ky, kx + 1)) * g(ky + 1, kx) +
     & dconjg(g(ky + 1, kx)) * g(ky, kx + 1))

c	>>> vacuum fluctuations needed for \eta
	x1 = 2.0d0 * dimag(g(kx, kx + 1))
	x3 = g(kx, kx) - g(kx + 1, kx + 1)
	y1 = - 2.0d0 * dimag(g(ky, ky + 1))
	y3 = g(ky, ky) - g(ky + 1, ky + 1)

	do ip = 0, NSPACE / 2

	  phase = dcos(idx * (dble(pi2 * ip) / dble(NSPACE)))

c         >>> Tr(sigma_1 G_{xy} sigma_1 G_{yx})
	    dsp(it, ip, 1) = dsp(it, ip, 1) + (mixdg - mixnd) * phase
c         >>> x1 * y1
	    dsp(it, ip, 2) = dsp(it, ip, 2) + (x1 * y1) * phase
c         >>> Tr(sigma_3 G_{xy} sigma_3 G_{yx})
	    dsp(it, ip, 3) = dsp(it, ip, 3) + (absdg + absnd) * phase
c         >>> x3 * y3
	    dsp(it, ip, 4) = dsp(it, ip, 4) + (x3 * y3) * phase

c         >>> x1
	    con(ip, 1) = con(ip, 1) + x1 * phase
c         >>> x3
	    con(ip, 2) = con(ip, 2) + x3 * phase
	end do

	  end do
	end do
      end do

      end do

c	>>> normalize everything to number of sites
	do ic = 1, 4
          do ip = 0, NSPACE / 2
	    do it = 1, NTIME - 1
	      dsp(it, ip, ic) = dsp(it, ip, ic) / dble(nsite)
	    end do
	  end do
	end do

	do ic = 1, 2
          do ip = 0, NSPACE / 2
	    con(ip, ic) = con(ip, ic) / dble(nsite)
	  end do
	end do

      return
      end

