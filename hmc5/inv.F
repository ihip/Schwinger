C approved to v3 - 29 Dec 96 - ivh (without imassp!!!)

C     ******************************************************************
      character*64 function v_source()
      v_source = 'point source'
      return
      end
C     ******************************************************************
      integer*4 function inv()
C     ******************************************************************
C     **                    **                                        **
C     **  INV               **  Ivan Hip, 15 May 96                   **
C     **                    **  Last modified: 21 Feb 97              **
C     ******************************************************************
C     Computes full inverse of the fermion matrix
C     ------------------------------------------------------------------
C     Part of package:   HMC2DU1
C     ------------------------------------------------------------------
C     Input parameters:  
C     ------------------------------------------------------------------
C     Output parameters: 
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), s(nferm)
      integer*4 nexcp, nexcp_b, nexcp1, nexcp1_b
C     ------------------------------------------------------------------
      common/invers/ g
      common/excpcomm/ nexcp, nexcp_b, nexcp1, nexcp1_b
C     ------------------------------------------------------------------
      max_niter1 = 0
      
      do k = 1, nferm

c	do i = 1, nferm
c	  s(i) = 0.0d0
c	end do
c	s(k) = 1.0d0

	call source(nferm, k, s)
c	call source2(nferm, k, s)

      niter1 = invert1(g(1, k), s)
c      print *, 'inv: niter1 = ', niter1

	if(niter1 .lt. 0) then
          nexcp1 = nexcp1 + 1
	  do i = 1, nferm
 	    g(i, k) = 1.0d0
	  end do
          niter1 = invert1_b(g(1, k), s)
          if(niter1 .lt. 0) then
            nexcp1_b = nexcp1_b + 1
            call save_gfield(nexcp1_b, 1)
            if(niter1 .eq. -1) niter1 = -3
            if(niter1 .eq. -2) niter1 = -4
	    inv = niter1
	    return
          end if
        end if
	if(niter1 .gt. max_niter1) max_niter1 = niter1
      end do
      
      inv = max_niter1
      return
      end


      subroutine source(n, k, s)
      complex*16 s(*)

	do i = 1, n
	  s(i) = 0.0d0
	end do
	s(k) = 1.0d0

      return
      end


      subroutine source2(nferm, k, s)
      complex*16 s(*)
	real*8 r

	ir = 1
c	r = 1.0d0 / dsqrt(dble(2 * ir + 1))
c	r = 1.0d0 / dble(2 * ir + 1)
	r = 1.0d0
	do i = 1, nferm
	  s(i) = 0.0d0
	end do
	n = 2 * NSPACE
	do i = -2 * ir, 2 * ir, 2
	  s(n * ((k - 1) / n) + mod(k + i + n - 1, n) + 1) = r
	end do

      return
      end


      complex*16 function op_pbp_e()
C     ------------------------------------------------------------------
C     computes chiral condensate (ivh, 16 May 96)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------
      common/invers/ g
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm
	sum = sum + g(i, i) 
      end do
      op_pbp_e = sum

      return
      end

      complex*16 function op_pbg5p_e()
C     ------------------------------------------------------------------
C     computes \bar{\Psi} \gamma_5 \Psi (ivh, 18 May 96)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------
      common/invers/ g
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
	sum = sum + (g(i, i) - g(i + 1, i + 1))
      end do
      op_pbg5p_e = sum

      return
      end

      complex*16 function op_pbg1p_e()
C     ------------------------------------------------------------------
C     computes \bar{\Psi} \gamma_1 \Psi (ivh, 18 May 96)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------
      common/invers/ g
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
c	sum = sum + g(i, i + 1) + dconjg(g(i + 1, i))
	sum = sum + g(i, i + 1) + g(i + 1, i)
      end do
      op_pbg1p_e = sum

      return
      end

      complex*16 function op_pbg2p_e()
C     ------------------------------------------------------------------
C     computes \bar{\Psi} \gamma_2 \Psi (ivh, 18 May 96) LM: 21 Mar 97
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------
      common/invers/ g
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
c        sum = sum - dcmplx(0.0d0, 1.0d0) * (dconjg(g(i, i + 1)) +
c     &    g(i + 1, i))
        sum = sum - dcmplx(0.0d0, 1.0d0) * (g(i, i + 1) +
     &    g(i + 1, i))
      end do
      op_pbg2p_e = sum

      return
      end

      complex*16 function op_ubudbd_e()
C     ------------------------------------------------------------------
C     computes $ \bar{u} u \bar{d} d $ (ivh, 29 May 96)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------
      common/invers/ g
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
	sum = sum + g(i, i)**2 + 2.0d0 * g(i, i) * g(i + 1, i + 1) +
     &    g(i + 1, i + 1)**2
      end do
      op_ubudbd_e = sum

      return
      end


      complex*16 function op_ubg5udbg5d_e()
C     ------------------------------------------------------------------
C     computes $ \bar{u} u \bar{d} d $ (ivh, 16 Jun 97)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------
      common/invers/ g
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
	sum = sum + g(i, i)**2 - 2.0d0 * g(i, i) * g(i + 1, i + 1) +
     &    g(i + 1, i + 1)**2
      end do
      op_ubg5udbg5d_e = sum

      return
      end


      complex*16 function op_pbppbp_e()
C     ------------------------------------------------------------------
C     computes $ \bar{\Psi} \Psi \bar{\Psi} \Psi $ (ivh, 29 May 96)
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm), sum
C     ------------------------------------------------------------------
      common/invers/ g
C     ------------------------------------------------------------------

      sum = dcmplx(0.0d0, 0.0d0)
      do i = 1, nferm, 2
	sum = sum + g(i, i)**2 + 4.0d0 * g(i, i) * g(i + 1, i + 1) +
     &    g(i + 1, i + 1)**2
      end do
      op_pbppbp_e = sum

      return
      end

      real*8 function op_linkc_e()
C     ------------------------------------------------------------------
C     computes "link condensate" (ivh, 16 May 96): \sum_{x,y}
C       (\bar{\Psi} (1 + \sigma_{\mu}) U_{x \mu} \Psi_{x + \mu} + h.c.)
C     WARNING! In this subroutine is assumed that ndrc = 2
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      complex*16 g(nferm, nferm)
      complex*16 u(nsite, ndim)
      integer ind(4, nsite)
      real*8 sum
C     ------------------------------------------------------------------
      common/gauge_fields/ u
      common/index_arrays/ ind
      common/invers/ g
C     ------------------------------------------------------------------
      sum = 0.0d0

      do ix = 1, nsite
	kx = ndrc * ix
	
	ixm = ind(1, ix)
	kxm = ndrc * ixm
	sum = sum + dreal( u(ix, 1) * ( g(kxm - 1, kx - 1) + 
     &    g(kxm - 1, kx) + g(kxm, kx - 1) + g(kxm, kx) ) )
	ixm = ind(2, ix)
	kxm = ndrc * ixm
	sum = sum + dreal(u(ix, 2) * (g(kxm - 1, kx - 1) + g(kxm, kx) + 
     &    dcmplx(0, 1) * (g(kxm - 1, kx) - g(kxm, kx - 1))))
      end do

      op_linkc_e = sum
      return
      end


C     ------------------------------------------------------------------
      character*64 function v_direction()
      v_direction = 'time'
      return
      end
C     ------------------------------------------------------------------
      subroutine imass(f4)
C     ------------------------------------------------------------------
C     v3 - ivh - 29 Dec 96 - Last modified: 07 Aug 97
C     ------------------------------------------------------------------
C     Computes current correlations and PCAC
C     - propagation between time-slices -
C     - summation over all time slices - because of symmetry, only 
C       NTIME / 2 (instead of NTIME - 1) entries should be keeped;
C       however, PCAC isn't symmetric (?), and we also want to be
C       compatible with old (or any other unsymmetric) measurements
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (mcomp=20)
C     ------------------------------------------------------------------
C INOUT real*4 array f4:
      real*4 f4(NTIME - 1, mcomp)

C GLOBAL gauge field (in common /gauge_fields/):
      complex*16 u(nsite, ndim)

C GLOBAL inverted matrix - result of inv (in common /invers/)
      complex*16 g(nferm, nferm)

C GLOBAL index array (in common /index_arrays/)
      integer ind(4, nsite)

C double precision equivalent of INOUT parameter f
      real*8 f(NTIME - 1, mcomp)

C t-s are triplet correlations Tr(sigma_i G_{yx} sigma_i G_{yx}):
      real*8 t0, t1, t2, t3

C v-s are vacuum fluctuations Tr(sigma_i G_{xx}) Tr(sigma_i G_{yy}):
      real*8 v0, v1, v2, v3

C x-s are Tr(sigma_i G_{xx}):
      real*8 x0, x1, x2, x3, cx0, cx1, cx2, cx3

C y-s are Tr(sigma_i G_{yy}):
      real*8 y0, y1, y2, y3, cy0, cy1, cy2, cy3

C PCAC (sp - sigma1 PCAC, sr - sigma1 r-corrections, t - triplet):
      real*8 sp, sr, tp, tr

C          | a   b |
C G_{xy} = |       |
C          | c   d |
      complex*16 a, b, c, d

C needed for triplet correlations computations:
      real*8 absdg, absnd, mixdg, mixnd

C needed for PCAC computations
      real*8 s1, s2, p1, p2, r
      real*8 s11, s12, s13, s14, p11, p12, p13, p14

c needed for PCAC with local currents
      complex*16 a1, b1, c1, d1, a2, b2, c2, d2
      complex*16 a3, b3, c3, d3, a4, b4, c4, d4
      real*8 spa, spb, spc

C     ------------------------------------------------------------------
      common/gauge_fields/ u
      common/index_arrays/ ind
      common/invers/ g
C     ------------------------------------------------------------------

C Wilson parameter r is set to 1
      r = 1.0d0

C loop through all (time) distances:
      do it = 1, NTIME - 1

C initialization with zero:
        t0 = 0.0d0
        t1 = 0.0d0
        t2 = 0.0d0
        t3 = 0.0d0
        v0 = 0.0d0
        v1 = 0.0d0
        v2 = 0.0d0
        v3 = 0.0d0
        cx0 = 0.0d0
        cx1 = 0.0d0
        cx2 = 0.0d0
        cx3 = 0.0d0
        cy0 = 0.0d0
        cy1 = 0.0d0
        cy2 = 0.0d0
        cy3 = 0.0d0
        sp = 0.0d0
        sr = 0.0d0
        tp = 0.0d0
        tr = 0.0d0

	spa = 0.0d0
	spb = 0.0d0
	spc = 0.0d0

	s11 = 0.0d0
	s12 = 0.0d0
	s13 = 0.0d0
	s14 = 0.0d0
	p11 = 0.0d0
	p12 = 0.0d0
	p13 = 0.0d0
	p14 = 0.0d0

C loop through all slices
      do it0 = 0, NTIME - 1

C   loop through the "source" time slice:
        do ix = 1, NSPACE

C     loop through the "sink" time slice:
          do j = 1, NSPACE

        ixt0 = it0 * NSPACE + ix
        iyt0 = mod(it0 + it, NTIME) * NSPACE + j
 
        kx = ndrc * ixt0
        ky = ndrc * iyt0

	ix1 = ind(1, ixt0)
	ix2 = ind(2, ixt0)
	ix3 = ind(3, ixt0)
	ix4 = ind(4, ixt0)

	kx1 = ndrc * ix1
	kx2 = ndrc * ix2
	kx3 = ndrc * ix3
	kx4 = ndrc * ix4

 	a = g(kx - 1, ky - 1)
	b = g(kx - 1, ky)
	c = g(kx, ky - 1)
	d = g(kx, ky)

C computation of triplet correlations (what is computed is:
C -Tr(\gamma_i G_yx \gamma_i G_xy); triplet propagator: -2 * t#
	absdg = dreal(dconjg(a) * a + dconjg(d) * d)
	absnd = dreal(dconjg(b) * b + dconjg(c) * c)
	mixdg = dreal(dconjg(a) * d + a * dconjg(d))
        mixnd = dreal(dconjg(b) * c + b * dconjg(c))
	t0 = t0 + (absdg - absnd)
	t1 = t1 + (mixdg - mixnd)
	t2 = t2 + (mixdg + mixnd)
	t3 = t3 + (absdg + absnd)

C computation of vacuum fluctuations (what is computed is:
C Tr(\gamma_i G_xx) * Tr(\gamma_i G_yy); to compute singlet propagator
C following should be used: 4 * v# - 2 * t# - 4 * cx# * cy#
C WARNING: x1 = 2 * I * Im(b_xx) and y1 = 2 * I * Im(b_yy) but instead
C   of I-s in both expressions we take x1 positive and y1 negative,
C   that means that x1 * y1 is OK (although x1 and y1 alone are not!)
C   - the same is true for x2 and y2 - 
	x0 = g(kx - 1, kx - 1) + g(kx, kx)
	x1 = 2.0d0 * dimag(g(kx - 1, kx))
	x2 = 2.0d0 * dreal(g(kx - 1, kx))
	x3 = g(kx - 1, kx - 1) - g(kx, kx)

	cx0 = cx0 + x0
	cx1 = cx1 + x1
	cx2 = cx2 + x2
	cx3 = cx3 + x3

	y0 = g(ky - 1, ky - 1) + g(ky, ky)
	y1 = - 2.0d0 * dimag(g(ky - 1, ky))
	y2 = - 2.0d0 * dreal(g(ky - 1, ky))
	y3 = g(ky - 1, ky - 1) - g(ky, ky)

        cy0 = cy0 + y0
	cy1 = cy1 + y1
	cy2 = cy2 + y2
	cy3 = cy3 + y3

	v0 = v0 + x0 * y0
	v1 = v1 + x1 * y1
	v2 = v2 + x2 * y2
	v3 = v3 + x3 * y3

C PCAC for triplet current:
C - tp is PCAC
C - tr is r-correction

c this was experiment with \sigma_1 current instead \sigma_3...
c        s11 = s11 + dreal( u(ixt0, 1) * (
c     &  - dconjg(d) * g(kx1 - 1, ky - 1) - dconjg(b) * g(kx1, ky - 1) +
c     &    dconjg(c) * g(kx1 - 1, ky) + dconjg(a) * g(kx1, ky) ) )
c        s12 = s12 + dreal( u(ixt0, 2) * dcmplx(0.0d0, 1.0d0) * (
c     &  - dconjg(d) * g(kx2 - 1, ky - 1) + dconjg(b) * g(kx2, ky - 1) +
c     &    dconjg(c) * g(kx2 - 1, ky) - dconjg(a) * g(kx2, ky) ) ) 
c        s13 = s13 + dreal( u(ix3, 1) * (
c     &    d * dconjg(g(kx3 - 1, ky - 1)) + b * dconjg(g(kx3, ky - 1)) -
c     &    c * dconjg(g(kx3 - 1, ky)) - a * dconjg(g(kx3, ky)) ) ) 
c        s14 = s14 + dreal( u(ix4, 2) * dcmplx(0.0d0, 1.0d0) * (
c     &  - d * dconjg(g(kx4 - 1, ky - 1)) + b * dconjg(g(kx4, ky - 1)) +
c     &    c * dconjg(g(kx4 - 1, ky)) - a * dconjg(g(kx4, ky)) ) ) 

	p1 =
     &	dreal( u(ixt0, 1) * (
     &    dconjg(c) * g(kx1 - 1, ky - 1) + dconjg(a) * g(kx1, ky - 1) +
     &    dconjg(d) * g(kx1 - 1, ky) + dconjg(b) * g(kx1, ky) ) ) +
     &  dreal( u(ixt0, 2) * dcmplx(0.0d0, 1.0d0) * (
     &    dconjg(c) * g(kx2 - 1, ky - 1) - dconjg(a) * g(kx2, ky - 1) +
     &    dconjg(d) * g(kx2 - 1, ky) - dconjg(b) * g(kx2, ky) ) ) -
     &  dreal( u(ix3, 1) * (
     &    c * dconjg(g(kx3 - 1, ky - 1)) + a * dconjg(g(kx3, ky - 1)) +
     &    d * dconjg(g(kx3 - 1, ky)) + b * dconjg(g(kx3, ky)) ) ) -
     &  dreal( u(ix4, 2) * dcmplx(0.0d0, 1.0d0) * (
     &    a * dconjg(g(kx4, ky - 1)) - c * dconjg(g(kx4 - 1, ky - 1)) +
     &    b * dconjg(g(kx4, ky)) - d * dconjg(g(kx4 - 1, ky)) ) )

c unclear what was the purpose of this part of the code...
c      p11 = p11 + dreal( u(ixt0, 1) * (
c     &    dconjg(c) * g(kx1 - 1, ky - 1) + dconjg(a) * g(kx1, ky - 1) +
c     &    dconjg(d) * g(kx1 - 1, ky) + dconjg(b) * g(kx1, ky) ) ) 
c      p12 = p12 + dreal( u(ixt0, 2) * dcmplx(0.0d0, 1.0d0) * (
c     &    dconjg(c) * g(kx2 - 1, ky - 1) - dconjg(a) * g(kx2, ky - 1) +
c     &    dconjg(d) * g(kx2 - 1, ky) - dconjg(b) * g(kx2, ky) ) ) 
c      p13 = p13 + dreal( u(ix3, 1) * (
c     &    c * dconjg(g(kx3 - 1, ky - 1)) + a * dconjg(g(kx3, ky - 1)) +
c     &    d * dconjg(g(kx3 - 1, ky)) + b * dconjg(g(kx3, ky)) ) ) 
c      p14 = p14 + dreal( u(ix4, 2) * dcmplx(0.0d0, 1.0d0) * (
c     &    a * dconjg(g(kx4, ky - 1)) - c * dconjg(g(kx4 - 1, ky - 1)) +
c     &    b * dconjg(g(kx4, ky)) - d * dconjg(g(kx4 - 1, ky)) ) )

	p2 =
     &  dreal( u(ixt0, 1) * (
     &    dconjg(a) * g(kx1 - 1, ky - 1) + dconjg(c) * g(kx1, ky - 1) +
     &    dconjg(b) * g(kx1 - 1, ky) + dconjg(d) * g(kx1, ky) ) ) +
     &  dreal( u(ixt0, 2) * (
     &    dconjg(a) * g(kx2 - 1, ky - 1) + dconjg(c) * g(kx2, ky - 1) +
     &    dconjg(b) * g(kx2 - 1, ky) + dconjg(d) * g(kx2, ky) ) ) +
     &  dreal( u(ix3, 1) * (
     &    a * dconjg(g(kx3 - 1, ky - 1)) + c * dconjg(g(kx3, ky - 1)) +
     &    b * dconjg(g(kx3 - 1, ky)) + d * dconjg(g(kx3, ky)) ) ) +
     &  dreal( u(ix4, 2) * (
     &    a * dconjg(g(kx4 - 1, ky - 1)) + c * dconjg(g(kx4, ky - 1)) +
     &    b * dconjg(g(kx4 - 1, ky)) + d * dconjg(g(kx4, ky)) ) )

c also \sigma_1 experiment...
c	s1 = s11 + s12 - s13 - s14
c	sp = sp + s1 / 2.0d0
C >>> future expansion: r-correction for sigma1 current
C >>> temporary set to zero
c	sr = 0.0d0
	
	tp = tp + p1 / 2.0d0
	tr = tr + r * p2 / 2.0d0

c experiment with local currents

 	a1 = g(kx1 - 1, ky - 1)
	b1 = g(kx1 - 1, ky)
	c1 = g(kx1, ky - 1)
	d1 = g(kx1, ky)

 	a2 = g(kx2 - 1, ky - 1)
	b2 = g(kx2 - 1, ky)
	c2 = g(kx2, ky - 1)
	d2 = g(kx2, ky)

 	a3 = g(kx3 - 1, ky - 1)
	b3 = g(kx3 - 1, ky)
	c3 = g(kx3, ky - 1)
	d3 = g(kx3, ky)

 	a4 = g(kx4 - 1, ky - 1)
	b4 = g(kx4 - 1, ky)
	c4 = g(kx4, ky - 1)
	d4 = g(kx4, ky)

c discretization: (f(x) - f(x - a)) / a
	spa = spa + (dconjg(a) * c + dconjg(b) * d +
     &         a * dconjg(c) + b * dconjg(d) +
     &       dcmplx(0.0d0, 1.0d0) * (
     &         a * dconjg(c) - dconjg(a) * c +
     &         b * dconjg(d) - dconjg(b) * d) -
     &       dconjg(a3) * c3 - dconjg(b3) * d3 -
     &         a3 * dconjg(c3) - b3 * dconjg(d3) -
     &       dcmplx(0.0d0, 1.0d0) * (
     &         a4 * dconjg(c4) - dconjg(a4) * c4 +
     &         b4 * dconjg(d4) - dconjg(b4) * d4)) / 2.0d0

c discretization: (f(x + a) - f(x)) / a
 	spb = spb + (dconjg(a1) * c1 + dconjg(b1) * d1 +
     &         a1 * dconjg(c1) + b1 * dconjg(d1) +
     &       dcmplx(0.0d0, 1.0d0) * (
     &         a2 * dconjg(c2) - dconjg(a2) * c2 +
     &         b2 * dconjg(d2) - dconjg(b2) * d2) -
     &       dconjg(a) * c - dconjg(b) * d -
     &         a * dconjg(c) - b * dconjg(d) -
     &       dcmplx(0.0d0, 1.0d0) * (
     &         a * dconjg(c) - dconjg(a) * c +
     &         b * dconjg(d) - dconjg(b) * d)) / 2.0d0

c discretization: (f(x + a) - f(x - a)) / (2 a)
 	spc = spc + (dconjg(a1) * c1 + dconjg(b1) * d1 +
     &         a1 * dconjg(c1) + b1 * dconjg(d1) +
     &       dcmplx(0.0d0, 1.0d0) * (
     &         a2 * dconjg(c2) - dconjg(a2) * c2 +
     &         b2 * dconjg(d2) - dconjg(b2) * d2) -
     &       dconjg(a3) * c3 - dconjg(b3) * d3 -
     &         a3 * dconjg(c3) - b3 * dconjg(d3) -
     &       dcmplx(0.0d0, 1.0d0) * (
     &         a4 * dconjg(c4) - dconjg(a4) * c4 +
     &         b4 * dconjg(d4) - dconjg(b4) * d4)) / 4.0d0

c discretization: (f(x) - f(x - a)) / a
c        s1 = dconjg(a) * c + dconjg(b) * d +
c     &         a * dconjg(c) + b * dconjg(d)
c        s2 = dcmplx(0.0d0, 1.0d0) * (
c     &         a * dconjg(c) - dconjg(a) * c +
c     &         b * dconjg(d) - dconjg(b) * d)
c        s3 = dconjg(a3) * c3 + dconjg(b3) * d3 +
c     &         a3 * dconjg(c3) + b3 * dconjg(d3)
c        s4 = dcmplx(0.0d0, 1.0d0) * (
c     &         a4 * dconjg(c4) - dconjg(a4) * c4 +
c     &         b4 * dconjg(d4) - dconjg(b4) * d4)
c
c        sp = sp + (s1 + s2 - s3 - s4) / 2.0d0

	sr = 0.0d0

            end do
          end do
        end do

	f(it, 1)  = t0
	f(it, 2)  = t1
	f(it, 3)  = t2
	f(it, 4)  = t3

	f(it, 5)  = v0
	f(it, 6)  = v1
	f(it, 7)  = v2
	f(it, 8)  = v3

	f(it, 9)  = cx0
	f(it, 10) = cx1
	f(it, 11) = cx2
	f(it, 12) = cx3

	f(it, 13) = cy0
	f(it, 14) = cy1
	f(it, 15) = cy2
	f(it, 16) = cy3

c	write(*, '(1x, i3, 4e12.3)') it, s11, s12, -s13, -s14
c	write(*, '(1x, i3, 4e12.3)') it, p11, p12, -p13, -p14
c	f(it, 17) = sp
c	f(it, 18) = sr
c	f(it, 19) = tp
c	f(it, 20) = tr

	f(it, 17) = spa
	f(it, 18) = spb
	f(it, 19) = tp
	f(it, 20) = spc

c	write(*, *) t3 / (tp + tr)

      end do
c	write(*, *)

C everything was summed once over all time and twice over all space;
C that should be normalized!
      do it = 1, NTIME - 1
        do ic = 1, mcomp
          f4(it, ic) = sngl(f(it, ic) / dble(NTIME * NSPACE**2))
        end do
      end do

      return
      end


C     ------------------------------------------------------------------
      subroutine imass_pp(f4)
C     ------------------------------------------------------------------
C     v3 - ivh - 06 Oct 97 - Last modified: 06 Oct 97
C     ------------------------------------------------------------------
C     Computes current correlations and PCAC
C     - POINT TO POINT propagation
C     - summation over all time slices - because of symmetry, only 
C       NTIME / 2 (instead of NTIME - 1) entries should be keeped;
C       however, PCAC isn't symmetric (?), and we also want to be
C       compatible with old (or any other unsymmetric) measurements
C     ------------------------------------------------------------------
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
      parameter (mcomp=20)
C     ------------------------------------------------------------------
C INOUT real*4 array f4:
      real*4 f4(NTIME - 1, mcomp)

C GLOBAL gauge field (in common /gauge_fields/):
      complex*16 u(nsite, ndim)

C GLOBAL inverted matrix - result of inv (in common /invers/)
      complex*16 g(nferm, nferm)

C GLOBAL index array (in common /index_arrays/)
      integer ind(4, nsite)

C double precision equivalent of INOUT parameter f
      real*8 f(NTIME - 1, mcomp)

C t-s are triplet correlations Tr(sigma_i G_{yx} sigma_i G_{yx}):
      real*8 t0, t1, t2, t3

C v-s are vacuum fluctuations Tr(sigma_i G_{xx}) Tr(sigma_i G_{yy}):
      real*8 v0, v1, v2, v3

C x-s are Tr(sigma_i G_{xx}):
      real*8 x0, x1, x2, x3, cx0, cx1, cx2, cx3

C y-s are Tr(sigma_i G_{yy}):
      real*8 y0, y1, y2, y3, cy0, cy1, cy2, cy3

C PCAC (sp - sigma1 PCAC, sr - sigma1 r-corrections, t - triplet):
      real*8 sp, sr, tp, tr

C          | a   b |
C G_{xy} = |       |
C          | c   d |
      complex*16 a, b, c, d

C needed for triplet correlations computations:
      real*8 absdg, absnd, mixdg, mixnd

C needed for PCAC computations
      real*8 s1, s2, p1, p2, r
      real*8 s11, s12, s13, s14, p11, p12, p13, p14

c needed for PCAC with local currents
      complex*16 a3, b3, c3, d3, a4, b4, c4, d4

C     ------------------------------------------------------------------
      common/gauge_fields/ u
      common/index_arrays/ ind
      common/invers/ g
C     ------------------------------------------------------------------

C Wilson parameter r is set to 1
      r = 1.0d0

C loop through all (time) distances:
      do it = 1, NTIME - 1

C initialization with zero:
        t0 = 0.0d0
        t1 = 0.0d0
        t2 = 0.0d0
        t3 = 0.0d0
        v0 = 0.0d0
        v1 = 0.0d0
        v2 = 0.0d0
        v3 = 0.0d0
        cx0 = 0.0d0
        cx1 = 0.0d0
        cx2 = 0.0d0
        cx3 = 0.0d0
        cy0 = 0.0d0
        cy1 = 0.0d0
        cy2 = 0.0d0
        cy3 = 0.0d0
        sp = 0.0d0
        sr = 0.0d0
        tp = 0.0d0
        tr = 0.0d0

	s11 = 0.0d0
	s12 = 0.0d0
	s13 = 0.0d0
	s14 = 0.0d0
	p11 = 0.0d0
	p12 = 0.0d0
	p13 = 0.0d0
	p14 = 0.0d0

C loop through all slices
      do it0 = 0, NTIME - 1

C   loop through the "source" and "sink" time slice:
        do ix = 1, NSPACE

        ixt0 = it0 * NSPACE + ix
        iyt0 = mod(it0 + it, NTIME) * NSPACE + ix
 
        kx = ndrc * ixt0
        ky = ndrc * iyt0

	ix1 = ind(1, ixt0)
	ix2 = ind(2, ixt0)
	ix3 = ind(3, ixt0)
	ix4 = ind(4, ixt0)

	kx1 = ndrc * ix1
	kx2 = ndrc * ix2
	kx3 = ndrc * ix3
	kx4 = ndrc * ix4

 	a = g(kx - 1, ky - 1)
	b = g(kx - 1, ky)
	c = g(kx, ky - 1)
	d = g(kx, ky)

C computation of triplet correlations (what is computed is:
C -Tr(\gamma_i G_yx \gamma_i G_xy); triplet propagator: -2 * t#
	absdg = dreal(dconjg(a) * a + dconjg(d) * d)
	absnd = dreal(dconjg(b) * b + dconjg(c) * c)
	mixdg = dreal(dconjg(a) * d + a * dconjg(d))
        mixnd = dreal(dconjg(b) * c + b * dconjg(c))
	t0 = t0 + (absdg - absnd)
	t1 = t1 + (mixdg - mixnd)
	t2 = t2 + (mixdg + mixnd)
	t3 = t3 + (absdg + absnd)

C computation of vacuum fluctuations (what is computed is:
C Tr(\gamma_i G_xx) * Tr(\gamma_i G_yy); to compute singlet propagator
C following should be used: 4 * v# - 2 * t# - 4 * cx# * cy#
C WARNING: x1 = 2 * I * Im(b_xx) and y1 = 2 * I * Im(b_yy) but instead
C   of I-s in both expressions we take x1 positive and y1 negative,
C   that means that x1 * y1 is OK (although x1 and y1 alone are not!)
C   - the same is true for x2 and y2 - 
	x0 = g(kx - 1, kx - 1) + g(kx, kx)
	x1 = 2.0d0 * dimag(g(kx - 1, kx))
	x2 = 2.0d0 * dreal(g(kx - 1, kx))
	x3 = g(kx - 1, kx - 1) - g(kx, kx)

	cx0 = cx0 + x0
	cx1 = cx1 + x1
	cx2 = cx2 + x2
	cx3 = cx3 + x3

	y0 = g(ky - 1, ky - 1) + g(ky, ky)
	y1 = - 2.0d0 * dimag(g(ky - 1, ky))
	y2 = - 2.0d0 * dreal(g(ky - 1, ky))
	y3 = g(ky - 1, ky - 1) - g(ky, ky)

        cy0 = cy0 + y0
	cy1 = cy1 + y1
	cy2 = cy2 + y2
	cy3 = cy3 + y3

	v0 = v0 + x0 * y0
	v1 = v1 + x1 * y1
	v2 = v2 + x2 * y2
	v3 = v3 + x3 * y3

C PCAC for triplet current:
C - tp is PCAC
C - tr is r-correction

c this was experiment with \sigma_1 current instead \sigma_3...
c        s11 = s11 + dreal( u(ixt0, 1) * (
c     &  - dconjg(d) * g(kx1 - 1, ky - 1) - dconjg(b) * g(kx1, ky - 1) +
c     &    dconjg(c) * g(kx1 - 1, ky) + dconjg(a) * g(kx1, ky) ) )
c        s12 = s12 + dreal( u(ixt0, 2) * dcmplx(0.0d0, 1.0d0) * (
c     &  - dconjg(d) * g(kx2 - 1, ky - 1) + dconjg(b) * g(kx2, ky - 1) +
c     &    dconjg(c) * g(kx2 - 1, ky) - dconjg(a) * g(kx2, ky) ) ) 
c        s13 = s13 + dreal( u(ix3, 1) * (
c     &    d * dconjg(g(kx3 - 1, ky - 1)) + b * dconjg(g(kx3, ky - 1)) -
c     &    c * dconjg(g(kx3 - 1, ky)) - a * dconjg(g(kx3, ky)) ) ) 
c        s14 = s14 + dreal( u(ix4, 2) * dcmplx(0.0d0, 1.0d0) * (
c     &  - d * dconjg(g(kx4 - 1, ky - 1)) + b * dconjg(g(kx4, ky - 1)) +
c     &    c * dconjg(g(kx4 - 1, ky)) - a * dconjg(g(kx4, ky)) ) ) 

	p1 =
     &	dreal( u(ixt0, 1) * (
     &    dconjg(c) * g(kx1 - 1, ky - 1) + dconjg(a) * g(kx1, ky - 1) +
     &    dconjg(d) * g(kx1 - 1, ky) + dconjg(b) * g(kx1, ky) ) ) +
     &  dreal( u(ixt0, 2) * dcmplx(0.0d0, 1.0d0) * (
     &    dconjg(c) * g(kx2 - 1, ky - 1) - dconjg(a) * g(kx2, ky - 1) +
     &    dconjg(d) * g(kx2 - 1, ky) - dconjg(b) * g(kx2, ky) ) ) -
     &  dreal( u(ix3, 1) * (
     &    c * dconjg(g(kx3 - 1, ky - 1)) + a * dconjg(g(kx3, ky - 1)) +
     &    d * dconjg(g(kx3 - 1, ky)) + b * dconjg(g(kx3, ky)) ) ) -
     &  dreal( u(ix4, 2) * dcmplx(0.0d0, 1.0d0) * (
     &    a * dconjg(g(kx4, ky - 1)) - c * dconjg(g(kx4 - 1, ky - 1)) +
     &    b * dconjg(g(kx4, ky)) - d * dconjg(g(kx4 - 1, ky)) ) )

c unclear what was the purpose of this part of the code...
c      p11 = p11 + dreal( u(ixt0, 1) * (
c     &    dconjg(c) * g(kx1 - 1, ky - 1) + dconjg(a) * g(kx1, ky - 1) +
c     &    dconjg(d) * g(kx1 - 1, ky) + dconjg(b) * g(kx1, ky) ) ) 
c      p12 = p12 + dreal( u(ixt0, 2) * dcmplx(0.0d0, 1.0d0) * (
c     &    dconjg(c) * g(kx2 - 1, ky - 1) - dconjg(a) * g(kx2, ky - 1) +
c     &    dconjg(d) * g(kx2 - 1, ky) - dconjg(b) * g(kx2, ky) ) ) 
c      p13 = p13 + dreal( u(ix3, 1) * (
c     &    c * dconjg(g(kx3 - 1, ky - 1)) + a * dconjg(g(kx3, ky - 1)) +
c     &    d * dconjg(g(kx3 - 1, ky)) + b * dconjg(g(kx3, ky)) ) ) 
c      p14 = p14 + dreal( u(ix4, 2) * dcmplx(0.0d0, 1.0d0) * (
c     &    a * dconjg(g(kx4, ky - 1)) - c * dconjg(g(kx4 - 1, ky - 1)) +
c     &    b * dconjg(g(kx4, ky)) - d * dconjg(g(kx4 - 1, ky)) ) )

	p2 =
     &  dreal( u(ixt0, 1) * (
     &    dconjg(a) * g(kx1 - 1, ky - 1) + dconjg(c) * g(kx1, ky - 1) +
     &    dconjg(b) * g(kx1 - 1, ky) + dconjg(d) * g(kx1, ky) ) ) +
     &  dreal( u(ixt0, 2) * (
     &    dconjg(a) * g(kx2 - 1, ky - 1) + dconjg(c) * g(kx2, ky - 1) +
     &    dconjg(b) * g(kx2 - 1, ky) + dconjg(d) * g(kx2, ky) ) ) +
     &  dreal( u(ix3, 1) * (
     &    a * dconjg(g(kx3 - 1, ky - 1)) + c * dconjg(g(kx3, ky - 1)) +
     &    b * dconjg(g(kx3 - 1, ky)) + d * dconjg(g(kx3, ky)) ) ) +
     &  dreal( u(ix4, 2) * (
     &    a * dconjg(g(kx4 - 1, ky - 1)) + c * dconjg(g(kx4, ky - 1)) +
     &    b * dconjg(g(kx4 - 1, ky)) + d * dconjg(g(kx4, ky)) ) )

c also \sigma_1 experiment...
c	s1 = s11 + s12 - s13 - s14
c	sp = sp + s1 / 2.0d0
C >>> future expansion: r-correction for sigma1 current
C >>> temporary set to zero
c	sr = 0.0d0
	
	tp = tp + p1 / 2.0d0
	tr = tr + r * p2 / 2.0d0

c experiment with local currents

 	a3 = g(kx3 - 1, ky - 1)
	b3 = g(kx3 - 1, ky)
	c3 = g(kx3, ky - 1)
	d3 = g(kx3, ky)

 	a4 = g(kx4 - 1, ky - 1)
	b4 = g(kx4 - 1, ky)
	c4 = g(kx4, ky - 1)
	d4 = g(kx4, ky)

	sp = sp + dconjg(a) * c + dconjg(b) * d +
     &         a * dconjg(c) + b * dconjg(d) +
     &       dcmplx(0.0d0, 1.0d0) * (
     &         a * dconjg(c) - dconjg(a) * c +
     &         b * dconjg(d) - dconjg(b) * d) -
     &       dconjg(a3) * c3 - dconjg(b3) * d3 -
     &         a3 * dconjg(c3) - b3 * dconjg(d3) -
     &       dcmplx(0.0d0, 1.0d0) * (
     &         a4 * dconjg(c4) - dconjg(a4) * c4 +
     &         b4 * dconjg(d4) - dconjg(b4) * d4)

	sr = 0.0d0

          end do
        end do

	f(it, 1)  = t0
	f(it, 2)  = t1
	f(it, 3)  = t2
	f(it, 4)  = t3

	f(it, 5)  = v0
	f(it, 6)  = v1
	f(it, 7)  = v2
	f(it, 8)  = v3

	f(it, 9)  = cx0
	f(it, 10) = cx1
	f(it, 11) = cx2
	f(it, 12) = cx3

	f(it, 13) = cy0
	f(it, 14) = cy1
	f(it, 15) = cy2
	f(it, 16) = cy3

c	write(*, '(1x, i3, 4e12.3)') it, s11, s12, -s13, -s14
c	write(*, '(1x, i3, 4e12.3)') it, p11, p12, -p13, -p14
	f(it, 17) = sp
	f(it, 18) = sr
	f(it, 19) = tp
	f(it, 20) = tr

c	write(*, *) t3 / (tp + tr)

      end do
c	write(*, *)

C everything was summed once over all time and once (POINT TO
C POINT!) over all space; that should be normalized!
      do it = 1, NTIME - 1
        do ic = 1, mcomp
          f4(it, ic) = sngl(f(it, ic) / dble(NTIME * NSPACE))
        end do
      end do

      return
      end


