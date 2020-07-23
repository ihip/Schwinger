C     ******************************************************************
      function invert1(x, b)
C     ******************************************************************
C     **                     **                                       **
C     **     INVERT1EO       **     Program by I. Hip, 13 Feb 97      **
C     **                     **     Last modified: 21 Feb 97          **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program computes x which satisfies Q^+ Q x = b
C     even-odd preconditioning (ASSUMING: NSPACE = NTIME)
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input var: vector b
C     ------------------------------------------------------------------
C     Output variables: vector x
C     ------------------------------------------------------------------
C     Remarks:
C       - the following index arrays have to be initialized:
C           indeo (call mk_index_eo(indeo))
C           neo (call mkneo(NSPACE, neo))
C           eon (call mkeon(NSPACE, eon))
C       - common block /eo/ should be present
C       - subroutines INVERT1EO, INVERT1TEO, MDEO, MTDEO, MDOE and
C         MTDOE should be available
C     -----------------------------------------------------------------
C     On exit x contains (Q^+ Q)^(-1) b   (i.e. the inverse of Q^+ Q)
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
C     ------------------------------------------------------------------
      complex*16 x(ndrc, nsite), b(ndrc, nsite)
      complex*16 beo(ndrc, nsite), xeo(ndrc, nsite)
      complex*16 u(nsite, ndim)
      complex*16 ueo(nsite, ndim)
      integer indeo(4, nsite), neo(nsite), eon(nsite)
C     ------------------------------------------------------------------
      common /gauge_fields/u
c      common /ufix/u
      common /eo/ indeo, neo, eon, ueo
C     ------------------------------------------------------------------

      do i = 1, nsite
	ueo(i, 1) = u(eon(i), 1)
	ueo(i, 2) = u(eon(i), 2)
	beo(1, i) = b(1, eon(i))
	beo(2, i) = b(2, eon(i))
      end do

      n1 = invert1eo1(xeo, beo)
      if(n1 .lt. 0) then
	invert1 = n1
	return
      end if

      do i = 1, nsite
	x(1, i) = xeo(1, neo(i))
	x(2, i) = xeo(2, neo(i))
      end do

      invert1 = n1
      return
      end

C     >>> dummy definition

C      function invert1eo1(xx, bb)
C      real*8 xx(*), bb(*)
C      write(*, *) 'invert1eo1'
C      invert1eo1 = 0
C      end

