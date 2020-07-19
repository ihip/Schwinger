C package inverteo consists of:
C - invert(x, b) - invert routine compatible with original cgr(x, b)
C   inversion of (Q^+ Q) matrix
C - mk_index_eo(indeo) - make index for even-odd preconditioned matrix
C - mkeon(n, eon) & mkneo(n, neo) - make indices for transformations 
C   between normal and vectors which are compatible with even-odd 
C   preconditioning
C - mdeo(a, b), mtdeo(a, b), mdoe(a, b) & mtdoe(a, b) multiplication
C   routines compatible with even odd preconditioning
C
C two additional subroutines should be provided:
C - invert1eo(x, b) & invert1teo(x, b) - routines which solve separate
C   inversions for Q and Q^+ matrices

C     ******************************************************************
      function invert(x, b)
C     ******************************************************************
C     **                     **                                       **
C     **      INVERTEO       **     Program by I. Hip, 05 Sep 95      **
C     **                     **     Last modified: 13 Feb 97          **
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
      complex*16 beo(ndrc, nsite), xeo(ndrc, nsite), tmp(ndrc, nsite)
      complex*16 u(nsite, ndim)
C common /eo/: (razmisliti o indeksiranju ueo polja)      
      complex*16 ueo(nsite, ndim)
      integer indeo(4, nsite), neo(nsite), eon(nsite)
C     ------------------------------------------------------------------
      common /gauge_fields/u
      common /eo/ indeo, neo, eon, ueo
C     ------------------------------------------------------------------
      
      do i = 1, nsite
	ueo(i, 1) = u(eon(i), 1)
	ueo(i, 2) = u(eon(i), 2)
	beo(1, i) = b(1, eon(i))
	beo(2, i) = b(2, eon(i))
      end do
      
      n1t = invert1teo(tmp, beo)
      if(n1t .lt. 0) then
	invert = n1t
	return
      end if
 
      n1 = invert1eo(xeo, tmp)
      if(n1 .lt. 0) then
	invert = n1
	return
      end if

      do i = 1, nsite
	x(1, i) = xeo(1, neo(i))
	x(2, i) = xeo(2, neo(i))
      end do

      invert = n1t + n1
      return
      end

C     ******************************************************************
      subroutine mk_index_eo(ind)
C     ******************************************************************
C     **                    **                                        **
C     **  MK_INDEX_EO       **    Program by Ivan Hip,  19-09-1995    **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for construction of d=2 index vectors
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Output variables:   ind(1:4,1:nsite)
C     ------------------------------------------------------------------
C     Remarks: index (1:nsite) = (1:NTIME, 1:NSPACE)
C                            i = (it-1)*NTIME+ix
C     ------------------------------------------------------------------
C     WARNING: NSPACE and NTIME should be equal and even !!!
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      integer ind(4,nsite)
C     ------------------------------------------------------------------
      integer h, s, i, n

      if(NSPACE .ne. NTIME) stop 'NSPACE and NTIME are not equal'
      if(mod(NSPACE, 2) .ne. 0) stop 'NSPACE and NTIME are not even'

      n = NSPACE
      h = n / 2
      s = n**2 / 2
      do i = 1, s
	if(mod((i - 1) / h, 2) .eq. 1) then
	  if(mod(i, n) .eq. 0) then
	    ind(1, i) = s + i + 1 - h
	  else
	    ind(1, i) = s + i + 1
	  end if
	  if(i .gt. (s - h)) then
	    ind(2, i) = i + h
	    ind(2, s + i) = i + h - s
	  else
	    ind(2, i) = s + i + h
	    ind(2, s + i) = i + h
	  end if
	  ind(3, i) = s + i
	  ind(4, i) = s + i - h
	  ind(1, s + i) = i
	  if(mod(i, n) .eq. (h + 1)) then
	    ind(3, s + i) = i - 1 + h
	  else
	    ind(3, s + i) = i - 1
	  end if
	  ind(4, s + i) = i - h
	else
	  ind(1, i) = s + i
	  ind(2, i) = s + i + h
	  if(mod(i, n) .eq. 1) then
	    ind(3, i) = s + i - 1 + h
	  else
	    ind(3, i) = s + i - 1
	  end if
	  if(i .le. h) then
	    ind(4, i) = 2 * s + i - h
	    ind(4, s + i) = s + i - h
	  else
	    ind(4, i) = s + i - h
	    ind(4, s + i) = i - h
	  end if
	  if(mod(i, n) .eq. h) then
	    ind(1, s + i) = i + 1 - h
	  else
	    ind(1, s + i) = i + 1
	  end if
	  ind(2, s + i) = i + h
	  ind(3, s + i) = i
	end if
      end do
      end

	subroutine mkneo(n, neo)
	integer n, neo(*)

	integer h, s, i

	h = n / 2
	s = n**2 / 2
	do i = 1, n**2
	  if(mod((i - 1) / n, 2) .eq. 0) then
	    if(mod(i, 2) .eq. 1) then
	      neo(i) = (i + 1) / 2
	    else
	      neo(i) = s + i / 2
	    end if
	  else
	    if(mod(i, 2) .eq. 1) then
	      neo(i) = s + (i + 1) / 2
	    else
	      neo(i) = i / 2
	    end if
	  end if
	end do
	end

	subroutine mkeon(n, eon)
	integer n, eon(*)

	integer h, s, i

	h = n / 2
	s = n**2 / 2
	do i = 1, s
	  if(mod((i - 1) / h, 2) .eq. 0) then
	    eon(i) = 2 * i - 1
	    eon(s + i) = 2 * i
	  else
	    eon(i) = 2 * i
	    eon(s + i) = 2 * i - 1
	  end if
	end do
	end

C     ******************************************************************
      subroutine mdeo(a,b)
C     ******************************************************************
C     **                    **                                        **
C     **  MDEO              **    Program by Ivan Hip, 21 Sep 1995    **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for multiplication b := Q a
C     d=2
C     gauge group: U(1)
C     Wilson fermions
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input variables:   vector a
C     ------------------------------------------------------------------
C     Output variables:  vector b
C     ------------------------------------------------------------------
C     Remarks: index arrays have to be initialized (call mk_index(ind))
C     ------------------------------------------------------------------
C     index = (1:nall) = (1:ndrc,1:ngrp,1:NTIME,1:NSPACE)
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nh=nsite/2)
C     ------------------------------------------------------------------
      real*8 a(2,ndrc,nh), b(2,ndrc,nh)
      real*8 akap
      real*8 xr, xi, c1r, c1i, c2r, c2i, c3r, c3i, c4r, c4i
C common /eo/: (razmisliti o indeksiranju ueo polja)      
      real*8 ueo(2, nsite, ndim)
      integer indeo(4, nsite), neo(nsite), eon(nsite)
C     ------------------------------------------------------------------
      common /kappa/ akap
      common /eo/ indeo, neo, eon, ueo
C     ------------------------------------------------------------------
      do i = 1, nh

	i1=indeo(1,i)
	i1a = i1 - nh
	xr= a(1,1,i1a) + a(1,2,i1a)
	xi= a(2,1,i1a) + a(2,2,i1a)
	c1r =  ueo(1,i,1)*xr - ueo(2,i,1)*xi
	c1i =  ueo(1,i,1)*xi + ueo(2,i,1)*xr

	i2=indeo(2,i)
	i2a = i2 - nh
	xr= a(1,1,i2a) + a(2,2,i2a)
	xi= a(2,1,i2a) - a(1,2,i2a)
	c2r =  ueo(1,i,2)*xr - ueo(2,i,2)*xi
	c2i =  ueo(1,i,2)*xi + ueo(2,i,2)*xr

	i3=indeo(3,i)
	i3a = i3 - nh
	xr= a(1,1,i3a) - a(1,2,i3a)
	xi= a(2,1,i3a) - a(2,2,i3a)
	c3r=  ueo(1,i3,1)*xr + ueo(2,i3,1)*xi
	c3i=  ueo(1,i3,1)*xi - ueo(2,i3,1)*xr

	i4=indeo(4,i)
	i4a = i4 - nh
	xr= a(1,1,i4a) - a(2,2,i4a)
	xi= a(2,1,i4a) + a(1,2,i4a)
	c4r=  ueo(1,i4,2)*xr + ueo(2,i4,2)*xi
	c4i=  ueo(1,i4,2)*xi - ueo(2,i4,2)*xr

	b(1,1,i) = akap * ( c1r + c2r + c3r + c4r)
	b(2,1,i) = akap * ( c1i + c2i + c3i + c4i)
	b(1,2,i) = akap * ( c1r - c3r - c2i + c4i)
	b(2,2,i) = akap * ( c1i - c3i + c2r - c4r)

      enddo
      end

C     ******************************************************************
      subroutine mtdeo(a,b)
C     ******************************************************************
C     **                    **                                        **
C     **  MTDEO             **    Program by Ivan Hip, 21 Sep 1995    **
C     **                    **                                        **
C     ******************************************************************

C     ------------------------------------------------------------------
C     Program for multiplication b := (Q^+) a
C     d=2
C     gauge group: U(1)
C     Wilson fermions
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input variables:   vector a
C     ------------------------------------------------------------------
C     Output variables:  vector b
C     ------------------------------------------------------------------

C     Remarks: index arrays have to be initialized (call mk_index(ind))
C     ------------------------------------------------------------------
C     index = (1:nall) = (1:ndrc,1:ngrp,1:NTIME,1:NSPACE)
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nh=nsite/2)
C     ------------------------------------------------------------------
      real*8 a(2,ndrc,nh), b(2,ndrc,nh)
      real*8 akap
      real*8 xr, xi, c1r, c1i, c2r, c2i, c3r, c3i, c4r, c4i
C common /eo/: (razmisliti o indeksiranju ueo polja)      
      real*8 ueo(2, nsite, ndim)
      integer indeo(4, nsite), neo(nsite), eon(nsite)
C     ------------------------------------------------------------------
      common /kappa/akap
      common /eo/ indeo, neo, eon, ueo
C     ------------------------------------------------------------------
      do i = 1, nh

	i1=indeo(1,i)
	i1a = i1 - nh
	xr= a(1,1,i1a) - a(1,2,i1a)
	xi= a(2,1,i1a) - a(2,2,i1a)
	c1r =  ueo(1,i,1)*xr - ueo(2,i,1)*xi
	c1i =  ueo(1,i,1)*xi + ueo(2,i,1)*xr

	i2=indeo(2,i)
	i2a = i2 - nh
	xr= a(1,1,i2a) - a(2,2,i2a)
	xi= a(2,1,i2a) + a(1,2,i2a)
	c2r =  ueo(1,i,2)*xr - ueo(2,i,2)*xi
	c2i =  ueo(1,i,2)*xi + ueo(2,i,2)*xr

	i3 =indeo(3,i)
	i3a = i3 - nh
	xr= a(1,1,i3a) + a(1,2,i3a)
	xi= a(2,1,i3a) + a(2,2,i3a)
	c3r=  ueo(1,i3,1)*xr + ueo(2,i3,1)*xi
	c3i=  ueo(1,i3,1)*xi - ueo(2,i3,1)*xr

	i4 =indeo(4,i)
	i4a = i4 - nh
	xr= a(1,1,i4a) + a(2,2,i4a)
	xi= a(2,1,i4a) - a(1,2,i4a)
	c4r=  ueo(1,i4,2)*xr + ueo(2,i4,2)*xi
	c4i=  ueo(1,i4,2)*xi - ueo(2,i4,2)*xr

	b(1,1,i) = akap * ( c1r + c2r + c3r + c4r)
	b(2,1,i) = akap * ( c1i + c2i + c3i + c4i)
	b(1,2,i) = akap * (-c1r + c3r + c2i - c4i)
	b(2,2,i) = akap * (-c1i + c3i - c2r + c4r)

      enddo
      end

C     ******************************************************************
      subroutine mdoe(a,b)
C     ******************************************************************
C     **                    **                                        **
C     **  MDOE              **    Program by Ivan Hip, 21 Sep 1995    **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for multiplication b := Q a
C     d=2
C     gauge group: U(1)
C     Wilson fermions
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input variables:   vector a
C     ------------------------------------------------------------------
C     Output variables:  vector b
C     ------------------------------------------------------------------
C     Remarks: index arrays have to be initialized (call mk_index(ind))
C     ------------------------------------------------------------------
C     index = (1:nall) = (1:ndrc,1:ngrp,1:NTIME,1:NSPACE)
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nh=nsite/2)
C     ------------------------------------------------------------------
      real*8 a(2,ndrc,nh), b(2,ndrc,nh)
      real*8 akap
      real*8 xr, xi, c1r, c1i, c2r, c2i, c3r, c3i, c4r, c4i
C common /eo/: (razmisliti o indeksiranju ueo polja)      
      real*8 ueo(2, nsite, ndim)
      integer indeo(4, nsite), neo(nsite), eon(nsite)
C     ------------------------------------------------------------------
      common /kappa/akap
      common /eo/ indeo, neo, eon, ueo
C     ------------------------------------------------------------------
      do k = 1, nsite / 2
	i = k + nsite / 2

	i1=indeo(1,i)
	xr= a(1,1,i1) + a(1,2,i1)
	xi= a(2,1,i1) + a(2,2,i1)
	c1r =  ueo(1,i,1)*xr - ueo(2,i,1)*xi
	c1i =  ueo(1,i,1)*xi + ueo(2,i,1)*xr

	i2=indeo(2,i)
	xr= a(1,1,i2) + a(2,2,i2)
	xi= a(2,1,i2) - a(1,2,i2)
	c2r =  ueo(1,i,2)*xr - ueo(2,i,2)*xi
	c2i =  ueo(1,i,2)*xi + ueo(2,i,2)*xr

	i3 =indeo(3,i)
	xr= a(1,1,i3) - a(1,2,i3)
	xi= a(2,1,i3) - a(2,2,i3)
	c3r=  ueo(1,i3,1)*xr + ueo(2,i3,1)*xi
	c3i=  ueo(1,i3,1)*xi - ueo(2,i3,1)*xr

	i4=indeo(4,i)
	xr= a(1,1,i4) - a(2,2,i4)
	xi= a(2,1,i4) + a(1,2,i4)
	c4r=  ueo(1,i4,2)*xr + ueo(2,i4,2)*xi
	c4i=  ueo(1,i4,2)*xi - ueo(2,i4,2)*xr

	b(1,1,k) = akap * ( c1r + c2r + c3r + c4r)
	b(2,1,k) = akap * ( c1i + c2i + c3i + c4i)
	b(1,2,k) = akap * ( c1r - c3r - c2i + c4i)
	b(2,2,k) = akap * ( c1i - c3i + c2r - c4r)

      enddo
      end

C     ******************************************************************
      subroutine mtdoe(a,b)
C     ******************************************************************
C     **                    **                                        **
C     **  MTDOE             **    Program by Ivan Hip, 21 Sep 1995    **
C     **                    **                                        **
C     ******************************************************************

C     ------------------------------------------------------------------
C     Program for multiplication b := (Q^+) a
C     d=2
C     gauge group: U(1)
C     Wilson fermions
C     double precision version
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Input variables:   vector a
C     ------------------------------------------------------------------
C     Output variables:  vector b
C     ------------------------------------------------------------------

C     Remarks: index arrays have to be initialized (call mk_index(ind))
C     ------------------------------------------------------------------
C     index = (1:nall) = (1:ndrc,1:ngrp,1:NTIME,1:NSPACE)
C     ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nh=nsite/2)
C     ------------------------------------------------------------------
      real*8 a(2,ndrc,nh), b(2,ndrc,nh)
      real*8 akap
      real*8 xr, xi, c1r, c1i, c2r, c2i, c3r, c3i, c4r, c4i
C common /eo/: (razmisliti o indeksiranju ueo polja)      
      real*8 ueo(2, nsite, ndim)
      integer indeo(4, nsite), neo(nsite), eon(nsite)
C     ------------------------------------------------------------------
      common /kappa/akap
      common /eo/ indeo, neo, eon, ueo
C     ------------------------------------------------------------------
      do k = 1, nsite / 2
	i = k + nsite / 2

	i1=indeo(1,i)
	xr= a(1,1,i1) - a(1,2,i1)
	xi= a(2,1,i1) - a(2,2,i1)
	c1r =  ueo(1,i,1)*xr - ueo(2,i,1)*xi
	c1i =  ueo(1,i,1)*xi + ueo(2,i,1)*xr

	i2=indeo(2,i)
	xr= a(1,1,i2) - a(2,2,i2)
	xi= a(2,1,i2) + a(1,2,i2)
	c2r =  ueo(1,i,2)*xr - ueo(2,i,2)*xi
	c2i =  ueo(1,i,2)*xi + ueo(2,i,2)*xr

	i3 =indeo(3,i)
	xr= a(1,1,i3) + a(1,2,i3)
	xi= a(2,1,i3) + a(2,2,i3)
	c3r=  ueo(1,i3,1)*xr + ueo(2,i3,1)*xi
	c3i=  ueo(1,i3,1)*xi - ueo(2,i3,1)*xr

	i4=indeo(4,i)
	xr= a(1,1,i4) + a(2,2,i4)
	xi= a(2,1,i4) - a(1,2,i4)
	c4r=  ueo(1,i4,2)*xr + ueo(2,i4,2)*xi
	c4i=  ueo(1,i4,2)*xi - ueo(2,i4,2)*xr

	b(1,1,k) = akap * ( c1r + c2r + c3r + c4r)
	b(2,1,k) = akap * ( c1i + c2i + c3i + c4i)
	b(1,2,k) = akap * (-c1r + c3r + c2i - c4i)
	b(2,2,k) = akap * (-c1i + c3i - c2r + c4r)

      enddo
      end


C     >>> dummy definitions

C      function invert1eo(xx, bb)
C      real*8 xx(*), bb(*)
C      write(*, *) 'invert1eo'
C      invert1eo = 0
C      end

C      function invert1teo(xx, bb)
C      real*8 xx(*), bb(*)
C      write(*, *) 'invert1teo'
C      invert1teo = 0
C      end
