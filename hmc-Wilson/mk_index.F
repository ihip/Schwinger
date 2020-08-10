C     ******************************************************************
      subroutine mk_index(ind)
C     ******************************************************************
C     **                    **                                        **
C     **  MK_INDEX          **    Program by C.B. Lang, 02-03-1993    **
C     **                    **    Debugged by I. Hip, 26 Nov 96       **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Program for construction of d=2 index vectors
C     ------------------------------------------------------------------
C     Part of package:   hmc2du1
C     ------------------------------------------------------------------
C     Output variables:   ind(1:4,1:nsite)
C     ------------------------------------------------------------------
C     Remarks: index (1:nsite) = (1:NSPACE,1:NTIME)
C                            i = (it-1)*NSPACE+ix
C    ------------------------------------------------------------------
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nall=2*ndrc*ngrp*nsite)
C     ------------------------------------------------------------------
      integer ind(4,nsite)
C     ------------------------------------------------------------------
C     first the n.n. index arrays
C     ------------------------------------------------------------------
      i=0
      do it=1,NTIME
	itp=mod(it,NTIME)+1
	itm=mod(it+NTIME-2,NTIME)+1
	do ix=1,NSPACE
	  i=i+1
	  ixp=mod(ix,NSPACE)+1
	  ixm=mod(ix+NSPACE-2,NSPACE)+1
	  ind(1,i)=(it-1)*NSPACE+ixp
	  ind(2,i)=(itp-1)*NSPACE+ix
	  ind(3,i)=(it-1)*NSPACE+ixm
	  ind(4,i)=(itm-1)*NSPACE+ix
	enddo
      enddo
      end
      