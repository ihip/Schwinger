C slight modifications by hip, 29 Sep 98 (hopefully not essential!)

      character*32 function v_integrate()
      v_integrate = 'metropolis'
      end

C     ******************************************************************
      subroutine gauge_met_c(beta)
C     ------------------------------------------------------------------
C	IN (real*8) beta
C	IN COMMON/gauge_index_arrays/ ...  !!! should be initialized !!!
C	INOUT COMMON/gauge_fields/u
C	  complex*16 u(nspace, ntime, ndim)
C	USES cbl_rcarry, cbl_rc_log FROM mylib
C     ------------------------------------------------------------------
C	CALL gauge_mk_index before calling gauge_met_c !!!
C     ------------------------------------------------------------------
C     Part of package:  gauge_met_c
C     ------------------------------------------------------------------
C     Subroutine for:   Metropolis updating 
C                       compact U(1) gauge field
C     ------------------------------------------------------------------
C     Remarks:
C     CPU-usage on SGI R5000 180MHz for 16^2: 0.0059 sec/ update
C     ------------------------------------------------------------------
C     ******************************************************************

      parameter (ndrc=2,nspace=NSPACE,ntime=NTIME,ndim=2)
      parameter (nsite=ntime*nspace,nall=ndrc*nsite)
      parameter (nspa6=nspace+6,ntim6=ntime+6)
      parameter (nmet=3)

      complex*16 u(nspace,ntime,ndim)
      complex*16 usum,ulink,uold,unew,umul
      real*8 beta
      real ran1(nsite*nmet),ran2(nsite*nmet)
      integer iper(nspa6,ntim6),iperx(nspa6),ipert(ntim6),
     &         ixp(nspace),itp(ntime),ixm(nspace),itm(ntime)
      save aim,eps
      data aim,eps/-.300162274,0.5/ 
C                  =ALOG10(0.501)
C     ------------------------------------------------------------------
      common /gauge_index_arrays/ixp,ixm,itp,itm,iper,iperx,ipert
      common /gauge_fields/u
C     ------------------------------------------------------------------
c
c     update loop
c
      nacc=0
c
c     update direction 1
c
      call cbl_rcarry(ran1,nsite*nmet,-0.5)
      call cbl_rc_log(ran2,nsite*nmet,0.)
      nr=0
      do ix=1,nspace
        do it=1,ntime
          usum=u(ixp(ix),it,2)*
     &       dconjg(u(ix,itp(it),1))*
     &       dconjg(u(ix,it,2))
     &    +  dconjg(u(ixp(ix),itm(it),2))*
     &       dconjg(u(ix,itm(it),1))*
     &       u(ix,itm(it),2)
c
c    action = beta*Re(Ulink*usum)
c
          usum=beta*usum
          ulink=u(ix,it,1)
          uold=ulink*usum
          sold=real(uold)

          do i=1,nmet
            nr=nr+1
            umul=zexp(dcmplx(0.,eps*ran1(nr)))
            unew=uold*umul
            snew=real(unew)
            if(snew-sold.ge.ran2(nr))then
              nacc=nacc+1
              ulink=ulink*umul
              uold=unew
              sold=snew
            endif
          enddo


          u(ix,it,1)=ulink
        enddo
      enddo
c
c     update direction 2
c
      call cbl_rcarry(ran1,nsite*nmet,-0.5)
      call cbl_rc_log(ran2,nsite*nmet,0.)
      nr=0
      do ix=1,nspace
        do it=1,ntime
          usum=dconjg(u(ixm(ix),itp(it),1))*
     &       dconjg(u(ixm(ix),it,2))*
     &       u(ixm(ix),it,1)
     &    +  u(ix,itp(it),1)*
     &       dconjg(u(ixp(ix),it,2))*
     &       dconjg(u(ix,it,1))
c
c    action = beta*Re(Ulink*usum)
c
          usum=beta*usum
          ulink=u(ix,it,2)
          uold=ulink*usum
          sold=real(uold)

          do i=1,nmet
            nr=nr+1
            umul=zexp(dcmplx(0.,eps*ran1(nr)))
            unew=uold*umul
            snew=real(unew)
            if(snew-sold.ge.ran2(nr))then
              nacc=nacc+1
              ulink=ulink*umul
              uold=unew
              sold=snew
            endif
          enddo
          u(ix,it,2)=ulink
        enddo
      enddo
c
      acc=nacc/float(ntime*nspace*nmet*2)
      epsnew=eps*2.**(1.-alog10(acc)/aim) 
      eps=amax1(epsnew,0.01) 
      eps=amin1(epsnew,6.28) 
c      print *,'Acceptance = ',beta,acc,eps
      end


C     ******************************************************************
      subroutine gauge_mk_index()
C	OUT COMMON/gauge_index_arrays/ ...
C     ------------------------------------------------------------------
C     Part of package:  gauge_met_c
C     ------------------------------------------------------------------
C     ******************************************************************

      integer ndrc,nsite,nspace,ntime,ndim,nall,nspa6,ntim6,i,j

      parameter (ndrc=2,nspace=NSPACE,ntime=NTIME,ndim=2)
      parameter (nsite=ntime*nspace,nall=ndrc*nsite)
      parameter (nspa6=nspace+6,ntim6=ntime+6)

      integer iper(nspa6,ntim6),iperx(nspa6),ipert(ntim6),
     &         ixp(nspace),itp(ntime),ixm(nspace),itm(ntime)
C     ------------------------------------------------------------------
      common /gauge_index_arrays/ixp,ixm,itp,itm,iper,iperx,ipert
C     ------------------------------------------------------------------
      do i=1,nspa6
         iperx(i)=mod(i-4+nspace,nspace)+1      
      enddo
      do i=1,ntim6
         ipert(i)=mod(i-4+ntime,ntime)+1      
      enddo
      do i=1,nspa6
        do j=1,ntim6
          iper(i,j)=iperx(i)+(ipert(j)-1)*nspace
        enddo
      enddo
      do ix=1,nspace
        ixp(ix)=mod(ix,nspace)+1
        ixm(ix)=mod(ix+nspace-2,nspace)+1
      enddo
      do it=1,ntime
        itp(it)=mod(it,ntime)+1
        itm(it)=mod(it+ntime-2,ntime)+1
      enddo
      return
      end
