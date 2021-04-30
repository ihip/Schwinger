C     ******************************************************************
C     **                    **   BASED ON:                            **
C     **  gauge_met_c_chiT  **   Program by C.B. Lang, DATE: 09/29/98 **  
C     **                    **   MODIFIED BY I. HIP / 2021-04-16      **
C     ******************************************************************
#define	NTIME	16
#define	NSPACE  16
#define KTIME	(NTIME)
#define KSPACE	(NSPACE)
C     ------------------------------------------------------------------
C     Program for testing d=2 compact U(1) gauge updates with 
C     the Metropolis method
C     ------------------------------------------------------------------
C     Input variables:
C
C     ------------------------------------------------------------------
C     Output variables:
C
C     ------------------------------------------------------------------
C     Remarks:
C     CPU-usage on SGI R500 180MHz for 16^2: 0.0059 sec/ update
C     ------------------------------------------------------------------
C     Compile with:
C     /usr/bin/f77 -O3 -mips2 -o gauge_met_c_test.x gauge_met_c_test.f \
C                  gauge_met_c.o $HOME/mylib/mylib.a
C     maybe later also:
C	 -L/usr/local/lib -lnag -lnagblas -lcomplib.slatec
C     or equivalently, using the script:
C     f77 gauge_met_c.o gauge_met_c_test.f 
C     ******************************************************************
      parameter (nspace=KSPACE,ntime=KTIME)
      parameter (nsite=nspace*ntime)
      parameter (nspa6=nspace+6,ntim6=ntime+6)
      parameter (nhist=41)
      complex*16 u(nspace,ntime,2)
      real*8 apl,top,const,plsum,beta
      integer ih(nhist)
      integer iper(nspa6,ntim6),iperx(nspa6),ipert(ntim6),
     &         ixp(nspace),itp(ntime),ixm(nspace),itm(ntime)

      real*8 qav, qsum, q2sum, q2av, chiT
C     ------------------------------------------------------------------
      common /index_fields/ixp,ixm,itp,itm,iper,iperx,ipert
      common /gauge_fields/u
C     ******************************************************************
C     ------------------------------------------------------------------
      call cbl_header(' gauge_met_c_test ') 
C     ------------------------------------------------------------------
c
c     initialization
c
C     ------------------------------------------------------------------
      npre=10000
      nmeas=100000
      call mkindex
C
      do ix=1,nspace
        do it=1,ntime
          u(ix,it,1)=dcmplx(1.,0.)
          u(ix,it,2)=dcmplx(1.,0.)
        enddo
      enddo
      print *,"nspace, ntime=", nspace,ntime
      print *,"Cold start."
C
      call cbl_rtim(iseed)
C     for testing use fixed value like e.g.
C      iseed=99727569

      print *,"cbl_rcarry initialized with iseed=",iseed
      call cbl_rc_seed(iseed)
      const=8.*datan(1.d0)
C     ------------------------------------------------------------------

55    print *,"enter next beta or 0:"
      read *,beta
      if(beta.eq.0)goto 99
c
c     pre-iterations
c
        iunit=10+int(beta+0.0001)
        call cbl_cputime(1,t1)
        do ipre=1,npre
          call gauge_met_c(beta)
        enddo
        call cbl_cputime(2,t1)
        print *, t1/npre
c
c     measurement-iterations
c
        plsum=0
c
c     initialize histogram
c
        do i=1,nhist
          ih(i)=0
        enddo
c
        do niter=1,nmeas
          call gauge_met_c(beta)
          call meas(apl,top)
          jtop=int(0.001+1000+top/const)-1000
          jt=jtop+(nhist+1)/2
          jt=max0(1,jt)
          jt=min0(nhist,jt)
          ih(jt)=ih(jt)+1
          apl=apl/nsite
          plsum=plsum+apl
          write(iunit,*)apl,top,jtop
        enddo
        close(iunit)
        plsum=plsum/nmeas
        print *,"beta = ",beta,"  plaquette = ",plsum

        qsum = 0.0d0
        q2sum = 0.0d0
        print *,"histogram"
        do i = 1, nhist
          jtop = i - (nhist + 1) / 2
          qsum = qsum + jtop * ih(i)          
          q2sum = q2sum + jtop**2 * ih(i)
          print *, jtop, ih(i)
        end do
        qav = qsum / dble(nmeas)
        q2av = q2sum / dble(nmeas)
        chiT = (q2av - qav**2) / dble(nsite)
        print *
        print *, "<Q> = ", qav
        print *, "<Q^2> = ", q2av
        write(*, '(" \chi_T         = ", f9.6)') chiT
        write(*, '(" \chi_T * \beta = ", f9.6)') chiT * beta
        print *
      goto 55
C     ------------------------------------------------------------------
99    end

