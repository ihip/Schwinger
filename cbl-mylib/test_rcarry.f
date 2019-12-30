c     f77 test_rcarry.f    
c
c     50K: 0.0152 each call (asterix)
c
      parameter(nr=50000)
      real a(nr),b(nr)
      call cbl_rc_seed(1526565)
      call cbl_cputime(0,t1)
      do i=1,10
        call cbl_rgauss(a,nr,1.001)
      enddo
      call cbl_cputime(1,tdiff)
      sum=0
      do i=1,nr
        sum=sum+a(i)
      enddo
      print *,sum
      print *,'cbl_rgauss (50K) :',tdiff/(10.*nr)
      call cbl_cputime(0,t1)
      do i=1,10
        call cbl_rgauss1(a,nr,1.001,b)
      enddo
      call cbl_cputime(1,tdiff)
      sum=0
      do i=1,nr
        sum=sum+a(i)
      enddo
      print *,sum
      print *,'cbl_rgauss1(50K) :',tdiff/(10.*nr)
      end

C     ******************************************************************
      subroutine cbl_rgauss1(x,n,sigma,ws)
C     ******************************************************************
C     ------------------------------------------------------------------
C     Part of package:        ctu1
C     ------------------------------------------------------------------
C     Subroutine for: GAUSSIAN RANDOM NUMBER GENERATOR
C                     a la Box-Muller
C     ------------------------------------------------------------------
C     Input variables:
C         sigma   :(related to width of distribution)
C         n       :length of random vector x 
C     ------------------------------------------------------------------
C     Output variables: 
C         x(1..n) :vector of random numbers with gaussian 
C                  distribution: exp( -x**2/sigma ) 
C     ------------------------------------------------------------------
C     Remarks:
C               Following Press (p 203) one could avoid the sin and cos,
C               call, but has to call more random numbers and check
C               them and eventually call more. This does not vectorize.
C               Eventually one should do this in a part cbl_rgauss
C     ------------------------------------------------------------------
C     CPU-Requirements (in 10**(-6) sec):  
C     vectorlength       MIPS4000
C     50000              2.42
C     gibt falsche Ergebnisse??? 
C     ******************************************************************
C     ------------------------------------------------------------------
      real x(n),xx(2),ws(n)
C     ------------------------------------------------------------------
      pi2=atan(1.0)*8.0
      n2=n/2

      call cbl_rcarry(ws,n,0.)
      do i=1,n
        ws(i)=2.*ws(i)-1.
      enddo

      do i=1,n2
        r2=ws(i)**2+ws(n2+i)**2
        if((r2.lt.1.0).and.(r2.ne.0.))then
          a=sqrt(-sigma*alog(r2)/r2)
          x(i)     = a*ws(i)
          x(n2+i)  = a*ws(n2+i)
        else
          a=sqrt(-sigma*alog(0.5*ws(n2+i)+0.5))
          x(i)     = a*cos(pi2*ws(i))
          x(n2+i)  = a*sin(pi2*ws(i))
        endif
      enddo

      if(n2*2.ne.n)then
        print *,'padding'
        call cbl_rcarry(xx,2,0.0)
        x(n)=sqrt(-sigma*alog(xx(1)))*cos(pi2*xx(2))
      endif

      end
C     ******************************************************************
      subroutine cbl_rgauss(x,n,sigma)
C     ******************************************************************
C     ------------------------------------------------------------------
C     Part of package:        ctu1
C     ------------------------------------------------------------------
C     Subroutine for: GAUSSIAN RANDOM NUMBER GENERATOR
C                     a la Box-Muller
C     ------------------------------------------------------------------
C     Input variables:
C         sigma   :(related to width of distribution)
C         n       :length of random vector x and workspace ws
C     ------------------------------------------------------------------
C     Output variables: 
C         x(1..n) :vector of random numbers with gaussian 
C                  distribution: exp( -x**2/sigma ) 
C     ------------------------------------------------------------------
C     Remarks:
C               Following Press (p 203) one couls avoid the sin and cos,
C               call, but has to call more random numbers and check
C               them and eventually call more. This does not vectorize.
C               Eventually one should do this in a part cbl_rgauss
C     ------------------------------------------------------------------
C     CPU-Requirements (in 10**(-6) sec):  
C     vectorlength       MIPS4000
C     50000              3.4 
C     ------------------------------------------------------------------
C     ******************************************************************
C     ------------------------------------------------------------------
      real x(n),xx(2)
C     ------------------------------------------------------------------
      pi2=atan(1.0)*8.0
      n2=n/2

      call cbl_rc_log(x,n2,0.0)
      call cbl_rcarry(x(n2+1),n2,0.0)

      do i=1,n2
        a        = sqrt(-sigma*x(i))
        b        = pi2*x(n2+i)
        x(i)     = a*cos(b)
        x(n2+i)  = a*sin(b)
      enddo

      if(n2*2.ne.n)then
        call cbl_rcarry(xx,2,0.0)
        x(n)=sqrt(-sigma*alog(xx(1)))*cos(pi2*xx(2))
      endif

      end
