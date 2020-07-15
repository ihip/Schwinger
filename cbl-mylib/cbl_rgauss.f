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
C               Following Press (p 203) one could avoid the sin and cos,
C               call, but has to call more random numbers and check
C               them and eventually call more. This does not vectorize.
C               Eventually one should do this in a part cbl_rgauss
C     ------------------------------------------------------------------
C     CPU-Requirements (in 10**(-6) sec):  
C     vectorlength       MIPS4000
C     50000              2.94 
C     ------------------------------------------------------------------
C     ******************************************************************
C     ------------------------------------------------------------------
      real x(n),xx(2)
C     ------------------------------------------------------------------
      pi2=atan(1.0)*8.0
      n2=n/2

      call cbl_rc_log(x,n2,0.0)
      call cbl_rcarry(x(n2+1),n/2,0.0)

      do i=1,n2
        a        = sqrt(-sigma*x(i))
        x(i)     = a*cos(pi2*x(n2+i))
        x(n2+i)  = a*sin(pi2*x(n2+i))
      enddo

      if(n2*2.ne.n)then  ! if n is odd last entry is computed separately
        ! print *,'padding'
        call cbl_rcarry(xx,2,0.0)
        x(n)=sqrt(-sigma*alog(xx(1)))*cos(pi2*xx(2))
      endif

      end

