C     ******************************************************************
      complex*16 function escalp(a, b, m, length, ndrc)
C     ******************************************************************
C     **                    **                                        **
C     **  extended scalar   **  by HCH , 21-12-1995                   **
C     **  product routine   **                                        **
C     **                    **                                        **
C     ******************************************************************
C     ------------------------------------------------------------------
C     Scalar Product (a,mb) where m is acting on the dirac indices only.
C     ------------------------------------------------------------------
C     Part of package:   HMC2DFU1
C     ------------------------------------------------------------------
C     Input parameters:  a,b,length,m
C     ------------------------------------------------------------------
C     Function result: <a|m*1|b>
C     ------------------------------------------------------------------
C     ******************************************************************

      INTEGER length,ndrc
      COMPLEX*16 a(*), b(*)
      COMPLEX*16 m(ndrc,ndrc)
      
      COMPLEX*16 sum
      
      INTEGER i,j,k,olen
      
C     BEGIN
       olen=length/ndrc
       sum=dcmplx(0.0, 0.0)
       DO i=0,olen-1
        DO j=1,ndrc
         DO k=1,ndrc
 	  sum=sum+dconjg(a(i*ndrc+j))*m(j,k)*b(i*ndrc+k)
          
         
	 ENDDO
        ENDDO
       ENDDO
       escalp=sum
      END

