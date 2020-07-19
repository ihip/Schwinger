C     ******************************************************************
      subroutine jackknife_s(nmeas, jkblocks, c, ac, ab, var_naiv, var,
     &  sigma2s, var_sigma2s)
C     ******************************************************************
C     **                   **                                         **
C     ** JACKKNIFE_S       **   I. Hip, 15 Mai 97                     **
C     **                   **   Last modified: 15 Mai 97              **
C     ******************************************************************
C     - computes unbiased mean and variance using jackknife analysis
C     - based on M.C.K. Yang, D.H. Robinson: "Understanding and learning
C       statistics by computer", World Scientific, Singapore, 1986 
C     - tested against jackknife.f and tint.f in CBL/mylib/ftn/gen
C     ******************************************************************
C     IN nmeas                  - number of measurements
C     IN jkblocks               - number of jk subsamples to be used
C     IN (real*4) c             - array with data
C     OUT (real*8) ac           - simple arithmetic mean
C     OUT (real*8) ab           - arithmetic mean adjusted for bias
C     OUT (real*8) var_naiv     - naive variance of ac
C     OUT (real*8) var          - variance of ac (jk)
C     OUT (real*8) sigma2s      - sigma^2_sample
C     OUT (real*8) var_sigma2s  - variance of sigma^2_sample (jk)
C     ******************************************************************
	parameter(max_jkblocks=100)
	real*4 c(*)
	real*8 psum(max_jkblocks), psum2(max_jkblocks)
	real*8 a(max_jkblocks), a2(max_jkblocks)
	real*8 sum_c, sum_c2, sum_a, sum_a2
	real*8 psigma2s(max_jkblocks)
	real*8 aa, aa2, ab, ac, bias, var_naiv, var
	real*8 sum_sigma2s, sum_sigma2s2
	real*8 aa_sigma2s, aa2_sigma2s, ab_sigma2s, sigma2s
	real*8 bias_sigma2s, var_naiv_sigma2s, var_sigma2s
	
	if(jkblocks .gt. max_jkblocks) then
	  jkblocks = max_jkblocks
	  write(*, *) 'Jackknife WARNING: number of parts to big, ',
     &      'set to ', max_jkblocks, '.'
	end if

	if(nmeas .lt. jkblocks) then 
	  write(*, *) 'Jackknife ERROR: to few measurements.'
	  stop
	end if

	irest = mod(nmeas, jkblocks)
	if(irest .ne. 0) then
c         write(*, *) 'Jackknife WARNING: last ', irest,
c     &      ' measurements will be ignored!'
	end if

C       >>> prepare partial sums
	npart = nmeas / jkblocks
	jkmeas = npart * jkblocks
	sum_c = 0.0d0
	sum_c2 = 0.0d0
	do i = 1, jkblocks
	  psum(i) = 0.0d0
	  psum2(i) = 0.0d0
	  do j = (i - 1) * npart + 1, i * npart
	    psum(i) = psum(i) + c(j)
	    psum2(i) = psum2(i) + c(j)**2
	  end do
	  sum_c = sum_c + psum(i)
	  sum_c2 = sum_c2 + psum2(i)
	end do
	ac = sum_c / dble(jkmeas)
	ac2 = sum_c2 / dble(jkmeas)
	sigma2s = ac2 - ac**2
	var_naiv = sigma2s / dble(jkmeas)
	
C     >>> compute averages
	do i = 1, jkblocks
	  a(i) = (sum_c - psum(i)) / dble(jkmeas - npart)
	  a2(i) = (sum_c2 - psum2(i)) / dble(jkmeas - npart)
	  psigma2s(i) = (a2(i) - a(i)**2)
	end do
	
	sum_a = 0.0d0
	sum_a2 = 0.0d0
	sum_sigma2s = 0.0d0
	sum_sigma2s2 = 0.0d0
	do i = 1, jkblocks
	  sum_a = sum_a + a(i)
	  sum_a2 = sum_a2 + a(i)**2
	  sum_sigma2s = sum_sigma2s + psigma2s(i)
	  sum_sigma2s2 = sum_sigma2s2 + psigma2s(i)**2
	end do

	aa = sum_a / dble(jkblocks)
	aa2 = sum_a2 / dble(jkblocks)
	bias = dble(jkblocks - 1) * (aa - ac)
	ab = ac - bias
	var = dble(jkblocks - 1) * (aa2 - aa**2)
	
	aa_sigma2s = sum_sigma2s / dble(jkblocks)
	aa2_sigma2s = sum_sigma2s2 / dble(jkblocks)
	bias_sigma2s = dble(jkblocks - 1) * (aa_sigma2s - sigma2s)
	ab_sigma2s = sigma2s - bias_sigma2s
	var_sigma2s = dble(jkblocks - 1) * (aa2_sigma2s - aa_sigma2s**2)
	
	return
	end


