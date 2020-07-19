C     ******************************************************************
      subroutine jackknife(nmeas, jkblocks, c, ac, ab, var_naiv, var)
C     ******************************************************************
C     **                   **                                         **
C     ** JACKKNIFE         **   I. Hip, 20 Jun 96                     **
C     **                   **   Last modified: 04 Jul 96              **
C     ******************************************************************
C     - computes unbiased mean and variance using jackknife analysis
C     - based on M.C.K. Yang, D.H. Robinson: "Understanding and learning
C       statistics by computer", World Scientific, Singapore, 1986 
C     - tested against jackknife.f and tint.f in CBL/mylib/ftn/gen
C     ******************************************************************
C     IN nmeas               - number of measurements
C     IN jkblocks            - number of jackknife subsamples to be used
C     IN (real*4) c          - array with data
C     OUT ac (real*8)        - simple arithmetic mean
C     OUT ab (real*8)        - arithmetic mean adjusted for bias
C     OUT var_naiv (real*8)  - naive variance
C     OUT var (real*8)       - variance (result of jackknife)
C     ******************************************************************
	parameter(max_jkblocks=100)
	real*4 c(*)
	real*8 psum(max_jkblocks), a(max_jkblocks)
	real*8 sum_c, sum_c2, sum_a, sum_a2
	real*8 aa, aa2, ab, ac, ac2, bias, var_naiv, var

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
	  do j = (i - 1) * npart + 1, i * npart
	    psum(i) = psum(i) + c(j)
	    sum_c2 = sum_c2 + c(j)**2
	  end do
	  sum_c = sum_c + psum(i)
	end do
	ac = sum_c / dble(jkmeas)
	ac2 = sum_c2 / dble(jkmeas)
	var_naiv = (ac2 - ac**2) / dble(jkmeas)
	
C     >>> compute averages
	do i = 1, jkblocks
	  a(i) = (sum_c - psum(i)) / dble(jkmeas - npart)
	end do
	
	sum_a = 0.0d0
	sum_a2 = 0.0d0
	do i = 1, jkblocks
	  sum_a = sum_a + a(i)
	  sum_a2 = sum_a2 + a(i)**2
	end do
	aa = sum_a / dble(jkblocks)
	aa2 = sum_a2 / dble(jkblocks)
	bias = dble(jkblocks - 1) * (aa - ac)
	ab = ac - bias
	var = dble(jkblocks - 1) * (aa2 - aa**2)
	
	return
	end


