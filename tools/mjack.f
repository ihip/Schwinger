c     ----------------------------------------------------------------
	subroutine mjack(nmeas, jkblocks, ntime, mcomp, c, t1, t1_var,
     &    t3, t3_var, s1, s1_var, mode, nplat)
c     ----------------------------------------------------------------
c	  IN integer*4 nmeas - number of measurements
c	  IN integer*4 jkblocks - number of jackknife blocks
c	  IN integer*4 ntime - lattice size in time direction
c	  IN integer*4 mcomp - number of mass computables
c	  IN real*4 c - measurement data
c	  OUT real*8 t1 - triplet sigma1 mass
c	  OUT real*8 t1_var - variance of t1
c	  OUT real*8 t3 - triplet sigma3 mass
c	  OUT real*8 t3_var - variance of t3
c	  OUT real*8 s1 - singlet sigma1 mass
c	  OUT real*8 s1_var - variance of s1
c	  IN integer*4 mode - choose plateau according to:
c	    1: minimal variance
c	    2: minimal chi^2
c	  IN integer*4 nplat - search for plateau through nplat points
c     ----------------------------------------------------------------
c	  ihip / 18 Sep 97 / last modified: 25 Aug 20
c     ----------------------------------------------------------------

	parameter(max_jkblocks=25)

	real*4 c(MAX_NMEAS, MAX_NTIME, MAX_NCOMP)

	real*8 mean(MAX_NTIME, MAX_NCOMP)
	real*8 var_naiv(MAX_NTIME, MAX_NCOMP)
	real*8 mean_p(MAX_NTIME, MAX_NCOMP, max_jkblocks)
	real*8 triplet(MAX_NTIME, 4), singlet(MAX_NTIME, 4)
	real*8 triplet_p(MAX_NTIME, max_jkblocks, 4)
	real*8 singlet_p(MAX_NTIME, max_jkblocks, 4)
	real*8 triplet_efm_s1(MAX_NTIME - 2)
	real*8 triplet_efm_s1_var(MAX_NTIME - 2)
	real*8 triplet_efm_s3(MAX_NTIME - 2)
	real*8 triplet_efm_s3_var(MAX_NTIME - 2)
	real*8 singlet_efm_s1(MAX_NTIME - 2)
	real*8 singlet_efm_s1_var(MAX_NTIME - 2)
	real*8 triplet_p_efm_s1(MAX_NTIME - 2, max_jkblocks)
	real*8 triplet_p_efm_s3(MAX_NTIME - 2, max_jkblocks)
	real*8 singlet_p_efm_s1(MAX_NTIME - 2, max_jkblocks)
	real*8 varJK
	real*8 av, var
	real*8 t1, t1_var, t3, t3_var, s1, s1_var

c	>>> number of jackknife blocks
c	jkblocks = 10

c	>>> choose plateau according to: (1) minimal var; (2) minimal chi^2
c	mode = 1

c   >>> search for plateau through nplat points (3 or 4 or more?)
c	nplat = 4

c   >>> loop through time slices
	do it = 1, ntime - 1
	  do ic = 1, mcomp
	    call jack(nmeas, jkblocks, c(1, it, ic), mean(it, ic),
     &        var_naiv(it, ic), mean_p(1, it, ic))
	  end do
          do ig = 1, 4
            triplet(it, ig) = 2.0d0 * mean(it, ig)
 	    singlet(it, ig) = 2.0d0 * (mean(it, ig) +
     &        2.0d0 * (mean(it, 8 + ig) * mean(it, 12 + ig) -
     &        mean(it, 4 + ig)))
	    do j = 1, jkblocks
              triplet_p(it, j, ig) = 2.0d0 * mean_p(j, it, ig)
 	      singlet_p(it, j, ig) = 2.0d0 * (mean_p(j, it, ig) +
     &          2.0d0 * (mean_p(j, it, 8 + ig) *
     &          mean_p(j, it, 12 + ig) - mean_p(j, it, 4 + ig)))
	    end do
          end do        
	end do

	call massb(ntime, triplet(1, 2), triplet_efm_s1(1))
	call massb(ntime, triplet(1, 4), triplet_efm_s3(1))
	call massb(ntime, singlet(1, 2), singlet_efm_s1(1))
	do j = 1, jkblocks
	  call massb(ntime, triplet_p(1, j, 2), triplet_p_efm_s1(1, j))
	  call massb(ntime, triplet_p(1, j, 4), triplet_p_efm_s3(1, j))
	  call massb(ntime, singlet_p(1, j, 2), singlet_p_efm_s1(1, j))
	end do

	do it = 1, ntime / 2 - 1
	  triplet_efm_s1_var(it) = 
     &      varJK(jkblocks, MAX_NTIME - 2, triplet_p_efm_s1(it, 1))
	  triplet_efm_s3_var(it) = 
     &      varJK(jkblocks, MAX_NTIME - 2, triplet_p_efm_s3(it, 1))
	  singlet_efm_s1_var(it) = 
     &      varJK(jkblocks, MAX_NTIME - 2, singlet_p_efm_s1(it, 1))
	end do

	write(*, *)
	write(*, *) '# triplet_efm_s1'
	do it = 1, ntime / 2 - 1
	  write(*, *) it, triplet_efm_s1(it),
     & dsqrt(triplet_efm_s1_var(it))
	end do
	call weinfit(mode, nplat, ntime / 2 - 1, triplet_efm_s1,
     & triplet_efm_s1_var, av, var, k_min)
	write(*, *) '#', av, dsqrt(var), k_min
	write(*, *) '&'
	t1 = av
	t1_var = var

	write(*, *)
	write(*, *) '# triplet_efm_s3'
	do it = 1, ntime / 2 - 1
	  write(*, *) it, triplet_efm_s3(it),
     &   dsqrt(triplet_efm_s3_var(it))
	end do
	call weinfit(mode, nplat, ntime / 2 - 1, triplet_efm_s3,
     & triplet_efm_s3_var, av, var, k_min)
	write(*, *) '#', av, dsqrt(var), k_min
	write(*, *) '&'
	t3 = av
	t3_var = var

	write(*, *)
	write(*, *) '# singlet_efm_s1'
	do it = 1, ntime / 2 - 1
	  write(*, *) it, singlet_efm_s1(it),
     &  dsqrt(singlet_efm_s1_var(it))
	end do
	call weinfit(mode, nplat, ntime / 2 - 1, singlet_efm_s1,
     & singlet_efm_s1_var, av, var, k_min)
	write(*, *) '#', av, dsqrt(var), k_min
	write(*, *)
	s1 = av
	s1_var = var

	return
	end


	subroutine jack(nmeas, jkblocks, c, ac, var_naiv, a)
	parameter(max_jkblocks=25)
	real*4 c(*)
	real*8 psum(max_jkblocks), a(max_jkblocks)
	real*8 sum_c, sum_c2, sum_a, sum_a2
	real*8 ac, ac2, var_naiv

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

	return
	end


	real*8 function varJK(jkblocks, n, a)
	real*8 a(n, *)

	real*8 sum_a, sum_a2, aa, aa2

	sum_a = 0.0d0
	sum_a2 = 0.0d0
	do i = 1, jkblocks
	  sum_a = sum_a + a(1, i)
	  sum_a2 = sum_a2 + a(1, i)**2
	end do
	aa = sum_a / dble(jkblocks)
	aa2 = sum_a2 / dble(jkblocks)
	varJK = dble(jkblocks - 1) * (aa2 - aa**2)

	return
	end
	
	
	    
