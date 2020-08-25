c       ----------------------------------------------------------------
	subroutine mpjack(nmeas, jkblocks, ntime, np, trip, vac, con,
     & t, t_var, s, s_var, mode, nplat, wm_t, wm_t_var, wm_s, wm_s_var)
c       ----------------------------------------------------------------
c	  IN integer*4 nmeas - number of measurements
c	  IN integer*4 jkblocks - number of jackknife blocks
c	  IN integer*4 ntime - lattice size in time direction
c	  IN integer*4 np - maximal impuls (0:np)
c	  IN real*8 trip(nmeas, ntime, np) - triplet current correlation
c	  IN real*8 vac(nmeas, ntime, np) - vacuum fluctuations
c	  IN real*8 con(nmeas, np) - connected correcture
c	  OUT real*8 t(ntime - 1, np) - triplet masses
c	  OUT real*8 t_var(ntime - 1, np) - variances of t
c	  OUT real*8 s(ntime - 1, np) - singlet masses
c	  OUT real*8 s_var(ntime - a, np) - variances of s
c	  IN integer*4 mode - choose plateau according to:
c	    1: minimal variance
c	    2: minimal chi^2
c	  IN integer*4 nplat - search for plateau through nplat points
c	  OUT real*8 wm_t(np) - triplet mass by Weingarten fit
c	  OUT real*8 wm_t_var(np) - wm_t variance
c	  OUT real*8 wm_s(np) - singlet mass by Weingarten fit
c	  OUT real*8 wm_s_var(np) - wm_s variance
c       ----------------------------------------------------------------
c	hip, 18 Aug 98; last modified: 26 Aug 98
c       ----------------------------------------------------------------

	parameter(max_jkblocks=25)
	parameter(max_np=MAX_NSPACE / 2)

	real*8 trip(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 vac(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 con(MAX_NMEAS, 0:max_np)

	real*8 mean_t(MAX_NTIME, 0:max_np)
	real*8 var_t(MAX_NTIME, 0:max_np)
	real*8 mean_v(MAX_NTIME, 0:max_np)
	real*8 var_v(MAX_NTIME, 0:max_np)
	real*8 mean_s(MAX_NTIME, 0:max_np)
	real*8 mean_c(0:max_np)
	real*8 var_c(0:max_np)
c	>>> necessary temporary variable to interchange indices
	real*8 mean_tpt(max_jkblocks, MAX_NTIME, 0:max_np)
	real*8 mean_vp(max_jkblocks, MAX_NTIME, 0:max_np)
	real*8 mean_tp(MAX_NTIME, max_jkblocks, 0:max_np)
	real*8 mean_sp(MAX_NTIME, max_jkblocks, 0:max_np)
	real*8 mean_cp(max_jkblocks, 0:max_np)
	real*8 t(MAX_NTIME - 2, 0:max_np)
	real*8 t_var(MAX_NTIME - 2, 0:max_np)
	real*8 s(MAX_NTIME - 2, 0:max_np)
	real*8 s_var(MAX_NTIME - 2, 0:max_np)
	real*8 tp(MAX_NTIME - 2, max_jkblocks, 0:max_np)
	real*8 sp(MAX_NTIME - 2, max_jkblocks, 0:max_np)
	real*8 wm_t(0:max_np), wm_t_var(0:max_np)
	real*8 wm_s(0:max_np), wm_s_var(0:max_np)
	real*8 varpJK
	real*8 av, var

c       >>> loop through impulses and time slices
	do ip = 0, np
	  do it = 1, ntime - 1
	    call pjack(nmeas, jkblocks, trip(1, it, ip), mean_t(it, ip),
     &        var_t(it, ip), mean_tpt(1, it, ip))
	    call pjack(nmeas, jkblocks, vac(1, it, ip), mean_v(it, ip),
     &        var_v(it, ip), mean_vp(1, it, ip))
	  end do
	  call pjack(nmeas, jkblocks, con(1, ip), mean_c(ip),
     &      var_c(ip), mean_cp(1, ip))
          do it = 1, ntime - 1
 	    mean_s(it, ip) = mean_t(it, ip) -
     &        2.0d0 * (mean_c(ip)**2 + mean_v(it, ip))
        write(*, *) 'rrr = ', ip, it, mean_t(it, ip), mean_v(it, ip)
 	    do j = 1, jkblocks
	      mean_tp(it, j, ip) = mean_tpt(j, it, ip)
 	      mean_sp(it, j, ip) = mean_tpt(j, it, ip) -
     &          2.0d0 * (mean_cp(j, ip)**2 + mean_vp(j, it, ip))
	    end do
          end do        

	  call massb(ntime, mean_t(1, ip), t(1, ip))
	  call massb(ntime, mean_s(1, ip), s(1, ip))
	  do j = 1, jkblocks
	    call massb(ntime, mean_tp(1, j, ip), tp(1, j, ip))
	    call massb(ntime, mean_sp(1, j, ip), sp(1, j, ip))
	  end do

	  do it = 1, ntime / 2 - 1
	    t_var(it, ip) = varpJK(jkblocks, MAX_NTIME - 2, tp(it, 1, ip))
	    s_var(it, ip) = varpJK(jkblocks, MAX_NTIME - 2, sp(it, 1, ip))
	  end do

	end do

	write(*, *)
	write(*, *) '#', ' triplet sigma_1 current'

	do ip = 0, np

	write(*, *)
	do it = 1, ntime / 2 - 1
	  write(*, *) it, t(it, ip), dsqrt(t_var(it, ip))
	end do
	call weinfit(mode, nplat, ntime / 2 - 1, t(1, ip),
     &    t_var(1, ip), av, var, k_min)
	write(*, *) '#', av, dsqrt(var), k_min
	write(*, *) '&'
	wm_t(ip) = av
	wm_t_var(ip) = var

	end do

	write(*, *)
	write(*, *) '#', ' singlet sigma_1 current'

	do ip = 0, np

	write(*, *)
	do it = 1, ntime / 2 - 1
	  write(*, *) it, s(it, ip), dsqrt(s_var(it, ip))
	end do
	call weinfit(mode, nplat, ntime / 2 - 1, s(1, ip),
     &    s_var(1, ip), av, var, k_min)
	write(*, *) '#', av, dsqrt(var), k_min
	write(*, *) '&'
	wm_s(ip) = av
	wm_s_var(ip) = var

	end do


	return
	end


c       -----------------------------------------------------
	subroutine pjack(nmeas, jkblocks, c, ac, var_naiv, a)
c       -----------------------------------------------------
c	  integer*4 nmeas - number of measurements
c	  integer*4 jkblocks - number of jackknife blocks
c	  real*8 c(nmeas) - measurements
c	  real*8 ac - average of c
c	  real*8 var_naiv - variance of c
c	  real*8 a(jkblocks) - averages of jkblocks
c       -----------------------------------------------------
c	hip, Aug 98
c       -----------------------------------------------------
	parameter(max_jkblocks=25)
	real*8 c(*)
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


	real*8 function varpJK(jkblocks, n, a)
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
	varpJK = dble(jkblocks - 1) * (aa2 - aa**2)

	return
	end
	
	
	    
