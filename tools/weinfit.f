C	----------------------------------------------------------------
	subroutine weinfit(mode, nplat, n, a, a_var, av, var, k_min)
C	----------------------------------------------------------------
C	Weingarten fit - ihip, 10 Sep 97; last modified: 25 Aug 20
C	----------------------------------------------------------------
C	IN mode - choose plateau according to:
C		1: minimal var
C		2: minimal chi^2
C	IN nplat	- number of points in fit-plateau
C	IN n - number of points/masses to be fitted
C	IN (real*8) a(n) - vector with effective masses
C	IN (real*8) a_var(n) - vector with variances
C	OUT (real*8) av - best fit
C	OUT (real*8) var - variance of av
C	OUT k_min - the point where the fit-plateau begins
C	----------------------------------------------------------------
	real*8 a(*), a_var(*), av, var
	real*8 y, ss, chi2
	real*8 chi2_min, y_min

	real*8 var_min
	real*8 var1, var2

	do k = 0, n - nplat

	  y = 0.0d0
	  ss = 0.0d0
	  do i = 1, nplat
	    y = y + a(k + i) / a_var(k + i)
	    ss = ss + 1.0d0 / a_var(k + i)
	  end do
	  y = y / ss

	  if(mode .eq. 1) then
	    var1 = 0.0d0
	    var2 = 0.0d0
	    do i = 1, nplat
	      var1 = var1 + a_var(k + i)
	      var2 = var2 + (a(k + i) - y)**2
	    end do
	    var1 = var1 / dble(nplat)
	    var2 = var2 / dble(nplat - 1)
	    var = var1 + var2
		write(*, '(1x, ''# ['', i3, ''..'', i3, '']'', $)') k + 1,
     &    k + nplat
		write(*, *) y, dsqrt(var)

	    if((k .eq. 0) .or. (var .lt. var_min)) then
	      var_min = var
	      k_min = k
	      y_min = y
	    end if
	  else if(mode .eq. 2) then
	    chi2 = 0.0d0
	    do i = 1, nplat
	      chi2 = chi2 + (a(k + i) - y)**2 / a_var(k + i)
	    end do

	    if((k .eq. 0) .or. (chi2 .lt. chi2_min)) then
	      chi2_min = chi2
	      k_min = k
	      y_min = y
	    end if
	  else
	    write(*, *) 'weinfit: mode ', mode, ' not implemented!'
	    stop
	  end if

	end do

	var1 = 0.0d0
	var2 = 0.0d0
	do i = 1, nplat
	  var1 = var1 + a_var(k_min + i)
	  var2 = var2 + (a(k_min + i) - y_min)**2
	end do
	var1 = var1 / dble(nplat)
	var2 = var2 / dble(nplat - 1)
	var = var1 + var2

	av = y_min
	k_min = k_min + 1
	return
	end


	
