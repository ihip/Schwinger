C     ******************************************************************
      program gcollector
C     ******************************************************************
C     **                   **                                         **
C     ** DATA COLLECTOR    **   I. Hip, 07 May 96                     **
C     ** v4.1              **   Adapted for gnuplot: 31 Oct 20        **
C     **                   **   Last modified: 13 Nov 20              **
C     ******************************************************************
C     Collects statistics for different .data files:
C     - the program reads the name of the file which contains a list
C       with data files which should be processed
C     - asks for the computable which should be collected and for the 
C       name of the output file
C     - opens data files, computes arithmetic mean and stat. error
C     - creates file with entries: kappa, ar. mean, stat. error, of
C       chosen computable (added are also number of valid measurements
C       and exceptional configurations)
C     - writes on the screen the simulation details for
C       every data file as a check
C     ******************************************************************

C       >>> needed to read .data file
	real*8 beta, akap
	integer*4 nspace, ntime, nmeas, ncomp
	character*64 c_computable(MAX_NCOMP)

C       >>> other variables which are needed
	character*1 cjack, cch
	character*64 data_file_name
	character*5 check
	character*32 datname, outname

	write(*, *)
	write(*, *) 'Data collector v4.1 [gnuplot] (ihip, 13 Nov 20)'
	write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)
	write(*, '(1x, ''File list file name: '', $)')
	read(*, '(a)') datname
	write(*, '(1x, ''Collect which computable: '', $)')
	read(*, *) kcomp

	write(*, *)
	write(*, *) ' 1. no sort - naive error'
	write(*, *) ' 2. no sort - jackknife'
	write(*, *) ' 3. sort - all channels - naive'
	write(*, *) ' 4. sort - all channels - jackknife'
	write(*, *) ' 5. sort - given channel - jackknife'
	write(*, *) ' 6. S_eff (sort, all, jk)'
	write(*, *) ' 7. A, B (no sort, jk)'
	write(*, *) ' 8. chi_A, chi_B (no sort, jk)'
	write(*, *) ' 9. A, B (sort, all, jk)'
	write(*, *) '10. A, B (sort, given channel, jk)'
	write(*, *) '11. chi_A, chi_B (sort, given channel, jk)'
	write(*, *) '12. W'
	write(*, *) '13. autocorrelation [NEW in 4.1]'
	write(*, *) '99. Quit'
	write(*, *)

	write(*, '(''> '', $)')
	read(*, *) ichoice

	if(ichoice .ne. 13) then
	  write(*, '(1x, ''Number of control comp. (usually 2): '', $)')
	  read(*, *) kcon
	end if

	write(*, *)
	write(*, '(1x, ''Output file name: '', $)')
	read(*, '(a)') outname
	write(*, *)

	if((ichoice .eq. 3) .or. (ichoice .eq. 4) .or. (ichoice .eq. 5)
     &    .or. (ichoice .eq. 6) .or. (ichoice .eq. 9)) then
	  write(*, '(1x, ''Sort acc. to which comp. (0 no sort,'')')
          write(*, '(1x,''  negative means acc. to abs. value): '', $)')
	  read(*, *) is
	end if

	if((ichoice .eq. 5)) then
	  write(*, '(1x, ''Which channel: '', $)')
	  read(*, *) iws
	end if	
	if((ichoice .eq. 10) .or. (ichoice .eq. 11)) then
	  iws = 1
	end if

	open(3, file = datname, form = 'formatted', status = 'old')
	open(2, file = outname, form = 'formatted', status = 'unknown')

c	if(is .ne. 0) then
c	  write(2, '(''#'',
c     &''-----------------------------------------------------------'')')
c	end if

10      read(3, '(a)', end = 99) data_file_name

	open(1, file = data_file_name, form = 'unformatted',
     &    status = 'old')

        call load_header(1, nspace, ntime, nmeas, beta, akap)

	if(nmeas .gt. MAX_NMEAS) then
	  write(*, *) 'nmeas is greater than MAX_NMEAS'
	  write(*, *) 'nmeas = ', nmeas, '  MAX_NMEAS = ', MAX_NMEAS
	  stop
	end if

        call load_dheader(1, ncomp, c_computable)
	if(ncomp .gt. MAX_NCOMP) then
	  write(*, *) 'ncomp is greater than MAX_NCOMP'
	  stop
	end if

	write(*, *) '  ', c_computable(kcomp)
	write(2, '(''#'', ''  '', a)') c_computable(kcomp)
	write(2, '(2i4, f6.2, f8.4, $)') nspace, ntime, beta, akap

	if(ichoice .eq. 1) then
	  call coll_stat(akap, ncomp, nmeas, kcomp, kcon)
	else if(ichoice .eq. 2) then
	  call coll_stat_jk(beta, akap, ncomp, nmeas, kcomp, kcon)
	else if(ichoice .eq. 3) then
	  call coll_sort(ncomp, nmeas, kcomp, is)
	else if(ichoice .eq. 4) then
	  call coll_sort_jk(ncomp, nmeas, kcomp, is, c_computable,
     &      kcon)
	else if(ichoice .eq. 5) then
	   call coll_sort_jk1(beta, akap, ncomp, nmeas, kcomp, is, iws,
     &       c_computable, kcon)
	else if(ichoice .eq. 6) then
	   call coll_sort_jk_seff(ntime, nspace, ncomp, nmeas,
     &       kcomp, is, c_computable, kcon)
	else if(ichoice .eq. 7) then
	   call coll4_stat_jk(beta, akap, ncomp, nmeas, kcomp, kcon)
	else if(ichoice .eq. 8) then
	   call coll4abs_stat_jk(beta, akap, ncomp, nmeas, kcomp, kcon)
 	else if(ichoice .eq. 9) then
	  call coll4_sort_jk(ncomp, nmeas, kcomp, is, c_computable,
     &      kcon)
	else if(ichoice .eq. 10) then
	  call coll4_sort1_jk(beta, akap, ncomp, nmeas, kcomp, iws,
     &      kcon)
	else if(ichoice .eq. 11) then
	  call coll4abs_sort1_jk(beta, akap, ncomp, nmeas, kcomp,
     &      iws, kcon)
 	else if(ichoice .eq. 12) then
	  call wformel(beta, akap, ntime, nspace, ncomp, nmeas,
     &      kcomp, kcon)
	else if(ichoice .eq. 13) then
	  call autocorr_analysis(ncomp, nmeas, kcomp)
	else
	  goto 99
	end if

c	write(2, '(''&'',
c     &''-----------------------------------------------------------'')')

	call load_tail(1)
	close(1)
	goto 10

99      close(3)
	close(2)
	write(*, *) 'ok!'
	end

C     ******************************************************************
      subroutine coll_stat(akap, ncomp, nmeas, kcomp, kniter)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL_STAT         **   I. Hip, 09 Mai 96                     **
C     ** v3                **   Last modified: 21 Mar 97              **
C     **                   **                                         **
C     ******************************************************************
C     - this routine computes some basic statistical values:
C         arithmetic mean
C         naive statistical error
C     ******************************************************************
C     IN (real*8) akap - actual kappa (used only for output)
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN kcomp - number of computable which should be analysed
C     IN knit - number of niter computable 
C     ******************************************************************
	real*8 akap
	real*4 computable(MAX_NCOMP)
	real*8 sum_x, sum_x2
	real*8 mean_x, mean_x2, sd_naiv
	real*4 c_min, c_max

C       >>> nexcp - number of exceptional configurations
	nexcp = 0

	sum_x = 0.0d0
	sum_x2 = 0.0d0

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(kniter) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	    if(i .eq. 1) then
	      c_min = computable(kcomp)
	      c_max = computable(kcomp)
	    else
	      if(computable(kcomp) .lt. c_min) c_min = computable(kcomp)
	      if(computable(kcomp) .gt. c_max) c_max = computable(kcomp)
	    end if
	    sum_x = sum_x + computable(kcomp)
	    sum_x2 = sum_x2 + computable(kcomp)**2
	  end if
	end do

	nmeas = nmeas - nexcp

	mean_x = sum_x / dble(nmeas)
	mean_x2 = sum_x2 / dble(nmeas)
	sd_naiv = dsqrt((mean_x2 - mean_x**2) / dble(nmeas))

	write(2, '(1x, e16.9, e11.3, i8, i6)')
     & mean_x, sd_naiv, nmeas, nexcp
	write(*, *) 'valid meas.: ', nmeas, '  except.: ', nexcp
	write(*, *)
	
	end


C     ******************************************************************
      subroutine coll_stat_jk(beta, akap, ncomp, nmeas, kcomp, knit)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL_STAT_JK      **   I. Hip, 04 Jul 96                     **
C     ** v3                **   Last modified: 24 May 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine prepares for the call of jackknife subroutine
C     ******************************************************************
C     IN (real*8) akap - actual kappa (used only for output)
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN kcomp - number of computable which should be analysed
C     IN knit - number of niter computable 
C     ******************************************************************
	real*8 beta, akap, ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 aa, sum_x
	real*4 computable(MAX_NCOMP)
	real*4 c(MAX_NMEAS)

C       >>> nexcp - number of exceptional configurations
	nexcp = 0

	k = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(knit) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	    k = k + 1
	    c(k) = computable(kcomp)
	  end if
	end do

	sum_x = 0.0d0
	do i = 1, k
	  sum_x = sum_x + c(i)
	end do
	aa = sum_x / dble(k)

c	>>> number of jk blocks: from 10 to 25 and then averaged!
	var_all = 0.0d0
	var_sigma2s_all = 0.0d0
	do i = 0, 15
	  call jackknife_s(k, 10 + i, c, ac, ab, var_naiv, var_jk,
     & sigma2s_jk, var_sigma2s) 
	  var_all = var_all + var_jk
	  var_sigma2s_all = var_sigma2s_all + var_sigma2s
	end do
	var_all = var_all / dble(16)
	var_sigma2s_all = var_sigma2s_all / dble(16)
	
	write(2, '(1x, e16.9, e11.3, i7, i4)')
     & aa, dsqrt(var_all), k, nexcp

	write(*, *) 'valid meas.: ', k, '  except.: ', nexcp
	write(*, *)

	return
	end


C     ******************************************************************
      subroutine coll4_stat_jk(beta, akap, ncomp, nmeas, kcomp, knit)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL4_STAT_JK     **   I. Hip, 04 Jul 96                     **
C     ** v3                **   Last modified: 02 Jul 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine prepares for the call of jackknife subroutine
C     ******************************************************************
C     IN (real*8) akap - actual kappa (used only for output)
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN kcomp - number of computable which should be analysed
C     IN knit - number of niter computable 
C     ******************************************************************
	real*8 beta, akap, ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 aa, sum_x
	real*4 computable(MAX_NCOMP)
	real*4 c(4 * MAX_NMEAS)

C       >>> nexcp - number of exceptional configurations
	nexcp = 0

	ktop = 5

	k = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(knit) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	    itop = iabs(nint(computable(ktop)))
	    k = k + 1
	    c(k) = computable(kcomp)
	    if(itop .gt. 1) then
	      k = k + 1
	      c(k) = computable(kcomp + 1)
	      if(itop .gt. 2) then
	        k = k + 1
	        c(k) = computable(kcomp + 2)
	        if(itop .gt. 3) then
	          k = k + 1
	          c(k) = computable(kcomp + 3)
		end if
	      end if
	    end if
	  end if
	end do

	sum_x = 0.0d0
	do i = 1, k
	  sum_x = sum_x + c(i)
	end do
	aa = sum_x / dble(k)

c	>>> number of jk blocks: from 10 to 25 and then averaged!
	var_all = 0.0d0
	var_sigma2s_all = 0.0d0
	do i = 0, 15
	  call jackknife_s(k, 10 + i, c, ac, ab, var_naiv, var_jk,
     & sigma2s_jk, var_sigma2s) 
	  var_all = var_all + var_jk
	  var_sigma2s_all = var_sigma2s_all + var_sigma2s
	end do
	var_all = var_all / dble(16)
	var_sigma2s_all = var_sigma2s_all / dble(16)
	
c	if(akap .eq. 0.0) then
	  write(2, '(1x, f10.6, e16.8, e11.3, e12.4, e10.2, i8, i4)')
     & beta, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
     & dsqrt(var_sigma2s_all), k, nexcp
c	else
c	  write(2, '(1x, f10.7, e16.8, e11.3, e12.4, e10.2, i8, i4)')
c     & akap, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
c     & dsqrt(var_sigma2s_all), k, nexcp
c	end if

	write(*, *) 'valid meas.: ', k, '  except.: ', nexcp
	write(*, *)

	return
	end


C     ******************************************************************
      subroutine coll4_sort1_jk(beta, akap, ncomp, nmeas, kcomp, iws, 
     & knit)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL4_SORT1_JK    **   I. Hip, 04 Jul 96                     **
C     ** v3                **   Last modified: 02 Jul 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine prepares for the call of jackknife subroutine
C     ******************************************************************
C     IN (real*8) akap - actual kappa (used only for output)
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN kcomp - number of computable which should be analysed
C     IN knit - number of niter computable 
C     ******************************************************************
	real*8 beta, akap, ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 aa, sum_x
	real*4 computable(MAX_NCOMP)
	real*4 c(4 * MAX_NMEAS)

C       >>> nexcp - number of exceptional configurations
	nexcp = 0

	ktop = 5

	k = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(knit) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	    itop = nint(computable(ktop))
	    if(itop .eq. iws) then
	      k = k + 1
	      c(k) = computable(kcomp)
	    end if
	  end if
	end do

	sum_x = 0.0d0
	do i = 1, k
	  sum_x = sum_x + c(i)
	end do
	aa = sum_x / dble(k)

c	>>> number of jk blocks: from 10 to 25 and then averaged!
	var_all = 0.0d0
	var_sigma2s_all = 0.0d0
	do i = 0, 15
	  call jackknife_s(k, 10 + i, c, ac, ab, var_naiv, var_jk,
     & sigma2s_jk, var_sigma2s) 
	  var_all = var_all + var_jk
	  var_sigma2s_all = var_sigma2s_all + var_sigma2s
	end do
	var_all = var_all / dble(16)
	var_sigma2s_all = var_sigma2s_all / dble(16)
	
c	if(akap .eq. 0.0) then
	  write(2, '(1x, f10.6, e16.8, e11.3, e12.4, e10.2, i8, i4)')
     & beta, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
     & dsqrt(var_sigma2s_all), k, nexcp
c	else
c	  write(2, '(1x, f10.7, e16.8, e11.3, e12.4, e10.2, i8, i4)')
c     & akap, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
c     & dsqrt(var_sigma2s_all), k, nexcp
c	end if

	write(*, *) 'valid meas.: ', k, '  except.: ', nexcp
	write(*, *)

	return
	end

C     ******************************************************************
      subroutine wformel(beta, akap, ntime, nspace,ncomp, nmeas,
     &  kcomp, knit)
C     ******************************************************************
C     **                   **                                         **
C     ** WFORMEL           **   I. Hip, 09 Jul 96                     **
C     ** v3                **   Last modified: 09 Jul 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine prepares for the call of jackknife subroutine
C     ******************************************************************
C     IN (real*8) akap - actual kappa (used only for output)
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN kcomp - number of computable which should be analysed
C     IN knit - number of niter computable 
C     ******************************************************************
	real*8 beta, akap, ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 aa, sum_x
	real*4 computable(MAX_NCOMP)
	real*4 c(4 * MAX_NMEAS)

	real*8 a1, a2

C       >>> nexcp - number of exceptional configurations
	nexcp = 0

	ktop = 5

	k = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(knit) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	    itop = iabs(nint(computable(ktop)))
	    k = k + 1
	    c(k) = 2.0 / (abs(computable(14)) *
     & (1.0 - (akap * computable(10))**2)) - 
     & 2.0 / (abs(computable(22)) *
     & (1.0 - (akap * computable(18))**2))

	a1 = c(k) / dble(ntime * nspace)
	a2 = abs(computable(28))

	write(*, '(1x, i3, 3e12.3, f7.2)') itop, a1, a2, a1 - a2, a1 / a2

	    if(itop .gt. 1) then
	      k = k + 1
	    c(k) = 2.0 / (abs(computable(15)) *
     & (1.0 - (akap * computable(11))**2)) - 
     & 2.0 / (abs(computable(23)) *
     & (1.0 - (akap * computable(19))**2))
	      if(itop .gt. 2) then
	        k = k + 1
	    c(k) = 2.0 / (abs(computable(16)) *
     & (1.0 - (akap * computable(12))**2)) - 
     & 2.0 / (abs(computable(24)) *
     & (1.0 - (akap * computable(20))**2))
	        if(itop .gt. 3) then
	          k = k + 1
	    c(k) = 2.0 / (abs(computable(17)) *
     & (1.0 - (akap * computable(13))**2)) - 
     & 2.0 / (abs(computable(25)) *
     & (1.0 - (akap * computable(21))**2))
		end if
	      end if
	    end if
	  end if
	end do

	sum_x = 0.0d0
	do i = 1, k
	  sum_x = sum_x + c(i)
	end do
	aa = sum_x / dble(k)

c	>>> number of jk blocks: from 10 to 25 and then averaged!
	var_all = 0.0d0
	var_sigma2s_all = 0.0d0
	do i = 0, 15
	  call jackknife_s(k, 10 + i, c, ac, ab, var_naiv, var_jk,
     & sigma2s_jk, var_sigma2s) 
	  var_all = var_all + var_jk
	  var_sigma2s_all = var_sigma2s_all + var_sigma2s
	end do
	var_all = var_all / dble(16)
	var_sigma2s_all = var_sigma2s_all / dble(16)
	
c	if(akap .eq. 0.0) then
	  write(2, '(1x, f10.6, e16.8, e11.3, e12.4, e10.2, i8, i4)')
     & beta, aa / dble(ntime * nspace),
     & dsqrt(var_all) / dble(ntime * nspace),
     & dsqrt(sigma2s_jk) / dble(ntime * nspace),
     & dsqrt(var_sigma2s_all) / dble(ntime * nspace),
     & k, nexcp
c	else
c	  write(2, '(1x, f10.7, e16.8, e11.3, e12.4, e10.2, i8, i4)')
c     & akap, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
c     & dsqrt(var_sigma2s_all), k, nexcp
c	end if

	write(*, *) 'valid meas.: ', k, '  except.: ', nexcp
	write(*, *)

	return
	end


C     ******************************************************************
      subroutine coll4abs_stat_jk(beta, akap, ncomp, nmeas, kcomp, knit)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL4ABS_STAT_JK  **   I. Hip, 04 Jul 96                     **
C     ** v3                **   Last modified: 02 Jul 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine prepares for the call of jackknife subroutine
C     ******************************************************************
C     IN (real*8) akap - actual kappa (used only for output)
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN kcomp - number of computable which should be analysed
C     IN knit - number of niter computable 
C     ******************************************************************
	real*8 beta, akap, ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 aa, sum_x
	real*4 computable(MAX_NCOMP)
	real*4 c(4 * MAX_NMEAS)

C       >>> nexcp - number of exceptional configurations
	nexcp = 0

	ktop = 5

	k = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(knit) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	    itop = iabs(nint(computable(ktop)))
	    k = k + 1
	    c(k) = abs(computable(kcomp))
	    if(itop .gt. 1) then
	      k = k + 1
	      c(k) = abs(computable(kcomp + 1))
	      if(itop .gt. 2) then
	        k = k + 1
	        c(k) = abs(computable(kcomp + 2))
	        if(itop .gt. 3) then
	          k = k + 1
	          c(k) = abs(computable(kcomp + 3))
		end if
	      end if
	    end if
	  end if
	end do

	sum_x = 0.0d0
	do i = 1, k
	  sum_x = sum_x + c(i)
	end do
	aa = sum_x / dble(k)

c	>>> number of jk blocks: from 10 to 25 and then averaged!
	var_all = 0.0d0
	var_sigma2s_all = 0.0d0
	do i = 0, 15
	  call jackknife_s(k, 10 + i, c, ac, ab, var_naiv, var_jk,
     & sigma2s_jk, var_sigma2s) 
	  var_all = var_all + var_jk
	  var_sigma2s_all = var_sigma2s_all + var_sigma2s
	end do
	var_all = var_all / dble(16)
	var_sigma2s_all = var_sigma2s_all / dble(16)
	
c	if(akap .eq. 0.0) then
	  write(2, '(1x, f10.6, e16.8, e11.3, e12.4, e10.2, i8, i4)')
     & beta, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
     & dsqrt(var_sigma2s_all), k, nexcp
c	else
c	  write(2, '(1x, f10.7, e16.8, e11.3, e12.4, e10.2, i8, i4)')
c     & akap, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
c     & dsqrt(var_sigma2s_all), k, nexcp
c	end if

	write(*, *) 'valid meas.: ', k, '  except.: ', nexcp
	write(*, *)

	return
	end


C     ******************************************************************
      subroutine coll4abs_sort1_jk(beta, akap, ncomp, nmeas, kcomp, iws,
     & knit)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL4ABS_STAT_JK  **   I. Hip, 04 Jul 96                     **
C     ** v3                **   Last modified: 02 Jul 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine prepares for the call of jackknife subroutine
C     ******************************************************************
C     IN (real*8) akap - actual kappa (used only for output)
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN kcomp - number of computable which should be analysed
C     IN knit - number of niter computable 
C     ******************************************************************
	real*8 beta, akap, ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 aa, sum_x
	real*4 computable(MAX_NCOMP)
	real*4 c(4 * MAX_NMEAS)

C       >>> nexcp - number of exceptional configurations
	nexcp = 0

	ktop = 5

	k = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(knit) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	    itop = nint(computable(ktop))
	    if(itop .eq. iws) then
	      k = k + 1
	      c(k) = computable(kcomp)
	    end if
	  end if
	end do

	sum_x = 0.0d0
	do i = 1, k
	  sum_x = sum_x + c(i)
	end do
	aa = sum_x / dble(k)

c	>>> number of jk blocks: from 10 to 25 and then averaged!
	var_all = 0.0d0
	var_sigma2s_all = 0.0d0
	do i = 0, 15
	  call jackknife_s(k, 10 + i, c, ac, ab, var_naiv, var_jk,
     & sigma2s_jk, var_sigma2s) 
	  var_all = var_all + var_jk
	  var_sigma2s_all = var_sigma2s_all + var_sigma2s
	end do
	var_all = var_all / dble(16)
	var_sigma2s_all = var_sigma2s_all / dble(16)
	
c	if(akap .eq. 0.0) then
	  write(2, '(1x, f10.6, e16.8, e11.3, e12.4, e10.2, i8, i4)')
     & beta, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
     & dsqrt(var_sigma2s_all), k, nexcp
c	else
c	  write(2, '(1x, f10.7, e16.8, e11.3, e12.4, e10.2, i8, i4)')
c     & akap, aa, dsqrt(var_all), dsqrt(sigma2s_jk),
c     & dsqrt(var_sigma2s_all), k, nexcp
c	end if

	write(*, *) 'valid meas.: ', k, '  except.: ', nexcp
	write(*, *)

	return
	end

C     ******************************************************************
      subroutine coll_sort(ncomp, nmeas, kcomp, is)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL_SORT         **   I. Hip, 26 Jun 96                     **
C     ** v3                **   Last modified: 21 Mar 97              **
C     ******************************************************************        
C     - this routine prepares for the call of sort_data subroutine
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN kcomp - number of computable which should be sorted
C     IN is - number of comp. according to which is sorting taking place  
C     ******************************************************************
	real*4 computable(MAX_NCOMP)
	real*4 c(MAX_NMEAS), cs(MAX_NMEAS)

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  c(i) = computable(kcomp)
	  cs(i) = computable(is)
	end do

	call sort_data(nmeas, cs, c, is) 
	
	return
	end


C     ******************************************************************
      subroutine sort_data(nmeas, cs, c, iflag)
C     ******************************************************************
C     **                   **                                         **
C     ** SORT_DATA         **   I. Hip, 25 Jun 96                     **
C     ** v3                **   Last modified: 21 Mar 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine computes some basic statistical values for data
C       sorted according some discrete variable (e.g. top. charge)
C     ******************************************************************
C     IN nmeas - number of measurements
C     IN (real*4) cs - key variable (according to which is sorted)
C     IN (real*4) c - variable which should be sorted
C     IN iflag - if iflag < 0 then the absolute value of cs is taken
C     ******************************************************************
	parameter(max_k=16)

	real*4 cs(*), c(*)
	real*8 sum_x(-max_k:max_k) 
	real*8 sum_x2(-max_k:max_k)
	real*8 mean_x, mean_x2, var
	integer n(-max_k:max_k)

C       >>> counts how many measurements are not included in statistics
C       >>> (because they are out of allowed interval)
	ni = 0

	do k = -max_k, max_k
	  n(k) = 0
	  sum_x(k) = 0.0d0
	  sum_x2(k) = 0.0d0
	end do
	
	do i = 1, nmeas
	  if(iflag .lt. 0) then
	    k = nint(abs(cs(i)))
	  else
	    k = nint(cs(i))
	  end if
	  if(iabs(k) .gt. max_k) then
	    write(*, *) 'k = ', k, ' Not included in statistics.'
	    ni = ni + 1
	  else
	    n(k) = n(k) + 1
	    sum_x(k) = sum_x(k) + c(i) 
	    sum_x2(k) = sum_x2(k) + c(i)**2
	  end if
	end do

	write(*, *)
	write(2, *)
	do k = -max_k, max_k
	  if(n(k) .ne. 0) then
	    mean_x = sum_x(k) / dble(n(k))
	    mean_x2 = sum_x2(k) / dble(n(k))
	    var = dsqrt(dabs(mean_x2 - mean_x**2) / dble(n(k)))
C           >>> human readable
C            write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
C     &          e13.5)') k, n(k), mean_x, var
C           >>> xmgr readable
	    write(2, '(1x, i4, e13.5, e13.5, i6)')
     &        k, mean_x, var, n(k)
	  end if
	end do
	if(ni .gt. 0) then
	  write(*, *) ni, ' measurements not included in statistics.'        
	  write(2, *) ni, ' measurements not included in statistics.'        
	end if

	return
	end


C     ******************************************************************
      subroutine coll_sort_jk(ncomp, nmeas, ic, is, c_computable, kcon)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL_SORT_JK      **   I. Hip, 24 May 97                     **
C     ** v3                **   Last modified: 02 Jun 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine computes some basic statistical values for data
C       sorted according some discrete variable (e.g. top. charge)
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     ******************************************************************
	parameter(max_k=16)

	character*64 c_computable(MAX_NCOMP)
	real*4 cs, c, computable(MAX_NCOMP)
	real*4 ck(MAX_NMEAS, -max_k:max_k)
	real*4 c_min(-max_k:max_k), c_max(-max_k:max_k) 
	real*8 sum_x(-max_k:max_k) 
	real*8 sum_x2(-max_k:max_k)
	real*8 mean_x(-max_k:max_k), mean_x2, var
	real*8 sigma2s(-max_k:max_k)
	integer n(-max_k:max_k)
	real*8 ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 dble, dreal, dabs

C       >>> counts how many measurements are not included in statistics
C       >>> (because they are out of allowed interval)
	ni = 0

	do k = -max_k, max_k
	  n(k) = 0
	  sum_x(k) = 0.0d0
	  sum_x2(k) = 0.0d0
	end do
	
	nexcp = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(kcon) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	  c = computable(ic)
	  cs = computable(is)
	  if(iflag .lt. 0) then
	    k = nint(abs(cs))
	  else
	    k = nint(cs)
	  end if
	  if(iabs(k) .gt. max_k) then
	    write(*, *) 'k = ', k, ' Not included in statistics.'
	    ni = ni + 1
	  else
	    n(k) = n(k) + 1
	    ck(n(k), k) = c
	    sum_x(k) = sum_x(k) + c
	    sum_x2(k) = sum_x2(k) + dble(c)**2
	    if(n(k) .eq. 1) then
	      c_min(k) = c
	      c_max(k) = c
	    else
	      if(c .lt. c_min(k)) c_min(k) = c
	      if(c .gt. c_max(k)) c_max(k) = c
	    end if
	  end if
	  end if
	end do

	do k = -max_k, max_k
	  if(n(k) .gt. 0) then
	    mean_x(k) = sum_x(k) / dble(n(k))
	    mean_x2 = sum_x2(k) / dble(n(k))
	    sigma2s(k) = dabs(mean_x2 - mean_x(k)**2)
	    var = sigma2s(k) / dble(n(k))
	  end if
	end do

	write(2, *)
	do k = -max_k, max_k
	  if(n(k) .ge. 25) then
	    var_all = 0.0d0
	    var_sigma2s_all = 0.0d0
	    do i = 0, 15
	      call jackknife_s(n(k), 10 + i, ck(1, k), ac, ab, var_naiv,
     &          var_jk, sigma2s_jk, var_sigma2s) 
	      var_all = var_all + var_jk
	      var_sigma2s_all = var_sigma2s_all + var_sigma2s
	    end do
	    var_all = var_all / dble(16)
	    var_sigma2s_all = var_sigma2s_all / dble(16)
	
            write(2, '(1x, i6, e13.5, e11.3, e11.3,  e11.3, i6)')
     &        k, mean_x(k), dsqrt(var_all), dsqrt(sigma2s(k)),
     &        dsqrt(var_sigma2s_all), n(k)
	  else if(n(k) .ge. 10) then
	    call jackknife_s(n(k), n(k), ck(1, k), ac, ab, var_naiv,
     &        var_jk, sigma2s_jk, var_sigma2s) 

            write(2, '(1x, i6, e13.5, e11.3, e11.3,  e11.3, i6)')
     &        k, mean_x(k), dsqrt(var_jk), dsqrt(sigma2s(k)),
     &        dsqrt(var_sigma2s), n(k)
	  else if(n(k) .gt. 0) then
            write(2, '(''#'', i6, e13.5, ''  jk failed'', e11.3,
     & ''  jk failed'', i6)') k, mean_x(k), dsqrt(sigma2s(k)), n(k)
	  end if
	end do
	if(ni .gt. 0) then
	  write(2, '(''#'', i6,
     &      '' measurement(s) out of range.'')') ni        
	end if
	if(nexcp .gt. 0) then
	  write(2, '(''#'', i6,
     &      '' excp. measurement(s) not included in stat. '')') nexcp        
	end if

	return
	end


C     ******************************************************************
      subroutine coll4_sort_jk(ncomp, nmeas, ic, is, c_computable, kcon)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL4_SORT_JK     **   I. Hip, 24 May 97                     **
C     ** v3                **   Last modified: 02 Jul 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine computes some basic statistical values for data
C       sorted according some discrete variable (e.g. top. charge)
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     ******************************************************************
	parameter(max_k=16)

	character*64 c_computable(MAX_NCOMP)
	real*4 cs, c, computable(MAX_NCOMP)
	real*4 ck(MAX_NMEAS, -max_k:max_k)
	real*4 c_min(-max_k:max_k), c_max(-max_k:max_k) 
	real*8 sum_x(-max_k:max_k) 
	real*8 sum_x2(-max_k:max_k)
	real*8 mean_x(-max_k:max_k), mean_x2, var
	real*8 sigma2s(-max_k:max_k)
	integer n(-max_k:max_k)
	real*8 ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 dble, dreal, dabs

C       >>> counts how many measurements are not included in statistics
C       >>> (because they are out of allowed interval)
	ni = 0

	do k = -max_k, max_k
	  n(k) = 0
	  sum_x(k) = 0.0d0
	  sum_x2(k) = 0.0d0
	end do
	
	nexcp = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(kcon) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	  c = computable(ic)
	  cs = computable(is)
	  if(iflag .lt. 0) then
	    k = nint(abs(cs))
	  else
	    k = nint(cs)
	  end if
	  if(iabs(k) .gt. max_k) then
	    write(*, *) 'k = ', k, ' Not included in statistics.'
	    ni = ni + 1
	  else
	    n(k) = n(k) + 1
	    ck(n(k), k) = c
	    if(iabs(k) .gt. 1) then
	      n(k) = n(k) + 1
	      ck(n(k), k) = computable(ic + 1)
	      if(iabs(k) .gt. 2) then
	        n(k) = n(k) + 1
	        ck(n(k), k) = computable(ic + 2)
	        if(iabs(k) .gt. 3) then
	          n(k) = n(k) + 1
	          ck(n(k), k) = computable(ic + 3)
		end if
	      end if
	    end if
	  end if
	  end if
	end do

	do k = -max_k, max_k
	  if(n(k) .gt. 0) then
	    sum_x(k) = 0.0d0
	    sum_x2(k) = 0.0d0
	    do i = 1, n(k)
	      c = ck(i, k)
	      sum_x(k) = sum_x(k) + c
	      sum_x2(k) = sum_x2(k) + dble(c)**2
	      if(i .eq. 1) then
	        c_min(k) = c
	        c_max(k) = c
	      else
	        if(c .lt. c_min(k)) c_min(k) = c
	        if(c .gt. c_max(k)) c_max(k) = c
	      end if
	    end do
	    mean_x(k) = sum_x(k) / dble(n(k))
	    mean_x2 = sum_x2(k) / dble(n(k))
	    sigma2s(k) = dabs(mean_x2 - mean_x(k)**2)
	    var = sigma2s(k) / dble(n(k))
	  end if
	end do

	do k = -max_k, max_k
	  if(n(k) .ge. 25) then
	    var_all = 0.0d0
	    var_sigma2s_all = 0.0d0
	    do i = 0, 15
	      call jackknife_s(n(k), 10 + i, ck(1, k), ac, ab, var_naiv,
     &          var_jk, sigma2s_jk, var_sigma2s) 
	      var_all = var_all + var_jk
	      var_sigma2s_all = var_sigma2s_all + var_sigma2s
	    end do
	    var_all = var_all / dble(16)
	    var_sigma2s_all = var_sigma2s_all / dble(16)
	
            write(2, '(1x, i6, e13.5, e11.3, e11.3,  e11.3, i6)')
     &        k, mean_x(k), dsqrt(var_all), dsqrt(sigma2s(k)),
     &        dsqrt(var_sigma2s_all), n(k)
	  else if(n(k) .ge. 10) then
	    call jackknife_s(n(k), n(k), ck(1, k), ac, ab, var_naiv,
     &        var_jk, sigma2s_jk, var_sigma2s) 

            write(2, '(1x, i6, e13.5, e11.3, e11.3,  e11.3, i6)')
     &        k, mean_x(k), dsqrt(var_jk), dsqrt(sigma2s(k)),
     &        dsqrt(var_sigma2s), n(k)
	  else if(n(k) .gt. 0) then
            write(2, '(''#'', i6, e13.5, ''  jk failed'', e11.3,
     & ''  jk failed'', i6)') k, mean_x(k), dsqrt(sigma2s(k)), n(k)
	  end if
	end do
	if(ni .gt. 0) then
	  write(2, '(''#'', i6,
     &      '' measurement(s) out of range.'')') ni        
	end if
	if(nexcp .gt. 0) then
	  write(2, '(''#'', i6,
     &      '' excp. measurement(s) not included in stat. '')') nexcp        
	end if

	return
	end


C     ******************************************************************
      subroutine coll_sort_jk_seff(ntime, nspace, ncomp, nmeas, ic, is,
     & c_computable, kcon)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL_SORT_JK_SEFF **   I. Hip, 02 Jul 97                     **
C     ** v3                **   Last modified: 02 Jul 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine computes some basic statistical values for data
C       sorted according some discrete variable (e.g. top. charge)
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     ******************************************************************
	parameter(max_k=16)

	character*64 c_computable(MAX_NCOMP)
	real*4 cs, c, computable(MAX_NCOMP)
	real*4 ck(MAX_NMEAS, -max_k:max_k)
	real*4 c_min(-max_k:max_k), c_max(-max_k:max_k) 
	real*8 sum_x(-max_k:max_k) 
	real*8 sum_x2(-max_k:max_k)
	real*8 mean_x(-max_k:max_k), mean_x2, var
	real*8 sigma2s(-max_k:max_k)
	integer n(-max_k:max_k)
	real*8 ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 dble, dreal, dabs

C       >>> counts how many measurements are not included in statistics
C       >>> (because they are out of allowed interval)
	ni = 0

	do k = -max_k, max_k
	  n(k) = 0
	  sum_x(k) = 0.0d0
	  sum_x2(k) = 0.0d0
	end do
	
	nexcp = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(kcon) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	  c = -2.0 * log(computable(ic)) / (ntime * nspace)
	  cs = computable(is)
	  if(iflag .lt. 0) then
	    k = nint(abs(cs))
	  else
	    k = nint(cs)
	  end if
	  if(iabs(k) .gt. max_k) then
	    write(*, *) 'k = ', k, ' Not included in statistics.'
	    ni = ni + 1
	  else
	    n(k) = n(k) + 1
	    ck(n(k), k) = c
	    sum_x(k) = sum_x(k) + c
	    sum_x2(k) = sum_x2(k) + dble(c)**2
	    if(n(k) .eq. 1) then
	      c_min(k) = c
	      c_max(k) = c
	    else
	      if(c .lt. c_min(k)) c_min(k) = c
	      if(c .gt. c_max(k)) c_max(k) = c
	    end if
	  end if
	  end if
	end do

	do k = -max_k, max_k
	  if(n(k) .gt. 0) then
	    mean_x(k) = sum_x(k) / dble(n(k))
	    mean_x2 = sum_x2(k) / dble(n(k))
	    sigma2s(k) = dabs(mean_x2 - mean_x(k)**2)
	    var = sigma2s(k) / dble(n(k))
	  end if
	end do

	do k = -max_k, max_k
	  if(n(k) .ge. 25) then
	    var_all = 0.0d0
	    var_sigma2s_all = 0.0d0
	    do i = 0, 15
	      call jackknife_s(n(k), 10 + i, ck(1, k), ac, ab, var_naiv,
     &          var_jk, sigma2s_jk, var_sigma2s) 
	      var_all = var_all + var_jk
	      var_sigma2s_all = var_sigma2s_all + var_sigma2s
	    end do
	    var_all = var_all / dble(16)
	    var_sigma2s_all = var_sigma2s_all / dble(16)
	
            write(2, '(1x, i6, e13.5, e11.3, e11.3,  e11.3, i6)')
     &        k, mean_x(k), dsqrt(var_all), dsqrt(sigma2s(k)),
     &        dsqrt(var_sigma2s_all), n(k)
	  else if(n(k) .ge. 10) then
	    call jackknife_s(n(k), n(k), ck(1, k), ac, ab, var_naiv,
     &        var_jk, sigma2s_jk, var_sigma2s) 

            write(2, '(1x, i6, e13.5, e11.3, e11.3,  e11.3, i6)')
     &        k, mean_x(k), dsqrt(var_jk), dsqrt(sigma2s(k)),
     &        dsqrt(var_sigma2s), n(k)
	  else if(n(k) .gt. 0) then
            write(2, '(''#'', i6, e13.5, ''  jk failed'', e11.3,
     & ''  jk failed'', i6)') k, mean_x(k), dsqrt(sigma2s(k)), n(k)
	  end if
	end do
	if(ni .gt. 0) then
	  write(2, '(''#'', i6,
     &      '' measurement(s) out of range.'')') ni        
	end if
	if(nexcp .gt. 0) then
	  write(2, '(''#'', i6,
     &      '' excp. measurement(s) not included in stat. '')') nexcp        
	end if

	return
	end

C     ******************************************************************
      subroutine coll_sort_jk1(beta, akap, ncomp, nmeas, ic, is, iws,
     & c_computable, kcon)
C     ******************************************************************
C     **                   **                                         **
C     ** COLL_SORT_JK1     **   I. Hip, 02 Jun 97                     **
C     ** v3                **   Last modified: 02 Jun 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine computes some basic statistical values for data
C       sorted according some discrete variable (e.g. top. charge)
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     ******************************************************************
	parameter(max_k=16)

	real*8 beta, akap
	character*64 c_computable(MAX_NCOMP)
	real*4 cs, c, computable(MAX_NCOMP)
	real*4 ck(MAX_NMEAS, -max_k:max_k)
	real*4 c_min(-max_k:max_k), c_max(-max_k:max_k) 
	real*8 sum_x(-max_k:max_k) 
	real*8 sum_x2(-max_k:max_k)
	real*8 mean_x(-max_k:max_k), mean_x2, var
	real*8 sigma2s(-max_k:max_k)
	integer n(-max_k:max_k)
	real*8 ac, ab, var_naiv, var_jk, var_all
	real*8 sigma2s_jk, var_sigma2s, var_sigma2s_all
	real*8 dble, dreal, dabs

C       >>> counts how many measurements are not included in statistics
C       >>> (because they are out of allowed interval)
	ni = 0

	do k = -max_k, max_k
	  n(k) = 0
	  sum_x(k) = 0.0d0
	  sum_x2(k) = 0.0d0
	end do
	
	nexcp = 0
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(computable(kcon) .lt. 0) then
	    nexcp = nexcp + 1
	  else
	  c = computable(ic)
	  cs = computable(iabs(is))
	  if(is .lt. 0) then
	    k = nint(abs(cs))
	  else
	    k = nint(cs)
	  end if
	  if(iabs(k) .gt. max_k) then
	    write(*, *) 'k = ', k, ' Not included in statistics.'
	    ni = ni + 1
	  else
	    n(k) = n(k) + 1
	    ck(n(k), k) = c
	    sum_x(k) = sum_x(k) + c
	    sum_x2(k) = sum_x2(k) + dble(c)**2
	    if(n(k) .eq. 1) then
	      c_min(k) = c
	      c_max(k) = c
	    else
	      if(c .lt. c_min(k)) c_min(k) = c
	      if(c .gt. c_max(k)) c_max(k) = c
	    end if
	  end if
	  end if
	end do

	k = iws
	  if(n(k) .gt. 0) then
	    mean_x(k) = sum_x(k) / dble(n(k))
	    mean_x2 = sum_x2(k) / dble(n(k))
	    sigma2s(k) = dabs(mean_x2 - mean_x(k)**2)
	    var = sigma2s(k) / dble(n(k))
	  end if
 
	k = iws 
	  if(n(k) .ge. 25) then
	    var_all = 0.0d0
	    var_sigma2s_all = 0.0d0
	    do i = 0, 15
	      call jackknife_s(n(k), 10 + i, ck(1, k), ac, ab, var_naiv,
     &          var_jk, sigma2s_jk, var_sigma2s) 
	      var_all = var_all + var_jk
	      var_sigma2s_all = var_sigma2s_all + var_sigma2s
	    end do
	    var_all = var_all / dble(16)
	    var_sigma2s_all = var_sigma2s_all / dble(16)
	
          write(2, '(1x, 2f10.4, 2i6, e13.5, e11.3, e11.3,  e11.3, i7)')
     &        beta, akap, k, n(k), mean_x(k), dsqrt(var_all),
     &        dsqrt(sigma2s(k)), dsqrt(var_sigma2s_all), nmeas
	  else if(n(k) .ge. 10) then
	    call jackknife_s(n(k), n(k), ck(1, k), ac, ab, var_naiv,
     &        var_jk, sigma2s_jk, var_sigma2s) 

          write(2, '(1x, 2f10.4, 2i6, e13.5, e11.3, e11.3,  e11.3, i7)')
     &        beta, akap, k, n(k), mean_x(k), dsqrt(var_jk),
     &        dsqrt(sigma2s(k)), dsqrt(var_sigma2s), nmeas
	  else
            write(2, '(1x, 2f10.4, 2i6, e13.5, ''  jk failed'', e11.3,
     &        ''  jk failed'', i7)') beta, akap, k, n(k), mean_x(k),
     &        dsqrt(sigma2s(k)), nmeas
	  end if

	if(ni .gt. 0) then
	  write(2, '(''#'', i6,
     &      '' measurement(s) out of range.'')') ni        
	end if
	if(nexcp .gt. 0) then
	  write(2, '(''#'', i6,
     &      '' excp. measurement(s) not included in stat. '')') nexcp        
	end if

	return
	end


c >>> ADDED IN VERSION 4.1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C     ******************************************************************
      subroutine autocorr_analysis(ncomp, nmeas, kcomp)
C     ******************************************************************
C     **                   **                                         **
C     ** AUTOCORR_ANALYSIS **   I. Hip, 07 Apr 96                     **
C     **                   **   Last modified: 13 Nov 20 (ihip)       **
C     ******************************************************************
C     - computes autocorrelation for given computable
C     ******************************************************************
	parameter(max_ac=1000)
	real c(MAX_NMEAS), tcpu, computable(MAX_NCOMP)
	real*8 x(max_ac), acln(max_ac), dlsq_k
	real*8 sum_x, sum_x2, sum_y, sum_y2, sum_xy
	real*8 ac, ac_last, sum_ac
	real*8 a, ab, var_naiv, var, tau1, tau2, tau3

	if(nmeas .gt. MAX_NMEAS) then
	  write(*, *) 'ERROR: nmeas greater than MAX_NMEAS'
	  stop
	end if

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  c(i) = computable(kcomp)
	end do

	write(*, *)
	ac_last = 1.0d0
	sum_ac = 0.5d0
	k = 1
10        n = nmeas - k
	  sum_x = 0.0d0
	  sum_x2 = 0.0d0
	  sum_y = 0.0d0
	  sum_y2 = 0.0d0
	  sum_xy = 0.0d0
	  do i = 1, n
	    sum_x = sum_x + c(i)
	    sum_x2 = sum_x2 + c(i)**2
	    sum_y = sum_y + c(i + k)
	    sum_y2 = sum_y2 + c(i + k)**2
	    sum_xy = sum_xy + c(i) * c(i + k)
	  end do
	  ac = (n * sum_xy - sum_x * sum_y) / 
     &      dsqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))
	  write(*, *) k, ac
	  if((ac .le. 0.0d0) .or. (ac .ge. ac_last)) goto 99
	  acln(k) = dlog(ac)
	  sum_ac = sum_ac + ac
	  ac_last = ac
	  k = k + 1
	  if(k .gt. max_ac) then
	    write(*, *) 'Autoccorelation WARNING: ',
     &        'ac longer than max_ac.'
	    goto 99
	  end if
	goto 10

99      write(*, *)
	
	do i = 1, k
	  x(i) = dble(i)
	end do
	tau1 = -1.0d0 / dlsq_k(k - 1, x, acln)
	write(*, *) 'Least squares fit:     tau = ', tau1

	tau2 = sum_ac
	write(*, *) 'Cumulative method:     tau = ', tau2

	call jackknife(nmeas, 20, c, a, ab, var_naiv, var)
	tau3 = var / (2.0d0 * var_naiv)
	write(*, *) 'Jackknife (20 blocks): tau = ', tau3

	write(2, '(3f7.3)') tau1, tau2, tau3

	return
	end

C     ******************************************************************
      real*8 function dlsq_k(n, x, y)
C     ******************************************************************
C     **                   **                                         **
C     ** DLSQ_K            **   I. Hip, 23 Jun 96                     **
C     **                   **   Last modified: 13 Nov 20              **
C     ******************************************************************
C     - computes least-squares approximation ...
C     ******************************************************************
C     IN n - number of points
C     IN x - x-coordinates of points
C     IN y - y-coordinates of points
C     ******************************************************************
	real*8 x(*), y(*)
	real*8 sum_x, sum_y, sum_xy, sum_x2 

	  sum_x = 0.0d0
	  sum_y = 0.0d0
	  sum_xy = 0.0d0
	  sum_x2 = 0.0d0

	do i = 1, n
c          write(*, *) i, x(i), y(i)
	  sum_x = sum_x + x(i)
	  sum_y = sum_y + y(i)
	  sum_xy = sum_xy + x(i) * y(i)
	  sum_x2 = sum_x2 + x(i)**2
	end do

	dlsq_k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x**2)
	return
	end

