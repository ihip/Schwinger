C     ******************************************************************
      program datan
C     ******************************************************************
C     **                   **                                         **
C     ** DATAN             **   I. Hip, Mar 96                        **
C     **                   **   Last modified: 02 Jun 97              **
C     ******************************************************************        
C     - DATa ANalysis program - 
C     ******************************************************************

C       >>> needed to read .data file
	character*64 line, c_computable(MAX_NCOMP)
	character*20 v_integrate
	character*32 v_invert
	character*8 date, time
	integer*4 ncomp, nspace, ntime, nsteps, ntherm, nmeas,
     &    istart, iseed
	real*8 beta, akap, eps

C       >>> other variables which are needed
	character*64 data_file_name
	character*5 check
	character*1 cyn

	write(*, *)
	write(*, *) 'DATa ANalysis program (ivh, 02 Jul 97)'
	write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)

	write(*, '(1x, ''Data file name: '', $)')
	read(*, '(a)') data_file_name
	data_file_name(64:64) = ' '
	i = 1
	do while((data_file_name(i:i) .ne. ' ') .and.
     &    (data_file_name(i:i) .ne. '.'))
	  i = i + 1
	end do

	if(data_file_name(i:i) .eq. '.') then
	  check = data_file_name(i:i+4)
c          if(check .ne. '.data') then
c            stop 'Improper input file name!'
c          end if
	else
	  if(i .gt. 58) then
	    stop 'ERROR: Name too long!'
	  else
	    data_file_name = data_file_name(1:i-1)//'.data'
	  end if
	end if
	write(*, *)

c       check = data_file_name(13:17)
c       if(check .ne. '.data') then
c         write(*, *) 'Improper input file name'
c         stop
c       end if

	open(1, file = data_file_name, form = 'unformatted',
     &    status = 'old')

10	call load_header(1, nspace, ntime, nmeas, beta, akap)

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

	write(*, *)
	write(*, *) 'Enter a choice: '
	write(*, *) '  1. Standard statistical analysis'
	write(*, *) '  2. Standard stat. analysis of sorted data'
	write(*, *) '  3. Jackknife analysis'
	write(*, *) '  4. Jackknife analysis of sorted data'
	write(*, *) '  5. Autocorrelation analysis'
	write(*, *) '  6. Data correlation analysis'
	write(*, *) '  7. Create .xmgr file'
	write(*, *) '  8. pbp left & right - create .xmgr file'
	write(*, *) '  9. pbp left & right - sort analysis'
	write(*, *) ' 10. <AB> - <A><B>'
	write(*, *) ' 11. Jackknife analysis (old)'
	write(*, *) ' 12. S^{eff}_f = -2 ln(det M)'
	write(*, *) ' 13. ergo'
	write(*, *) ' 99. Quit'
	read(*, *) ichoice

	if(ichoice .eq. 99) stop

	write(*, *)
	write(*, *) 'The following computables are available:'
	write(*, *)
	do i = 1, ncomp
	  write(*, *) i, '. ', c_computable(i)
	end do
	write(*, *)

	if(ichoice .eq. 1) then
	  call stat_analysis(ncomp, nmeas, c_computable)
	else if(ichoice .eq. 2) then
	  call sort_analysis(ncomp, nmeas, c_computable)
	else if(ichoice .eq. 3) then
	  call jack_analysis(ncomp, nmeas, c_computable)
	else if(ichoice .eq. 4) then
	  call sort_data_jk(ncomp, nmeas, c_computable)
	else if(ichoice .eq. 5) then
	  call autocorr_analysis(ncomp, nmeas)
	else if(ichoice .eq. 6) then
	  call correlation(ncomp, nmeas, c_computable)
	else if(ichoice .eq. 7) then
	  call create_xmgr_file(ncomp, nmeas, data_file_name)
	else if(ichoice .eq. 8) then
	  call pbp_left_right(ncomp, nmeas, data_file_name)
	else if(ichoice .eq. 9) then
	  call sort_pbp_lr(ncomp, nmeas)
	else if(ichoice .eq. 10) then
	  call abminusab(ncomp, nmeas)
	else if(ichoice .eq. 11) then
	  call jack_analysis_old(ncomp, nmeas)
	else if(ichoice .eq. 12) then
	  call s_eff(ncomp, nmeas, c_computable)
	else if(ichoice .eq. 13) then
	  call ergo(ncomp, nmeas, ntime, nspace, c_computable)
	else if(ichoice .eq. 99) then
	  goto 99
	end if

	call load_tail(1)

	write(*,
     &   '(1x, ''Continue analysis of the same data file (y/n)? '', $)')
	read(*, '(a)') cyn
	write(*, *)

	if(cyn .eq. 'y') then
	  rewind(1)
	  go to 10
	end if 

99	close(1)
	end


C     ******************************************************************
      subroutine jack_analysis(ncomp, nmeas, c_computable)
C     ******************************************************************
C     **                   **                                         **
C     ** JACK_ANALYSIS     **   I. Hip, 23 May 97                     **
C     **                   **   Last modified: 02 Jul 97              **
C     ******************************************************************        
C     - this routine computes some basic statistical values:
C
C         x  +/-  J(x)  [ s(x),  J(s)] 
C
C       x    - arithmetic mean
C       J(x) - sigma_(x_mean) computed using jackknife
C       s(x) - sigma_x, standard deviation of sample
C       J(s) - standard dev of s(x) computed using jackknife
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN (character*64) c_computable - names of computables
C     ******************************************************************
	character*64 c_computable(MAX_NCOMP)
	real*4 computable(MAX_NCOMP), cc(MAX_NMEAS), c
	real*8 sum_x, sum_x2, c_min, c_max
	real*8 mean_x, mean_x2, sigma2s, var
	real*8 var_all, var_sigma2s_all, ac, ab, var_naiv
	real*8 var_jk, sigma2s_jk, var_sigma2s

	write(*, '(1x, ''Analyse which computable: '', $)')
	read(*, *) ic
	
	sum_x = 0.0d0
	sum_x2 = 0.0d0

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  c = computable(ic)
	  cc(i) = c
	  sum_x = sum_x + c
	  sum_x2 = sum_x2 + dble(c)**2
	  if(i .eq. 1) then
	    c_min = c
	    c_max = c
	  else
	    if(c .lt. c_min) c_min = c
	    if(c .gt. c_max) c_max = c
	  end if
	end do

	write(*, *)
	write(*, '('' ---------------------------------------------'',
     &    ''-----------------------'')')
	write(*, *) c_computable(ic)
	write(*, *)

	mean_x = sum_x / dble(nmeas)
	mean_x2 = sum_x2 / dble(nmeas)
	sigma2s = dabs(mean_x2 - mean_x**2)
	var = sigma2s / dble(nmeas)

	var_all = 0.0d0
	var_sigma2s_all = 0.0d0
	do i = 0, 15
	  call jackknife_s(nmeas, 10 + i, cc, ac, ab, var_naiv,
     &      var_jk, sigma2s_jk, var_sigma2s) 
	  var_all = var_all + var_jk
	  var_sigma2s_all = var_sigma2s_all + var_sigma2s
	end do
	var_all = var_all / dble(16)
	var_sigma2s_all = var_sigma2s_all / dble(16)
	
c       >>> human readable
        write(*, '(e13.5, '' +/- '', e10.3,
     &    ''  ['', e11.4, '' +/- '', e9.2, '']'')')
     &    mean_x, dsqrt(var_all),
     &    dsqrt(sigma2s), dsqrt(var_sigma2s_all)

	write(*, '('' ---------------------------------------------'',
     &    ''-----------------------'')')

c	write(*, *) dsqrt(sigma2s_jk)

	return
	end


C     ******************************************************************
      subroutine s_eff(ncomp, nmeas, c_computable)
C     ******************************************************************
C     **                   **                                         **
C     ** S_EFF             **   I. Hip, 27 Jun 97                     **
C     **                   **   Last modified: 02 Jul 97              **
C     ******************************************************************        
C     - this routine computes some basic statistical values:
C
C         x  +/-  J(x)  [ s(x),  J(s)] 
C
C       x    - arithmetic mean
C       J(x) - sigma_(x_mean) computed using jackknife
C       s(x) - sigma_x, standard deviation of sample
C       J(s) - standard dev of s(x) computed using jackknife
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN (character*64) c_computable - names of computables
C     ******************************************************************
	character*64 c_computable(MAX_NCOMP)
	real*4 computable(MAX_NCOMP), cc(MAX_NMEAS)
	real*8 sum_x, sum_x2, c_min, c_max, c
	real*8 mean_x, mean_x2, sigma2s, var
	real*8 var_all, var_sigma2s_all, ac, ab, var_naiv
	real*8 var_jk, sigma2s_jk, var_sigma2s

	write(*, '(1x, ''Which computable is det M: '', $)')
	read(*, *) ic
	
	sum_x = 0.0d0
	sum_x2 = 0.0d0

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  c = - 2.0d0 * dlog(dble(computable(ic)))
	  cc(i) = real(c)
	  sum_x = sum_x + c
	  sum_x2 = sum_x2 + c**2
	  if(i .eq. 1) then
	    c_min = c
	    c_max = c
	  else
	    if(c .lt. c_min) c_min = c
	    if(c .gt. c_max) c_max = c
	  end if
	end do

	write(*, *)
	write(*, '('' ---------------------------------------------'',
     &    ''-----------------------'')')
	write(*, *) c_computable(ic)
	write(*, *)

	mean_x = sum_x / dble(nmeas)
	mean_x2 = sum_x2 / dble(nmeas)
	sigma2s = dabs(mean_x2 - mean_x**2)
	var = sigma2s / dble(nmeas)

	var_all = 0.0d0
	var_sigma2s_all = 0.0d0
	do i = 0, 15
	  call jackknife_s(nmeas, 10 + i, cc, ac, ab, var_naiv,
     &      var_jk, sigma2s_jk, var_sigma2s) 
	  var_all = var_all + var_jk
	  var_sigma2s_all = var_sigma2s_all + var_sigma2s
	end do
	var_all = var_all / dble(16)
	var_sigma2s_all = var_sigma2s_all / dble(16)
	
c       >>> human readable
        write(*, '(e13.5, '' +/- '', e10.3,
     &    ''  ['', e11.4, '' +/- '', e9.2, '']'')')
     &    mean_x, dsqrt(var_all),
     &    dsqrt(sigma2s), dsqrt(var_sigma2s_all)

	write(*, '('' ---------------------------------------------'',
     &    ''-----------------------'')')

c	write(*, *) dsqrt(sigma2s_jk)

	return
	end

C     ******************************************************************
      subroutine abminusab(ncomp, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** ABMINUSAB         **   I. Hip, 18 May 96                     **
C     **                   **   Last modified: 21 Feb 97              **
C     ******************************************************************        
C     - this routine computes <AB> - <A><B>
C     ******************************************************************
	
	real*4 computable(MAX_NCOMP), tcpu
	real*8 sum_a, sum_b, sum_ab, sum_a2, sum_b2, sum_ab2
	real*8 mean_a, mean_b, mean_ab
	
	write(*, '(1x, ''Number of A operator: '', $)')
	read(*, *) n_a
	write(*, '(1x, ''Number of B operator: '', $)')
	read(*, *) n_b

	sum_a = 0.0d0
	sum_a2 = 0.0d0
	sum_b = 0.0d0
	sum_b2 = 0.0d0
	sum_ab = 0.0d0
	sum_ab2 = 0.0d0

	write(*, *) 'Reading data...'
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  sum_a = sum_a + computable(n_a)
	  sum_a2 = sum_a2 + computable(n_a)**2
	  sum_b = sum_b + computable(n_b)
	  sum_b2 = sum_b2 + computable(n_b)**2
	  sum_ab = sum_ab + computable(n_a) * computable(n_b)
	  sum_a2 = sum_a2 + (computable(n_a) * computable(n_b))**2
	end do
	write(*, *) '...done.'
	write(*, *)

	mean_a = sum_a / dble(nmeas)
	mean_b = sum_b / dble(nmeas)
	mean_ab = sum_ab / dble(nmeas)
	
c tu bi trebalo promisliti o pogreshci...
c        do j = 1, ncomp
c          var(j) = dsqrt(
c     &      dabs( (mean_x2(j) - mean_x(j)**2) / dble(nmeas)) )
c        end do

	write(*, *) '<AB> - <A><B> = ', mean_ab - mean_a * mean_b
	
	end


C     ******************************************************************
      subroutine stat_analysis(ncomp, nmeas, c_computable)
C     ******************************************************************
C     **                   **                                         **
C     ** STAT_ANALYSIS     **   I. Hip, Mar 96                        **
C     **                   **   Last modified: 25 Jan 97              **
C     ******************************************************************        
C     - this routine computes some basic statistical values:
C         arithmetic mean
C         naive statistical error
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN (character*64) c_computable - names of computables
C     ******************************************************************

	character*32 hostname
	character*64 c_computable(MAX_NCOMP)
	real*4 computable(MAX_NCOMP), tcpu
	real*8 sum_x(MAX_NCOMP), sum_x2(MAX_NCOMP)
	real*8 mean_x(MAX_NCOMP), mean_x2(MAX_NCOMP), var(MAX_NCOMP)
	real*4 c_min(MAX_NCOMP), c_max(MAX_NCOMP)

	do j = 1, ncomp
	  sum_x(j) = 0.0d0
	  sum_x2(j) = 0.0d0
	end do

	write(*, *) 'Reading data...'
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  do j = 1, ncomp
	    if(i .eq. 1) then
	      c_min(j) = computable(j)
	      c_max(j) = computable(j)
	    else
	      if(computable(j) .lt. c_min(j)) c_min(j) = computable(j)
	      if(computable(j) .gt. c_max(j)) c_max(j) = computable(j)
	    end if
	    sum_x(j) = sum_x(j) + computable(j)
	    sum_x2(j) = sum_x2(j) + dble(computable(j))**2
	  end do
	end do

	write(*, *) '...done.'
	write(*, *)

	do j = 1, ncomp
	  mean_x(j) = sum_x(j) / dble(nmeas)
	  mean_x2(j) = sum_x2(j) / dble(nmeas)
	end do

C       >>> dabs is added because of limited fp-precision when var
C           should be exactly zero 
	do j = 1, ncomp
	  var(j) = dsqrt(dabs(mean_x2(j) - mean_x(j)**2) / dble(nmeas))
	end do

	write(*, *)
     &'----------------------------------------------------------------'
	do j = 1, ncomp
	  write(*, '(1x, i3, ''. Computable: '', a)') j, c_computable(j)
	  write(*, *)
	  write(*, '(1x, e13.5, '' +/- '', e10.3, ''  (min: '',
     &      e10.3, ''  max: '', e10.3, '')'')')
     &      mean_x(j), var(j), c_min(j), c_max(j)
	  write(*, *)
     &'----------------------------------------------------------------'
	end do
	
	end


C     ******************************************************************
      subroutine sort_analysis(ncomp, nmeas, c_computable)
C     ******************************************************************
C     **                   **                                         **
C     ** SORT_ANALYSIS     **   I. Hip, 07 Apr 96                     **
C     **                   **   Last modified: 21 Feb 97              **
C     ******************************************************************        
C     - this routine computes some basic statistical values for data
C       sorted according some discrete variable (e.g. top. charge)
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN (character*64) c_computable - names of computables
C     ******************************************************************
	parameter(max_k=16)

	character*64 c_computable(MAX_NCOMP)
	real*4 computable(MAX_NCOMP), tcpu
	real*8 sum_x(MAX_NCOMP, -max_k:max_k) 
	real*8 sum_x2(MAX_NCOMP, -max_k:max_k)
	real*8 mean_x, mean_x2, var
	integer n(-max_k:max_k)

	write(*, *) 'Sort according to which computable (negative'
	write(*, *) 'value means according to absolute value):'
	read(*, *) is

	ni = 0
	do k = -max_k, max_k
	  n(k) = 0
	  do i = 1, ncomp
	    sum_x(i, k) = 0.0d0
	    sum_x2(i, k) = 0.0d0
	  end do
	end do
	
	write(*, *)
	write(*, *) 'Reading data...'

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(is .lt. 0) then
	    k = nint(abs(computable(-is)))
	  else
	    k = nint(computable(is))
	  end if
	  if(iabs(k) .gt. max_k) then
	    write(*, *) 'k = ', k, ' Not included in statistics.'
	    ni = ni + 1
	  else
	    n(k) = n(k) + 1
	    do j = 1, ncomp
	      sum_x(j, k) = sum_x(j, k) + computable(j) 
	      sum_x2(j, k) = sum_x2(j, k) + dble(computable(j))**2
	    end do
	  end if
	end do

	write(*, *) '...done.'
	write(*, *)
	do j = 1, ncomp
	  write(*, *) '---------------------------------------------'
	  write(*, *) c_computable(j)
	  do k = -max_k, max_k
	    if(n(k) .ne. 0) then
	      mean_x = sum_x(j, k) / dble(n(k))
	      mean_x2 = sum_x2(j, k) / dble(n(k))
	      var = dsqrt(dabs(mean_x2 - mean_x**2) / dble(n(k)))
c             write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
c     &          e13.5)') k, n(k), mean_x, var
	      write(*, '(1x, i4, i6, e13.5, e13.5)')
     &          k, n(k), mean_x, var
	    end if
	  end do
	end do
	write(*, *) '---------------------------------------------'
	if(ni .gt. 0) then
	  write(*, *) ni, ' measurements not included in statistics.'        
	end if

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
	    sum_x2(k) = sum_x2(k) + dble(c(i))**2
	  end if
	end do

	write(*, *)
	write(*, *) '---------------------------------------------'
	do k = -max_k, max_k
	  if(n(k) .ne. 0) then
	    mean_x = sum_x(k) / dble(n(k))
	    mean_x2 = sum_x2(k) / dble(n(k))
	    var = dsqrt(dabs(mean_x2 - mean_x**2) / dble(n(k)))
c           >>> human readable
            write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
     &          e13.5)') k, n(k), mean_x, var
c           >>> xmgr readable 
c	    write(2, '(1x, i4, e13.5, e13.5, i6)')
c     &        k, mean_x, var, n(k)
	  end if
	end do
	write(*, *) '---------------------------------------------'
	if(ni .gt. 0) then
	  write(*, *) ni, ' measurements not included in statistics.'        
	end if

	return
	end


C     ******************************************************************
      subroutine sort_data_jk(ncomp, nmeas, c_computable)
C     ******************************************************************
C     **                   **                                         **
C     ** SORT_DATA_JK      **   I. Hip, 13 May 97                     **
C     ** v3                **   Last modified: 15 May 97              **
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

	write(*, '(1x, ''Sort which computable: '', $)')
	read(*, *) ic
	write(*, *) 'Sort according to which computable (negative'
	write(*,
     &    '(1x, ''value means according to absolute value): '', $)')
	read(*, *) is

	do k = -max_k, max_k
	  n(k) = 0
	  sum_x(k) = 0.0d0
	  sum_x2(k) = 0.0d0
	end do
	
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
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
	end do

	write(*, *)
	write(*, '('' ---------------------------------------------'',
     &    ''-----------------------'')')
	write(*, *) c_computable(ic)
	write(*, *)
	do k = -max_k, max_k
	  if(n(k) .gt. 0) then
	    mean_x(k) = sum_x(k) / dble(n(k))
	    mean_x2 = sum_x2(k) / dble(n(k))
	    sigma2s(k) = dabs(mean_x2 - mean_x(k)**2)
	    var = sigma2s(k) / dble(n(k))
c           >>> human readable
            write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
     &        e10.3, ''  ('', e10.3, '', '', e10.3, '')'')')
     &        k, n(k), mean_x(k), dsqrt(var), c_min(k), c_max(k)
c           >>> xmgr readable 
c	    write(2, '(1x, i4, e13.5, e13.5, i6)')
c     &        k, mean_x(k), var, n(k)
	  end if
	end do
	write(*, '('' ---------------------------------------------'',
     &    ''-----------------------'')')
	write(*, *) '                                 jackknife'
	write(*, *)
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
	
c           >>> human readable
c            write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
c     &        e10.3, ''  ('', e10.3, '', '', e10.3, '')'')')
c     &        k, n(k), mean_x(k), dsqrt(var_all), c_min(k), c_max(k)
             write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
     &         e10.3, ''  ['', e10.3, '' +/- '', e9.2, '']'')')
     &         k, n(k), mean_x(k), dsqrt(var_all), dsqrt(sigma2s(k)),
     &         dsqrt(var_sigma2s_all)
c           >>> xmgr readable 
c	    write(2, '(1x, i4, e13.5, e13.5, i6)')
c     &        k, mean_x(k), dsqrt(var_all), n(k)
	  else if(n(k) .ge. 10) then
	    call jackknife_s(n(k), n(k), ck(1, k), ac, ab, var_naiv,
     &        var_jk, sigma2s_jk, var_sigma2s) 

c           >>> human readable
c            write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
c     &        e10.3, ''  ('', e10.3, '', '', e10.3, '')'')')
c     &        k, n(k), mean_x(k), dsqrt(var_jk), c_min(k), c_max(k)
             write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
     &         e10.3, ''  ['', e10.3, '' +/- '', e9.2, '']'')')
     &         k, n(k), mean_x(k), dsqrt(var_jk), dsqrt(sigma2s(k)),
     &         dsqrt(var_sigma2s)
c           >>> xmgr readable 
c	    write(2, '(1x, i4, e13.5, e13.5, i6)')
c     &        k, mean_x(k), dsqrt(var_jk), n(k)	    
	  else if(n(k) .gt. 0) then
c           >>> human readable
c            write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
c     &        '' jk failed'', ''  ('', e10.3, '', '', e10.3, '')'')')
            write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
     &        '' jk failed'', ''  ['', e10.3, '' +/- jk failed]'')')
     &        k, n(k), mean_x(k), dsqrt(sigma2s(k))
c           >>> xmgr readable 
c	    write(2, '(1x, i4, e13.5, e13.5, i6)')
c     &        k, mean_x(k), dsqrt(var_all), n(k)	    
	  end if
	end do
	write(*, '('' ---------------------------------------------'',
     &    ''-----------------------'')')
	if(ni .gt. 0) then
	  write(*, *) ni, ' measurements not included in statistics.'        
	end if

	return
	end


C     ******************************************************************
      subroutine autocorr_analysis(ncomp, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** AUTOCORR_ANALYSIS **   I. Hip, 07 Apr 96                     **
C     **                   **   Last modified: 21 Feb 97              **
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
	write(*, *) 'Compute autocorrelation for which computable: '
	read(*, *) ic

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  c(i) = computable(ic)
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

	return
	end

C     ******************************************************************
      real*8 function dlsq_k(n, x, y)
C     ******************************************************************
C     **                   **                                         **
C     ** DLSQ_K            **   I. Hip, 23 Jun 96                     **
C     **                   **                                         **
C     ******************************************************************
C     - computes least-squares approximation ...
C     ******************************************************************
C     IN n - number of points
C     IN x - x-coordinates of points
C     IN y - y-coordinates of points
C     ******************************************************************
	real*8 x(*), y(*)
	real*8 sum_x, sum_y, sum_xy, sum_x2 
	data sum_x, sum_y, sum_xy, sum_x2/4*0.0d0/

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


C     ******************************************************************
      subroutine jack_analysis_old(ncomp, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** JACK_ANALYSIS_OLD **   I. Hip, 23 Jun 96                     **
C     **                   **   Last modified: 21 Feb 97              **
C     ******************************************************************        
C     - statistical analysis using jackknife
C     ******************************************************************
	real c(MAX_NMEAS), tcpu, computable(MAX_NCOMP)
	real*8 a, ab, var_naiv, var, tau3

	if(nmeas .gt. MAX_NMEAS) then
	  write(*, *) 'ERROR: nmeas greater than MAX_NMEAS'
	  stop
	end if

	write(*, *) 'Analyse which computable: '
	read(*, *) ic
c       write(*, *) 'How many blocks (e.g. 20): '
c       read(*, *) jkblocks
	
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  c(i) = computable(ic)
	end do
     
        do jkblocks = 5, 25
	  call jackknife(nmeas, jkblocks, c, a, ab, var_naiv, var)
	  write(*, *)
	  write(*, *) 'Blocks: ', jkblocks
          write(*, *) 'Naively:   ', a, ' +/- ', dsqrt(var_naiv)
          write(*, *) 'Jackknife: ', ab, ' +/- ', dsqrt(var)
	  tau3 = var / (2.0d0 * var_naiv)
	  write(*, *) 'which means that tau = ', tau3
        end do

	return
	end


C     ******************************************************************
      subroutine correlation(ncomp, nmeas, c_computable)
C     ******************************************************************
C     **                   **                                         **
C     ** CORRELATION       ** I. Hip, 05 Apr 96, Last mod.: 21 Feb 97 **
C     **                   **                                         **
C     ******************************************************************
C     - this routine computes (Pearson) correlation between two
C       computables ix and iy
C     ******************************************************************

	character*64 c_computable(MAX_NCOMP)
	real*4 computable(MAX_NCOMP), tcpu
	real*8 sum_x, sum_x2, sum_y, sum_y2, sum_xy
	real*8 corr

	write(*, *) 'Enter the numbers of two computables whose'
	write(*, *) '(Pearson) correlation should be computed'
	write(*, *) '(use negative sign for absolute values):'
	write(*, '(1x, ''1. '', $)')
	read(*, *) ix
	write(*, '(1x, ''2. '', $)')
	read(*, *) iy
	
	sum_x = 0.0d0
	sum_x2 = 0.0d0
	sum_y = 0.0d0
	sum_y2 = 0.0d0
	sum_xy = 0.0d0

	write(*, *)
	write(*, *) 'Reading data...'
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(ix .lt. 0) then
	    sum_x = sum_x + abs(computable(-ix))
	  else
	    sum_x = sum_x + computable(ix)
	  end if
	  sum_x2 = sum_x2 + computable(iabs(ix))**2
	  if(iy .lt. 0) then
	    sum_y = sum_y + abs(computable(-iy))
	  else
	    sum_y = sum_y + computable(iy)
	  end if
	  sum_y2 = sum_y2 + computable(iabs(iy))**2
	  if(ix .lt. 0) then
	    if(iy .lt. 0) then
	      sum_xy = sum_xy + abs(computable(-ix)) * 
     &          abs(computable(-iy))
	    else
	      sum_xy = sum_xy + abs(computable(-ix)) * computable(iy)
	    end if
	  else
	    if(iy .lt. 0) then
	      sum_xy = sum_xy + computable(ix) * abs(computable(-iy))
	    else
	      sum_xy = sum_xy + computable(ix) * computable(iy)
	    end if
	  end if
	end do

	write(*, *) '...done.'
c       write(*, *) 'CPU time used for simulation: ', tcpu, ' sec'
	write(*, *)

	corr = (nmeas * sum_xy - sum_x * sum_y) /
     & dsqrt((nmeas * sum_x2 - sum_x**2) * (nmeas * sum_y2 - sum_y**2))

	write(*, *) '=============================================='
	write(*, *)
	write(*, *) 'Correlation between computables'
	write(*, *) '  ', c_computable(iabs(ix))
	write(*, *) '  ', c_computable(iabs(iy))
	write(*, *) 'is ', corr
	write(*, *)
	write(*, *) '=============================================='
	
	end


C     ******************************************************************
      subroutine create_xmgr_file(ncomp, nmeas, data_file_name)
C     ******************************************************************
C     **                   **                                         **
C     ** CREATE_XMGR_FILE  **   I. Hip, Mar 96, Last mod.: 21 Feb 97  **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine should create file which can be read by xmgr
C     ******************************************************************
C       >>> max_out = maximal number of data columns in output file
	parameter(max_out=32)

	character*32 data_file_name
	character*17 xmgr_file_name
	real*4 computable(MAX_NCOMP), tcpu
	integer*4 icomp(max_out)

	write(*, *) 'Enter the numbers of computables which should'
	write(*, *) 'be written to .xmgr file (0 = numbers 1 to nmeas,'
	write(*, *) 'negative sign gives absolute value, 99 = end)'

	i = 1
	read(*, *) n
	do while (n .ne. 99)
	  if(abs(n) .gt. ncomp) then
	    write(*, *) 'This computable does not exist (ignored).'
	  else
	    icomp(i) = n
	    i = i + 1
	    if(i .gt. max_out) then
	       write(*, *) 'Maximal number of output columns exceeded'
	       stop
	    end if
	  end if
	  read(*, *) n
	end do
	ncol = i - 1

	write(*, *)
	write(*, *) '.xmgr will contain computables: '
	write(*, *) (icomp(j), j = 1, ncol)
	write(*, *)
	
	if(ncol .eq. 0) then
	  write(*, *) 'Nothing to do!'
	  stop
	end if

	print *, data_file_name

	xmgr_file_name = data_file_name(1:12)//'.xmgr'
	write(*, *) 'Creating file: ', xmgr_file_name, '... '

	open(2, file = xmgr_file_name, form = 'formatted',
     &    status = 'unknown')

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  do j = 1, ncol
	    if(icomp(j) .eq. 0) then
	      write(2, '(i11, $)') i
	    else if(icomp(j) .lt. 0) then
	      write(2, '(e11.3, $)') abs(computable(iabs(icomp(j))))
	    else
	      write(2, '(e11.3, $)') computable(icomp(j))
	    end if
	  end do
	  write(2, *)
	end do

	close(2)
	write(*, *) '...file created!'

	end


C     ******************************************************************
      subroutine pbp_left_right(ncomp, nmeas, data_file_name)
C     ******************************************************************
C     **                   **                                         **
C     ** PBP_LEFT_RIGHT    ** I. Hip, 09 Apr 96, Last mod.: 21 Feb 97 **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine should create file which contains pbp_left and
C       pbp_right and can be read by xmgr
C     ******************************************************************

	character*32 data_file_name
	character*16 pbp_file_name
	real*4 computable(MAX_NCOMP), tcpu

	write(*, *) 'Enter the numbers of computables pbp and pbg5p'
	write(*, '(1x, ''pbp = '', $)')
	read(*, *) ipbp
	write(*, '(1x, ''pbg5p = '', $)')
	read(*, *) ipbg5p

	pbp_file_name = data_file_name(1:12)//'.pbp'
	write(*, *) 'Creating file: ', pbp_file_name, '... '

	open(2, file = pbp_file_name, form = 'formatted',
     &    status = 'unknown')

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  pbp_left = 0.5 * (computable(ipbp) - computable(ipbg5p))
	  pbp_right = 0.5 * (computable(ipbp) + computable(ipbg5p))
	  write(2, '(i11, 2e12.4)') i, pbp_left, pbp_right
	end do

	close(2)

	write(*, *) '...file created!'

	end

C     ******************************************************************
      subroutine sort_pbp_lr(ncomp, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** SORT_PBP_LR       **   I. Hip, 25 Jun 96                     **
C     **                   **   Last modified: 21 Feb 97              **
C     ******************************************************************        
C     - this routine should compute pbp_left and pbp_right and makes
C       all kinds of statistical analysis on them
C     ******************************************************************

	real*4 computable(MAX_NCOMP), tcpu
	real*4 pbp_left(MAX_NMEAS), pbp_right(MAX_NMEAS), ck(MAX_NMEAS)

	write(*, *) 'Enter the numbers of computables pbp and pbg5p'
	write(*, '(1x, ''pbp = '', $)')
	read(*, *) ipbp
	write(*, '(1x, ''pbg5p = '', $)')
	read(*, *) ipbg5p
	write(*, *) 'Sort according to which computable (negative'
	write(*, *) 'value means according to absolute value):'
	read(*, *) is

	write(*, *) 'Reading data...'
	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  pbp_left(i) = 0.5 * (computable(ipbp) - computable(ipbg5p))
	  pbp_right(i) = 0.5 * (computable(ipbp) + computable(ipbg5p))
	  ck(i) = computable(is)
	end do
	write(*, *) '...done!'

	write(*, *) 'pbp_left'
c        call jackknife(...)
	call sort_data(nmeas, ck, pbp_left, is) 
	
	write(*, *) 'pbp_right'
c        call jackknife(...)
	call sort_data(nmeas, ck, pbp_right, is) 
	
	return
	end


C     ******************************************************************
      subroutine ergo(ncomp, nmeas, ntime, nspace, c_computable)
C     ******************************************************************
C     **                   **                                         **
C     ** SORT_ANALYSIS     **   I. Hip, 07 Apr 96                     **
C     **                   **   Last modified: 21 Feb 97              **
C     ******************************************************************        
C     - this routine computes some basic statistical values for data
C       sorted according some discrete variable (e.g. top. charge)
C     ******************************************************************
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN (character*64) c_computable - names of computables
C     ******************************************************************
	parameter(max_k=16)

	character*64 c_computable(MAX_NCOMP)
	real*4 computable(MAX_NCOMP), tcpu
	real*8 sum_x(-max_k:max_k) 
	real*8 sum_x2(-max_k:max_k)
	real*8 mean_x, mean_x2, var
	integer n(-max_k:max_k)
	real*8 c

	write(*, *) 'Sort according to which computable (negative'
	write(*, *) 'value means according to absolute value):'
	read(*, *) is

	ni = 0
	do k = -max_k, max_k
	  n(k) = 0
	  sum_x(k) = 0.0d0
	  sum_x2(k) = 0.0d0
	end do
	
	write(*, *)
	write(*, *) 'Reading data...'

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ncomp)
	  if(is .lt. 0) then
	    k = nint(abs(computable(-is)))
	  else
	    k = nint(computable(is))
	  end if
	  if(iabs(k) .gt. max_k) then
	    write(*, *) 'k = ', k, ' Not included in statistics.'
	    ni = ni + 1
	  else
	    n(k) = n(k) + 1
	    c = dexp(- dble(computable(4)) * dble(ntime * nspace)) *
     &        dble(computable(26))**2
	    sum_x(k) = sum_x(k) + c 
	    sum_x2(k) = sum_x2(k) + c**2
	  end if
	end do

	write(*, *) '...done.'
	write(*, *)
	write(*, *) '---------------------------------------------'
	write(*, *) c_computable(4)
	write(*, *) c_computable(26)
	do k = -max_k, max_k
	  if(n(k) .ne. 0) then
	    mean_x = sum_x(k) / dble(n(k))
	    mean_x2 = sum_x2(k) / dble(n(k))
	    var = dsqrt(dabs(mean_x2 - mean_x**2) / dble(n(k)))
c           write(*, '(1x, i4, '' ('', i6, '') '', e13.5, '' +/- '',
c     &        e13.5)') k, n(k), mean_x, var
	    write(*, '(1x, i4, i6, e13.5, e13.5)')
     &        k, n(k), mean_x, var
	  end if
	end do
	write(*, *) '---------------------------------------------'
	if(ni .gt. 0) then
	  write(*, *) ni, ' measurements not included in statistics.'        
	end if

	end



