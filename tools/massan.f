C     ******************************************************************
      program massan
C     ******************************************************************
C     **                   **                                         **
C     ** MASSAN            **   I. Hip, Sep 96                        **
C     ** v4                **   Last modified: 25 Aug 20              **
C     **                   **                                         **
C     ******************************************************************        
C     - MASS ANalysis program - 
C     ******************************************************************

C       >>> needed to read .mass file
	real*8 beta, akap
	character*64 c_mcomputable(MAX_MCOMP)
	integer*4 nspace, ntime, nmeas, ncomp

C       >>> other variables which are needed
	character*64 data_file_name
	character*5 check

	write(*, *)
	write(*, *) 'MASS ANalysis program (ihip, 25 Aug 20)'
	write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)

	write(*, '(1x, ''Mass file name: '', $)')
	read(*, '(a)') data_file_name
	write(*, *)

c       >>> check of the file name - not really necessary
c	data_file_name(64:64) = ' '
c	i = 1
c	do while((data_file_name(i:i) .ne. ' ') .and.
c    &    (data_file_name(i:i) .ne. '.'))
c	  i = i + 1
c	end do
c        if(data_file_name(i:i) .eq. '.') then
c          check = data_file_name(i:i+4)
c          if(check .ne. '.mass') then
c            stop 'Improper input file name!'
c          end if
c        else
c          if(i .gt. 58) then
c           stop 'ERROR: Name too long!'
c          else
c            data_file_name = data_file_name(1:i-1)//'.mass'
c          end if
c        end if
c        write(*, *)

	open(1, file = data_file_name, form = 'unformatted',
     &    status = 'old')

        call load_header(1, nspace, ntime, nmeas, beta, akap)

c	write(*, *) '!!! load header ok !!!'

        call load_mheader(1, mcomp)
	  write(*, *) 'mcomp = ', mcomp
	if(mcomp .gt. MAX_MCOMP) then
	  write(*, *) 'mcomp is greater than MAX_MCOMP'
	  stop
	end if

	write(*, *)
	write(*, *) 'Enter a choice: '
	write(*, *) '  1. Mass analysis (naive error)'
	write(*, *) '  2. Mass analysis (jackknife 10-25)'
	write(*, *) '  3. Massp analysis (naive error)'
	write(*, *) '  4. Eta-mass analysis'
	write(*, *) '  5. Eta-mass(2) analysis'
	write(*, *) '  6. PCAC analysis'
	write(*, *) '  7. Effective mass jack-analysis'
        write(*, *) '  8. p2f analysis'
	write(*, *) '  9. Effective mass .masp JK analysis'
	write(*, *) ' 99. Quit'
	read(*, *) ichoice

	if(ichoice .eq. 99) stop

	if(ichoice .eq. 1) then
	  call mass_analysis(ntime, nspace, mcomp, nmeas, 0)
	else if(ichoice .eq. 2) then
	  call mass_analysis(ntime, nspace, mcomp, nmeas, 1)
	else if(ichoice .eq. 3) then
	  call massp_analysis(ncomp, nmeas, ntime, nspace)
	else if(ichoice .eq. 4) then
	  call etamass_analysis(ncomp, nmeas, ntime, nspace)
	else if(ichoice .eq. 5) then
	  call etamass2_analysis(ncomp, nmeas, ntime, nspace)
	else if(ichoice .eq. 6) then
	  call pcac_analysis(ncomp, nmeas, ntime)
	else if(ichoice .eq. 7) then
	  call meffjack(ntime, mcomp, nmeas)
	else if(ichoice .eq. 8) then
	  call p2f_analysis(ntime, nspace, mcomp, nmeas)
	else if(ichoice .eq. 9) then
	  write(*, *) 'nmeas = ', nmeas
	  write(*, *) 'how many of them to analyse?'
	  read(*, *) nmeas
	  call masp_eff_jk(ntime, nspace, nmeas)
	end if

	call load_tail(1)
        close(1)

	end


C     ******************************************************************
      subroutine masp_eff_jk(ntime, nspace, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** MASP_EFF_JK       **   I. Hip, 18 Aug 98                     **
C     ** v3                **   Last modified: 18 Aug 98              **
C     **                   **                                         **
C     ******************************************************************
C     IN ntime - lattice size in time dimension
C     IN nspace - lattice size in space dimension
C     IN nmeas - number of measurements
C     ******************************************************************

	parameter(max_np=MAX_NSPACE / 2)

	real*8 dsp(MAX_NTIME, 0:max_np, 4)
	real*8 conn(0:max_np, 2)

	real*8 trip(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 vac(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 con(MAX_NMEAS, 0:max_np)

	real*8 t(MAX_NTIME - 2, 0:max_np)
	real*8 t_var(MAX_NTIME - 2, 0:max_np)
	real*8 s(MAX_NTIME - 2, 0:max_np)
	real*8 s_var(MAX_NTIME - 2, 0:max_np)
	real*8 wm_t(0:max_np), wm_t_var(0:max_np)
	real*8 wm_s(0:max_np), wm_s_var(0:max_np)

	real*8 rtol, detr, edetr, pbp, pbg5p, uudd, ug5udg5d

	real*8 det, det2, detsum, det2sum

	real*4 tcpu

c	>>> read some user parameter
	write(*, '(1x, ''Top. charge (0-selected, 1-all): '', $)')
	read(*, *) iall
	if(iall .eq. 0) then
	  write(*, '(1x, ''  Which |nu| to collect: '', $)')
	  read(*, *) inu
	end if
	write(*, *)

	write(*, *)
	write(*, *) '0 quenched'
	write(*, *) '1 1-flavour'
	write(*, *) '2 2-flavours'
	read(*, *) iq 

	write(*, *) 'Reading data...'

	detsum = 0.0d0
	det2sum = 0.0d0

	do i = 1, nmeas
	  read(1) nu, nuf, rtol
	  read(1) detr, edetr, pbp, pbg5p, uudd, ug5udg5d

	  read(1)      
     &      (((dsp(j, ip, k), j = 1, ntime - 1),
     &      ip = 0, nspace / 2), k = 1, 4)
	  read(1) ((conn(ip, k), ip = 0, nspace / 2), k = 1, 2) 

	  det = edetr
	  det2 = det**2
	  detsum = detsum + det
	  det2sum = det2sum + det2

	  if(iq .eq. 0) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      trip(i, it, ip) = dsp(it, ip, 1)
	      vac(i, it, ip) = dsp(it, ip, 2)
	      con(i, ip) = conn(ip, 1)
	    end do
	  end do
	  else if(iq .eq. 1) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      trip(i, it, ip) = dsp(it, ip, 1) * det
	      vac(i, it, ip) = dsp(it, ip, 2) * det
	      con(i, ip) = conn(ip, 1) * det
	    end do
	  end do
	  else if(iq .eq. 2) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      trip(i, it, ip) = dsp(it, ip, 1) * det2
	      vac(i, it, ip) = dsp(it, ip, 2) * det2
	      con(i, ip) = conn(ip, 1) * det2
	    end do
	  end do
	  else
	    stop 'unallowed iq'
	  end if  
	end do

	if(iq .eq. 1) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      do i = 1, nmeas
	        trip(i, it, ip) = trip(i, it, ip) / detsum
	        vac(i, it, ip) = vac(i, it, ip) / detsum
	        con(i, ip) = con(i, ip) / detsum
	      end do
	    end do
	  end do
	else if(iq .eq. 2) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      do i = 1, nmeas
	        trip(i, it, ip) = trip(i, it, ip) / det2sum
	        vac(i, it, ip) = vac(i, it, ip) / det2sum
	        con(i, ip) = con(i, ip) / det2sum
	      end do
	    end do
	  end do
	end if

	write(*, *) '...done.'
	write(*, *)

	write(*, *) 'mode (1-var, 2-chi^2) = '
	read(*, *) mode

	write(*, *) 'nplat = '
	read(*, *) nplat

	write(*, *) 'jkblocks = '
	read(*, *) jkblocks

	np = nspace / 2

	call mpjack(nmeas, jkblocks, ntime, np, trip, vac, con, 
     &    t, t_var, s, s_var, 
     &    mode, nplat, wm_t, wm_t_var, wm_s, wm_s_var)

	write(*, *)
	do ip = 0, nspace / 2
	  write(*, *) ip, wm_t(ip), dsqrt(wm_t_var(ip))
	end do

	write(*, *)
	do ip = 0, nspace / 2
	  write(*, *) ip, wm_s(ip), dsqrt(wm_s_var(ip))
	end do

	write(*, *)

	write(*, *)
	do ip = 0, nspace / 2
	  write(*, *) ip * 3.14159265358979d0 / 8.0d0,
     &      wm_t(ip), dsqrt(wm_t_var(ip))
	end do

	write(*, *)
	do ip = 0, nspace / 2
	  write(*, *) ip * 3.14159265358979d0 / 8.0d0,
     &      wm_s(ip), dsqrt(wm_s_var(ip))
	end do

	write(*, *)
        return
        end


C     ******************************************************************
      subroutine meffjack(ntime, mcomp, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** MEFFJACK          **   I. Hip, 17 Jun 97                     **
C     ** v3                **   Last modified: 09 Sep 97              **
C     **                   **                                         **
C     ******************************************************************
C     IN ntime - lattice size in time dimension
C     IN mcomp - number of mass computables
C     IN nmeas - number of measurements
C     ******************************************************************
	real*4 c(MAX_NMEAS, MAX_NTIME, MAX_MCOMP)
	real*8 t1, t1_var, t3, t3_var, s1, s1_var

	write(*, *) 'Reading data...'

	do i = 1, nmeas
	  read(1) ((c(i, it, ic), it = 1, ntime - 1), ic = 1, mcomp)
	end do
	
	write(*, *) '...done.'
	write(*, *)

	write(*, '(1x, ''mode (1-var, 2-chi^2) = '', $)')
	read(*, *) mode

	write(*, '(1x, ''nplat [>= 3] = '', $)')
	read(*, *) nplat

	write(*, '(1x, ''jkblocks [10..25] = '', $)')
	read(*, *) jkblocks
	write(*, *)

	call mjack(nmeas, jkblocks, ntime, mcomp, c, t1, t1_var,
     & t3, t3_var, s1, s1_var, mode, nplat)

        return
        end


C     ******************************************************************
      subroutine jk_parsumav(nmeas, jkblocks, c, a)
C     ******************************************************************
C     - gives the averages of partial sums of c
C     ******************************************************************
C     IN nmeas - number of measurements
C     IN jkblocks - number of jk-blocks i.e. partial sums
C     IN (real*4) c - array with data
C     OUT (real*8) a - array with averages of partial sums
C     ******************************************************************
	real*4 c(*)
        real*8 a(*)
        real*8 sum_c, sum_c2

	irest = mod(nmeas, jkblocks)
	if(irest .ne. 0) then
          write(*, *) 'jk_parsumav WARNING: last ', irest,
     &      ' measurements will be ignored!'
	end if

C       >>> prepare partial sums
	npart = nmeas / jkblocks
	jkmeas = npart * jkblocks
	sum_c = 0.0d0
	sum_c2 = 0.0d0
	do i = 1, jkblocks
	  a(i) = 0.0d0
	  do j = (i - 1) * npart + 1, i * npart
	    a(i) = a(i) + c(j)
	    sum_c2 = sum_c2 + c(j)**2
	  end do
	  sum_c = sum_c + a(i)
	end do
	
C     >>> compute averages
	do i = 1, jkblocks
	  a(i) = (sum_c - a(i)) / dble(jkmeas - npart)
	end do

        return
        end


C     ******************************************************************
      subroutine jk_final(jkblocks, a, aa, var)
C     ******************************************************************
C     IN jkblocks - number of jk-blocks i.e. partial averages
C     IN (real*8) a(*) - array with partial averages
C     OUT (real*8) aa - arithmetic mean
C     OUT (real*8) var - jackknife variance
C     ******************************************************************
	real*8 a(*), aa, var
	real*8 sum_a, sum_a2, aa2

	sum_a = 0.0d0
	sum_a2 = 0.0d0
	do i = 1, jkblocks
	  sum_a = sum_a + a(i)
	  sum_a2 = sum_a2 + a(i)**2
	end do
	aa = sum_a / dble(jkblocks)
	aa2 = sum_a2 / dble(jkblocks)
	var = dble(jkblocks - 1) * (aa2 - aa**2)
	
	return
	end


C     ******************************************************************
      subroutine p2f_analysis(ntime, nspace, mcomp, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** P2F_ANALYSIS      **   I. Hip, 20 Jan 98                     **
C     ** v3                **   Last modified: 20 Jan 98              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine computes the quantities for which there are som
C       analytical predictions by Christof
C     ******************************************************************
C     IN ntime - lattice size in time dimension
C     IN nspace - lattice size in space dimension
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     ******************************************************************
C       >>> number of jk blocks
        parameter(jk=10)

	real*4 c(MAX_NMEAS, MAX_NTIME, MAX_MCOMP)

	real*8 p0(jk, MAX_NTIME), p3(jk, MAX_NTIME)
	real*8 v0(jk, MAX_NTIME), v3(jk, MAX_NTIME)
	real*8 x0(jk, MAX_NTIME), x3(jk, MAX_NTIME)
	real*8 y0(jk, MAX_NTIME), y3(jk, MAX_NTIME)

	real*8 t0(jk, MAX_NTIME), t3(jk, MAX_NTIME)
	real*8 s0(jk, MAX_NTIME), s3(jk, MAX_NTIME)

	real*8 at0(MAX_NTIME), at3(MAX_NTIME)
	real*8 as0(MAX_NTIME), as3(MAX_NTIME)
	real*8 vt0(MAX_NTIME), vt3(MAX_NTIME)
	real*8 vs0(MAX_NTIME), vs3(MAX_NTIME)
C     ******************************************************************

	write(*, *) 'Reading data...'

	do i = 1, nmeas
	  read(1) ((c(i, it, ic), it = 1, ntime - 1), ic = 1, mcomp)
	end do

	write(*, *) '...done.'
	write(*, *)

c       >>> prepare partial averages
        do it = 1, ntime - 1
	  call jk_parsumav(nmeas, jk, c(1, it,  1), p0(1, it))
	  call jk_parsumav(nmeas, jk, c(1, it,  4), p3(1, it))
	  call jk_parsumav(nmeas, jk, c(1, it,  5), v0(1, it))
	  call jk_parsumav(nmeas, jk, c(1, it,  8), v3(1, it))
	  call jk_parsumav(nmeas, jk, c(1, it,  9), x0(1, it))
	  call jk_parsumav(nmeas, jk, c(1, it, 12), x3(1, it))
	  call jk_parsumav(nmeas, jk, c(1, it, 13), y0(1, it))
	  call jk_parsumav(nmeas, jk, c(1, it, 16), y3(1, it))
        end do

c       >>> compute the necessary values
        do it = 1, ntime - 1
          do i = 1, 10
            t0(i, it) = v0(i, it) - x0(i, it) * y0(i, it)
	    s0(i, it) = t0(i, it) - p0(i, it)
            t3(i, it) = v3(i, it) - x3(i, it) * y3(i, it)
	    s3(i, it) = t3(i, it) - p3(i, it)
          end do
          call jk_final(jk, t0(1, it), at0(it), vt0(it))
          call jk_final(jk, s0(1, it), as0(it), vs0(it))
          call jk_final(jk, t3(1, it), at3(it), vt3(it))
          call jk_final(jk, s3(1, it), as3(it), vs3(it))
    	end do

	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 's0 - scalar current; p2f (7)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, as0(it),
     &      dsqrt(vs0(it)) 
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 's3 - pseudoscalar current; p2f (5)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, as3(it),
     &      dsqrt(vs3(it)) 
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 't0 - flavor mixed scalar density; p2f (12)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, at0(it),
     &      dsqrt(vt0(it))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 't3 - flavor mixed pseudoscalar density; p2f (10)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, at3(it),
     &      dsqrt(vt3(it))
	end do
c	write(*, *)
c     &'----------------------------------------------------------------'
c	write(*, *) 'test1 [p2f (5)]: s0 = -s3  =>  s0 / (-s3)'
c	do it = 1, ntime - 1
c	  write(*, '(1x, i4, e13.5, e11.3)') it, - as0(it) / as3(it)
c	end do
c	write(*, *)
c     &'----------------------------------------------------------------'
c	write(*, *) 'test2 [p2f (10)]: t0 = t3  =>  t0 / t3'
c	do it = 1, ntime - 1
c	  write(*, '(1x, i4, e13.5, e11.3)') it, at0(it) / at3(it)
c	end do
	write(*, *)
     &'----------------------------------------------------------------'

	return
	end


C     ******************************************************************
      subroutine mass_analysis(ntime, nspace, mcomp, nmeas, ijack)
C     ******************************************************************
C     **                   **                                         **
C     ** MASS_ANALYSIS     **   I. Hip, 25 Apr 96, v3: Jan 96         **
C     ** v3                **   Last modified: 09 Apr 97              **
C     **                   **                                         **
C     ******************************************************************        
C     - this routine computes arithmetic means and (naive) standard
C       deviations for all sorts of mass computables
C     ******************************************************************
C     IN ntime - lattice size in time dimension
C     IN nspace - lattice size in space dimension
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     IN ijack - flag: 0 - naive stat. analysis, 1 - jackknife
C     ******************************************************************
	real*8 connect
	real*8 mean_x(MAX_NTIME, MAX_MCOMP) 
	real*8 var_x(MAX_NTIME, MAX_MCOMP)
	real*8 var_s(MAX_NTIME, MAX_MCOMP)
	real*8 triplet(MAX_NTIME, 4), singlet(MAX_NTIME, 4) 
	real*8 tpcac(MAX_NTIME), tpcac_var(MAX_NTIME)
        real*8 trpcac(MAX_NTIME), tallpcac(MAX_NTIME)
	real*8 spapcac(MAX_NTIME), spapcac_var(MAX_NTIME)
	real*8 spbpcac(MAX_NTIME), spbpcac_var(MAX_NTIME)
	real*8 spcpcac(MAX_NTIME), spcpcac_var(MAX_NTIME)
	real*8 qmass, qmass2, vqmass, qmass_var
	real*8 flmix(MAX_NTIME, 4), density(MAX_NTIME, 4) 
	real*8 flmix_var(MAX_NTIME, 4), density_var(MAX_NTIME, 4) 

	real*8 teffm(MAX_NTIME - 2), t3effm(MAX_NTIME - 2)
	real*8 seffm(MAX_NTIME - 2)

        if(ijack .eq. 0) then
          call mass_an_naive(ntime, mcomp, nmeas, mean_x, var_x)
        else
          call mass_an_jack(ntime, mcomp, nmeas, mean_x, var_x)
        end if

C       >>> the sign is opposite to the sign in CBL Notes!
        do ig = 1, 4
	  do it = 1, ntime - 1
            triplet(it, ig) = 2.0d0 * mean_x(it, ig)
	    singlet(it, ig) = 2.0d0 * (mean_x(it, ig) +
     &        2.0d0 * (mean_x(it, 8 + ig) * mean_x(it, 12 + ig) -
     &        mean_x(it, 4 + ig)))
	    var_s(it, ig) = 4.0d0 * (var_x(it, ig) + 4.0d0 * (
     &        mean_x(it, 12 + ig)**2 * var_x(it,  8 + ig) +
     &        mean_x(it,  8 + ig)**2 * var_x(it, 12 + ig) +
     &        var_x(it, 4 + ig)))
C         >>> provjeriti proraccun varijanci !!!
          end do
        end do        
	
C	>>> check of (5), (7) and (12) in 2pf
        do ig = 1, 4
	  do it = 1, ntime - 1
C	    >>> flavor-mixed correlator (12)	
            flmix(it, ig) = mean_x(it, 4 + ig) -
     &        mean_x(it, 8 + ig) * mean_x(it, 12 + ig)
C           >>> density (5)       
	    density(it, ig) = flmix(it, ig) - mean_x(it, ig)
	    flmix_var(it, ig) = var_x(it, 4 + ig) +
     &        mean_x(it, 12 + ig)**2 * var_x(it,  8 + ig) +
     &        mean_x(it,  8 + ig)**2 * var_x(it, 12 + ig)
	    density_var(it, ig) = flmix_var(it, ig) + var_x(it, ig)
          end do
        end do        

        do it = 1, ntime - 1
	  tpcac(it) = mean_x(it, 19) / mean_x(it, 4)
	  tpcac_var(it) = (var_x(it, 19) +
     &      (mean_x(it, 19) / mean_x(it, 4))**2 * var_x(it, 4)) /
     &      mean_x(it, 4)**2
	  trpcac(it) = mean_x(it, 20) / mean_x(it, 4) - 2.0d0
	  tallpcac(it) = tpcac(it) + trpcac(it)
C         >>> fali proraccun varijanci za pcac!!!
	end do        

        do it = 2, ntime - 2
c	  spapcac(it) = mean_x(it, 17) / mean_x(it, 4)

	  spapcac(it) = 2.0d0 * mean_x(it, 17) /
     &      (mean_x(it, 4) + mean_x(it + 1, 4))

	  spapcac_var(it) = (var_x(it, 17) +
     &      (mean_x(it, 17) / mean_x(it, 4))**2 * var_x(it, 4)) /
     &      mean_x(it, 4)**2

c	  spbpcac(it) = mean_x(it, 18) / mean_x(it, 4)

	  spbpcac(it) = 2.0d0 * mean_x(it, 18) /
     &      (mean_x(it, 4) + mean_x(it - 1, 4))

	  spbpcac_var(it) = (var_x(it, 18) +
     &      (mean_x(it, 18) / mean_x(it, 4))**2 * var_x(it, 4)) /
     &      mean_x(it, 4)**2

	  spcpcac(it) = mean_x(it, 20) / mean_x(it, 4)
	  spcpcac_var(it) = (var_x(it, 20) +
     &      (mean_x(it, 20) / mean_x(it, 4))**2 * var_x(it, 4)) /
     &      mean_x(it, 4)**2
	end do        

	qmass = 0.0d0
	qmass2 = 0.0d0
	qmass_var = 0.0d0
c >>> tu eventualno odbaciti josh koju toccku kod vechih reshetki...
	do it = 3, ntime - 3
	  qmass = qmass + tpcac(it)
	  qmass2 = qmass2 + tpcac(it)**2
	  qmass_var = qmass_var + tpcac_var(it)
	end do
	qmass = qmass / dble(ntime - 5)
	qmass2 = qmass2 / dble(ntime - 5)
c >>> ovo josh provjeriti/razmisliti...
c	write(*, *) dsqrt(qmass_var / dble(ntime - 5)),
c     &    dsqrt(qmass2 - qmass**2)
	qmass_var = qmass_var / dble(ntime - 5) +
     &    (qmass2 - qmass**2)

	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'triplet scalar density correlation'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, triplet(it, 1), 
     &      dsqrt(var_x(it, 1))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'triplet j_1 current correlation'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, triplet(it, 2), 
     &      dsqrt(var_x(it, 2))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'triplet j_2 current correlation'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, triplet(it, 3), 
     &      dsqrt(var_x(it, 3))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'triplet pseudoscalar density correlation'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, triplet(it, 4), 
     &      dsqrt(var_x(it, 4))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'singlet scalar density correlation'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, singlet(it, 1), 
     &      dsqrt(var_s(it, 1))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'scalar density; p2f (7)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, density(it, 1), 
     &      dsqrt(density_var(it, 1))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'flavor mixed s density; p2f (12)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, flmix(it, 1), 
     &      dsqrt(flmix_var(it, 1))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'j0 vacuum loops'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it,
     &      4.0d0 * (mean_x(it, 9) * mean_x(it, 13) - mean_x(it, 5)),
     &      4.0d0 * dsqrt(var_x(it, 5))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'singlet j_1 current correlation'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, singlet(it, 2), 
     &      dsqrt(var_s(it, 2))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'j_1 vacuum loops'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it,
     &      4.0d0 * mean_x(it, 6), 4.0d0 * dsqrt(var_x(it, 6))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'singlet j_2 current correlation'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, singlet(it, 3), 
     &      dsqrt(var_s(it, 3))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'singlet pseudoscalar density correlation'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, singlet(it, 4), 
     &      dsqrt(var_s(it, 4))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'j3 vacuum loops'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it,
     &      4.0d0 * (mean_x(it, 12) * mean_x(it, 16) - mean_x(it, 8)),
     &      4.0d0 * dsqrt(var_x(it, 8))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'pseudoscalar density; p2f (5)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, density(it, 4), 
     &      dsqrt(density_var(it, 4))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'flavor mixed ps density; p2f (10)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, flmix(it, 4), 
     &      dsqrt(flmix_var(it, 4))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'raw local currents PCAC'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, mean_x(it, 17), 
     &      dsqrt(var_x(it, 17))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'raw PCAC (\pi J)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, mean_x(it, 19), 
     &      dsqrt(var_x(it, 19))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'local currents PCAC (naiv) a) f(x) - f(x - a)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, spapcac(it), 
     &      dsqrt(spapcac_var(it))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'local currents PCAC (naiv) b) f(x + a) - f(x)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, spbpcac(it), 
     &      dsqrt(spbpcac_var(it))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'local currents PCAC (naiv) c) f(x + a) - f(x - a)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, spcpcac(it), 
     &      dsqrt(spcpcac_var(it))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'triplet PCAC (naiv)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, tpcac(it), 
     &      dsqrt(tpcac_var(it))
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'triplet PCAC (all)'
	do it = 1, ntime - 1
	  write(*, '(1x, i4, e13.5, e11.3)') it, tallpcac(it), 
     &      0.0
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'qmass = ', qmass, '   ', dsqrt(qmass_var)
	write(*, *) 
     &'----------------------------------------------------------------'
	

	call massb(ntime, triplet(1, 2), teffm)
	call massb(ntime, triplet(1, 4), t3effm)
	call massb(ntime, singlet(1, 2), seffm)

	write(*, *) 'effective pion mass (from j1 current)'
        do i = 1, ntime / 2 - 1
          write(*, *) i, teffm(i)
        end do
	write(*, *) 
     &'----------------------------------------------------------------'

	write(*, *) 'effective pion mass (from j3 current)'
        do i = 1, ntime / 2 - 1
          write(*, *) i, t3effm(i)
        end do
	write(*, *) 
     &'----------------------------------------------------------------'

	write(*, *) 'effective eta mass (from j1 current)'
        do i = 1, ntime / 2 - 1
          write(*, *) i, seffm(i)
        end do
	write(*, *) 
     &'----------------------------------------------------------------'

	return
	end


C     ******************************************************************
      subroutine mass_an_naive(ntime, mcomp, nmeas, mean_x, var_x)
C     ******************************************************************
C     **                   **                                         **
C     ** MASS_AN_NAIVE     **   I. Hip, 22 Mar 97                     **
C     ** v3                **   Last modified: 22 Mar 97              **
C     **                   **                                         **
C     ******************************************************************
C     IN ntime - lattice size in time dimension
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     OUT (real*8) mean_x - arithmetic means of the computables
C     OUT (real*8) var_x - naive variances of computables
C     ******************************************************************
	real*8 sum_x(MAX_NTIME, MAX_MCOMP)
	real*8 sum_x2(MAX_NTIME, MAX_MCOMP)
	real*8 mean_x(MAX_NTIME, MAX_MCOMP) 
	real*8 mean_x2(MAX_NTIME, MAX_MCOMP)
	real*8 var_x(MAX_NTIME, MAX_MCOMP)

	real*4 f(MAX_NTIME, MAX_MCOMP)
	real*8 ff

	do it = 1, ntime - 1
	  do ic = 1, mcomp
	    sum_x(it, ic) = 0.0d0
	    sum_x2(it, ic) = 0.0d0
	  end do
	end do

	write(*, *) 'Reading data...'

	do i = 1, nmeas
	  read(1) ((f(it, ic), it = 1, ntime - 1), ic = 1, mcomp)
	  do ic = 1, mcomp
	    do it = 1, ntime - 1
	      ff = dble(f(it, ic))
	      sum_x(it, ic) = sum_x(it, ic) + ff
	      sum_x2(it, ic) = sum_x2(it, ic) + ff**2
	    end do
	  end do
	end do

	write(*, *) '...done.'
	write(*, *)

	do ic = 1, mcomp
	  do it = 1, ntime - 1
	    mean_x(it, ic) = sum_x(it, ic) / dble(nmeas)
	    mean_x2(it, ic) = sum_x2(it, ic) / dble(nmeas)
	    var_x(it, ic) = (mean_x2(it, ic) - mean_x(it, ic)**2) /
     &        dble(nmeas)
	  end do
	end do

        return
        end


C     ******************************************************************
      subroutine mass_an_jack(ntime, mcomp, nmeas, mean_x, var_x)
C     ******************************************************************
C     **                   **                                         **
C     ** MASS_AN_JACK      **   I. Hip, 22 Mar 97                     **
C     ** v3                **   Last modified: 22 Mar 97              **
C     **                   **                                         **
C     ******************************************************************
C     IN ntime - lattice size in time dimension
C     IN ncomp - number of computables
C     IN nmeas - number of measurements
C     OUT (real*8) mean_x - biased arithmetic means of the computables
C     OUT (real*8) var_x - jackknife variances of computables
C     ******************************************************************
	real*8 mean_x(MAX_NTIME, MAX_MCOMP) 
	real*8 var_x(MAX_NTIME, MAX_MCOMP)
	real*4 c(MAX_NMEAS, MAX_NTIME, MAX_MCOMP)
	real*8 ac, ab, var_naiv, var_jack, var_all

	write(*, *) 'Reading data...'

	do i = 1, nmeas
	  read(1) ((c(i, it, ic), it = 1, ntime - 1), ic = 1, mcomp)
	end do

	write(*, *) '...done.'
	write(*, *)

        do ic = 1, mcomp
          do it = 1, ntime - 1
            var_all = 0.0d0
            do j = 0, 15
              call jackknife(nmeas, 10 + j, c(1, it, ic), ac, ab,
     &          var_naiv, var_jack)
              var_all = var_all + var_jack
            end do
            mean_x(it, ic) = ab
            var_x(it, ic) = var_all / dble(16)
          end do
	end do

        return
        end


C     ******************************************************************
      subroutine massp_analysis(ncomp, nmeas, ntime, nspace)
C     ******************************************************************
C     **                   **                                         **
C     ** MASSP_ANALYSIS    **   I. Hip, 03 Nov 96                     **
C     **                   **   Last modified: 31 Aug 97              **
C     ******************************************************************        
C     - this routine computes 
C     ******************************************************************
	parameter(max_np = MAX_NSPACE / 2)

	real*8 dsp(MAX_NTIME, 0:max_np, 4)
	real*8 con(0:max_np, 2)
	real*8 d, c

c	>>> needed for dsp analyse
	real*8 dsum_x(MAX_NTIME, 0:max_np, 4)
	real*8 dsum_x2(MAX_NTIME, 0:max_np, 4)
	real*8 dmean_x(MAX_NTIME, 0:max_np, 4)
	real*8 dmean_x2(MAX_NTIME, 0:max_np, 4)

c	>>> needed for con analyse
	real*8 csum_x(0:max_np, 4)
	real*8 csum_x2(0:max_np, 4)
	real*8 cmean_x(0:max_np, 4)
	real*8 cmean_x2(0:max_np, 4)

	real*8 var(MAX_NTIME, 0:max_np, 4)
	real*8 teffm(MAX_NTIME - 2)

	real*8 singlet1(MAX_NTIME, 0:max_np)
	real*8 singlet3(MAX_NTIME, 0:max_np)

	real*8 rtol, detr, edetr, pbp, pbg5p, uudd, ug5udg5d

	real*8 det, det2, detsum, det2sum

	real*4 tcpu

c	>>> read some user parameter
	write(*, '(1x, ''Top. charge (0-selected, 1-all): '', $)')
	read(*, *) iall
	if(iall .eq. 0) then
	  write(*, '(1x, ''  Which |nu| to collect: '', $)')
	  read(*, *) inu
	end if
	write(*, *)

	write(*, *)
	write(*, *) '0 quenched'
	write(*, *) '1 1-flavour'
	write(*, *) '2 2-flavours'
	read(*, *) iq 

c	>>> initialize everything to zero
	do ip = 0, nspace / 2
	  do k = 1, 4
	    do j = 1, ntime
	      dsum_x(j, ip, k) = 0.0d0
	      dsum_x2(j, ip, k) = 0.0d0
	    end do
	  end do
	  do k = 1, 2
	    csum_x(ip, k) = 0.0d0
	    csum_x2(ip, k) = 0.0d0
	  end do
	end do

	detsum = 0.0d0
	det2sum = 0.0d0

c	>>> main loop
	write(*, *) 'Reading data...'

	kmeas = 0
	do i = 1, nmeas
	  read(1) nu, nuf, rtol
	  read(1) detr, edetr, pbp, pbg5p, uudd, ug5udg5d

	read(1)      
     & (((dsp(j, ip, k), j = 1, ntime - 1),
     & ip = 0, nspace / 2), k = 1, 4)
	read(1) ((con(ip, k), ip = 0, nspace / 2), k = 1, 2) 

c	write(*, *) kmeas
c	write(*, *) (dsp(j, 0, 2), j = 1, 15)
c	write(*, *)

c	write(*, *) kmeas, nu, con(0, 1), con(0, 2)

	if((iall .eq. 1) .or. (iabs(nu) .eq. inu)) then
	  kmeas = kmeas + 1

	  det = edetr
	  det2 = det**2
	  detsum = detsum + det
	  det2sum = det2sum + det2

C      dsp(index)
C       1 : Tr(sigma_1 G_{xy} sigma_1 G_{yx}) - sigma_1 pion
C       2 : Tr(sigma_1 G_{xx}) Tr(sigma_1 G_{yy}) - vacuum fluctuation
C       3 : Tr(sigma_3 G_{xy} sigma_3 G_{yx}) = Tr(G_{yx}^+ G_{yx})
C	4 : Tr(sigma_3 G_{xx}) Tr(sigma_3 G_{yy})
C      con(index)
C	1 : Tr(sigma_1 G_{xx})
C	2 : Tr(sigma_3 G_{xx})

	if(iq .eq. 0) then
	  do ip = 0, nspace / 2
	    do k = 1, 4
	      do j = 1, ntime - 1
		d = dsp(j, ip, k)
		dsum_x(j, ip, k) = dsum_x(j, ip, k) + d
		dsum_x2(j, ip, k) = dsum_x2(j, ip, k) + d**2
	      end do
	    end do
	    do k = 1, 2
	      c = con(ip, k)
	      csum_x(ip, k) = csum_x(ip, k) + c
	      csum_x2(ip, k) = csum_x2(ip, k) + c**2
	    end do
	  end do
	else if(iq .eq. 1) then
	  do ip = 0, nspace / 2
	    do k = 1, 4
	      do j = 1, ntime - 1
		d = dsp(j, ip, k)
		dsum_x(j, ip, k) = dsum_x(j, ip, k) + d * det
		dsum_x2(j, ip, k) = dsum_x2(j, ip, k) + d**2
	      end do
	    end do
	    do k = 1, 2
	      c = con(ip, k)
	      csum_x(ip, k) = csum_x(ip, k) + c * det
	      csum_x2(ip, k) = csum_x2(ip, k) + c**2
	    end do
	  end do
	else if(iq .eq. 2) then
	  do ip = 0, nspace / 2
	    do k = 1, 4
	      do j = 1, ntime - 1
		d = dsp(j, ip, k)
		dsum_x(j, ip, k) = dsum_x(j, ip, k) + d * det2
		dsum_x2(j, ip, k) = dsum_x2(j, ip, k) + d**2
	      end do
	    end do
	    do k = 1, 2
	      c = con(ip, k)
	      csum_x(ip, k) = csum_x(ip, k) + c * det2
	      csum_x2(ip, k) = csum_x2(ip, k) + c**2
	    end do
	  end do
	else
	  stop 'unallowed iq'
	end if

	end if

	end do

	write(*, *) '...done.'
	write(*, *)

	if(iq .eq. 0) then
	do ip = 0, nspace / 2
	  do k = 1, 2
	    cmean_x(ip, k) = csum_x(ip, k) / dble(kmeas)
	    cmean_x2(ip, k) = csum_x2(ip, k) / dble(kmeas)
	  end do
	  do k = 1, 4
	    do j = 1, ntime - 1
	      dmean_x(j, ip, k) = dsum_x(j, ip, k) / dble(kmeas)
	      dmean_x2(j, ip, k) = dsum_x2(j, ip, k) / dble(kmeas)
	      var(j, ip, k) = dsqrt((dmean_x2(j, ip, k) -
     &          dmean_x(j, ip, k)**2) / dble(kmeas))

c	    >>> minus at cmean is because of i^2 in x1 and y1
	    singlet1(j, ip) = 2.0d0 * (dmean_x(j, ip, 1) -
     &        2.0d0 * (cmean_x(ip, 1)**2 + dmean_x(j, ip, 2)))

	    singlet3(j, ip) = 2.0d0 * (dmean_x(j, ip, 3) +
     &        2.0d0 * (cmean_x(ip, 2)**2 - dmean_x(j, ip, 4)))

c	    var_s(it, ig) = 4.0d0 * (var_x(it, ig) + 4.0d0 * (
c     &        mean_x(it, 12 + ig)**2 * var_x(it,  8 + ig) +
c     &        mean_x(it,  8 + ig)**2 * var_x(it, 12 + ig) +
c     &        var_x(it, 4 + ig)))
	    end do
	  end do
	end do
	else if(iq .eq. 1) then
	do ip = 0, nspace / 2
	  do k = 1, 2
	    cmean_x(ip, k) = csum_x(ip, k) / detsum / dble(kmeas)
	    cmean_x2(ip, k) = csum_x2(ip, k) / dble(kmeas)
	  end do
	  do k = 1, 4
	    do j = 1, ntime - 1
	      dmean_x(j, ip, k) = dsum_x(j, ip, k) /
     &          detsum / dble(kmeas)
	      dmean_x2(j, ip, k) = dsum_x2(j, ip, k) / dble(kmeas)
	      var(j, ip, k) = dsqrt((dmean_x2(j, ip, k) -
     &          dmean_x(j, ip, k)**2) / dble(kmeas))

c	    >>> minus at cmean is because of i^2 in x1 and y1
	    singlet1(j, ip) = 2.0d0 * (dmean_x(j, ip, 1) -
     &        2.0d0 * (cmean_x(ip, 1)**2 + dmean_x(j, ip, 2)))

	    singlet3(j, ip) = 2.0d0 * (dmean_x(j, ip, 3) +
     &        2.0d0 * (cmean_x(ip, 2)**2 - dmean_x(j, ip, 4)))

c	    var_s(it, ig) = 4.0d0 * (var_x(it, ig) + 4.0d0 * (
c     &        mean_x(it, 12 + ig)**2 * var_x(it,  8 + ig) +
c     &        mean_x(it,  8 + ig)**2 * var_x(it, 12 + ig) +
c     &        var_x(it, 4 + ig)))
	    end do
	  end do
	end do
	else if(iq .eq. 2) then
	do ip = 0, nspace / 2
	  do k = 1, 2
	    cmean_x(ip, k) = csum_x(ip, k) / det2sum / dble(kmeas)
	    cmean_x2(ip, k) = csum_x2(ip, k) / dble(kmeas)
	  end do
	  do k = 1, 4
	    do j = 1, ntime - 1
	      dmean_x(j, ip, k) = dsum_x(j, ip, k) /
     &          det2sum / dble(kmeas)
	      dmean_x2(j, ip, k) = dsum_x2(j, ip, k) / dble(kmeas)
	      var(j, ip, k) = dsqrt((dmean_x2(j, ip, k) -
     &          dmean_x(j, ip, k)**2) / dble(kmeas))

c	    >>> minus at cmean is because of i^2 in x1 and y1
	    singlet1(j, ip) = 2.0d0 * (dmean_x(j, ip, 1) -
     &        2.0d0 * (cmean_x(ip, 1)**2 + dmean_x(j, ip, 2)))

	    singlet3(j, ip) = 2.0d0 * (dmean_x(j, ip, 3) +
     &        2.0d0 * (cmean_x(ip, 2)**2 - dmean_x(j, ip, 4)))

c	    var_s(it, ig) = 4.0d0 * (var_x(it, ig) + 4.0d0 * (
c     &        mean_x(it, 12 + ig)**2 * var_x(it,  8 + ig) +
c     &        mean_x(it,  8 + ig)**2 * var_x(it, 12 + ig) +
c     &        var_x(it, 4 + ig)))
	    end do
	  end do
	end do
	end if

	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) '\pion with sigma_1 current'
	do ip = 0, nspace / 2
c         write(*, *) ip
	  do j = 1, ntime - 1
	    write(*, '(1x, i4, e13.5, e10.3)') j, dmean_x(j, ip, 1),
     &        var(j, ip, 1)
	  end do
	  write(*, *)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) '\eta with sigma_1 current'
	write(*, *)
	do ip = 0, nspace / 2
c         write(*, *) ip
	  do j = 1, ntime - 1
c	    write(*, '(1x, i4, e13.5, e10.3)') j, singlet1(j, ip),
c     &        cmean_x(ip, 1)
	    write(*, '(1x, i4, e13.5)') j, singlet1(j, ip)
	  end do
	  write(*, *)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) '\pion with sigma_3 current'
	write(*, *)
	do ip = 0, nspace / 2
c         write(*, *) ip
	  do j = 1, ntime - 1
	    write(*, '(1x, i4, e13.5, e10.3)') j, dmean_x(j, ip, 3),
     &        var(j, ip, 2)
	  end do
	  write(*, *)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) '\eta with sigma_3 current'
	write(*, *)
	do ip = 0, nspace / 2
c         write(*, *) ip
	  do j = 1, ntime - 1
	    write(*, '(1x, i4, e13.5, e10.3)') j, singlet3(j, ip),
     &        cmean_x(ip, 2)
	  end do
	  write(*, *)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) '\pion with sigma_1 current'
	write(*, *)
	write(*, *) 'nspace / 2 = ', nspace / 2
	write(*, *) 'upper bound: '
	read(*, *) mip
c	do ip = 0, nspace / 2
	do ip = 0, mip
c	  write(*, *) 'ip = ', ip
	  call massb(ntime, dmean_x(1, ip, 1), teffm)
          do i = 1, ntime / 2 - 1
            write(*, *) i, teffm(i)
          end do
	  write(*, *)
	end do
	write(*, *) 
     &'----------------------------------------------------------------'
	write(*, *) '\eta with sigma_1 current'
	write(*, *)
c	do ip = 0, nspace / 2
	do ip = 0, mip
c	  write(*, *) 'ip = ', ip
	  call massb(ntime, singlet1(1, ip), teffm)
          do i = 1, ntime / 2 - 1
            write(*, *) i, teffm(i)
          end do
	  write(*, *)
	end do
	write(*, *) 
     &'----------------------------------------------------------------'
	write(*, *) 'kmeas = ', kmeas, ' out of ', nmeas, ' ( ',
     &    dble(kmeas * 100) / dble(nmeas), ' % )'
	write(*, *) 
     &'----------------------------------------------------------------'

	return
	end


C     ******************************************************************
      subroutine etamass_analysis(ncomp, nmeas, ntime, nspace)
C     ******************************************************************
C     **                   **                                         **
C     ** ETAMASS_ANALYSIS  **   I. Hip, 04 Nov 96                     **
C     **                   **   Last modified: 04 Nov 96              **
C     ******************************************************************        
C     - this routine computes 
C     ******************************************************************
	parameter(max_ncomp=144)

 	real*4 computable(max_ncomp), tcpu
	real*8 com, sum_x(4, MAX_NTIME), sum_x2(4, MAX_NTIME)
	real*8 mean_x(4, MAX_NTIME), mean_x2(4, MAX_NTIME)
	real*8 var(4, MAX_NTIME)
	real*8 connected(MAX_NTIME), vacuum(MAX_NTIME), eta(MAX_NTIME)

	ns = 4 * (ntime - 1) 

	do k = 1, 4
	  do j = 1, ntime
	    sum_x(k, j) = 0.0d0
	    sum_x2(k, j) = 0.0d0
	  end do
	end do

	write(*, *) 'Reading data...'

c       write(*, *) nmeas, ns

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ns)
c       write(*, *) 'i = ', i
c         >>> k = 1 : absdg - absnd (scalar density)
c         >>> k = 2 : mixdg - mixnd (j_1 current)
c         >>> k = 3 : mixdg + mixnd (j_2 current)
c         >>> k = 4 : absdg + absnd (pseudoscalar density)
	  do k = 1, 4
	    do j = 1, ntime - 1
	      com = computable(j + (k - 1) * (ntime - 1))
	      sum_x(k, j) = sum_x(k, j) + com
	      sum_x2(k, j) = sum_x2(k, j) + com**2
	    end do
	  end do
	end do

	close(1)
	write(*, *) '...done.'
	write(*, *)

	do k = 1, 4
	  do j = 1, ntime - 1
	    mean_x(k, j) = sum_x(k, j) / dble(nmeas)
	    mean_x2(k, j) = sum_x2(k, j) / dble(nmeas)
	    var(k, j) = dsqrt((mean_x2(k, j) - mean_x(k, j)**2) /
     &        dble(nmeas))
	  end do
	end do
	
	do j = 1, ntime - 1
	  connected(j) = (mean_x(3, j) * mean_x(4, j)) / dble(nspace**2)
	  vacuum(j) = 2.0d0 * (mean_x(2, j) - connected(j))
	  eta(j) = mean_x(1, j) - vacuum(j) 
	end do

	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'pion'
	do j = 1, ntime - 1
c         write(*, '(1x, i4, e13.5, e10.3)') j, mean_x(1, j), var(1, j)
	  write(*, '(1x, i4, e13.5)') j, mean_x(1, j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'vacuum'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, vacuum(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) '"connected"'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, connected(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'eta'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, eta(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	return
	end


C     ******************************************************************
      subroutine etamass2_analysis(ncomp, nmeas, ntime, nspace)
C     ******************************************************************
C     **                   **                                         **
C     ** ETAMASS2_ANALYSIS **   I. Hip, 19 Nov 96                     **
C     **                   **   Last modified: 19 Nov 96              **
C     ******************************************************************        
C     - this routine computes 
C     ******************************************************************
	parameter(max_ncomp=144)
	real*4 computable(max_ncomp), tcpu
	real*4 computable2(max_ncomp)
	real*8 com, sum_x(4, MAX_NTIME), sum_x2(4, MAX_NTIME)
	real*8 mean_x(4, MAX_NTIME), mean_x2(4, MAX_NTIME)
	real*8 var(4, MAX_NTIME)
	real*8 connected(MAX_NTIME), vacuum(MAX_NTIME), eta(MAX_NTIME)
	real*8 com2, sum2_x(4, MAX_NTIME), sum2_x2(4, MAX_NTIME)
	real*8 mean2_x(4, MAX_NTIME), mean2_x2(4, MAX_NTIME)
	real*8 var2(4, MAX_NTIME)
	real*8 connected2(MAX_NTIME), vacuum2(MAX_NTIME), eta2(MAX_NTIME)

	ns = 4 * (ntime - 1) 

	do k = 1, 4
	  do j = 1, ntime
	    sum_x(k, j) = 0.0d0
	    sum_x2(k, j) = 0.0d0
	    sum2_x(k, j) = 0.0d0
	    sum2_x2(k, j) = 0.0d0
	  end do
	end do

	write(*, *) 'Reading data...'

c       write(*, *) nmeas, ns

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, 2*ns)
c         read(1) (computable2(j), j = 1, ns)
	  do j = 1, ns
	    computable2(j) = computable(j + ns)
	  end do
c       write(*, *) 'i = ', i
	  do k = 1, 4
	    do j = 1, ntime - 1
	      com = computable(j + (k - 1) * (ntime - 1))
	      sum_x(k, j) = sum_x(k, j) + com
	      sum_x2(k, j) = sum_x2(k, j) + com**2
	      com2 = computable2(j + (k - 1) * (ntime - 1))
	      sum2_x(k, j) = sum2_x(k, j) + com2
	      sum2_x2(k, j) = sum2_x2(k, j) + com2**2
	    end do
	  end do
	end do

	close(1)
	write(*, *) '...done.'
	write(*, *)

	do k = 1, 4
	  do j = 1, ntime - 1
	    mean_x(k, j) = sum_x(k, j) / dble(nmeas)
	    mean_x2(k, j) = sum_x2(k, j) / dble(nmeas)
	    var(k, j) = dsqrt((mean_x2(k, j) - mean_x(k, j)**2) /
     &        dble(nmeas))
	    mean2_x(k, j) = sum2_x(k, j) / dble(nmeas)
	    mean2_x2(k, j) = sum2_x2(k, j) / dble(nmeas)
	    var2(k, j) = dsqrt((mean2_x2(k, j) - mean2_x(k, j)**2) /
     &        dble(nmeas))
	  end do
	end do
	
	do j = 1, ntime - 1
	  connected(j) = (mean_x(3, j) * mean_x(4, j)) / dble(nspace**2)
	  vacuum(j) = 2.0d0 * (mean_x(2, j) - connected(j))
	  eta(j) = mean_x(1, j) - vacuum(j) 
	  connected2(j) = (mean2_x(3, j) * mean2_x(4, j)) / 
     &      dble(nspace**2)
	  vacuum2(j) = 2.0d0 * (mean2_x(2, j) - connected2(j))
	  eta2(j) = mean2_x(1, j) - vacuum2(j) 
	end do

	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'pion'
	do j = 1, ntime - 1
c         write(*, '(1x, i4, e13.5, e10.3)') j, mean_x(1, j), var(1, j)
	  write(*, '(1x, i4, e13.5)') j, mean_x(1, j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'vacuum'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, vacuum(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) '"connected"'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, connected(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'eta'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, eta(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'

	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'pion'
	do j = 1, ntime - 1
c         write(*, '(1x, i4, e13.5, e10.3)') j, mean_x(1, j), var(1, j)
	  write(*, '(1x, i4, e13.5)') j, mean2_x(1, j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'vacuum'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, vacuum2(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) '"connected"'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, connected2(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'eta'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, eta2(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	return
	end

C     ******************************************************************
      subroutine pcac_analysis(ncomp, nmeas, ntime)
C     ******************************************************************
C     **                   **                                         **
C     ** PCAC_ANALYSIS     **   I. Hip, 05 Nov 96                     **
C     **                   **   Last modified: 05 Nov 96              **
C     ******************************************************************        
C     - this routine computes 
C     ******************************************************************
c       >>> max_ncomp = maximal number of computables in .data file
	parameter(max_ncomp=144)
	real*4 computable(max_ncomp), tcpu
	real*8 com, sum_x(4, MAX_NTIME), sum_x2(4, MAX_NTIME)
	real*8 mean_x(4, MAX_NTIME), mean_x2(4, MAX_NTIME)
	real*8 var(4, MAX_NTIME)
	real*8 naiv(MAX_NTIME), correctur(MAX_NTIME), all(MAX_NTIME)
	real*8 qmass, qmass2, vqmass

	ns = 4 * (ntime - 1) 

	do k = 1, 4
	  do j = 1, ntime
	    sum_x(k, j) = 0.0d0
	    sum_x2(k, j) = 0.0d0
	  end do
	end do

	write(*, *) 'Reading data...'

c       write(*, *) nmeas, ns

	do i = 1, nmeas
	  read(1) (computable(j), j = 1, ns)
c       write(*, *) 'i = ', i
c         >>> k = 1 : absdg - absnd (scalar density)
c         >>> k = 2 : mixdg - mixnd (j_1 current)
c         >>> k = 3 : mixdg + mixnd (j_2 current)
c         >>> k = 4 : absdg + absnd (pseudoscalar density)
	  do k = 1, 4
	    do j = 1, ntime - 1
	      com = computable(j + (k - 1) * (ntime - 1))
	      sum_x(k, j) = sum_x(k, j) + com
	      sum_x2(k, j) = sum_x2(k, j) + com**2
	    end do
	  end do
	end do

	close(1)
	write(*, *) '...done.'
	write(*, *)

	do k = 1, 4
	  do j = 1, ntime - 1
	    mean_x(k, j) = sum_x(k, j) / dble(nmeas)
	    mean_x2(k, j) = sum_x2(k, j) / dble(nmeas)
	    var(k, j) = dsqrt((mean_x2(k, j) - mean_x(k, j)**2) /
     &        dble(nmeas))
	  end do
	end do

	do j = 1, ntime - 1
	  naiv(j) = mean_x(2, j) / mean_x(4, j)
	  correctur(j) = mean_x(3, j) / mean_x(4, j) - 2.0
	  all(j) = naiv(j) + correctur(j)
	end do
	
	qmass = 0.0d0
	qmass2 = 0.0d0
	do j = 2, ntime - 2
	  qmass = qmass + naiv(j)
	  qmass2 = qmass2 + naiv(j)**2
	end do
	qmass = qmass / dble(ntime - 3)
	qmass2 = qmass2 / dble(ntime - 3)
	vqmass = dsqrt((qmass2 - qmass**2) / dble(ntime - 3))

c       write(*, *)
c     &'----------------------------------------------------------------'
c       write(*, *) 'scalar density (absdg - absnd)'
c       do j = 1, ntime - 1
c         write(*, '(1x, i4, e13.5') j, 0.0
c       end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'naiv PCAC'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, naiv(j)
	end do

	write(*, *)
	write(*, *) 'qmass = ', qmass, '   ', vqmass
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'r corrections'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, correctur(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	write(*, *) 'all together'
	do j = 1, ntime - 1
	  write(*, '(1x, i4, e13.5)') j, all(j)
	end do
	write(*, *)
     &'----------------------------------------------------------------'
	return
	end

