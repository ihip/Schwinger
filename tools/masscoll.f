C     ******************************************************************
      program masscoll
C     ******************************************************************
C     **                   **                                         **
C     ** MASSCOLL          **   I. Hip, 09 Apr 97                     **
C     ** v4                **   Last modified: 25 Aug 20              **
C     **                   **                                         **
C     ******************************************************************
C     ...
C     ******************************************************************
	parameter(max_file=100)

C     >>> needed to read .mass file
	real*8 beta, akap
	integer*4 nspace, ntime, nmeas, ncomp
	character*64 c_computable(MAX_MCOMP)

C     >>> other variables which are needed
	character*1 cjack
	character*64 mass_file_name
	character*5 check
	character*32 filelistname, outname

	real*8 t1(max_file), t1_var(max_file)
	real*8 t3(max_file), t3_var(max_file)
	real*8 s1(max_file), s1_var(max_file)
	real*8 beta_f(max_file), akap_f(max_file)
	integer*4 nspace_f(max_file), ntime_f(max_file)

	write(*, *)
	write(*, *) 'Mass collector v4 (ihip, 30 Oct 20)'
	write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)
	write(*, '(1x, ''File list file name: '', $)')
	read(*, '(a)') filelistname
	write(*, '(1x, ''Output file name (e.g. mass.col): '', $)')
	read(*, '(a)') outname
	write(*, *)

c	write(*, *) 'Output file name: mass.col'
c	outname = 'mass.col'

	open(3, file = filelistname, form = 'formatted', status = 'old')
	open(2, file = outname, form = 'formatted', status = 'unknown')

c   >>> read fit parameters
	write(*, '(1x, ''mode (1-var, 2-chi^2) = '', $)')
	read(*, *) mode

	write(*, '(1x, ''nplat [>= 3] = '', $)')
	read(*, *) nplat

	write(*, '(1x, ''jkblocks [10..25] = '', $)')
	read(*, *) jkblocks
	write(*, *)

c   >>> loop over files in the list
	kf = 0
10    read(3, '(a)', end = 99) mass_file_name

	  open(1, file = mass_file_name, form = 'unformatted',
     &    status = 'old')

	  call load_header(1, nspace, ntime, nmeas, beta, akap)

      call load_mheader(1, mcomp)
	  if(mcomp .gt. MAX_MCOMP) then
	    stop 'mcomp is greater than MAX_MCOMP'
	  end if

	  kf = kf + 1
	  nspace_f(kf) = nspace
	  ntime_f(kf) = ntime
	  beta_f(kf) = beta
	  akap_f(kf) = akap
	  call meffjack_coll(ntime, jkblocks, mcomp, nmeas, t1(kf),
     & t1_var(kf), t3(kf), t3_var(kf), s1(kf), s1_var(kf), mode, nplat)

	  call load_tail(1)
	  close(1)
	  goto 10

99      close(3)

	write(2, '(a)') '# pion mass (\sigma_1)'
	write(2, '(''"'', i1, ''"'')') nplat
	do i = 1, kf
	  write(2, '(1x, 2i4, f6.2, f8.4, 2e24.16, i7)') 
     & nspace_f(i), ntime_f(i), beta_f(i), akap_f(i),
     & t1(i), dsqrt(t1_var(i)), nmeas
	end do
	write(2, *)
	write(2, *)

	write(2, '(a)') '# pseudoscalar mass (\sigma_3)'
	write(2, '(''"'', i1, ''"'')') nplat
	do i = 1, kf
	  write(2, '(1x, 2i4, f6.2, f8.4, 2e24.16, i7)') 
     & nspace_f(i), ntime_f(i), beta_f(i), akap_f(i),
     & t3(i), dsqrt(t3_var(i)), nmeas
	end do
	write(2, *)
	write(2, *)

	write(2, '(a)') '# eta mass (\sigma_1)'
	write(2, '(''"'', i1, ''"'')') nplat
	do i = 1, kf
	  write(2, '(1x, 2i4, f6.2, f8.4, 2e24.16, i7)') 
     & nspace_f(i), ntime_f(i), beta_f(i), akap_f(i),
     & s1(i), dsqrt(s1_var(i)), nmeas
	end do

	close(2)
	write(*, *) 'ok!'
	end


C     ******************************************************************
      subroutine meffjack_coll(ntime, jkblocks, mcomp, nmeas,
     &  t1, t1_var, t3, t3_var, s1, s1_var, mode, nplat)
C     ******************************************************************
C     **                   **                                         **
C     ** MEFFJACK_COLL     **   I. Hip, 17 Jun 97                     **
C     ** v3                **   Last modified: 17 Sep 97              **
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

	call mjack(nmeas, jkblocks, ntime, mcomp, c, t1, t1_var,
     & t3, t3_var, s1, s1_var, mode, nplat)

	return
	end


C     ******************************************************************
      subroutine pcac_analysis(akap, ntime, nspace, mcomp, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** PCAC_ANALYSIS     **   I. Hip, 08 Apr 97                     **
C     ** v3                **   Last modified: 08 Apr 97              **
C     **                   **                                         **
C     ******************************************************************        
C     ...
C     ******************************************************************
C     IN ntime - lattice size in time dimension
C     IN nspace - lattice size in space dimension
C     IN mcomp - number of mcomputables
C     IN nmeas - number of measurements
C     ******************************************************************
	real*8 akap
	real*8 connect
	real*8 mean_x(MAX_NTIME, MAX_MCOMP) 
	real*8 var_x(MAX_NTIME, MAX_MCOMP)
	real*8 var_s(MAX_NTIME, MAX_MCOMP)
	real*8 triplet(MAX_NTIME, 4), singlet(MAX_NTIME, 4) 
	real*8 tpcac(MAX_NTIME), trpcac(MAX_NTIME), tallpcac(MAX_NTIME)
        real*8 tpcac_var(MAX_NTIME)
	real*8 qmass, qmass2, vqmass, qmass_var

        call mass_an_jack(ntime, mcomp, nmeas, mean_x, var_x)

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
	
        do it = 1, ntime - 1
	  tpcac(it) = mean_x(it, 19) / mean_x(it, 4)
	  tpcac_var(it) = (var_x(it, 19) +
     &      (mean_x(it, 19) / mean_x(it, 4))**2 * var_x(it, 4)) /
     &      mean_x(it, 4)**2
	  trpcac(it) = mean_x(it, 20) / mean_x(it, 4) - 2.0d0
	  tallpcac(it) = tpcac(it) + trpcac(it)
C         >>> fali proraccun varijanci za pcac!!!
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

	do it = 1, ntime - 1
	  write(2, '(1x, i4, e13.5, e11.3)') it, triplet(it, 2), 
     &      dsqrt(var_x(it, 2))
	end do
	write(2, *)

	do it = 1, ntime - 1
	  write(2, '(1x, i4, e13.5, e11.3)') it, singlet(it, 2), 
     &      dsqrt(var_s(it, 2))
	end do
	write(2, *)
	
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

