C     ******************************************************************
      program pcaccoll
C     ******************************************************************
C     **                   **                                         **
C     ** PCACCOLL          **   I. Hip, 08 Apr 97                     **
C     ** v3                **   Last modified: 08 Sep 97              **
C     **                   **                                         **
C     ******************************************************************
C     ...
C     ******************************************************************

C       >>> needed to read .mass file
	real*8 beta, akap
	integer*4 nspace, ntime, nmeas, ncomp
	character*64 c_computable(MAX_MCOMP)

C       >>> other variables which are needed
	character*1 cjack
	character*64 mass_file_name
	character*5 check
	character*32 listfilename, outname

	write(*, *)
	write(*, *) 'PCAC collector v3 (ivh, 08 Sep 97)'
	write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)
	write(*, '(1x, ''File list file name: '', $)')
	read(*, '(a)') listfilename
	write(*, '(1x, ''Output file name (e.g. pcac.col): '', $)')
	read(*, '(a)') outname
	write(*, *)

c	write(*, *) 'Output file name: pcac.col'
c	outname = 'pcac.col'

	open(3, file = listfilename, form = 'formatted', status = 'old')
	open(2, file = outname, form = 'formatted', status = 'unknown')

c     >>> main loop over files in the list
10    read(3, '(a)', end = 99) mass_file_name

	  open(1, file = mass_file_name, form = 'unformatted',
     &    status = 'old')

	  call load_header(1, nspace, ntime, nmeas, beta, akap)

	  call load_mheader(1, mcomp)
	  if(mcomp .gt. MAX_MCOMP) then
	    write(*, *) 'mcomp is greater than MAX_MCOMP'
	    stop
	  end if

	  call pcac_analysis(akap, ntime, nspace, mcomp, nmeas)

	  call load_tail(1)
	  close(1)
	  goto 10

99      close(3)
	close(2)
	write(*, *) 'ok!'
	end


C     ******************************************************************
      subroutine pcac_analysis(akap, ntime, nspace, mcomp, nmeas)
C     ******************************************************************
C     **                   **                                         **
C     ** PCAC_ANALYSIS     **   I. Hip, 08 Apr 97                     **
C     ** v3                **   Last modified: 08 Sep 97              **
C     **                   **                                         **
C     ******************************************************************        
C     ...
C     ******************************************************************
C     IN (real*8) akap - kappa
C     IN ntime - lattice size in time dimension
C     IN nspace - lattice size in space dimension
C     IN mcomp - number of mcomputables
C     IN nmeas - number of measurements
C     ******************************************************************
	real*8 akap
	real*8 mean_x(MAX_NTIME, MAX_MCOMP), var_x(MAX_NTIME, MAX_MCOMP)
	real*8 tpcac(MAX_NTIME), trpcac(MAX_NTIME), tallpcac(MAX_NTIME)
        real*8 tpcac_var(MAX_NTIME)
	real*8 qmass, ss, qmass_var, smm

        call mass_an_jack(ntime, mcomp, nmeas, mean_x, var_x)

        do it = 1, ntime - 1
	  tpcac(it) = mean_x(it, 19) / mean_x(it, 4)
	  tpcac_var(it) = (var_x(it, 19) +
     &      (mean_x(it, 19) / mean_x(it, 4))**2 * var_x(it, 4)) /
     &      mean_x(it, 4)**2
	  trpcac(it) = mean_x(it, 20) / mean_x(it, 4) - 2.0d0
	  tallpcac(it) = tpcac(it) + trpcac(it)
	end do        

c >>> for 8x8 lattice 2, for bigger lattices 3 points from both
c     sides are excluded from the computation of the average
	if(ntime .le. 8) then
	  leave = 2
	else 
          leave = 3
	end if

c >>> computation of the plateau (weighted) average: qmass
	qmass = 0.0d0
	ss = 0.0d0
	do it = 1 + leave, ntime - 1 - leave
	  qmass = qmass + tpcac(it) / tpcac_var(it)
	  ss = ss + 1 / tpcac_var(it)
	end do
	qmass = qmass / ss

c >>> computation of final error (errors + scattering): qmass_var
	qmass_var = 0.0d0
	smm = 0.0d0
	do it = 1 + leave, ntime - 1 - leave
	  qmass_var = qmass_var + tpcac_var(it)
	  smm = smm + (tpcac(it) - qmass)**2
	end do
	qmass_var = (qmass_var + smm) / dble(ntime - 1 - 2 * leave)

	write(2, '(1x, f8.5, 2e24.16, i7)') akap, qmass,
     &    dsqrt(qmass_var), nmeas
	
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

