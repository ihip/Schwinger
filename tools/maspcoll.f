C     ******************************************************************
      program maspcoll
C     ******************************************************************
C     **                   **                                         **
C     ** MASPCOLL          **   I. Hip, 2022-04-14                    **
C     ** v6                **   Last modified: 2023-02-02             **
C     **                   **                                         **
C     ******************************************************************
C     >>> new in v6:
C			- j1 masses before fitting are written to file name.Mj1
C     ******************************************************************
C     >>> new in v5:
C			- eta mass, both with sigma1 and sigma3 currents
C     ******************************************************************
C     >>> new in v4:
C	        - analyze subset of configurations
C			- "proper" jackknife for mass error bars
C			- jackknife error for GMOR
C     ******************************************************************

	character*64 masplistname, cname, outname, Mj1name, masp_file_name

	real*8 beta, fmass

	write(*, *)
	write(*, *) 'Masp collector v6 (Hip, 2023-02-02)'
	write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)
	write(*, '(1x, ''.masp list file name: '', $)')
	read(*, '(a)') masplistname
	write(*, *)
	write(*, '(1x, ''name for .col and .Mj1 output files: '', $)')
	read(*, '(a)') cname
	write(*, *)
	outname = cname(1:lnblnk(cname))//'.col'
	Mj1name = cname(1:lnblnk(cname))//'.Mj1'

c	>>> not yet implemented
c	write(*, '(1x, ''Top. charge (0-selected, 1-all): '', $)')
c	read(*, *) iall
c	if(iall .eq. 0) then
c	  write(*, '(1x, ''  Which |nu| to collect: '', $)')
c	  read(*, *) inu
c	end if
c	write(*, *)

	write(*, '(1x, ''number of flavors: '', $)')
	read(*, *) nf
	write(*, *)

c   >>> read fit parameters
	write(*, '(1x, ''mode (1-var, 2-chi^2) = '', $)')
	read(*, *) mode

	write(*, '(1x, ''nplat [>= 3] = '', $)')
	read(*, *) nplat

	write(*, '(1x, ''jkblocks [10..25] = '', $)')
	read(*, *) jkblocks
	write(*, *)

c	>>> read start configuration and the number of measurements
	write(*, '(1x, ''start with conf. No. = '', $)')
	read(*, *) kcfstart

	write(*, '(1x, ''nmeas = '', $)')
	read(*, *) nmeas
	write(*, *)

c	>>> open file with the list of .masp input files
	open(3, file = masplistname, form = 'formatted', status = 'old')

c	>>> create temporary output files for pion and eta
	open(14, file = 'p1.tmp', form = 'formatted', status = 'unknown')
	write(14, '(''# triplet j1 / nf'', i2, $)') nf
	write(14,
     & '('' / mode'', i2, '' / nplat'', i3, '' / jkb'', i4, $)')
     & mode, nplat, jkblocks
	write(14,
     & '('' / start'', i6, '' / nmeas'', i6)')
     & kcfstart, nmeas
	open(15, file = 'p3.tmp', form = 'formatted', status = 'unknown')
	write(15, '(''# triplet j3 / nf'', i2, $)') nf
	write(15,
     & '('' / mode'', i2, '' / nplat'', i3, '' / jkb'', i4, $)')
     & mode, nplat, jkblocks
	write(15,
     & '('' / start'', i6, '' / nmeas'', i6)')
     & kcfstart, nmeas
	open(16, file = 'e1.tmp', form = 'formatted', status = 'unknown')
	write(16, '(''# singlet j1 / nf'', i2, $)') nf
	write(16,
     & '('' / mode'', i2, '' / nplat'', i3, '' / jkb'', i4, $)')
     & mode, nplat, jkblocks
	write(16,
     & '('' / start'', i6, '' / nmeas'', i6)')
     & kcfstart, nmeas
	open(17, file = 'e3.tmp', form = 'formatted', status = 'unknown')
	write(17, '(''# singlet j3 / nf'', i2, $)') nf
	write(17,
     & '('' / mode'', i2, '' / nplat'', i3, '' / jkb'', i4, $)')
     & mode, nplat, jkblocks
	write(17,
     & '('' / start'', i6, '' / nmeas'', i6)')
     & kcfstart, nmeas

c	>>> create output file for Mj1
	open(18, file = Mj1name, form = 'formatted', status = 'unknown')

c   >>> loop over files in the list
	ifile = 0
10    read(3, '(a)', end = 99) masp_file_name

	  open(1, file = masp_file_name, form = 'unformatted',
     &    status = 'old')

	  call load_header(1, nspace, ntime, nconf, beta, fmass)

      call load_mheader(1, mcomp)

	  ifile = ifile + 1
      call masp_eff_jk(ntime, nspace, nconf, kcfstart, nmeas, fmass,	   
     &  nf, mode, nplat, jkblocks)

	  call load_tail(1)
	  close(1)
	  goto 10

99      close(3)

c	>>> close output file p1.tmp
	write(14, *)
	write(14, *)
	close(14)

c	>>> close output file p3.tmp
	write(15, *)
	write(15, *)
	close(15)

c	>>> close output file e1.tmp
	write(16, *)
	write(16, *)
	close(16)

c	>>> close output file e3.tmp
	write(17, *)
	write(17, *)
	close(17)

c	>>> close output file .Mj1
	close(18)

c	>>> create final output file
	call system('cat p1.tmp p3.tmp e1.tmp e3.tmp > '//outname)
	call system('rm p1.tmp')
	call system('rm p3.tmp')
	call system('rm e1.tmp')
	call system('rm e3.tmp')

	write(*, *) '...done.'
	write(*, *)
	write(*, *) '-> results are written to files: '
	write(*, *) outname
	write(*, *) Mj1name
	write(*, *)
	end


C     ******************************************************************
      subroutine masp_eff_jk(ntime, nspace, nconf, kcfstart, nmeas,
     &   fmass, nf, mode, nplat, jkblocks)
C     ******************************************************************
C     **                   **                                         **
C     ** MASP_EFF_JK       **   I. Hip, 18 Aug 98                     **
C     ** v4                **   Last modified: 2022-06-25             **
C     **                   **                                         **
C     ******************************************************************
C     IN integer*4 ntime - lattice size in time dimension
C     IN integer*4 nspace - lattice size in space dimension
C     IN integer*4 nconf - number of configurations in .masp file
C     IN integer*4 kcfstart - start with this configuration
C     IN integer*4 nmeas - number of measurements
C     IN real*8 fmass - fermion (quark) mass
C     IN integer*4 nf - number of flavors
C     IN integer*4 mode - mass fit parameter (1-var, 2-chi^2)
C     IN integer*4 nplat - mass fit parameter: plateau size
C     IN integer*4 jkblocks - number of jackknife subsamples
C     ******************************************************************
	real*8 fmass

	parameter(max_np=MAX_NSPACE / 2)

	real*8 dsp(MAX_NTIME, 0:max_np, 4)
	real*8 conn(0:max_np, 2)

c	>>> j1 (sigma1) current
	real*8 trip(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 vac(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 con(MAX_NMEAS, 0:max_np)

c	>>> j3 (sigma3) current
	real*8 trip3(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 vac3(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 con3(MAX_NMEAS, 0:max_np)

	real*8 sigma_list(MAX_NMEAS)
	real*8 sigma_av, sigma_var
	real*8 det_list(MAX_NMEAS)

c	>>> j1 (sigma1) current
	real*8 t(MAX_NTIME - 2, 0:max_np)
	real*8 t_var(MAX_NTIME - 2, 0:max_np)
	real*8 s(MAX_NTIME - 2, 0:max_np)
	real*8 s_var(MAX_NTIME - 2, 0:max_np)
	real*8 wm_t(0:max_np), wm_t_var(0:max_np)
	real*8 wm_s(0:max_np), wm_s_var(0:max_np)

c	>>> j3 (sigma3) current
	real*8 t3(MAX_NTIME - 2, 0:max_np)
	real*8 t3_var(MAX_NTIME - 2, 0:max_np)
	real*8 s3(MAX_NTIME - 2, 0:max_np)
	real*8 s3_var(MAX_NTIME - 2, 0:max_np)
	real*8 wm_t3(0:max_np), wm_t3_var(0:max_np)
	real*8 wm_s3(0:max_np), wm_s3_var(0:max_np)

	real*8 rtol, sigma, edetr, pbp, pbg5p, uudd, ug5udg5d
	real*8 det, det_max, det_nf, det_nf_sum, det_nf_av
	real*8 gmor, gmor_var, gmor3, gmor3_var

	real*4 tcpu

c	>>> loop over all configurations in .masp file
	i = 0
	do icf = 1, nconf
c	  >>> read all data for one configuration
	  read(1) nu, nuf, rtol
	  read(1) sigma, edetr, pbp, pbg5p, uudd, ug5udg5d
	  read(1)  
     &  (((dsp(j, ip, k), j = 1, ntime - 1),
     &  ip = 0, nspace / 2), k = 1, 4)
	  read(1) ((conn(ip, k), ip = 0, nspace / 2), k = 1, 2) 

c	  >>> start with kcfstart and make nmeas measurements
	  if((icf .ge. kcfstart) .and. (i .lt. nmeas)) then
		i = i + 1
	    det = edetr
	    sigma_list(i) = sigma
	    det_list(i) = det

c	    >>> find det_max
	    if(i .eq. 1) then
	      det_max = det
	    else
	      if(det .gt. det_max) det_max = det
	    end if

	    do ip = 0, nspace / 2
	      do it = 1, ntime
c		    >>> j1 (sigma1) current
	        trip(i, it, ip) = dsp(it, ip, 1)
	        vac(i, it, ip) = dsp(it, ip, 2)
	        con(i, ip) = conn(ip, 1)

c		    >>> j3 (sigma3) current
	        trip3(i, it, ip) = dsp(it, ip, 3)
	        vac3(i, it, ip) = dsp(it, ip, 4)
	        con3(i, ip) = conn(ip, 2)
	      end do
	    end do
	  end if
	end do

c	>>> for number of flavors != 0 reweighting is necessary
	if(nf .ne. 0) then

c	  >>> normalize determinant (det_max -> 1) and compute average
	  det_nf_sum = 0.0d0
	  do i = 1, nmeas
	    det_list(i) = det_list(i) / det_max
	  	det_nf = det_list(i)**nf
		det_nf_sum = det_nf_sum + det_nf
	  end do
	  det_nf_av = det_nf_sum / dble(nmeas)

c	  >>> reweighting
      do i = 1, nmeas
	  	det_nf = det_list(i)**nf / det_nf_av
		do ip = 0, nspace / 2
	      do it = 1, ntime
c			>>> j1 (sigma1) current
	    	trip(i, it, ip) = trip(i, it, ip) * det_nf
	    	vac(i, it, ip) = vac(i, it, ip) * det_nf
	    	con(i, ip) = con(i, ip) * det_nf

c			>>> j3 (sigma3) current
	    	trip3(i, it, ip) = trip3(i, it, ip) * det_nf
	    	vac3(i, it, ip) = vac3(i, it, ip) * det_nf
	    	con3(i, ip) = con3(i, ip) * det_nf
		  end do
		end do
	  end do  

	end if

	np = nspace / 2

c	>>> j1 (sigma1) current
	call mpjack(nmeas, jkblocks, ntime, np, trip, vac, con, 
     &  t, t_var, s, s_var, 
     &  mode, nplat, wm_t, wm_t_var, wm_s, wm_s_var)

c	>>> j3 (sigma3) current
	call mpjack(nmeas, jkblocks, ntime, np, trip3, vac3, con3,
     &  t3, t3_var, s3, s3_var, 
     &  mode, nplat, wm_t3, wm_t3_var, wm_s3, wm_s3_var)

c	>>> Sigma (chiral condensate)
	call rw_jack(nmeas, jkblocks, sigma_list, det_list, nf,
     &  sigma_av, sigma_var)

c	>>> Gell-Mann--Oakes--Renner relation (j1 current)
	gmor = dsqrt(2.0d0 * fmass * sigma_av) / wm_t(0)
	gmor_var = fmass * sigma_var / (2.0d0 * sigma_av * wm_t(0)**2) +
     &  2.0d0 * fmass * sigma_av * wm_t_var(0) / wm_t(0)**4

c	>>> Gell-Mann--Oakes--Renner relation (j3 current)
	gmor3 = dsqrt(2.0d0 * fmass * sigma_av) / wm_t3(0)
	gmor3_var = fmass * sigma_var / (2.0d0 * sigma_av * wm_t3(0)**2) +
     &  2.0d0 * fmass * sigma_av * wm_t3_var(0) / wm_t3(0)**4

c	>>> write to output file p1.tmp
	write(14, '(1x, f6.4, $)') fmass
	write(14, *) wm_t(0), dsqrt(wm_t_var(0)),
     &  sigma_av, dsqrt(sigma_var), gmor, dsqrt(gmor_var)

c	>>> write to output file p3.tmp
	write(15, '(1x, f6.4, $)') fmass
	write(15, *) wm_t3(0), dsqrt(wm_t3_var(0)),
     &  sigma_av, dsqrt(sigma_var), gmor3, dsqrt(gmor3_var)

c	>>> write to output file j1.tmp
	write(16, '(1x, f6.4, $)') fmass
	write(16, *) wm_t(0), dsqrt(wm_t_var(0)),
     &  wm_s(0), dsqrt(wm_s_var(0))

c	>>> write to output file j3.tmp
	write(17, '(1x, f6.4, $)') fmass
	write(17, *) wm_t3(0), dsqrt(wm_t3_var(0)),
     &  wm_s3(0), dsqrt(wm_s3_var(0))

c	>>> write to .Mj1 file
	write(18, '(A, I3)') '# j1 / nplat = ', nplat
	write(18, *) fmass, wm_t(0), dsqrt(wm_t_var(0))
	do it = 1, ntime / 2 - 1
		write(18, *) it, t(it, 0), dsqrt(t_var(it, 0))
	end do
	write(18, *)

	write(*, *)
	return
	end


C     ******************************************************************
      subroutine rw_jack(nconf, jkblocks, a, det, nf, av, var)
C     ******************************************************************
C     **                   **                                         **
C     ** RW_JACK           **   I. Hip, 2022-04-17                    **
C     **                   **   Last modified: 2022-05-05             **
C     ******************************************************************
C     IN integer*4 nconf - number of configurations
C     IN integer*4 jkblocks - number of jackknife subsamples
C     IN real*8 a(nconf) - value of a for each configuration
C     IN real*8 det(nconf) - determinant for each configuration
C     IN integer*4 nf - number of flavors
C     OUT real*8 av - average of jkblocks
C     OUT real*8 var - variance (result of jackknife)
C     ******************************************************************
C     jackknife for reweighted variable a
C     ******************************************************************
	real*8 a(nconf), det(nconf), av, var

	parameter(max_jkblocks=100)

	real*8 weight, sum, wsum 
	real*8 e(max_jkblocks), sum_e, sum_e2
	real*8 av2

	if(jkblocks .gt. max_jkblocks) then
	  jkblocks = max_jkblocks
	  write(*, *) 'Jackknife WARNING: number of blocks to big, ',
     &      'set to ', max_jkblocks, '.'
	end if

	if(mod(nconf, jkblocks) .ne. 0) then
	  write(*, *) 'mod(nconf, jkblocks) != 0'
	  write(*, *) 'Jackknife ERROR: not all conf. included'
	  stop
	end if

C   >>> prepare jkblocks different ensembles
	npart = nconf / jkblocks

	do i = 1, jkblocks
	  sum = 0.0d0
	  wsum = 0.0d0
	  do j = 1, nconf
	    if((j .lt. (i - 1) * npart + 1) .or. (j .gt. i * npart)) then
		  weight = det(j)**nf
		  sum = sum + a(j) * weight
		  wsum = wsum + weight
		end if
	  end do
	  e(i) = sum / wsum
	end do

	sum_e = 0.0d0
	sum_e2 = 0.0d0
	do i = 1, jkblocks
	  sum_e = sum_e + e(i)
	  sum_e2 = sum_e2 + e(i)**2
	end do
c	>>> biased average ???
	av = sum_e / dble(jkblocks)
	av2 = sum_e2 / dble(jkblocks)
	var = dble(jkblocks - 1) * (av2 - av**2)

c	>>> straightforward reweighted average
	sum = 0.0d0
	wsum = 0.0d0
	do j = 1, nconf
	  weight = det(j)**nf
	  sum = sum + a(j) * weight
	  wsum = wsum + weight
	end do
	av = sum / wsum

	return
	end

