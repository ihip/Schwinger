C     ******************************************************************
      program maspcoll
C     ******************************************************************
C     **                   **                                         **
C     ** MASPCOLL          **   I. Hip, 2022-04-14                    **
C     ** v2                **   Last modified: 2022-05-01             **
C     **                   **                                         **
C     ******************************************************************
C     ...
C     ******************************************************************

	character*64 masplistname, outname, masp_file_name

	real*8 beta, fmass

	write(*, *)
	write(*, *) 'Masp collector v2 (Hip, 2022-05-01)'
	write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*, *)
	write(*, '(1x, ''.masp list file name: '', $)')
	read(*, '(a)') masplistname
	write(*, *)
	write(*, '(1x, ''output file name (e.g. masp.col): '', $)')
	read(*, '(a)') outname
	write(*, *)

	open(3, file = masplistname, form = 'formatted', status = 'old')
	open(2, file = outname, form = 'formatted', status = 'unknown')

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

c   >>> loop over files in the list
	kf = 0
10    read(3, '(a)', end = 99) masp_file_name

	  open(1, file = masp_file_name, form = 'unformatted',
     &    status = 'old')

	  call load_header(1, nspace, ntime, nmeas, beta, fmass)

      call load_mheader(1, mcomp)

	  kf = kf + 1
      call masp_eff_jk(ntime, nspace, nmeas, fmass, nf, mode, nplat,
     &  jkblocks)

	  call load_tail(1)
	  close(1)
	  goto 10

99      close(3)

c	>>> close output file
	write(2, *)
	close(2)

	write(*, *) '...done. -> results are written to file: ', outname
	write(*, *)
	end


C     ******************************************************************
      subroutine masp_eff_jk(ntime, nspace, nmeas, fmass, nf, mode,
     &   nplat, jkblocks)
C     ******************************************************************
C     **                   **                                         **
C     ** MASP_EFF_JK       **   I. Hip, 18 Aug 98                     **
C     ** v4                **   Last modified: 2022-04-20             **
C     **                   **                                         **
C     ******************************************************************
C     IN integer*4 ntime - lattice size in time dimension
C     IN integer*4 nspace - lattice size in space dimension
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

	real*8 trip(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 vac(MAX_NMEAS, MAX_NTIME, 0:max_np)
	real*8 con(MAX_NMEAS, 0:max_np)

	real*8 sigma_list(MAX_NMEAS)
	real*8 sigma_av, sigma_var
	real*8 det_list(MAX_NMEAS)

	real*8 t(MAX_NTIME - 2, 0:max_np)
	real*8 t_var(MAX_NTIME - 2, 0:max_np)
	real*8 s(MAX_NTIME - 2, 0:max_np)
	real*8 s_var(MAX_NTIME - 2, 0:max_np)
	real*8 wm_t(0:max_np), wm_t_var(0:max_np)
	real*8 wm_s(0:max_np), wm_s_var(0:max_np)

	real*8 rtol, sigma, edetr, pbp, pbg5p, uudd, ug5udg5d
	real*8 det, det_max, det_nf, det_nf_sum, det_nf_av
	real*8 gmor, gmor_var

	real*4 tcpu

c	>>> loop over all measurements (configurations)
	do i = 1, nmeas
	  read(1) nu, nuf, rtol
	  read(1) sigma, edetr, pbp, pbg5p, uudd, ug5udg5d

	  det = edetr
	  sigma_list(i) = sigma
	  det_list(i) = det

c	  >>> find det_max
	  if(i .eq. 1) then
	    det_max = det
	  else
	    if(det .gt. det_max) det_max = det
	  end if

	  read(1)  
     &  (((dsp(j, ip, k), j = 1, ntime - 1),
     &  ip = 0, nspace / 2), k = 1, 4)
	  read(1) ((conn(ip, k), ip = 0, nspace / 2), k = 1, 2) 

	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      trip(i, it, ip) = dsp(it, ip, 1)
	      vac(i, it, ip) = dsp(it, ip, 2)
	      con(i, ip) = conn(ip, 1)
	    end do
	  end do

	end do

c	>>> for number of flavors != 0 reweighting is necessary
	if(nf .ne. 0) then

c	  >>> normalize determinant (det_max -> 1)
      do i = 1, nmeas
	    det_list(i) = det_list(i) / det_max
c		write(*, *) i, det_list(i)
	  end do

c	  >>> reweighting
	  det_nf_sum = 0.0d0
      do i = 1, nmeas
	  	det_nf = det_list(i)**nf
		det_nf_sum = det_nf_sum + det_nf
		do ip = 0, nspace / 2
	      do it = 1, ntime
	    	trip(i, it, ip) = trip(i, it, ip) * det_nf
	    	vac(i, it, ip) = vac(i, it, ip) * det_nf
	    	con(i, ip) = con(i, ip) * det_nf
		  end do
		end do
	  end do  

	  det_nf_av = det_nf_sum / dble(nmeas)
	  do i = 1, nmeas
		do ip = 0, nspace / 2
	      do it = 1, ntime
	        trip(i, it, ip) = trip(i, it, ip) / det_nf_av
	        vac(i, it, ip) = vac(i, it, ip) / det_nf_av
	        con(i, ip) = con(i, ip) / det_nf_av
	      end do
	    end do
	  end do
	end if

	np = nspace / 2

	call mpjack(nmeas, jkblocks, ntime, np, trip, vac, con, 
     &  t, t_var, s, s_var, 
     &  mode, nplat, wm_t, wm_t_var, wm_s, wm_s_var)

	call rw_jack(nmeas, jkblocks, sigma_list, det_list, nf,
     &  sigma_av, sigma_var)

c	>>> Gell-Mann--Oakes--Renner relation
	gmor = dsqrt(2.0d0 * fmass * sigma_av) / wm_t(0)
	gmor_var = fmass * sigma_var / (2.0d0 * sigma_av * wm_t(0)**2) +
     &  2.0d0 * fmass * sigma_av * wm_t_var(0) / wm_t(0)**4

c	>>> write to output file
	write(2, '(1x, f6.4, $)') fmass
	write(2, *) wm_t(0), dsqrt(wm_t_var(0)),
     &  sigma_av, dsqrt(sigma_var), gmor, dsqrt(gmor_var)

	write(*, *)
	return
	end


C     ******************************************************************
      subroutine rw_jack(nconf, jkblocks, a, det, nf, av, var)
C     ******************************************************************
C     **                   **                                         **
C     ** RW_JACK           **   I. Hip, 2022-04-17                    **
C     **                   **   Last modified: 2022-04-17             **
C     ******************************************************************
C     IN integer*4 nmeas - number of configurations
C     IN integer*4 jkblocks - number of jackknife subsamples
C     IN real*8 a(nmeas) - value of a for each configuration
C     IN real*8 det(nmeas) - determinant for each configuration
C     IN integer*4 nf - number of flavors
C     OUT real*8 av - average of jkblocks
C     OUT real*8 var - variance (result of jackknife)
C     ******************************************************************
C     jackknife for reweighted variable a
C     ******************************************************************
	real*8 a(nmeas), det(nmeas), av, var

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

