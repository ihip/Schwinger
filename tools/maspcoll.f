C     ******************************************************************
      program maspcoll
C     ******************************************************************
C     **                   **                                         **
C     ** MASPCOLL          **   I. Hip, 2022-04-14                    **
C     ** v1                **   Last modified: 2022-04-14             **
C     **                   **                                         **
C     ******************************************************************
C     ...
C     ******************************************************************

	character*64 masplistname, outname, masp_file_name

	real*8 beta, fmass

	write(*, *)
	write(*, *) 'Masp collector v1 (Hip, 2022-04-14)'
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

	write(*, *) '0 quenched'
	write(*, *) '2 2-flavours'
	write(*, *) '4 4-flavours'
	read(*, *) iq 
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
      call masp_eff_jk(ntime, nspace, nmeas, fmass, iq, mode, nplat,
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
      subroutine masp_eff_jk(ntime, nspace, nmeas, fmass, iq, mode,
     &   nplat, jkblocks)
C     ******************************************************************
C     **                   **                                         **
C     ** MASP_EFF_JK       **   I. Hip, 18 Aug 98                     **
C     ** v4                **   Last modified: 2022-04-14             **
C     **                   **                                         **
C     ******************************************************************
C     IN ntime - lattice size in time dimension
C     IN nspace - lattice size in space dimension
C     IN nmeas - number of measurements
C     ******************************************************************
	real*8 fmass

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

	real*8 rtol, sigma, edetr, pbp, pbg5p, uudd, ug5udg5d

	real*8 det, det2, det2sum, det4, det4sum

	real*4 tcpu

	detsum = 0.0d0
	det2sum = 0.0d0
	det4sum = 0.0d0

	do i = 1, nmeas
	  read(1) nu, nuf, rtol
	  read(1) sigma, edetr, pbp, pbg5p, uudd, ug5udg5d

	  read(1)      
     &      (((dsp(j, ip, k), j = 1, ntime - 1),
     &      ip = 0, nspace / 2), k = 1, 4)
	  read(1) ((conn(ip, k), ip = 0, nspace / 2), k = 1, 2) 

	  det = edetr
	  det2 = det**2
	  det2sum = det2sum + det2
	  det4 = det**4
	  det4sum = det4sum + det4

	  if(iq .eq. 0) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      trip(i, it, ip) = dsp(it, ip, 1)
	      vac(i, it, ip) = dsp(it, ip, 2)
	      con(i, ip) = conn(ip, 1)
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
	  else if(iq .eq. 4) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      trip(i, it, ip) = dsp(it, ip, 1) * det4
	      vac(i, it, ip) = dsp(it, ip, 2) * det4
	      con(i, ip) = conn(ip, 1) * det4
	    end do
	  end do
	  else
	    stop 'unallowed iq'
	  end if  
	end do

	if(iq .eq. 2) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      do i = 1, nmeas
	        trip(i, it, ip) = trip(i, it, ip) / det2sum
	        vac(i, it, ip) = vac(i, it, ip) / det2sum
	        con(i, ip) = con(i, ip) / det2sum
	      end do
	    end do
	  end do
	else if(iq .eq. 4) then
	  do ip = 0, nspace / 2
	    do it = 1, ntime
	      do i = 1, nmeas
	        trip(i, it, ip) = trip(i, it, ip) / det4sum
	        vac(i, it, ip) = vac(i, it, ip) / det4sum
	        con(i, ip) = con(i, ip) / det4sum
	      end do
	    end do
	  end do
	end if

	write(*, *) '...done.'
	write(*, *)

	np = nspace / 2

	call mpjack(nmeas, jkblocks, ntime, np, trip, vac, con, 
     &    t, t_var, s, s_var, 
     &    mode, nplat, wm_t, wm_t_var, wm_s, wm_s_var)

	write(*, *)

c	>>> write to output file
	write(2, '(1x, f6.4, $)') fmass
	write(2, *) wm_t(0), dsqrt(wm_t_var(0))

	return
	end
