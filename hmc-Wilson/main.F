C     ******************************************************************
      program hmc2dfu1
C     ******************************************************************
C     **                    **                                        **
C     **  HMC2DFU1          **    Based on program by                 **
C     **                    **    H. Gausterer and C.B. Lang          **
C     **                    **    Modified by I. Hip, Nov 95          **
C     **                    **    Major revisions: Mar, Dec 96 (ihip) **
C     **  v4.2.1            **                                        **
C     **                    **    Last modified: 30 Jun 25 (ihip)     **
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)
C     >>> ncomp, mcomp - number of computables & mass computables
      parameter (ncomp=27,mcomp=20,nmass=mcomp*(NTIME-1))
C     >>> nbufsize - size of data buffer
      parameter (nbufsize=10000)
C     ------------------------------------------------------------------

C     >>> gauge fields
      complex*16 u(nsite, ndim)

C     >>> simulation parameters
      real*8 beta, akap, akap_start, akap_end, akap_inc, eps

C     >>> needed to call update
      integer*4 update
      real*8 s_final, h_delta

C     >>> variables used for some measurements - not necessary for HMC
      real*4 eig(3), eign(3)
      real*8 cumtime
      real*8 detm

C     >>> declarations for functions in inv.f, etc.
      integer*4 inv
      real*8 op_linkc_e 
      complex*16 op_pbp_e, op_pbg5p_e, op_pbg1p_e, op_pbg2p_e
      complex*16 op_ubudbd_e, op_pbppbp_e, op_ubg5udbg5d_e
      complex*16 pbg1p_e, pbg2p_e
      real*8 sgauge

C     >>> file names
      character*32 start_file
      character*27 conf_name, data_name, mass_name
      character*8 c_ident

C     >>> file name preparation
      character*4 c_nspace, c_ntime
      character*10 c_size
      character*6 c_beta, c_kappa

      character*7 c_iconf
      character*12 c_iconf_name

C     >>> declarations needed for .data file
      character*64 c_computable(ncomp)
      real*4 computable(ncomp, nbufsize)

C     >>> declarations needed for .mass file
      character*64 c_source, c_direction, c_mfuture
      character*64 v_source, v_direction
      character*64 c_mcomputable(mcomp)
      real*4 f(nmass, nbufsize)

C     >>> declarations needed for tail
      integer*4 nexcp, nexcp_b, nexcp1, nexcp1_b
      integer*4 itfut(4)
      real*4 tcpu

C     ------------------------------------------------------------------
      common/kappa/ akap
      common/constants/ beta
      common/gauge_fields/ u
      common/excpcomm/ nexcp, nexcp_b, nexcp1, nexcp1_b
C     ------------------------------------------------------------------

C     >>> number of exceptional configurations
      nexcp = 0
      nexcp_b = 0
      nexcp1 = 0
      nexcp1_b = 0

C     >>> interactive parameter input
      write(*, *)
      write(*, *) 'HMC Wilson - v4.2.1 (30 Jun 25)'
      write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*, *)
      write(*, '(1x, ''beta = '', $)')
      read(*, *) beta
      write(*, '(1x, ''start kappa     = '', $)')
      read(*, *) akap_start
      write(*, '(1x, ''end kappa       = '', $)')
      read(*, *) akap_end
      write(*, '(1x, ''increment kappa = '', $)')
      read(*, *) akap_inc
      write(*, '(1x, ''epsilon = '', $)')
      read(*, *) eps
      write(*, '(1x, ''trajectory steps = '', $)')
      read(*, *) nsteps
      write(*, '(1x, ''ntherm = '', $)')
      read(*, *) ntherm
      write(*, '(1x, ''nmeas = '', $)')
      read(*, *) nmeas
      write(*, '(1x, ''measurement step (neg. -> save)  = '', $)')
      read(*, *) mstep
      write(*, '(1x, ''nsave = '', $)')
      read(*, *) nsave
      if(nsave .gt. nbufsize) stop 'Maximum buffer size exceeded'
      write(*, '(1x, ''istart (0=COLD, 1=HOT, 2=FILE) = '', $)')
      read(*,*) istart
      if(istart .eq. 2) then
	write(*, '(1x, ''start file = '', $)')
	read(*, '(a)') start_file
      end if
5     write(*, '(1x, ''1=triplet, 2=all, 3=dispersion meas. = '', $)')
      read(*, *) imode
      if(imode .ne. 2) stop 'Sorry, not yet implemented!'
      if((imode .lt. 1) .or. (imode .gt. 3)) goto 5

      write(*, *)

C     >>> initialize random number generator (seed is wall-time and host
C         dependent) - be sure to initialize it properly: there is no
C         internal check in rcarry anymore
      call cbl_rtim(iseed)
c      write(*, *) 'WARNING: seed is set to 38768509' 
c      iseed = 38768509
      call cbl_rc_seed(iseed)

C     >>> initialize EXPLICIT DEFINITION of gamma matrices (not needed!)
C      call init_gamma()

C     >>> initialize index arrays
      call mk_index()

C     >>> initialization of gauge fields
      call create_u(istart, u, start_file)

C     >>> here starts akap loop: akap = akap_start, akap_end, akap_inc
      akap = akap_start
10    continue

C     >>> prepare the names for configuration and data files
C     >>> - cname = L[nspace]x[ntime]b[beta]k[kappa]
C     >>> - gauge field configuration is saved in cname.conf
C     >>> - results of measurements in cname.data
C     >>> - mass measurements in cname.mass
      write(c_nspace, '(i4)') NSPACE + 1000
      write(c_ntime, '(i4)') NTIME + 1000
      c_size = 'L'//c_nspace(2:4)//'x'//c_ntime(2:4)//'_b'
      ival = int((beta + 100.0000001d0) * 1000.0d0)
      write(c_beta, '(i6)') ival
      ival = int((akap + 1.000000001d0) * 100000.0d0)
      write(c_kappa, '(i6)') ival
      conf_name = c_size//c_beta(2:6)//'_k'//c_kappa(2:6)//'.conf'
      data_name = c_size//c_beta(2:6)//'_k'//c_kappa(2:6)//'.data'
      mass_name = c_size//c_beta(2:6)//'_k'//c_kappa(2:6)//'.mass'

C     >>> prepare textual data description
      c_computable(1) = 'acc'
      c_computable(2) = 'niter'
      c_computable(3) = 'eps'
      c_computable(4) = 'sgauge / nsites'
      c_computable(5) = 'topological charge'
      c_computable(6) = 'exp(-h_delta)'
      c_computable(7) = 'niter1'
      c_computable(8) = 'Re(op_pbp_e)'
      c_computable(9) = 'link condensate'
      c_computable(10) = 'Re(op_pbg5p_e)'
      c_computable(11) = 'Re(op_ubudbd_e)'
      c_computable(12) = 'Re(op_pbppbp_e)'
      c_computable(13) = 'Re(op_pbg1p_e)'
      c_computable(14) = 'Im(op_pbg1p_e)'
      c_computable(15) = 'Re(op_pbg2p_e)'
      c_computable(16) = 'Im(op_pbg2p_e)'
      c_computable(17) = 'smallest eigenvalue (M^+ M) AP'
      c_computable(18) = '2nd smallest eigenvalue'
      c_computable(19) = '3rd smallest eigenvalue'
      c_computable(20) = 'energy (s_final)'
      c_computable(21) = 'cummulative time'
      c_computable(22) = 'num. of real ev-s'
      c_computable(23) = 'smallest abs(ev(M)) AP'
      c_computable(24) = '2nd smallest'
      c_computable(25) = '3rd smallest'
      c_computable(26) = 'det M'
      c_computable(27) = 'Re(op_ubg5udbg5d)'

      open(unit=7,file=conf_name,form='unformatted',status='unknown')
      open(unit=8,file=data_name,form='unformatted',status='unknown')
      if(imode .eq. 2) then
        open(unit=9,file=mass_name,form='unformatted',status='unknown')
      end if

      write(c_ident, '(a8)') '%mass-v3'

C for the first kappa, write the initialization data on the screen
      if(akap .eq. akap_start) then
        nw = 1
      else
        nw = 0
      end if
C save header in .data file
      call save_header(8, nw, c_ident, NSPACE, NTIME, beta,
     &  akap, eps, nsteps, ntherm, nmeas, mstep, istart, iseed)

C save header in .mass file
      call save_header(9, 0, c_ident, NSPACE, NTIME, beta,
     &  akap, eps, nsteps, ntherm, nmeas, mstep, istart, iseed)

C .data specific part - saving description of computables
      write(8) ncomp
      write(8) (c_computable(i), i = 1, ncomp)

C .mass specific part - saving desc. of mass computables
      c_source = v_source()
      c_direction = v_direction()
      write(9) mcomp, c_source, c_direction, c_mfuture
      write(9) (c_mcomputable(i), i = 1, mcomp)

C     ------------------------------------------------------------------
C     Equilibration - thermalization
C     ------------------------------------------------------------------
      do 100 it = 1, ntherm
	niter = update(eps, nsteps, iacc, s_final, h_delta)
	call adjust_eps(eps, iacc, it)
	if(mod(it, 100) .eq. 0) call unitarize(ngaug, u)
	if(mod(it, nsave) .eq. 0) then
	  rewind 7
	  write(7) NTIME, NSPACE
          write(7) beta, akap, eps
          write(7) it
	  write(7) u
	  endfile 7
	end if
100   continue

C     ------------------------------------------------------------------
C     Measuring  updates
C     ------------------------------------------------------------------

c      call initwindow()


      call cbl_cputime(0, tcpu)

      ns = 1
      it = 0
      cumtime = 0.0d0
200   it = it + 1

C       >>> update
	do i = 1, iabs(mstep)
	  ii = (it - 1) * iabs(mstep) + i
	  if(mod(ii, 100) .eq. 0) call unitarize(ngaug, u)
	  niter = update(eps, nsteps, iacc, s_final, h_delta)
          cumtime = cumtime + iacc * nsteps * eps
	  call adjust_eps(eps, iacc, ii)
	end do

C     >>> added in version 4.2 / 2025-05-08 / ihip
C     >>> if mstep is negative save the configuration
      if(mstep .lt. 0) then

        write(c_iconf, '(i7)') it + 1000000
        c_iconf_name = 'i'//c_iconf(2:7)//'.conf'

        open(unit=77, file=c_iconf_name, form='unformatted',
     &status='unknown')
	  write(77) NTIME, NSPACE
          write(77) beta, akap, eps
          write(77) ntherm + it * iabs(mstep)
	  write(77) u
        close(77)
      end if
	
C       >>> measurements
	computable(1, ns) = real(iacc)
	computable(2, ns) = real(niter)
	computable(3, ns) = real(eps)
	computable(4, ns) = real(sgauge(beta, u) / dble(nsite))
	computable(5, ns) = real(itopch())
	computable(6, ns) = real(dexp(-h_delta))

c	call eigen(eig)
c	computable(17, ns) = eig(1)
c	computable(18, ns) = eig(2)
c	computable(19, ns) = eig(3)

        computable(20, ns) = real(s_final / dble(nsite))
        computable(21, ns) = real(cumtime)

C >>> just temporary for test:
c	write(*, *) 'topch = ', computable(5, ns)
c	write(*, *)
c	do i = 1, 3
c	  write(*, *) i, eig(i)
c	end do
c	computable(22, ns) = real(nageig(akap, detm, eign))
c	computable(23, ns) = eign(1)
c	computable(24, ns) = eign(2)
c	computable(25, ns) = eign(3)
c	computable(26, ns) = real(detm)

C >>> end of test

	if(imode .eq. 2) then
        niter1 = inv()
        
C        print *, 'update: ', it, niter, niter1

C         >>> if both inversion routines fail, forget the measurement!
          if(niter1 .lt. -2) then
            it = it - 1
            go to 200
          end if
	  computable(7, ns) = real(niter1)
	  if(niter1 .lt. 0) then
C           >>> in the case of inv failure, fill computables with zeroes
	    do i = 7, 16
	      computable(i, ns) = 0.0
	    end do
	  else
	  computable(8, ns) = real(dreal(op_pbp_e() / dble(nsite)))
	  computable(9, ns) = - real(2 * op_linkc_e() / dble(nsite))
	  computable(10, ns) = real(dreal(op_pbg5p_e() / dble(nsite)))
	  computable(11, ns) = real(dreal(op_ubudbd_e() / dble(nsite)))
	  computable(27, ns) = real(dreal(op_ubg5udbg5d_e() /
     & dble(nsite)))
	  computable(12, ns) = real(dreal(op_pbppbp_e() / dble(nsite)))

	  pbg1p_e = op_pbg1p_e()
	  computable(13, ns) = real(dreal(pbg1p_e) / dble(nsite))
	  computable(14, ns) = real(dimag(pbg1p_e) / dble(nsite))

	  pbg2p_e = op_pbg2p_e()
	  computable(15, ns) = real(dreal(pbg2p_e) / dble(nsite))
	  computable(16, ns) = real(dimag(pbg2p_e) / dble(nsite))
	  end if
	
	  call imass(f(1, ns))
c	  call imass_sp(f(1, ns))
	
	else if(imode .eq. 1) then
	  stop 'Sorry, not yet implemented!'
c         m = mass(...)
	else if(imode .eq. 3) then
	  stop 'Sorry, not yet implemented!'
c         m = imassp(...)
	end if

	if((ns .eq. nsave) .or. (it .eq. nmeas)) then
	  rewind 7
	  write(7) NTIME, NSPACE
          write(7) beta, akap, eps
          write(7) ntherm + it * iabs(mstep)
	  write(7) u
	  endfile 7
	  do i = 1, nsave
	    write(8) (computable(j, i), j = 1, ncomp)
	  end do
c	  endfile 8
	  if(imode .eq. 2) then
	    do i = 1, nsave
	      write(9) (f(j, i), j = 1, nmass)
	    end do
	  end if
c	  endfile 9
	  ns = 0
	end if
	ns = ns + 1
      if(it .lt. nmeas) go to 200

      call cbl_cputime(1, tcpu)

      call save_tail(8, nexcp, nexcp_b, nexcp1, nexcp1_b, tcpu)
      call save_tail(9, nexcp, nexcp_b, nexcp1, nexcp1_b, tcpu)

      close(8)
      if(imode .eq. 2) then
        close(9)
      end if

      akap = akap + akap_inc
      if(akap .gt. akap_end + 0.00000001) goto 999
      goto 10

999   continue

      end

