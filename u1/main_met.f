C     ******************************************************************
      program main_metropolis
C     ******************************************************************
C     **                    **                                        **
C     **  METROPOLIS        **    I. Hip, 29 Sep 98                   **
C     **                    **    Last modified: 29 Sep 98 (hip)      **
C     **                    **                                        **
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)

C     >>> ncomp - number of computables
      parameter (ncomp=8)
C     >>> nbufsize - size of data buffer
      parameter (nbufsize=10000)
C     ------------------------------------------------------------------
C     >>> gauge fields
      complex*16 u(nsite, ndim)

C     >>> indexes needed for sgauge & itopch
      integer ind(4, nsite)

C     >>> simulation parameters
      real*8 beta, akap, eps

C     >>> declarations for functions
      real*8 sgauge, sintopch

C     >>> file names handling
      character*32 start_file
      character*6 c_beta, c_kappa
      character*17 conf_name, data_name
      character*8 c_ident

C     >>> declarations needed for .data file
      character*64 c_computable(ncomp)
      real*4 computable(ncomp, nbufsize)

C     >>> declarations needed for tail
      integer*4 nexcp, nexcp_b, nexcp1, nexcp1_b
      integer*4 itfut(4)
      real*4 tcpu

C     ------------------------------------------------------------------
      common/constants/ beta
      common/gauge_fields/ u
      common /index_arrays/ind
C     ------------------------------------------------------------------
C     Initialization
C     ------------------------------------------------------------------
C     >>> pure gauge means \kappa = 0
      akap = 0.0d0

C     >>> interactive parameter input
      write(*, *)
      write(*, *) 'Metropolis (29 Sep 98 / 12 Jul 21)'
      write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*, *)
      write(*, '(1x, ''start beta     = '', $)')
      read(*, *) beta_start
      write(*, '(1x, ''end beta       = '', $)')
      read(*, *) beta_end
      write(*, '(1x, ''increment beta = '', $)')
      read(*, *) beta_inc
      write(*, '(1x, ''ntherm = '', $)')
      read(*, *) ntherm
      write(*, '(1x, ''nmeas = '', $)')
      read(*, *) nmeas
      write(*, '(1x, ''measurement step = '', $)')
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
      write(*, *)

C     >>> initialize random number generator (seed is wall-time and host
C         dependent) - be sure to initialize it properly: there is no
C         internal check in rcarry anymore
      call cbl_rtim(iseed)
C      write(*, *) 'WARNING: seed is set to 38768509' 
C      iseed = 38768509
      call cbl_rc_seed(iseed)

C     >>> initialize gauge index arrays
      call gauge_mk_index

C     >>> initialize indexes for sgauge & itopch
      call mk_index(ind)

C     >>> initialization of gauge fields
      call create_u(istart, u, start_file)

C     >>> here starts beta loop: beta = beta_start, beta_end, beta_inc
      beta = beta_start
10    continue

C     >>> prepare the names for configuration and data files
C     >>> - gauge field configuration is saved in b[beta]k[kappa].conf
C     >>> - results of measurements in b[beta]k[kappa].data
      ival = int((beta + 100.0000001d0) * 1000.0d0)
      write(c_beta, '(i6)') ival
      ival = int((akap + 1.000000001d0) * 100000.0d0)
      write(c_kappa, '(i6)') ival
      conf_name = 'b'//c_beta(2:6)//'k'//c_kappa(2:6)//'.conf'
      data_name = 'b'//c_beta(2:6)//'k'//c_kappa(2:6)//'.data'

C     >>> prepare textual data description
      c_computable(1) = 'acc'
      c_computable(2) = 'niter'
      c_computable(3) = 'eps'
      c_computable(4) = 'sgauge / nsites'
      c_computable(5) = 'topological charge'
      c_computable(6) = 'topological charge (sin)'
      c_computable(7) = 'beta * chiT'
      c_computable(8) = 'beta * chiT (sin)'

      open(unit=7,file=conf_name,form='unformatted',status='unknown')
      open(unit=8,file=data_name,form='unformatted',status='unknown')

      write(c_ident, '(a8)') '%data-v3'

C for the first beta, write the initialization data on the screen
      if(beta .eq. beta_start) then
        nw = 1
      else
        nw = 0
      end if
C save header in .data file
      call save_header(8, nw, c_ident, NSPACE, NTIME, beta,
     &  akap, eps, nsteps, ntherm, nmeas, mstep, istart, iseed)

C .data specific part - saving description of computables
      write(8) ncomp
      write(8) (c_computable(i), i = 1, ncomp)

C     ------------------------------------------------------------------
C     Equilibration - thermalization
C     ------------------------------------------------------------------
      do 100 it = 1, ntherm

	call gauge_met_c(beta)

c	>>> unclear if unitarization is necessary...
c	if(mod(it, 100) .eq. 0) call unitarize(u, ngaug)

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

      call cbl_cputime(0, tcpu)

      ns = 1
      it = 0
      cumtime = 0.0d0
200   it = it + 1

C       >>> update
	do i = 1, mstep
	  ii = (it - 1) * mstep + i

c	  >>> unclear if unitarization is necessary...
c	  if(mod(ii, 100) .eq. 0) call unitarize(u, ngaug)

	  call gauge_met_c(beta)

	end do
	
C       >>> measurements
	computable(1, ns) = float(iacc)
	computable(2, ns) = float(niter)
	computable(3, ns) = sngl(eps)
      computable(4, ns) = sngl(sgauge(beta, u) / dble(nsite))
C      print *, 'sgauge: ', sngl(sgauge(beta, u) / dble(nsite))
      computable(5, ns) = float(itopch())
C      print *, 'itopch: ', float(itopch())
      computable(6, ns) = sngl(sintopch())
      computable(7, ns) = sngl(beta * itopch()**2 / dble(nsite))
      computable(8, ns) = sngl(beta * sintopch()**2 / dble(nsite))
      
	if((ns .eq. nsave) .or. (it .eq. nmeas)) then
	  rewind 7
	  write(7) NTIME, NSPACE
          write(7) beta, akap, eps
          write(7) ntherm + it * mstep
	  write(7) u
	  endfile 7
	  do i = 1, nsave
	    write(8) (computable(j, i), j = 1, ncomp)
	  end do
	  ns = 0
	end if
	ns = ns + 1
      if(it .lt. nmeas) go to 200

      call cbl_cputime(1, tcpu)

      call save_tail(8, 0, 0, 0, 0, tcpu)

      close(8)

      beta = beta + beta_inc
      if(beta .le. beta_end + 0.00000001) goto 10

999   continue

      end

