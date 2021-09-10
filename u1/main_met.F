C     ******************************************************************
      program main_metropolis
C     ******************************************************************
C     **                    **                                        **
C     **  METROPOLIS        **    I. Hip, 29 Sep 98                   **
C     **                    **    Last modified: 25 Aug 21 (Hip)      **
C     **                    **                                        **
C     ******************************************************************
      parameter (ndrc=2,ngrp=1,ndim=2)
      parameter (nsite=NTIME*NSPACE)
      parameter (nferm=ndrc*nsite,ngaug=ndim*nsite)
      parameter (nall=2*ndrc*ngrp*nsite)

C     >>> ncomp - number of computables
      parameter (ncomp=5)
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
      character*4 c_nspace, c_ntime
      character*10 c_size
      character*6 c_beta, c_iconf
      character*27 conf_name, data_name
      character*8 c_ident
      character*1 c_saveconf

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
      write(*, *) 'Metropolis (29 Sep 98 / 25 Aug 21)'
      write(*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*, *)
      write(*, '(1x, ''beta = '', $)')
      read(*, *) beta
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
      if(istart .eq. 2) stop 'option 2 not yet implemented!'
c	  write(*, '(1x, ''start file = '', $)')
c	  read(*, '(a)') start_file
      write(*, '(1x, ''save configurations? [Y/n] '', $)')
      read(*, '(a)') c_saveconf
      write(*, *)

C     >>> initialize random number generator (seed is wall-time and host
C         dependent) - be sure to initialize it properly: there is no
C         internal check in rcarry anymore
      call cbl_rtim(iseed)
C      write(*, *) 'WARNING: seed is set to 38768509' 
C      iseed = 38768509
      call cbl_rc_seed(iseed)

C     >>> initialize gauge index arrays
      call gauge_mk_index()

C     >>> initialize indexes for sgauge & itopch
      call mk_index(ind)

C     >>> initialization of gauge fields
      call create_u(istart, u, start_file)

C     >>> configuration counter
      iconf = 1

C     >>> prepare the file name L[nspace]x[ntime]_b[beta]k[kappa].data
      write(c_nspace, '(i4)') NSPACE + 1000
      write(c_ntime, '(i4)') NTIME + 1000
      c_size = 'L'//c_nspace(2:4)//'x'//c_ntime(2:4)//'_b'
      ival = int((beta + 100.0000001d0) * 1000.0d0)
      write(c_beta, '(i6)') ival
      data_name = c_size//c_beta(2:6)//'k00000.data'

C     >>> prepare textual data description
      c_computable(1) = 'sgauge / nsites'
      c_computable(2) = 'topological charge'
      c_computable(3) = 'topological charge (sin)'
      c_computable(4) = 'beta * chiT'
      c_computable(5) = 'beta * chiT (sin)'

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
      do it = 1, ntherm
	  call gauge_met_c(beta)
c	  >>> unclear if unitarization is necessary...
c	  if(mod(it, 100) .eq. 0) call unitarize(u, ngaug)
      end do

C     ------------------------------------------------------------------
C     Measuring  updates
C     ------------------------------------------------------------------

      call cbl_cputime(0, tcpu)

c     >>> measurement counter
      it = 0
c     >>> counter for nsave
      ns = 0
c     >>> configuration counter
      iconf = 0
c     >>> timer
      cumtime = 0.0d0

c     loop from 1 to nmeas
200   it = it + 1

C     >>> update
	do i = 1, mstep  
c	  >>> unclear if unitarization is necessary...
c       ii = (it - 1) * mstep + i
c	  if(mod(ii, 100) .eq. 0) call unitarize(u, ngaug)
	  call gauge_met_c(beta)
	end do
      ns = ns + 1

C     >>> measurements
      computable(1, ns) = sngl(sgauge(beta, u) / dble(nsite))
      computable(2, ns) = float(itopch())
      computable(3, ns) = sngl(sintopch())
      computable(4, ns) = sngl(beta * itopch()**2 / dble(nsite))
      computable(5, ns) = sngl(beta * sintopch()**2 / dble(nsite))
 
	if(ns .eq. nsave) then
c       >>> save the computables in the .data file
        do i = 1, nsave
          write(8) (computable(j, i), j = 1, ncomp)
        end do
        if(c_saveconf .eq. 'Y') then
c         >>> save the gauge configuration in .conf file
          iconf = iconf + 1
          write(c_iconf, '(i6)') 100000 + iconf 
          conf_name = c_size//c_beta(2:6)//'i'//c_iconf(2:6)//'.conf'
          open(unit=7,file=conf_name,form='unformatted',status='new')
	      write(7) NTIME, NSPACE
            write(7) beta, akap, eps
            write(7) ntherm + it * mstep
	      write(7) u
	    close(7)
        end if

	  ns = 0
	end if

      if(it .lt. nmeas) go to 200

      call cbl_cputime(1, tcpu)

      call save_tail(8, 0, 0, 0, 0, tcpu)

      close(8)

999   continue

      end

