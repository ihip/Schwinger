C v3: 07 Jan 97, ivh - Last modified: Hip, 2022-04-12

      subroutine save_header(ifile, nw, c_ident, nspace, ntime, beta,
     &  fmass, eps, nsteps, ntherm, nmeas, mstep, istart, iseed)
C     ------------------------------------------------------------------
      character*8 c_ident
      real*8 beta, fmass, eps
C     ------------------------------------------------------------------
      character*64 c_line
      character*32 c_integrate, c_invert, c_invert_b, c_invert1,
     &  c_invert1_b, c_future, hostname
      character*8  c_date, c_time
      integer*4 ifut(8)
C     ------------------------------------------------------------------
      character*32 v_integrate, v_invert, v_invert_b,
     &  v_invert1, v_invert1_b
C     ------------------------------------------------------------------

C     >>> prepare human-readable line
      write(c_line, 1000) nspace, ntime, beta, fmass, nmeas
1000  format(1x, i3, ' * ', i3, '  beta = ', f6.3, ' kappa = ', f6.3,
     &  '  nmeas = ', i7)

      c_integrate = v_integrate()
      c_invert = v_invert()
      c_invert_b = v_invert_b()
      c_invert1 = v_invert1()
      c_invert1_b = v_invert1_b()

      call cbl_datim(c_date, c_time)

C     >>> WARNING: nonportable code !!!
      ierr = gethostname(hostname, 32)

      write(ifile) c_ident, c_line,
     &  c_integrate, c_invert, c_invert_b, c_invert1,
     &  c_invert1_b, c_future, c_date, c_time, hostname,
     &  nspace, ntime, beta, fmass, eps, nsteps, ntherm, nmeas,
     &  mstep, istart, iseed, ifut

      if(nw .eq. 1) then
        write(*, 888)
        write(*, *)
        write(*, *) '                           H M C 2 D F U 1 (v4)'
        write(*, *)
        write(*, 888)
        write(*, 1001) nspace, ntime
        write(*, 888)
        write(*, 1002) c_integrate
        write(*, 1003) c_invert, c_invert_b
        write(*, 1004) c_invert1, c_invert1_b
        write(*, 888)
        write(*, 1005) beta
        write(*, 1006) fmass
        write(*, 1007) eps
        write(*, 1008) nsteps
        write(*, 1009) ntherm
        write(*, 1010) nmeas
        write(*, 1011) mstep
        write(*, 1012) istart
        write(*, 1013) iseed
        write(*, 888)
888     format(1x, 13('******'))
1001    format(' 2 Flavor Schwinger Model on a ',i3,' x',i3,' Lattice')
1002    format(' Integration: ', a32)
1003    format(' Inversion: ', a32, ' (', a32, ')')
1004    format(' Invers(1): ', a32, ' (', a32, ')')
1005    format(' beta                                = ', f10.5)
1006    format(' fmass                               = ', f10.5)
1007    format(' eps                                 = ', f10.5)
1008    format(' Length of Traj.                     = ', i10)
1009    format(' Thermalization                      = ', i10)
1010    format(' Measurements                        = ', i10)
1011    format(' Measurement step                    = ', i10)
1012    format(' Start (0 = COLD, 1 = HOT, 2 = FILE) = ', i10)
1013    format(' Seed                                = ', i10)
        write(*, *) 'Now you can press "CTRL+Z" and type "bg <RET>"...'
      end if

      return
      end

      subroutine save_tail(ifile, nexcp, nexcp_b, nexcp1, nexcp1_b,
     &  tcpu)
      integer*4 nexcp, nexcp_b, nexcp1, nexcp1_b, itfut(4)
      real*4 tcpu
      write(ifile) nexcp, nexcp_b, nexcp1, nexcp1_b, itfut, tcpu
      return
      end

c >>> gauge configurations generated using Metropolis
      character*32 function v_integrate()
      v_integrate = 'Metropolis'
      return
      end
      
c >>> dummy definition of v_invert
      character*32 function v_invert()
      v_invert = 'none'
      return
      end

c >>> dummy definition of v_invert_b
      character*32 function v_invert_b()
      v_invert_b = 'none'
      return
      end

c >>> LAPACK inversion subroutine(s)
      character*32 function v_invert1()
      v_invert1 = 'LAPACK: zgetrf + zgetri'
      return
      end

c >>> dummy definition of v_invert1_b
      character*32 function v_invert1_b()
      v_invert1_b = 'none'
      return
      end



