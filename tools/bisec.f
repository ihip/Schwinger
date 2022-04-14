        subroutine massb(ntime, corr, efmass)

c	IN ntime - lattice dimension in time direction
c	IN (real*8) corr - correlation data
c       OUT (real*8) efmass - effective masses

        real*8 corr(ntime - 1), efmass(ntime - 2)

        real*8 r12, bisectm
        
c	>>> make symmetrization
	do i = 1, ntime / 2 - 1
	  corr(i) = (corr(i) + corr(ntime - i)) / 2.0d0
	end do

        do i = 1, ntime - 2
	  if((corr(i) .lt. 0) .or. (corr(i + 1) .lt. 0)) then
c	    write(*, *) 'neg. corr'
	    efmass(i) = -100.0d0
	  else

          r12 = corr(i) / corr(i + 1)
          nt1 = i - ntime / 2
          nt2 = i + 1 - ntime / 2
	  if(((i .lt. ntime / 2) .and. (r12 .ge. 1.0d0)) .or.
     &       ((i .ge. ntime / 2) .and. (r12 .le. 1.0d0))) then
            efmass(i) =
     &      - dlog(bisectm(nt1, nt2, r12, 10.0d0, 1.0d-6, 1.0d-9))
	  else
            efmass(i) =
     & dlog(bisectm(nt1, nt2, 1.0d0 / r12, 10.0d0, 1.0d-6, 1.0d-9))
c	write(*, *) '### ', efmass(i)
	  end if
	  end if	    
        end do

        return
        end

        
        real*8 function bisectm(nt1, nt2, r12, aa, bb, epsilon)
        real*8 r12, aa, bb, epsilon
        
        real*8 x, ff, a, b, massfunc

	a = dexp(-aa)
	b = dexp(-bb)
        if( massfunc(a, nt1, nt2, r12) * 
     &    massfunc(b, nt1, nt2, r12) .gt. 0.0d0 ) then
          write(*, *) 'Unproper starting interval'
	  write(*, *) nt1, nt2, r12
          bisectm = -1.0d0
          return
        end if

        x = (a + b) / 2.0d0
        ff = massfunc(x, nt1, nt2, r12)
        
	i = 0
        do while( (dabs(ff) .gt. epsilon) .and. (i .lt. 200) )
          if( massfunc(a, nt1, nt2, r12) * ff .gt. 0.0d0 ) then
            a = x
          else
            b = x
          end if
          x = (a + b) / 2.0d0
          ff = massfunc(x, nt1, nt2, r12)
	  i = i + 1
        end do
c	if(i .eq. 200) write(*, *) r12, i
        
        bisectm = x
        return
        end

        real*8 function massfunc(x, nt1, nt2, r12)
        real*8 x, r12

        massfunc = x**nt1 + x**(-nt1) - r12 * (x**nt2 + x**(-nt2))
        return
        end


