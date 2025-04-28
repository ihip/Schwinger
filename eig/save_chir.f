	  subroutine save_chir(n, eig, chir, conf_file, dirop)
c     IN integer n - number of eigenvalues
c     IN complex*16 eig(n) - eigenvalues
c     IN character*64 conf_file
c     IN character*2 abbreviation for Dirac operator (e.g.'HF')
c     >>> ihip, 05 Aug 20; last modified: 06 Aug 20
	  complex*16 eig(n)
	  character*64 conf_file, eigname
	  character*2 dirop
	  eigname = conf_file(1:lnblnk(conf_file)-5)//'-'//dirop//'.eig'
	  open(21, file = eigname, status = 'unknown')
	    do i = 1, n
	      write(21, '(2g24.15)') dreal(eig(i)), dimag(eig(i)), chir(i)
	    end do
	  close(21)
	  write(*, '(1x, ''Eigenvalues written to: '', a)') eigname
	  write(*, *)
	  end
