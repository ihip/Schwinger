#define VLEN 21007

      real ran(VLEN)

      sum = 0.0
      call cbl_cputime(0,td)
      call cbl_rcarry(ran,VLEN,0.)
      do i=1,VLEN
        sum=sum+ran(i)
      enddo
      call cbl_cputime(1,td)
      write(*, *) 'rcarry - vector'
      print *,sum,td/VLEN

      sum = 0.0
      call cbl_cputime(0,td)
      do i=1,VLEN
        sum=sum+cbl_ranf()
      enddo
      call cbl_cputime(1,td)
      write(*, *) 'rcarry - one by one'
      print *,sum,td/VLEN
      end
      