#define MAX_POINTS 100

	character*64 infile
	real x(MAX_POINTS), y(MAX_POINTS)
	real x_sig(MAX_POINTS), y_sig(MAX_POINTS)
	real al, bk, sig_al, sig_bk, chi2, q, sig

	write(*, *) 'PCACFIT X'
	write(*, *)
	write(*, '(1x, ''Input file (e.g. pcac.col): '', $)')
	read(*, '(a)') infile
	write(*, *)

	open(1, file=infile, status='old')

	do i = 1, MAX_POINTS
	  y_sig(i) = 0.0
	end do

	i = 1
10	read(1, *, end = 99) y(i), x(i), x_sig(i)
	i = i + 1
	if(i .gt. MAX_POINTS) stop 'error: too many data points!'
	goto 10

99	close(1)
	n = i - 1
	write(*, *) 'Data points read: ', n
	write(*, *)
	write(*, *) 'Linear fit: y = k x + l'

	call fitexy(x, y, n, x_sig, y_sig, al, bk, sig_al, sig_bk,
     & chi2, q)

	write(*, *) '  k = ', bk, ' +/-', sig_bk
	write(*, *) '  l = ', al, ' +/-', sig_al
	write(*, *) 
	write(*, *) '  chi2 = ', chi2
	write(*, *) '  q = ', q, '   (q > 0.1 good, q < 0.001 poor)'
	write(*, *)
	write(*, *) 'K_c = ', al, ' +/- ', sig_al
	write(*, *)
	write(*, *) al, sig_al

	end



	
