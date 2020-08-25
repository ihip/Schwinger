#define MAX_POINTS 100

	character*64 infile
	real x(MAX_POINTS), y(MAX_POINTS), y_sig(MAX_POINTS)
	real al, bk, sig_al, sig_bk, chi2, q, sig

	write(*, *) 'PCACFIT'
	write(*, *)
	write(*, '(1x, ''Input file (e.g. pcac.col): '', $)')
	read(*, '(a)') infile
	write(*, *)

	open(1, file=infile, status='old')

	i = 1
10	read(1, *, end = 99) x(i), y(i), y_sig(i)
	i = i + 1
	if(i .gt. MAX_POINTS) stop 'error: too many data points!'
	goto 10

99	close(1)
	n = i - 1
	write(*, *) 'Data points read: ', n
	write(*, *)
	write(*, *) 'Linear fit: y = k x + l'

	call fit(x, y, n, y_sig, 1, al, bk, sig_al, sig_bk, chi2, q)

	write(*, *) '  k = ', bk, ' +/-', sig_bk
	write(*, *) '  l = ', al, ' +/-', sig_al
	write(*, *) 
	write(*, *) '  chi2 = ', chi2
	write(*, *) '  q = ', q, '   (q > 0.1 good, q < 0.001 poor)'
	write(*, *)
	sig = sqrt((al * sig_bk / bk)**2 / bk**2 +
     &   (sig_al / bk)**2)
	write(*, *) 'K_c = ', -al / bk, ' +/- ', sig
	write(*, *)
	write(*, *) -al / bk, sig

	end



	
