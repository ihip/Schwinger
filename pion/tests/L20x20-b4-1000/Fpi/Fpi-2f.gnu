set terminal pdfcairo
set size ratio 0.75
set output 'Fpi-2f-L20x20-b4-1000-j1-chi2-pl4-jk20.pdf'

set title "F_{/Symbol p} / 2-flavor / j1 / L20x20 / {/Symbol b}=4 / 1000 meas."

set xlabel 'm'
set ylabel 'F_{/Symbol p}'

set key right bottom

plot [0:0.51][0:0.61] \
  0.398942 lc "black" dashtype "." title "F_{/Symbol p} = 1 / sqrt(2 * {/Symbol p})", \
  "../pion-L20x20-b4-1000-j1-chi2-pl4-jk20-2f-HOR.col" u 1:6:7 w yerr pt 13 lc "black" \
    ps 0.4 lw 0.3 title "reweighted hypercube overlap (N_f=2)"

quit
