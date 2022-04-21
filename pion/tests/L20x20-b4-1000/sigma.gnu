set terminal pdfcairo
set size ratio 0.75
set output 'sigma-L20x20-b4-1000-jk20.pdf'

set title "chiral condensate / L20x20 / {/Symbol b}=4 / 1000 meas."

set xlabel 'm'
set ylabel '{/Symbol S}'

set key right top

plot [0:0.51][0:0.36] \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-0f-HOM.col" u 1:4:5 w yerr pt 7 lc "orange" ps 0.3 lw 0.3 \
     title "quenched hypercube overlap (N_f = 0)", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-2f-HOR.col" u 1:4:5 w yerr pt 13 lc "red" ps 0.4 lw 0.3 \
     title "reweighted hypercube overlap (N_f=2)", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-4f-HOR.col" u 1:4:5 w yerr pt 13 lc "green" ps 0.4 lw 0.3 \
     title "reweighted hypercube overlap (N_f=4)"

quit
